// GPL-3.0-or-later. Compact FP64 graph derived from the frozen compact_math.h.
// New namespace and explicit unit map; old diagnostic implementation untouched.
#ifndef MECHSYS_BENCHMARK_MATH_H
#define MECHSYS_BENCHMARK_MATH_H
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#ifdef __CUDACC__
#define CB_HD __host__ __device__
#else
#define CB_HD
#endif
namespace MechSysBenchmarkMath {
struct Units {
    double dx,dt,conversion,C,mass,inertia,impulse,angular_impulse,force,torque;
    Units(double x,double t,double f):dx(x),dt(t),conversion(f) {
        if(!(std::isfinite(x)&&x>0&&std::isfinite(t)&&t>0&&std::isfinite(f)&&f>0))
            throw std::runtime_error("CB: positive finite dx/dt/Fconv required");
        C=x/t;mass=f*x*x*x;inertia=mass*x*x;impulse=mass*C;
        angular_impulse=impulse*x;force=impulse/t;torque=angular_impulse/t;
        for(double v:{C,mass,inertia,impulse,angular_impulse,force,torque})
            if(!(std::isfinite(v)&&v>0))throw std::runtime_error("CB: unit scale overflow/underflow");
    }
};
struct Trace {
    unsigned n=0;int error=0;
    CB_HD void reset(){n=0;error=0;}
    CB_HD bool range(double x){return x==0||(::fabs(x)>=0x1p-450&&::fabs(x)<=0x1p450);}
    CB_HD double put(double z,double a,double b){if(!range(a)||!range(b)||!range(z))error=1;if(n<2048)++n;else error=2;return z;}
};
CB_HD inline double add(Trace &t,double a,double b){
#ifdef __CUDA_ARCH__
    double z=__dadd_rn(a,b);
#else
    volatile double z=a+b;
#endif
    return t.put(z,a,b);
}
CB_HD inline double sub(Trace &t,double a,double b){return add(t,a,-b);}
CB_HD inline double mul(Trace &t,double a,double b){
#ifdef __CUDA_ARCH__
    double z=__dmul_rn(a,b);
#else
    volatile double z=a*b;
#endif
    if(a!=0&&b!=0&&z==0)t.error=3;return t.put(z,a,b);
}
CB_HD inline double div(Trace &t,double a,double b){
    if(!(b>0)){t.error=4;return t.put(0,a,b);}
#ifdef __CUDA_ARCH__
    double z=__ddiv_rn(a,b);
#else
    volatile double z=a/b;
#endif
    if(a!=0&&z==0)t.error=3;return t.put(z,a,b);
}
CB_HD inline double root(Trace &t,double a){
    if(!(a>0)){t.error=4;return t.put(0,a,0);}
#ifdef __CUDA_ARCH__
    double z=__dsqrt_rn(a);
#else
    volatile double z=::sqrt(a);
#endif
    return t.put(z,a,0);
}
CB_HD inline void cross(Trace &t,double const *a,double const *b,double *out){for(int i=0;i<3;i++){int j=(i+1)%3,k=(i+2)%3;double x=mul(t,a[j],b[k]),y=mul(t,a[k],b[j]);out[i]=sub(t,x,y);}}
CB_HD inline void quat(Trace &t,double const *a,double const *b,double *out){
    double dot=0;for(int i=1;i<4;i++){double p=mul(t,a[i],b[i]);dot=add(t,dot,p);}double p=mul(t,a[0],b[0]);out[0]=sub(t,p,dot);
    for(int i=0;i<3;i++){int j=(i+1)%3,k=(i+2)%3;double x=mul(t,a[0],b[i+1]),y=mul(t,b[0],a[i+1]);double s=add(t,x,y);x=mul(t,a[j+1],b[k+1]);y=mul(t,a[k+1],b[j+1]);double c=sub(t,x,y);out[i+1]=add(t,s,c);}
}
CB_HD inline void rotate(Trace &t,double const *q,double const *v,double *out,bool inverse=false){
    double a[4]={q[0],inverse?-q[1]:q[1],inverse?-q[2]:q[2],inverse?-q[3]:q[3]},b[4]={0,v[0],v[1],v[2]},conj[4]={a[0],-a[1],-a[2],-a[3]},tmp[4],r[4];
    quat(t,a,b,tmp);quat(t,tmp,conj,r);for(int i=0;i<3;i++)out[i]=r[i+1];
}
struct MathCell {
    int owner=-1;double A[18]{},p[3]{},rho=0,D[36]{},b[6]{},j[3]{},kw[3]{},kb[3]{};Trace trace;
    CB_HD void clear(){owner=-1;rho=0;trace.reset();for(auto &x:A)x=0;for(auto &x:D)x=0;for(auto &x:b)x=0;for(int i=0;i<3;i++)p[i]=j[i]=kw[i]=kb[i]=0;}
};
struct Body {
    double M[6]{},z[6]{},Q[4]{},D[36]{},C[36]{},b[6]{},L[36]{},y[6]{},delta[6]{},t[6]{};unsigned owned=0;Trace trace;
    CB_HD void clear(){owned=0;trace.reset();for(int i=0;i<36;i++)D[i]=C[i]=L[i]=0;for(int i=0;i<6;i++)M[i]=z[i]=b[i]=y[i]=delta[i]=t[i]=0;for(auto &x:Q)x=0;}
};
CB_HD inline void assemble(MathCell &m,Body const &s,double const *r,double B,double const *f,double rho,int const c[15][3]){
    Trace &t=m.trace;m.rho=rho;for(int i=0;i<3;i++)for(int j=0;j<3;j++)m.A[6*i+j]=(i==j);
    for(int j=0;j<3;j++){double e[3]={0,0,0},w[3],v[3];e[j]=1;rotate(t,s.Q,e,w);cross(t,w,r,v);for(int i=0;i<3;i++)m.A[6*i+3+j]=v[i];}
    for(int i=0;i<3;i++){m.p[i]=0;for(int k=0;k<15;k++){double v=c[k][i]==0?0:c[k][i]>0?f[k]:-f[k];m.p[i]=add(t,m.p[i],v);}}
    double a=2*B,w=mul(t,a,rho);
    for(int j=0;j<6;j++)for(int k=j;k<6;k++){double v=0;for(int i=0;i<3;i++){double p=mul(t,m.A[6*i+j],m.A[6*i+k]);v=add(t,v,p);}m.D[6*j+k]=mul(t,w,v);m.D[6*k+j]=m.D[6*j+k];}
    double u[3]={},rhs[3];for(int i=0;i<3;i++)for(int j=0;j<6;j++){double p=mul(t,m.A[6*i+j],s.z[j]);u[i]=add(t,u[i],p);}
    for(int i=0;i<3;i++){double p=mul(t,rho,u[i]);rhs[i]=sub(t,m.p[i],p);}
    for(int j=0;j<6;j++){double v=0;for(int i=0;i<3;i++){double p=mul(t,m.A[6*i+j],rhs[i]);v=add(t,v,p);}m.b[j]=mul(t,a,v);}
}
CB_HD inline void solve(Body &s){
    Trace &t=s.trace;for(int j=0;j<6;j++)for(int k=j;k<6;k++){if(j!=k)s.D[6*k+j]=s.D[6*j+k];s.C[6*j+k]=s.C[6*k+j]=s.D[6*j+k];}
    for(int j=0;j<6;j++)s.C[6*j+j]=add(t,s.M[j],s.D[6*j+j]);
    for(int j=0;j<6;j++)for(int i=j;i<6;i++){double v=s.C[6*i+j];for(int k=0;k<j;k++){double p=mul(t,s.L[6*i+k],s.L[6*j+k]);v=sub(t,v,p);}s.L[6*i+j]=i==j?root(t,v):div(t,v,s.L[6*j+j]);if(t.error)return;}
    for(int i=0;i<6;i++){double v=s.b[i];for(int k=0;k<i;k++){double p=mul(t,s.L[6*i+k],s.y[k]);v=sub(t,v,p);}s.y[i]=div(t,v,s.L[6*i+i]);}
    for(int i=5;i>=0;i--){double v=s.y[i];for(int k=i+1;k<6;k++){double p=mul(t,s.L[6*k+i],s.delta[k]);v=sub(t,v,p);}s.delta[i]=div(t,v,s.L[6*i+i]);}
    for(int i=0;i<6;i++)s.t[i]=add(t,s.z[i],s.delta[i]);
}
CB_HD inline double equilibrium(Trace &t,int k,double rho,double const *u,double W,int const c[15][3]){
    double cu=0,uu=0;for(int i=0;i<3;i++){double p=mul(t,c[k][i],u[i]);cu=add(t,cu,p);p=mul(t,u[i],u[i]);uu=add(t,uu,p);}
    double l=mul(t,3,cu),sq=mul(t,cu,cu),q=mul(t,4.5,sq),v=mul(t,1.5,uu),a=add(t,1,l);a=add(t,a,q);a=sub(t,a,v);double wr=mul(t,W,rho);return mul(t,wr,a);
}
CB_HD inline void surface_omega(MathCell &m,Body const &s,double const *r,double *surface,double *omega,double const *f,double rho,double const *W,int const c[15][3],int const *op){
    double w[3],v[3];rotate(m.trace,s.Q,s.t+3,w);cross(m.trace,w,r,v);for(int i=0;i<3;i++)surface[i]=add(m.trace,s.t[i],v[i]);
    for(int k=0;k<15;k++){double a=equilibrium(m.trace,op[k],rho,surface,W[op[k]],c),b=equilibrium(m.trace,k,rho,surface,W[k],c);a=sub(m.trace,f[op[k]],a);b=sub(m.trace,f[k],b);omega[k]=sub(m.trace,a,b);}
}
CB_HD inline void reaction(MathCell &m,Body const &s,double const *r,double const *solid,int const c[15][3]){
    for(auto &x:m.j)x=0;for(int k=0;k<15;k++)for(int i=0;i<3;i++){double p=mul(m.trace,solid[k],c[k][i]);m.j[i]=sub(m.trace,m.j[i],p);}
    cross(m.trace,r,m.j,m.kw);rotate(m.trace,s.Q,m.kw,m.kb,true);
}
} // namespace MechSysBenchmarkMath
#endif
