// FP64 exact expansions for a bounded sphere/cube geometry.
// Rounding-error compensation helpers used by the experimental coupled benchmark
// path; the bound they enforce is documented at their point of use.
#ifndef MECHSYS_SPHERE_CUBE_ROBUST_H
#define MECHSYS_SPHERE_CUBE_ROBUST_H
#include <cmath>
#include <cfloat>
#ifdef __CUDACC__
#define MECHSYS_EDGE_HD __host__ __device__
#else
#define MECHSYS_EDGE_HD
#endif
#define MECHSYS_ROBUST_EDGE_PROXY_AVAILABLE 1
namespace MechSysRobustEdge {
MECHSYS_EDGE_HD inline double add(double a,double b) {
#ifdef __CUDA_ARCH__
    return __dadd_rn(a,b);
#else
    volatile double r=a+b;return r;
#endif
}
MECHSYS_EDGE_HD inline double mul(double a,double b) {
#ifdef __CUDA_ARCH__
    return __dmul_rn(a,b);
#else
    volatile double r=a*b;return r;
#endif
}
MECHSYS_EDGE_HD inline double div(double a,double b) {
#ifdef __CUDA_ARCH__
    return __ddiv_rn(a,b);
#else
    volatile double r=a/b;return r;
#endif
}
MECHSYS_EDGE_HD inline double fused(double a,double b,double c) {
#ifdef __CUDA_ARCH__
    return __fma_rn(a,b,c);
#else
    return std::fma(a,b,c);
#endif
}
MECHSYS_EDGE_HD inline void two_sum(double a,double b,double &h,double &l) {
    h=add(a,b);double bv=add(h,-a),av=add(h,-bv);
    l=add(add(a,-av),add(b,-bv));
}
// Grow-expansion-zero-elim (general TwoSum, no magnitude-order assumption).
// Shewchuk's expansion arithmetic; independent implementation for this small
// fixed-degree expression. Max capacity proven from the number of terms.
struct Expansion {
    double v[64];int n;
    MECHSYS_EDGE_HD Expansion():n(0) {}
    MECHSYS_EDGE_HD void grow(double b) {
        if(b==0.0)return;
        double q=b;int m=0;
        for(int i=0;i<n;i++) {
            double h,l;two_sum(q,v[i],h,l);if(l!=0.0)v[m++]=l;q=h;
        }
        if(q!=0.0)v[m++]=q;n=m;
    }
    MECHSYS_EDGE_HD int sign() const {return n?(v[n-1]>0.0?1:-1):0;}
    MECHSYS_EDGE_HD double value() const {
        double q=0.0;for(int i=0;i<n;i++)q=add(q,v[i]);return q;
    }
};
MECHSYS_EDGE_HD inline Expansion coordinate(double c,double x,double h) {
    Expansion q;q.grow(c);q.grow(-x);q.grow(h);return q;
}
MECHSYS_EDGE_HD inline void product(Expansion &s,double a,double b,double sign) {
    double hi=mul(a,b),lo=fused(a,b,-hi);s.grow(sign*lo);s.grow(sign*hi);
}
MECHSYS_EDGE_HD inline void square(Expansion &s,Expansion const &q,double sign) {
    for(int i=0;i<q.n;i++)for(int j=0;j<q.n;j++)product(s,q.v[i],q.v[j],sign);
}
MECHSYS_EDGE_HD inline Expansion endpoint(Expansion m,Expansion const &q) {square(m,q,-1.0);return m;}
struct Trace {double margin_hi=0.0,margin_lo=0.0,root=0.0;int sign=0,components=0,scale=0;};
// This conservative exponent envelope guarantees exact transforms and normal
// products/residuals. Outside it, NaN/false is explicit unsupported input, NOT
// a zero intersection. No ambiguous sign is replaced by zero.
MECHSYS_EDGE_HD inline bool normalize(double const *x,double const *c,double R,double dx,
                                       double *xn,double *cn,double &rn,double &dn,int &exponent) {
    if(!(R>=0.0) || !(dx>0.0))return false;
    double maximum=R>dx?R:dx;
    for(int k=0;k<3;k++) {
        double a=::fabs(x[k]),b=::fabs(c[k]);
        if(!(a<=DBL_MAX) || !(b<=DBL_MAX))return false;
        if(a>maximum)maximum=a;if(b>maximum)maximum=b;
    }
    if(!(maximum<=DBL_MAX))return false;
    exponent=::ilogb(maximum);if(exponent< -400 || exponent>400)return false;
    rn=::scalbn(R,-exponent);dn=::scalbn(dx,-exponent);
    // Least input bit >=2^-452 (half dx >=2^-453), hence every
    // exact quadratic product term >=2^-906: safely normal binary64.
    const double minimum=0x1p-400;
    if((rn!=0.0 && rn<minimum)||dn<minimum)return false;
    for(int k=0;k<3;k++) {
        xn[k]=::scalbn(x[k],-exponent);cn[k]=::scalbn(c[k],-exponent);
        if((xn[k]!=0.0&&::fabs(xn[k])<minimum)||(cn[k]!=0.0&&::fabs(cn[k])<minimum))return false;
    }
    return true;
}
MECHSYS_EDGE_HD inline bool supported(double const *x,double const *c,double R,double dx) {
    double xn[3],cn[3],rn,dn;int e;return normalize(x,c,R,dx,xn,cn,rn,dn,e);
}
MECHSYS_EDGE_HD inline double length(double const *x,double const *c,double R,double dx,
                                     double *edges=nullptr,Trace *trace=nullptr) {
    double xn[3],cn[3],rn,dn;int exponent;
    if(!normalize(x,c,R,dx,xn,cn,rn,dn,exponent))return NAN;
    const int signs[8][3]={{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1},
                           {-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1}};
    Expansion total;
    for(int e=0;e<12;e++) {
        int j=e%4,vertex=e<4?j:(e<8?j+4:j),axis=e<8?j%2:2;
        Expansion m;product(m,rn,rn,1.0);
        for(int k=0;k<3;k++)if(k!=axis) {
            Expansion q=coordinate(cn[k],xn[k],mul(0.5*dn,signs[vertex][k]));square(m,q,-1.0);
        }
        double L=0.0,root=0.0;
        if(m.sign()>0) {
            root=::sqrt(m.value());
            Expansion lo=coordinate(cn[axis],xn[axis],-0.5*dn),hi=coordinate(cn[axis],xn[axis],0.5*dn);
            Expansion ml=endpoint(m,lo),mh=endpoint(m,hi);
            if(!((lo.sign()>=0&&ml.sign()<=0)||(hi.sign()<=0&&mh.sign()<=0))) {
                bool inside_lo=lo.sign()>=0||ml.sign()>=0,inside_hi=hi.sign()<=0||mh.sign()>=0;
                if(inside_lo&&inside_hi)L=dn;
                else if(!inside_lo&&!inside_hi)L=mul(2.0,root);
                else if(inside_lo) {
                    if(lo.sign()>=0)L=div(ml.value(),add(root,lo.value()));
                    else L=add(root,-lo.value());
                } else {
                    if(hi.sign()<=0)L=div(mh.value(),add(root,-hi.value()));
                    else L=add(root,hi.value());
                }
                // Only a last-bit overshoot of a mathematically bounded length
                // can reach this bound. No positive minimum and no zero clamp.
                if(L>dn)L=dn;
                if(!(L>0.0))return NAN;
            }
        }
        double scaled=::scalbn(L,exponent);
        if(L>0.0 && scaled==0.0)return NAN; // unrepresentable, never silently zero
        L=scaled;
        if(edges)edges[e]=L;total.grow(L);
        if(trace) {
            trace[e].sign=m.sign();trace[e].components=m.n;trace[e].scale=exponent;trace[e].root=root;
            trace[e].margin_hi=m.n?m.v[m.n-1]:0.0;
            Expansion tail;for(int k=0;k<m.n-1;k++)tail.grow(m.v[k]);trace[e].margin_lo=tail.value();
        }
    }
    double result=total.value(),bound=12.0*dx;return result>bound?bound:result;
}
} // namespace MechSysRobustEdge
#undef MECHSYS_EDGE_HD
#endif
