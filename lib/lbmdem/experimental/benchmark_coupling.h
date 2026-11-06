// GPL-3.0-or-later. Default-off single-fluid benchmark API, NOT production Solve.
#ifndef MECHSYS_EXPERIMENTAL_BENCHMARK_COUPLING_H
#define MECHSYS_EXPERIMENTAL_BENCHMARK_COUPLING_H
#ifdef MECHSYS_EXPERIMENTAL_BENCHMARK_COUPLING
#if !defined(USE_CUDA) || !defined(__CUDACC__)
#error "Benchmark coupling requires CUDA compilation and USE_CUDA"
#endif
#include <mechsys/lbmdem/Domain.h>
#include <mechsys/dem/sphere_cube_robust.h>
#include "detail/benchmark_math.h"
#include "detail/density_sum.h"
#include <array>
#include <algorithm>
#include <climits>
#include <cstring>
#include <functional>
#include <set>
#include <vector>
namespace LBMDEM { namespace Experimental { namespace Benchmark {
namespace M=MechSysBenchmarkMath;
inline void need(bool b,char const *s){if(!b)throw std::runtime_error(s);}
inline void ck(cudaError_t e){if(e!=cudaSuccess)throw std::runtime_error(cudaGetErrorString(e));}
template<class T>inline T *ptr(thrust::device_vector<T> &v){return thrust::raw_pointer_cast(v.data());}
template<class T>inline std::vector<T> read(T const *p,size_t n){std::vector<T> v(n);if(n)ck(cudaMemcpy(v.data(),p,n*sizeof(T),cudaMemcpyDeviceToHost));return v;}
// Keep recorder files cell-contiguous while using upstream's device layout.
inline std::vector<double> read_cell_populations(double const *p,size_t cells){
    auto raw=read(p,15*cells);std::vector<double> out(raw.size());
    for(size_t i=0;i<cells;i++)for(size_t k=0;k<15;k++)
        out[15*i+k]=raw[FLBM::FIDX(cells,15,0,k,i)];
    return out;
}
enum class Selection {Disabled,ExplicitPV,ImplicitVV};
enum class DensitySum {SequentialV1,CompensatedV1};
struct Options {
    Selection selection=Selection::Disabled;
    DensitySum density_sum=DensitySum::SequentialV1;
    size_t max_intervals=0,max_cells=16000000,max_claim_bytes=size_t(4)*1024*1024*1024;
    Vec3_t fluid_acceleration=OrthoSys::O;
};
struct Claim {
    int owner=-1,cell=-1,next=-1;
    double gamma=0,B=0,lever[3]{},surface[3]{},omega[15]{},solid[15]{};
    M::MathCell math;
};
struct Cell {int head=-1,count=0;double gamma_sum=0,Gamma=0,B=0;};
struct Impulse {double J[3]{},Kw[3]{},Kb[3]{};};
namespace Detail {
// The small algebra routines consume one cell's 15 populations in k order.
// Gathering changes storage addressing only; arithmetic and ordering are retained.
__device__ inline void load_populations(double const *f,size_t cell,
                                       FLBM::lbm_aux const *la,double (&out)[15]){
    for(int k=0;k<15;k++)out[k]=f[FLBM::FIDX(la->Ncells,15,0,k,cell)];
}

__device__ inline bool prescribed(DEM::ParticleCU const &p){return p.vxf&&p.vyf&&p.vzf&&p.wxf&&p.wyf&&p.wzf;}
__global__ void claims(Claim *a,Cell *cells,size_t cap,int reach,DEM::ParticleCU const *p,
    DEM::DynParticleCU const *d,DEM::dem_aux const *da,FLBM::lbm_aux const *la,bool const *mask,int *error){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=cap)return;
    int side=2*reach+1;size_t per=size_t(side)*side*side,ip=i/per,t=i%per;
    auto &c=a[i];c.owner=-1;c.cell=-1;c.math.clear();
    double xx[3]={d[ip].x.x/la->dx,d[ip].x.y/la->dx,d[ip].x.z/la->dx};
    for(double x:xx)if(!(x>double(INT_MIN)+side&&x<double(INT_MAX)-side)){atomicOr(error,1);return;}
    int q[3]={int(floor(xx[0]))+int(t%side)-reach,int(floor(xx[1]))+int((t/side)%side)-reach,int(floor(xx[2]))+int(t/(side*side))-reach};
    int dim[3]={int(la->Nx),int(la->Ny),int(la->Nz)};double peraxis[3]={da->Per.x,da->Per.y,da->Per.z};
    for(int j=0;j<3;j++){if(q[j]<0||q[j]>=dim[j]){if(!(peraxis[j]>0))return;q[j]=(q[j]%dim[j]+dim[j])%dim[j];}}
    real3 center=make_real3(q[0]*la->dx,q[1]*la->dx,q[2]*la->dx),r;
    DEM::BranchVec(d[ip].x,center,r,da->Per);real3 image=center-r;
    double x[3]={image.x,image.y,image.z},z[3]={center.x,center.y,center.z};
    if(!MechSysRobustEdge::supported(x,z,p[ip].R,la->dx)){atomicOr(error,2);return;}
    double length=MechSysRobustEdge::length(x,z,p[ip].R,la->dx);
    if(!(length>=0&&length<=12*la->dx)){atomicOr(error,4);return;}
    if(fabs(length/la->dx)<1e-12)return;
    c.cell=q[0]+la->Nx*(q[1]+la->Ny*q[2]);
    if(mask[c.cell]){atomicOr(error,128);return;} // no silent particle/wall-cell overlap
    c.owner=int(ip);c.gamma=length/(12*la->dx);
    c.lever[0]=r.x/la->dx;c.lever[1]=r.y/la->dx;c.lever[2]=r.z/la->dx;
    c.next=atomicExch(&cells[c.cell].head,int(i));atomicAdd(&cells[c.cell].count,1);
}
__global__ void partition(Claim *a,Cell *cells,size_t nc,double tau,int *error){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=nc)return;auto &c=cells[i];double g=0;int seen=0;
    for(int j=c.head;j!=-1;j=a[j].next){g=__dadd_rn(g,a[j].gamma);seen++;}
    if(seen!=c.count||!(g>=0&&g<=DBL_MAX)){atomicOr(error,8);return;}
    c.gamma_sum=g;c.Gamma=fmin(1.,g);c.B=c.Gamma*(tau-.5)/(1-c.Gamma+tau-.5);
    for(int j=c.head;j!=-1;j=a[j].next)a[j].B=__dmul_rn(c.B,__ddiv_rn(a[j].gamma,g));
}
__global__ void init(M::Body *b,DEM::ParticleCU const *p,DEM::DynParticleCU const *d,size_t n,double C,double dt,double mass,double inertia){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;auto &v=b[i];v.clear();
    double mm[6]={p[i].m/mass,p[i].m/mass,p[i].m/mass,p[i].I.x/inertia,p[i].I.y/inertia,p[i].I.z/inertia};
    double zz[6]={d[i].v.x/C,d[i].v.y/C,d[i].v.z/C,d[i].w.x*dt,d[i].w.y*dt,d[i].w.z*dt};
    for(int j=0;j<6;j++){v.M[j]=mm[j];v.z[j]=v.t[j]=zz[j];}
    v.Q[0]=d[i].Q.w;v.Q[1]=d[i].Q.x;v.Q[2]=d[i].Q.y;v.Q[3]=d[i].Q.z;
}
__global__ void assemble(Claim *a,size_t cap,M::Body *b,double const *f,double const *rho,FLBM::lbm_aux const *la){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=cap||a[i].owner<0)return;auto &c=a[i];c.math.owner=c.owner;
    int cv[15][3];for(int k=0;k<15;k++){cv[k][0]=int(la->C[k].x);cv[k][1]=int(la->C[k].y);cv[k][2]=int(la->C[k].z);}
    double cell_f[15];load_populations(f,c.cell,la,cell_f);
    M::assemble(c.math,b[c.owner],c.lever,c.B,cell_f,rho[c.cell],cv);
    for(int j=0;j<6;j++){atomicAdd(&b[c.owner].b[j],c.math.b[j]);for(int k=j;k<6;k++)atomicAdd(&b[c.owner].D[6*j+k],c.math.D[6*j+k]);}atomicAdd(&b[c.owner].owned,1u);
}
__global__ void solve(M::Body *b,DEM::ParticleCU const *p,size_t n){size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n&&!prescribed(p[i]))M::solve(b[i]);}
__global__ void surfaces(Claim *a,size_t cap,M::Body const *b,double const *f,double const *rho,FLBM::lbm_aux const *la){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=cap||a[i].owner<0)return;auto &c=a[i];int cv[15][3],op[15];
    for(int k=0;k<15;k++){cv[k][0]=int(la->C[k].x);cv[k][1]=int(la->C[k].y);cv[k][2]=int(la->C[k].z);op[k]=la->Op[k];}
    double cell_f[15];load_populations(f,c.cell,la,cell_f);
    M::surface_omega(c.math,b[c.owner],c.lever,c.surface,c.omega,cell_f,rho[c.cell],la->W,cv,op);
    for(int k=0;k<15;k++)c.solid[k]=M::mul(c.math.trace,c.B,c.omega[k]);
}
__global__ void acceleration(real3 *force,double const *rho,size_t nc,real3 acc){size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i<nc)force[i]=rho[i]*acc;}
__global__ void collision(Cell const *cells,Claim const *claims,size_t nc,bool const *mask,double const *f,double *post,
    double const *rho,real3 const *vel,real3 const *force,FLBM::lbm_aux const *la,int *error){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=nc)return;
    if(mask[i]){for(int k=0;k<15;k++)post[FLBM::FIDX(nc,15,0,k,i)]=f[FLBM::FIDX(nc,15,0,la->Op[k],i)];return;}
    auto const &c=cells[i];double tau=la->Tau[0];real3 shifted=vel[i]+(la->dt*tau/rho[i])*force[i];
    for(int k=0;k<15;k++){
        double solid=0,absolute=0;for(int j=c.head;j!=-1;j=claims[j].next){solid=__dadd_rn(solid,claims[j].solid[k]);absolute+=fabs(claims[j].solid[k]);}
        double before=f[FLBM::FIDX(nc,15,0,k,i)],bulk=MechSysRobustEdge::div(MechSysRobustEdge::mul(1-c.B,before-FLBM::FeqFluid(k,rho[i],shifted,la)),tau);
        double after=__dsub_rn(before,__dsub_rn(bulk,solid));double pad=64*DBL_EPSILON*(1+fabs(before)+fabs(bulk)+absolute);
        if(!(before>fabs(bulk)+absolute+pad)||!(after>0&&after<=DBL_MAX))atomicOr(error,16);post[FLBM::FIDX(nc,15,0,k,i)]=after;
    }
}
__global__ void reactions(Claim *a,size_t cap,M::Body const *b,Impulse *imp,FLBM::lbm_aux const *la){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=cap||a[i].owner<0)return;auto &c=a[i];int cv[15][3];
    for(int k=0;k<15;k++){cv[k][0]=int(la->C[k].x);cv[k][1]=int(la->C[k].y);cv[k][2]=int(la->C[k].z);}
    M::reaction(c.math,b[c.owner],c.lever,c.solid,cv);for(int j=0;j<3;j++){atomicAdd(&imp[c.owner].J[j],c.math.j[j]);atomicAdd(&imp[c.owner].Kw[j],c.math.kw[j]);atomicAdd(&imp[c.owner].Kb[j],c.math.kb[j]);}
}
__global__ void moments(bool const *mask,double const *f,double *rho,real3 *vel,real3 *force,FLBM::lbm_aux const *la,int *error,bool compensated=false){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=la->Ncells)return;double r=0;real3 v=make_real3(0,0,0);
    double cell_f[15];load_populations(f,i,la,cell_f);
    for(int k=0;k<15;k++){double z=cell_f[k];if(!(z>0&&z<=DBL_MAX))atomicOr(error,256);r+=z;v=v+z*la->C[k];}
    if(compensated){bool valid=true;r=MechSysDensitySum::sum15(cell_f,valid);if(!valid)atomicOr(error,256);}
    if(!(r>0&&r<=DBL_MAX))atomicOr(error,256);
    rho[i]=mask[i]?0:r;vel[i]=mask[i]?make_real3(0,0,0):(la->Cs/r)*v;force[i]=make_real3(0,0,0);
}
__global__ void apply(Impulse const *a,DEM::ParticleCU *p,DEM::DynParticleCU *d,size_t n,bool vv,double Junit,double Kunit,double dt){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    real3 J=Junit*make_real3(a[i].J[0],a[i].J[1],a[i].J[2]),K=Kunit*make_real3(a[i].Kb[0],a[i].Kb[1],a[i].Kb[2]);
    if(vv&&!prescribed(p[i])){d[i].v=d[i].v+(1./p[i].m)*J;d[i].w=d[i].w+make_real3(K.x/p[i].I.x,K.y/p[i].I.y,K.z/p[i].I.z);}
    else if(!vv){d[i].F=d[i].F+(1./dt)*J;p[i].T=p[i].T+(1./dt)*K;}
    d[i].Flbm=(1./dt)*J;
}
__global__ void status(Claim const *a,size_t cap,M::Body const *b,DEM::DynParticleCU const *d,size_t n,int *error){
    size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i<cap&&a[i].owner>=0&&a[i].math.trace.error)atomicOr(error,32);
    if(i<n){if(b[i].trace.error)atomicOr(error,64);auto const &v=d[i];
        double q=v.Q.w*v.Q.w+v.Q.x*v.Q.x+v.Q.y*v.Q.y+v.Q.z*v.Q.z;
        if(!(isfinite(v.x.x)&&isfinite(v.x.y)&&isfinite(v.x.z)&&isfinite(v.v.x)&&isfinite(v.v.y)&&isfinite(v.v.z)&&
             isfinite(v.w.x)&&isfinite(v.w.y)&&isfinite(v.w.z)&&fabs(q-1)<3e-9))atomicOr(error,512);}
}
__global__ void displacement(real3 const *v,real3 const *old,double *out,size_t n){size_t i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)out[i]=norm(v[i]-old[i]);}
__global__ void clock(DEM::dem_aux *da,FLBM::lbm_aux *la,size_t n,double t){da->iter=n;da->Time=t;la->iter=n;la->Time=t;}
} // Detail

class Session {
    Domain &d_;Options options_;M::Units units_;bool vv_,ready_=false,failed_=false,any_free_=false;size_t step_=0,capacity_=0;int reach_=0;
    thrust::device_vector<Cell> cells_;thrust::device_vector<Claim> claims_;
    thrust::device_vector<M::Body> bodies_;thrust::device_vector<Impulse> impulses_;thrust::device_vector<int> error_;
    std::vector<uint64_t> ids_;
    void checked(){ck(cudaGetLastError());Detail::status<<<(std::max(capacity_,ids_.size())+127)/128,128>>>(ptr(claims_),capacity_,ptr(bodies_),d_.DEMDOM.pDynParticlesCU,ids_.size(),ptr(error_));ck(cudaDeviceSynchronize());int e=read(ptr(error_),1)[0];if(e)throw std::runtime_error("CB: geometry/collision/graph/state failure bits="+std::to_string(e));}
    void reset(){auto &m=d_.DEMDOM;m.pReset<<<(m.demaux.nparts+m.demaux.ncoint+127)/128,128>>>(m.pParticlesCU,m.pDynParticlesCU,m.pInteractons,m.pComInteractons,m.pdemaux,m.pExtraParams);ck(cudaDeviceSynchronize());}
    void forces(){auto &m=d_.DEMDOM;if(m.demaux.nvvint)m.pForceVV<<<(m.demaux.nvvint+127)/128,128>>>(m.pInteractons,m.pComInteractons,m.pDynInteractonsVV,m.pParticlesCU,m.pDynParticlesCU,m.pdemaux,m.pExtraParams);ck(cudaDeviceSynchronize());}
    void rebuild(bool pv){if(!any_free_)return;auto &m=d_.DEMDOM;Detail::displacement<<<(m.demaux.nverts+127)/128,128>>>(m.pVertsCU,m.pVertsoCU,m.pMaxDCU,m.demaux.nverts);ck(cudaDeviceSynchronize());auto v=read(m.pMaxDCU,m.demaux.nverts);if(*std::max_element(v.begin(),v.end())>=m.Alpha){if(pv)m.UpdateContactsDevice();else m.UpdateContactsDevicePreForce();}}
public:
    // Observational hook for benchmark audit packets. Must not mutate domain.
    std::function<void(char const *)> hook;
    Session(Domain &d,Options const &o):d_(d),options_(o),units_(d.dx,d.dt,d.Fconv),vv_(o.selection==Selection::ImplicitVV){
        need((o.selection==Selection::ExplicitPV||o.selection==Selection::ImplicitVV)&&o.max_intervals>0,"CB: explicit selection and finite interval limit required");
        need(d.DEMDOM.UseVelocityVerlet==vv_,"CB: select and verify matching mechanical integrator");
    }
    Session(Session const&)=delete;Session &operator=(Session const&)=delete;
    static void validate(Domain const &d,Options const &o){
        auto const &m=d.DEMDOM;auto const &l=d.LBMDOM;M::Units u(d.dx,d.dt,d.Fconv);
        need((o.selection==Selection::ExplicitPV||o.selection==Selection::ImplicitVV)&&o.max_intervals>0,"CB: default-off selection/window");
        need(o.density_sum==DensitySum::SequentialV1||o.density_sum==DensitySum::CompensatedV1,"CB: unknown density summation version");
        need(!m.Initialized&&!m.CheckpointReady&&!m.VVStateReady&&m.iter==0&&d.iter==0&&d.Time==0,"CB: fresh domain; no coupled restart");
        need(sizeof(real)==8&&l.Nl==1&&l.Nneigh==15&&l.Sc==0&&l.G[0]==0&&l.Gmix==0&&l.Tau[0]>.5,"CB: FP64 smooth single-fluid D3Q15 only");
        need(l.Ncells>0&&l.Ncells<=o.max_cells&&m.Particles.Size()>0&&m.Particles.Size()<=512,"CB: cell/body resource bound");
        need(m.ContactLaw==0&&m.SphereCoulombMode==0&&!m.SphereTensileCutoff&&!m.SphereFirstContactCorrection&&!m.MostlySpheres&&!m.BInteractons.Size(),"CB: original sphere contact-law subset");
        for(int j=0;j<3;j++)need(std::isfinite(o.fluid_acceleration(j)),"CB: finite imposed acceleration");
        for(auto *p:m.Particles){
            bool none=!p->vxf&&!p->vyf&&!p->vzf&&!p->wxf&&!p->wyf&&!p->wzf;
            bool all=p->vxf&&p->vyf&&p->vzf&&p->wxf&&p->wyf&&p->wzf;
            need((none||all)&&!p->FixFree,"CB: free or fully prescribed only; no partial constraints");
            need(p->Verts.Size()==1&&!p->Edges.Size()&&!p->Faces.Size()&&!p->Bdry&&std::isfinite(p->Props.R)&&p->Props.R>0&&std::isfinite(p->Props.m)&&p->Props.m>0&&p->Props.Gv==0&&p->Props.Gm==0,"CB: positive physical canonical spheres, no artificial damping");
            need(!all||norm(p->w)==0,"CB: prescribed spin not supported");
            need(norm(p->Flbmf)==0,"CB: prescribed Flbm diagnostic offset unsupported");
            for(int j=0;j<3;j++)need(std::isfinite(p->x(j))&&std::isfinite(p->v(j))&&std::isfinite(p->w(j))&&std::isfinite(p->Ff(j))&&std::isfinite(p->Tf(j))&&std::isfinite(p->I(j))&&p->I(j)>0,"CB: finite state/force and positive inertia");
        }
    }
    void initialize(){
        validate(d_,options_);auto &m=d_.DEMDOM;auto &l=d_.LBMDOM;double maxr=0;bool any_free=false;
        m.Alpha=d_.Alpha;need(m.Alpha>0,"CB: positive Verlet margin");
        m.Xmin=m.Ymin=m.Zmin=0;m.Xmax=d_.PeriodicX?l.Ndim(0)*d_.dx:0;m.Ymax=d_.PeriodicY?l.Ndim(1)*d_.dx:0;m.Zmax=d_.PeriodicZ?l.Ndim(2)*d_.dx:0;
        for(auto *p:m.Particles){maxr=std::max(maxr,p->Props.R);any_free=any_free||p->IsFree();}
        need(maxr/d_.dx<double(INT_MAX/4),"CB: radius/cell-index bound");
        reach_=int(ceil(maxr/d_.dx+.5))+1;size_t side=2*reach_+1;
        need(side<=l.Ndim(0)&&side<=l.Ndim(1)&&side<=l.Ndim(2),"CB: periodic alias/capacity guard");
        capacity_=m.Particles.Size()*side*side*side;need(capacity_<size_t(INT_MAX)&&capacity_<=options_.max_claim_bytes/sizeof(Claim),"CB: claim capacity byte bound");
        any_free_=any_free;
        if(any_free)m.Solve(0,d_.dt,d_.dt,NULL,NULL,NULL,false,1);
        else {
            // No free bodies means no physical contact-force pairs. Avoid the
            // standalone solver's free-body skin/grid heuristics (MaxDmax=0).
            m.Dt=d_.dt;m.Nproc=1;m.Initialize(d_.dt);m.AssignPersistentIds();m.Per=Vec3_t(m.Xmax,m.Ymax,m.Zmax);
            m.A.Resize(m.Particles.Size());m.Wdot.Resize(m.Particles.Size());m.FreePar.Resize(0);m.NoFreePar.Resize(0);
            for(size_t i=0;i<m.Particles.Size();i++){m.Particles[i]->Index=i;m.NoFreePar.Push(i);m.A[i]=m.Wdot[i]=OrthoSys::O;m.Particles[i]->F=m.Particles[i]->Ff;m.Particles[i]->T=m.Particles[i]->Tf;}
            m.MTD=new DEM::MtData[1];m.MTD[0].N_Proc=1;m.MTD[0].ProcRank=0;m.MTD[0].Dom=&m;m.MTD[0].Dmx=0;m.UpLoadDevice(1,true,true);
        }
        m.VVStateReady=false;l.UpLoadDevice(1);std::set<uint64_t> unique;
        for(auto *p:m.Particles){need(p->PersistentId&&unique.insert(p->PersistentId).second,"CB: unique persistent IDs");ids_.push_back(p->PersistentId);}
        cells_.resize(l.Ncells);claims_.resize(capacity_);bodies_.resize(ids_.size());impulses_.resize(ids_.size());error_.resize(1);ready_=true;
    }
    void advance_one(){
        need(ready_&&!failed_&&step_<options_.max_intervals,"CB: terminal state/interval limit");auto &m=d_.DEMDOM;auto &l=d_.LBMDOM;
        try {
            need(d_.dx==units_.dx&&d_.dt==units_.dt&&d_.Fconv==units_.conversion&&m.UseVelocityVerlet==vv_&&m.Particles.Size()==ids_.size(),"CB: frozen scales/mode/body extent");
            for(size_t i=0;i<ids_.size();i++)need(ids_[i]==m.Particles[i]->PersistentId,"CB: identity changed");
            if(hook)hook("before");if(!vv_){reset();if(hook)hook("pre_force");forces();if(hook)hook("force");}
            auto history=read(m.pDynInteractonsVV,m.demaux.nvvint);
            cells_.assign(l.Ncells,Cell{});error_.assign(1,0);impulses_.assign(ids_.size(),Impulse{});
            Detail::claims<<<(capacity_+127)/128,128>>>(ptr(claims_),ptr(cells_),capacity_,reach_,m.pParticlesCU,m.pDynParticlesCU,m.pdemaux,l.plbmaux,l.pIsSolid,ptr(error_));
            Detail::partition<<<(l.Ncells+127)/128,128>>>(ptr(claims_),ptr(cells_),l.Ncells,l.Tau[0],ptr(error_));
            Detail::init<<<(ids_.size()+127)/128,128>>>(ptr(bodies_),m.pParticlesCU,m.pDynParticlesCU,ids_.size(),units_.C,units_.dt,units_.mass,units_.inertia);
            Detail::assemble<<<(capacity_+127)/128,128>>>(ptr(claims_),capacity_,ptr(bodies_),l.pF,l.pRho,l.plbmaux);
            if(vv_)Detail::solve<<<(ids_.size()+127)/128,128>>>(ptr(bodies_),m.pParticlesCU,ids_.size());
            Detail::surfaces<<<(capacity_+127)/128,128>>>(ptr(claims_),capacity_,ptr(bodies_),l.pF,l.pRho,l.plbmaux);
            real3 acc=make_real3(options_.fluid_acceleration(0),options_.fluid_acceleration(1),options_.fluid_acceleration(2));
            Detail::acceleration<<<(l.Ncells+127)/128,128>>>(l.pBForce,l.pRho,l.Ncells,acc);
            Detail::collision<<<(l.Ncells+127)/128,128>>>(ptr(cells_),ptr(claims_),l.Ncells,l.pIsSolid,l.pF,l.pFtemp,l.pRho,l.pVel,l.pBForce,l.plbmaux,ptr(error_));
            Detail::reactions<<<(capacity_+127)/128,128>>>(ptr(claims_),capacity_,ptr(bodies_),ptr(impulses_),l.plbmaux);checked();if(hook)hook("exchange");
            std::swap(l.pF,l.pFtemp);FLBM::cudaStream1<<<(l.Ncells+127)/128,128>>>(l.pF,l.pFtemp,l.pBForce,l.plbmaux);std::swap(l.pF,l.pFtemp);
            Detail::moments<<<(l.Ncells+127)/128,128>>>(l.pIsSolid,l.pF,l.pRho,l.pVel,l.pBForce,l.plbmaux,ptr(error_),options_.density_sum==DensitySum::CompensatedV1);
            if(vv_)reset();Detail::apply<<<(ids_.size()+127)/128,128>>>(ptr(impulses_),m.pParticlesCU,m.pDynParticlesCU,ids_.size(),vv_,units_.impulse,units_.angular_impulse,d_.dt);checked();
            auto after=read(m.pDynInteractonsVV,m.demaux.nvvint);need(after.size()==history.size()&&(after.empty()||!memcmp(after.data(),history.data(),after.size()*sizeof(after[0]))),"CB: hydro changed contact history");
            if(hook)hook("hydro");
            if(vv_){m.pVerletStep1<<<(ids_.size()+127)/128,128>>>(m.pVertsCU,m.pParticlesCU,m.pDynParticlesCU,m.pA,m.pdemaux);m.pOrientationUpdate<<<(ids_.size()+127)/128,128>>>(m.pVertsCU,m.pParticlesCU,m.pDynParticlesCU,m.pWdot,m.pdemaux);rebuild(false);if(hook)hook("pre_force");forces();if(hook)hook("force");m.pFinalizeVelocity<<<(ids_.size()+127)/128,128>>>(m.pParticlesCU,m.pDynParticlesCU,m.pA,m.pdemaux);m.pFinalizeRotation<<<(ids_.size()+127)/128,128>>>(m.pParticlesCU,m.pDynParticlesCU,m.pWdot,m.pdemaux);}
            else {m.pTranslate<<<(ids_.size()+127)/128,128>>>(m.pVertsCU,m.pParticlesCU,m.pDynParticlesCU,m.pdemaux,m.pExtraParams);m.pRotate<<<(ids_.size()+127)/128,128>>>(m.pVertsCU,m.pParticlesCU,m.pDynParticlesCU,m.pdemaux,m.pExtraParams);rebuild(true);}
            checked();step_++;d_.iter=m.iter=step_;d_.Time=m.Time=l.Time=step_*d_.dt;m.demaux.iter=step_;m.demaux.Time=d_.Time;
            Detail::clock<<<1,1>>>(m.pdemaux,l.plbmaux,step_,d_.Time);ck(cudaDeviceSynchronize());m.VVStateReady=false;if(hook)hook("after");
        }catch(...){failed_=true;ready_=false;m.VVStateReady=false;throw;}
    }
    void download(){
        need(ready_&&!failed_,"CB: no snapshot after terminal failure");auto &m=d_.DEMDOM;auto &l=d_.LBMDOM;m.DnLoadDevice(1,false);l.DnLoadDevice(1);
        auto p=read(m.pParticlesCU,ids_.size());auto impulses=this->impulses();auto cells=this->cells();
        for(size_t i=0;i<ids_.size();i++){auto *v=m.Particles[i];v->T=Vec3_t(p[i].T.x,p[i].T.y,p[i].T.z);if(vv_){v->F+=v->Flbm;v->T+=units_.torque*Vec3_t(impulses[i].Kb[0],impulses[i].Kb[1],impulses[i].Kb[2]);}}
        for(size_t z=0;z<l.Ndim(2);z++)for(size_t y=0;y<l.Ndim(1);y++)for(size_t x=0;x<l.Ndim(0);x++)l.Gamma[x][y][z]=cells[x+l.Ndim(0)*(y+l.Ndim(1)*z)].Gamma;
    }
    size_t completed()const{return step_;}size_t capacity()const{return capacity_;}
    M::Units const &units()const{return units_;}
    std::vector<Claim> claims(){return read(ptr(claims_),capacity_);}std::vector<Cell> cells(){return read(ptr(cells_),cells_.size());}
    std::vector<M::Body> bodies(){return read(ptr(bodies_),ids_.size());}std::vector<Impulse> impulses(){return read(ptr(impulses_),ids_.size());}
};
}}} // LBMDEM::Experimental::Benchmark
#endif // default-off explicit feature
#endif
