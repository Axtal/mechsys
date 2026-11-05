// GPL-3.0-or-later. Analytic contact-law and exported-row regression.
#include <mechsys/dem/domain.h>
#include <cmath>
#include <stdexcept>
#include <string>

static void check(bool value, char const *message)
{
    if (!value) throw std::runtime_error(message);
}
static void close(double a, double b, char const *message)
{
    check(std::isfinite(a) && std::abs(a-b)<=1e-11*(1+std::abs(b)),message);
}
#ifdef USE_CUDA
static size_t cudaComparisons=0;
template<class T> static T *ptr(thrust::device_vector<T> &v)
{ return thrust::raw_pointer_cast(v.data()); }

static void device_case(DEM::Domain &d, DEM::CInteractonSphere const &c,
                        int law, int mode, bool cutoff, bool onset, double dt,
                        Vec3_t const &initialHistory=OrthoSys::O, bool persistent=false)
{
    check(sizeof(real)==8,"FP64 required");
    thrust::host_vector<DEM::ParticleCU> p(2);
    thrust::host_vector<DEM::DynParticleCU> q(2);
    for (int i=0;i<2;i++) {
        d.Particles[i]->UpdateRotCacheIfChanged();
        DEM::UploadParticle(q[i],p[i],*d.Particles[i]);
    }
    DEM::InteractonCU interaction={};
    interaction.BothFree=true; interaction.Kn=c.Kn; interaction.Kt=c.Kt;
    interaction.Gn=c.Gn; interaction.Gt=c.Gt; interaction.Mu=c.Mu;
    DEM::ComInteractonCU common={}; common.I1=0; common.I2=1;
    DEM::DynInteractonCU history={}; history.Idx=0;
    double const delta=2-norm(d.Particles[1]->x-d.Particles[0]->x);
    history.PrevDelta=persistent?delta:-.1; history.InContact=persistent;
    double const stiffness=c.Kt*(law?std::sqrt(delta):1);
    history.Ft=make_real3(stiffness*initialHistory(0),stiffness*initialHistory(1),
                         stiffness*initialHistory(2));
    DEM::dem_aux aux={}; aux.nparts=2;aux.ncoint=1;aux.nvvint=1;
    aux.dt=dt;aux.sphereCoulombMode=mode;aux.sphereTensileCutoff=cutoff;
    aux.sphereFirstContactCorrection=onset;
    thrust::device_vector<DEM::ParticleCU> dp=p;
    thrust::device_vector<DEM::DynParticleCU> dq=q;
    thrust::device_vector<DEM::InteractonCU> di(1,interaction);
    thrust::device_vector<DEM::ComInteractonCU> dc(1,common);
    thrust::device_vector<DEM::DynInteractonCU> dh(1,history);
    thrust::device_vector<DEM::dem_aux> da(1,aux);
    DEM::Reset<<<1,32>>>(ptr(dp),ptr(dq),ptr(di),ptr(dc),ptr(da),nullptr);
    if(law) DEM::CalcForceVV_Hertz<<<1,32>>>(ptr(di),ptr(dc),ptr(dh),ptr(dp),ptr(dq),ptr(da),nullptr);
    else DEM::CalcForceVV<<<1,32>>>(ptr(di),ptr(dc),ptr(dh),ptr(dp),ptr(dq),ptr(da),nullptr);
    check(cudaDeviceSynchronize()==cudaSuccess,"contact kernel failed");
    common=dc[0]; auto body=DEM::DynParticleCU(dq[0]);
    close(common.Fnnet.x,c.Fnet(0),"CPU/GPU normal elastic");
    close(common.Ftnet.y,c.Ftnet(1),"CPU/GPU tangential elastic y");
    close(common.Ftnet.z,c.Ftnet(2),"CPU/GPU tangential elastic z");
    close(common.Fndpot.x,c.Fndpot(0),"CPU/GPU normal dashpot");
    close(common.Ftdpot.y,c.Ftdpot(1),"CPU/GPU tangential dashpot y");
    close(common.Ftdpot.z,c.Ftdpot(2),"CPU/GPU tangential dashpot z");
    close(body.F.x,c.F1(0),"CPU/GPU force x");
    close(body.F.y,c.F1(1),"CPU/GPU force y");
    close(body.F.z,c.F1(2),"CPU/GPU force z");
    check(common.Nc==1,"GPU geometric contact count");
    close(common.NormalOverlap,.1,"GPU normal overlap");
    ++cudaComparisons;
}
#endif

static void contact_case(int law,int mode,bool cutoff,bool onset)
{
    DEM::Domain d;
    d.AddSphere(-1,Vec3_t(0,0,0),1,1);
    d.AddSphere(-2,Vec3_t(1.9,0,0),1,1);
    d.Particles[0]->v=Vec3_t(cutoff?-20:1,4,0);
    DEM::CInteractonSphere c(d.Particles[0],d.Particles[1]);
    c.Kn=100;c.Kt=20;c.Gn=2;c.Gt=onset?0:3;c.Mu=onset?100:.5;
    c.beta=0;c.eta=0;c.PrevDelta=-.1;c.InContact=false;
    double dt=onset?.02:.2;
    c.CalcForce(dt,OrthoSys::O,1,law,mode,cutoff,onset);
    double scale=law?std::sqrt(.1):1;
    double damping=law?std::sqrt(scale):1;
    double elastic=100*.1*scale;
    double dash=2*damping*(cutoff?-20:1);
    if(cutoff && elastic+dash<0) dash=-elastic;
    double trial=20*scale*dt*4*(onset?.5:1);
    double tDash=(onset?0:3)*damping*4;
    // Mode 1 triggers on total tangential force; modes 0, 2 and 3 on
    // elastic force. Modes 1/3 project the total onto the Coulomb ball;
    // modes 0/2 cap elastic force, and mode 2 drops the sliding dashpot.
    double const normalLimit=((mode==1||mode==3)?std::abs(elastic+dash):elastic);
    double const limit=(onset?100:.5)*normalLimit;
    bool const triggered = mode==1 ? (trial+tDash>limit) : (trial>limit);
    double tangential=trial;
    if(mode==1||mode==3) {
        if(triggered) tangential=std::min(trial+tDash,limit)-tDash;
    }
    else if(triggered) {tangential=limit;if(mode==2)tDash=0;}
    close(c.Fnet(0),elastic,"analytic normal elastic");
    close(c.Fndpot(0),dash,"analytic normal dashpot");
    close(c.Ftnet(1),tangential,"analytic tangential elastic");
    close(c.Ftdpot(1),tDash,"analytic tangential dashpot");
    close(c.F1(0),-(elastic+dash),"analytic total normal");
    close(c.F1(1),-(tangential+tDash),"analytic total tangential");
    close(norm(c.F1+c.F2),0,"action/reaction");
    check(c.Nc==1,"CPU geometric contact count");
    close(c.NormalOverlap,.1,"CPU normal overlap");
#ifdef USE_CUDA
    // The CUDA kernels implement only modes 0 (legacy) and 1 (LAMMPS); modes 2
    // (PFC) and 3 (Yade) are CPU only and are clamped to 0 on upload.
    if (mode<=1) device_case(d,c,law,mode,cutoff,onset,dt);
#endif
}

// Independent prescribed-force cases: effective coefficients are Kn=100,
// Gn=2, Kt=20, Gt=3 for both laws, giving normal limits 5 (elastic) and 6
// (total). Expected totals are explicit, including opposed and 2-D forces.
struct TangentialCase {
    double elasticY,elasticZ,dashY,dashZ,dt;
    double expectedY[4],expectedZ[4];
    bool sliding[4];
};
static TangentialCase const tangentialCases[] = {
    {4,0,12,0,.01, {16,6,16,16}, {0,0,0,0}, {false,true,false,false}},
    {8,0,-5,0,.01, {0,3,5,3}, {0,0,0,0}, {true,false,true,true}},
    {8,0,-8,0,0, {-3,0,5,0}, {0,0,0,0}, {true,false,true,true}},
    {8,0,-8+1e-9,0,.01, {-3+1e-9,1e-9,5,1e-9}, {0,0,0,0}, {true,false,true,true}},
    {8,0,-8-1e-9,0,.01, {-3-1e-9,-1e-9,5,-1e-9}, {0,0,0,0}, {true,false,true,true}},
    {8,0,0,6,.01, {5,4.8,5,4.8}, {6,3.6,0,3.6}, {true,true,true,true}},
    {0,0,12,0,.01, {12,6,12,12}, {0,0,0,0}, {false,true,false,false}},
    {8,0,-20,0,.01, {-15,-6,5,-6}, {0,0,0,0}, {true,true,true,true}}
};

static void tangential_case(int law,int mode,TangentialCase const &test)
{
    DEM::Domain d;
    d.AddSphere(-1,Vec3_t(0,0,0),1,1);
    d.AddSphere(-2,Vec3_t(1.9,0,0),1,1);
    Vec3_t const vt(0,test.dashY/3,test.dashZ/3);
    d.Particles[0]->v=Vec3_t(1,vt(1),vt(2));
    double const delta=2-norm(d.Particles[1]->x-d.Particles[0]->x);
    double const scale=law?std::sqrt(delta):1;
    double const damping=std::sqrt(scale);
    DEM::CInteractonSphere c(d.Particles[0],d.Particles[1]);
    c.Kn=100/scale;c.Gn=2/damping;c.Kt=20/scale;c.Gt=3/damping;
    c.Mu=.5;c.beta=0;c.eta=0;c.PrevDelta=delta;c.InContact=true;
    Vec3_t const initialHistory=Vec3_t(0,test.elasticY/20,test.elasticZ/20)-vt*test.dt;
    c.Fdvv=initialHistory;
    c.CalcForce(test.dt,OrthoSys::O,1,law,mode,false,false);
    Vec3_t const expected(0,test.expectedY[mode],test.expectedZ[mode]);
    Vec3_t const dash=(mode==2 && test.sliding[mode])?OrthoSys::O:Vec3_t(0,test.dashY,test.dashZ);
    for(int axis=1;axis<3;axis++) {
        close(c.Ftnet(axis)+c.Ftdpot(axis),expected(axis),"prescribed tangential total");
        close(c.Ftdpot(axis),dash(axis),"prescribed tangential dashpot");
        close(c.Fdvv(axis),(expected(axis)-dash(axis))/20,"prescribed elastic history");
        close(c.F1(axis),-expected(axis),"prescribed applied force");
    }
    check(c.Nsc==size_t(test.sliding[mode]),"elastic versus total sliding trigger");
    close(norm(c.F1+c.F2),0,"prescribed action/reaction");
#ifdef USE_CUDA
    if(mode<=1) device_case(d,c,law,mode,false,false,test.dt,initialHistory,true);
#endif
}

static void fast_tangential_case(int mode,TangentialCase const &test)
{
    DEM::Domain d;
    d.Nproc=1;d.SphereCoulombMode=mode;
    d.AddSphere(-1,Vec3_t(0,0,0),1,1);
    d.AddSphere(-2,Vec3_t(1.9,0,0),1,1);
    for(size_t i=0;i<2;i++) {
        auto *p=d.Particles[i];
        p->Props.m=2;p->Props.Kn=200;p->Props.Kt=20;
        p->Props.Gn=2;p->Props.Gt=3;p->Props.Mu=.5;
        p->Props.Beta=0;p->Props.Eta=0;p->F=p->T=OrthoSys::O;
    }
    Vec3_t const vt(0,test.dashY/3,test.dashZ/3);
    d.Particles[0]->v=Vec3_t(1,vt(1),vt(2));
    d.ListPosPairs.Push(std::make_pair(size_t(0),size_t(1)));
    size_t const key=HashFunction(0,1);
    d.FricSpheres[key]=Vec3_t(0,test.elasticY/20,test.elasticZ/20)-vt*test.dt;
    d.CalcForceSphere(test.dt);
    Vec3_t const expected(0,test.expectedY[mode],test.expectedZ[mode]);
    Vec3_t const dash=(mode==2 && test.sliding[mode])?OrthoSys::O:Vec3_t(0,test.dashY,test.dashZ);
    for(int axis=1;axis<3;axis++) {
        close(d.Particles[0]->F(axis),-expected(axis),"fast prescribed applied force");
        close(d.FricSpheres[key](axis),(expected(axis)-dash(axis))/20,"fast elastic history");
    }
    close(norm(d.Particles[0]->F+d.Particles[1]->F),0,"fast action/reaction");
}

static void output_cases(std::string const &prefix)
{
    DEM::Domain d;
    d.AddSphere(-1,Vec3_t(0,0,0),1,1);
    d.AddSphere(-2,Vec3_t(1.9,0,0),1,1);
    DEM::CInteractonSphere c(d.Particles[0],d.Particles[1]);
    d.PairtoCInt[HashFunction(0,1)]=&c;
    c.Fnet=c.Ftnet=c.Fndpot=c.Ftdpot=c.Fther=OrthoSys::O;
    // Rolling must use its own force history, even if the normal scratch field
    // has never been assigned by a GPU force evaluation.
    c.Fn=Vec3_t(NAN,NAN,NAN);
    c.Fdr=Vec3_t(3,4,5);
    c.Nc=1;c.NormalOverlap=.1;
    d.WriteBF((prefix+"-overlap-zero").c_str());
    c.Nc=0;c.NormalOverlap=0;c.Fnet=Vec3_t(1,0,0);
    d.WriteBF((prefix+"-force-only").c_str());
    c.Fnet=OrthoSys::O;c.Ftnet=Vec3_t(0,1e-15,0);
    d.WriteBF((prefix+"-tangential-only").c_str());
    c.Ftnet=OrthoSys::O;c.Ftdpot=Vec3_t(0,2,0);
    d.WriteBF((prefix+"-dashpot-only").c_str());
    c.Ftdpot=OrthoSys::O;
    check(!DEM::HasContactOutput(c),"empty candidate must not be exported");
    d.WriteBF((prefix+"-empty").c_str());
    d.PairtoCInt.clear();
}

int main(int argc,char **argv) try
{
    check(argc==2,"usage: test_contact OUTPUT_PREFIX");
    // Modes 0-3: legacy elastic cap with permanent dashpot, LAMMPS total-force
    // trigger and cap, PFC elastic cap with the dashpot off while sliding, and
    // the Yade hybrid (total normal limit, elastic trigger, total cap).
    size_t cpuAnalyticCases=0,cpuFastCases=0;
    for(int law=0;law<2;law++)for(int mode=0;mode<4;mode++) {
        contact_case(law,mode,false,false);++cpuAnalyticCases;
        contact_case(law,mode,true,false);++cpuAnalyticCases;
        contact_case(law,mode,false,true);++cpuAnalyticCases;
    }
    for(auto const &test:tangentialCases) {
        for(int mode=0;mode<4;mode++) {
            for(int law=0;law<2;law++) {
                tangential_case(law,mode,test);++cpuAnalyticCases;
            }
            fast_tangential_case(mode,test);++cpuFastCases;
        }
    }
    output_cases(argv[1]);
    std::cout<<"CONTACT_REGRESSION_PASS cpu_analytic_cases="<<cpuAnalyticCases
             <<" cpu_fast_cases="<<cpuFastCases
#ifdef USE_CUDA
             <<" cuda_comparisons="<<cudaComparisons
#endif
             <<std::endl;
    return 0;
} MECHSYS_CATCH
