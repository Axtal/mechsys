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
template<class T> static T *ptr(thrust::device_vector<T> &v)
{ return thrust::raw_pointer_cast(v.data()); }

static void device_case(DEM::Domain &d, DEM::CInteractonSphere const &c,
                        int law, int mode, bool cutoff, bool onset, double dt)
{
    check(sizeof(real)==8,"FP64 required");
    thrust::host_vector<DEM::ParticleCU> p(2);
    thrust::host_vector<DEM::DynParticleCU> q(2);
    for (int i=0;i<2;i++) {
        d.Particles[i]->UpdateRotCacheIfChanged();
        DEM::UploadParticle(q[i],p[i],*d.Particles[i]);
    }
    DEM::InteractonCU interaction={};
    interaction.BothFree=true; interaction.Kn=100; interaction.Kt=20;
    interaction.Gn=2; interaction.Gt=onset?0:3; interaction.Mu=onset?100:.5;
    DEM::ComInteractonCU common={}; common.I1=0; common.I2=1;
    DEM::DynInteractonCU history={}; history.Idx=0; history.PrevDelta=-.1;
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
    close(common.Ftnet.y,c.Ftnet(1),"CPU/GPU tangential elastic");
    close(common.Fndpot.x,c.Fndpot(0),"CPU/GPU normal dashpot");
    close(common.Ftdpot.y,c.Ftdpot(1),"CPU/GPU tangential dashpot");
    close(body.F.x,c.F1(0),"CPU/GPU force x");
    close(body.F.y,c.F1(1),"CPU/GPU force y");
    check(common.Nc==1,"GPU geometric contact count");
    close(common.NormalOverlap,.1,"GPU normal overlap");
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
    // Modes 0 (legacy) and 2 (PFC) limit on the elastic normal force; modes 1
    // (LAMMPS) and 3 (Yade) on the total normal force. Modes 1 and 3 trigger on
    // and cap the total tangential force, so they apply the same sliding force
    // and keep the dashpot; modes 0 and 2 cap the elastic part, and mode 2
    // alone drops the dashpot once sliding.
    double const normalLimit=((mode==1||mode==3)?std::abs(elastic+dash):elastic);
    double const limit=(onset?100:.5)*normalLimit;
    bool const triggered = (mode==1||mode==3) ? (trial+tDash>limit) : (trial>limit);
    double tangential=trial;
    if(mode==1||mode==3) { if(triggered) tangential=limit-tDash; }
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
    for(int law=0;law<2;law++)for(int mode=0;mode<4;mode++) {
        contact_case(law,mode,false,false);
        contact_case(law,mode,true,false);
        contact_case(law,mode,false,true);
    }
    output_cases(argv[1]);
    std::cout<<"CONTACT_REGRESSION_PASS cases="<<(2*4*3)<<std::endl;
    return 0;
} MECHSYS_CATCH
