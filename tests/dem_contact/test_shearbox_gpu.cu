// Exercise the actual custom kernels and physical-step helper used by the example.
#define main shearbox_example_main
#define LE_REBUILD_VERIFY
#include "../../tdem/shearbox_LE_const_vol.cu"
#undef main
#include <stdexcept>

static void require(bool ok,char const *message) { if(!ok) throw std::runtime_error(message); }
static void close(double a,double b,char const *message) {
    require(std::isfinite(a) && fabs(a-b)<1e-11*(1+fabs(b)),message);
}
template<class T> static T *ptr(thrust::device_vector<T> &v) { return thrust::raw_pointer_cast(v.data()); }

static void contact_outputs(int law,bool cutoff,bool separated)
{
    DEM::Domain d;
    d.AddSphere(-1,Vec3_t(0,0,0),1,1);
    d.AddSphere(-2,Vec3_t(separated?2.1:1.9,0,0),1,1);
    d.Particles[0]->v=Vec3_t(cutoff?-20:1,cutoff?0:4,0);
    thrust::host_vector<DEM::ParticleCU> hp(2);
    thrust::host_vector<DEM::DynParticleCU> hq(2);
    for(int j=0;j<2;++j) { d.Particles[j]->UpdateRotCacheIfChanged();DEM::UploadParticle(hq[j],hp[j],*d.Particles[j]); }
    DEM::InteractonCU i={};i.BothFree=true;i.Kn=100;i.Kt=20;i.Gn=2;i.Gt=3;i.Mu=100;
    DEM::ComInteractonCU c={};c.I1=0;c.I2=1;c.Nc=19;c.NormalOverlap=.9;
    DEM::DynInteractonCU h={};h.Ft=make_real3(0,.2,0);h.Fr=make_real3(0,.1,0);
    DEM::dem_aux a={};a.nparts=2;a.ncoint=1;a.nvvint=1;a.dt=.02;a.Per=make_real3(10,10,10);
    LeExtra extra={};
    thrust::device_vector<DEM::ParticleCU> dp=hp;
    thrust::device_vector<DEM::DynParticleCU> dq=hq;
    thrust::device_vector<DEM::InteractonCU> di(1,i);
    thrust::device_vector<DEM::ComInteractonCU> dc(1,c);
    thrust::device_vector<DEM::DynInteractonCU> dh(1,h);
    thrust::device_vector<DEM::dem_aux> da(1,a);
    thrust::device_vector<LeExtra> de(1,extra);
    DEM::ResetLE<<<1,32>>>(ptr(dp),ptr(dq),ptr(di),ptr(dc),ptr(da),ptr(de));
    if(law) DEM::CalcForceVV_Hertz_LE<<<1,32>>>(ptr(di),ptr(dc),ptr(dh),ptr(dp),ptr(dq),ptr(da),ptr(de));
    else DEM::CalcForceVV_LE<<<1,32>>>(ptr(di),ptr(dc),ptr(dh),ptr(dp),ptr(dq),ptr(da),ptr(de));
    require(cudaDeviceSynchronize()==cudaSuccess,"LE contact kernel failed");
    c=dc[0];h=dh[0];
    close(c.NormalOverlap,separated?0:.1,"wrong LE overlap");
    require(c.Nc==size_t(!separated),"stale LE contact count");
    real3 sum=c.Fnnet+c.Ftnet+c.Fndpot+c.Ftdpot;
    close(sum.x,h.F.x,"LE normal force decomposition");
    close(sum.y,h.F.y,"LE tangential force decomposition");
    auto q0=DEM::DynParticleCU(dq[0]);auto q1=DEM::DynParticleCU(dq[1]);
    close(q0.F.x+q1.F.x,0,"LE action reaction");
    require(h.InContact==!separated,"LE activation state");
    close(h.PrevDelta,separated?-.1:.1,"LE signed overlap history");
    if(separated) {
        close(h.Ft.y,.2,"temporary gap changed custom tangential history");
        close(h.Fr.y,.1,"temporary gap changed custom rolling history");
    }
    else if(cutoff) close(h.F.x,0,"LE tensile clamp");
    else close(h.F.x,law?100*pow(.1,1.5)+2*pow(.1,.25):12,"LE analytic normal force");
}

static void prepare(DEM::Domain &d,UserData &data,double dt)
{
    d.Alpha=.01;d.Nthread=32;d.UseVelocityVerlet=true;
    for(size_t j=0;j<d.Particles.Size();++j) {
        auto *p=d.Particles[j];p->Props.Kn=1000;p->Props.Kt=200;
        p->Props.Gn=0;p->Props.Gt=0;p->Props.Mu=100;
        p->Props.Beta=0;p->Props.Eta=0;
    }
    data.Lx=data.Ly=10;data.base_alpha=d.Alpha;
    SetPeriodicBox(d,10,10,-5,5);
    PrepareShearSolve(d,dt,"le_test",1);
    SetLeState(0,data.shear_rate,0,0);
    d.pForceVV=d.ContactLaw?DEM::CalcForceVV_Hertz_LE:DEM::CalcForceVV_LE;
    d.pReset=DEM::ResetLE;
    LeExtra extra={};extra.state=g_le_state;
    require(cudaMalloc(&d.pExtraParams,sizeof(extra))==cudaSuccess,"LE state allocation");
    cudaMemcpy(d.pExtraParams,&extra,sizeof(extra),cudaMemcpyHostToDevice);
    ResetContactsLeesEdwardsSpheres(d,0);
    d.UpLoadDevice(1,true,true);
}

static void history_rebuild(int law)
{
    UserData data={};DEM::Domain d(&data,law);
    d.AddSphere(-1,Vec3_t(0,0,0),1,1);
    d.AddSphere(-2,Vec3_t(1.9,0,0),1,1);
    d.Particles[0]->v=Vec3_t(0,1,0);
    prepare(d,data,.001);
    DEM::DynInteractonCU h=d.bDynInteractonsVV[0];
    h.Ft=make_real3(0,.7,0);h.Fr=make_real3(0,0,.3);h.PrevDelta=.1;h.InContact=true;
    d.bDynInteractonsVV[0]=h;
    RefreshDeviceAccelerations(d);
    auto after=DEM::DynInteractonCU(d.bDynInteractonsVV[0]);
    close(after.Ft.y,h.Ft.y,"bootstrap changed tangential history");
    close(after.Fr.z,h.Fr.z,"bootstrap changed rolling history");
    close(after.PrevDelta,h.PrevDelta,"bootstrap changed signed overlap history");
    close(d.demaux.dt,.001,"bootstrap did not restore timestep");
    d.bA[0]=make_real3(3,4,5);d.bWdot[0]=make_real3(6,7,8);
    for(int repeat=0;repeat<2;++repeat) {
        RebuildLeesEdwardsContacts(d,data);
        after=d.bDynInteractonsVV[0];
        close(after.Ft.y,.7,"rebuild changed tangential history");
        close(after.Fr.z,.3,"rebuild changed rolling history");
        close(DEM::DynParticleCU(d.bDynParticlesCU[0]).v.y,1,"rebuild changed half velocity");
        close(real3(d.bA[0]).x,3,"rebuild changed translation acceleration");
        close(real3(d.bWdot[0]).z,8,"rebuild changed rotation acceleration");
    }
    // A gap must not lose a retained force history through a zero Hertz
    // stiffness during the host/device rebuild round trip.
    auto body=DEM::DynParticleCU(d.bDynParticlesCU[1]);body.x.x=2.01;
    d.bDynParticlesCU[1]=body;d.bVertsCU[1]=body.x;
    DEM::ResetLE<<<1,32>>>(d.pParticlesCU,d.pDynParticlesCU,d.pInteractons,d.pComInteractons,d.pdemaux,d.pExtraParams);
    d.pForceVV<<<1,32>>>(d.pInteractons,d.pComInteractons,d.pDynInteractonsVV,d.pParticlesCU,d.pDynParticlesCU,d.pdemaux,d.pExtraParams);
    RebuildLeesEdwardsContacts(d,data);
    after=d.bDynInteractonsVV[0];
    close(after.Ft.y,.7,"gap rebuild lost custom tangential force history");
    close(after.Fr.z,.3,"gap rebuild lost custom rolling history");
    require(!after.InContact,"gap rebuild changed activation state");
    cudaFree(d.pExtraParams);d.pExtraParams=nullptr;
}

static void fast_contact(int kind)
{
    UserData data={};DEM::Domain d(&data);
    if(kind==0) {
        d.AddSphere(-1,Vec3_t(0,0,0),.1,1000);
        d.AddSphere(-2,Vec3_t(.24,0,0),.1,1000);
        d.Particles[0]->v=Vec3_t(300,0,0);d.Particles[1]->v=Vec3_t(-300,0,0);
    } else if(kind==1) {
        d.AddSphere(-1,Vec3_t(-.12,0,-4.95),.1,1000);
        d.AddSphere(-2,Vec3_t(.12,0,4.95),.1,1000);
        data.shear_rate=120;
    } else {
        d.AddSphere(-1,Vec3_t(-.11005,0,-4.9995),.1,1000);
        d.AddSphere(-2,Vec3_t(.11005,0,4.9995),.1,1000);
        d.Particles[0]->v=Vec3_t(85,0,0);d.Particles[1]->v=Vec3_t(-85,0,0);
        data.shear_rate=4;
    }
    prepare(d,data,1e-4);
    require(d.demaux.nvvint==0,"fast contact must start outside the list");
    RefreshDeviceAccelerations(d);
    require(AdvanceLeesEdwardsStep(d,data),"post-drift list was not rebuilt before force");
    require(cudaDeviceSynchronize()==cudaSuccess,"LE step failed");
    require(d.demaux.nvvint==1,"post-drift contact is missing");
    auto c=DEM::ComInteractonCU(d.bComInteractons[0]);
    double distance=kind==0?.18:(kind==1?sqrt(.12*.12+.1*.1):sqrt(.1991*.1991+.001*.001));
    double overlap=.2-distance;
    close(c.NormalOverlap,overlap,"post-drift oracle overlap");
    close(norm(c.Fnnet),1000*overlap,"first contact received no force");
    auto q0=DEM::DynParticleCU(d.bDynParticlesCU[0]);auto q1=DEM::DynParticleCU(d.bDynParticlesCU[1]);
    close(norm(q0.F+q1.F),0,"fast-contact action reaction");
    close(d.Time,1e-4,"step clock advanced twice");
    cudaFree(d.pExtraParams);d.pExtraParams=nullptr;
}

static void wrap_reference()
{
    UserData data={};DEM::Domain d(&data);
    d.AddSphere(-1,Vec3_t(0,0,4.95),.1,1000);
    d.Particles[0]->v=Vec3_t(0,0,2000);
    prepare(d,data,1e-4);
    DEM::VerletStep1<<<1,32>>>(d.pVertsCU,d.pParticlesCU,d.pDynParticlesCU,d.pA,d.pdemaux);
    DEM::WrapZ<<<1,32>>>(d.pDynParticlesCU,d.pParticlesCU,d.pVertsCU,d.pVertsoCU,d.pdemaux,
                         static_cast<LeExtra const *>(d.pExtraParams));
    auto q=DEM::DynParticleCU(d.bDynParticlesCU[0]);
    real3 v=d.bVertsCU[0],ref=d.bVertsoCU[0];
    close(q.x.z,-4.85,"LE wrapped center");close(v.z,q.x.z,"LE vertices disagree with center");
    close(v.z-ref.z,.2,"wrap erased physical displacement");
    cudaFree(d.pExtraParams);d.pExtraParams=nullptr;
}
int main() {
    for(int law=0;law<2;++law) { contact_outputs(law,false,false);contact_outputs(law,false,true);history_rebuild(law); }
    contact_outputs(0,true,false);
    for(int kind=0;kind<3;++kind) fast_contact(kind);
    wrap_reference();
    printf("SHEARBOX_REGRESSION_PASS outputs=5 histories=2 fast_contacts=3 wraps=1\n");
}
