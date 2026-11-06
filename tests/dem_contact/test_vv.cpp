// Repeated VV rebuild/reference-state regression; compile both with and without USE_CUDA.
#include <mechsys/dem/domain.h>
#include <cstring>
#include <stdexcept>
#include <vector>

static void require(bool ok, char const * message)
{
    if (!ok) throw std::runtime_error(message);
}

static void material(DEM::Domain & d)
{
    for (size_t i=0;i<d.Particles.Size();i++)
    {
        d.Particles[i]->Props.Kn=1000.0;
        d.Particles[i]->Props.Kt=4.0;
        d.Particles[i]->Props.Gn=0.0;
        d.Particles[i]->Props.Gt=0.0;
        d.Particles[i]->Props.Beta=0.5;
        d.Particles[i]->Props.Mu=0.4;
    }
}

#ifdef USE_CUDA
template<class T>
static bool same(thrust::host_vector<T> const & a, thrust::host_vector<T> const & b)
{
    return a.size()==b.size() && (a.empty() || std::memcmp(a.data(),b.data(),a.size()*sizeof(T))==0);
}

static void narrow_rebuild(bool expect_stale)
{
    DEM::Domain d;
    d.UseVelocityVerlet=true;
    d.Alpha=0.025;
    d.AddSphere(-1,Vec3_t(0,0,0),0.1,1000.0);
    d.AddSphere(-2,Vec3_t(0.199,0,0),0.1,1000.0);
    material(d);
    d.Solve(1e-8,1e-8,1.0,NULL,NULL,NULL,false,1);
    require(d.bDynInteractonsVV.size()==1,"Expected one persistent contact");
    thrust::host_vector<DEM::DynParticleCU> dyn=d.bDynParticlesCU;
    thrust::host_vector<real3> vertices=d.bVertsCU;
    for (size_t i=0;i<dyn.size();i++)
    {
        dyn[i].x.x+=0.04;
        dyn[i].xb.x+=0.04;
        dyn[i].v=make_real3(1,2,3);
        dyn[i].w=make_real3(4,5,6);
        vertices[i].x+=0.04;
    }
    d.bDynParticlesCU=dyn;
    d.bVertsCU=vertices;
    thrust::host_vector<real3> accel(d.bA.size(),make_real3(0.25,0.5,0.75));
    thrust::host_vector<real3> angular(d.bWdot.size(),make_real3(1.25,1.5,1.75));
    d.bA=accel;
    d.bWdot=angular;
    thrust::host_vector<DEM::DynInteractonCU> history=d.bDynInteractonsVV;
    history[0].Ft=make_real3(0,8,16);
    history[0].Fr=make_real3(0,2,4);
    history[0].PrevDelta=0.001;
    history[0].InContact=true;
    d.bDynInteractonsVV=history;
    auto before=d.PairtoCInt.begin()->second;
    d.UpdateContactsDevicePreForce();
    thrust::host_vector<DEM::DynParticleCU> dyn_after=d.bDynParticlesCU;
    thrust::host_vector<real3> a_after=d.bA,w_after=d.bWdot,ref=d.bVertsoCU;
    thrust::host_vector<DEM::DynInteractonCU> h_after=d.bDynInteractonsVV;
    thrust::host_vector<DEM::ComInteractonCU> common=d.bComInteractons;
    require(same(dyn,dyn_after),"Pre-force rebuild changed particle state");
    require(same(accel,a_after)&&same(angular,w_after),"Pre-force rebuild changed A/Wdot");
    require(before==d.PairtoCInt.begin()->second,"Persistent physical contact identity changed");
    require(h_after.size()==1 && std::memcmp(&history[0].Ft,&h_after[0].Ft,sizeof(real3))==0 &&
            std::memcmp(&history[0].Fr,&h_after[0].Fr,sizeof(real3))==0 &&
            history[0].PrevDelta==h_after[0].PrevDelta && history[0].InContact==h_after[0].InContact,
            "Pre-force rebuild changed stored tangential/rolling/contact history");
    bool clear=norm(common[0].Fnnet)==0 && norm(common[0].Ftnet)==0 &&
               norm(common[0].Fndpot)==0 && norm(common[0].Ftdpot)==0 && norm(common[0].Fther)==0 &&
               common[0].Epot==0 && common[0].dEvis==0 && common[0].dEfric==0;
    printf("V1_NEW_CONTACT_ACCUMULATORS zero_initialized=%d\n",int(clear));
    if(!expect_stale) require(clear,"New contact accumulators were not initialized after pre-force rebuild");
    double residual=0;
    for(size_t i=0;i<vertices.size();i++) residual=std::max(residual,double(norm(vertices[i]-ref[i])));
    printf("V1_REFERENCE_DIRECT residual=%.17g state_unchanged=1 acceleration_unchanged=1 history_unchanged=1\n",residual);
    require(expect_stale ? residual>d.Alpha : residual==0.0,"Unexpected reference-reset result");
    if (!expect_stale)
    {
        thrust::host_vector<real> maxd=d.bMaxDCU;
        for(size_t i=0;i<maxd.size();i++) require(maxd[i]==0,"Displacement cache was not cleared");
        d.UpdateContactsDevicePreForce();
        thrust::host_vector<real3> ref2=d.bVertsoCU;
        require(same(ref,ref2),"Repeated reference reset is not idempotent");
    }
}

static void vertex_rotation()
{
    std::vector<real3> old,current;
    for(int x : {-1,1}) for(int y : {-1,1}) for(int z : {-1,1})
    {
        old.push_back(make_real3(0.5*x,0.5*y,0.5*z));
        current.push_back(make_real3(0.5*x*cos(0.1)-0.5*y*sin(0.1),
                                    0.5*x*sin(0.1)+0.5*y*cos(0.1),0.5*z));
    }
    thrust::device_vector<real3> base(old.begin(),old.end()),v(current.begin(),current.end());
    thrust::device_vector<real> delta(old.size());
    DEM::dem_aux aux={}; aux.nverts=old.size(); aux.dt=1e-4;
    DEM::dem_aux * device_aux;
    require(cudaMalloc(&device_aux,sizeof(aux))==cudaSuccess,"Device allocation failed");
    cudaMemcpy(device_aux,&aux,sizeof(aux),cudaMemcpyHostToDevice);
    DEM::MaxD<<<1,32>>>(thrust::raw_pointer_cast(v.data()),thrust::raw_pointer_cast(base.data()),
                       thrust::raw_pointer_cast(delta.data()),device_aux);
    thrust::host_vector<real> observed=delta;
    double maximum=*std::max_element(observed.begin(),observed.end());
    double expected=sqrt(2.0)*sin(0.05);
    cudaFree(device_aux);
    require(fabs(maximum-expected)<1e-14,"Vertex rotation displacement was not measured");
    printf("V1_VERTEX_ROTATION center_motion=0 maximum_vertex_motion=%.17g expected=%.17g\n",maximum,expected);
}
#endif

static void zero_step_history()
{
    DEM::Domain d;
    d.UseVelocityVerlet=true;d.Alpha=.025;
    d.AddSphere(-1,Vec3_t(0,0,0),.1,1000);
    d.AddSphere(-2,Vec3_t(.199,0,0),.1,1000);
    material(d);d.Particles[0]->v=Vec3_t(0,.1,0);
    d.Particles[0]->w=Vec3_t(0,0,1);
    d.Solve(0,1e-4,1,NULL,NULL,NULL,false,1);
    auto *c=dynamic_cast<DEM::CInteractonSphere *>(d.PairtoCInt.begin()->second);
    require(c && d.iter==0,"Expected zero physical startup steps");
    require(norm(c->Fdvv)==0 && norm(c->Fdr)==0,
            "Zero-step bootstrap advanced tangential or rolling history");
    require(c->Nc==1 && fabs(c->NormalOverlap-.001)<1e-14,
            "Entry-state overlap output was lost during upload/download");
    printf("VV_ZERO_STEP_HISTORY_PASS\n");
}

int main(int argc,char ** argv) try
{
    zero_step_history();
    bool expect_stale=argc>1&&std::strcmp(argv[1],"--expect-stale")==0;
#ifdef USE_CUDA
    narrow_rebuild(expect_stale);
    vertex_rotation();
#endif
    DEM::Domain d;
    d.UseVelocityVerlet=true;
    d.Alpha=0.025;
#ifdef USE_CUDA
#endif
    d.AddSphere(-1,Vec3_t(0,0,0),0.1,1000.0);
    d.AddSphere(-2,Vec3_t(0,0.6,0),0.1,1000.0);
    material(d);
    for(size_t i=0;i<d.Particles.Size();i++) d.Particles[i]->v=Vec3_t(100,0,0);
    d.Solve(9.5e-4,1e-4,1.0,NULL,NULL,NULL,false,1);
    double residual=0,reference=0;
#ifdef USE_CUDA
    thrust::host_vector<real3> v=d.bVertsCU,ref=d.bVertsoCU;
    for(size_t i=0;i<v.size();i++) residual=std::max(residual,double(norm(v[i]-ref[i])));
    reference=ref[0].x;
#else
    for(size_t i=0;i<d.Particles.Size();i++) residual=std::max(residual,d.Particles[i]->MaxDisplacement());
    reference=(*d.Particles[0]->Vertso[0])(0);
#endif
    printf("V1_REFERENCE_REPEAT steps=%zu position=%.17g reference=%.17g residual=%.17g\n",
           d.iter,d.Particles[0]->x(0),reference,residual);
    require(d.iter==10,"Unexpected integration step count");
#ifdef USE_CUDA
    if(expect_stale) require(residual>d.Alpha && fabs(reference)<1e-14,"Original stale-reference defect was not reproduced");
    else
#endif
        require(fabs(reference-0.09)<1e-13 && fabs(residual-0.01)<1e-13,"References did not track the latest completed rebuild");

    DEM::Domain collision;
    collision.UseVelocityVerlet=true;
    collision.Alpha=0.001;
#ifdef USE_CUDA
#endif
    collision.AddSphere(-1,Vec3_t(0,0,0),0.1,1000.0);
    collision.AddSphere(-2,Vec3_t(0.24,0,0),0.1,1000.0);
    material(collision);
    collision.Particles[0]->v=Vec3_t(300,0,0);
    collision.Particles[1]->v=Vec3_t(-300,0,0);
    collision.Solve(1e-4,1e-4,1.0,NULL,NULL,NULL,false,1);
    double f1=collision.Particles[0]->F(0),f2=collision.Particles[1]->F(0);
    require(fabs(f1+20.0)<1e-10 && fabs(f2-20.0)<1e-10,"Wrong force in fast VV contact");
    printf("V1_FAST_CONTACT F1=%.17g F2=%.17g expected_abs_force=20 action_reaction=%.17g\n",f1,f2,f1+f2);
#ifdef USE_CUDA
    DEM::ComInteractonCU c=collision.bComInteractons[0];
    printf("V1_COMMON_FORCE normal_x=%.17g expected=20\n",double(c.Fnnet.x));
    if(!expect_stale) require(fabs(c.Fnnet.x-20.0)<1e-10,"Wrong accumulated contact normal force");
#endif
    return 0;
}
MECHSYS_CATCH
