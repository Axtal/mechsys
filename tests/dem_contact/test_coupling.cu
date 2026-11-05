// GPL-3.0-or-later. Tiny endpoint regressions; not a new physics certificate.
#include <mechsys/lbmdem/experimental/benchmark_coupling.h>
#include <fstream>
#include <iomanip>
#include <array>
namespace B=LBMDEM::Experimental::Benchmark;
static std::array<long double,4> budget(LBMDEM::Domain &d,B::M::Units const &u){
    auto f=B::read_cell_populations(d.LBMDOM.pF,d.LBMDOM.Ncells);
    std::array<long double,4> out{};
    for(size_t i=0;i<d.LBMDOM.Ncells;i++)for(int k=0;k<15;k++){
        double v=f[15*i+k];B::need(std::isfinite(v)&&v>0,"positive populations");
        out[0]+=static_cast<long double>(u.mass)*v;
        for(int j=0;j<3;j++)out[1+j]+=static_cast<long double>(u.impulse)*v*d.LBMDOM.C[k](j);
    }
    auto dyn=B::read(d.DEMDOM.pDynParticlesCU,1);auto *p=d.DEMDOM.Particles[0];
    out[1]+=static_cast<long double>(p->Props.m)*dyn[0].v.x;
    out[2]+=static_cast<long double>(p->Props.m)*dyn[0].v.y;
    out[3]+=static_cast<long double>(p->Props.m)*dyn[0].v.z;
    return out;
}
int main(int argc,char **argv)try {
    B::need(argc==2&&sizeof(real)==8,"usage: coupling_smoke SNAPSHOT_PREFIX; FP64 required");
    bool const preflight=std::string(argv[1])=="--preflight";
    for(int test=0;test<3;test++){
        double a=test? .4:1,h=test? .016:1,nu=.1*a*a/h;
        LBMDEM::Domain d(D3Q15,nu,iVec3_t(12,12,12),a,h,0);
        d.Fconv=test==2?2.5:1;d.Alpha=.1*a;
        d.PeriodicX=d.PeriodicY=d.PeriodicZ=true;
        auto &m=d.DEMDOM;auto &l=d.LBMDOM;m.UseVelocityVerlet=true;
        // Match the declared smooth single-fluid fixture and the qualified
        // benchmark driver, not FLBM's stock Sc=0.17 constructor default.
        l.Nproc=1;l.Sc=0;l.G[0]=l.Gmix=0;
        m.AddSphere(-1,a*Vec3_t(6.13,6.19,6.23),2*a,1);
        auto *p=m.Particles[0];p->PersistentId=101;
        if(test){p->v=(a/h)*Vec3_t(.0003,-.0001,.0002);p->Ff=Vec3_t(.001,0,0);}
        Vec3_t rest(0,0,0);
        for(int z=0;z<12;z++)for(int y=0;y<12;y++)for(int x=0;x<12;x++)
            l.Initialize(0,iVec3_t(x,y,z),1,rest);
        B::Options o;o.selection=B::Selection::ImplicitVV;o.max_intervals=2;
        if(test==2)o.density_sum=B::DensitySum::CompensatedV1;
        if(preflight){B::Session::validate(d,o);std::cout<<"TB_PREFLIGHT_PASS case="<<test<<" solver_steps=0"<<std::endl;continue;}
        B::Session s(d,o);s.initialize();B::M::Units u(a,h,d.Fconv);
        long double maxmass=0,maxmomentum=0;
        for(int n=1;n<=2;n++){
            auto before=budget(d,u);s.advance_one();auto after=budget(d,u);
            maxmass=std::max(maxmass,std::abs(after[0]-before[0]));
            for(int j=0;j<3;j++)maxmomentum=std::max(maxmomentum,
                std::abs(after[j+1]-before[j+1]-static_cast<long double>(h)*p->Ff(j)));
            B::need(maxmass<=1e-10L&&maxmomentum<=1e-10L,"unchanged physical mass/momentum cap");
            B::need(m.iter==size_t(n)&&d.iter==size_t(n)&&d.Time==n*h,"common clocks");
        }
        s.download();auto f=B::read_cell_populations(l.pF,l.Ncells);auto rho=B::read(l.pRho,l.Ncells);
        for(double r:rho)B::need(std::isfinite(r)&&r>0,"positive density");
        auto dyn=B::read(m.pDynParticlesCU,1)[0];
        double particle[]={dyn.x.x,dyn.x.y,dyn.x.z,dyn.v.x,dyn.v.y,dyn.v.z,
                           dyn.w.x,dyn.w.y,dyn.w.z,dyn.Q.w,dyn.Q.x,dyn.Q.y,dyn.Q.z};
        for(double x:particle)B::need(std::isfinite(x),"finite particle state");
        std::string name=std::string(argv[1])+"-"+std::to_string(test);
        std::ofstream fluid(name+"-fluid.bin",std::ios::binary),body(name+"-particle.bin",std::ios::binary);
        fluid.write(reinterpret_cast<char const*>(f.data()),f.size()*sizeof(double));
        body.write(reinterpret_cast<char const*>(particle),sizeof(particle));
        B::need(bool(fluid)&&bool(body),"endpoint output");
        std::cout<<std::setprecision(17)<<"TB_COUPLING_PASS {\"case\":"<<test
                 <<",\"vv\":1,\"steps\":2,\"mass\":"<<double(maxmass)
                 <<",\"momentum\":"<<double(maxmomentum)<<"}"<<std::endl;
    }
    return 0;
}MECHSYS_CATCH
