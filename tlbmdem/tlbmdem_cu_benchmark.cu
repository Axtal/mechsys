// GPL-3.0-or-later. Explicit experimental API; existing examples remain intact.
#include <mechsys/lbmdem/experimental/benchmark_coupling.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
namespace B=LBMDEM::Experimental::Benchmark;
using Tree=boost::property_tree::ptree;
static std::string output;
static size_t interval=0;
static void number(Tree &t,char const *k,double x){DEM::Checkpoint::number(t,k,x);}
static Tree array(double const *x,size_t n){Tree t;for(size_t i=0;i<n;i++){Tree a;number(a,"",x[i]);t.push_back({"",a});}return t;}
static void vec(Tree &t,char const *key,real3 x){double a[3]={x.x,x.y,x.z};t.add_child(key,array(a,3));}
static void write(std::string const &name,Tree const &t){boost::property_tree::write_json(output+"/"+name,t);}
template<class T>static void binary(std::string const &name,std::vector<T> const &v){std::ofstream f(output+"/"+name,std::ios::binary);f.write(reinterpret_cast<char const*>(v.data()),v.size()*sizeof(T));B::need(bool(f),"CB: output write failed");}
static Tree state(DEM::Domain &m){
    Tree root,parts;auto d=B::read(m.pDynParticlesCU,m.Particles.Size());auto p=B::read(m.pParticlesCU,m.Particles.Size());
    auto a=B::read(m.pA,m.Particles.Size()),w=B::read(m.pWdot,m.Particles.Size());auto da=B::read(m.pdemaux,1)[0];
    for(size_t i=0;i<d.size();i++){Tree t;t.put("id",m.Particles[i]->PersistentId);t.put("device_index",p[i].Index);
        number(t,"mass",p[i].m);number(t,"radius",p[i].R);vec(t,"I",p[i].I);vec(t,"x",d[i].x);vec(t,"xb",d[i].xb);vec(t,"v",d[i].v);vec(t,"w",d[i].w);
        vec(t,"wa",d[i].wa);vec(t,"wb",d[i].wb);vec(t,"F",d[i].F);vec(t,"Flbm",d[i].Flbm);vec(t,"T",p[i].T);vec(t,"A",a[i]);vec(t,"Wdot",w[i]);
        vec(t,"Ff",p[i].Ff);vec(t,"Tf",p[i].Tf);t.put("fixed",int(p[i].vxf&&p[i].vyf&&p[i].vzf&&p[i].wxf&&p[i].wyf&&p[i].wzf));
        double q[4]={d[i].Q.w,d[i].Q.x,d[i].Q.y,d[i].Q.z};t.add_child("Q",array(q,4));parts.push_back({"",t});}
    root.add_child("particles",parts);root.put("contacts",m.demaux.nvvint);root.put("step",m.iter);number(root,"time",m.Time);
    root.put("device_step",da.iter);number(root,"device_time",da.Time);root.put("vv",int(m.UseVelocityVerlet));root.put("checkpoint_ready",int(m.CheckpointReady));return root;
}
static void fluid(LBMDEM::Domain &d,std::string const &label){
    auto &l=d.LBMDOM;binary(label+"-f.bin",B::read_cell_populations(l.pF,l.Ncells));auto r=B::read(l.pRho,l.Ncells);auto u=B::read(l.pVel,l.Ncells);std::vector<double> m;
    for(size_t i=0;i<r.size();i++){m.push_back(r[i]);m.push_back(u[i].x);m.push_back(u[i].y);m.push_back(u[i].z);}binary(label+"-moments.bin",m);
}
static void exchange(B::Session &s,LBMDEM::Domain &d,std::string const &label){
    auto claims=s.claims();auto cells=s.cells();auto bodies=s.bodies();auto impulses=s.impulses();Tree root,bs,js;
    std::ofstream f(output+"/"+label+"-claims.csv");f<<std::setprecision(17)<<"cell,owner,gamma,B";
    for(auto const &n:std::vector<std::pair<std::string,int>>{{"r",3},{"surface",3},{"omega",15},{"solid",15},{"D",36},{"b",6},{"J",3},{"Kw",3},{"Kb",3}})for(int j=0;j<n.second;j++)f<<','<<n.first<<j;f<<",trace_error\n";
    size_t active=0;for(auto const &c:claims)if(c.owner>=0){active++;f<<c.cell<<','<<c.owner<<','<<c.gamma<<','<<c.B;
        for(auto const &block:std::vector<std::pair<double const*,int>>{{c.lever,3},{c.surface,3},{c.omega,15},{c.solid,15},{c.math.D,36},{c.math.b,6},{c.math.j,3},{c.math.kw,3},{c.math.kb,3}})for(int j=0;j<block.second;j++)f<<','<<block.first[j];f<<','<<c.math.trace.error<<'\n';}B::need(bool(f),"CB: claim output");
    std::vector<double> data;for(auto const &c:cells){data.push_back(c.gamma_sum);data.push_back(c.Gamma);data.push_back(c.B);data.push_back(c.count);}binary(label+"-cells.bin",data);
    for(auto const &b:bodies){Tree x;x.add_child("M",array(b.M,6));x.add_child("z",array(b.z,6));x.add_child("Q",array(b.Q,4));x.add_child("D",array(b.D,36));x.add_child("b",array(b.b,6));x.add_child("t",array(b.t,6));x.put("owned",b.owned);x.put("error",b.trace.error);bs.push_back({"",x});}
    for(auto const &j:impulses){Tree x;x.add_child("J",array(j.J,3));x.add_child("Kw",array(j.Kw,3));x.add_child("Kb",array(j.Kb,3));js.push_back({"",x});}
    root.add_child("bodies",bs);root.add_child("impulses",js);root.put("active_claims",active);root.put("candidate_capacity",s.capacity());write(label+"-exchange.json",root);
    binary(label+"-collision-f.bin",B::read_cell_populations(d.LBMDOM.pFtemp,d.LBMDOM.Ncells));
}
static size_t positive(char const *s){std::string z(s);B::need(!z.empty()&&z.find_first_not_of("0123456789")==std::string::npos,"CB: unsigned integer required");size_t n=std::stoull(z);B::need(n>0,"CB: positive integer required");return n;}
static int guards(){
    LBMDEM::Domain d(D3Q15,.1,iVec3_t(16,16,16),1,1,0);auto &m=d.DEMDOM;m.UseVelocityVerlet=true;
    m.AddSphere(-1,Vec3_t(8,8,8),2,1);B::Options o;o.selection=B::Selection::ImplicitVV;o.max_intervals=2;
    int n=0;auto reject=[&](auto action){bool rejected=false;try{action();}catch(std::runtime_error const&){rejected=true;}B::need(rejected,"CB: guard not rejected");n++;};
    reject([&]{B::Session s(d,B::Options{});});o.max_intervals=0;reject([&]{B::Session s(d,o);});o.max_intervals=2;
    o.selection=B::Selection::ExplicitPV;reject([&]{B::Session s(d,o);});o.selection=B::Selection::ImplicitVV;
    m.Particles[0]->vxf=true;reject([&]{B::Session::validate(d,o);});m.Particles[0]->vxf=false;
    d.LBMDOM.Sc=.1;reject([&]{B::Session::validate(d,o);});d.LBMDOM.Sc=0;
    d.Fconv=0;reject([&]{B::Session::validate(d,o);});d.Fconv=1;
    m.Particles[0]->FixVeloc();m.Particles[0]->w=Vec3_t(0,0,.1);reject([&]{B::Session::validate(d,o);});m.Particles[0]->w=OrthoSys::O;
    o.max_claim_bytes=1;reject([&]{B::Session s(d,o);s.initialize();});
    B::need(!m.Initialized&&m.iter==0,"CB: guards must precede bootstrap");std::cout<<"CB_GUARDS_PASS "<<n<<std::endl;return 0;
}
int main(int argc,char **argv)try{
    if(argc==2&&std::string(argv[1])=="--guards")return guards();
    B::need(argc==7||argc==8,"usage: tlbmdem_cu_benchmark vv|pv 01|02|03|tiny:FIXTURE STEPS REPORT_EVERY NEW_OUTPUT [audit]");
    // argv[6] is the explicit opt-in token, not a silently enabled default.
    B::need(std::string(argv[6])=="--experimental","CB: --experimental required");
    std::string method=argv[1],fixture=argv[2];size_t steps=positive(argv[3]),cadence=positive(argv[4]);output=argv[5];
    bool tiny=fixture.rfind("tiny:",0)==0,audit=argc==8&&std::string(argv[7])=="audit";std::string test=tiny?fixture.substr(5):fixture;
    std::set<std::string> tests={"unit","drag","moving","buoyancy","prescribed","conversion"};
    B::need((method=="vv"||method=="pv")&&(tiny?tests.count(test):(test=="01"||test=="02"||test=="03")),"CB: supported mode/fixture");
    B::need(argc==7||audit,"CB: unknown output mode");B::need(!audit||(tiny&&steps==2&&cadence==1),"CB: full audit restricted to registered two-interval fixtures");
    B::need(steps<=100000000,"CB: finite clock bound");B::need(!std::filesystem::exists(output),"CB: fresh output only");std::filesystem::create_directories(output);
    bool drag=test=="drag"||test=="01",moving=test=="moving"||test=="02",buoy=test=="buoyancy"||test=="03",prescribed=test=="prescribed";
    double dx=test=="unit"?1:.4,dt=test=="unit"?1:.016,C=dx/dt;
    int nx=tiny?16:test=="01"?301:test=="02"?200:100,ny=tiny?16:test=="01"?61:nx,nz=ny;
    double nu=tiny?(buoy?.45:.3)*dx*dx/(3*dt):test=="01"?1:test=="02"?1.6:1.5;
    double R=tiny?2*dx:test=="01"?3.6:test=="02"?1.8:7.2;
    LBMDEM::Domain d(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt,0);auto &m=d.DEMDOM;auto &l=d.LBMDOM;
    m.UseVelocityVerlet=method=="vv";l.Nproc=1;l.Sc=0;l.G[0]=l.Gmix=0;d.Fconv=test=="conversion"?2.5:1;d.Alpha=2*dx;
    d.PeriodicX=d.PeriodicY=d.PeriodicZ=!(drag||buoy);
    Vec3_t center(.5*dx*nx,.5*dx*ny,.5*dx*nz);if(tiny)center+=dx*Vec3_t(.13,.19,.23);
    m.AddSphere(-1,center,R,1);auto *p=m.Particles[0];p->PersistentId=101;
    if(drag||buoy)p->FixVeloc();if(prescribed)p->FixVeloc(.0002*C,0,0);
    if(moving){p->Ff=p->Props.m*Vec3_t(.1,0,0);if(tiny){p->Tf=Vec3_t(1e-5,-2e-5,3e-5);p->w=Vec3_t(.0001,-.0002,.0003)/dt;p->Q=Quaternion_t(cos(.17),sin(.17),0,0);}}
    if(tiny&&!drag&&!buoy&&!prescribed)p->v=C*Vec3_t(.0003,-.0001,.0002);
    double rhof=buoy?2:1,Cs2=C*C/3,g=buoy?1e-5*Cs2/(nz*dx):0;
    B::Options options;options.selection=m.UseVelocityVerlet?B::Selection::ImplicitVV:B::Selection::ExplicitPV;options.max_intervals=steps;
    options.fluid_acceleration=drag?Vec3_t(.03,0,0):Vec3_t(0,0,-g);
    double rho0=buoy?rhof*g*(nz*dx)/(Cs2*(exp(nz*dx*g/Cs2)-1)):rhof;
    std::vector<unsigned char> mask;
    for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
        double rho=buoy?rho0*exp((nz*dx-dx*(z+1))*g/Cs2):rhof;
        if(buoy&&norm(Vec3_t(x*dx,y*dx,z*dx)-center)<R)rho=rhof;
        l.Initialize(0,iVec3_t(x,y,z),rho,OrthoSys::O);bool wall=(moving&&(y==0||y==ny-1||z==0||z==nz-1))||(buoy&&z==0);
        l.IsSolid[0][x][y][z]=wall;mask.push_back(wall);
    }
    Tree input;input.put("method",method);input.put("fixture",fixture);input.put("steps",steps);input.put("report_every",cadence);input.put("audit",int(audit));
    input.put("nx",nx);input.put("ny",ny);input.put("nz",nz);number(input,"dx",dx);number(input,"dt",dt);number(input,"Fconv",d.Fconv);number(input,"nu",nu);number(input,"tau",l.Tau[0]);number(input,"rhof",rhof);number(input,"gravity",g);
    double acc[3]={options.fluid_acceleration(0),options.fluid_acceleration(1),options.fluid_acceleration(2)},per[3]={double(d.PeriodicX),double(d.PeriodicY),double(d.PeriodicZ)};
    input.add_child("fluid_acceleration",array(acc,3));input.add_child("periodic_particles",array(per,3));input.put("classification","EXPERIMENTAL_FEATURE_VALIDATION_NOT_STEADY_ACCURACY_OR_PERFORMANCE");write("INPUT.json",input);binary("mask.bin",mask);
    B::Session s(d,options);s.initialize();if(audit){write("0-after.json",state(m));fluid(d,"0-after");}
    std::ofstream report(output+"/report.csv");report<<std::setprecision(17)<<"interval,time,x,y,z,vx,vy,vz,Fx,Fy,Fz,Flbmx,Flbmy,Flbmz,Tx,Ty,Tz,reference\n";
    s.hook=[&](char const *stage){if(!audit)return;std::string label=std::to_string(interval)+"-"+stage;
        if(std::string(stage)=="exchange"){exchange(s,d,std::to_string(interval));return;}write(label+".json",state(m));
        if(std::string(stage)=="before"||std::string(stage)=="after")fluid(d,label);
        if(std::string(stage)=="after"){std::cout<<"CB_FRAME "<<interval<<std::endl;std::string ack;B::need(bool(std::getline(std::cin,ack))&&ack=="ACCEPT","CB: independent audit ACK required");}};
    for(interval=1;interval<=steps;interval++){s.advance_one();if(interval%cadence==0||interval==steps){s.download();report<<interval<<','<<d.Time;
        for(auto const &a:{p->x,p->v,p->F,p->Flbm,p->T})for(int j=0;j<3;j++)report<<','<<a(j);
        double reference=buoy?p->Props.V*rhof*g:moving?6*M_PI*nu*R*norm(p->v):0;report<<','<<reference<<'\n';B::need(bool(report),"CB: report write failed");}}
    Tree done;done.put("completed",s.completed());done.put("retries",0);done.put("steady_accuracy_certified",false);write("COMPLETED.json",done);return 0;
}MECHSYS_CATCH
