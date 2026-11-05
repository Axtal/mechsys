// Explicit full-step VV restart. Legacy HDF5 Save/Load deliberately unchanged.
#ifndef MECHSYS_DEM_CHECKPOINT_H
#define MECHSYS_DEM_CHECKPOINT_H
#include <boost/property_tree/json_parser.hpp>
#include <iomanip>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <fcntl.h>
#include <unistd.h>

namespace DEM { namespace Checkpoint {
using Tree = boost::property_tree::ptree;
inline void need(bool ok, char const * what) {
    if (!ok) throw std::runtime_error(std::string("VV checkpoint: ")+what);
}
inline void number(Tree & t,char const * key,double x) {
    need(std::isfinite(x),"non-finite mechanical state");
    std::ostringstream s; s<<std::setprecision(17)<<x; t.put(key,s.str());
}
inline double number(Tree const & t,char const * key) {
    double x=t.get<double>(key); need(std::isfinite(x),"non-finite number"); return x;
}
inline bool flag(Tree const & t,char const * key) {
    unsigned x=t.get<unsigned>(key); need(x<=1,"invalid boolean"); return x!=0;
}
inline void vector(Tree & t,char const * key,Vec3_t const & v) {
    Tree a; for(size_t i=0;i<3;i++){Tree x; number(x,"",v(i)); a.push_back({"",x});} t.add_child(key,a);
}
inline Vec3_t vector(Tree const & t,char const * key) {
    auto const & a=t.get_child(key); need(a.size()==3,"invalid vector length");
    Vec3_t v; size_t i=0; for(auto const & x:a) v(i++)=number(x.second,""); return v;
}
inline void named_keys(Tree const & t) {
    std::set<std::string> keys;
    for(auto const & x:t) {
        need(x.first.empty() || keys.insert(x.first).second,"duplicate object field");
        named_keys(x.second);
    }
}
inline std::pair<size_t,size_t> unhash(size_t key) {
    // Exact inverse of the baseline Cantor local-feature key, no FP sqrt.
    uint64_t lo=0,hi=uint64_t(1)<<33;
    while(lo+1<hi) { uint64_t m=(lo+hi)/2;
        if((__uint128_t)m*(m+1)/2<=key) lo=m; else hi=m; }
    size_t j=key-size_t((__uint128_t)lo*(lo+1)/2);
    need(j<=lo,"invalid feature key"); return {lo-j,j};
}
inline void features(Tree & t,char const * key,ListContacts_t const & list,FrictionMap_t const & map) {
    Tree data,active,hist;
    for(auto const & f:list) {Tree x; x.put("a",f.first); x.put("b",f.second); active.push_back({"",x});}
    std::vector<size_t> keys; for(auto const & f:map) keys.push_back(f.first); std::sort(keys.begin(),keys.end());
    for(auto k:keys) {auto f=unhash(k); Tree x; x.put("a",f.first); x.put("b",f.second); vector(x,"displacement",map.at(k)); hist.push_back({"",x});}
    data.put("active_count",list.Size()); data.put("history_count",map.size());
    data.add_child("active",active); data.add_child("history",hist); t.add_child(key,data);
}
inline void features(Tree const & t,char const * key,ListContacts_t & list,FrictionMap_t & map,size_t na,size_t nb) {
    auto const & d=t.get_child(key); auto const & active=d.get_child("active"); auto const & hist=d.get_child("history");
    need(active.size()==d.get<size_t>("active_count") && hist.size()==d.get<size_t>("history_count"),"feature count mismatch");
    std::set<size_t> seen;
    for(auto const & f:active) {
        size_t a=f.second.get<size_t>("a"),b=f.second.get<size_t>("b");
        need(a<na && b<nb && seen.insert(HashFunction(a,b)).second,"invalid/duplicate active feature"); list.Push({int(a),int(b)});
    }
    for(auto const & f:hist) {
        size_t a=f.second.get<size_t>("a"),b=f.second.get<size_t>("b");
        need(a<na && b<nb && map.emplace(HashFunction(a,b),vector(f.second,"displacement")).second,"invalid/duplicate feature history");
    }
}
inline void write_new(char const * filename,Tree const & t) {
    std::ostringstream stream; boost::property_tree::write_json(stream,t);
    std::string bytes=stream.str(),temp=std::string(filename)+".tmp";
    int fd=::open(temp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);
    need(fd>=0,"cannot create unique temporary file (existing evidence is never overwritten)");
    bool ok=true; size_t offset=0;
    while(offset<bytes.size()) {ssize_t n=::write(fd,bytes.data()+offset,bytes.size()-offset); if(n<=0){ok=false;break;} offset+=size_t(n);}
    ok=(::fsync(fd)==0)&&ok; ok=(::close(fd)==0)&&ok;
    // link publishes atomically and refuses an existing target, unlike rename.
    if(ok) ok=(::link(temp.c_str(),filename)==0);
    ::unlink(temp.c_str()); need(ok,"checkpoint write/publication failed; target not overwritten");
}
#ifdef USE_CUDA
inline void gpu_history(Tree & t,char const * key,thrust::device_vector<DynInteractonCU> const & buffer,bool sphere) {
    thrust::host_vector<DynInteractonCU> h=buffer; Tree a;
    for(auto const & f:h) {Tree x; x.put("contact",f.Idx);
        if(!sphere){x.put("feature_a",f.IF1);x.put("feature_b",f.IF2);}
        vector(x,"Ft",Vec3_t(f.Ft.x,f.Ft.y,f.Ft.z));
        if(sphere){vector(x,"Fr",Vec3_t(f.Fr.x,f.Fr.y,f.Fr.z)); number(x,"PrevDelta",f.PrevDelta);x.put("InContact",int(f.InContact));}
        a.push_back({"",x});}
    t.add_child(key,a);
}
inline void gpu_history(Tree const & t,char const * key,std::vector<DynInteractonCU> & h,bool sphere,size_t nc) {
    for(auto const & row:t.get_child(key)) {auto const & x=row.second; DynInteractonCU f={};
        f.Idx=x.get<size_t>("contact"); need(f.Idx<nc,"invalid GPU contact reference");
        if(!sphere){f.IF1=x.get<size_t>("feature_a");f.IF2=x.get<size_t>("feature_b");}
        Vec3_t v=vector(x,"Ft"); f.Ft=make_real3(v(0),v(1),v(2));
        if(sphere){v=vector(x,"Fr"); f.Fr=make_real3(v(0),v(1),v(2)); f.PrevDelta=number(x,"PrevDelta"); f.InContact=flag(x,"InContact");}
        h.push_back(f);}
}
inline void restore_gpu_history(thrust::device_vector<DynInteractonCU> & buffer,std::vector<DynInteractonCU> const & saved,bool sphere) {
    need(buffer.size()==saved.size(),"GPU history count mismatch"); thrust::host_vector<DynInteractonCU> h=buffer;
    for(size_t i=0;i<h.size();i++) {
        need(h[i].Idx==saved[i].Idx && (sphere || (h[i].IF1==saved[i].IF1 && h[i].IF2==saved[i].IF2)),"GPU feature identity mismatch");
        h[i].Ft=saved[i].Ft;
        if(sphere){h[i].Fr=saved[i].Fr;h[i].PrevDelta=saved[i].PrevDelta;h[i].InContact=saved[i].InContact;}
    }
    buffer=h;
}
#endif
} // Checkpoint

inline void Domain::AssignPersistentIds() {
    std::set<uint64_t> ids;
    for(auto * p:Particles) if(p->PersistentId) {
        Checkpoint::need(ids.insert(p->PersistentId).second,"duplicate particle lifetime ID");
        Checkpoint::need(p->PersistentId<UINT64_MAX,"particle ID exhausted");
        NextPersistentId=std::max(NextPersistentId,p->PersistentId+1);
    }
    for(auto * p:Particles) if(!p->PersistentId) {
        Checkpoint::need(NextPersistentId<UINT64_MAX,"particle ID exhausted"); p->PersistentId=NextPersistentId++;
    }
}

#define MS_CP_DOMAIN_D(X) X(Time) X(Dt) X(Alpha) X(Beta) X(Evis) X(Efric) X(Wext) X(Xmin) X(Xmax) X(Ymin) X(Ymax) X(Zmin) X(Zmax) X(ThermostatInteractionRange)
#define MS_CP_PROPS(X) X(Kn) X(Kt) X(Bn) X(Bt) X(Bm) X(Gn) X(Gt) X(Gv) X(Gm) X(Mu) X(eps) X(Beta) X(Eta) X(R) X(rho) X(V) X(m)
#define MS_CP_PV(X) X(x) X(xb) X(v) X(w) X(wa) X(wb) X(F) X(T) X(Ff) X(Tf) X(Flbm) X(Flbmf) X(I) X(A) X(Wdot) X(v_half) X(w_half)
#define MS_CP_PB(X) X(PropsReady) X(Eroded) X(Bdry) X(Closed) X(vxf) X(vyf) X(vzf) X(wxf) X(wyf) X(wzf) X(FixFree)
#define MS_CP_CD(X) X(NormalOverlap) X(Kn) X(Kt) X(Gn) X(Gt) X(Mu) X(Epot) X(dEvis) X(dEfric)
#define MS_CP_CV(X) X(F1) X(F2) X(T1) X(T2) X(Fnet) X(Ftnet) X(Fndpot) X(Ftdpot) X(Fther)
#define MS_CP_SD(X) X(beta) X(eta) X(Bn) X(Bt) X(eps) X(PrevDelta)
#define MS_CP_SV(X) X(Fdvv) X(Fdr) X(vt_prev)

inline void Domain::SaveCheckpoint(char const * filename) {
    using namespace Checkpoint;
    need(UseVelocityVerlet && VVStateReady,"save requires a completed full-step VV boundary");
    need(std::numeric_limits<double>::is_iec559 && sizeof(double)==8,"requires IEEE754 binary64");
    need(ContactLaw<=1 && !MostlySpheres && BInteractons.Size()==0 && UserData==NULL,
         "unsupported shortcut/bond/custom callback-owned state; nothing was silently omitted");
    need(PairtoCInt.size()==Interactons.Size(),"inactive pair history cannot be omitted");
#ifdef USE_CUDA
    need(!CheckpointReady,"resume the loaded device state before saving another GPU checkpoint");
#endif
    for(auto * p:Particles) need(p->Tori.Size()==0 && p->Cylinders.Size()==0,"torus/cylinder checkpoint not supported");
    AssignPersistentIds(); Tree out;
    out.put("magic","MechSys full-state VV checkpoint"); out.put("version",1);out.put("precision","IEEE754-binary64-17digits");out.put("phase","completed-full-step");
#ifdef USE_CUDA
    need(sizeof(real)==8,"CUDA checkpoint requires FP64");
    need(pForceVV==(ContactLaw==0?CalcForceVV:CalcForceVV_Hertz) && pForceEE==CalcForceEEPre && pForceVF==CalcForceVFPre && pForceFV==CalcForceFVPre &&
         pReset==Reset && pExtraParams==NULL && pVerletStep1==VerletStep1 && pFinalizeVelocity==FinalizeVelocity &&
         pOrientationUpdate==OrientationUpdate && pFinalizeRotation==FinalizeRotation,"custom CUDA/coupled state is unsupported");
    need(cudaDeviceSynchronize()==cudaSuccess,"device synchronization failed"); DnLoadDevice(Nproc,true);
    thrust::host_vector<ParticleCU> ps=bParticlesCU;
    thrust::host_vector<real3> aa=bA,ww=bWdot,refs=bVertsoCU;
    for(size_t i=0;i<Particles.Size();i++) {
        Particles[i]->T=Vec3_t(ps[i].T.x,ps[i].T.y,ps[i].T.z);
        A[i]=Vec3_t(aa[i].x,aa[i].y,aa[i].z);Wdot[i]=Vec3_t(ww[i].x,ww[i].y,ww[i].z);
        for(size_t v=0;v<Particles[i]->Verts.Size();v++) {auto r=refs[Particles[i]->Nvi+v]; *Particles[i]->Vertso[v]=Vec3_t(r.x,r.y,r.z);}
    }
    out.put("backend","cuda");
    gpu_history(out,"gpu_vv",bDynInteractonsVV,true);gpu_history(out,"gpu_ee",bDynInteractonsEE,false);
    gpu_history(out,"gpu_vf",bDynInteractonsVF,false);gpu_history(out,"gpu_fv",bDynInteractonsFV,false);
#else
    out.put("backend","cpu");
#endif
    out.put("particle_count",Particles.Size());out.put("contact_count",Interactons.Size());out.put("candidate_count",ListPosPairs.Size());
    out.put("next_particle_id",NextPersistentId);out.put("iter",iter);out.put("law",ContactLaw);out.put("coulomb",SphereCoulombMode);
    out.put("RotPar",int(RotPar));out.put("tensile",int(SphereTensileCutoff));out.put("activation",int(SphereFirstContactCorrection));
#define WRITE_D(n) number(out,#n,n);
    MS_CP_DOMAIN_D(WRITE_D)
#undef WRITE_D
    Tree particles,contacts,pairs,friction;
    for(size_t i=0;i<Particles.Size();i++) {
        auto * p=Particles[i];Tree t,verts,edges,faces,q;
        t.put("id",p->PersistentId);t.put("tag",p->Tag);t.put("cluster",p->Cluster);
#define WRITE_D(n) number(t,"props." #n,p->Props.n);
        MS_CP_PROPS(WRITE_D)
#undef WRITE_D
#define WRITE_V(n) vector(t,#n,p->n);
        MS_CP_PV(WRITE_V)
#undef WRITE_V
#define WRITE_B(n) t.put(#n,int(p->n));
        MS_CP_PB(WRITE_B)
#undef WRITE_B
        number(t,"Dmax",p->Dmax);number(t,"Diam",p->Diam);number(t,"Ekin",p->Ekin);number(t,"Erot",p->Erot);
        vector(t,"domain_A",A[i]);vector(t,"domain_Wdot",Wdot[i]);
        for(size_t k=0;k<4;k++){Tree x;number(x,"",p->Q(k));q.push_back({"",x});}t.add_child("Q",q);
        for(size_t k=0;k<p->Verts.Size();k++){Tree x;vector(x,"current",*p->Verts[k]);vector(x,"reference",*p->Vertso[k]);verts.push_back({"",x});}
        for(auto const & e:p->EdgeCon){Tree x;x.put("a",e[0]);x.put("b",e[1]);edges.push_back({"",x});}
        for(auto const & f:p->FaceCon){Tree x;for(int j:f){Tree v;v.put("",j);x.push_back({"",v});}faces.push_back({"",x});}
        t.put("vertices_count",verts.size());t.put("edges_count",edges.size());t.put("faces_count",faces.size());
        t.add_child("vertices",verts);t.add_child("edges",edges);t.add_child("faces",faces);particles.push_back({"",t});
    }
    std::set<std::pair<uint64_t,uint64_t>> seen;
    for(auto * interaction:Interactons) {
        auto * c=dynamic_cast<CInteracton *>(interaction);need(c!=NULL,"unsupported interaction class");Tree t;
        auto id1=c->P1->PersistentId,id2=c->P2->PersistentId;
        need(id1!=id2 && seen.insert(std::minmax(id1,id2)).second,"duplicate physical contact");
        t.put("particle_a",id1);t.put("particle_b",id2);t.put("BothFree",int(c->BothFree));t.put("First",int(c->First));
        t.put("Nc",c->Nc);t.put("Nsc",c->Nsc);t.put("Nr",c->Nr);
#define WRITE_D(n) number(t,#n,c->n);
        MS_CP_CD(WRITE_D)
#undef WRITE_D
#define WRITE_V(n) vector(t,#n,c->n);
        MS_CP_CV(WRITE_V)
#undef WRITE_V
        if(auto * s=dynamic_cast<CInteractonSphere *>(c)) {
            t.put("kind","sphere");t.put("InContact",int(s->InContact));t.put("cohesion",int(s->cohesion));
#define WRITE_D(n) number(t,#n,s->n);
            MS_CP_SD(WRITE_D)
#undef WRITE_D
#define WRITE_V(n) vector(t,#n,s->n);
            MS_CP_SV(WRITE_V)
#undef WRITE_V
        } else {
            t.put("kind","features");
            features(t,"EE",c->Lee,c->Fdee);features(t,"VF",c->Lvf,c->Fdvf);features(t,"FV",c->Lfv,c->Fdfv);
            features(t,"VT",c->Lvt,c->Fdvt);features(t,"TV",c->Ltv,c->Fdtv);features(t,"VC",c->Lvc,c->Fdvc);features(t,"CV",c->Lcv,c->Fdcv);
        }
        contacts.push_back({"",t});
    }
    for(auto const & p:ListPosPairs){Tree x;x.put("a",Particles[p.first]->PersistentId);x.put("b",Particles[p.second]->PersistentId);pairs.push_back({"",x});}
    for(auto const & f:FricCoeff){Tree x;x.put("tag_a",f.first.first);x.put("tag_b",f.first.second);number(x,"mu",f.second);friction.push_back({"",x});}
    out.add_child("particles",particles);out.add_child("contacts",contacts);out.add_child("candidates",pairs);out.add_child("friction_overrides",friction);
    write_new(filename,out);
}

inline void Domain::LoadCheckpoint(char const * filename) {
    using namespace Checkpoint;
    need(Particles.Size()==0 && !Initialized,"load requires a fresh, empty Domain");
    Tree in;
    try {boost::property_tree::read_json(filename,in);}
    catch(boost::property_tree::json_parser_error const &) {
        throw std::runtime_error("VV checkpoint: invalid JSON/full-state file; legacy HDF5 Save/Load is not a checkpoint");
    }
    named_keys(in);
    need(in.get<std::string>("magic")=="MechSys full-state VV checkpoint" && in.get<unsigned>("version")==1,"incompatible magic/version (legacy Save is not a checkpoint)");
    need(in.get<std::string>("precision")=="IEEE754-binary64-17digits" && std::numeric_limits<double>::is_iec559 && sizeof(double)==8,"incompatible precision");
    need(in.get<std::string>("phase")=="completed-full-step","incompatible integration phase");
#ifdef USE_CUDA
    need(sizeof(real)==8 && in.get<std::string>("backend")=="cuda","incompatible backend/precision");
#else
    need(in.get<std::string>("backend")=="cpu","incompatible backend");
#endif
    auto const & rows=in.get_child("particles");auto const & contact_rows=in.get_child("contacts");
    need(rows.size()>0 && rows.size()==in.get<size_t>("particle_count") && contact_rows.size()==in.get<size_t>("contact_count") &&
         in.get_child("candidates").size()==in.get<size_t>("candidate_count"),"record count mismatch");
    Domain staged;std::vector<std::unique_ptr<Particle>> owners;std::vector<std::unique_ptr<CInteracton>> history_owners;
#define READ_D(n) staged.n=number(in,#n);
    MS_CP_DOMAIN_D(READ_D)
#undef READ_D
    need(staged.Dt>0 && staged.Alpha>0 && staged.Beta>0,"invalid timestep/skin/bin factor");
    staged.ContactLaw=in.get<size_t>("law");staged.SphereCoulombMode=in.get<size_t>("coulomb");
    need(staged.ContactLaw<=1 && staged.SphereCoulombMode<=3,"unsupported law/Coulomb mode");
    staged.Time=number(in,"Time");staged.iter=in.get<size_t>("iter");staged.NextPersistentId=in.get<uint64_t>("next_particle_id");
    staged.RotPar=flag(in,"RotPar");staged.SphereTensileCutoff=flag(in,"tensile");staged.SphereFirstContactCorrection=flag(in,"activation");
    staged.Per=Vec3_t(staged.Xmax-staged.Xmin,staged.Ymax-staged.Ymin,staged.Zmax-staged.Zmin);
    std::map<uint64_t,size_t> ids;
    staged.A.Resize(rows.size());staged.Wdot.Resize(rows.size());
    for(auto const & row:rows) {
        auto const & t=row.second;uint64_t id=t.get<uint64_t>("id");size_t i=owners.size();
        need(id>0 && id<staged.NextPersistentId && ids.emplace(id,i).second,"invalid/duplicate particle ID");
        Array<Vec3_t> v,refs;Array<Array<int>> e,f;
        for(auto const & x:t.get_child("vertices")){v.Push(vector(x.second,"current"));refs.Push(vector(x.second,"reference"));}
        need(v.Size()>0 && v.Size()==t.get<size_t>("vertices_count"),"vertex count mismatch");
        for(auto const & x:t.get_child("edges")){int a=x.second.get<int>("a"),b=x.second.get<int>("b");need(a>=0&&b>=0&&size_t(a)<v.Size()&&size_t(b)<v.Size()&&a!=b,"bad edge connectivity");Array<int> pair;pair.Push(a);pair.Push(b);e.Push(pair);}
        for(auto const & x:t.get_child("faces")){Array<int> face;for(auto const & j:x.second){int k=j.second.get_value<int>();need(k>=0&&size_t(k)<v.Size(),"bad face connectivity");face.Push(k);}need(face.Size()>=3,"short face");f.Push(face);}
        need(e.Size()==t.get<size_t>("edges_count") && f.Size()==t.get<size_t>("faces_count"),"topology count mismatch");
        auto p=std::unique_ptr<Particle>(new Particle(t.get<int>("tag"),v,e,f,OrthoSys::O,OrthoSys::O,number(t,"props.R"),number(t,"props.rho")));
        p->PersistentId=id;p->Index=i;p->Cluster=t.get<int>("cluster");
#define READ_D(n) p->Props.n=number(t,"props." #n);
        MS_CP_PROPS(READ_D)
#undef READ_D
#define READ_V(n) p->n=vector(t,#n);
        MS_CP_PV(READ_V)
#undef READ_V
#define READ_B(n) p->n=flag(t,#n);
        MS_CP_PB(READ_B)
#undef READ_B
        need(p->PropsReady && p->Props.m>0 && p->Props.R>=0 && p->I(0)>0 && p->I(1)>0 && p->I(2)>0,"invalid particle material/inertia");
        p->Dmax=number(t,"Dmax");p->Diam=number(t,"Diam");p->Ekin=number(t,"Ekin");p->Erot=number(t,"Erot");
        auto const & q=t.get_child("Q");need(q.size()==4,"invalid quaternion length");size_t k=0;double qnorm=0;
        for(auto const & x:q){p->Q(k)=number(x.second,"");qnorm+=p->Q(k)*p->Q(k);k++;}need(fabs(qnorm-1)<1e-10,"non-unit quaternion");
        for(size_t j=0;j<refs.Size();j++) *p->Vertso[j]=refs[j];
        staged.A[i]=vector(t,"domain_A");staged.Wdot[i]=vector(t,"domain_Wdot");
        staged.Particles.Push(p.get());owners.push_back(std::move(p));
    }
    std::set<std::pair<uint64_t,uint64_t>> seen;
    for(auto const & row:contact_rows) {
        auto const & t=row.second;auto a=t.get<uint64_t>("particle_a"),b=t.get<uint64_t>("particle_b");
        need(a!=b && ids.count(a) && ids.count(b) && seen.insert(std::minmax(a,b)).second,"invalid/duplicate contact identity");
        Particle * p=staged.Particles[ids.at(a)],* q=staged.Particles[ids.at(b)];
        bool sphere=t.get<std::string>("kind")=="sphere";
        need(sphere==(p->Verts.Size()==1 && q->Verts.Size()==1) && (sphere || t.get<std::string>("kind")=="features"),"contact kind/geometry mismatch");
        std::unique_ptr<CInteracton> c(sphere?static_cast<CInteracton *>(new CInteractonSphere(p,q,staged.ContactLaw)):new CInteracton(p,q,staged.ContactLaw));
        c->BothFree=flag(t,"BothFree");c->First=flag(t,"First");need(c->BothFree==(p->IsFree()&&q->IsFree()),"periodicity/fixity mismatch");
        c->Nc=t.get<size_t>("Nc");c->Nsc=t.get<size_t>("Nsc");c->Nr=t.get<size_t>("Nr");
#define READ_D(n) c->n=number(t,#n);
        MS_CP_CD(READ_D)
#undef READ_D
#define READ_V(n) c->n=vector(t,#n);
        MS_CP_CV(READ_V)
#undef READ_V
        if(sphere){auto * s=static_cast<CInteractonSphere *>(c.get());s->InContact=flag(t,"InContact");s->cohesion=flag(t,"cohesion");
#define READ_D(n) s->n=number(t,#n);
            MS_CP_SD(READ_D)
#undef READ_D
#define READ_V(n) s->n=vector(t,#n);
            MS_CP_SV(READ_V)
#undef READ_V
        } else {
            features(t,"EE",c->Lee,c->Fdee,p->Edges.Size(),q->Edges.Size());features(t,"VF",c->Lvf,c->Fdvf,p->Verts.Size(),q->Faces.Size());features(t,"FV",c->Lfv,c->Fdfv,p->Faces.Size(),q->Verts.Size());
            features(t,"VT",c->Lvt,c->Fdvt,p->Verts.Size(),q->Tori.Size());features(t,"TV",c->Ltv,c->Fdtv,p->Tori.Size(),q->Verts.Size());
            features(t,"VC",c->Lvc,c->Fdvc,p->Verts.Size(),q->Cylinders.Size());features(t,"CV",c->Lcv,c->Fdcv,p->Cylinders.Size(),q->Verts.Size());
        }
        staged.PairtoCInt.emplace(HashFunction(p->Index,q->Index),c.get());staged.Interactons.Push(c.get());history_owners.push_back(std::move(c));
    }
    for(auto const & row:in.get_child("candidates")){auto const & x=row.second;auto a=x.get<uint64_t>("a"),b=x.get<uint64_t>("b");need(a!=b&&ids.count(a)&&ids.count(b),"bad candidate identity");staged.ListPosPairs.Push({ids.at(a),ids.at(b)});}
    for(auto const & row:in.get_child("friction_overrides")){auto const & x=row.second;need(staged.FricCoeff.emplace(std::make_pair(x.get<int>("tag_a"),x.get<int>("tag_b")),number(x,"mu")).second,"duplicate friction override");}
#ifdef USE_CUDA
    gpu_history(in,"gpu_vv",staged.CheckpointVV,true,contact_rows.size());gpu_history(in,"gpu_ee",staged.CheckpointEE,false,contact_rows.size());
    gpu_history(in,"gpu_vf",staged.CheckpointVF,false,contact_rows.size());gpu_history(in,"gpu_fv",staged.CheckpointFV,false,contact_rows.size());
    // Validate all transient GPU references against the independently restored
    // physical contact/feature records BEFORE committing the loaded domain.
    std::vector<DynInteractonCU> evv,eee,evf,efv;
    std::vector<size_t> vo(staged.Particles.Size()),eo(vo.size()),fo(vo.size());
    size_t nv=0,ne=0,nf=0;
    for(size_t i=0;i<vo.size();i++){vo[i]=nv;eo[i]=ne;fo[i]=nf;nv+=staged.Particles[i]->Verts.Size();ne+=staged.Particles[i]->Edges.Size();nf+=staged.Particles[i]->Faces.Size();}
    for(size_t i=0;i<staged.Interactons.Size();i++) {
        auto * c=static_cast<CInteracton *>(staged.Interactons[i]);DynInteractonCU f={};f.Idx=i;
        if(dynamic_cast<CInteractonSphere *>(c)) evv.push_back(f);
        else {
            for(auto const & p:c->Lee){f.IF1=eo[c->I1]+p.first;f.IF2=eo[c->I2]+p.second;eee.push_back(f);}
            for(auto const & p:c->Lvf){f.IF1=vo[c->I1]+p.first;f.IF2=fo[c->I2]+p.second;evf.push_back(f);}
            for(auto const & p:c->Lfv){f.IF1=fo[c->I1]+p.first;f.IF2=vo[c->I2]+p.second;efv.push_back(f);}
        }
    }
    auto check=[](std::vector<DynInteractonCU> const & expected,std::vector<DynInteractonCU> const & saved,bool sphere) {
        need(expected.size()==saved.size(),"GPU history record count mismatch");
        for(size_t i=0;i<expected.size();i++) need(expected[i].Idx==saved[i].Idx && (sphere ||
            (expected[i].IF1==saved[i].IF1 && expected[i].IF2==saved[i].IF2)),"GPU physical-feature reference mismatch");
    };
    check(evv,staged.CheckpointVV,true);check(eee,staged.CheckpointEE,false);check(evf,staged.CheckpointVF,false);check(efv,staged.CheckpointFV,false);
#endif
    // Commit only after parsing/validating into RAII-owned temporary objects.
#define COPY_D(n) n=staged.n;
    MS_CP_DOMAIN_D(COPY_D)
#undef COPY_D
    Time=staged.Time;iter=staged.iter;NextPersistentId=staged.NextPersistentId;ContactLaw=staged.ContactLaw;SphereCoulombMode=staged.SphereCoulombMode;
    RotPar=staged.RotPar;SphereTensileCutoff=staged.SphereTensileCutoff;SphereFirstContactCorrection=staged.SphereFirstContactCorrection;Per=staged.Per;
    Particles=staged.Particles;Interactons=staged.Interactons;PairtoCInt=staged.PairtoCInt;ListPosPairs=staged.ListPosPairs;FricCoeff=staged.FricCoeff;A=staged.A;Wdot=staged.Wdot;
    for(auto & p:owners) p.release();for(auto & c:history_owners) c.release();
#ifdef USE_CUDA
    CheckpointVV=std::move(staged.CheckpointVV);CheckpointEE=std::move(staged.CheckpointEE);CheckpointVF=std::move(staged.CheckpointVF);CheckpointFV=std::move(staged.CheckpointFV);
    pForceVV=ContactLaw==0?CalcForceVV:CalcForceVV_Hertz;
#endif
    UseVelocityVerlet=true;Initialized=true;Finished=false;CheckpointReady=true;VVStateReady=true;
}

#ifdef USE_CUDA
inline void Domain::RestoreCheckpointDevice() {
    using namespace Checkpoint;
    restore_gpu_history(bDynInteractonsVV,CheckpointVV,true);restore_gpu_history(bDynInteractonsEE,CheckpointEE,false);
    restore_gpu_history(bDynInteractonsVF,CheckpointVF,false);restore_gpu_history(bDynInteractonsFV,CheckpointFV,false);
    thrust::host_vector<real3> refs(bVertsoCU.size());
    for(auto * p:Particles) for(size_t i=0;i<p->Vertso.Size();i++){auto const & v=*p->Vertso[i];refs[p->Nvi+i]=make_real3(v(0),v(1),v(2));}bVertsoCU=refs;
    thrust::host_vector<ComInteractonCU> common=bComInteractons;
    for(size_t i=0;i<common.size();i++){auto * c=static_cast<CInteracton *>(Interactons[i]);
        auto pack=[](Vec3_t const & v){return make_real3(v(0),v(1),v(2));};
        common[i].Fnnet=pack(c->Fnet);common[i].Ftnet=pack(c->Ftnet);common[i].Fndpot=pack(c->Fndpot);common[i].Ftdpot=pack(c->Ftdpot);common[i].Fther=pack(c->Fther);
        common[i].Nc=c->Nc;common[i].NormalOverlap=c->NormalOverlap;
        common[i].Epot=c->Epot;common[i].dEvis=c->dEvis;common[i].dEfric=c->dEfric;
    }bComInteractons=common;
    CheckpointVV.clear();CheckpointEE.clear();CheckpointVF.clear();CheckpointFV.clear();
}
#endif
#undef MS_CP_DOMAIN_D
#undef MS_CP_PROPS
#undef MS_CP_PV
#undef MS_CP_PB
#undef MS_CP_CD
#undef MS_CP_CV
#undef MS_CP_SD
#undef MS_CP_SV
} // DEM
#endif
