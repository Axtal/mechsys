// Independent bounded convex-polygon/segment reference; no DEM Distance/Overlap.
#include <map>
#include <set>
#include <tuple>
#include <stdexcept>
namespace V2Geometry {
using Key=std::tuple<int,uint64_t,uint64_t,size_t,size_t>;
static size_t checks=0,required[3]={0,0,0},missing=0,duplicates=0;
static bool drop_feature=false;
static double clamp(double x){return std::max(0.0,std::min(1.0,x));}
static double point_segment(Vec3_t const & x,Vec3_t const & a,Vec3_t const & b) {
    Vec3_t u=b-a;return norm(x-(a+clamp(dot(x-a,u)/dot(u,u))*u));
}
static double segment_segment(Vec3_t const & a,Vec3_t const & b,Vec3_t const & c,Vec3_t const & d) {
    Vec3_t u=b-a,v=d-c,w=a-c;double aa=dot(u,u),bb=dot(u,v),cc=dot(v,v),dd=dot(u,w),ee=dot(v,w);
    double best=std::min(std::min(point_segment(a,c,d),point_segment(b,c,d)),std::min(point_segment(c,a,b),point_segment(d,a,b)));
    double determinant=aa*cc-bb*bb;
    if(determinant>1e-24) {double s=(bb*ee-cc*dd)/determinant,t=(aa*ee-bb*dd)/determinant;
        if(s>=0&&s<=1&&t>=0&&t<=1)best=std::min(best,norm(w+s*u-t*v));}
    return best;
}
static double vertex_face(Vec3_t const & x,std::vector<Vec3_t> const & face) {
    Vec3_t n=cross(Vec3_t(face[1]-face[0]),Vec3_t(face[2]-face[0]));n/=norm(n);
    Vec3_t projected=x-dot(x-face[0],n)*n;double best=1e300;bool pos=false,neg=false;
    for(size_t i=0;i<face.size();i++) {auto const & a=face[i];auto const & b=face[(i+1)%face.size()];
        double side=dot(cross(Vec3_t(b-a),Vec3_t(projected-a)),n);pos|=side>1e-14;neg|=side< -1e-14;
        best=std::min(best,point_segment(x,a,b));}
    if(!(pos&&neg))best=std::min(best,fabs(dot(x-face[0],n)));return best;
}
static Key key(int family,DEM::Particle * a,DEM::Particle * b,size_t i,size_t j) {
    uint64_t x=a->PersistentId,y=b->PersistentId;
    if(x>y){std::swap(x,y);std::swap(i,j);if(family)family=3-family;}
    return {family,x,y,i,j};
}
static void check(DEM::Domain & d,void *) {
#ifdef USE_CUDA
    d.DnLoadDevice(1,true);
#endif
    std::set<Key> reference;std::map<Key,size_t> uploaded;
    for(size_t i=0;i<d.Particles.Size();i++)for(size_t j=i+1;j<d.Particles.Size();j++) {
        auto * a=d.Particles[i];auto * b=d.Particles[j];
        if(!a->IsFree()&&!b->IsFree())continue;
        Vec3_t shift(0,0,0),branch=b->x-a->x;
        if(a->IsFree()&&b->IsFree())for(size_t k=0;k<3;k++)if(d.Per(k)>0)shift(k)=-std::round(branch(k)/d.Per(k))*d.Per(k);
        auto av=[&](size_t k){return *a->Verts[k];};auto bv=[&](size_t k){return Vec3_t(*b->Verts[k]+shift);};
        double cutoff=a->Props.R+b->Props.R;
        for(size_t p=0;p<a->EdgeCon.Size();p++)for(size_t q=0;q<b->EdgeCon.Size();q++) {
            auto const & e=a->EdgeCon[p];auto const & f=b->EdgeCon[q];
            if(segment_segment(av(e[0]),av(e[1]),bv(f[0]),bv(f[1]))<cutoff)reference.insert(key(0,a,b,p,q));
        }
        for(size_t p=0;p<a->Verts.Size();p++)for(size_t q=0;q<b->FaceCon.Size();q++) {
            std::vector<Vec3_t> face;for(int k:b->FaceCon[q])face.push_back(bv(k));
            if(vertex_face(av(p),face)<cutoff)reference.insert(key(1,a,b,p,q));
        }
        for(size_t p=0;p<a->FaceCon.Size();p++)for(size_t q=0;q<b->Verts.Size();q++) {
            std::vector<Vec3_t> face;for(int k:a->FaceCon[p])face.push_back(av(k));
            if(vertex_face(bv(q),face)<cutoff)reference.insert(key(2,a,b,p,q));
        }
    }
#ifdef USE_CUDA
    auto inspect=[&](thrust::device_vector<DEM::DynInteractonCU> & buffer,int family) {
        thrust::host_vector<DEM::DynInteractonCU> h=buffer;
        for(size_t i=0;i<h.size();) {
            auto f=h[i];if(f.Idx>=d.Interactons.Size())throw std::runtime_error("V2 feature oracle invalid common identity");
            auto * c=static_cast<DEM::CInteracton *>(d.Interactons[f.Idx]);auto * a=c->P1;auto * b=c->P2;
            size_t x=f.IF1-(family==0?a->Nei:family==1?a->Nvi:a->Nfi),y=f.IF2-(family==0?b->Nei:family==1?b->Nfi:b->Nvi);
            Key k=key(family,a,b,x,y);
            if(drop_feature && reference.count(k)) {
                // Negative control deletes an ACTUAL uploaded record and extent.
                h.erase(h.begin()+i);buffer=h;
                if(family==0){d.demaux.neeint--;d.pDynInteractonsEE=thrust::raw_pointer_cast(buffer.data());}
                if(family==1){d.demaux.nvfint--;d.pDynInteractonsVF=thrust::raw_pointer_cast(buffer.data());}
                if(family==2){d.demaux.nfvint--;d.pDynInteractonsFV=thrust::raw_pointer_cast(buffer.data());}
                cudaMemcpy(d.pdemaux,&d.demaux,sizeof(d.demaux),cudaMemcpyHostToDevice);drop_feature=false;
                continue;
            }
            uploaded[k]++;
            i++;
        }
    };
    inspect(d.bDynInteractonsEE,0);inspect(d.bDynInteractonsVF,1);inspect(d.bDynInteractonsFV,2);
#else
    for(auto * interaction:d.Interactons) {
        auto * c=static_cast<DEM::CInteracton *>(interaction);
        DEM::ListContacts_t * lists[]={&c->Lee,&c->Lvf,&c->Lfv};
        for(int family=0;family<3;family++)for(auto const & f:*lists[family])uploaded[key(family,c->P1,c->P2,f.first,f.second)]++;
    }
#endif
    size_t miss=0,dup=0;
    for(auto const & k:reference){required[std::get<0>(k)]++;if(!uploaded.count(k))miss++;}
    for(auto const & k:uploaded)if(k.second!=1)dup+=k.second-1;
    checks++;missing+=miss;duplicates+=dup;
    printf("V2_FEATURE_ORACLE step=%zu reference=%zu missing=%zu duplicates=%zu\n",d.iter,reference.size(),miss,dup);
    if(miss||dup)throw std::runtime_error("V2 independent feature-completeness oracle FAIL");
}
static void summary(){printf("V2_FEATURE_SUMMARY checks=%zu EE=%zu VF=%zu FV=%zu missing=%zu duplicates=%zu\n",checks,required[0],required[1],required[2],missing,duplicates);}
}
