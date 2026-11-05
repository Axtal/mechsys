#include <mechsys/dem/domain.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include "feature_oracle.h"

static int scenario=0;
static void vec(char const * name,Vec3_t const & v) {
    printf(",\"%s\":[%.17g,%.17g,%.17g]",name,v(0),v(1),v(2));
}
static void setup(DEM::Domain & d,void *) {
    if(scenario==15 && (d.iter==200 || d.iter==400)) {
#ifdef USE_CUDA
        d.DnLoadDevice(1,true);
#endif
        Vec3_t v=d.Particles[1]->v;v(2)=d.iter==200?0.5:-0.5;
#ifdef USE_CUDA
        d.SetVelocityDevice(1,v(0),v(1),v(2));
#endif
        d.Particles[1]->v=v;
    }
    if(scenario==1 && d.iter==75) {
#ifdef USE_CUDA
        // Setup precedes the solver's output download: read the CURRENT device
        // velocity before imposing a matched reversal, not the prior host report.
        d.DnLoadDevice(1,true);
#endif
        for(size_t i=0;i<2;i++) {
            Vec3_t v=d.Particles[i]->v;v(1)=-v(1);
#ifdef USE_CUDA
            d.SetVelocityDevice(i,v(0),v(1),v(2));
#endif
            d.Particles[i]->v=v;
        }
    }
}
static void report(DEM::Domain & d,void *) {
    printf("V2_STATE {\"step\":%zu,\"time\":%.17g,\"particles\":[",d.iter,d.Time);
    double kinetic=0,potential=0;
    Vec3_t balance(0,0,0),momentum(0,0,0);
    for(size_t i=0;i<d.Particles.Size();i++) {
        if(i) printf(",");auto * p=d.Particles[i];
        if(p->Tag>=0)throw std::runtime_error("V2 fixture exceeded baseline overlap guard; mutated diagnostic tag is not accepted physics");
        printf("{\"id\":%llu,\"tag\":%d",(unsigned long long)p->PersistentId,p->Tag);
        vec("x",p->x);vec("v",p->v);vec("w",p->w);vec("F",p->F);vec("Ff",p->Ff);
        Vec3_t torque=p->T;
#ifdef USE_CUDA
        if(!d.CheckpointReady) {DEM::ParticleCU cu=d.bParticlesCU[i];torque=Vec3_t(cu.T.x,cu.T.y,cu.T.z);}
#endif
        vec("T",torque);
        printf(",\"Q\":[%.17g,%.17g,%.17g,%.17g]}",p->Q(0),p->Q(1),p->Q(2),p->Q(3));
        kinetic+=0.5*p->Props.m*dot(p->v,p->v)+0.5*(p->I(0)*p->w(0)*p->w(0)+p->I(1)*p->w(1)*p->w(1)+p->I(2)*p->w(2)*p->w(2));
        momentum+=p->Props.m*p->v;balance+=p->F-p->Ff;
    }
    printf("],\"contacts\":[");
    for(size_t i=0;i<d.Interactons.Size();i++) {
        if(i)printf(",");auto * c=static_cast<DEM::CInteracton *>(d.Interactons[i]);
        printf("{\"a\":%llu,\"b\":%llu",(unsigned long long)c->P1->PersistentId,(unsigned long long)c->P2->PersistentId);
        vec("normal",c->Fnet);vec("tangent",c->Ftnet);vec("normal_damping",c->Fndpot);vec("tangent_damping",c->Ftdpot);
        printf(",\"Epot\":%.17g,\"dEvis\":%.17g,\"dEfric\":%.17g",c->Epot,c->dEvis,c->dEfric);
        printf(",\"Nc\":%zu,\"max_overlap\":%.17g",c->Nc,c->NormalOverlap);
        if(auto * s=dynamic_cast<DEM::CInteractonSphere *>(c)) {
            vec("displacement",s->Fdvv);vec("rolling",s->Fdr);
            printf(",\"contact\":%d,\"overlap\":%.17g",int(s->InContact),s->PrevDelta);
        } else {
            printf(",\"features\":[");bool first=true;
            DEM::FrictionMap_t * maps[]={&c->Fdee,&c->Fdvf,&c->Fdfv};
            for(int family=0;family<3;family++) {
                std::vector<size_t> keys;for(auto const & f:*maps[family]) keys.push_back(f.first);std::sort(keys.begin(),keys.end());
                for(auto key:keys) {if(!first)printf(",");first=false;printf("{\"family\":%d,\"key\":%zu",family,key);vec("displacement",maps[family]->at(key));printf("}");}
            }
            printf("]");
        }
        printf("}");potential+=c->Epot;
    }
    printf("],\"kinetic\":%.17g,\"potential\":%.17g",kinetic,potential);
    vec("force_balance",balance);vec("momentum",momentum);printf("}\n");
}
static void create(DEM::Domain & d,int law,bool damping) {
    d.ContactLaw=law;d.UseVelocityVerlet=true;d.Alpha=0.001;
    if(scenario>=11) {
        if(scenario==11 || scenario==12 || scenario==16) {
            if(scenario!=12)d.AddSphere(-1,Vec3_t(0.074,0,0),0.02,1000);
            d.AddCube(-1,Vec3_t(0,0,0),0.005,0.1,1000,0,&OrthoSys::e0);
            if(scenario==12)d.AddSphere(-1,Vec3_t(0.074,0,0),0.02,1000);
        } else {
            d.AddRecBox(-1,Vec3_t(0,0,0),Vec3_t(0.1,0.01,0.01),0.005,1000,0,&OrthoSys::e0);
            d.AddRecBox(-1,Vec3_t(0,0,scenario==15?0.0209:0.019),Vec3_t(0.01,0.1,0.01),0.005,1000,0,&OrthoSys::e0);
        }
        d.Alpha=0.0001;
    } else if(scenario==6) {
        d.AddSphere(-1,Vec3_t(0.109,0,0),0.1,1);
        d.AddPlane(-2,Vec3_t(0,0,0),0.01,1,1,1,M_PI/2,&OrthoSys::e1);
        d.Particles[1]->FixVeloc(0,0.05,0);
    } else {
        d.AddSphere(-1,Vec3_t(0,0,0),0.1,1);
        d.AddSphere(-1,Vec3_t(scenario>=9?0.202:0.199,0,0),0.1,1); // Repeated tags exercise independent IDs.
        if(scenario==7) d.AddSphere(-1,Vec3_t(0.398,0,0),0.1,1);
    }
    double kn=1000,kt=law?400:20;
    double nu=0.5*kn/kt-1;
    double effective_kn=law?(4.0/3.0)*sqrt(0.05)*(kn/(2*(1-nu*nu))):kn;
    double holding=0.75*effective_kn*(law?pow(0.001,1.5):0.001);
    for(size_t i=0;i<d.Particles.Size();i++) {
        auto * p=d.Particles[i];p->Props.Kn=kn;p->Props.Kt=kt;p->Props.Gn=damping?3:0;p->Props.Gt=damping?2:0;
        p->Props.Mu=0.4;p->Props.Beta=0.3;p->Props.Eta=0.2;
        if(scenario!=6) {
            p->v=Vec3_t(10,i?-0.1:0.1,0);p->w=Vec3_t(0,0,i?-1:1);
            p->Ff=Vec3_t(i?-holding:holding,0,0);
        }
        if(scenario==2) {p->v(1)=i?-5:5;p->Ff=OrthoSys::O;}
        if(scenario==4) {p->v(0)=10+(i?0.3:-0.3);p->Ff=Vec3_t(i?-0.2:0.2,0,0);}
        if(scenario==9 || scenario==10) {p->v=Vec3_t(i?-0.3:0.3,scenario==10?(i?-0.1:0.1):0,0);p->w=OrthoSys::O;p->Ff=OrthoSys::O;}
        if(scenario>=11) {
            p->FixVeloc(0,i?0.01:0,0);p->FixFree=true;p->Ff=OrthoSys::O;
            p->w=Vec3_t(0,0,scenario==14?(i?-10:10):0);
            if(scenario==15 && i==1)p->v(2)=-0.5;
            if(scenario==16)p->v(0)=10;
        }
    }
    if(scenario==5) {
        d.Xmin=-0.5;d.Xmax=0.5;d.Particles[0]->Translate(Vec3_t(0.49,0,0));d.Particles[1]->Translate(Vec3_t(-0.51,0,0));
    }
    if(scenario==6) {d.Particles[0]->v=Vec3_t(0,0.1,0);d.Particles[0]->Ff=Vec3_t(-0.375,0,0);}
    if(scenario==7) {d.Particles[0]->Ff=Vec3_t(holding,0,0);d.Particles[1]->Ff=OrthoSys::O;d.Particles[2]->Ff=Vec3_t(-holding,0,0);}
    if(scenario==16){d.Xmin=-0.5;d.Xmax=0.5;for(auto * p:d.Particles)p->Translate(Vec3_t(0.49,0,0));}
#ifdef USE_CUDA
    d.pForceVV=law?DEM::CalcForceVV_Hertz:DEM::CalcForceVV;
#endif
}
int main(int argc,char ** argv) try {
    if(argc<7) throw std::runtime_error("usage: probe mode scenario law damping steps checkpoint_file");
    char const * mode=argv[1];scenario=atoi(argv[2]);int law=atoi(argv[3]);bool damping=atoi(argv[4]);int steps=atoi(argv[5]);
    DEM::Domain d;
    if(!strcmp(mode,"continue") || !strcmp(mode,"validate")) {
        d.LoadCheckpoint(argv[6]);report(d,NULL);
        if(!strcmp(mode,"validate"))return 0;
#ifdef USE_CUDA
#endif
    } else create(d,law,damping);
    if(scenario>=11) {
        if(law)throw std::runtime_error("V2 bounded feature fixtures are registered for the linear feature law only");
        if(std::getenv("MECHSYS_V2_FEATURE_ORACLE"))d.ForceTimeCheck=V2Geometry::check;
#ifdef USE_CUDA
#endif
        V2Geometry::drop_feature=!strcmp(mode,"negative");
    }
    printf("V2_MODE vv=%d scenario=%d law=%zu damping=%d\n",int(d.UseVelocityVerlet),scenario,d.ContactLaw,int(damping));
    d.Solve((steps-0.5)*2e-5,2e-5,2e-5,setup,report,NULL,false,1);
    if(!strcmp(mode,"split")) d.SaveCheckpoint(argv[6]);
    if(scenario>=11)V2Geometry::summary();
    return 0;
}
MECHSYS_CATCH
