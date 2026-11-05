#include <mechsys/dem/domain.h>
#include <stdexcept>
#include <type_traits>

struct OldInteraction : DEM::Interacton {
    int calls=0;
    bool fail=false;
    bool UpdateContacts(double,Vec3_t const &,size_t) override { return true; }
    void UpdateParameters(size_t) override {}
    bool CalcForce(double,Vec3_t const &,size_t,size_t) override {
        ++calls;
        if (fail) throw std::runtime_error("custom contact failure");
        return sphereOptions.mode==2 && sphereOptions.cutoff && sphereOptions.onset;
    }
};
struct OldSphere : DEM::CInteractonSphere {
    int calls=0;
    OldSphere(DEM::Particle *a,DEM::Particle *b):CInteractonSphere(a,b) {}
    bool CalcForce(double,Vec3_t const &,size_t,size_t) override { ++calls; return true; }
};
static_assert(!std::is_abstract<OldInteraction>::value,"old interaction interface must remain concrete");
int main() {
    OldInteraction custom;
    DEM::Interacton *p=&custom;
    if (!p->CalcForce(.01,OrthoSys::O,0,0,2,true,true) || custom.calls!=1)
        throw std::runtime_error("sphere options bypassed old virtual override");
    if (p->CalcForce(.01,OrthoSys::O,0,0))
        throw std::runtime_error("sphere options leaked into legacy call");
    custom.fail=true;
    try { p->CalcForce(.01,OrthoSys::O,0,0,2,true,true); }
    catch (std::runtime_error const &) { custom.fail=false; }
    if (p->CalcForce(.01,OrthoSys::O,0,0))
        throw std::runtime_error("sphere options leaked after exception");
    DEM::Domain d;
    d.AddSphere(-1,Vec3_t(0,0,0),1,1);
    d.AddSphere(-2,Vec3_t(1.9,0,0),1,1);
    OldSphere sphere(d.Particles[0],d.Particles[1]);
    p=&sphere;
    if (!p->CalcForce(.01,OrthoSys::O,0,0,2,true,true) || sphere.calls!=1)
        throw std::runtime_error("old concrete sphere override was lost");
    printf("FORCE_API_COMPATIBILITY_PASS\n");
}
