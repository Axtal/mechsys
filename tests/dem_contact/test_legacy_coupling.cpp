// GPL-3.0-or-later. Ordinary coupled Solve compatibility, two tiny steps.
#include <mechsys/lbmdem/Domain.h>
#include <stdexcept>

static void need(bool ok,char const *message)
{
    if (!ok) throw std::runtime_error(message);
}
int main() try
{
    LBMDEM::Domain d(D3Q15,.1,iVec3_t(12,12,12),1,1);
    d.Alpha=.1;d.DEMDOM.UseVelocityVerlet=false;
    d.PeriodicX=d.PeriodicY=d.PeriodicZ=true;
    d.DEMDOM.AddSphere(-1,Vec3_t(6.13,6.19,6.23),2,1);
    d.LBMDOM.Sc=0;d.LBMDOM.G[0]=d.LBMDOM.Gmix=0;
    Vec3_t rest(0,0,0);
    for(int z=0;z<12;z++)for(int y=0;y<12;y++)for(int x=0;x<12;x++)
        d.LBMDOM.Initialize(0,iVec3_t(x,y,z),1,rest);
    need(d.DEMDOM.A.Size()==0 && d.DEMDOM.Wdot.Size()==0,
         "Fixture must exercise coupled entry without standalone VV bootstrap");
    d.Solve(2,1,NULL,NULL,NULL,false,1);
#ifdef USE_CUDA
    d.DnLoadDevice(1);
    need(d.DEMDOM.bA.size()==1 && d.DEMDOM.bWdot.size()==1,
         "Shared legacy upload must initialize acceleration storage");
    need(cudaDeviceSynchronize()==cudaSuccess,"legacy coupled GPU failure");
#endif
    need(d.iter==2 && d.Time==2,"legacy coupled clocks");
    auto *p=d.DEMDOM.Particles[0];
    for(int j=0;j<3;j++)need(std::isfinite(p->x(j))&&std::isfinite(p->v(j)),"finite legacy particle state");
    need(norm(p->x-Vec3_t(6.13,6.19,6.23))<1e-10,"resting legacy sphere moved");
    for(int z=0;z<12;z++)for(int y=0;y<12;y++)for(int x=0;x<12;x++)
        need(std::isfinite(d.LBMDOM.Rho[0][x][y][z])&&d.LBMDOM.Rho[0][x][y][z]>0,"positive legacy density");
    std::cout<<"LEGACY_COUPLING_PASS steps=2"<<std::endl;
    return 0;
} MECHSYS_CATCH
