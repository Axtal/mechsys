#define main periodic_pid_example_main
#include "../../tdem/TX_PID_periodic.cpp"
#undef main
#include <stdexcept>
int main() {
    DEM::Domain d;
    d.Xmin=d.Ymin=d.Zmin=-5; d.Xmax=d.Ymax=d.Zmax=5;
    d.Per=Vec3_t(10,10,10);
    d.AddSphere(-1,Vec3_t(0,0,0),1,1);
    d.AddSphere(-2,Vec3_t(1.9,.2,.3),1,1);
    auto *c=new DEM::CInteractonSphere(d.Particles[0],d.Particles[1]);
    d.PairtoCInt[HashFunction(0,1)]=c;
    for (int elastic=0;elastic<2;++elastic) {
        c->Fnet=Vec3_t(10*elastic,0,0); c->Ftnet=Vec3_t(0,elastic,0);
        c->Fndpot=Vec3_t(2,0,0); c->Ftdpot=Vec3_t(0,3,0); c->Fther=Vec3_t(0,0,4);
        Mat3_t stress=ContactStress(d);
        double force[3]={double(10*elastic+2),double(elastic+3),4};
        double branch[3]={-1.9,-.2,-.3};
        for(int i=0;i<3;++i) for(int j=0;j<3;++j)
            if(fabs(stress(i,j)-force[i]*branch[j]/1000.)>1e-13)
                throw std::runtime_error("PID virial must include every applied force component");
    }
    printf("PERIODIC_PID_STRESS_PASS cases=2\n");
}
