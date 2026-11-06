#define main periodic_pid_example_main
#include "../../tdem/TX_PID_periodic.cu"
#undef main
#include <stdexcept>
template<class T> static T *ptr(thrust::device_vector<T> &v) { return thrust::raw_pointer_cast(v.data()); }
int main() {
    for(int elastic=0;elastic<2;++elastic) {
        DEM::InteractonCU i={};i.BothFree=true;
        DEM::ComInteractonCU c={};c.I1=0;c.I2=1;
        c.Fnnet=make_real3(10*elastic,0,0);c.Ftnet=make_real3(0,elastic,0);
        c.Fndpot=make_real3(2,0,0);c.Ftdpot=make_real3(0,3,0);c.Fther=make_real3(0,0,4);
        DEM::DynInteractonCU h={};DEM::dem_aux a={};a.nvvint=1;a.Per=make_real3(10,10,10);
        DEM::DynParticleCU q[2]={};q[1].x=make_real3(1.9,.2,.3);
        thrust::device_vector<DEM::InteractonCU> di(1,i);
        thrust::device_vector<DEM::ComInteractonCU> dc(1,c);
        thrust::device_vector<DEM::DynInteractonCU> dh(1,h);
        thrust::device_vector<DEM::DynParticleCU> dq(q,q+2);
        thrust::device_vector<DEM::dem_aux> da(1,a);thrust::device_vector<real> ds(6,0);
        DEM::CalcStressVV<<<1,32>>>(ptr(di),ptr(dc),ptr(dh),ptr(dq),ptr(da),ptr(ds));
        if(cudaDeviceSynchronize()!=cudaSuccess) throw std::runtime_error("PID stress kernel failed");
        thrust::host_vector<real> s=ds;
        double expected[6]={-1.9*(10*elastic+2),-.2*(elastic+3),-1.2,
                            -.2*(10*elastic+2),-.3*(10*elastic+2),-.3*(elastic+3)};
        for(int j=0;j<6;++j) if(fabs(s[j]-expected[j])>1e-12)
            throw std::runtime_error("GPU PID virial omitted applied force components");
    }
    printf("PERIODIC_PID_STRESS_PASS cases=2\n");
}
