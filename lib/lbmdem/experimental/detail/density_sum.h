// GPL-3.0-or-later. FP64 Neumaier sum; no population modification.
#ifndef MECHSYS_EXPERIMENTAL_DENSITY_SUM_H
#define MECHSYS_EXPERIMENTAL_DENSITY_SUM_H
#include <cmath>
#ifdef __CUDACC__
#define MS_DENSITY_HD __host__ __device__
#else
#define MS_DENSITY_HD
#endif
namespace MechSysDensitySum {
MS_DENSITY_HD inline double add(double a,double b){
#ifdef __CUDA_ARCH__
    return __dadd_rn(a,b);
#else
    volatile double r=a+b;return r;
#endif
}
// Positive-population, conservative normal-arithmetic envelope. Small residuals
// are multiples of at least 2^-502; partials are below 2^455. No subnormal
// arithmetic/overflow is needed for this graph inside the supported envelope.
MS_DENSITY_HD inline bool supported(double x){return x>=0x1p-450&&x<=0x1p450;}
MS_DENSITY_HD inline double sum15(double const *f,bool &valid){
    double s=0,c=0;
    for(int k=0;k<15;k++){
        double x=f[k];valid=valid&&supported(x);double t=add(s,x);
        double e=fabs(s)>=fabs(x)?add(add(s,-t),x):add(add(x,-t),s);
        c=add(c,e);s=t;
    }
    return add(s,c);
}
}
#undef MS_DENSITY_HD
#endif
