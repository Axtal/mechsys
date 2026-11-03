/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 ************************************************************************/

#ifndef MECHSYS_MPM_EOS_H
#define MECHSYS_MPM_EOS_H

#include <cmath>

namespace MPM
{

enum FluidEOS
{
    EOS_TAIT   = 0,
    EOS_LINEAR = 1
};

// Weakly-compressible Tait/Cole EOS, usually gamma=7 for water.
inline double TaitEOS(double rho, double rho0, double c0, double gamma)
{
    return rho0*c0*c0/gamma*(std::pow(rho/rho0,gamma) - 1.0);
}

// Linear acoustic EOS: first-order form of the pressure-density relation.
inline double LinearEOS(double rho, double rho0, double c0)
{
    return c0*c0*(rho-rho0);
}

inline double EOS(int type, double rho, double rho0, double c0, double gamma)
{
    return (type == EOS_LINEAR) ? LinearEOS(rho,rho0,c0)
                                : TaitEOS(rho,rho0,c0,gamma);
}

inline const char * EOSName(int type)
{
    return (type == EOS_LINEAR) ? "Linear EOS" : "Tait-Cole EOS";
}

#ifdef USE_CUDA
__host__ __device__ inline real TaitEOSCU(real rho, real rho0, real c0, real gamma)
{
    return rho0*c0*c0/gamma*(pow(rho/rho0,gamma) - 1.0);
}

__host__ __device__ inline real LinearEOSCU(real rho, real rho0, real c0)
{
    return c0*c0*(rho-rho0);
}

__host__ __device__ inline real EOSCU(int type, real rho, real rho0, real c0, real gamma)
{
    return (type == EOS_LINEAR) ? LinearEOSCU(rho,rho0,c0)
                                : TaitEOSCU(rho,rho0,c0,gamma);
}
#endif

}

#endif // MECHSYS_MPM_EOS_H
