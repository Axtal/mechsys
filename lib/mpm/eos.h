/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 ************************************************************************/

#ifndef MECHSYS_MPM_EOS_H
#define MECHSYS_MPM_EOS_H

#include <cmath>

namespace MPM
{

// Weakly-compressible Tait/Cole EOS for free-surface water:
// p = rho0*c0^2/gamma*((rho/rho0)^gamma - 1), usually gamma=7.
inline double TaitPressure(double rho, double rho0, double c0, double gamma)
{
    const double tol = 1.0e-12;
    if (rho <= tol || rho0 <= tol || c0 <= tol || gamma <= tol) return 0.0;
    return rho0*c0*c0/gamma*(std::pow(rho/rho0,gamma) - 1.0);
}

#ifdef USE_CUDA
__host__ __device__ inline real TaitPressureCU(real rho, real rho0, real c0, real gamma)
{
    const real tol = 1.0e-12;
    if (rho <= tol || rho0 <= tol || c0 <= tol || gamma <= tol) return 0.0;
    return rho0*c0*c0/gamma*(pow(rho/rho0,gamma) - 1.0);
}
#endif

}

#endif // MECHSYS_MPM_EOS_H
