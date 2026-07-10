/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 ************************************************************************/

#ifndef MECHSYS_MPM_FLUID_FILTER_H
#define MECHSYS_MPM_FLUID_FILTER_H

#include <algorithm>
#include <cmath>

namespace MPM
{

// vNRL artificial viscosity from Telikicherla and Moutsanidis (2025),
// Eq. (8): added pressure only during local compression, tr(D) < 0.
inline double VNRLPressure(double rho, double divv, double dx,
                                     double c0, double cL, double cQ)
{
    if (cL <= 0.0 && cQ <= 0.0) return 0.0;

    const double comp = std::max(-divv,0.0);
    double q = 0.0;
    if (cL > 0.0) q += cL*dx*c0*comp;
    if (cQ > 0.0) q += cQ*dx*dx*comp*comp;
    return rho*q;
}

#ifdef USE_CUDA
__host__ __device__ inline real VNRLPressureCU(real rho, real divv, real dx,
                                                         real c0, real cL, real cQ)
{
    if (cL <= 0.0 && cQ <= 0.0) return 0.0;

    real comp = fmax(-divv,(real)0.0);
    real q = 0.0;
    if (cL > 0.0) q += cL*dx*c0*comp;
    if (cQ > 0.0) q += cQ*dx*dx*comp*comp;
    return rho*q;
}
#endif

}

#endif // MECHSYS_MPM_FLUID_FILTER_H
