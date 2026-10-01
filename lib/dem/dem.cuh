/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2023 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

/////////////////////////////DEM CUDA implementation////////////////////

#ifndef MECHSYS_DEM_CUH
#define MECHSYS_DEM_CUH



//Mechsys
#include <mechsys/dem/interacton.h>
namespace DEM
{

struct dem_aux
{
    size_t ncoint = 0;   ///< number of common interactons
    size_t nverts = 0;   ///< number of total vertices
    size_t nfacid = 0;   ///< number of vertices list for faces
    size_t nfaces = 0;   ///< number of total faces
    size_t nedges = 0;   ///< number of total edges
    size_t nfvint = 0;   ///< number of total face vertex interactions
    size_t nvfint = 0;   ///< number of total vertex face interactions
    size_t neeint = 0;   ///< number of total edge edge interactions
    size_t nvvint = 0;   ///< number of total vertex vertex interactions
    size_t nparts = 0;   ///< number of particles
    real   dt;           ///< Time step
    real   Time;         ///< Time clock
    size_t iter   = 0;   ///< Iteration clock;
    real3  Per;          ///< Vector with the periodic boundary condition information
    real   Xmin   = 0.0;   
    real   Ymin   = 0.0;
    real   Zmin   = 0.0;
    real   Xmax   = 0.0;
    real   Ymax   = 0.0;
    real   Zmax   = 0.0;
    bool   px     = false;
    bool   py     = false;
    bool   pz     = false;
    real   strain = 0.0;   ///< Accumulated shear strain (gamma = gamma_dot * t)
    real   shear_rate = 0.0; ///< Shear rate gamma_dot (1/T)
    bool   use_wrapZ = false; ///< If true, pWrapZ handles LE shift; ResetMaxD skips it to avoid double-shifting
    real   vZmin = 0.0;    ///< Velocity of Zmin boundary (for contact-wrapping compression)
    real   vZmax = 0.0;    ///< Velocity of Zmax boundary (for contact-wrapping compression)
    size_t sphereCoulombMode = 0;       ///< Sphere Coulomb mode: 0 elastic cap with dashpot, 1 elastic cap dashpot off sliding, 2 total force cap
    bool   sphereTensileCutoff = true;  ///< Clamp tensile total normal force for sphere contacts
    bool   sphereFirstContactCorrection = false; ///< Scale first active tangential increment by contact activation fraction

};

typedef void (*ForceVV_ptr_t)(InteractonCU const *, ComInteractonCU *, DynInteractonCU *, ParticleCU *
        , DynParticleCU *, dem_aux const *, void *);

typedef void (*ForceEE_ptr_t)(size_t const *, real3 const *, InteractonCU const *, ComInteractonCU *, DynInteractonCU *, ParticleCU *
        , DynParticleCU *, dem_aux const *, void *);

typedef void (*ForceVF_ptr_t)(size_t const *, size_t const *, real3 const *, InteractonCU const *, ComInteractonCU *, DynInteractonCU *, ParticleCU *
        , DynParticleCU *, dem_aux const *, void *);

typedef void (*ForceFV_ptr_t)(size_t const *, size_t const *, real3 const *, InteractonCU const *, ComInteractonCU *, DynInteractonCU *, ParticleCU *
        , DynParticleCU *, dem_aux const *, void *);

typedef void (*Translate_ptr_t)(real3 *, ParticleCU const *, DynParticleCU *, dem_aux const *, void *);

typedef void (*Rotate_ptr_t)   (real3 *, ParticleCU const *, DynParticleCU *, dem_aux const *, void *);

typedef void (*Reset_ptr_t)    (ParticleCU *, DynParticleCU *, InteractonCU const * , ComInteractonCU *, dem_aux const *, void *);
typedef void (*VerletStep1_ptr_t)(real3 *, ParticleCU const *, DynParticleCU *, const real3  * , dem_aux const *);
typedef void (*FinalizeVelocity_ptr_t)(ParticleCU const *, DynParticleCU *, real3 * A, dem_aux const *);
typedef void (*OrientationUpdate_ptr_t)(real3 *, ParticleCU const *, DynParticleCU *,
                                        real3 const *, dem_aux const *);
typedef void (*FinalizeRotation_ptr_t)(ParticleCU const *, DynParticleCU *, real3 *,
                                       dem_aux const *);
typedef void (*ZeroFixedTorque_ptr_t)(ParticleCU *, dem_aux const *);

typedef void (*EnforceAngularFixity_ptr_t)(ParticleCU const *,
                                           DynParticleCU *,
                                           dem_aux const *);
typedef void (*RecomputeAccelerations_ptr_t)(ParticleCU const *, DynParticleCU const *,
                                             real3 *, real3 *, dem_aux const *);
typedef void (*WrapZ_ptr_t)(DynParticleCU *, ParticleCU const *, dem_aux const *);
// Quaternion exponential map (body‑frame angular velocity → rotation increment)
__device__ inline real4 Exp_GPU(real3 w, real dt)
{
    real theta = sqrt(w.x*w.x + w.y*w.y + w.z*w.z) * dt;
    real4 q;
    if (theta < 1e-12) {
        q = MakeQuat(1.0,make_real3(0.0,0.0,0.0));
    } else {
        real s = sin(theta * 0.5);
        real c = cos(theta * 0.5);
        real ax = w.x / sqrt(w.x*w.x + w.y*w.y + w.z*w.z);
        real ay = w.y / sqrt(w.x*w.x + w.y*w.y + w.z*w.z);
        real az = w.z / sqrt(w.x*w.x + w.y*w.y + w.z*w.z);
        q = MakeQuat(c,make_real3(s*ax,s*ay,s*az));
    }
    return q;
}

__device__ inline real4 QMult(real4 q1, real4 q2)
{
    real s1 = QScalar(q1);
    real s2 = QScalar(q2);
    real3 v1 = QVector(q1);
    real3 v2 = QVector(q2);
    real scalar = s1*s2 - dotreal3(v1,v2);
    real3 vector = s1*v2 + s2*v1 + cross(v1,v2);
    return MakeQuat(scalar,vector);
}



// =========================================================================
// Force calculation kernels
// =========================================================================

__global__ void CalcForceVV(InteractonCU const * Int, ComInteractonCU * CInt, DynInteractonCU * DIntVV, ParticleCU * Par,
        DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nvvint) return;
    size_t id = DIntVV[ic].Idx;
    size_t i1 = CInt [id].I1;
    size_t i2 = CInt [id].I2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    real3  Branch;
    if (Int[id].BothFree) BranchVec(xf, xi, Branch, demaux[0].Per);
    else Branch = xi - xf;

    real dist  = norm(Branch);
    real delta = r1 + r2 - dist;

    DIntVV[ic].Fn = make_real3(0.0, 0.0, 0.0);
    if (delta <= 0.0) {
        DIntVV[ic].Ft = make_real3(0.0, 0.0, 0.0);
        DIntVV[ic].InContact = false;
        DIntVV[ic].PrevDelta = delta;
    }

    if (delta > 0.0)
    {
        bool newContact = !DIntVV[ic].InContact;
        real prevDelta = DIntVV[ic].PrevDelta;
        DIntVV[ic].InContact = true;

        real3  n   = -1.0 * Branch / dist;
        real   d   = (r1*r1 - r2*r2 + dist*dist) / (2.0 * dist);
        real3  x1c = xi + d * n;
        real3  x2c = xf - (dist - d) * n;

        real3 t1, t2, x1, x2;
        Rotation(DPar[i1].w, DPar[i1].Q, t1);
        Rotation(DPar[i2].w, DPar[i2].Q, t2);
        x1 = x1c - xi;
        x2 = x2c - xf;
        real3 vrel = (DPar[i1].v + cross(t1, x1)) - (DPar[i2].v + cross(t2, x2));
        real3 vt   = vrel - dotreal3(n, vrel) * n;

        real3 Fn_elastic = Int[id].Kn * delta * n;
        real3 Fn_dashpot = Int[id].Gn * dotreal3(n, vrel) * n;
        real3 Fn_total = Fn_elastic + Fn_dashpot;
        if (demaux[0].sphereTensileCutoff && dotreal3(Fn_total, n) < 0.0) {
            Fn_total = make_real3(0.0, 0.0, 0.0);
            Fn_dashpot = -1.0 * Fn_elastic;
        }
        DIntVV[ic].Fn = Fn_total;

        real dtTangential = demaux[0].dt;
        if (demaux[0].sphereFirstContactCorrection && newContact && prevDelta < 0.0 && delta > prevDelta) {
            real activationFraction = delta/(delta - prevDelta);
            if (activationFraction < 0.0) activationFraction = 0.0;
            if (activationFraction > 1.0) activationFraction = 1.0;
            dtTangential *= activationFraction;
        }
        DIntVV[ic].Ft = DIntVV[ic].Ft + (Int[id].Kt * dtTangential) * vt;
        DIntVV[ic].Ft = DIntVV[ic].Ft - dotreal3(DIntVV[ic].Ft, n) * n;

        real3 Ft_elastic_trial = DIntVV[ic].Ft;
        real3 Ft_elastic = Ft_elastic_trial;
        real3 Ft_dashpot = Int[id].Gt * vt;
        real3 Ft_total = Ft_elastic + Ft_dashpot;
        size_t coulombMode = demaux[0].sphereCoulombMode;
        if (coulombMode > 2) coulombMode = 0;
        real friction_limit = Int[id].Mu
            * (coulombMode == 2 ? norm(Fn_total) : norm(Fn_elastic));
        bool sliding = false;
        bool tangentialDashpotApplied = true;
        if (coulombMode == 0 || coulombMode == 1) {
            real ftElasticNorm = norm(Ft_elastic);
            if (ftElasticNorm > friction_limit) {
                sliding = true;
                Ft_elastic = friction_limit * Ft_elastic / ftElasticNorm;
                DIntVV[ic].Ft = Ft_elastic;
                if (coulombMode == 1) {
                    Ft_dashpot = make_real3(0.0, 0.0, 0.0);
                    tangentialDashpotApplied = false;
                }
            }
            Ft_total = Ft_elastic + Ft_dashpot;
        } else {
            real ftTotalNorm = norm(Ft_total);
            if (ftTotalNorm > friction_limit) {
                sliding = true;
                real3 Ft_cap = friction_limit * Ft_total / ftTotalNorm;
                Ft_elastic = Ft_cap - Ft_dashpot;
                DIntVV[ic].Ft = Ft_elastic;
                Ft_total = Ft_cap;
            }
        }
        if (sliding && Int[id].Kt > 0.0) {
            real3 plastic_disp = (Ft_elastic_trial - Ft_elastic) / Int[id].Kt;
            CInt[id].dEfric = friction_limit * norm(plastic_disp);
        }

        CInt[id].Epot = 0.5 * Int[id].Kn * delta * delta
                        + (Int[id].Kt > 0.0
                           ? 0.5 * dotreal3(Ft_elastic, Ft_elastic) / Int[id].Kt
                           : 0.0);
        CInt[id].dEvis = (Int[id].Gn * dotreal3(vrel - vt, vrel - vt)
                          + (tangentialDashpotApplied ? Int[id].Gt * dotreal3(vt, vt) : 0.0))
                         * demaux[0].dt;

        real3 vr = r1 * r2 * cross((t1 - t2), n) / (r1 + r2);
        DIntVV[ic].Fr = DIntVV[ic].Fr + (Int[id].Beta * Int[id].Kt * demaux[0].dt) * vr;
        DIntVV[ic].Fr = DIntVV[ic].Fr - dotreal3(DIntVV[ic].Fr, n) * n;

        real3 tan = DIntVV[ic].Fr;
        if (norm(tan) > 0.0) tan = tan / norm(tan);
        if (norm(DIntVV[ic].Fr) > Int[id].Eta * Int[id].Mu * norm(Fn_elastic)) {
            DIntVV[ic].Fr = Int[id].Eta * Int[id].Mu * norm(Fn_elastic) * tan;
        }

        // total normal and tangential contact force
        DIntVV[ic].F = Fn_total + Ft_total;

        real3 T1, T2, T, Tt;
        Tt = cross(x1, DIntVV[ic].F) + r1 * cross(n, DIntVV[ic].Fr);
        real4 q;
        Conjugate(DPar[i1].Q, q);
        Rotation(Tt, q, T);
        T1 = -1.0 * T;
        Tt = cross(x2, DIntVV[ic].F) + r2 * cross(n, DIntVV[ic].Fr);
        Conjugate(DPar[i2].Q, q);
        Rotation(Tt, q, T);
        T2 = T;

        // accumulate ELASTIC into contact net for output
        atomicAdd(&CInt[id].Fnnet.x, Fn_elastic.x);
        atomicAdd(&CInt[id].Fnnet.y, Fn_elastic.y);
        atomicAdd(&CInt[id].Fnnet.z, Fn_elastic.z);
        atomicAdd(&CInt[id].Ftnet.x, Ft_elastic.x);
        atomicAdd(&CInt[id].Ftnet.y, Ft_elastic.y);
        atomicAdd(&CInt[id].Ftnet.z, Ft_elastic.z);

        // accumulate DASHPOT into contact dpot for output
        atomicAdd(&CInt[id].Fndpot.x, Fn_dashpot.x);
        atomicAdd(&CInt[id].Fndpot.y, Fn_dashpot.y);
        atomicAdd(&CInt[id].Fndpot.z, Fn_dashpot.z);
        atomicAdd(&CInt[id].Ftdpot.x, Ft_dashpot.x);
        atomicAdd(&CInt[id].Ftdpot.y, Ft_dashpot.y);
        atomicAdd(&CInt[id].Ftdpot.z, Ft_dashpot.z);


        // ignore the Fthermostat component, only updated if there is a thermostat active that overloads this function

        atomicAdd(&DPar[i1].F.x, -DIntVV[ic].F.x);
        atomicAdd(&DPar[i1].F.y, -DIntVV[ic].F.y);
        atomicAdd(&DPar[i1].F.z, -DIntVV[ic].F.z);
        atomicAdd(&DPar[i2].F.x,  DIntVV[ic].F.x);
        atomicAdd(&DPar[i2].F.y,  DIntVV[ic].F.y);
        atomicAdd(&DPar[i2].F.z,  DIntVV[ic].F.z);
        atomicAdd(&Par[i1].T.x, T1.x);
        atomicAdd(&Par[i1].T.y, T1.y);
        atomicAdd(&Par[i1].T.z, T1.z);
        atomicAdd(&Par[i2].T.x, T2.x);
        atomicAdd(&Par[i2].T.y, T2.y);
        atomicAdd(&Par[i2].T.z, T2.z);
        DIntVV[ic].PrevDelta = delta;
    }
}

__global__ void CalcForceVV_Hertz(InteractonCU const * Int, ComInteractonCU * CInt, DynInteractonCU * DIntVV, ParticleCU * Par,
        DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nvvint) return;
    size_t id = DIntVV[ic].Idx;
    size_t i1 = CInt [id].I1;
    size_t i2 = CInt [id].I2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    real3  Branch;
    if (Int[id].BothFree) BranchVec(xf, xi, Branch, demaux[0].Per);
    else Branch = xi - xf;

    real dist  = norm(Branch);
    real delta = r1 + r2 - dist;

    DIntVV[ic].Fn = make_real3(0.0, 0.0, 0.0);
    if (delta <= 0.0) {
        DIntVV[ic].Ft = make_real3(0.0, 0.0, 0.0);
        DIntVV[ic].InContact = false;
        DIntVV[ic].PrevDelta = delta;
    }

    if (delta > 0.0)
    {
        bool newContact = !DIntVV[ic].InContact;
        real prevDelta = DIntVV[ic].PrevDelta;
        DIntVV[ic].InContact = true;

        real3  n   = -1.0 * Branch / dist;
        real   d   = (r1*r1 - r2*r2 + dist*dist) / (2.0 * dist);
        real3  x1c = xi + d * n;
        real3  x2c = xf - (dist - d) * n;

        real3 t1, t2, x1, x2;
        Rotation(DPar[i1].w, DPar[i1].Q, t1);
        Rotation(DPar[i2].w, DPar[i2].Q, t2);
        x1 = x1c - xi;
        x2 = x2c - xf;
        real3 vrel = (DPar[i1].v + cross(t1, x1)) - (DPar[i2].v + cross(t2, x2));
        real3 vt   = vrel - dotreal3(n, vrel) * n;

        real sqrtDelta = sqrt(delta);
        real tangentialStiffness = Int[id].Kt * sqrtDelta;
        real normalDamping = Int[id].Gn * sqrt(sqrtDelta);
        real tangentialDamping = Int[id].Gt * sqrt(sqrtDelta);
        real3 Fn_elastic = Int[id].Kn * sqrtDelta * delta * n;
        real3 Fn_dashpot = normalDamping * dotreal3(n, vrel) * n;
        real3 Fn_total = Fn_elastic + Fn_dashpot;
        if (demaux[0].sphereTensileCutoff && dotreal3(Fn_total, n) < 0.0) {
            Fn_total = make_real3(0.0, 0.0, 0.0);
            Fn_dashpot = -1.0 * Fn_elastic;
        }
        DIntVV[ic].Fn = Fn_total;

        real dtTangential = demaux[0].dt;
        if (demaux[0].sphereFirstContactCorrection && newContact && prevDelta < 0.0 && delta > prevDelta) {
            real activationFraction = delta/(delta - prevDelta);
            if (activationFraction < 0.0) activationFraction = 0.0;
            if (activationFraction > 1.0) activationFraction = 1.0;
            dtTangential *= activationFraction;
        }
        // The CPU stores tangential displacement; the device stores elastic force.
        // Rescale the old force when Hertz stiffness changes with overlap.
        if (prevDelta > 0.0 && Int[id].Kt > 0.0)
            DIntVV[ic].Ft = (sqrtDelta / sqrt(prevDelta)) * DIntVV[ic].Ft;
        else
            DIntVV[ic].Ft = make_real3(0.0, 0.0, 0.0);
        DIntVV[ic].Ft = DIntVV[ic].Ft + (tangentialStiffness * dtTangential) * vt;
        DIntVV[ic].Ft = DIntVV[ic].Ft - dotreal3(DIntVV[ic].Ft, n) * n;

        real3 Ft_elastic_trial = DIntVV[ic].Ft;
        real3 Ft_elastic = Ft_elastic_trial;
        real3 Ft_dashpot = tangentialDamping * vt;
        real3 Ft_total = Ft_elastic + Ft_dashpot;
        size_t coulombMode = demaux[0].sphereCoulombMode;
        if (coulombMode > 2) coulombMode = 0;
        real friction_limit = Int[id].Mu
            * (coulombMode == 2 ? norm(Fn_total) : norm(Fn_elastic));
        bool sliding = false;
        bool tangentialDashpotApplied = true;
        if (coulombMode == 0 || coulombMode == 1) {
            real ftElasticNorm = norm(Ft_elastic);
            if (ftElasticNorm > friction_limit) {
                sliding = true;
                Ft_elastic = friction_limit * Ft_elastic / ftElasticNorm;
                DIntVV[ic].Ft = Ft_elastic;
                if (coulombMode == 1) {
                    Ft_dashpot = make_real3(0.0, 0.0, 0.0);
                    tangentialDashpotApplied = false;
                }
            }
            Ft_total = Ft_elastic + Ft_dashpot;
        } else {
            real ftTotalNorm = norm(Ft_total);
            if (ftTotalNorm > friction_limit) {
                sliding = true;
                real3 Ft_cap = friction_limit * Ft_total / ftTotalNorm;
                Ft_elastic = Ft_cap - Ft_dashpot;
                DIntVV[ic].Ft = Ft_elastic;
                Ft_total = Ft_cap;
            }
        }
        if (sliding && tangentialStiffness > 0.0) {
            real3 plastic_disp = (Ft_elastic_trial - Ft_elastic) / tangentialStiffness;
            CInt[id].dEfric = friction_limit * norm(plastic_disp);
        }

        CInt[id].Epot = 0.4 * Int[id].Kn * delta * delta * sqrtDelta
                        + (tangentialStiffness > 0.0
                           ? 0.5 * dotreal3(Ft_elastic, Ft_elastic) / tangentialStiffness
                           : 0.0);
        CInt[id].dEvis = (normalDamping * dotreal3(vrel - vt, vrel - vt)
                          + (tangentialDashpotApplied ? tangentialDamping * dotreal3(vt, vt) : 0.0))
                         * demaux[0].dt;

        real3 vr = r1 * r2 * cross((t1 - t2), n) / (r1 + r2);
        DIntVV[ic].Fr = DIntVV[ic].Fr + (Int[id].Beta * tangentialStiffness * demaux[0].dt) * vr;
        DIntVV[ic].Fr = DIntVV[ic].Fr - dotreal3(DIntVV[ic].Fr, n) * n;

        real3 tan = DIntVV[ic].Fr;
        if (norm(tan) > 0.0) tan = tan / norm(tan);
        if (norm(DIntVV[ic].Fr) > Int[id].Eta * Int[id].Mu * norm(Fn_elastic)) {
            DIntVV[ic].Fr = Int[id].Eta * Int[id].Mu * norm(Fn_elastic) * tan;
        }

        // total normal and tangential contact force
        DIntVV[ic].F = Fn_total + Ft_total;

        real3 T1, T2, T, Tt;
        Tt = cross(x1, DIntVV[ic].F) + r1 * cross(n, DIntVV[ic].Fr);
        real4 q;
        Conjugate(DPar[i1].Q, q);
        Rotation(Tt, q, T);
        T1 = -1.0 * T;
        Tt = cross(x2, DIntVV[ic].F) + r2 * cross(n, DIntVV[ic].Fr);
        Conjugate(DPar[i2].Q, q);
        Rotation(Tt, q, T);
        T2 = T;

        // accumulate ELASTIC into contact net for output
        atomicAdd(&CInt[id].Fnnet.x, Fn_elastic.x);
        atomicAdd(&CInt[id].Fnnet.y, Fn_elastic.y);
        atomicAdd(&CInt[id].Fnnet.z, Fn_elastic.z);
        atomicAdd(&CInt[id].Ftnet.x, Ft_elastic.x);
        atomicAdd(&CInt[id].Ftnet.y, Ft_elastic.y);
        atomicAdd(&CInt[id].Ftnet.z, Ft_elastic.z);

        // accumulate DASHPOT into contact dpot for output
        atomicAdd(&CInt[id].Fndpot.x, Fn_dashpot.x);
        atomicAdd(&CInt[id].Fndpot.y, Fn_dashpot.y);
        atomicAdd(&CInt[id].Fndpot.z, Fn_dashpot.z);
        atomicAdd(&CInt[id].Ftdpot.x, Ft_dashpot.x);
        atomicAdd(&CInt[id].Ftdpot.y, Ft_dashpot.y);
        atomicAdd(&CInt[id].Ftdpot.z, Ft_dashpot.z);


        // ignore the Fthermostat component, only updated if there is a thermostat active that overloads this function

        atomicAdd(&DPar[i1].F.x, -DIntVV[ic].F.x);
        atomicAdd(&DPar[i1].F.y, -DIntVV[ic].F.y);
        atomicAdd(&DPar[i1].F.z, -DIntVV[ic].F.z);
        atomicAdd(&DPar[i2].F.x,  DIntVV[ic].F.x);
        atomicAdd(&DPar[i2].F.y,  DIntVV[ic].F.y);
        atomicAdd(&DPar[i2].F.z,  DIntVV[ic].F.z);
        atomicAdd(&Par[i1].T.x, T1.x);
        atomicAdd(&Par[i1].T.y, T1.y);
        atomicAdd(&Par[i1].T.z, T1.z);
        atomicAdd(&Par[i2].T.x, T2.x);
        atomicAdd(&Par[i2].T.y, T2.y);
        atomicAdd(&Par[i2].T.z, T2.z);
        DIntVV[ic].PrevDelta = delta;
    }
}

__global__ void CalcForceEE(size_t const * Edges, real3 const * Verts, InteractonCU const * Int, ComInteractonCU * CInt, DynInteractonCU * DIntEE,
        ParticleCU * Par, DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].neeint) return;
    size_t id = DIntEE[ic].Idx;
    size_t i1 = CInt  [id].I1;
    size_t i2 = CInt  [id].I2;
    size_t f1 = DIntEE[ic].IF1;
    size_t f2 = DIntEE[ic].IF2;
    real  dm1 = DIntEE[ic].Dmax1;
    real  dm2 = DIntEE[ic].Dmax2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    
    DIntEE[ic].Fn = make_real3(0.0,0.0,0.0);
    
    real3 s;
    real3 Pert = make_real3(0.0,0.0,0.0);
    if (Int[id].BothFree) Pert = demaux[0].Per;
    if (!OverlapEE(Edges,Verts,f1,f2,dm1,dm2,Pert)) return;
    DistanceEE(Edges,Verts,f1,f2,xi,xf,s,Pert);
    real dist  = norm(s);
    real delta = r1 + r2 - dist;
    if (delta>0)
    {
        real3  n = s/dist;
        real   d = (r1*r1-r2*r2+dist*dist)/(2*dist);
        real3  x1c = xi+d*n;
        real3  x2c = xf-(dist-d)*n;
        real3 t1,t2,x1,x2;
        Rotation(DPar[i1].w,DPar[i1].Q,t1);
        Rotation(DPar[i2].w,DPar[i2].Q,t2);
        x1 = x1c - DPar [i1].x;
        x2 = x2c - DPar [i2].x;
        real3 vrel = (DPar[i1].v+cross(t1,x1))-(DPar[i2].v+cross(t2,x2));
        real3 vt   = vrel - dotreal3(n,vrel)*n;

        DIntEE[ic].Fn  = Int[id].Kn*delta*n;
        DIntEE[ic].Ft  = DIntEE[ic].Ft + (Int[id].Kt*demaux[0].dt)*vt;
        DIntEE[ic].Ft  = DIntEE[ic].Ft - dotreal3(DIntEE[ic].Ft,n)*n;

        real3 tan = DIntEE[ic].Ft;
        if (norm(tan)>0.0) tan = tan/norm(tan);
        if (norm(DIntEE[ic].Ft)>Int[id].Mu*norm(DIntEE[ic].Fn))
        {
            DIntEE[ic].Ft = Int[id].Mu*norm(DIntEE[ic].Fn)*tan;
        }

        DIntEE[ic].F = DIntEE[ic].Fn + DIntEE[ic].Ft + Int[id].Gn*dotreal3(n,vrel)*n + Int[id].Gt*vt;

        real3 T1,T2,T, Tt;
        Tt = cross (x1,DIntEE[ic].F);
        real4 q;
        Conjugate (DPar[i1].Q,q);
        Rotation  (Tt,q,T);
        T1 = -1.0*T;
        Tt = cross (x2,DIntEE[ic].F);
        Conjugate (DPar[i2].Q,q);
        Rotation  (Tt,q,T);
        T2 =      T;

        atomicAdd(&CInt[id].Fnnet.x, DIntEE[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntEE[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntEE[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntEE[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntEE[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntEE[ic].Ft.z);

        atomicAdd(&DPar[i1].F.x,-DIntEE[ic].F .x);
        atomicAdd(&DPar[i1].F.y,-DIntEE[ic].F .y);
        atomicAdd(&DPar[i1].F.z,-DIntEE[ic].F .z);
        atomicAdd(&DPar[i2].F.x, DIntEE[ic].F .x);
        atomicAdd(&DPar[i2].F.y, DIntEE[ic].F .y);
        atomicAdd(&DPar[i2].F.z, DIntEE[ic].F .z);
        atomicAdd(& Par[i1].T.x,            T1.x);
        atomicAdd(& Par[i1].T.y,            T1.y);
        atomicAdd(& Par[i1].T.z,            T1.z);
        atomicAdd(& Par[i2].T.x,            T2.x);
        atomicAdd(& Par[i2].T.y,            T2.y);
        atomicAdd(& Par[i2].T.z,            T2.z);
    }
}

__global__ void CalcForceVF(size_t const * Faces, size_t const * Facid, real3 const * Verts, InteractonCU const * Int, ComInteractonCU * CInt,
        DynInteractonCU * DIntVF, ParticleCU * Par, DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nvfint) return;
    size_t id = DIntVF[ic].Idx;
    size_t i1 = CInt  [id].I1;
    size_t i2 = CInt  [id].I2;
    size_t f1 = DIntVF[ic].IF1;
    size_t f2 = DIntVF[ic].IF2;
    real  dm1 = DIntVF[ic].Dmax1;
    real  dm2 = DIntVF[ic].Dmax2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    
    DIntVF[ic].Fn = make_real3(0.0,0.0,0.0);

    real3 s;
    xi = Verts[f1];
    real3 Pert = make_real3(0.0,0.0,0.0);
    if (Int[id].BothFree) Pert = demaux[0].Per;
    if (!OverlapVF(Faces,Facid,Verts,xi,f2,dm1,dm2,Pert)) return;
    DistanceVF(Faces,Facid,Verts,xi,f2,xf,s,Pert);
    real dist  = norm(s);
    real delta = r1 + r2 - dist;
    if (delta>0)
    {
        real3  n = s/dist;
        real   d = (r1*r1-r2*r2+dist*dist)/(2*dist);
        real3  x1c = xi+d*n;
        real3  x2c = xf-(dist-d)*n;
        real3 t1,t2,x1,x2;
        Rotation(DPar[i1].w,DPar[i1].Q,t1);
        Rotation(DPar[i2].w,DPar[i2].Q,t2);

        x1 = x1c - DPar [i1].x;
        x2 = x2c - DPar [i2].x;
        real3 vrel = (DPar[i1].v+cross(t1,x1))-(DPar[i2].v+cross(t2,x2));
        real3 vt   = vrel - dotreal3(n,vrel)*n;

        DIntVF[ic].Fn  = Int[id].Kn*delta*n;
        DIntVF[ic].Ft  = DIntVF[ic].Ft + (Int[id].Kt*demaux[0].dt)*vt;
        DIntVF[ic].Ft  = DIntVF[ic].Ft - dotreal3(DIntVF[ic].Ft,n)*n;

        real3 tan = DIntVF[ic].Ft;
        if (norm(tan)>0.0) tan = tan/norm(tan);
        if (norm(DIntVF[ic].Ft)>Int[id].Mu*norm(DIntVF[ic].Fn))
        {
            DIntVF[ic].Ft = Int[id].Mu*norm(DIntVF[ic].Fn)*tan;
        }

        DIntVF[ic].F = DIntVF[ic].Fn + DIntVF[ic].Ft + Int[id].Gn*dotreal3(n,vrel)*n + Int[id].Gt*vt;

        real3 T1,T2,T, Tt;
        Tt = cross (x1,DIntVF[ic].F);
        real4 q;
        Conjugate (DPar[i1].Q,q);
        Rotation  (Tt,q,T);
        T1 = -1.0*T;
        Tt = cross (x2,DIntVF[ic].F);
        Conjugate (DPar[i2].Q,q);
        Rotation  (Tt,q,T);
        T2 =      T;

        atomicAdd(&CInt[id].Fnnet.x, DIntVF[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntVF[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntVF[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntVF[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntVF[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntVF[ic].Ft.z);

        atomicAdd(&DPar[i1].F.x,-DIntVF[ic].F .x);
        atomicAdd(&DPar[i1].F.y,-DIntVF[ic].F .y);
        atomicAdd(&DPar[i1].F.z,-DIntVF[ic].F .z);
        atomicAdd(&DPar[i2].F.x, DIntVF[ic].F .x);
        atomicAdd(&DPar[i2].F.y, DIntVF[ic].F .y);
        atomicAdd(&DPar[i2].F.z, DIntVF[ic].F .z);
        atomicAdd(& Par[i1].T.x,            T1.x);
        atomicAdd(& Par[i1].T.y,            T1.y);
        atomicAdd(& Par[i1].T.z,            T1.z);
        atomicAdd(& Par[i2].T.x,            T2.x);
        atomicAdd(& Par[i2].T.y,            T2.y);
        atomicAdd(& Par[i2].T.z,            T2.z);
    }
}

__global__ void CalcForceFV(size_t const * Faces, size_t const * Facid, real3 const * Verts, InteractonCU const * Int, ComInteractonCU * CInt,
        DynInteractonCU * DIntFV, ParticleCU * Par, DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nfvint) return;
    size_t id = DIntFV[ic].Idx;
    size_t i1 = CInt  [id].I1;
    size_t i2 = CInt  [id].I2;
    size_t f1 = DIntFV[ic].IF1;
    size_t f2 = DIntFV[ic].IF2;
    real  dm1 = DIntFV[ic].Dmax1;
    real  dm2 = DIntFV[ic].Dmax2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    
    DIntFV[ic].Fn = make_real3(0.0,0.0,0.0);

    real3 s;
    xf = Verts[f2];
    real3 Pert = make_real3(0.0,0.0,0.0);
    if (Int[id].BothFree) Pert = demaux[0].Per;
    if (!OverlapFV(Faces,Facid,Verts,f1,xf,dm1,dm2,Pert)) return;
    DistanceFV(Faces,Facid,Verts,f1,xf,xi,s,Pert);
    real dist  = norm(s);
    real delta = r1 + r2 - dist;
    if (delta>0)
    {
        real3  n = s/dist;
        real   d = (r1*r1-r2*r2+dist*dist)/(2*dist);
        real3  x1c = xi+d*n;
        real3  x2c = xf-(dist-d)*n;
        real3 t1,t2,x1,x2;
        Rotation(DPar[i1].w,DPar[i1].Q,t1);
        Rotation(DPar[i2].w,DPar[i2].Q,t2);
        x1 = x1c - DPar [i1].x;
        x2 = x2c - DPar [i2].x;
        real3 vrel = (DPar[i1].v+cross(t1,x1))-(DPar[i2].v+cross(t2,x2));
        real3 vt   = vrel - dotreal3(n,vrel)*n;

        DIntFV[ic].Fn  = Int[id].Kn*delta*n;
        DIntFV[ic].Ft  = DIntFV[ic].Ft + (Int[id].Kt*demaux[0].dt)*vt;
        DIntFV[ic].Ft  = DIntFV[ic].Ft - dotreal3(DIntFV[ic].Ft,n)*n;

        real3 tan = DIntFV[ic].Ft;
        if (norm(tan)>0.0) tan = tan/norm(tan);
        if (norm(DIntFV[ic].Ft)>Int[id].Mu*norm(DIntFV[ic].Fn))
        {
            DIntFV[ic].Ft = Int[id].Mu*norm(DIntFV[ic].Fn)*tan;
        }

        DIntFV[ic].F = DIntFV[ic].Fn + DIntFV[ic].Ft + Int[id].Gn*dotreal3(n,vrel)*n + Int[id].Gt*vt;

 
        
        real3 T1,T2,T, Tt;
        Tt = cross (x1,DIntFV[ic].F);
        real4 q;
        Conjugate (DPar[i1].Q,q);
        Rotation  (Tt,q,T);
        T1 = -1.0*T;
        Tt = cross (x2,DIntFV[ic].F);
        Conjugate (DPar[i2].Q,q);
        Rotation  (Tt,q,T);
        T2 =      T;

        atomicAdd(&CInt[id].Fnnet.x, DIntFV[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntFV[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntFV[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntFV[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntFV[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntFV[ic].Ft.z);

        atomicAdd(&DPar[i1].F.x,-DIntFV[ic].F .x);
        atomicAdd(&DPar[i1].F.y,-DIntFV[ic].F .y);
        atomicAdd(&DPar[i1].F.z,-DIntFV[ic].F .z);
        atomicAdd(&DPar[i2].F.x, DIntFV[ic].F .x);
        atomicAdd(&DPar[i2].F.y, DIntFV[ic].F .y);
        atomicAdd(&DPar[i2].F.z, DIntFV[ic].F .z);
        atomicAdd(& Par[i1].T.x,            T1.x);
        atomicAdd(& Par[i1].T.y,            T1.y);
        atomicAdd(& Par[i1].T.z,            T1.z);
        atomicAdd(& Par[i2].T.x,            T2.x);
        atomicAdd(& Par[i2].T.y,            T2.y);
        atomicAdd(& Par[i2].T.z,            T2.z);
    }
}


/// Set the velocity of a single particle directly on the device.
/// Needed for prescribed-motion boundary conditions: MechSys's only host->device
/// state path is the full re-upload, which also overwrites the accelerations, so
/// a per-step velocity prescription requires this one-line kernel.
__global__ void SetVelocity (DynParticleCU * DPar, size_t idx, real3 v)
{
    DPar[idx].v = v;
}

__global__ void VerletStep1(real3 * Verts, ParticleCU const * Par, DynParticleCU * DPar,
                            real3 const * A, dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nparts) return;



    // Use the acceleration from the previous step (already includes damping)
    real3 a_prev = A[ic];
    if (Par[ic].vxf) a_prev.x = 0.0f;
    if (Par[ic].vyf) a_prev.y = 0.0f;
    if (Par[ic].vzf) a_prev.z = 0.0f;

    real3 v = DPar[ic].v;
    real3 v_half = make_real3(v.x + 0.5f * demaux[0].dt * a_prev.x,
                              v.y + 0.5f * demaux[0].dt * a_prev.y,
                              v.z + 0.5f * demaux[0].dt * a_prev.z);

    real3 x_new = make_real3(DPar[ic].x.x + v_half.x * demaux[0].dt,
                             DPar[ic].x.y + v_half.y * demaux[0].dt,
                             DPar[ic].x.z + v_half.z * demaux[0].dt);

    real3 dis = make_real3(x_new.x - DPar[ic].x.x,
                           x_new.y - DPar[ic].x.y,
                           x_new.z - DPar[ic].x.z);

    DPar[ic].x = x_new;
    for (size_t iv = Par[ic].Nvi; iv < Par[ic].Nvf; iv++) {
        Verts[iv].x += dis.x;
        Verts[iv].y += dis.y;
        Verts[iv].z += dis.z;
    }

    DPar[ic].v = v_half;
}
__global__ void OrientationUpdate(real3 * Verts, ParticleCU const * Par,
                                  DynParticleCU * DPar, real3 const * Wdot,
                                  dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nparts) return;

    // Use the angular acceleration from the previous step
    real3 wdot_prev = Wdot[ic];
    if (Par[ic].wxf) wdot_prev.x = 0.0f;
    if (Par[ic].wyf) wdot_prev.y = 0.0f;
    if (Par[ic].wzf) wdot_prev.z = 0.0f;

    real3 w = DPar[ic].w;

    real3 w_half;
    w_half.x = w.x + 0.5f * demaux[0].dt * wdot_prev.x;
    w_half.y = w.y + 0.5f * demaux[0].dt * wdot_prev.y;
    w_half.z = w.z + 0.5f * demaux[0].dt * wdot_prev.z;

    // exponential map to update quaternion
    real4 dq = Exp_GPU(w_half, demaux[0].dt);
    real4 Q_old = DPar[ic].Q;
    real4 Q_new = QMult(Q_old, dq);
    real inv_norm = rsqrt(Q_new.x*Q_new.x + Q_new.y*Q_new.y +
                          Q_new.z*Q_new.z + Q_new.w*Q_new.w);
    Q_new.x *= inv_norm; Q_new.y *= inv_norm;
    Q_new.z *= inv_norm; Q_new.w *= inv_norm;

    // rotate vertices only if angular velocity is non-zero
    if (w_half.x != 0.0f || w_half.y != 0.0f || w_half.z != 0.0f)
    {
    real4 Q_old_conj;
    Conjugate(Q_old,Q_old_conj);
    for (size_t iv = Par[ic].Nvi; iv < Par[ic].Nvf; iv++) {
        real3 xt = make_real3(Verts[iv].x - DPar[ic].x.x,
                              Verts[iv].y - DPar[ic].x.y,
                              Verts[iv].z - DPar[ic].x.z);
        real3 xt_body;
        Rotation(xt, Q_old_conj, xt_body);
        real3 xt_new;
        Rotation(xt_body, Q_new, xt_new);
        Verts[iv].x = xt_new.x + DPar[ic].x.x;
        Verts[iv].y = xt_new.y + DPar[ic].x.y;
        Verts[iv].z = xt_new.z + DPar[ic].x.z;
    }
    }

    DPar[ic].Q = Q_new;
    DPar[ic].w = w_half;
}
__global__ void FinalizeVelocity(ParticleCU const * Par, DynParticleCU * DPar,
                                 real3 * A, dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nparts) return;

    real3 F = DPar[ic].F;               // force from current step (already computed)
    if (Par[ic].vxf) F.x = 0.0f;
    if (Par[ic].vyf) F.y = 0.0f;
    if (Par[ic].vzf) F.z = 0.0f;

    // damping uses half‑step velocity (stored in DPar[ic].v)
    real3 v_half = DPar[ic].v;
    if (Par[ic].Gv>0){
        F.x -= Par[ic].Gv * Par[ic].m * v_half.x;
        F.y -= Par[ic].Gv * Par[ic].m * v_half.y;
        F.z -= Par[ic].Gv * Par[ic].m * v_half.z;
    }

    real3 a_new = make_real3(F.x / Par[ic].m,
                             F.y / Par[ic].m,
                             F.z / Par[ic].m);

    // full velocity: v = v_half + 0.5*dt*a_new
    DPar[ic].v = make_real3(v_half.x + 0.5f * demaux[0].dt * a_new.x,
                            v_half.y + 0.5f * demaux[0].dt * a_new.y,
                            v_half.z + 0.5f * demaux[0].dt * a_new.z);

    // Store the same damped acceleration used in the full-step update so the
    // next half-step matches the CPU velocity-Verlet path.
    A[ic] = a_new;
}

__global__ void FinalizeRotation(ParticleCU const * Par, DynParticleCU * DPar,
                                 real3 * Wdot, dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nparts) return;

    real3 T = Par[ic].T;               // torque from current step
    if (Par[ic].wxf) T.x = 0.0f;
    if (Par[ic].wyf) T.y = 0.0f;
    if (Par[ic].wzf) T.z = 0.0f;

    real3 w_half = DPar[ic].w;         // half‑step angular velocity
    // damping
    if (Par[ic].Gm>0){
        T.x -= Par[ic].Gm * Par[ic].I.x * w_half.x;
        T.y -= Par[ic].Gm * Par[ic].I.y * w_half.y;
        T.z -= Par[ic].Gm * Par[ic].I.z * w_half.z;
    }

    real3 wdot_new;
    wdot_new.x = (T.x + (Par[ic].I.y - Par[ic].I.z) * w_half.y * w_half.z) / Par[ic].I.x;
    wdot_new.y = (T.y + (Par[ic].I.z - Par[ic].I.x) * w_half.z * w_half.x) / Par[ic].I.y;
    wdot_new.z = (T.z + (Par[ic].I.x - Par[ic].I.y) * w_half.x * w_half.y) / Par[ic].I.z;

    // full angular velocity
    DPar[ic].w = make_real3(
        w_half.x + 0.5f * demaux[0].dt * wdot_new.x,
        w_half.y + 0.5f * demaux[0].dt * wdot_new.y,
        w_half.z + 0.5f * demaux[0].dt * wdot_new.z);

    // Store the same damped angular acceleration used in the full-step update.
    Wdot[ic] = wdot_new;
}

__global__ void EnforceAngularFixity(ParticleCU const * Par, DynParticleCU * DPar,
                                     dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nparts) return;
    if (Par[ic].wxf) DPar[ic].w.x = 0.0f;
    if (Par[ic].wyf) DPar[ic].w.y = 0.0f;
    if (Par[ic].wzf) DPar[ic].w.z = 0.0f;
}

// OLD TRANSLATE AND ROTATE METHODS
__global__ void Translate(real3 * Verts, ParticleCU const * Par, DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nparts) return;
    real3 Ft = DPar[ic].F;
    if (Par[ic].vxf) Ft.x = 0.0;
    if (Par[ic].vyf) Ft.y = 0.0;
    if (Par[ic].vzf) Ft.z = 0.0;
    //if (ic==0) printf("Fz = %g \n",DPar[ic].Flbm.z);
    real3 temp,xa;
    xa    = 2.0*DPar[ic].x - DPar[ic].xb + (demaux[0].dt*demaux[0].dt/Par[ic].m)*Ft;
    temp  = xa - DPar[ic].x;
    DPar[ic].v    = 0.5*(xa - DPar[ic].xb)/demaux[0].dt;
    DPar [ic].xb  = DPar[ic].x;
    DPar[ic].x    = xa;


    //if (isnan(norm(DPar[ic].x))) printf("ic %lu it %lu \n",ic,demaux[0].iter);

    for (size_t iv=Par[ic].Nvi;iv<Par[ic].Nvf;iv++)
    {
        Verts [iv] = Verts [iv] + temp;
    }
}

__global__ void Rotate(real3 * Verts, ParticleCU const * Par, DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nparts) return;
    real q0,q1,q2,q3,wx,wy,wz;
    real3 qv;

    q0 = 0.5*QScalar(DPar[ic].Q);
    qv = QVector(DPar[ic].Q);
    q1 = 0.5*qv.x;
    q2 = 0.5*qv.y;
    q3 = 0.5*qv.z;

    real3 Tt = Par[ic].T;

    if (Par[ic].wxf) Tt.x = 0.0;
    if (Par[ic].wyf) Tt.y = 0.0;
    if (Par[ic].wzf) Tt.z = 0.0;

    //if (isnan(norm(Tt)))
    //{
        //printf("ic %d \n",ic);
        //printf("it %d \n",demaux[0].iter);
    //}

    DPar[ic].wa.x=(Tt.x+(Par[ic].I.y-Par[ic].I.z)*DPar[ic].wb.y*DPar[ic].wb.z)/Par[ic].I.x;
    DPar[ic].wa.y=(Tt.y+(Par[ic].I.z-Par[ic].I.x)*DPar[ic].wb.x*DPar[ic].wb.z)/Par[ic].I.y;
    DPar[ic].wa.z=(Tt.z+(Par[ic].I.x-Par[ic].I.y)*DPar[ic].wb.y*DPar[ic].wb.x)/Par[ic].I.z;
    DPar[ic].w = DPar[ic].wb+0.5*demaux[0].dt*DPar[ic].wa;

    wx = DPar[ic].w.x;
    wy = DPar[ic].w.y;
    wz = DPar[ic].w.z;
    real4 dq,qm;
    SetQuaternion (-(q1*wx+q2*wy+q3*wz),
                   make_real3(q0*wx-q3*wy+q2*wz,
                              q3*wx+q0*wy-q1*wz,
                              -q2*wx+q1*wy+q0*wz),
                   dq);

    DPar[ic].wb  = DPar[ic].wb+demaux[0].dt*DPar[ic].wa;
    qm  = DPar[ic].Q+(0.5*demaux[0].dt)*dq;

    q0 = 0.5*QScalar(qm);
    qv = QVector(qm);
    q1 = 0.5*qv.x;
    q2 = 0.5*qv.y;
    q3 = 0.5*qv.z;

    wx  = DPar[ic].wb.x;
    wy  = DPar[ic].wb.y;
    wz  = DPar[ic].wb.z;
    
    SetQuaternion (-(q1*wx+q2*wy+q3*wz),
                   make_real3(q0*wx-q3*wy+q2*wz,
                              q3*wx+q0*wy-q1*wz,
                              -q2*wx+q1*wy+q0*wz),
                   dq);

    real4 Qd = qm+0.5*demaux[0].dt*dq,temp;
    Conjugate(DPar[ic].Q,temp);

    for (size_t iv=Par[ic].Nvi;iv<Par[ic].Nvf;iv++)
    {
        real3 xt = Verts[iv] - DPar[ic].x;
        Rotation(xt,temp,Verts[iv]);
        Verts[iv] = Verts[iv] + DPar[ic].x;
    }

    DPar[ic].Q = Qd/norm(Qd);

    for (size_t iv=Par[ic].Nvi;iv<Par[ic].Nvf;iv++)
    {
        real3 xt = Verts[iv] - DPar[ic].x;
        Rotation(xt,DPar[ic].Q,Verts[iv]);
        Verts[iv] = Verts[iv] + DPar[ic].x;
    }
}


__global__ void Reset (ParticleCU * Par, DynParticleCU * DPar, InteractonCU const * Int, ComInteractonCU * CInt, dem_aux const * demaux, void *
        extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic<demaux[0].nparts)
    {
        DPar[ic].F    = Par[ic].Ff;
        DPar[ic].Flbm = Par[ic].Flbmf;
        Par [ic].T    = Par[ic].Tf;
    }
    else if (ic<demaux[0].nparts+demaux[0].ncoint)
    {
        size_t id = ic-demaux[0].nparts;
        CInt[id].Fnnet = Int[id].Fnf;
        CInt[id].Ftnet = Int[id].Ftf;
        CInt[id].Fndpot = make_real3(0.0,0.0,0.0);
        CInt[id].Ftdpot = make_real3(0.0,0.0,0.0);
        CInt[id].Fther = make_real3(0.0,0.0,0.0);
        CInt[id].Epot = 0.0;
        CInt[id].dEvis = 0.0;
        CInt[id].dEfric = 0.0;
    }
    else return;
}

__global__ void MaxD(real3 const * Verts, real3 const * Vertso, real * maxd, dem_aux * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nverts) return;
    if (ic==0) 
    {
        demaux[0].Time += demaux[0].dt;
        demaux[0].iter++;
    }
    maxd[ic] = norm(Vertso[ic]-Verts[ic]);

    /*
    if (maxd[ic] > 1.0)
        printf("DEBUG MaxD_BIG: ic=%d Vso=(%g,%g,%g) Vs=(%g,%g,%g) d=%g\n",
               (int)ic, Vertso[ic].x,Vertso[ic].y,Vertso[ic].z,
               Verts[ic].x,Verts[ic].y,Verts[ic].z, maxd[ic]);
    if (ic == 0)
        printf("DEBUG MaxD: it=%lu Vso=(%g,%g,%g) Vs=(%g,%g,%g) d=%g\n",
               demaux[0].iter, Vertso[ic].x,Vertso[ic].y,Vertso[ic].z,
               Verts[ic].x,Verts[ic].y,Verts[ic].z, maxd[ic]);
               */
}
__global__ void ResetMaxD(real3 * Verts, real3 * Vertso, real * maxd, ParticleCU const * Par, DynParticleCU * DPar, dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nparts) return;

    bool isfree = ((!Par[ic].vxf&&!Par[ic].vyf&&!Par[ic].vzf&&!Par[ic].wxf&&!Par[ic].wyf&&!Par[ic].wzf)||Par[ic].FixFree);

    real3 dis = make_real3(0.0,0.0,0.0);
    bool wrapped = false;

    if (isfree)
    {
        if (demaux[0].px)
        {
            if (DPar[ic].x.x< demaux[0].Xmin) { dis.x = demaux[0].Xmax - demaux[0].Xmin; wrapped = true; }
            if (DPar[ic].x.x>=demaux[0].Xmax) { dis.x = demaux[0].Xmin - demaux[0].Xmax; wrapped = true; }
        }
        if (demaux[0].py)
        {
            if (DPar[ic].x.y< demaux[0].Ymin) { dis.y = demaux[0].Ymax - demaux[0].Ymin; wrapped = true; }
            if (DPar[ic].x.y>=demaux[0].Ymax) { dis.y = demaux[0].Ymin - demaux[0].Ymax; wrapped = true; }
        }
        if (demaux[0].pz && fabs(demaux[0].Per.z) > 0.0)
        {
            real shear_shift = demaux[0].strain * (demaux[0].Zmax - demaux[0].Zmin);
            if (DPar[ic].x.z< demaux[0].Zmin)
            {
                dis.z = demaux[0].Zmax - demaux[0].Zmin;
                dis.x += shear_shift;
                wrapped = true;
            }
            if (DPar[ic].x.z>=demaux[0].Zmax)
            {
                dis.z = demaux[0].Zmin - demaux[0].Zmax;
                dis.x -= shear_shift;
                wrapped = true;
            }
        }
    }

    DPar[ic].x  = DPar[ic].x  + dis;
    DPar[ic].xb = DPar[ic].xb + dis;

    for (size_t iv=Par[ic].Nvi;iv<Par[ic].Nvf;iv++)
    {
        Verts [iv] = Verts [iv] + dis;
        // Only reset displacement reference for wrapped particles.
        // Non-wrapped particles retain their accumulated displacement so
        // that MaxD can trigger a Verlet list rebuild when needed.
        if (wrapped)
        {
            Vertso[iv] = Verts [iv];
            maxd  [iv] = 0.0;
        }
    }

}

// Zero out accelerations AND velocities on device after contact rebuilding + wrapping.
// This prevents VerletStep1 from producing a large displacement from stale velocities
// computed from the pre-wrap force state.
__global__ void ZeroAccel(real3 * A, real3 * Wdot, DynParticleCU * DPar, dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nparts) return;
    A[ic]        = make_real3(0.0, 0.0, 0.0);
    Wdot[ic]     = make_real3(0.0, 0.0, 0.0);
    //DPar[ic].v   = make_real3(0.0, 0.0, 0.0);
    //DPar[ic].w   = make_real3(0.0, 0.0, 0.0);
}

// Recompute translational and angular accelerations from current forces/torques
// after a contact list rebuild. This prevents stale accelerations from being used
// in the next Verlet half‑step.
__global__ void RecomputeAccelerations(ParticleCU const * Par,
                                       DynParticleCU const * DPar,
                                       real3 * A,
                                       real3 * Wdot,
                                       dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nparts) return;

    // Translational acceleration from current net force
    real3 F = DPar[ic].F;
    if (Par[ic].vxf) F.x = 0.0f;
    if (Par[ic].vyf) F.y = 0.0f;
    if (Par[ic].vzf) F.z = 0.0f;

    real3 v = DPar[ic].v;
    if (Par[ic].Gv>0){
        F.x -= Par[ic].Gv * Par[ic].m * v.x;
        F.y -= Par[ic].Gv * Par[ic].m * v.y;
        F.z -= Par[ic].Gv * Par[ic].m * v.z;
    }
    A[ic] = make_real3(F.x / Par[ic].m,
                       F.y / Par[ic].m,
                       F.z / Par[ic].m);

    // Angular acceleration from current torque (body‑frame)
    real3 T = Par[ic].T;
    if (Par[ic].wxf) T.x = 0.0f;
    if (Par[ic].wyf) T.y = 0.0f;
    if (Par[ic].wzf) T.z = 0.0f;

    // Use current angular velocity (DPar[ic].w) for the inertia cross terms
    real3 w = DPar[ic].w;
    if (Par[ic].Gm>0){
        T.x -= Par[ic].Gm * Par[ic].I.x * w.x;
        T.y -= Par[ic].Gm * Par[ic].I.y * w.y;
        T.z -= Par[ic].Gm * Par[ic].I.z * w.z;
    }
    real Ix = Par[ic].I.x, Iy = Par[ic].I.y, Iz = Par[ic].I.z;

    real3 wdot_new;
    wdot_new.x = (T.x + (Iy - Iz) * w.y * w.z) / Ix;
    wdot_new.y = (T.y + (Iz - Ix) * w.z * w.x) / Iy;
    wdot_new.z = (T.z + (Ix - Iy) * w.x * w.y) / Iz;

    Wdot[ic] = wdot_new;
}


}
#endif //MECHSYS_DEM_CUH
