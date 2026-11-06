/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
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
 * ************************************************************************/

// -------------------------------------------------------------------------
// BUILD NOTE
//   This solver is the constant-volume Lees-Edwards shearbox.  It is built by the
//   MechSys CMake project as the tdem target of the same name.  It builds against
//   an unmodified MechSys; see the note below for exactly which library hooks it
//   uses and which library calls it deliberately avoids.
// -------------------------------------------------------------------------
//
// Three-phase periodic shearbox, CONSTANT VOLUME, Lees-Edwards kinematics.
//
// Ensemble: the cell is a fixed-size periodic box (Lx, Ly, Lz all prescribed),
// so the simulation is run at constant volume / constant solid fraction.  The
// pressure reported in <key>_shear.res is a *measured* observable (contact
// virial + kinetic term); there is no barostat, pressure servo or stress
// feedback anywhere, and the z-boundaries are held fixed throughout Phase 3.
// Phase 2 is a prescribed-rate densification (the two z-periodic boundaries
// are driven inward over a fixed comp_time to reach a target solid fraction),
// which is still a volume-controlled protocol, not a pressure-controlled one.
//
// Phases:
//   Phase 1: generate particles with MechSys dompargen helpers in a fully
//            periodic box.
//   Phase 2: compress to a target solid fraction by moving both z-periodic
//            boundaries inward, then translate the compressed cell to z=[0,Lz].
//            Driven by SolveCompression() below (see that function for why the
//            driver is local).
//   Phase 3: shear with Lees-Edwards shifted z-periodic boundaries, driven by
//            SolveLeesEdwardsShear() below.
//
// =========================================================================
// How this file couples to MechSys: what can and what cannot be replaced
// -------------------------------------------------------------------------
// MechSys' DEM is HEADER-ONLY (there is no lib/dem/*.cpp), so every library
// definition is compiled into this translation unit.  Two consequences drive
// the whole design of this file:
//
//   * A library definition CANNOT be replaced by re-defining it here: the
//     compiler sees the header definition first and reports a redefinition
//     (e.g. 'error: function "DEM::ResetMaxD" has already been defined').
//     There is also no library object file, so there is no symbol to interpose
//     at link time.  "Overloading" does not help either: overload resolution
//     happens at the call site by name and signature, and the library's call
//     sites are already written.
//
//   * The ONLY supported extension mechanism is Domain's table of function
//     pointers.  Anything launched through one of those pointers can be swapped
//     from here:
//
//       pForceVV  pForceEE  pForceVF  pForceFV    contact laws
//       pTranslate  pRotate                       integration
//       pReset                                    per-step reset
//       pVerletStep1  pFinalizeVelocity
//       pOrientationUpdate  pFinalizeRotation     velocity Verlet
//       pRecomputeAccelerations
//       pExtraParams                              (void*, shared state)
//
//   What this file installs, and nothing more:
//       pForceVV     -> CalcForceVV_LE / CalcForceVV_Hertz_LE
//       pReset       -> ResetLE
//       pExtraParams -> LeExtra   (Lees-Edwards state + contact virial sums)
//
// What this file CANNOT install, and what it does instead:
//
//   1. ResetMaxD and MaxD are NOT in the table.  They are free __global__
//      functions defined in dem.cuh and launched *by name* from library code
//      (Domain::Solve, Domain::UpdateContactsDevice), so both the definition and
//      every call site belong to the library.  The only options are a local
//      kernel under a different name, and calling it from a loop we own -- which
//      is exactly why SolveCompression() and SolveLeesEdwardsShear() exist.
//      (This file actually keeps calling the library ResetMaxD, which is correct
//      here: our WrapZ has already applied the Lees-Edwards x-shift before
//      ResetMaxD runs, so ResetMaxD only has to do the plain periodic wrap.)
//
//   2. dem_aux is a concrete object owned by MechSys.  Domain::demaux is a
//      member of that exact type, and Domain::pdemaux is a dem_aux* that the
//      library allocates (cudaMalloc inside Domain::UpLoadDevice), frees and
//      re-uploads wholesale (cudaMemcpy of sizeof(dem_aux)); every library launch
//      passes pdemaux.  Deriving from it (struct my_aux : dem_aux) does compile,
//      but MechSys would still allocate and upload the base-sized object, so an
//      added member would sit outside the storage MechSys manages and would be
//      invisible to library kernels.  Hence the Lees-Edwards state travels in a
//      separate block reached through the pre-existing void* pExtraParams.
//      Consequence to respect: pdemaux does not exist before the first
//      Domain::UpLoadDevice, so dem_aux must never be uploaded earlier than that
//      (see the NOTE in SolveCompression()).
//
//   3. The per-step wrap inside Domain::Solve is a hard-wired <<<>>> launch, not
//      a table entry, so it cannot be replaced from here at all.  This file
//      therefore owns the step loop and launches its own WrapZ directly, which
//      also lets the Lees-Edwards state travel in pExtraParams instead of having
//      to be added to dem_aux.
//
//   4. Library methods (Domain::Initialize, UpLoadDevice, DnLoadDevice,
//      UpdateContactsDevice, WriteXDMF, WriteBF, GenSpheresBox, SetProps, ...)
//      are not virtual and are not in the table: they can only be CALLED, never
//      replaced.  This file calls them and re-uses their data structures
//      (PairtoCInt, CInteractonSphere, FricCoeff, ContactLaw, ...).
//
// Net effect: this solver runs against an unmodified MechSys.  It needs no
// Lees-Edwards support inside the library -- no dem_aux strain/shear-rate
// fields, no extra wrap hook, no shear offset in the wrap, and no
// BranchVec(...,Per,strain) in distance.h.  All of that lives in this file.
// =========================================================================


#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/linalg/matvec.h>
#include <cuda_runtime_api.h>
#include <thrust/extrema.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <random>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <vector>

using std::cout;
using std::endl;

// =========================================================================
// UserData structure for compression + shear
// =========================================================================
struct UserData
{
    double         Lx;               // x-period
    double         Ly;               // y-period
    double         Lz_target;        // target box height after compression
    double         Lz_initial;       // initial box height before compression
    double         Zmin_initial;     // starting lower periodic boundary
    double         Zmin_target;      // final lower periodic boundary
    double         Zmax_initial;     // starting upper periodic boundary
    double         Zmax_target;      // final upper periodic boundary
    double         last_contact_Lz;  // last compression box height used for contacts
    double         comp_start_time;  // time at which compression began
    double         comp_duration;    // total duration of compression phase
    double         cum_strain;       // cumulative shear strain (host)
    double         reduced_strain;   // strain modulo x-period, report-only diagnostic
    double         shear_rate;       // current shear rate (0 during compression)
    double         base_alpha;       // input Verlet distance
    size_t         comp_rebuilds;    // number of boundary-driven compression rebuilds
    double         le_contact_strain;// continuous LE strain used for the last shear contact rebuild
    size_t         le_rebuilds;      // number of LE-aware shear contact rebuilds
    bool           compressing;      // true during compression phase
    bool           shearing;         // true during shear phase
    double         dt;              // time step
    double         dtOut;           // output interval
    std::ofstream  oss_comp;        // compression output
    std::ofstream  oss_shear;       // shear output
};

// -------------------------------------------------------------------------
// Device-resident Lees-Edwards state.
//
// MechSys owns the dem_aux structure and uploads it wholesale, so its layout
// cannot be extended from here.  The Lees-Edwards state therefore travels
// through Domain::pExtraParams, which is a generic void* that MechSys already
// hands to every kernel.  The first member is the part that is re-uploaded
// every step; the virial accumulators come after it so a partial upload never
// clobbers the device-side sums.
// -------------------------------------------------------------------------
struct LeState
{
    double strain;      ///< accumulated shear strain gamma = gamma_dot * t
    double shear_rate;  ///< shear rate gamma_dot (1/T)
    double vZmin;       ///< velocity of the lower z-periodic boundary
    double vZmax;       ///< velocity of the upper z-periodic boundary
};

// LeExtra -- the complete device-side block reached through Domain::pExtraParams.
//   Layout matters: the re-uploaded LeState prefix comes first, so UploadLeState()
//   can copy just sizeof(LeState) and leave the accumulate-only virial sums below it
//   untouched.  Sized once in main() before any kernel runs.
struct LeExtra
{
    LeState state;      ///< LE state, re-uploaded every step
    double  contact_xx; ///< contact virial accumulators (device-side sums)
    double  contact_yy;
    double  contact_zz;
    double  contact_xz;
};

// Host shadow of the device LE state. It is the single source of truth for the
// host-side wrapping helpers and for the per-step upload.
static LeState g_le_state = {0.0, 0.0, 0.0, 0.0};

// SetLeState -- update the host shadow of the Lees-Edwards state.  The device copy
//   is pushed separately by UploadLeState(), once per step, so that the four fields
//   always change together and cannot drift apart between call sites.
static void SetLeState(double strain, double shear_rate,
                       double vZmin, double vZmax)
{
    g_le_state.strain     = strain;
    g_le_state.shear_rate = shear_rate;
    g_le_state.vZmin      = vZmin;
    g_le_state.vZmax      = vZmax;
}

// =========================================================================
// Device kernels (namespace DEM)
// -------------------------------------------------------------------------
// Everything here is either installed into MechSys through Domain's function
// pointer table (pForceVV, pReset) or launched directly by this file's own
// drivers (WrapZ).  Nothing the library launches by name (ResetMaxD, MaxD,
// Translate, Rotate, ...) is touched -- see the header note
// "How this file couples to MechSys" for why that is not possible.
// =========================================================================
namespace DEM {

// BranchVecLE -- Lees-Edwards branch vector (minimum image) between two points.
//   Branch = V1 - V0, wrapped in y, then in z recording the integer image nz,
//   with the x component shifted by -nz*strain*Per.z to undo the Lees-Edwards
//   skew of that z-image, and finally wrapped in x.  nz is returned so the
//   caller can build the matching image velocity.
//   Why local: MechSys' device BranchVec() (dem.cuh) takes only Per and cannot
//   express the shear offset.  This is a __device__ inline helper, not a kernel,
//   so it is simply a new name in namespace DEM -- it replaces nothing.
__device__ inline void BranchVecLE(real3 const & V0, real3 const & V1, real3 & Branch,
                                   real3 const & Per, real strain, int & nz)
{
    Branch = V1 - V0;
    nz = 0;

    if (fabs(Per.y) > 0.0) Branch.y -= round(Branch.y / Per.y) * Per.y;
    if (fabs(Per.z) > 0.0)
    {
        real zround = round(Branch.z / Per.z);
        nz = (int)zround;
        Branch.z -= zround * Per.z;
        Branch.x -= zround * strain * Per.z;
    }
    if (fabs(Per.x) > 0.0) Branch.x -= round(Branch.x / Per.x) * Per.x;
}

// ImageVelocityLE -- velocity of the periodic z-image selected by BranchVecLE:
//   it slides in x at shear_rate*Lz, and in z at the relative velocity of the two
//   z-boundaries (vZmax - vZmin), which is non-zero only while the box height is
//   changing.  It is added to particle 2's velocity when forming the relative
//   velocity, so a contact across the sliding z-boundary sees the correct slip.
//   Why local: MechSys has no concept of a sliding periodic image.
__device__ inline real3 ImageVelocityLE(int nz, dem_aux const * demaux,
                                        LeExtra const * le)
{
    real Lz = demaux[0].Zmax - demaux[0].Zmin;
    return make_real3(nz * le[0].state.shear_rate * Lz,
                      0.0,
                      nz * (le[0].state.vZmax - le[0].state.vZmin));
}

// CalcForceVV_LE -- sphere-sphere contact force, linear spring-dashpot law,
//   Lees-Edwards aware.  Installed as Domain::pForceVV for the "linear" contact
//   law (.inp line 2 = "linear").
//   Differences from MechSys' CalcForceVV:
//     * the branch vector comes from BranchVecLE (sheared minimum image);
//     * the relative velocity includes ImageVelocityLE, so a contact across the
//       sliding z-boundary sees the correct slip velocity;
//     * the contact virial F (x) branch is accumulated in pExtraParams, which is
//       what ShearReport() turns into pressure and shear stress.
//   Why this CAN be overloaded: Domain exposes pForceVV and the library launches
//   it through that pointer from both integrator branches.
__global__ void CalcForceVV_LE(InteractonCU const * Int, ComInteractonCU * CInt,
                               DynInteractonCU * DIntVV, ParticleCU * Par,
                               DynParticleCU * DPar, dem_aux const * demaux,
                               void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nvvint) return;
    LeExtra * diagnostics =
        static_cast<LeExtra *>(extraparams);
    LeExtra const * le = diagnostics;   // LE state shares this device block
    size_t id = DIntVV[ic].Idx;
    size_t i1 = CInt[id].I1;
    size_t i2 = CInt[id].I2;
    real   r1 = Par[i1].R;
    real   r2 = Par[i2].R;
    real3  xi = DPar[i1].x;
    real3  xf = DPar[i2].x;
    real3  Branch;
    real3  image_vel = make_real3(0.0, 0.0, 0.0);

    if (Int[id].BothFree)
    {
        int nz = 0;
        BranchVecLE(xf, xi, Branch, demaux[0].Per, le[0].state.strain, nz);
        image_vel = ImageVelocityLE(nz, demaux, le);
    }
    else Branch = xi - xf;

    real dist  = norm(Branch);
    real delta = r1 + r2 - dist;

    DIntVV[ic].Fn = make_real3(0.0, 0.0, 0.0);
    DIntVV[ic].F = make_real3(0.0, 0.0, 0.0);
    DIntVV[ic].PrevDelta = delta;
    DIntVV[ic].InContact = delta > 0.0;
    // Retain the custom law's existing tangential/rolling state during a
    // temporary gap; inactive candidates contribute no force or overlap row.

    if (delta > 0.0)
    {
        RecordContactOverlap(CInt[id],delta);
        real3 n = -1.0 * Branch / dist;
        real  d = (r1 * r1 - r2 * r2 + dist * dist) / (2.0 * dist);
        real3 x1c = xi + d * n;
        real3 x2c = xf - (dist - d) * n;

        real3 t1, t2, x1, x2;
        Rotation(DPar[i1].w, DPar[i1].Q, t1);
        Rotation(DPar[i2].w, DPar[i2].Q, t2);
        x1 = x1c - xi;
        x2 = x2c - xf;
        real3 vrel = (DPar[i1].v + cross(t1, x1)) -
                     (DPar[i2].v + image_vel + cross(t2, x2));
        real3 vt = vrel - dotreal3(n, vrel) * n;

        real3 Fn_elastic = Int[id].Kn * delta * n;
        real3 Fn_dashpot = Int[id].Gn * dotreal3(n, vrel) * n;
        real3 Fn_total = Fn_elastic + Fn_dashpot;

        if (dotreal3(Fn_total, n) < 0.0)
        {
            Fn_total = make_real3(0.0, 0.0, 0.0);
            Fn_dashpot = -1.0 * Fn_elastic;
        }

        DIntVV[ic].Fn = Fn_total;

        DIntVV[ic].Ft = DIntVV[ic].Ft + (Int[id].Kt * demaux[0].dt) * vt;
        DIntVV[ic].Ft = DIntVV[ic].Ft - dotreal3(DIntVV[ic].Ft, n) * n;

        real3 Ft_elastic = DIntVV[ic].Ft;
        real3 Ft_dashpot = Int[id].Gt * vt;
        real3 Ft_total_tentative = Ft_elastic + Ft_dashpot;

        real friction_limit = Int[id].Mu * norm(Fn_total);
        real3 tan = make_real3(0.0, 0.0, 0.0);

        if (norm(Ft_total_tentative) > friction_limit)
        {
            tan = Ft_total_tentative / norm(Ft_total_tentative);
            Ft_elastic = friction_limit * tan;
            DIntVV[ic].Ft = Ft_elastic;
            Ft_dashpot = make_real3(0.0, 0.0, 0.0);
        }

        real3 Ft_total = Ft_elastic + Ft_dashpot;

        real3 vr = r1 * r2 * cross((t1 - t2), n) / (r1 + r2);
        DIntVV[ic].Fr = DIntVV[ic].Fr + (Int[id].Beta * Int[id].Kt * demaux[0].dt) * vr;
        DIntVV[ic].Fr = DIntVV[ic].Fr - dotreal3(DIntVV[ic].Fr, n) * n;

        tan = DIntVV[ic].Fr;
        if (norm(tan) > 0.0) tan = tan / norm(tan);
        if (norm(DIntVV[ic].Fr) > Int[id].Eta * Int[id].Mu * norm(DIntVV[ic].Fn))
        {
            DIntVV[ic].Fr = Int[id].Eta * Int[id].Mu * norm(DIntVV[ic].Fn) * tan;
        }

        DIntVV[ic].F = Fn_total + Ft_total;

        if (diagnostics != NULL)
        {
            atomicAdd(&diagnostics[0].contact_xx,
                      static_cast<double>((-DIntVV[ic].F.x) * Branch.x));
            atomicAdd(&diagnostics[0].contact_yy,
                      static_cast<double>((-DIntVV[ic].F.y) * Branch.y));
            atomicAdd(&diagnostics[0].contact_zz,
                      static_cast<double>((-DIntVV[ic].F.z) * Branch.z));
            atomicAdd(&diagnostics[0].contact_xz,
                      static_cast<double>((-DIntVV[ic].F.x) * Branch.z));
        }

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

        atomicAdd(&CInt[id].Fnnet.x, Fn_elastic.x);
        atomicAdd(&CInt[id].Fnnet.y, Fn_elastic.y);
        atomicAdd(&CInt[id].Fnnet.z, Fn_elastic.z);
        atomicAdd(&CInt[id].Ftnet.x, Ft_elastic.x);
        atomicAdd(&CInt[id].Ftnet.y, Ft_elastic.y);
        atomicAdd(&CInt[id].Ftnet.z, Ft_elastic.z);

        atomicAdd(&CInt[id].Fndpot.x, Fn_dashpot.x);
        atomicAdd(&CInt[id].Fndpot.y, Fn_dashpot.y);
        atomicAdd(&CInt[id].Fndpot.z, Fn_dashpot.z);
        atomicAdd(&CInt[id].Ftdpot.x, Ft_dashpot.x);
        atomicAdd(&CInt[id].Ftdpot.y, Ft_dashpot.y);
        atomicAdd(&CInt[id].Ftdpot.z, Ft_dashpot.z);

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
    }
}

// CalcForceVV_Hertz_LE -- same Lees-Edwards kinematics as CalcForceVV_LE, but
//   with the Hertz-Mindlin nonlinear stiffness (normal force Kn*delta^1.5, and a
//   sqrt(delta) scaling of both the tangential stiffness and the dashpot terms).
//   Selected when .inp line 2 is "hertz"; installed as Domain::pForceVV.
//   Why this CAN be overloaded: same pForceVV table entry as above, so the two
//   laws are simply two kernels chosen at startup.
__global__ void CalcForceVV_Hertz_LE(InteractonCU const * Int, ComInteractonCU * CInt,
                                     DynInteractonCU * DIntVV, ParticleCU * Par,
                                     DynParticleCU * DPar, dem_aux const * demaux,
                                     void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nvvint) return;
    LeExtra * diagnostics =
        static_cast<LeExtra *>(extraparams);
    LeExtra const * le = diagnostics;   // LE state shares this device block

    size_t id = DIntVV[ic].Idx;
    size_t i1 = CInt[id].I1;
    size_t i2 = CInt[id].I2;
    real   r1 = Par[i1].R;
    real   r2 = Par[i2].R;
    real3  xi = DPar[i1].x;
    real3  xf = DPar[i2].x;
    real3  Branch;
    real3  image_vel = make_real3(0.0, 0.0, 0.0);

    if (Int[id].BothFree)
    {
        int nz = 0;
        BranchVecLE(xf, xi, Branch, demaux[0].Per, le[0].state.strain, nz);
        image_vel = ImageVelocityLE(nz, demaux, le);
    }
    else Branch = xi - xf;

    real dist  = norm(Branch);
    real delta = r1 + r2 - dist;

    DIntVV[ic].Fn = make_real3(0.0, 0.0, 0.0);
    DIntVV[ic].F = make_real3(0.0, 0.0, 0.0);
    DIntVV[ic].PrevDelta = delta;
    DIntVV[ic].InContact = delta > 0.0;
    // Retain the custom law's existing tangential/rolling state during a
    // temporary gap; inactive candidates contribute no force or overlap row.

    if (delta > 0.0)
    {
        RecordContactOverlap(CInt[id],delta);
        real3 n = -1.0 * Branch / dist;
        real  d = (r1 * r1 - r2 * r2 + dist * dist) / (2.0 * dist);
        real3 x1c = xi + d * n;
        real3 x2c = xf - (dist - d) * n;

        real3 t1, t2, x1, x2;
        Rotation(DPar[i1].w, DPar[i1].Q, t1);
        Rotation(DPar[i2].w, DPar[i2].Q, t2);
        x1 = x1c - xi;
        x2 = x2c - xf;
        real3 vrel = (DPar[i1].v + cross(t1, x1)) -
                     (DPar[i2].v + image_vel + cross(t2, x2));
        real3 vt = vrel - dotreal3(n, vrel) * n;

        real sqrtdelta = sqrt(delta);
        DIntVV[ic].Fn = Int[id].Kn * sqrtdelta * delta * n;
        DIntVV[ic].Ft = DIntVV[ic].Ft + (Int[id].Kt * sqrtdelta * demaux[0].dt) * vt;
        DIntVV[ic].Ft = DIntVV[ic].Ft - dotreal3(DIntVV[ic].Ft, n) * n;

        real3 tan = DIntVV[ic].Ft;
        if (norm(tan) > 0.0) tan = tan / norm(tan);
        if (norm(DIntVV[ic].Ft) > Int[id].Mu * norm(DIntVV[ic].Fn))
        {
            DIntVV[ic].Ft = Int[id].Mu * norm(DIntVV[ic].Fn) * tan;
        }

        real3 vr = r1 * r2 * cross((t1 - t2), n) / (r1 + r2);
        DIntVV[ic].Fr = DIntVV[ic].Fr + (Int[id].Beta * Int[id].Kt * sqrtdelta * demaux[0].dt) * vr;
        DIntVV[ic].Fr = DIntVV[ic].Fr - dotreal3(DIntVV[ic].Fr, n) * n;

        tan = DIntVV[ic].Fr;
        if (norm(tan) > 0.0) tan = tan / norm(tan);
        if (norm(DIntVV[ic].Fr) > Int[id].Eta * Int[id].Mu * norm(DIntVV[ic].Fn))
        {
            DIntVV[ic].Fr = Int[id].Eta * Int[id].Mu * norm(DIntVV[ic].Fn) * tan;
        }

        real3 Fn_dashpot = Int[id].Gn * sqrt(sqrtdelta) * dotreal3(n, vrel) * n;
        real3 Ft_dashpot = Int[id].Gt * sqrt(sqrtdelta) * vt;
        DIntVV[ic].F = DIntVV[ic].Fn + DIntVV[ic].Ft + Fn_dashpot + Ft_dashpot;
        atomicAdd(&CInt[id].Fndpot.x, Fn_dashpot.x);
        atomicAdd(&CInt[id].Fndpot.y, Fn_dashpot.y);
        atomicAdd(&CInt[id].Fndpot.z, Fn_dashpot.z);
        atomicAdd(&CInt[id].Ftdpot.x, Ft_dashpot.x);
        atomicAdd(&CInt[id].Ftdpot.y, Ft_dashpot.y);
        atomicAdd(&CInt[id].Ftdpot.z, Ft_dashpot.z);

        if (diagnostics != NULL)
        {
            atomicAdd(&diagnostics[0].contact_xx,
                      static_cast<double>((-DIntVV[ic].F.x) * Branch.x));
            atomicAdd(&diagnostics[0].contact_yy,
                      static_cast<double>((-DIntVV[ic].F.y) * Branch.y));
            atomicAdd(&diagnostics[0].contact_zz,
                      static_cast<double>((-DIntVV[ic].F.z) * Branch.z));
            atomicAdd(&diagnostics[0].contact_xz,
                      static_cast<double>((-DIntVV[ic].F.x) * Branch.z));
        }

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

        atomicAdd(&CInt[id].Fnnet.x, DIntVV[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntVV[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntVV[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntVV[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntVV[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntVV[ic].Ft.z);

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
    }
}

// ResetLE -- installed as Domain::pReset, called once per step before the forces.
//   1) thread 0 zeroes the four contact-virial accumulators held in pExtraParams;
//   2) seeds the force/torque containers from the fixed (Ff, Tf) sources and
//      clears the per-step energy accumulators -- the same duties as the library
//      Reset in dem.cuh.
//   The Lees-Edwards state is deliberately NOT touched: it is re-uploaded once per
//   step by the host (UploadLeState) and must survive the reset.
//   Why this CAN be overloaded: pReset is a table entry and the library calls it
//   through the pointer.
__global__ void ResetLE(ParticleCU * Par, DynParticleCU * DPar,
                                            InteractonCU const * Int, ComInteractonCU * CInt,
                                            dem_aux const * demaux, void * extraparams)
{
    LeExtra * diagnostics =
        static_cast<LeExtra *>(extraparams);
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;

    if (ic == 0 && diagnostics != NULL)
    {
        diagnostics[0].contact_xx = 0.0;
        diagnostics[0].contact_yy = 0.0;
        diagnostics[0].contact_zz = 0.0;
        diagnostics[0].contact_xz = 0.0;
    }

    if (ic < demaux[0].nparts)
    {
        DPar[ic].F = Par[ic].Ff;
        DPar[ic].Flbm = Par[ic].Flbmf;
        Par[ic].T = Par[ic].Tf;
    }
    else if (ic < demaux[0].nparts + demaux[0].ncoint)
    {
        size_t id = ic - demaux[0].nparts;
        CInt[id].Fnnet = Int[id].Fnf;
        CInt[id].Ftnet = Int[id].Ftf;
        CInt[id].Fndpot = make_real3(0.0, 0.0, 0.0);
        CInt[id].Ftdpot = make_real3(0.0, 0.0, 0.0);
        CInt[id].Fther = make_real3(0.0, 0.0, 0.0);
        CInt[id].Nc = 0;
        CInt[id].NormalOverlap = 0.0;
        CInt[id].Epot = 0.0;
        CInt[id].dEvis = 0.0;
        CInt[id].dEfric = 0.0;
    }
}

// WrapZ -- particle wrapping for a sheared (Lees-Edwards) z-periodic cell.
//   A particle leaving through z is folded back by -/+Lz, its x is shifted by
//   -/+strain*Lz (the sliding image) and its x-velocity by -/+shear_rate*Lz, so
//   that the wrap relabels images instead of kicking the particle.  While the box
//   is being compressed (vZmin/vZmax non-zero) the crossing also receives the
//   relative velocity of the two boundaries.
//   The plain x/y wrap follows the LE z shift. Current/reference vertices move
//   together, preserving accumulated displacement for the pre-force skin check.
__global__ void WrapZ(DynParticleCU * DPar, ParticleCU const * Par,
                      real3 * Verts, real3 * Vertso,
                      dem_aux const * demaux, LeExtra const * le)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= demaux[0].nparts) return;

    bool isfree = ((!Par[ic].vxf && !Par[ic].vyf && !Par[ic].vzf &&
                    !Par[ic].wxf && !Par[ic].wyf && !Par[ic].wzf) || Par[ic].FixFree);
    if (!isfree) return;

    real Lz = demaux[0].Zmax - demaux[0].Zmin;
    if (Lz <= 0.0) return;

    real3 dis = make_real3(0.0, 0.0, 0.0);
    if (demaux[0].pz)
    {
        // A constant image-velocity shift commutes with the final VV kick.
        while (DPar[ic].x.z + dis.z < demaux[0].Zmin) {
            dis.z += Lz;
            dis.x += le[0].state.strain * Lz;
            DPar[ic].v.x += le[0].state.shear_rate * Lz;
            DPar[ic].v.z += le[0].state.vZmax - le[0].state.vZmin;
        }
        while (DPar[ic].x.z + dis.z >= demaux[0].Zmax) {
            dis.z -= Lz;
            dis.x -= le[0].state.strain * Lz;
            DPar[ic].v.x -= le[0].state.shear_rate * Lz;
            DPar[ic].v.z -= le[0].state.vZmax - le[0].state.vZmin;
        }
    }
    if (demaux[0].px) {
        while (DPar[ic].x.x + dis.x < demaux[0].Xmin) dis.x += demaux[0].Per.x;
        while (DPar[ic].x.x + dis.x >= demaux[0].Xmax) dis.x -= demaux[0].Per.x;
    }
    if (demaux[0].py) {
        while (DPar[ic].x.y + dis.y < demaux[0].Ymin) dis.y += demaux[0].Per.y;
        while (DPar[ic].x.y + dis.y >= demaux[0].Ymax) dis.y -= demaux[0].Per.y;
    }
    DPar[ic].x = DPar[ic].x + dis;
    DPar[ic].xb = DPar[ic].xb + dis;
    for (size_t iv=Par[ic].Nvi;iv<Par[ic].Nvf;iv++) {
        Verts[iv] = Verts[iv] + dis;
        Vertso[iv] = Vertso[iv] + dis;
    }
}
} // namespace DEM

// SolidVolume -- total solid volume of the particles, used to convert a target
//   solid fraction into a target box height.  Pure host bookkeeping; no library
//   counterpart is needed or replaced.
static double SolidVolume(DEM::Domain const & dom)
{
    double volume = 0.0;
    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        volume += dom.Particles[i]->Props.V;
    }
    return volume;
}

// SolidFraction -- solid volume divided by the current periodic cell volume.
//   This is the constant-volume control variable: it is what Phase 2 drives to
//   the target and what Phase 3 holds fixed.
static double SolidFraction(DEM::Domain const & dom)
{
    double Lx = dom.Xmax - dom.Xmin;
    double Ly = dom.Ymax - dom.Ymin;
    double Lz = dom.Zmax - dom.Zmin;
    return SolidVolume(dom) / (Lx * Ly * Lz);
}

// ReducedLeesEdwardsStrain -- fold the accumulated strain into the fundamental
//   interval [-Lx/(2*Lz), +Lx/(2*Lz)].  Adding Lx/Lz to the strain is the same as
//   translating every z-image by one x-period, so the cell state is identical;
//   the folded value is therefore (a) the reported diagnostic and (b) the correct
//   strain to use when BINNING particles for the contact rebuild, where only the
//   sideways image position matters (see CollectLePairsCells).
static double ReducedLeesEdwardsStrain(double strain, double Lx, double Lz)
{
    if (Lx <= 0.0 || Lz <= 0.0) return strain;
    double period = Lx / Lz;
    if (period <= 0.0) return strain;
    double reduced = fmod(strain + 0.5 * period, period);
    if (reduced < 0.0) reduced += period;
    return reduced - 0.5 * period;
}

// InitializeProductionVelocities -- production start condition.
//   "affine" : vx = shear_rate*(z - Lz/2), i.e. the affine field with zero mean
//              flow and symmetric about the cell centre (a pure Galilean boost of
//              shear_rate*(z - Zmin)); vy = vz = 0.  The Lees-Edwards image
//              kinematics and every contact relative velocity are unaffected.
//   "random" : isotropic Gaussian per component, centre-of-mass velocity removed;
//              used only to study flow build-up from rest.
//   Local by nature: MechSys has no initial-condition model for a sheared cell.
static void InitializeProductionVelocities(DEM::Domain & dom,
                                           std::string const & mode,
                                           double shear_rate,
                                           double random_stddev,
                                           size_t seed)
{
    if (mode == "affine")
    {
        for (size_t i = 0; i < dom.Particles.Size(); ++i)
        {
            DEM::Particle * particle = dom.Particles[i];
            // Symmetric profile about the mid-plane: the affine field has zero
            // mass-weighted mean flow, so the exported velocity field is
            // symmetric with respect to the centre of the cell. This is a pure
            // Galilean boost of shear_rate*(z - Zmin) and leaves the Lees-Edwards
            // image kinematics and all contact relative velocities unchanged.
            particle->v(0) = shear_rate * (particle->x(2) - 0.5 * (dom.Zmin + dom.Zmax));
            particle->v(1) = 0.0;
            particle->v(2) = 0.0;
        }
        return;
    }

    if (mode != "random")
        throw new Fatal("initial_velocity_mode must be 'affine' or 'random', got <%s>",
                        mode.c_str());
    if (!(random_stddev > 0.0) || !std::isfinite(random_stddev))
        throw new Fatal("initial_random_velocity_stddev must be finite and positive for random initialization");

    std::mt19937_64 generator(static_cast<unsigned long long>(seed) ^ 0x9E3779B97F4A7C15ULL);
    std::normal_distribution<double> normal(0.0, random_stddev);
    Vec3_t momentum(0.0, 0.0, 0.0);
    double total_mass = 0.0;

    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        DEM::Particle * particle = dom.Particles[i];
        particle->v = Vec3_t(normal(generator), normal(generator), normal(generator));
        momentum += particle->Props.m * particle->v;
        total_mass += particle->Props.m;
    }
    if (!(total_mass > 0.0) || !std::isfinite(total_mass))
        throw new Fatal("cannot remove initial center-of-mass velocity: invalid particle mass");

    Vec3_t center_of_mass_velocity = momentum / total_mass;
    for (size_t i = 0; i < dom.Particles.Size(); ++i)
        dom.Particles[i]->v -= center_of_mass_velocity;
}

// BranchVecLEHost -- host twin of the device BranchVecLE, so that the host-side
//   contact rebuild measures exactly the distances the GPU kernels will use.
//   Kept byte-for-byte consistent with the device version on purpose: any change
//   must be made in both places.
static Vec3_t BranchVecLEHost(Vec3_t const & V0, Vec3_t const & V1,
                              Vec3_t const & Per, double strain)
{
    Vec3_t branch = V1 - V0;

    if (fabs(Per(1)) > 0.0) branch(1) -= round(branch(1) / Per(1)) * Per(1);
    if (fabs(Per(2)) > 0.0)
    {
        double zround = round(branch(2) / Per(2));
        branch(2) -= zround * Per(2);
        branch(0) -= zround * strain * Per(2);
    }
    if (fabs(Per(0)) > 0.0) branch(0) -= round(branch(0) / Per(0)) * Per(0);

    return branch;
}

// LeesEdwardsDistance -- |branch vector| between two particles under the
//   Lees-Edwards minimum image.  This is THE acceptance test for a contact: both
//   the brute-force scan and the linked-cell scan below funnel into it, which is
//   what makes the two produce identical contact sets.
//   Why local: MechSys' Distance(...,Per) is LE-unaware, and the host overload
//   BranchVec(...,Per,strain) that upstream added to distance.h has been removed
//   (it had no callers).
static double LeesEdwardsDistance(DEM::Particle const * p1, DEM::Particle const * p2,
                                  Vec3_t const & Per, double strain)
{
    Vec3_t branch = BranchVecLEHost(p2->x, p1->x, Per, strain);
    return norm(branch);
}

// WrapHostLeesEdwards -- host counterpart of the device WrapZ: fold free
//   particles back into the cell in x, y and z, applying the Lees-Edwards x-shift
//   and x-velocity shift when they cross z, and refresh every particle's
//   displacement reference (Particle::ResetDisplacements, which is what the
//   Verlet rebuild criterion later measures against).
//   Needed because DnLoadDevice / UpLoadDevice only copy state, they do not wrap,
//   and the host contact rebuild must see positions inside the cell.
static void WrapHostLeesEdwards(DEM::Domain & dom, double strain)
{
    double Lx = dom.Xmax - dom.Xmin;
    double Ly = dom.Ymax - dom.Ymin;
    double Lz = dom.Zmax - dom.Zmin;
    if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0) return;

    double shear_shift = strain * Lz;
    double vel_shift   = g_le_state.shear_rate * Lz;

    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        DEM::Particle * p = dom.Particles[i];
        Vec3_t dis = OrthoSys::O;

        if (p->IsFree())
        {
            while (p->x(2) + dis(2) < dom.Zmin)
            {
                dis(2) += Lz;
                dis(0) += shear_shift;
                p->v(0) += vel_shift;
            }
            while (p->x(2) + dis(2) >= dom.Zmax)
            {
                dis(2) -= Lz;
                dis(0) -= shear_shift;
                p->v(0) -= vel_shift;
            }

            while (p->x(0) + dis(0) < dom.Xmin) dis(0) += Lx;
            while (p->x(0) + dis(0) >= dom.Xmax) dis(0) -= Lx;
            while (p->x(1) + dis(1) < dom.Ymin) dis(1) += Ly;
            while (p->x(1) + dis(1) >= dom.Ymax) dis(1) -= Ly;

            if (norm(dis) > 0.0) p->Translate(dis);
        }

        p->ResetDisplacements();
    }
}

// SetPeriodicBox -- write the box bounds into BOTH the host Domain and its
//   dem_aux mirror (Xmin..Zmax, Per, px/py/pz).  Called whenever the box moves,
//   i.e. at phase boundaries and once per step while compressing.
//   Why it must touch dem_aux explicitly: Domain::demaux is what UpLoadDevice
//   copies to the device, and the box is a dem_aux field (nothing about the box
//   lives in the function-pointer table).
static void SetPeriodicBox(DEM::Domain & dom, double Lx, double Ly,
                           double Zmin, double Zmax)
{
    dom.Xmin = -0.5 * Lx;
    dom.Xmax =  0.5 * Lx;
    dom.Ymin = -0.5 * Ly;
    dom.Ymax =  0.5 * Ly;
    dom.Zmin = Zmin;
    dom.Zmax = Zmax;
    dom.Per = Vec3_t(dom.Xmax - dom.Xmin, dom.Ymax - dom.Ymin, dom.Zmax - dom.Zmin);

    dom.demaux.Xmin = dom.Xmin;
    dom.demaux.Xmax = dom.Xmax;
    dom.demaux.Ymin = dom.Ymin;
    dom.demaux.Ymax = dom.Ymax;
    dom.demaux.Zmin = dom.Zmin;
    dom.demaux.Zmax = dom.Zmax;
    dom.demaux.Per.x = dom.Per(0);
    dom.demaux.Per.y = dom.Per(1);
    dom.demaux.Per.z = dom.Per(2);
    dom.demaux.px = true;
    dom.demaux.py = true;
    dom.demaux.pz = true;
}

// UploadAux -- push the dem_aux mirror (box bounds, Per, periodicity flags) to
//   the device.
//   LIFETIME CONSTRAINT: dem_aux lives in a device buffer that Domain::UpLoadDevice
//   allocates (cudaMalloc of sizeof(dem_aux), inside that method).  Before the first
//   UpLoadDevice, pdemaux is uninitialised, so calling this any earlier writes
//   through a garbage pointer and corrupts the CUDA context (it then surfaces much
//   later as e.g. 'cudaErrorInvalidDevice' from a thrust call).  Every call site in
//   this file is after an UpLoadDevice for that reason -- see the NOTE in
//   SolveCompression().
static void UploadAux(DEM::Domain & dom)
{
    cudaMemcpy(dom.pdemaux, &dom.demaux, sizeof(DEM::dem_aux), cudaMemcpyHostToDevice);
}

// UploadLeState -- push the Lees-Edwards state into the head of the pExtraParams
//   block.  Only the LeState prefix is copied, so the device-side contact-virial
//   accumulators that follow it are left untouched (that is why LeExtra is split
//   into a re-uploaded prefix and an accumulate-only tail).
//   Safe at any time, unlike UploadAux: pExtraParams is allocated once in main()
//   and nothing in MechSys frees or reallocates it.
static void UploadLeState(DEM::Domain & dom)
{
    if (dom.pExtraParams == NULL) return;
    cudaMemcpy(dom.pExtraParams, &g_le_state, sizeof(LeState), cudaMemcpyHostToDevice);
}

// UploadAuxLe -- both blocks, for the places that change the box and the LE state
//   together (the compression ramp).  Only valid after the first UpLoadDevice.
static void UploadAuxLe(DEM::Domain & dom)
{
    UploadAux(dom);
    UploadLeState(dom);
}

// =========================================================================
// Lees-Edwards contact rebuild
// -------------------------------------------------------------------------
// The original version rescanned all N(N-1)/2 pairs. That is O(N^2), with a
// poor memory pattern, and it is entered ~10^4 times per production run, so it
// dominates the cost of the shear phase. This version uses the same strategy as
// MechSys' own contact update: a uniform linked-cell grid whose cell size is
// the interaction range, so a particle is only tested against the neighbouring
// cells.
//
// The grid cannot be built in the plain frame. Under Lees-Edwards a z-image of
// a particle is displaced in x by nz*strain*Lz, so two particles whose plain
// x-separation is large can still be in contact. The grid is therefore built in
// the AFFINE frame
//        a = x - gamma_r * z,    gamma_r = ReducedLeesEdwardsStrain(strain)
// which un-shears the image lattice: in (a, y, z) the Lees-Edwards minimum
// image is a plain minimum image. The link between the affine separation and
// the actual Lees-Edwards branch is, writing u = dx - nz*gamma_r*Lz and
// dz_w = dz - nz*Lz (|dz_w| <= Lz/2),
//        a_j - a_i = u - gamma_r*dz_w,
//        wrap_Lx(a_j - a_i) = wrap_Lx( LE_x - gamma_r*dz_w ),
// where LE_x = wrap_Lx(u) is the x component of BranchVecLE. Hence every
// interacting pair (|LE_x| <= cut) satisfies
//        |wrap_Lx(a_j - a_i)| <= cut + |gamma_r|*Lz/2  =:  W .
// Since |gamma_r| <= Lx/(2*Lz) the widening is at most Lx/4, and that is exactly
// the window the x-neighbourhood below covers. The candidate list is therefore
// a strict superset of the true contact set, and every candidate is still
// tested with the exact LeesEdwardsDistance(), so the accepted set is identical
// to the brute-force one.
//
// Accepted pairs are collected and then sorted by (i,j) before any interaction
// is created, so that dom.Interactons is filled in exactly the same order as
// the brute-force scan. That keeps the device-side interaction order -- and
// therefore the floating-point accumulation order -- unchanged.
//
// Define LE_REBUILD_VERIFY to cross-check the cell candidate set against the
// brute-force one on every rebuild.
// =========================================================================
typedef std::pair<size_t, size_t> LePair;

// CollectLePairsBruteForce -- the original all-pairs scan: every i<j tested with
//   LeesEdwardsDistance against Dmax_i + Dmax_j + 2*Alpha, wall-wall pairs skipped.
//   It is the REFERENCE implementation (any candidate generator must reproduce its
//   accepted set) and the runtime FALLBACK when the cell grid cannot be justified.
//   Why it cannot simply be replaced by the library's own contact update:
//   Domain::UpdateContacts()/Domain::UpdateContactsDevice() are library methods
//   (callable, not overloadable) and they are Lees-Edwards unaware -- both the
//   linked-cell binning and CInteracton::UpdateContacts(Alpha,Per,iter) use the
//   plain periodic image, so they would miss every contact that exists only
//   through the sheared z-image.
static void CollectLePairsBruteForce(DEM::Domain & dom, double strain,
                                     std::vector<LePair> & pairs)
{
    for (size_t i = 0; i + 1 < dom.Particles.Size(); ++i)
    {
        bool pi_has_vf = !dom.Particles[i]->IsFree();
        for (size_t j = i + 1; j < dom.Particles.Size(); ++j)
        {
            bool pj_has_vf = !dom.Particles[j]->IsFree();
            if (pi_has_vf && pj_has_vf) continue;

            double cutoff = dom.Particles[i]->Dmax + dom.Particles[j]->Dmax + 2.0 * dom.Alpha;
            if (LeesEdwardsDistance(dom.Particles[i], dom.Particles[j], dom.Per, strain) > cutoff) continue;

            pairs.push_back(LePair(i, j));
        }
    }
}

// LeCellRebuildUsable -- guard for the linked-cell path.  Every bound below is
//   conservative and parameter independent: if any of them fails we fall back to
//   the brute-force scan rather than risk a missed contact.  N < 64 also fails
//   (bookkeeping is not worth it for tiny systems).
//   Set MECHSYS_LE_BRUTE_FORCE=1 to force the brute-force scan.  That makes a
//   bit-for-bit A/B possible from one binary: the two paths produce the same
//   contact list in the same order, so all outputs must be identical.
static bool LeCellRebuildUsable(DEM::Domain const & dom)
{
    {
        static const bool brute_force_forced = (std::getenv("MECHSYS_LE_BRUTE_FORCE") != NULL);
        if (brute_force_forced) return false;
    }

    const size_t N = dom.Particles.Size();
    if (N < 64) return false;                       // not worth the bookkeeping

    const double Lx = dom.Xmax - dom.Xmin;
    const double Ly = dom.Ymax - dom.Ymin;
    const double Lz = dom.Zmax - dom.Zmin;
    if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0) return false;

    double dmax = 0.0;
    for (size_t i = 0; i < N; ++i)
        if (dom.Particles[i]->Dmax > dmax) dmax = dom.Particles[i]->Dmax;
    const double cut = 2.0 * dmax + 2.0 * dom.Alpha;

    // The affine window must not wrap around the x-period: W < Lx/2, and
    // W <= cut + Lx/4 because |gamma_r| <= Lx/(2*Lz).
    if (cut + 0.25 * Lx >= 0.5 * Lx) return false;
    // Only the nz = -1,0,1 z-images can be within range of a pair, so one cell
    // shell in y and z is enough only if a cell (>= cut) is smaller than half
    // the box in those directions.
    if (cut >= 0.5 * Ly || cut >= 0.5 * Lz) return false;
    return true;
}

// CollectLePairsCells -- linked-cell candidate generation (the derivation is in
//   the block comment above).  Builds a uniform grid in the affine frame, walks
//   the (rx,1,1)-shell neighbourhood, and hands every candidate to the same exact
//   LeesEdwardsDistance test the brute-force scan uses, so the accepted set is
//   identical while the number of distance evaluations drops by ~15x.
//   This is the local analogue of MechSys' UpdateLinkedCells(), which cannot be
//   reused because it bins and pairs in the plain (un-sheared) frame.
static void CollectLePairsCells(DEM::Domain & dom, double strain,
                                std::vector<LePair> & pairs)
{
    const size_t N  = dom.Particles.Size();
    const double Lx = dom.Xmax - dom.Xmin;
    const double Ly = dom.Ymax - dom.Ymin;
    const double Lz = dom.Zmax - dom.Zmin;
    const double gamma = ReducedLeesEdwardsStrain(strain, Lx, Lz);

    double dmax = 0.0;
    for (size_t i = 0; i < N; ++i)
        if (dom.Particles[i]->Dmax > dmax) dmax = dom.Particles[i]->Dmax;
    const double cut  = 2.0 * dmax + 2.0 * dom.Alpha;
    const double skew = 0.5 * std::fabs(gamma) * Lz;   // half of the Lees-Edwards widening

    // Flooring the cell count makes every cell slightly larger than the
    // interaction range, which is what the one-shell neighbourhood in y and z
    // requires.
    int nx = int(std::floor(Lx / cut)); if (nx < 1) nx = 1;
    int ny = int(std::floor(Ly / cut)); if (ny < 1) ny = 1;
    int nz = int(std::floor(Lz / cut)); if (nz < 1) nz = 1;
    const double dx = Lx / nx;
    const double dy = Ly / ny;
    const double dz = Lz / nz;
    const int rx = int(std::ceil((cut + skew) / dx));  // x shells for the affine window

    std::vector< std::vector<size_t> > cells(size_t(nx) * size_t(ny) * size_t(nz));
    for (size_t i = 0; i < N; ++i)
    {
        DEM::Particle const * p = dom.Particles[i];
        double ax = p->x(0) - gamma * p->x(2);         // affine x
        ax -= Lx * std::floor((ax - dom.Xmin) / Lx);   // wrapped into one x-period

        int ix = int(std::floor((ax      - dom.Xmin) / dx));
        int iy = int(std::floor((p->x(1) - dom.Ymin) / dy));
        int iz = int(std::floor((p->x(2) - dom.Zmin) / dz));
        if (ix < 0) ix = 0;  if (ix >= nx) ix = nx - 1;
        if (iy < 0) iy = 0;  if (iy >= ny) iy = ny - 1;
        if (iz < 0) iz = 0;  if (iz >= nz) iz = nz - 1;

        cells[size_t(ix) + size_t(nx) * (size_t(iy) + size_t(ny) * size_t(iz))].push_back(i);
    }

    // Exact Lees-Edwards test, identical to the brute-force one.
    auto consider = [&](size_t a, size_t b)
    {
        if (a > b) { size_t t = a; a = b; b = t; }
        if (!dom.Particles[a]->IsFree() && !dom.Particles[b]->IsFree()) return;
        const double cutoff = dom.Particles[a]->Dmax + dom.Particles[b]->Dmax + 2.0 * dom.Alpha;
        if (LeesEdwardsDistance(dom.Particles[a], dom.Particles[b], dom.Per, strain) > cutoff) return;
        pairs.push_back(LePair(a, b));
    };

    for (size_t idx = 0; idx < cells.size(); ++idx)
    {
        if (cells[idx].empty()) continue;
        const int cx = int(idx % size_t(nx));
        const int cy = int((idx / size_t(nx)) % size_t(ny));
        const int cz = int(idx / (size_t(nx) * size_t(ny)));

        std::unordered_set<size_t> visited;            // a wrapped cell can repeat
        for (int kz = -1; kz <= 1; ++kz)
        for (int ky = -1; ky <= 1; ++ky)
        for (int kx = -rx; kx <= rx; ++kx)
        {
            int jx = (cx + kx) % nx; if (jx < 0) jx += nx;
            int jy = (cy + ky) % ny; if (jy < 0) jy += ny;
            int jz = (cz + kz) % nz; if (jz < 0) jz += nz;
            const size_t jdx = size_t(jx) + size_t(nx) * (size_t(jy) + size_t(ny) * size_t(jz));
            if (!visited.insert(jdx).second) continue;

            if (jdx == idx)                            // internal pairs of this cell
            {
                for (size_t a = 0;     a + 1 < cells[idx].size(); ++a)
                for (size_t b = a + 1; b     < cells[idx].size(); ++b)
                    consider(cells[idx][a], cells[idx][b]);
            }
            else if (jdx > idx)                        // each unordered cell pair once
            {
                for (size_t a = 0; a < cells[idx].size(); ++a)
                for (size_t b = 0; b < cells[jdx].size(); ++b)
                    consider(cells[idx][a], cells[jdx][b]);
            }
        }
    }
}


// ResetContactsLeesEdwardsSpheres -- rebuild the whole contact list for the
//   sheared cell.  This is the function that replaces, for Phase 3, MechSys'
//   Domain::UpdateContacts():
//     * it cannot CALL UpdateContacts(), because that one pairs particles with the
//       plain periodic image and would drop every contact that exists only through
//       the sheared z-image;
//     * it cannot OVERRIDE UpdateContacts() either, because that is a non-virtual
//       library method, absent from Domain's function-pointer table.
//   So it is reimplemented here: wrap the host particles (WrapHostLeesEdwards),
//   collect the contact pairs (cells, or brute force), sort them back into the
//   brute-force (i,j) order, then create/reuse the interactions exactly as the
//   library does (PairtoCInt keyed by HashFunction(i,j), CInteractonSphere for
//   sphere-sphere, per-pair Mu from FricCoeff) and finally drop the pairs that
//   fell out of range.
//   The sort is what makes this a drop-in replacement: dom.Interactons ends up in
//   the same order as the original scan, so the device-side accumulation order --
//   and therefore every result -- is unchanged.
static void ResetContactsLeesEdwardsSpheres(DEM::Domain & dom, double strain)
{
    WrapHostLeesEdwards(dom, strain);
    dom.Interactons.Resize(0);

    // ---- 1. candidate generation (linked cells, or brute force) ----------
    std::vector<LePair> pairs;
    pairs.reserve(dom.Particles.Size() * 8);

    if (LeCellRebuildUsable(dom)) CollectLePairsCells(dom, strain, pairs);
    else                          CollectLePairsBruteForce(dom, strain, pairs);

#ifdef LE_REBUILD_VERIFY
    {
        std::vector<LePair> reference;
        CollectLePairsBruteForce(dom, strain, reference);
        std::vector<LePair> a(pairs), b(reference);
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        if (a != b)
            throw new Fatal("LE rebuild: cell candidates (%zd) differ from brute force (%zd)",
                            a.size(), b.size());
    }
#endif

    // ---- 2. restore the brute-force ordering -----------------------------
    // Sorting by (i,j) makes dom.Interactons exactly what the pair-by-pair scan
    // would have produced, so the device interaction order (and hence the
    // floating-point accumulation order) is unchanged.
    std::sort(pairs.begin(), pairs.end());

    // ---- 3. create / reuse the interactions ------------------------------
    std::unordered_set<size_t> active_pairs;
    active_pairs.reserve(pairs.size() * 2 + 16);

    for (size_t k = 0; k < pairs.size(); ++k)
    {
        const size_t i = pairs[k].first;
        const size_t j = pairs[k].second;
        const size_t hash = HashFunction(i, j);

        if (dom.PairtoCInt.count(hash) == 0)
        {
            if (dom.Particles[i]->Verts.Size() == 1 && dom.Particles[j]->Verts.Size() == 1)
            {
                dom.PairtoCInt[hash] = new DEM::CInteractonSphere(dom.Particles[i], dom.Particles[j], dom.ContactLaw);
            }
            else
            {
                dom.PairtoCInt[hash] = new DEM::CInteracton(dom.Particles[i], dom.Particles[j], dom.ContactLaw);
                dom.PairtoCInt[hash]->UpdateContacts(dom.Alpha, dom.Per, dom.iter);
            }

            std::pair<int, int> p(dom.Particles[i]->Tag, dom.Particles[j]->Tag);
            if (dom.FricCoeff.count(p) == 1) dom.PairtoCInt[hash]->Mu = dom.FricCoeff[p];
        }

        active_pairs.insert(hash);
        dom.Interactons.Push(dom.PairtoCInt[hash]);
    }

    for (auto it = dom.PairtoCInt.begin(); it != dom.PairtoCInt.end();)
    {
        if (active_pairs.count(it->first) == 0)
        {
            delete it->second;
            it = dom.PairtoCInt.erase(it);
        }
        else
        {
            ++it;
        }
    }
}

// UpdateContactsDeviceLeesEdwards -- one complete contact refresh:
//   download the device state to the host (DnLoadDevice), rebuild the list with the
//   Lees-Edwards-aware rebuild above, then upload it again (UpLoadDevice).  Single
//   entry point used by both phase drivers.
static void UpdateContactsDeviceLeesEdwards(DEM::Domain & dom, double strain)
{
    // The custom Hertz kernel stores tangential FORCE. Preserve it directly
    // by pair identity, including temporary gaps where a zero Hertz stiffness
    // makes the library's CPU displacement conversion non-invertible.
    thrust::host_vector<DEM::DynInteractonCU> old = dom.bDynInteractonsVV;
    thrust::host_vector<DEM::ComInteractonCU> common = dom.bComInteractons;
    std::unordered_map<size_t,DEM::DynInteractonCU> histories;
    for (auto const & h : old) {
        auto const & c = common[h.Idx];
        histories.emplace(HashFunction(c.I1,c.I2),h);
    }
    dom.DnLoadDevice(dom.Nproc, true);
    ResetContactsLeesEdwardsSpheres(dom, strain);
    dom.UpLoadDevice(dom.Nproc, true, false);
    thrust::host_vector<DEM::DynInteractonCU> current = dom.bDynInteractonsVV;
    common = dom.bComInteractons;
    for (auto & h : current) {
        auto const & c = common[h.Idx];
        auto const found = histories.find(HashFunction(c.I1,c.I2));
        if (found != histories.end()) {
            h.Ft=found->second.Ft;h.Fr=found->second.Fr;
            h.PrevDelta=found->second.PrevDelta;h.InContact=found->second.InContact;
        }
    }
    thrust::copy(current.begin(),current.end(),dom.bDynInteractonsVV.begin());
    cudaMemset(dom.pMaxDCU,0,dom.demaux.nverts*sizeof(real));
}

// ShiftZOriginToZero -- after compression the cell sits at [Zmin_target, Zmax_target].
//   This translates every particle so that Zmin becomes exactly 0, so that the affine
//   profile used at production start (and the z reported in the outputs) is measured
//   from the cell floor.  It then re-writes the periodic box for the new origin and
//   clears the Lees-Edwards state.  Returns the applied shift.
static double ShiftZOriginToZero(DEM::Domain & dom)
{
    double shift = -dom.Zmin;
    if (fabs(shift) <= 1.0e-12) return 0.0;

    Vec3_t dz(0.0, 0.0, shift);
    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        dom.Particles[i]->Translate(dz);
    }

    double Lz = dom.Zmax - dom.Zmin;
    SetPeriodicBox(dom, dom.Xmax - dom.Xmin, dom.Ymax - dom.Ymin, 0.0, Lz);
    SetLeState(0.0, 0.0, 0.0, 0.0);
    return shift;
}

// =========================================================================
// Phase 2 Setup callback: move z-boundaries inward for compression
// =========================================================================
// CompressionSetup -- per-step hook of Phase 2, called at the top of every step by
//   SolveCompression (it used to be the ptSetup callback passed to Domain::Solve).
//   Duties:
//     * move both z-periodic boundaries linearly from their initial to their target
//       position over comp_time, and publish the applied boundary velocities in the
//       LE state so a particle crossing z picks up the relative velocity of the two
//       moving images (see WrapZ);
//     * rebuild the contact list whenever the box height has changed by more than
//       half the Verlet distance since the last rebuild (this is the compression
//       analogue of the displacement criterion, and in practice it is the one that
//       fires: 89 rebuilds for the reference run);
//     * on the last step, clamp the box to its target and do a final rebuild.
//   Why a plain function and not a table entry: the library's ptSetup mechanism lives
//   inside Domain::Solve, which this file deliberately does not use.
void CompressionSetup(DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    if (!dat.compressing) return;
    if (dom.Finished) return;

    double elapsed = dom.Time - dat.comp_start_time;
    double bottom_velocity = (dat.Zmin_target - dat.Zmin_initial) / dat.comp_duration;
    double top_velocity = (dat.Zmax_target - dat.Zmax_initial) / dat.comp_duration;

    if (elapsed >= dat.comp_duration)
    {
        SetPeriodicBox(dom, dat.Lx, dat.Ly, dat.Zmin_target, dat.Zmax_target);
        SetLeState(0.0, 0.0, 0.0, 0.0);
        UploadAuxLe(dom);
        double Lz = dom.Zmax - dom.Zmin;
        if (fabs(Lz - dat.last_contact_Lz) > 1e-12)
        {
            dom.UpdateContactsDevice();
            dat.last_contact_Lz = Lz;
        }
        return;
    }

    double frac = elapsed / dat.comp_duration;
    double Zmin_cur = dat.Zmin_initial + (dat.Zmin_target - dat.Zmin_initial) * frac;
    double Zmax_cur = dat.Zmax_initial + (dat.Zmax_target - dat.Zmax_initial) * frac;

    SetPeriodicBox(dom, dat.Lx, dat.Ly, Zmin_cur, Zmax_cur);
    SetLeState(0.0, 0.0, bottom_velocity, top_velocity);
    UploadAuxLe(dom);

    double rebuild_dz = 0.5 * dat.base_alpha;
    if (rebuild_dz <= 0.0) rebuild_dz = 1e-12;
    double Lz = dom.Zmax - dom.Zmin;
    if (fabs(Lz - dat.last_contact_Lz) >= rebuild_dz)
    {
        dom.UpdateContactsDevice();
        dat.last_contact_Lz = dom.Zmax - dom.Zmin;
        dat.comp_rebuilds++;
    }
}

// =========================================================================
// Phase 2 Report callback
// =========================================================================
// CompressionReport -- writes <key>_comp.res (time, Lz, solid fraction) once per
//   output interval and closes the file when the phase ends.  Formerly the ptReport
//   callback of Domain::Solve; now called directly by SolveCompression.
void CompressionReport(DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    if (!dat.compressing) return;

    if (!dat.oss_comp.is_open())
    {
        String fs;
        fs.Printf("%s_comp.res", dom.FileKey.CStr());
        dat.oss_comp.open(fs.CStr());
        dat.oss_comp << Util::_10_6 << "Time"
                     << Util::_8s << "Lz"
                     << Util::_8s << "SolidFrac"
                     << Util::_8s << "Phase"
                     << std::endl;
    }

    if (!dom.Finished)
    {
        double Lz = dom.Zmax - dom.Zmin;
        double phi = SolidFraction(dom);

        dat.oss_comp << Util::_10_6 << dom.Time
                     << Util::_8s << Lz
                     << Util::_8s << phi
                     << Util::_8s << "COMPRESS"
                     << std::endl;
    }
    else
    {
        dat.oss_comp.close();
        dat.compressing = false;
    }
}

// =========================================================================
// Phase 3 Setup callback: Lees-Edwards shear
// =========================================================================
// AdvanceShearState -- one step of the shear bookkeeping: accumulate the strain
//   (gamma += shear_rate*dt), fold it into the diagnostic interval, re-assert the
//   periodic box, and publish strain/shear_rate into the LE state.  The boundary
//   velocities stay zero, which is what makes Phase 3 strictly constant volume.
//   The caller does the actual device upload, so the driver controls when the GPU
//   sees the new state.
static void AdvanceShearState(DEM::Domain & dom, UserData & dat)
{
    double dt = dom.Dt;
    double Lz = dom.Zmax - dom.Zmin;

    dat.cum_strain += dat.shear_rate * dt;
    dat.reduced_strain = ReducedLeesEdwardsStrain(dat.cum_strain, dat.Lx, Lz);

    SetPeriodicBox(dom, dat.Lx, dat.Ly, dom.Zmin, dom.Zmax);
    SetLeState(dat.cum_strain, dat.shear_rate, 0.0, 0.0);
}

// =========================================================================
// Phase 3 Report callback
// =========================================================================
// LeesEdwardsHostMeasurements -- host-side observables computed once per output
//   interval: temperature, and the contact and kinetic parts of the pressure and of
//   the shear stress.
struct LeesEdwardsHostMeasurements
{
    double temperature;
    double pressure_contact;
    double pressure_kinetic;
    double pressure_total;
    double shear_stress_contact;
    double shear_stress_kinetic;
    double shear_stress_total;
};

// MeasureLeesEdwardsState -- computes the observables reported by ShearReport.
//   Temperature and the kinetic stress come from the fluctuation field, i.e. the
//   particle velocity minus the affine field minus the mass-weighted mean
//   fluctuation.  The contact stress comes from the virial accumulators that
//   CalcForceVV_LE / CalcForceVV_Hertz_LE atomicAdd-ed into pExtraParams during the
//   last step.  The volume divisor is the constant cell volume.
static LeesEdwardsHostMeasurements MeasureLeesEdwardsState(DEM::Domain & dom,
                                                            UserData const & dat)
{
    LeesEdwardsHostMeasurements output = {};
    double volume = (dom.Xmax - dom.Xmin) * (dom.Ymax - dom.Ymin) *
                    (dom.Zmax - dom.Zmin);
    if (volume <= 0.0 || dom.Particles.Size() == 0) return output;

    double total_mass = 0.0;
    Vec3_t mean_fluctuation(0.0, 0.0, 0.0);
    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        DEM::Particle const * particle = dom.Particles[i];
        double affine_vx = dat.shear_rate * (particle->x(2) - dom.Zmin);
        Vec3_t fluctuation = particle->v;
        fluctuation(0) -= affine_vx;
        mean_fluctuation += particle->Props.m * fluctuation;
        total_mass += particle->Props.m;
    }
    if (total_mass > 0.0) mean_fluctuation /= total_mass;

    double kinetic_trace = 0.0;
    double kinetic_xz = 0.0;
    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        DEM::Particle const * particle = dom.Particles[i];
        double affine_vx = dat.shear_rate * (particle->x(2) - dom.Zmin);
        Vec3_t fluctuation = particle->v;
        fluctuation(0) -= affine_vx;
        fluctuation -= mean_fluctuation;
        kinetic_trace += particle->Props.m *
                         (fluctuation(0) * fluctuation(0) +
                          fluctuation(1) * fluctuation(1) +
                          fluctuation(2) * fluctuation(2));
        kinetic_xz += particle->Props.m * fluctuation(0) * fluctuation(2);
    }

    LeExtra diagnostics = {};
    if (dom.pExtraParams != NULL)
        cudaMemcpy(&diagnostics, dom.pExtraParams, sizeof(diagnostics), cudaMemcpyDeviceToHost);

    output.temperature = total_mass > 0.0 ? kinetic_trace / (3.0 * total_mass) : 0.0;
    output.pressure_contact =
        (diagnostics.contact_xx + diagnostics.contact_yy + diagnostics.contact_zz) /
        (3.0 * volume);
    output.pressure_kinetic = kinetic_trace / (3.0 * volume);
    output.pressure_total = output.pressure_contact + output.pressure_kinetic;
    output.shear_stress_contact = -diagnostics.contact_xz / volume;
    output.shear_stress_kinetic = -kinetic_xz / volume;
    output.shear_stress_total = output.shear_stress_contact + output.shear_stress_kinetic;
    return output;
}

// ShearReport -- writes <key>_shear.res once per output interval: time, strain,
//   reduced strain, Lz, solid fraction, temperature, and the total/contact/kinetic
//   parts of the pressure and of the shear stress.  Closes the file at the end of
//   the phase.
void ShearReport(DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    if (!dat.shearing) return;

    if (!dat.oss_shear.is_open())
    {
        String fs;
        fs.Printf("%s_shear.res", dom.FileKey.CStr());
        dat.oss_shear.open(fs.CStr());
        dat.oss_shear << std::setprecision(17)
                      << "time strain reduced_LE_strain Lz phi temperature "
                      << "pressure pressure_contact pressure_kinetic "
                      << "shear_stress shear_stress_contact shear_stress_kinetic"
                      << std::endl;
    }

    if (!dom.Finished)
    {
        double Lz = dom.Zmax - dom.Zmin;
        LeesEdwardsHostMeasurements measured = MeasureLeesEdwardsState(dom, dat);
        dat.oss_shear << dom.Time << " "
                      << dat.cum_strain << " "
                      << dat.reduced_strain << " "
                      << Lz << " "
                      << SolidFraction(dom) << " "
                      << measured.temperature << " "
                      << measured.pressure_total << " "
                      << measured.pressure_contact << " "
                      << measured.pressure_kinetic << " "
                      << measured.shear_stress_total << " "
                      << measured.shear_stress_contact << " "
                      << measured.shear_stress_kinetic
                      << std::endl;
    }
    else
    {
        dat.oss_shear.close();
        dat.shearing = false;
    }
}

#ifdef USE_CUDA
// PrepareShearSolve -- one-off setup shared by both phase drivers.  It reproduces the
//   part of Domain::Solve that runs BEFORE the step loop:
//     FileKey / idx_out / Dt / Nproc, Domain::Initialize (which allocates the device
//     buffers and uploads dem_aux, and therefore creates pdemaux), the free /
//     non-free partition with Vs, Ms and MaxDmax, zeroed A and Wdot arrays, one
//     MtData per OpenMP thread, the host Per, and Finished = false.
//   It deliberately does NOT touch the device: UpLoadDevice is left to the caller,
//   because the two phases initialise the acceleration state differently.
static void PrepareShearSolve(DEM::Domain & dom, double dt,
                              char const * filekey, size_t Nproc)
{
    dom.FileKey.Printf("%s", filekey);
    dom.idx_out = 0;
    dom.Dt = dt;
    dom.Nproc = Nproc;
    dom.Initialize(dt);

    dom.FreePar.Resize(0);
    dom.NoFreePar.Resize(0);
    dom.Vs = 0.0;
    dom.Ms = 0.0;
    dom.MaxDmax = 0.0;

    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        if (dom.Particles[i]->IsFree())
        {
            dom.Vs += dom.Particles[i]->Props.V;
            dom.Ms += dom.Particles[i]->Props.m;
            if (dom.Particles[i]->Dmax > dom.MaxDmax) dom.MaxDmax = dom.Particles[i]->Dmax;
            dom.FreePar.Push(i);
        }
        else
        {
            dom.NoFreePar.Push(i);
        }
        dom.Particles[i]->Index = i;
    }

    dom.A.Resize(dom.Particles.Size());
    dom.Wdot.Resize(dom.Particles.Size());
    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        dom.A[i] = OrthoSys::O;
        dom.Wdot[i] = OrthoSys::O;
    }

    dom.MTD = new DEM::MtData[Nproc];
    for (size_t i = 0; i < Nproc; ++i)
    {
        dom.MTD[i].N_Proc = Nproc;
        dom.MTD[i].ProcRank = i;
        dom.MTD[i].Dom = &dom;
        dom.MTD[i].Dmx = 0.0;
    }

    dom.Per = Vec3_t(dom.Xmax - dom.Xmin, dom.Ymax - dom.Ymin, dom.Zmax - dom.Zmin);
    dom.Finished = false;
}

// Sample initial forces at zero duration. Sampling must not consume a history
// step; list rebuilds never call this function and preserve A/Wdot instead.
static void RefreshDeviceAccelerations(DEM::Domain & dom)
{
    auto const history = dom.bDynInteractonsVV;
    double const dt = dom.demaux.dt;
    dom.demaux.dt = 0.0;
    UploadAux(dom);
    size_t nt = dom.Nthread;

    dom.pReset<<<(dom.demaux.nparts + dom.demaux.ncoint) / nt + 1, nt>>>
        (dom.pParticlesCU, dom.pDynParticlesCU, dom.pInteractons,
         dom.pComInteractons, dom.pdemaux, dom.pExtraParams);

    if (dom.demaux.nvvint > 0)
        dom.pForceVV<<<(dom.demaux.nvvint + nt - 1) / nt, nt>>>
            (dom.pInteractons, dom.pComInteractons, dom.pDynInteractonsVV,
             dom.pParticlesCU, dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);
    if (dom.demaux.neeint > 0)
        dom.pForceEE<<<(dom.demaux.neeint + nt - 1) / nt, nt>>>
            (dom.pEdgesCU, dom.pVertsCU, dom.pInteractons, dom.pComInteractons,
             dom.pDynInteractonsEE, dom.pParticlesCU, dom.pDynParticlesCU,
             dom.pdemaux, dom.pExtraParams);
    if (dom.demaux.nvfint > 0)
        dom.pForceVF<<<(dom.demaux.nvfint + nt - 1) / nt, nt>>>
            (dom.pFacesCU, dom.pFacidCU, dom.pVertsCU, dom.pInteractons,
             dom.pComInteractons, dom.pDynInteractonsVF, dom.pParticlesCU,
             dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);
    if (dom.demaux.nfvint > 0)
        dom.pForceFV<<<(dom.demaux.nfvint + nt - 1) / nt, nt>>>
            (dom.pFacesCU, dom.pFacidCU, dom.pVertsCU, dom.pInteractons,
             dom.pComInteractons, dom.pDynInteractonsFV, dom.pParticlesCU,
             dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);

    dom.pRecomputeAccelerations<<<(dom.demaux.nparts + nt - 1) / nt, nt>>>
        (dom.pParticlesCU, dom.pDynParticlesCU, dom.pA, dom.pWdot, dom.pdemaux);
    dom.bDynInteractonsVV = history;
    dom.pDynInteractonsVV = thrust::raw_pointer_cast(dom.bDynInteractonsVV.data());
    dom.demaux.dt = dt;
    UploadAux(dom);
}

// Rebuild only the candidate list, preserving post-drift state and histories.
static void RebuildLeesEdwardsContacts(DEM::Domain & dom, UserData & dat)
{
    UpdateContactsDeviceLeesEdwards(dom, dat.cum_strain);
    dat.le_contact_strain = dat.cum_strain;
    dat.le_rebuilds++;
}

// WriteShearFrame -- write one XDMF + BF frame pair for the current FileKey and bump
//   idx_out.  Used by both phases.  Written here because the library's own frame
//   writing happens inside Domain::Solve, which this file does not use.
static void WriteShearFrame(DEM::Domain & dom)
{
    String fn;
    String fnbf;
    fn.Printf("%s_%04d", dom.FileKey.CStr(), int(dom.idx_out));
    fnbf.Printf("%s_bf_%04d", dom.FileKey.CStr(), int(dom.idx_out));
    dom.WriteXDMF(fn.CStr());
    dom.WriteBF(fnbf.CStr());
    dom.idx_out++;
}

// One physical VV step. Particle motion and image motion share the 2*Alpha
// pair allowance: a list is valid only while 2*max_particle_motion + image_motion
// is smaller than that allowance. Test post-drift geometry before applying force.
static bool AdvanceLeesEdwardsStep(DEM::Domain & dom, UserData & dat)
{
    AdvanceShearState(dom,dat);
    UploadAuxLe(dom);
    size_t nt = dom.Nthread;
    dom.pVerletStep1<<<(dom.demaux.nparts+nt-1)/nt,nt>>>
        (dom.pVertsCU,dom.pParticlesCU,dom.pDynParticlesCU,dom.pA,dom.pdemaux);
    if (dom.RotPar) dom.pOrientationUpdate<<<(dom.demaux.nparts+nt-1)/nt,nt>>>
        (dom.pVertsCU,dom.pParticlesCU,dom.pDynParticlesCU,dom.pWdot,dom.pdemaux);
    DEM::WrapZ<<<(dom.demaux.nparts+nt-1)/nt,nt>>>
        (dom.pDynParticlesCU,dom.pParticlesCU,dom.pVertsCU,dom.pVertsoCU,
         dom.pdemaux,static_cast<LeExtra const *>(dom.pExtraParams));
    size_t const blocks = dom.demaux.nverts/nt+1;
    DEM::MaxD<<<blocks,nt>>>(dom.pVertsCU,dom.pVertsoCU,dom.pMaxDCU,dom.pdemaux);
    real const maxdis = *thrust::max_element(dom.bMaxDCU.begin(),dom.bMaxDCU.begin()+blocks);
    double const image_motion = fabs(dat.cum_strain-dat.le_contact_strain)*(dom.Zmax-dom.Zmin);
    bool const rebuilt = 2.0*maxdis + image_motion >= 2.0*dom.Alpha;
    if (rebuilt) RebuildLeesEdwardsContacts(dom,dat);
    dom.pReset<<<(dom.demaux.nparts+dom.demaux.ncoint)/nt+1,nt>>>
        (dom.pParticlesCU,dom.pDynParticlesCU,dom.pInteractons,
         dom.pComInteractons,dom.pdemaux,dom.pExtraParams);
    if (dom.demaux.nvvint > 0)
        dom.pForceVV<<<(dom.demaux.nvvint + nt - 1) / nt, nt>>>
            (dom.pInteractons, dom.pComInteractons, dom.pDynInteractonsVV,
             dom.pParticlesCU, dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);
    if (dom.demaux.neeint > 0)
        dom.pForceEE<<<(dom.demaux.neeint + nt - 1) / nt, nt>>>
            (dom.pEdgesCU, dom.pVertsCU, dom.pInteractons, dom.pComInteractons,
             dom.pDynInteractonsEE, dom.pParticlesCU, dom.pDynParticlesCU,
             dom.pdemaux, dom.pExtraParams);
    if (dom.demaux.nvfint > 0)
        dom.pForceVF<<<(dom.demaux.nvfint + nt - 1) / nt, nt>>>
            (dom.pFacesCU, dom.pFacidCU, dom.pVertsCU, dom.pInteractons,
             dom.pComInteractons, dom.pDynInteractonsVF, dom.pParticlesCU,
             dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);
    if (dom.demaux.nfvint > 0)
        dom.pForceFV<<<(dom.demaux.nfvint + nt - 1) / nt, nt>>>
            (dom.pFacesCU, dom.pFacidCU, dom.pVertsCU, dom.pInteractons,
             dom.pComInteractons, dom.pDynInteractonsFV, dom.pParticlesCU,
             dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);

    dom.pFinalizeVelocity<<<(dom.demaux.nparts + nt - 1) / nt, nt>>>
        (dom.pParticlesCU, dom.pDynParticlesCU, dom.pA, dom.pdemaux);
    if (dom.RotPar) dom.pFinalizeRotation<<<(dom.demaux.nparts + nt - 1) / nt, nt>>>
        (dom.pParticlesCU, dom.pDynParticlesCU, dom.pWdot, dom.pdemaux);


    dom.demaux.iter++;
    dom.demaux.Time += dom.Dt;
    dom.Time += dom.Dt;
    dom.iter++;
    return rebuilt;
}

// SolveLeesEdwardsShear -- Phase 3 driver (constant-volume simple shear).
//   Written here rather than using Domain::Solve because the step loop must launch
//   the local WrapZ, and the library's own wrap is a hard-wired launch with no
//   hook and no extraparams (see the header note, item 3).
//   The shared step helper checks post-drift particle/image motion before force
//   evaluation. Rebuilds migrate contact history without integrating it again.
static void SolveLeesEdwardsShear(DEM::Domain & dom, double shear_time,
                                  double dt, double dtOut, UserData & dat,
                                  bool Render, size_t Nproc)
{
    (void)Render;
    double tf = dom.Time + shear_time;
    double tout = dom.Time;
    size_t iter_b = dom.iter;
    size_t iter_t = 0;
    size_t numup = 0;

    dom.UseVelocityVerlet = true;
    PrepareShearSolve(dom, dt, "shearing", Nproc);
    SetPeriodicBox(dom, dat.Lx, dat.Ly, dom.Zmin, dom.Zmax);
    SetLeState(dat.cum_strain, dat.shear_rate, 0.0, 0.0);
    UploadLeState(dom);

    ResetContactsLeesEdwardsSpheres(dom, dat.cum_strain);
    dom.UpLoadDevice(Nproc, true, true);
    RefreshDeviceAccelerations(dom);
    dat.le_contact_strain = dat.cum_strain;
    dat.le_rebuilds++;

    while (dom.Time < tf)
    {
        if (dom.Time >= tout)
        {
            dom.DnLoadDevice(Nproc, true);
            ShearReport(dom, &dat);
            WriteShearFrame(dom);
            tout += dtOut;
        }

        if (AdvanceLeesEdwardsStep(dom,dat)) {
            numup++;
            iter_t += dom.iter-iter_b;
            iter_b = dom.iter;
        }
    }

    dom.Finished = true;
    dom.DnLoadDevice(Nproc, true);
    ShearReport(dom, &dat);

    printf("  Total shear iterations: %zd\n", dom.iter);
    if (numup > 0)
    {
        printf("  Average shear iterations between LE contact updates: %zd\n",
               iter_t / numup);
    }
}
// =========================================================================
// Phase 2 driver: prescribed-rate compression with moving z-boundaries
// -------------------------------------------------------------------------
// Local replacement for Domain::Solve.  The library wraps particles with a
// hard-wired plain periodic launch that this file cannot replace or parameterise,
// so the CUDA driver is reproduced here with the local WrapZ kernel launched
// directly.
//
// Why the loop is owned at all: the wrap has to happen inside the step loop, the
// library offers no hook for it, and giving it one would force the Lees-Edwards
// state into dem_aux.  Owning the loop keeps that state in
// pExtraParams and leaves MechSys untouched; the price is that this loop has to be
// kept in step with Domain::Solve by hand.
//
// It mirrors Domain::Solve for the configuration this solver uses
// (USE_CUDA, legacy position Verlet during compression):
//   * per-step order: setup callback -> dtOut output -> pReset ->
//     pForceVV/EE/VF/FV -> pTranslate -> pRotate -> wrap -> ResetMaxD -> MaxD,
//   * rebuild when maxdis > Alpha, via UpdateContactsDevice(),
//   * counters: demaux.iter/demaux.Time advance with the step; Time/iter advance
//     at the end of the step,
//   * same callbacks and same XDMF/BF frame naming.
// Known deviation: MechSys computes the initial accelerations on the host via
// Interacton::CalcForce (the CPU contact law). Here they are produced on the
// device by RefreshDeviceAccelerations, i.e. by the same GPU kernels used for
// the rest of the run. This is the only respect in which a compression started
// from this file can differ from one started from the MechSys driver; it keeps
// the whole integration on a single contact implementation. Compare Phase 2
// outputs before switching production runs over.
// =========================================================================
static void SolveCompression(DEM::Domain & dom, double comp_time,
                             double dt, double dtOut, UserData & dat,
                             bool Render, size_t Nproc)
{
    (void)Render;
    double tf     = dom.Time + comp_time;
    double tout   = dom.Time;
    size_t iter_b = dom.iter;
    size_t iter_t = 0;
    size_t numup  = 0;

    dom.UseVelocityVerlet = false;
    PrepareShearSolve(dom, dt, "compression", Nproc);
    SetPeriodicBox(dom, dat.Lx, dat.Ly, dom.Zmin, dom.Zmax);

    // NOTE: dom.pdemaux is allocated inside Domain::UpLoadDevice, so dem_aux must
    // NOT be uploaded before that call. The Lees-Edwards state lives in
    // pExtraParams (allocated in main), so it can be pushed here safely.
    SetLeState(0.0, 0.0, 0.0, 0.0);
    UploadLeState(dom);

    dom.UpLoadDevice(Nproc, true, true);
    RefreshDeviceAccelerations(dom);

    while (dom.Time < tf)
    {
        // Advance the z-boundaries for this step (and rebuild contacts if the
        // boundary moved by more than half the Verlet distance).
        CompressionSetup(dom, &dat);

        if (dom.Time >= tout)
        {
            dom.DnLoadDevice(Nproc, true);
            CompressionReport(dom, &dat);
            WriteShearFrame(dom);
            tout += dtOut;
        }

        size_t nt = dom.Nthread;
        dom.pReset<<<(dom.demaux.nparts + dom.demaux.ncoint) / nt + 1, nt>>>
            (dom.pParticlesCU, dom.pDynParticlesCU, dom.pInteractons,
             dom.pComInteractons, dom.pdemaux, dom.pExtraParams);

        // Position-Verlet + leapfrog path (UseVelocityVerlet == false), matching
        // the MechSys default integrator for this solver.
        if (dom.demaux.nvvint > 0)
            dom.pForceVV<<<(dom.demaux.nvvint + nt - 1) / nt, nt>>>
                (dom.pInteractons, dom.pComInteractons, dom.pDynInteractonsVV,
                 dom.pParticlesCU, dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);
        if (dom.demaux.neeint > 0)
            dom.pForceEE<<<(dom.demaux.neeint + nt - 1) / nt, nt>>>
                (dom.pEdgesCU, dom.pVertsCU, dom.pInteractons, dom.pComInteractons,
                 dom.pDynInteractonsEE, dom.pParticlesCU, dom.pDynParticlesCU,
                 dom.pdemaux, dom.pExtraParams);
        if (dom.demaux.nvfint > 0)
            dom.pForceVF<<<(dom.demaux.nvfint + nt - 1) / nt, nt>>>
                (dom.pFacesCU, dom.pFacidCU, dom.pVertsCU, dom.pInteractons,
                 dom.pComInteractons, dom.pDynInteractonsVF, dom.pParticlesCU,
                 dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);
        if (dom.demaux.nfvint > 0)
            dom.pForceFV<<<(dom.demaux.nfvint + nt - 1) / nt, nt>>>
                (dom.pFacesCU, dom.pFacidCU, dom.pVertsCU, dom.pInteractons,
                 dom.pComInteractons, dom.pDynInteractonsFV, dom.pParticlesCU,
                 dom.pDynParticlesCU, dom.pdemaux, dom.pExtraParams);

        dom.pTranslate<<<dom.demaux.nparts / nt + 1, nt>>>
            (dom.pVertsCU, dom.pParticlesCU, dom.pDynParticlesCU,
             dom.pdemaux, dom.pExtraParams);
        if (dom.RotPar)
            dom.pRotate<<<dom.demaux.nparts / nt + 1, nt>>>
                (dom.pVertsCU, dom.pParticlesCU, dom.pDynParticlesCU,
                 dom.pdemaux, dom.pExtraParams);

        // Lees-Edwards / moving-boundary z-wrap, launched directly from here.
        DEM::WrapZ<<<(dom.demaux.nparts + nt - 1) / nt, nt>>>
            (dom.pDynParticlesCU, dom.pParticlesCU, dom.pVertsCU, dom.pVertsoCU,
             dom.pdemaux, static_cast<LeExtra const *>(dom.pExtraParams));
        DEM::ResetMaxD<<<dom.demaux.nparts / nt + 1, nt>>>
            (dom.pVertsCU, dom.pVertsoCU, dom.pMaxDCU,
             dom.pParticlesCU, dom.pDynParticlesCU, dom.pdemaux);

        DEM::MaxD<<<dom.demaux.nverts / nt + 1, nt>>>
            (dom.pVertsCU, dom.pVertsoCU, dom.pMaxDCU, dom.pdemaux);
        thrust::device_vector<real>::iterator it =
            thrust::max_element(dom.bMaxDCU.begin(), dom.bMaxDCU.end());
        real maxdis = *it;

        dom.demaux.iter++;
        dom.demaux.Time += dt;

        if (maxdis > dom.Alpha)
        {
            numup++;
            iter_t += dom.iter - iter_b;
            iter_b = dom.iter;
            dom.UpdateContactsDevice();
        }

        dom.Time += dt;
        dom.iter++;
    }

    dom.Finished = true;
    dom.DnLoadDevice(Nproc, true);
    CompressionReport(dom, &dat);

    printf("  Total compression iterations: %zd\n", dom.iter);
    if (numup > 0)
    {
        printf("  Average compression iterations between contact updates: %zd\n",
               iter_t / numup);
    }
}



#endif

// =========================================================================
// Main
// =========================================================================
// main -- reads <filekey>.inp, runs the three phases and writes the outputs.
//   .inp fields, in order: Verlet distance; contact law ("linear" or "hertz");
//   render flag; Kn Kt Gn Gt Mu (rolling stiffness) Beta (rolling resistance) Eta;
//   eps; particle radius; RNG seed; dt; output interval; Lx; Ly; Lz of the
//   generation box; material density; fill fraction; target solid fraction;
//   compression duration; shear rate; production duration; then optionally
//   initial_velocity_mode ("affine"/"random"), the random-mode velocity stddev, and
//   the minimum acceptable compressed height in diameters.
//   Only three things are installed into MechSys (see the header note):
//   pForceVV, pReset, and pExtraParams -- which is sized for LeExtra here and must
//   be allocated before any kernel runs.
int main(int argc, char **argv) try
{
    if (argc < 2)
        throw new Fatal("Usage: %s filekey [Nproc]\n", argv[0]);

    size_t Nproc = (argc >= 3) ? atoi(argv[2]) : 1;

    String filekey(argv[1]);
    String filename = filekey + ".inp";
    if (!Util::FileExists(filename))
        throw new Fatal("File <%s> not found", filename.CStr());
    std::ifstream infile(filename.CStr());

    // ---- Read parameters ----
    double verlet;
    String contactlaw;
    size_t RenderVideo;
    double Kn, Kt, Gn, Gt, Mu0, Beta, Eta;
    double Eps, R;
    size_t seed;
    double dt, dtOut;
    double Lx, Ly, Lz_gen;
    double rho;
    double fraction;
    double target_phi;
    double comp_time;
    double shear_rate;
    double shear_time;
    std::string initial_velocity_mode = "affine";
    double initial_random_velocity_stddev = 0.01;
    double minimum_compressed_height_over_d = 1.0;

    {
        infile >> verlet;       infile.ignore(200, '\n');
        infile >> contactlaw;   infile.ignore(200, '\n');
        infile >> RenderVideo;  infile.ignore(200, '\n');
        infile >> Kn;           infile.ignore(200, '\n');
        infile >> Kt;           infile.ignore(200, '\n');
        infile >> Gn;           infile.ignore(200, '\n');
        infile >> Gt;           infile.ignore(200, '\n');
        infile >> Mu0;          infile.ignore(200, '\n');
        infile >> Beta;         infile.ignore(200, '\n');
        infile >> Eta;          infile.ignore(200, '\n');
        infile >> Eps;          infile.ignore(200, '\n');
        infile >> R;            infile.ignore(200, '\n');
        infile >> seed;         infile.ignore(200, '\n');
        infile >> dt;           infile.ignore(200, '\n');
        infile >> dtOut;        infile.ignore(200, '\n');
        infile >> Lx;           infile.ignore(200, '\n');
        infile >> Ly;           infile.ignore(200, '\n');
        infile >> Lz_gen;       infile.ignore(200, '\n');
        infile >> rho;          infile.ignore(200, '\n');
        infile >> fraction;     infile.ignore(200, '\n');
        infile >> target_phi;   infile.ignore(200, '\n');
        infile >> comp_time;    infile.ignore(200, '\n');
        infile >> shear_rate;   infile.ignore(200, '\n');
        infile >> shear_time;   infile.ignore(200, '\n');
        if (infile >> initial_velocity_mode)
        {
            infile.ignore(200, '\n');
            if (infile >> initial_random_velocity_stddev)
            {
                infile.ignore(200, '\n');
                if (infile >> minimum_compressed_height_over_d)
                    infile.ignore(200, '\n');
            }
        }
    }

    printf("=== Input parameters ===\n");
    printf("  R = %g  (diameter = %g)\n", R, 2.0*R);
    printf("  Box: %g x %g x %g (generation)\n", Lx, Ly, Lz_gen);
    printf("  Fill fraction: %g\n", fraction);
    printf("  Target solid fraction: %g\n", target_phi);
    printf("  Comp time: %g, Shear rate: %g, Shear time: %g\n", comp_time, shear_rate, shear_time);
    printf("  Initial production velocity mode: %s\n", initial_velocity_mode.c_str());
    if (initial_velocity_mode == "random")
        printf("  Initial random velocity component stddev: %g\n", initial_random_velocity_stddev);
    printf("  Minimum compressed height: %g particle diameters\n",
           minimum_compressed_height_over_d);

    if (initial_velocity_mode != "affine" && initial_velocity_mode != "random")
        throw new Fatal("initial_velocity_mode must be 'affine' or 'random', got <%s>",
                        initial_velocity_mode.c_str());
    if (!(minimum_compressed_height_over_d > 0.0)
        || !std::isfinite(minimum_compressed_height_over_d))
        throw new Fatal("minimum_compressed_height_over_d must be finite and positive");

    UserData dat;
    dat.Lx            = Lx;
    dat.Ly            = Ly;
    dat.Lz_initial    = Lz_gen;
    dat.Lz_target     = Lz_gen;
    dat.Zmin_initial  = 0.0;
    dat.Zmin_target   = 0.0;
    dat.Zmax_initial  = Lz_gen;
    dat.Zmax_target   = Lz_gen;
    dat.last_contact_Lz = Lz_gen;
    dat.comp_start_time = 0.0;
    dat.comp_duration  = comp_time;
    dat.cum_strain    = 0.0;
    dat.reduced_strain = 0.0;
    dat.shear_rate    = shear_rate;
    dat.base_alpha    = verlet;
    dat.comp_rebuilds = 0;
    dat.le_contact_strain = 0.0;
    dat.le_rebuilds   = 0;
    dat.compressing   = false;
    dat.shearing      = false;
    dat.dt            = dt;
    dat.dtOut         = dtOut;

    size_t cl = (contactlaw == "hertz") ? 1 : 0;

    // ---- Domain ----
    DEM::Domain dom(&dat, cl);
    dom.Alpha  = verlet;
    dom.Dilate = true;
    // Overload only pre-existing MechSys function pointers: the Lees-Edwards
    // contact branch vector, the LE force laws and the LE reset all live in
    // this file (see CalcForceVV_LE / CalcForceVV_Hertz_LE / ResetLE).
    dom.pForceVV = (cl == 1) ? DEM::CalcForceVV_Hertz_LE : DEM::CalcForceVV_LE;
    dom.pReset   = DEM::ResetLE;

    // =====================================================================
    // Phase 1: Generate particles in periodic domain (no gravity)
    // =====================================================================
    printf("\n=== Phase 1: Particle Generation ===\n");

    SetPeriodicBox(dom, Lx, Ly, dat.Zmin_initial, dat.Zmax_initial);

    Vec3_t Xmin(dom.Xmin, dom.Ymin, dom.Zmin);
    Vec3_t Xmax(dom.Xmax, dom.Ymax, dom.Zmax);
    dom.GenSpheresBox(-1, Xmin, Xmax, R, rho, "HCP", seed, fraction, Eps);

    if (dom.Particles.Size() == 0)
        throw new Fatal("Particle generation produced no particles; check Lx/Ly/Lz/R/fraction");

    printf("  Particles generated: %zd\n", dom.Particles.Size());
    printf("  Periodic box: x=[%g,%g], y=[%g,%g], z=[%g,%g]\n",
           dom.Xmin, dom.Xmax, dom.Ymin, dom.Ymax, dom.Zmin, dom.Zmax);

    Dict B1;
    B1.Set(-1, "Kn Kt Gn Gt Mu Beta Eta", Kn, Kt, Gn, Gt, Mu0, Beta, Eta);
    dom.SetProps(B1);

    LeExtra le_extra = {};
    cudaMalloc(&dom.pExtraParams, sizeof(LeExtra));
    cudaMemcpy(dom.pExtraParams, &le_extra, sizeof(le_extra), cudaMemcpyHostToDevice);

    SetPeriodicBox(dom, Lx, Ly, dat.Zmin_initial, dat.Zmax_initial);
    SetLeState(0.0, 0.0, 0.0, 0.0);

    printf("Phase 1 complete.\n\n");

    // =====================================================================
    // Phase 2: Compression via contact wrapping with moving boundaries
    // =====================================================================
    printf("=== Phase 2: Compression via contact wrapping ===\n");

    if (target_phi <= 0.0)
        throw new Fatal("target_phi must be positive");

    double V_solid = SolidVolume(dom);
    double A_xy = Lx * Ly;
    double Lz_current = dom.Zmax - dom.Zmin;  // = Lz_gen after Phase 1
    double phi_current = V_solid / (A_xy * Lz_current);
    double Lz_target = V_solid / (target_phi * A_xy);
    double d_max = 0.0;
    for (size_t i = 0; i < dom.Particles.Size(); ++i)
    {
        double diameter = 2.0 * dom.Particles[i]->Dmax;
        if (diameter > d_max) d_max = diameter;
    }
    double planned_final_Lz =
        (target_phi > phi_current && comp_time > 0.0) ? Lz_target : Lz_current;
    if (!(d_max > 0.0) || !std::isfinite(d_max))
        throw new Fatal("cannot validate compressed height: invalid particle diameter");
    if (planned_final_Lz + 1.0e-12 < minimum_compressed_height_over_d * d_max)
        throw new Fatal(
            "planned final Lz (%g) is below the required %g particle diameters (%g)",
            planned_final_Lz, minimum_compressed_height_over_d,
            minimum_compressed_height_over_d * d_max);
    double z_center = 0.5 * (dom.Zmin + dom.Zmax);
    double Zmin_target = z_center - 0.5 * Lz_target;
    double Zmax_target = z_center + 0.5 * Lz_target;

    printf("  Number of particles: %zd\n", dom.Particles.Size());
    printf("  Solid volume: %g\n", V_solid);
    printf("  Current box Lz: %g\n", Lz_current);
    printf("  Current solid fraction: %g\n", phi_current);
    printf("  Target solid fraction:  %g\n", target_phi);
    printf("  Target box Lz: %g\n", Lz_target);
    printf("  dLz: %g (%.1f%% compression)\n",
           Lz_current - Lz_target, 100.0*(Lz_current - Lz_target)/Lz_current);
    printf("  Lower boundary: %g -> %g\n", dom.Zmin, Zmin_target);
    printf("  Upper boundary: %g -> %g\n", dom.Zmax, Zmax_target);
    printf("  Post-compression origin shift: final Zmin -> 0\n");

    if (target_phi > phi_current && comp_time > 0.0)
    {
        dat.compressing   = true;
        dat.comp_start_time = dom.Time;
        dat.Lz_initial    = Lz_current;
        dat.Lz_target     = Lz_target;
        dat.Zmin_initial  = dom.Zmin;
        dat.Zmin_target   = Zmin_target;
        dat.Zmax_initial  = dom.Zmax;
        dat.Zmax_target   = Zmax_target;
        dat.last_contact_Lz = Lz_current;
        dat.comp_rebuilds = 0;

        printf("  Starting two-sided periodic-boundary compression...\n");
        printf("  Verlet alpha: %g, boundary rebuild dz: %g\n",
               dom.Alpha, 0.5 * dat.base_alpha);

        SolveCompression(dom, comp_time, dt, dtOut, dat, RenderVideo, Nproc);

        SetPeriodicBox(dom, Lx, Ly, dat.Zmin_target, dat.Zmax_target);
        SetLeState(0.0, 0.0, 0.0, 0.0);
        UploadLeState(dom);

        double z_shift = ShiftZOriginToZero(dom);

        dat.Lz_initial    = dom.Zmax - dom.Zmin;
        dat.Lz_target     = dat.Lz_initial;
        dat.Zmin_initial  = dom.Zmin;
        dat.Zmin_target   = dom.Zmin;
        dat.Zmax_initial  = dom.Zmax;
        dat.Zmax_target   = dom.Zmax;
        dat.last_contact_Lz = dat.Lz_initial;

        // Verify final state
        double final_Lz = dom.Zmax - dom.Zmin;
        double final_phi = V_solid / (A_xy * final_Lz);
        printf("  Shifted compressed system by dz = %.6f\n", z_shift);
        printf("  Final boundaries: z=[%.6f, %.6f]\n", dom.Zmin, dom.Zmax);
        printf("  Final Lz: %.6f, final phi: %.6f\n", final_Lz, final_phi);
        printf("  Boundary-driven compression rebuilds: %zd\n", dat.comp_rebuilds);

        if (dat.oss_comp.is_open()) dat.oss_comp.close();
        dat.compressing = false;
        printf("Phase 2 complete.\n\n");
    }
    else
    {
        printf("  Warning: target <= current or comp_time <= 0. Skipping compression.\n");
        printf("Phase 2 complete.\n\n");
    }

    // The random mode starts without an imposed velocity gradient so the
    // Lees-Edwards image kinematics establish the shear flow dynamically.
    InitializeProductionVelocities(dom, initial_velocity_mode, shear_rate,
                                   initial_random_velocity_stddev, seed);

    // =====================================================================
    // Phase 3: Lees-Edwards shear
    // =====================================================================
    printf("=== Phase 3: Lees-Edwards Shear ===\n");

    if (shear_time > 0.0 && fabs(shear_rate) > 1e-12)
    {
        SetPeriodicBox(dom, Lx, Ly, dom.Zmin, dom.Zmax);
        SetLeState(0.0, shear_rate, 0.0, 0.0);
        UploadLeState(dom);

        dat.shearing  = true;
        dat.cum_strain = 0.0;
        dat.reduced_strain = 0.0;
        dat.le_contact_strain = 0.0;
        dat.le_rebuilds = 0;

        printf("  Shear rate: %g\n", shear_rate);
        printf("  Shear duration: %g\n", shear_time);
        printf("  Box height Lz: %g\n", dom.Zmax - dom.Zmin);
        printf("  Verlet alpha: %g\n", dom.Alpha);
        printf("  Starting Lees-Edwards shear...\n");

        SolveLeesEdwardsShear(dom, shear_time, dt, dtOut,
                              dat, RenderVideo, Nproc);

        printf("  Shear phase complete. Final strain: %g\n", dat.cum_strain);
        printf("  Lees-Edwards contact rebuilds: %zd\n", dat.le_rebuilds);

        if (dat.oss_shear.is_open()) dat.oss_shear.close();
        dat.shearing = false;
        printf("Phase 3 complete.\n\n");
    }
    else
    {
        printf("  Skipping shear phase.\n");
    }

    printf("All phases completed successfully.\n");

    cudaFree(dom.pExtraParams);

    return 0;
}

MECHSYS_CATCH
