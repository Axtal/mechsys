/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// GPU (CUDA) version of tdem/TX_PID_periodic.cpp
//
// Same physics as the CPU driver: a PID servo that drives the *periodic cell*
// of a dompargen-generated sphere packing to a target confining pressure and
// then shears it at a constant axial strain rate, with velocity-Verlet
// integration and Coulomb mode 2.
//
// What runs on the GPU:
//   * contact detection lists (uploaded once per Verlet rebuild),
//   * all contact forces, both integration steps and the maximum-displacement
//     reduction (this is the stock DEM CUDA path, domain.h#ifdef USE_CUDA),
//   * the virial stress tensor, computed every step by a reduction kernel over
//     the device contact list (CalcStressVV below) so that no contact data is
//     ever transferred back to the host.
//
// The host only does the closed loop of the servo: it reads 6 numbers per step
// (the stress accumulator), runs the PID and writes the new periodic cell back
// to the device with a one-thread kernel (SetCellBox).

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/linalg/matvec.h>

#include <algorithm>
#include <cmath>
#include <fstream>

using std::cout;


//////////////////////////////////////////////////////////////////////////
// PID controller (identical formulation to tdem/PID_v2.cpp)
//////////////////////////////////////////////////////////////////////////

struct PIDController
{
    double integral;
    double previousError;
    bool   initialized;
};

double PID(PIDController & controller, double targetStress, double measuredStress,
           double pressureScale, double dt, double kp, double ki, double kd)
{
    const double maxRate = 0.10;
    const double error   = (targetStress - measuredStress) / pressureScale;
    if (!controller.initialized)
    {
        controller.previousError = error;
        controller.initialized   = true;
    }

    controller.integral += 0.5 * (error + controller.previousError) * dt;
    if (ki > 0.0)
    {
        const double integralLimit = maxRate / ki;
        controller.integral = std::max(
            -integralLimit, std::min(controller.integral, integralLimit));
    }

    const double derivative = (error - controller.previousError) / dt;
    controller.previousError = error;

    const double rate = kp * error + ki * controller.integral + kd * derivative;
    return std::max(-maxRate, std::min(rate, maxRate));
}


//////////////////////////////////////////////////////////////////////////
// User data shared with the Setup/Report callbacks
//////////////////////////////////////////////////////////////////////////

struct UserData
{
    bool          renderVideo;     ///< Write XDMF output ?
    double        dt;              ///< Time step
    double        pressure;        ///< Isotropic target pressure (positive)
    double        axialStrainRate; ///< Axial strain rate during shear (positive compresses)
    double        pidKp;           ///< Proportional gain
    double        pidKi;           ///< Integral gain
    double        pidKd;           ///< Derivative gain
    double        verletDist;      ///< Verlet distance (cell-change rebuild threshold)
    bool          shearStage;      ///< Is the deviatoric stage running ?
    Vec3_t        initialLength;   ///< Reference cell lengths
    Vec3_t        measured;        ///< Last measured (compressive positive) stress
    Vec3_t        lastVerletLength;///< Cell lengths at the last Verlet list rebuild
    bool          verletLengthSet; ///< Bookkeeping flag for lastVerletLength
    PIDController controller[3];   ///< One PID controller per axis
    std::ofstream stressStrain;    ///< Stress/strain output file
    real        * devSigma;        ///< Device buffer with the 6 stress components
};


//////////////////////////////////////////////////////////////////////////
// Geometric helpers operating on the periodic cell (host side)
//////////////////////////////////////////////////////////////////////////

double BoxLength(DEM::Domain const & dom, size_t axis)
{
    if (axis == 0) return dom.Xmax - dom.Xmin;
    if (axis == 1) return dom.Ymax - dom.Ymin;
    return dom.Zmax - dom.Zmin;
}

double BoxVolume(DEM::Domain const & dom)
{
    return BoxLength(dom,0) * BoxLength(dom,1) * BoxLength(dom,2);
}


//////////////////////////////////////////////////////////////////////////
// GPU kernels
//////////////////////////////////////////////////////////////////////////

namespace DEM
{

/// Virial stress over the active vertex-vertex (sphere) contact list.
/// For every active contact F is the total contact force added to particle i2
/// and Branch = x(i1) - x(i2) (minimum image), exactly as in CalcForceVV, so
/// the sign convention matches the CPU driver (compression is positive).
__global__ void CalcStressVV(InteractonCU const * Int, ComInteractonCU const * CInt,
                             DynInteractonCU const * DIntVV, DynParticleCU const * DPar,
                             dem_aux const * aux, real * sigma)
{
    __shared__ real s[6];
    const size_t tid = threadIdx.x;
    if (tid < 6) s[tid] = 0.0;
    __syncthreads();

    const size_t n    = aux[0].nvvint;
    const size_t step = (size_t)gridDim.x * blockDim.x;
    for (size_t ic = (size_t)blockIdx.x*blockDim.x + tid; ic < n; ic += step)
    {
        const size_t id = DIntVV[ic].Idx;
        const size_t i1 = CInt[id].I1;
        const size_t i2 = CInt[id].I2;

        real3 Branch;
        if (Int[id].BothFree) BranchVec(DPar[i2].x, DPar[i1].x, Branch, aux[0].Per);
        else                  Branch = DPar[i1].x - DPar[i2].x;

        const real3 F = make_real3(CInt[id].Fnnet.x + CInt[id].Ftnet.x +
                                   CInt[id].Fndpot.x + CInt[id].Ftdpot.x + CInt[id].Fther.x,
                                   CInt[id].Fnnet.y + CInt[id].Ftnet.y +
                                   CInt[id].Fndpot.y + CInt[id].Ftdpot.y + CInt[id].Fther.y,
                                   CInt[id].Fnnet.z + CInt[id].Ftnet.z +
                                   CInt[id].Fndpot.z + CInt[id].Ftdpot.z + CInt[id].Fther.z);

        atomicAdd(&s[0], F.x*Branch.x);
        atomicAdd(&s[1], F.y*Branch.y);
        atomicAdd(&s[2], F.z*Branch.z);
        atomicAdd(&s[3], F.x*Branch.y);
        atomicAdd(&s[4], F.x*Branch.z);
        atomicAdd(&s[5], F.y*Branch.z);
    }
    __syncthreads();
    if (tid < 6) atomicAdd(&sigma[tid], s[tid]);
}

/// One-thread kernel that updates the periodic cell on the device.  Only the
/// cell limits are touched, so the device's own time/iteration counters are
/// preserved.
__global__ void SetCellBox(dem_aux * aux, real xmin, real xmax,
                           real ymin, real ymax, real zmin, real zmax)
{
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
    aux[0].Xmin = xmin;  aux[0].Xmax = xmax;
    aux[0].Ymin = ymin;  aux[0].Ymax = ymax;
    aux[0].Zmin = zmin;  aux[0].Zmax = zmax;
    aux[0].Per  = make_real3(xmax-xmin, ymax-ymin, zmax-zmin);
}

} // namespace DEM


//////////////////////////////////////////////////////////////////////////
// Cell deformation (FCC shearbox style): only the cell limits and Per are
// rescaled about the cell centre; the stress develops through the
// over-the-boundary contacts.  The host copy of the cell is the reference
// that UpLoadDevice uses at every Verlet rebuild, so it is kept in sync.
//////////////////////////////////////////////////////////////////////////
void ScaleCell(DEM::Domain & dom, Vec3_t const & strainRate, double dt)
{
    const Vec3_t centre(0.5*(dom.Xmax + dom.Xmin),
                        0.5*(dom.Ymax + dom.Ymin),
                        0.5*(dom.Zmax + dom.Zmin));

    const double factor[3] =
    {
        1.0 - strainRate(0)*dt,
        1.0 - strainRate(1)*dt,
        1.0 - strainRate(2)*dt
    };

    const double length[3] =
    {
        (dom.Xmax - dom.Xmin) * factor[0],
        (dom.Ymax - dom.Ymin) * factor[1],
        (dom.Zmax - dom.Zmin) * factor[2]
    };

    dom.Xmin = centre(0) - 0.5*length[0];  dom.Xmax = centre(0) + 0.5*length[0];
    dom.Ymin = centre(1) - 0.5*length[1];  dom.Ymax = centre(1) + 0.5*length[1];
    dom.Zmin = centre(2) - 0.5*length[2];  dom.Zmax = centre(2) + 0.5*length[2];
    dom.Per  = Vec3_t(length[0], length[1], length[2]);

    // keep the host mirror of the auxiliary structure consistent
    dom.demaux.Xmin = dom.Xmin;  dom.demaux.Xmax = dom.Xmax;
    dom.demaux.Ymin = dom.Ymin;  dom.demaux.Ymax = dom.Ymax;
    dom.demaux.Zmin = dom.Zmin;  dom.demaux.Zmax = dom.Zmax;
    dom.demaux.Per  = make_real3(length[0], length[1], length[2]);
}



//////////////////////////////////////////////////////////////////////////
// Callbacks
//////////////////////////////////////////////////////////////////////////

void Setup(DEM::Domain & dom, void * userData)
{
    UserData & data = *static_cast<UserData *>(userData);

    // ---- 1. virial stress straight from the device contact list ----
    real hsig[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (dom.demaux.nvvint > 0)
    {
        const size_t block = dom.Nthread;
        size_t grid = (dom.demaux.nvvint + block - 1)/block;
        if (grid > 4096) grid = 4096;   // grid-stride loop handles the rest
        cudaMemset(data.devSigma, 0, 6*sizeof(real));
        DEM::CalcStressVV<<<grid, block>>>(dom.pInteractons, dom.pComInteractons,
                                           dom.pDynInteractonsVV, dom.pDynParticlesCU,
                                           dom.pdemaux, data.devSigma);
        cudaMemcpy(hsig, data.devSigma, 6*sizeof(real), cudaMemcpyDeviceToHost);
    }

    const double volume = BoxVolume(dom);
    if (volume <= 0.0) throw new Fatal("TX_PID_periodic: non-positive specimen volume");
    data.measured = Vec3_t(-(double)hsig[0]/volume,
                           -(double)hsig[1]/volume,
                           -(double)hsig[2]/volume);

    // ---- 2. PID control -> strain rates on the periodic cell ----
    const double targetStress  = data.pressure;
    const double pressureScale = std::max(data.pressure, 1.0e-12);
    const size_t controlledAxes = data.shearStage ? 2 : 3;

    Vec3_t strainRate(0.0, 0.0, 0.0);
    for (size_t axis = 0; axis < controlledAxes; ++axis)
        strainRate(axis) = PID(data.controller[axis], targetStress, data.measured(axis),
                               pressureScale, data.dt, data.pidKp, data.pidKi, data.pidKd);
    if (data.shearStage) strainRate(2) = data.axialStrainRate;

    // ---- 3. deform the cell on the host and push it to the device ----
    ScaleCell(dom, strainRate, data.dt);
    DEM::SetCellBox<<<1,1>>>(dom.pdemaux, dom.Xmin, dom.Xmax,
                                          dom.Ymin, dom.Ymax, dom.Zmin, dom.Zmax);

    // ---- 4. keep the Verlet list valid while the cell shrinks ----
    // The device rebuild criterion is grain displacement only and cannot see
    // the cell change, so refresh the list on the device when the cell has
    // changed by more than the Verlet distance (a handful of times per run).
    const Vec3_t length(BoxLength(dom,0), BoxLength(dom,1), BoxLength(dom,2));
    if (!data.verletLengthSet)
    {
        data.lastVerletLength = length;
        data.verletLengthSet  = true;
    }
    else
    {
        const double change = std::max(std::max(fabs(length(0) - data.lastVerletLength(0)),
                                                fabs(length(1) - data.lastVerletLength(1))),
                                       fabs(length(2) - data.lastVerletLength(2)));
        if (change > data.verletDist)
        {
            dom.UpdateContactsDevice();
            data.lastVerletLength = length;
        }
    }
}

void Report(DEM::Domain & dom, void * userData)
{
    UserData & data = *static_cast<UserData *>(userData);

    if (dom.idx_out == 0)
    {
        if (data.stressStrain.is_open()) data.stressStrain.close();
        String filename;
        filename.Printf("%s_walls.res", dom.FileKey.CStr());
        data.stressStrain.open(filename.CStr());
        if (!data.stressStrain)
            throw new Fatal("TX_PID_periodic: cannot open output file <%s>", filename.CStr());

        data.stressStrain << Util::_10_6 << "Time"
                          << Util::_8s << "sx" << Util::_8s << "sy" << Util::_8s << "sz"
                          << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez"
                          << Util::_8s << "ev" << "\n";
    }

    Vec3_t strain;
    for (size_t axis = 0; axis < 3; ++axis)
        strain(axis) = BoxLength(dom, axis) / data.initialLength(axis) - 1.0;

    data.stressStrain << Util::_10_6 << dom.Time
                      << Util::_8s << data.measured(0) << Util::_8s << data.measured(1) << Util::_8s << data.measured(2)
                      << Util::_8s << strain(0)       << Util::_8s << strain(1)       << Util::_8s << strain(2)
                      << Util::_8s << strain(0) + strain(1) + strain(2) << "\n";
}



//////////////////////////////////////////////////////////////////////////
// Input handling
//////////////////////////////////////////////////////////////////////////

template <typename T>
void ReadParameter(std::ifstream & input, T & value, char const * name)
{
    if (!(input >> value)) throw new Fatal("TX_PID_periodic: missing or invalid parameter <%s>", name);
    input.ignore(200, '\n');
}

int main(int argc, char ** argv) try
{
    if (argc < 2)
        throw new Fatal("Usage: %s <filekey> [Nproc] [Nthread]\n", argv[0]);

    const size_t nproc = (argc >= 3) ? atoi(argv[2]) : 1;
    const String fileKey(argv[1]);
    const String inputName(fileKey + ".inp");

    if (!Util::FileExists(inputName))
        throw new Fatal("File <%s> not found", inputName.CStr());

    std::ifstream input(inputName.CStr());
    double verlet, fraction, kn, kt, gn, gt, friction;
    double radius, rRatio, dt, dtOut, lx, ly, lz, density;
    double pressure, isotropicEnd, axialStrainRate, shearEnd, pidKp, pidKi, pidKd;
    size_t renderVideo, seed;
    ReadParameter(input, verlet, "verlet");
    ReadParameter(input, renderVideo, "renderVideo");
    ReadParameter(input, fraction, "fraction");
    ReadParameter(input, kn, "Kn");
    ReadParameter(input, kt, "Kt");
    ReadParameter(input, gn, "Gn");
    ReadParameter(input, gt, "Gt");
    ReadParameter(input, friction, "Mu");
    ReadParameter(input, radius, "R");
    ReadParameter(input, rRatio, "RRatio");
    ReadParameter(input, seed, "seed");
    ReadParameter(input, dt, "dt");
    ReadParameter(input, dtOut, "dtOut");
    ReadParameter(input, lx, "Lx");
    ReadParameter(input, ly, "Ly");
    ReadParameter(input, lz, "Lz");
    ReadParameter(input, density, "rho");
    ReadParameter(input, pressure, "pressure");
    ReadParameter(input, isotropicEnd, "isotropicEnd");
    ReadParameter(input, axialStrainRate, "axialStrainRate");
    ReadParameter(input, shearEnd, "shearEnd");
    ReadParameter(input, pidKp, "pidKp");
    ReadParameter(input, pidKi, "pidKi");
    ReadParameter(input, pidKd, "pidKd");

    UserData data;
    data.renderVideo     = (renderVideo != 0);
    data.dt              = dt;
    data.pressure        = pressure;
    data.axialStrainRate = axialStrainRate;
    data.pidKp           = pidKp;
    data.pidKi           = pidKi;
    data.pidKd           = pidKd;
    data.shearStage      = false;
    data.verletDist      = verlet;
    data.verletLengthSet = false;
    data.measured        = Vec3_t(0.0,0.0,0.0);

    int ndev = 0;
    cudaGetDeviceCount(&ndev);
    if (ndev < 1) throw new Fatal("TX_PID_periodic: no CUDA device found");
    cudaMalloc((void**)&data.devSigma, 6*sizeof(real));

    DEM::Domain dom(&data);
    dom.Alpha         = verlet;
    dom.Beta          = 2.0;
    dom.MostlySpheres = false;      // PairtoCInt must be populated for the CUDA upload
    dom.UseVelocityVerlet = true;   // velocity-Verlet integration
    dom.SphereCoulombMode = 1;      // LAMMPS: Coulomb cap on the total (elastic + dashpot) force
    if (argc >= 4) dom.Nthread = atoi(argv[3]);
    {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        printf("GPU: %s | contact mode %zu | velocity verlet: %s | block = %zu threads\n",
               prop.name, dom.SphereCoulombMode, dom.UseVelocityVerlet ? "on" : "off", dom.Nthread);
    }

    // dompargen generated (random) initial state.
    const Vec3_t minCorner(-0.5*lx, -0.5*ly, -0.5*lz);
    const Vec3_t maxCorner = -minCorner;
    dom.GenSpheresBox(-1, minCorner, maxCorner, radius, density, "HCP",
                      seed, fraction, rRatio);

    // Periodic boundary conditions (FCC style): the periodic cell is the box
    // that contains the generated packing and its size is the control input.
    dom.Xmin = minCorner(0);  dom.Xmax = maxCorner(0);
    dom.Ymin = minCorner(1);  dom.Ymax = maxCorner(1);
    dom.Zmin = minCorner(2);  dom.Zmax = maxCorner(2);
    dom.Per  = Vec3_t(dom.Xmax - dom.Xmin, dom.Ymax - dom.Ymin, dom.Zmax - dom.Zmin);

    Dict properties;
    properties.Set(-1, "Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",
                   kn, kt, gn, gt, friction, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    dom.SetProps(properties);

    // Phase 1: isotropic compression by shrinking the periodic cell.
    cout << "Phase 1: isotropic compression of the periodic cell.\n";
    data.shearStage = false;
    for (size_t axis = 0; axis < 3; ++axis)
    {
        data.controller[axis].integral      = 0.0;
        data.controller[axis].previousError = 0.0;
        data.controller[axis].initialized   = false;
        data.initialLength(axis) = BoxLength(dom, axis);
    }
    if (data.stressStrain.is_open()) data.stressStrain.close();

    const String isotropicKey(fileKey + "_iso");
    dom.Solve(isotropicEnd, dt, dtOut, &Setup, &Report,
              isotropicKey.CStr(), data.renderVideo, nproc);
    dom.Save(isotropicKey.CStr());

    // Phase 2: triaxial shear.  The lateral axes hold the confining pressure
    // through the PID servo; the axial axis is strained at a constant rate.
    cout << "Phase 2: constant-rate triaxial shear of the periodic cell.\n";
    data.shearStage = true;
    data.verletLengthSet = false;
    for (size_t axis = 0; axis < 3; ++axis)
    {
        data.controller[axis].integral      = 0.0;
        data.controller[axis].previousError = 0.0;
        data.controller[axis].initialized   = false;
        data.initialLength(axis) = BoxLength(dom, axis);
    }
    if (data.stressStrain.is_open()) data.stressStrain.close();

    const String shearKey(fileKey + "_shear");
    dom.Solve(shearEnd, dt, dtOut, &Setup, &Report,
              shearKey.CStr(), data.renderVideo, nproc);
    dom.Save(shearKey.CStr());

    cudaFree(data.devSigma);
    cout << "TX_PID_periodic (GPU) finished successfully.\n";
    return 0;
}
MECHSYS_CATCH
