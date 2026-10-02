// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

#include <algorithm>
#include <cmath>
#include <fstream>

using std::cout;


struct PIDController
{
    double integral;
    double previousError;
    bool   initialized;
};

double PID(PIDController & controller, double targetStress, double measuredStress,
           double pressureScale, double dt, double kp, double ki, double kd)
{
    // The normalized error makes these gains independent of the stress units.
    const double maxSpeed = 0.10;
    const double error = (targetStress - measuredStress) / pressureScale;
    if (!controller.initialized)
    {
        controller.previousError = error;
        controller.initialized   = true;
    }

    controller.integral += 0.5 * (error + controller.previousError) * dt;
    if (ki > 0.0)
    {
        const double integralLimit = maxSpeed / ki;
        controller.integral = std::max(
            -integralLimit, std::min(controller.integral, integralLimit));
    }

    const double derivative = (error - controller.previousError) / dt;
    controller.previousError = error;

    const double velocity = kp * error + ki * controller.integral + kd * derivative;
    return std::max(-maxSpeed, std::min(velocity, maxSpeed));
}

struct UserData
{
    bool          renderVideo;
    size_t        wallIndex;
    double        dt;
    double        pressure;
    double        axialVelocity;
    double        pidKp;
    double        pidKi;
    double        pidKd;
    bool          shearStage;
    Vec3_t        initialLength;
    PIDController controller[3];
    std::ofstream stressStrain;
};

double WallLength(DEM::Domain const & dom, UserData const & data, size_t axis)
{
    return dom.Particles[data.wallIndex + 2*axis]->x(axis)
         - dom.Particles[data.wallIndex + 2*axis + 1]->x(axis);
}

double WallArea(DEM::Domain const & dom, UserData const & data, size_t axis)
{
    return WallLength(dom, data, (axis + 1) % 3)
         * WallLength(dom, data, (axis + 2) % 3);
}

double MeasuredStress(DEM::Domain const & dom, UserData const & data, size_t axis)
{
    const double area = WallArea(dom, data, axis);
    if (area <= 0.0) throw new Fatal("PID_v2: non-positive specimen area");

    DEM::Particle const * positiveWall = dom.Particles[data.wallIndex + 2*axis];
    DEM::Particle const * negativeWall = dom.Particles[data.wallIndex + 2*axis + 1];
    return -0.5 * (positiveWall->F(axis) - negativeWall->F(axis)) / area;
}

void SetWallVelocity(DEM::Domain & dom, UserData const & data, size_t axis, double velocity)
{
    DEM::Particle * positiveWall = dom.Particles[data.wallIndex + 2*axis];
    DEM::Particle * negativeWall = dom.Particles[data.wallIndex + 2*axis + 1];
    positiveWall->v = Vec3_t(0.0, 0.0, 0.0);
    negativeWall->v = Vec3_t(0.0, 0.0, 0.0);
    positiveWall->v(axis) =  velocity;
    negativeWall->v(axis) = -velocity;
}

void Setup(DEM::Domain & dom, void * userData)
{
    UserData & data = *static_cast<UserData *>(userData);
    const double targetStress = -data.pressure;
    const size_t controlledAxes = data.shearStage ? 2 : 3;
    const double pressureScale = std::max(data.pressure, 1.0e-12);
    for (size_t axis = 0; axis < controlledAxes; ++axis)
    {
        const double velocity = PID(
            data.controller[axis],
            targetStress, MeasuredStress(dom, data, axis), pressureScale, data.dt,
            data.pidKp, data.pidKi, data.pidKd);
        SetWallVelocity(dom, data, axis, velocity);
    }

    if (data.shearStage)
    {
        SetWallVelocity(dom, data, 2, data.axialVelocity);
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
            throw new Fatal("PID_v2: cannot open output file <%s>", filename.CStr());

        data.stressStrain << Util::_10_6 << "Time"
                          << Util::_8s << "sx" << Util::_8s << "sy" << Util::_8s << "sz"
                          << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez"
                          << Util::_8s << "ev" << "\n";
    }

    Vec3_t stress;
    Vec3_t strain;
    for (size_t axis = 0; axis < 3; ++axis)
    {
        stress(axis) = MeasuredStress(dom, data, axis);
        strain(axis) = WallLength(dom, data, axis) / data.initialLength(axis) - 1.0;
    }

    data.stressStrain << Util::_10_6 << dom.Time
                      << Util::_8s << stress(0) << Util::_8s << stress(1) << Util::_8s << stress(2)
                      << Util::_8s << strain(0) << Util::_8s << strain(1) << Util::_8s << strain(2)
                      << Util::_8s << strain(0) + strain(1) + strain(2) << "\n";

    if (data.renderVideo) dom.WriteXDMF("pid");
}


template <typename T>
void ReadParameter(std::ifstream & input, T & value, char const * name)
{
    if (!(input >> value)) throw new Fatal("PID_v2: missing or invalid parameter <%s>", name);
    input.ignore(200, '\n');
}

int main(int argc, char ** argv) try
{
    if (argc < 2)
        throw new Fatal("Usage: %s <filekey> [Nproc]\n", argv[0]);

    const size_t nproc = (argc >= 3) ? atoi(argv[2]) : 1;
    const String fileKey(argv[1]);
    const String inputName(fileKey + ".inp");

    if (!Util::FileExists(inputName))
        throw new Fatal("File <%s> not found", inputName.CStr());

    std::ifstream input(inputName.CStr());
    double verlet, fraction, kn, kt, gn, gt, friction;
    double radius, rRatio, dt, dtOut, lx, ly, lz, density;
    double pressure, isotropicEnd, axialVelocity, shearEnd, pidKp, pidKi, pidKd;
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
    ReadParameter(input, axialVelocity, "axialVelocity");
    ReadParameter(input, shearEnd, "shearEnd");
    ReadParameter(input, pidKp, "pidKp");
    ReadParameter(input, pidKi, "pidKi");
    ReadParameter(input, pidKd, "pidKd");

    UserData data;
    data.renderVideo   = (renderVideo != 0);
    data.dt            = dt;
    data.pressure      = pressure;
    data.axialVelocity = axialVelocity;
    data.pidKp         = pidKp;
    data.pidKi         = pidKi;
    data.pidKd         = pidKd;

    DEM::Domain dom(&data);
    dom.Alpha         = verlet;
    dom.Beta          = 2.0;
    dom.MostlySpheres = false;

    const Vec3_t minCorner(-0.5*lx, -0.5*ly, -0.5*lz);
    const Vec3_t maxCorner = -minCorner;
    dom.GenSpheresBox(-1, minCorner, maxCorner, radius, density, "HCP",
                      seed, fraction, rRatio);

    data.wallIndex = dom.Particles.Size();

    dom.GenBoundingBox(-2, 0.02*radius, 1.3, false);

    Dict properties;
    properties.Set(-1, "Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",
                   kn, kt, gn, gt, friction, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for (int tag = -2; tag >= -7; --tag)
    {
        properties.Set(tag, "Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",
                       kn, kt, gn, gt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    dom.SetProps(properties);

    for (size_t i = 0; i < dom.Particles.Size(); ++i) dom.Particles[i]->Initialize(i);

    cout << "Phase 1: isotropic compression.\n";
    data.shearStage = false;
    for (size_t axis = 0; axis < 3; ++axis)
    {
        data.controller[axis].integral      = 0.0;
        data.controller[axis].previousError = 0.0;
        data.controller[axis].initialized   = false;
        data.initialLength(axis) = WallLength(dom, data, axis);

        DEM::Particle * positiveWall = dom.Particles[data.wallIndex + 2*axis];
        DEM::Particle * negativeWall = dom.Particles[data.wallIndex + 2*axis + 1];
        positiveWall->Ff = Vec3_t(0.0, 0.0, 0.0);
        negativeWall->Ff = Vec3_t(0.0, 0.0, 0.0);
        positiveWall->FixVeloc();
        negativeWall->FixVeloc();
    }
    if (data.stressStrain.is_open()) data.stressStrain.close();

    const String isotropicKey(fileKey + "_iso");
    dom.Solve(isotropicEnd, dt, dtOut, &Setup, &Report,
              isotropicKey.CStr(), data.renderVideo, nproc);
    dom.Save(isotropicKey.CStr());

    cout << "Phase 2: constant-velocity triaxial shear.\n";
    data.shearStage = true;
    for (size_t axis = 0; axis < 3; ++axis)
    {
        data.controller[axis].integral      = 0.0;
        data.controller[axis].previousError = 0.0;
        data.controller[axis].initialized   = false;
        data.initialLength(axis) = WallLength(dom, data, axis);

        DEM::Particle * positiveWall = dom.Particles[data.wallIndex + 2*axis];
        DEM::Particle * negativeWall = dom.Particles[data.wallIndex + 2*axis + 1];
        positiveWall->Ff = Vec3_t(0.0, 0.0, 0.0);
        negativeWall->Ff = Vec3_t(0.0, 0.0, 0.0);
        positiveWall->FixVeloc();
        negativeWall->FixVeloc();
    }
    if (data.stressStrain.is_open()) data.stressStrain.close();

    const String shearKey(fileKey + "_shear");
    dom.Solve(shearEnd, dt, dtOut, &Setup, &Report,
              shearKey.CStr(), data.renderVideo, nproc);
    dom.Save(shearKey.CStr());

    cout << "PID_v2 finished successfully.\n";
    return 0;
}
MECHSYS_CATCH
