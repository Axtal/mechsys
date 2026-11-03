/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 ************************************************************************/
// 3D periodic Poiseuille duct flow. The duct axis is x; y-z is the square section.

#include <mechsys/mpm/Domain.h>
#include <algorithm>
#include <cmath>
#include <iostream>

struct UserData
{
    double HalfWidth;
    double BodyAccelX;
};

void SetGravity(MPM::Domain & dom, double gx)
{
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        MPM::Particle * p = dom.Particles[ip];
        p->b(0) = p->m*gx;
        p->b(1) = 0.0;
        p->b(2) = 0.0;
    }
}

void Setup(MPM::Domain & dom, void * UD)
{
    UserData * dat = static_cast<UserData *>(UD);
    SetGravity(dom,dat->BodyAccelX);
}

double SquareDuctVmax(double rho, double gx, double half_width, double mu)
{
    const int nterms = 40;
    const double pi = M_PI;
    double sum = 0.0;

    for (int m=0; m<nterms; m++)
    {
        const int n = 2*m + 1;
        const double sign = (m%2 == 0) ? 1.0 : -1.0;
        const double a = 0.5*n*pi;
        sum += sign/(n*n*n)*(1.0 - 1.0/std::cosh(a));
    }

    return 16.0*rho*gx*half_width*half_width*sum/(mu*pi*pi*pi);
}

double MorrisChannelVelocity(double y, double time, double rho, double gx, double height, double mu)
{
    const int nterms = 80;
    const double pi = M_PI;
    const double nu = mu/rho;
    const double steady = 0.5*gx*y*(height - y)/nu;
    double transient = 0.0;

    for (int m=0; m<nterms; m++)
    {
        const int n = 2*m + 1;
        const double nn = (double)n;
        const double arg = nn*pi*y/height;
        const double decay = -time*nn*nn*pi*pi*nu/(height*height);
        transient += 4.0*gx*height*height*std::sin(arg)*std::exp(decay)/(nu*pi*pi*pi*nn*nn*nn);
    }

    return steady - transient;
}

void Report(MPM::Domain & dom, void * UD)
{
    UserData * dat = static_cast<UserData *>(UD);
    if (dom.idx_out % 10 != 0) return;

    double vmax = 0.0;

    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        MPM::Particle * p = dom.Particles[ip];
        vmax = std::max(vmax,p->v(0));
    }

    const double rho0 = 1000.0;
    const double mu = 1.0e-1;
    const double g = dat->BodyAccelX;
    const double height = 2.0*dat->HalfWidth;
    const double vmax_ana = SquareDuctVmax(rho0,g,dat->HalfWidth,mu);
    //const double vmax_ana = MorrisChannelVelocity(0.5*height,dom.Time,rho0,g,height,mu);
    const double rerr = fabs((vmax - vmax_ana)/vmax_ana);

    std::cout << "Time: " << dom.Time
              << ", Particles: " << dom.Particles.Size()
              << ", Vx_max: " << vmax
              << ", Vx_analytical: " << vmax_ana
              << ", Rerror: " << rerr
              << std::endl;
}

void TrimPipeParticles(MPM::Domain & dom, int tag, double radius)
{
    const int delete_tag = -999;
    const double r2 = radius*radius;

    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        MPM::Particle * p = dom.Particles[ip];
        if (p->Tag != tag) continue;

        const double y = p->x(1);
        const double z = p->x(2);
        if (y*y + z*z > r2) p->Tag = delete_tag;
    }

    dom.DelParticles(delete_tag);
}


void FixPipeWallNodes(MPM::Domain & dom, double radius, double dx)
{
    const double wall_radius = radius - 0.75*dx;
    const double wall_r2 = wall_radius*wall_radius;
    size_t fixed_count = 0;

    for (size_t in=0; in<dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        const double r2 = xn(1)*xn(1) + xn(2)*xn(2);
        if (r2 >= wall_r2)
        {
            dom.Nodes[in].FixVeloc();
            fixed_count++;
        }
    }

    std::cout << "No-slip circular-wall nodes: " << fixed_count << std::endl;
}

void FixSquarePipeWallNodes(MPM::Domain & dom, double half_width, double dx)
{
    const double wall_limit = half_width - 0.5*dx;
    size_t fixed_count = 0;

    for (size_t in=0; in<dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        if (fabs(xn(1)) >= wall_limit || fabs(xn(2)) >= wall_limit)
        {
            dom.Nodes[in].FixVeloc();
            fixed_count++;
        }
    }

    std::cout << "No-slip square-wall nodes: " << fixed_count << std::endl;
}


int main(int argc, char **argv) try
{
    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
    dom.ShapeFunction = MPM::SHAPE_BSPLINE_QUADRATIC;
    dom.Periodic[0] = true;

    size_t nproc = std::max<size_t>(1,(size_t)(0.75*omp_get_max_threads()));
    if (argc > 1) nproc = atoi(argv[1]);

    const double rho0 = 1000.0;
    const double mu = 1.0e-1;
    const double gamma = 7.0;
    const int fluid_eos = MPM::EOS_LINEAR; // EOS_TAIT

    const double pipe_diameter = 0.1;
    const double pipe_radius = 0.5*pipe_diameter;
    const double pipe_length = 0.1;
    const int particles_per_diameter = 60;
    const int particles_per_cell = 2;
    const double particle_spacing = pipe_diameter/particles_per_diameter;
    const double dx = particles_per_cell*particle_spacing;

    const double bulk_modulus = 4.9; // 2.0e5;
    const double cs = std::sqrt(bulk_modulus/rho0);
    const double tf = 200.0;
    dat.HalfWidth = pipe_radius;
    dat.BodyAccelX = 1.0e-4; // This is the "gravitational" acceleraction for this problem

    const double vmax_ana = SquareDuctVmax(rho0,dat.BodyAccelX,dat.HalfWidth,mu);
    //const double vmax_ana = MorrisChannelVelocity(dat.HalfWidth,tf,rho0,dat.BodyAccelX,pipe_diameter,mu);
    const double cfl = 0.1;
    const double dt = cfl*dx/(cs + vmax_ana);

    dom.AddRectangularBeam(-1,Vec3_t(0.0,-pipe_radius,-pipe_radius),
                              Vec3_t(pipe_length,pipe_radius,pipe_radius),
                              rho0,particles_per_diameter);
    // TrimPipeParticles(dom,-1,pipe_radius);

    const int wall_cells = 3;
    dom.ResizeDomain(Vec3_t(0.0,-pipe_radius - wall_cells*dx,-pipe_radius - wall_cells*dx),
                     Vec3_t(pipe_length, pipe_radius + wall_cells*dx, pipe_radius + wall_cells*dx),
                     dx);

    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        MPM::Particle * p = dom.Particles[ip];
        p->Constitutive = 4;
        p->Mu = mu;
        p->Cs = cs;
        p->EOSType = fluid_eos;
        p->EOSexp = gamma;
        //p->AVLin = 1.0;
        //p->AVQuad = 2.0;
        //p->FluidDx = dx;
    }

    // FixPipeWallNodes(dom,pipe_radius,dx);
    FixSquarePipeWallNodes(dom,pipe_radius,dx);
    

    std::cout << "********** 3D periodic square-pipe flow **********" << std::endl;
    std::cout << "Square duct side: " << pipe_diameter << " m" << std::endl;
    std::cout << "Pipe length: " << pipe_length << " m" << std::endl;
    std::cout << "Grid spacing: " << dx << " m" << std::endl;
    std::cout << "EOS: " << MPM::EOSName(fluid_eos) << std::endl;
    std::cout << "Sound speed: " << cs << " m/s" << std::endl;
    std::cout << "Body acceleration a: " << dat.BodyAccelX << " m/s2" << std::endl;
    std::cout << "Analytical Vmax at final time: " << vmax_ana << " m/s" << std::endl;
    std::cout << "Reynolds: " << rho0*pipe_diameter*vmax_ana/mu << std::endl;
    std::cout << "Periodic x: " << dom.Periodic[0] << std::endl;
    std::cout << "Time step: " << dt << " s" << std::endl;
    std::cout << "Final time: " << tf << " s" << std::endl;

    dom.Gn = 0.0;
    dom.Solve(tf,dt,tf/100.0,Setup,Report,"tmpm04",true,nproc);

    return 0;
}
MECHSYS_CATCH
