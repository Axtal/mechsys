/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Cylindrical water-column collapse in a closed 3D box.

#include <mechsys/mpm/Domain.h>
#include <algorithm>
#include <cmath>
#include <iostream>

struct UserData
{
    double Tf;
    double BoxLx;
    double BoxLy;
    double BoxLz;
    double Dx;
};

__global__ void SetGravity(MPM::ParticleCU * Par, real gx, real gy, real gz, MPM::mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= mpmaux[0].Npart) return;
    Par[ic].b.x = Par[ic].m*gx;
    Par[ic].b.y = Par[ic].m*gy;
    Par[ic].b.z = Par[ic].m*gz;
}

void Setup(MPM::Domain & dom, void * UD)
{
    (void)UD;
    const real gx = 0.0;
    const real gy = -9.81;
    const real gz = 0.0;
    SetGravity<<<dom.mpmaux.Npart/dom.Nthread + 1, dom.Nthread>>>(dom.pParticlesCU,gx,gy,gz,dom.pmpmaux);
}

void Report(MPM::Domain & dom, void * UD)
{
    (void)UD;
    if (dom.idx_out % 10 != 0) return;

    double vmax = 0.0;
    double pmax = 0.0;
    double pmin = 1.0e30;
    double jmin = 1.0e30;
    double jmax = -1.0e30;
    double rhomax = 0.0;
    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        vmax = std::max(vmax,norm(dom.Particles[ip]->v));
        pmax = std::max(pmax,dom.Particles[ip]->Pressure);
        pmin = std::min(pmin,dom.Particles[ip]->Pressure);
        if (dom.Particles[ip]->V0 > 0.0 && dom.Particles[ip]->V > 0.0)
        {
            double J = dom.Particles[ip]->V/dom.Particles[ip]->V0;
            jmin = std::min(jmin,J);
            jmax = std::max(jmax,J);
            rhomax = std::max(rhomax,dom.Particles[ip]->m/dom.Particles[ip]->V);
        }
    }

    std::cout << "Time: " << dom.Time
              << ", Particles: " << dom.Particles.Size()
              << ", vmax: " << vmax
              << ", pmax: " << pmax
              << ", pmin: " << pmin
              << ", J-range: [" << jmin << ", " << jmax << "]"
              << ", rhomax: " << rhomax << std::endl;
}

void TrimCylinderY(MPM::Domain & dom, int Tag, double radius, double center_x, double center_z)
{
    const int delete_tag = -999;

    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        if (dom.Particles[ip]->Tag != Tag) continue;

        double dx = dom.Particles[ip]->x(0) - center_x;
        double dz = dom.Particles[ip]->x(2) - center_z;
        if (std::sqrt(dx*dx + dz*dz) > radius)
        {
            dom.Particles[ip]->Tag = delete_tag;
        }
    }

    dom.DelParticles(delete_tag);
}

void InitializeHydrostaticWater(MPM::Domain & dom, int eos_type,
                                double rho0, double c0, double gamma,
                                double gravity, double free_surface_y)
{
    Mat3_t I = OrthoSys::I;
    const double bulk_factor = rho0*c0*c0;

    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        MPM::Particle * par = dom.Particles[ip];
        double depth = std::max(0.0,free_surface_y - par->x(1));
        double p = rho0*gravity*depth;
        set_to_zero(par->E);
        double density_ratio = (eos_type == MPM::EOS_LINEAR) ?
                               1.0 + p/bulk_factor :
                               std::pow(1.0 + gamma*p/bulk_factor,1.0/gamma);
        double J = 1.0/density_ratio;
        double stretch = std::pow(J,1.0/3.0);

        par->Pressure = p;
        par->S = -p*I;
        par->F = stretch*I;
        par->V = par->V0*J;
        par->Dx = Vec3_t(stretch*par->Dx0(0),stretch*par->Dx0(1),stretch*par->Dx0(2));
    }
}

void FixClosedBoxNodes(MPM::Domain & dom, double xmin, double xmax,
                       double ymin, double ymax, double zmin, double zmax,
                       double Dx)
{
    const double tol = 0.5*Dx;
    size_t fixed_count = 0;

    for (size_t in=0; in<dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        bool x_wall = (xn(0) <= xmin + tol) || (xn(0) >= xmax - tol);
        bool bottom_wall = (xn(1) <= ymin + tol);
        bool z_wall = (xn(2) <= zmin + tol) || (xn(2) >= zmax - tol);
        if (x_wall || bottom_wall || z_wall)
        {
            dom.Nodes[in].FixVeloc();
            fixed_count++;
        }
    }

    std::cout << "No-slip bottom/lateral boundary nodes: " << fixed_count << std::endl;
}

int main(int argc, char **argv) try
{
    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;

    // Coordinate convention: X = main horizontal spread direction,
    // Y = vertical gravity direction, Z = second horizontal direction.
    dom.ShapeFunction = MPM::SHAPE_BSPLINE_QUADRATIC;

    size_t Nproc = std::max<size_t>(1,(size_t)(0.75*omp_get_max_threads()));
    if (argc > 1) Nproc = atoi(argv[1]);

    const double rho0 = 1000.0;   // Water density [kg/m3]
    const double mu = 1.0e-3;     // Water dynamic viscosity [Pa s]
    const double gamma = 7.0;     // Tait/Cole EOS exponent for water
    const int fluid_eos = MPM::EOS_TAIT; // Switch to MPM::EOS_LINEAR for the acoustic EOS.
    const double gravity = 9.81;  // Gravity magnitude [m/s2]

    const double water_radius = 8.0e-2;
    const double water_height = 8.0e-2;
    const double box_length = 3.0e-1;
    const double box_lx = box_length;
    const double box_ly = box_length;
    const double box_lz = box_length;
    // Change these two values if you want to relocate the cylinder center.
    const double cylinder_center_x = 0.5*box_lx;
    const double cylinder_center_z = 0.5*box_lz;
    const double center_x = cylinder_center_x;
    const double center_z = cylinder_center_z;

    const int particles_per_diameter = 70;
    const int particles_per_cell_per_dim = 2;
    const double particle_spacing = 2.0*water_radius/particles_per_diameter;
    const double Dx = particles_per_cell_per_dim*particle_spacing;
    const size_t ndiv = particles_per_diameter;

    const double target_bulk_modulus = 2.0e5; // Telikicherla and Moutsanidis (2025) dam-break scale.
    const double sound_speed = std::sqrt(target_bulk_modulus/rho0);
    const double velocity_scale = std::sqrt(2.0*gravity*box_ly);
    const double cfl = 0.07;
    // Set both artificial-viscosity coefficients to zero to disable the filter.
    const double av_linear = 1.0;
    const double av_quadratic = 2.0;
    const double dt = cfl*Dx/(sound_speed + velocity_scale);
    const double Tf = 5.0;

    dat.Tf = Tf;
    dat.BoxLx = box_lx;
    dat.BoxLy = box_ly;
    dat.BoxLz = box_lz;
    dat.Dx = Dx;

    const double water_xmin = center_x - water_radius;
    const double water_xmax = center_x + water_radius;
    const double water_ymin = 0.0;
    const double water_ymax = water_height;
    const double water_zmin = center_z - water_radius;
    const double water_zmax = center_z + water_radius;

    dom.AddRectangularBeam(-1,Vec3_t(water_xmin,water_ymin,water_zmin),
                              Vec3_t(water_xmax,water_ymax,water_zmax),rho0,ndiv);
    TrimCylinderY(dom,-1,water_radius,center_x,center_z);

    const int ghost_cells = 3;
    const double xmin = 0.0;
    const double xmax = box_lx;
    const double ymin = 0.0;
    const double ymax = box_ly;
    const double zmin = 0.0;
    const double zmax = box_lz;

    dom.ResizeDomain(Vec3_t(xmin - ghost_cells*Dx,ymin - ghost_cells*Dx,zmin - ghost_cells*Dx),
                     Vec3_t(xmax + ghost_cells*Dx,ymax + ghost_cells*Dx,zmax + ghost_cells*Dx),Dx);

    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        dom.Particles[ip]->Constitutive = 4;
        dom.Particles[ip]->Mu = mu;
        dom.Particles[ip]->Cs = sound_speed;
        dom.Particles[ip]->EOSType = fluid_eos;
        dom.Particles[ip]->EOSexp = gamma;
        dom.Particles[ip]->AVLin = av_linear;
        dom.Particles[ip]->AVQuad = av_quadratic;
        dom.Particles[ip]->FluidDx = Dx;
    }

    InitializeHydrostaticWater(dom,fluid_eos,rho0,sound_speed,gamma,gravity,water_ymax);
    FixClosedBoxNodes(dom,xmin,xmax,ymin,ymax,zmin,zmax,Dx);

    std::cout << "3D closed-box cylindrical water collapse" << std::endl;
    std::cout << "Grid spacing: " << Dx << " m" << std::endl;
    std::cout << "EOS: " << MPM::EOSName(fluid_eos) << std::endl;
    std::cout << "Bulk modulus: " << target_bulk_modulus << " Pa" << std::endl;
    std::cout << "Sound speed: " << sound_speed << " m/s" << std::endl;
    std::cout << "Velocity scale for dt: " << velocity_scale << " m/s" << std::endl;
    std::cout << "CFL factor: " << cfl << std::endl;
    std::cout << "Pressure filter cL/cQ: " << av_linear << " / " << av_quadratic << std::endl;
    std::cout << "Time step: " << dt << " s" << std::endl;
    std::cout << "Final time: " << Tf << " s" << std::endl;
    std::cout << "Water particles: " << dom.Particles.Size() << std::endl;

    dom.Gn = 0.0;
    dom.Solve(Tf,dt,Tf/100.0,Setup,Report,"tmpm03",true,Nproc);

    return 0;
}
MECHSYS_CATCH
