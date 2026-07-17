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
// Cylindrical dry granular collapse with weakly-compressible mu(I) rheology


// MechSys
#include <mechsys/mpm/Domain.h>
#include <algorithm>
#include <math.h>
#include <iostream>

struct UserData
{
    double BodyAccelY;
};


// Set the gravity - body - force
__global__ void SetGravity(MPM::ParticleCU * Par, real gy, MPM::mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic >= mpmaux[0].Npart) return; 
    Par[ic].b.x = 0.0;
    Par[ic].b.y = Par[ic].m * gy;
    Par[ic].b.z = 0.0;
}


void Setup (MPM::Domain & dom, void * UD)
{
    UserData * dat = static_cast<UserData *>(UD);
    SetGravity<<<dom.mpmaux.Npart/dom.Nthread + 1, dom.Nthread>>>
              (dom.pParticlesCU,(real)dat->BodyAccelY,dom.pmpmaux);
}

void Report(MPM::Domain & dom, void * UD)
{
    (void)UD;
    if (dom.idx_out % 10 != 0) return;

    double vmax = 0.0;
    double pmax = 0.0;
    double syy_min = 1.0e100;
    double syy_max = -1.0e100;
    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        MPM::Particle * p = dom.Particles[ip];
        vmax = std::max(vmax,norm(p->v));
        pmax = std::max(pmax,p->Pressure);
        syy_min = std::min(syy_min,p->S(1,1));
        syy_max = std::max(syy_max,p->S(1,1));
    }

    std::cout << "Time: " << dom.Time
              << ", Particles: " << dom.Particles.Size()
              << ", vmax: " << vmax
              << ", pmax: " << pmax
              << ", Syy-range: [" << syy_min << ", " << syy_max << "]" << std::endl;
}


// Create a cylindrical column with its axis along Y.
void MakeCylinderY(MPM::Domain & dom, int Tag, double radius,double center_x, double center_z)
{
    int deleteTag = -999;
    
    for (size_t ip = 0; ip < dom.Particles.Size(); ip++)
    {
        if (dom.Particles[ip]->Tag == Tag)
        {
            double xb = dom.Particles[ip]->x(0) - center_x;
            double zb = dom.Particles[ip]->x(2) - center_z;
            double r = sqrt(xb*xb + zb*zb);
            
            // Delete only based on XZ plane distance (no Y condition)
            if (r > radius)
            {
                dom.Particles[ip]->Tag = deleteTag;
            }
        }
    }
    dom.DelParticles(deleteTag);
}


int main(int argc, char **argv) try
{

    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
    // Coordinate convention for this collapse case:
    // X = main horizontal spreading direction, Y = vertical gravity direction,
    // Z = second horizontal direction. The ground is the plane Y = 0.
    dom.ShapeFunction  = MPM::SHAPE_BSPLINE_QUADRATIC;
    const int constitutive_model = 5; // mu(I) rheology.
    const char * model_name = "mu(I)";
   
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);

    const double mui_K = 26000.0;       // Weakly-compressible granular bulk modulus [Pa].
    const double mui_mu_s = 0.3819;     // Quasi-static friction coefficient.
    const double mui_mu_2 = 0.5718;     // High-inertial-number friction coefficient.
    const double mui_I0 = 0.279;        // Characteristic inertial number.
    const double mui_d = 0.0053;        // Mean grain diameter [m].
    const double mui_rhop = 2000.0;     // Intrinsic grain density used by I = gamma_dot*d*sqrt(rho_p/p) [kg/m3].
    const double solid_volume_fraction = 0.65; // Initial packing fraction; rho_bulk = phi_s*rho_p.
    const double rho_bulk = solid_volume_fraction*mui_rhop;
    const double mui_reg = 1.0e-4;      // Shear-rate regularization [1/s].
    const double mui_cs = std::sqrt(mui_K/rho_bulk);
    double Lx   = 2.0e-1;
    double Ly   = 1.0e-1;
    double Lz   = 2.0e-1;

    int particles_per_cell_per_dim = 2;  // 2×2×2 = 8 particles per cell
    int num_particles_x = 40;
    double particle_spacing = Lx / num_particles_x;
    double Dx = particle_spacing * particles_per_cell_per_dim;
    size_t ndiv = num_particles_x;

    dat.BodyAccelY = -9.81;
    const double velocity_scale = std::sqrt(2.0*fabs(dat.BodyAccelY)*Ly);
    const double cfl = 0.01; // 1.7618550176e-3; // Gives dt = 3.0e-6 s with the current Dx, Cs, and velocity scale.
    const double dt = cfl*Dx/(mui_cs + velocity_scale);

    // Ground is the horizontal X-Z plane at Y=0; gravity acts in -Y.
    double ground_thickness = Dx;
    double ground_ymin = -ground_thickness;
    double ground_ymax = 0.0;

    // Domain bounds with padding (2 ghost cells)
    int ghost_cells = 2;
    double xmin = -Lx - ghost_cells * Dx;
    double xmax = 1.5*Lx + ghost_cells * Dx;
    double ymin = ground_ymin - ghost_cells * Dx;
    double ymax = Ly + ghost_cells * Dx;
    double zmin = -Lz - ghost_cells * Dx;
    double zmax = 1.5*Lz + ghost_cells * Dx;

    dom.AddRectangularBeam(-1, Vec3_t(0.0,0.0,0.0), Vec3_t(Lx,Ly,Lz), rho_bulk, ndiv);
    dom.ResizeDomain(Vec3_t(xmin, ymin, zmin),Vec3_t(xmax, ymax, zmax),Dx);
    double radius = std::min(Lx, Lz) / 2.0;

    //Setting properties for the material points
    for (size_t ip=0; ip < dom.Particles.Size(); ip++)
    {   
        MPM::Particle * p = dom.Particles[ip];
        p->Constitutive = constitutive_model;
        p->K = mui_K;
        p->Cs = mui_cs;
        p->EOSType = MPM::EOS_LINEAR;
        p->EOSexp = 1.0;
        p->MuIStatic = mui_mu_s;
        p->MuIDynamic = mui_mu_2;
        p->MuII0 = mui_I0;
        p->GrainDiameter = mui_d;
        p->GrainRho = mui_rhop;
        p->MuIReg = mui_reg;
    }

    MakeCylinderY(dom, -1, radius, Lx/2.0, Lz/2.0);

    //Setting properties for nodes
    for (size_t in=0; in < dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        if (xn(1) <= ground_ymax + 0.5*Dx)
        {
            dom.Nodes[in].FixVeloc();
        }
    }

    // After ResizeDomain, you can verify:
    size_t Nx = (size_t)((xmax - xmin) / Dx);  // Should be reasonable (~20-30)
    size_t Ny = (size_t)((ymax - ymin) / Dx);
    size_t Nz = (size_t)((zmax - zmin) / Dx);

    std::cout << "Grid: " << Nx << " × " << Ny << " × " << Nz << " = " << Nx*Ny*Nz << " nodes" << std::endl;
    std::cout << "Particles per cell (approx): " << (double)dom.Particles.Size() / (Nx*Ny*Nz) << std::endl;


    dom.Gn = 0.001; //5.0e0;
    double Tf = 1.5; // Final Time
    
    std::cout << "********** Granular dry collapse **********" << std::endl;
    std::cout << "Constitutive model: " << model_name << " (" << constitutive_model << ")" << std::endl;
    std::cout << "Grid spacing: " << Dx << " m" << std::endl;
    std::cout << "Solid volume fraction: " << solid_volume_fraction << std::endl;
    std::cout << "Continuum bulk density: " << rho_bulk << " kg/m3" << std::endl;
    std::cout << "Sound speed: " << mui_cs << " m/s" << std::endl;
    std::cout << "Velocity scale: " << velocity_scale << " m/s" << std::endl;
    std::cout << "CFL: " << cfl << std::endl;
    std::cout << "Parameter values of mu(I) --K/mu_s/mu_2/I0/d/rho_p/reg--: "
              << mui_K << " / " << mui_mu_s << " / " << mui_mu_2
              << " / " << mui_I0 << " / " << mui_d
              << " / " << mui_rhop << " / " << mui_reg << std::endl;
    std::cout << "Time step: " << dt << " s" << std::endl;
    double dtOut = std::min(Tf/100.0,5000.0*dt);
    std::cout << "Output interval: " << dtOut << " s" << std::endl;
    std::cout << "Final time: " << Tf << " s" << std::endl;

    dom.Solve(Tf,dt,dtOut,Setup,Report,"tmpm05",true,Nproc);

    return 0;
}
MECHSYS_CATCH

//
