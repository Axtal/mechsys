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
// Cylindrical dry granular collapse


// MechSys
#include <mechsys/mpm/Domain.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    thrust::device_vector<size_t> ForcePar;
    Array<MPM::Particle *> PlaneBeamPar0; //Particles at each plane of the beam
    Array<double >               x0;
    Array<double >               x2; 
    double                       Tf; //Final simulation time
    size_t                     ndiv;
    size_t *                 pnFpar;
    double                       Lx;
    double                       Ly;
    double                       Lz;
    double                       E;
    double                       BodyAccelY;
    real                         Fb; //Maximun applied body force
    real *                       pFb;
    std::ofstream                oss_ss; // The output file
};

enum GranularModel
{
    MODEL_DRUCKER_PRAGER = 1,
    MODEL_BUI_DRUCKER_PRAGER = 2,
    MODEL_HYPOPLASTIC  = 3
};

bool IsSoilConstitutive(int constitutive_model)
{
    return constitutive_model == MODEL_DRUCKER_PRAGER ||
           constitutive_model == MODEL_BUI_DRUCKER_PRAGER ||
           constitutive_model == MODEL_HYPOPLASTIC;
}

const char * SoilModelName(int constitutive_model)
{
    if (constitutive_model == MODEL_DRUCKER_PRAGER) return "Drucker-Prager";
    if (constitutive_model == MODEL_BUI_DRUCKER_PRAGER) return "Bui Drucker-Prager";
    if (constitutive_model == MODEL_HYPOPLASTIC) return "Hypoplastic";
    return "Unsupported soil model";
}


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
    double jmin = 1.0e100;
    double jmax = -1.0e100;
    double ymin = 1.0e100;
    double ymax = -1.0e100;
    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        MPM::Particle * p = dom.Particles[ip];
        vmax = std::max(vmax,norm(p->v));
        pmax = std::max(pmax,p->Pressure);
        syy_min = std::min(syy_min,p->S(1,1));
        syy_max = std::max(syy_max,p->S(1,1));
        if (p->V0 > 0.0)
        {
            double J = p->V/p->V0;
            jmin = std::min(jmin,J);
            jmax = std::max(jmax,J);
        }
        ymin = std::min(ymin,p->x(1));
        ymax = std::max(ymax,p->x(1));
    }

    std::cout << "Time: " << dom.Time
              << ", Particles: " << dom.Particles.Size()
              << ", vmax: " << vmax
              << ", pmax: " << pmax
              << ", Syy-range: [" << syy_min << ", " << syy_max << "]" << std::endl;
              //<< ", Jmin: " << jmin
              //<< ", Jmax: " << jmax
              //<< ", y-range: [" << ymin << ", " << ymax << "]" << std::endl;
}


// Función para crear cilindro con eje en Y
void MakeCylinderY(MPM::Domain & dom, int Tag, double radius,double center_x, double center_z)
{
    int deleteTag = -999;
    size_t deleted_count = 0;
    
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
                deleted_count++;
            }
        }
    }
    dom.DelParticles(deleteTag);
}

void InitializeGeostaticStress(MPM::Domain & dom, double rho, double gravity,
                               double free_surface_y, double k0)
{
    for (size_t ip=0; ip<dom.Particles.Size(); ip++)
    {
        MPM::Particle * p = dom.Particles[ip];
        if (!IsSoilConstitutive(p->Constitutive)) continue;

        double depth = std::max(0.0,free_surface_y - p->x(1));
        double sigma_v = -rho*gravity*depth;
        double sigma_h = k0*sigma_v;

        set_to_zero(p->S);
        p->S(0,0) = sigma_h;
        p->S(1,1) = sigma_v;
        p->S(2,2) = sigma_h;
        p->Pressure = -Trace(p->S)/3.0;
    }
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
    const int constitutive_model = MODEL_BUI_DRUCKER_PRAGER; // Soil options: MODEL_DRUCKER_PRAGER, MODEL_BUI_DRUCKER_PRAGER, MODEL_HYPOPLASTIC.
    const bool initialize_geostatic_stress = true;    // Writes initial vertical self-weight stress to the first output.
   
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);
    // size_t ndiv = 10; //number of divisions per x lenght - keep it one tenth of the lenght Lx
    double Nu	= 0.3;    // Poisson ratio
    double E	= 5.0e7; // Elastic modulus (Pa)
    double K    = E/(3.0*(1.0-2.0*Nu)); // Bulk modulus (Pa)
    double G	= E/(2.0*(1.0+Nu));     // Shear modulus (Pa)
    double c    = 0.0;      // Cohesion (Pa); dry cohesionless sand uses c = 0.
    double phi  = 28.0;    // Friction angle (°)
    double psi  = 0.0;     // Dilation angle (°)
    double rho  = 1800.0;  // Density of the material (kg/m3)
    // double Cs = sqrt(E/rho); // Legacy elastic sound speed; kept only for comparison.
    double Lx   = 2.0e-1;
    double Ly   = 1.0e-1;
    double Lz   = 2.0e-1;

    // === NEW CODE ===
    int particles_per_cell_per_dim = 2;  // 2×2×2 = 8 particles per cell
    // Or use 3×3×3 = 27 particles per cell for better quality

    // Desired number of particles across the main horizontal X direction.
    int num_particles_x = 40;

    // Calculate particle spacing
    double particle_spacing = Lx / num_particles_x;

    // Cell size = particle_spacing * particles_per_cell_per_dim
    double Dx = particle_spacing * particles_per_cell_per_dim;

    // AddRectangularBeam divides the X length by ndiv.
    size_t ndiv = num_particles_x;


    // The following 3 lines just help to keep the ratio of a mount of total material points
    // double Lx_ref = 1.0e-1;  // Lx de referencia = 0.1
    // size_t ndiv_ref = 50;     // ndiv para Lx_ref
    // size_t ndiv = (size_t)(ndiv_ref * (Lx / Lx_ref));

    double h    = Lx/ndiv;
    // double dt = 0.1*h/Cs; // Elastic CFL estimate; kept only for comparison.

    const double gravity = 9.81;
    const double tCFL = 0.1;
    const double vmax_est = std::max(std::sqrt(2.0*gravity*Ly),1.0e-12);
    const double elastic_wave_speed = std::sqrt((K + 4.0*G/3.0)/rho);
    double signal_speed = elastic_wave_speed;
    double dt = 0.1*h/signal_speed;
    const char * model_name = SoilModelName(constitutive_model);

    if (!IsSoilConstitutive(constitutive_model))
    {
        std::cout << "Unsupported constitutive_model in tmpm_cu_02.cu" << std::endl;
        std::cout << "1 => Drucker-Prager, 2 => Bui Drucker-Prager, 3 => Hypoplastic" << std::endl;
        return 1;
    }
    else if (constitutive_model == MODEL_HYPOPLASTIC)
    {
        // Wang et al. (2024) granular SPH CFL form:
        // dt = tCFL*min(h/(10*|v|max), h/cs, rho*h^2/eta_max). The current
        // MPM hypoplastic law has no eta_max state variable, so only the
        // convective, signal-speed, and initial gravity-acceleration limits are used.
        signal_speed = 100.0*vmax_est;
        const double dt_velocity = h/(10.0*vmax_est);
        const double dt_signal = h/signal_speed;
        const double dt_acceleration = std::sqrt(h/gravity);
        dt = tCFL*std::min(dt_velocity,std::min(dt_signal,dt_acceleration));
        model_name = "Hypoplastic";
    }

    // double dt = 0.1*h/Cs; // Legacy elastic CFL estimate kept only for comparison.
    // double Dx   = 2.0*Lx/ndiv; // cell grid size dx
    dat.ndiv    = ndiv;
    dat.Lx      = Lx;
    dat.Ly      = Ly;
    dat.Lz      = Lz;
    dat.E       = E;
    dat.BodyAccelY = -gravity;

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

    dom.AddRectangularBeam(-1, Vec3_t(0.0,0.0,0.0), Vec3_t(Lx,Ly,Lz), rho, ndiv);
    dom.ResizeDomain(Vec3_t(xmin, ymin, zmin),Vec3_t(xmax, ymax, zmax),Dx);
    double radius = std::min(Lx, Lz) / 2.0;

    const double geostatic_k0 = std::max(0.0,1.0 - std::sin(phi*M_PI/180.0)); // Jaky at-rest estimate for soil.
    // K0 represents the initial confined soil state before release.
    // It is not an active side wall; after t=0 those lateral stresses are free to relax.

    //Setting properties for the material points
    for (size_t ip=0; ip < dom.Particles.Size(); ip++)
    {   
        MPM::Particle * p = dom.Particles[ip];
        p->Constitutive = constitutive_model;
        p->Cohesion = c;

        if (constitutive_model == MODEL_DRUCKER_PRAGER || constitutive_model == MODEL_BUI_DRUCKER_PRAGER)
        {
            p->K = K;
            p->G = G;
            p->Friction = phi*M_PI/180.0;
            p->Dilation = psi*M_PI/180.0;
            p->Fail = 3;  // 3 = non-associated flow.
        }
    }

    MakeCylinderY(dom, -1, radius, Lx/2.0, Lz/2.0);
    if (initialize_geostatic_stress)
    {
        InitializeGeostaticStress(dom,rho,gravity,Ly,geostatic_k0);
    }

    dat.Fb = 0.0;
    dat.pFb = NULL;
    dat.pnFpar = NULL;

    //Setting properties for nodes
    for (size_t in=0; in < dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        //if (xn(1)=ground_ymax)
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
    dat.Tf = Tf;
    
    std::cout << "********** Granular dry collapse **********" << std::endl;
    std::cout << "Constitutive model: " << model_name << " (" << constitutive_model << ")" << std::endl;
    std::cout << "Grid spacing: " << Dx << " m" << std::endl;
    std::cout << "Signal speed for dt: " << signal_speed << " m/s" << std::endl;
    std::cout << "Time step: " << dt << " s" << std::endl;
    double dtOut = std::min(Tf/100.0,5000.0*dt);
    std::cout << "Output interval: " << dtOut << " s" << std::endl;
    std::cout << "Initial stress enabled: " << initialize_geostatic_stress
              << ", K0 used: " << geostatic_k0
              << " (Jaky)" << std::endl;
    std::cout << "Final time: " << Tf << " s" << std::endl;

    dom.Solve(Tf,dt,dtOut,Setup,Report,"tmpm02",true,Nproc);

    return 0;
    //Vec3_t a(1.0,2.0,3.0);
    //Vec3_t b(4.0,5.0,6.0);
    //real3    c = make_real3(1.0,2.0,3.0);
    //real3    d = make_real3(4.0,5.0,6.0);
    //Mat3    m3;
    //m3.print();
    //m3 = Dyad(c,d);
    //Mat3_t m3t;
    //Dyad(a,b,m3t);
    //m3.print();
    //std::cout << m3t << std::endl;
}
MECHSYS_CATCH

//
