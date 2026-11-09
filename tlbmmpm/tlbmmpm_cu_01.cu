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
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Elastic beam.


// CUDA LBMMPM.  The CMake target adds --fmad=false so the facet-edge
// classification matches the CPU SolAlt path exactly.

// MechSys
#include <mechsys/lbmmpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    double                         * Vel;
#ifdef USE_CUDA
    double                         * dVel;
#endif
    double                           rho;
    double                            Tf;
    double                            Ly;
    double                            Lz;
    double                            xc;
    double                            nu;
    double                            CD;
    double                            Re;
    double                          umax;
    Vec3_t                    fluidforce;
    Array<MPM::Particle *>    EndBeamPar;
    Array<MPM::Particle *>     middlePar;
    Array<double >                    x0;
    Array<double >                    x2; 
    std::ofstream                oss_ss1; // The output of Cd & Re & Dx & Dy changed with Time
    std::ofstream                oss_ss2; // The output of x&y position 
};

#ifdef USE_CUDA
__global__ void t01_inlet(real * F, real * Rho, real3 * Vel, real const * dVel,
                          real scale, FLBM::lbm_aux const * a)
{
    int t = threadIdx.x + blockIdx.x*blockDim.x;
    int Ny = (int)a[0].Ny, Nz = (int)a[0].Nz;
    if (t >= Ny*Nz) return;
    int i = t % Ny, j = t / Ny;
    size_t n = (size_t)i*a[0].Nx + (size_t)j*a[0].Nx*a[0].Ny;
    real f[15];
    for (int k=0;k<15;k++) f[k] = F[FLBM::FIDX(a[0].Ncells,a[0].Nneigh,0,k,n)];
    real rho = (f[0]+f[3]+f[4]+f[5]+f[6]+2.0*(f[2]+f[8]+f[10]+f[12]+f[14]))
             /(1.0-scale*dVel[j]/a[0].Cs);
    real u = scale*dVel[j]/a[0].Cs;
    f[1]  = f[2]  + (2.0/3.0 )*rho*u;
    f[7]  = f[8]  + (1.0/12.0)*rho*u;
    f[9]  = f[10] + (1.0/12.0)*rho*u;
    f[11] = f[12] + (1.0/12.0)*rho*u;
    f[13] = f[14] + (1.0/12.0)*rho*u;
    for (int k=0;k<15;k++) F[FLBM::FIDX(a[0].Ncells,a[0].Nneigh,0,k,n)] = f[k];
    real r=0.0; real3 v=make_real3(0.0,0.0,0.0);
    for (int k=0;k<15;k++) { r += f[k]; v = v + f[k]*a[0].C[k]; }
    Rho[n] = r; Vel[n] = (a[0].Cs/r)*v;
}

__global__ void t01_outlet(real * F, real * Rho, real3 * Vel, real rhof,
                           FLBM::lbm_aux const * a)
{
    int t = threadIdx.x + blockIdx.x*blockDim.x;
    int Nx=(int)a[0].Nx, Ny=(int)a[0].Ny, Nz=(int)a[0].Nz;
    if (t >= Ny*Nz) return;
    int i=t%Ny, j=t/Ny;
    size_t n = (size_t)(Nx-1) + (size_t)i*Nx + (size_t)j*Nx*Ny;
    real f[15];
    for (int k=0;k<15;k++) f[k] = F[FLBM::FIDX(a[0].Ncells,a[0].Nneigh,0,k,n)];
    real vx = (f[0]+f[3]+f[4]+f[5]+f[6]+2.0*(f[1]+f[7]+f[9]+f[11]+f[13]))/rhof - 1.0;
    f[2]  = f[1]  - (2.0/3.0 )*rhof*vx;
    f[8]  = f[7]  - (1.0/12.0)*rhof*vx;
    f[10] = f[9]  - (1.0/12.0)*rhof*vx;
    f[12] = f[11] - (1.0/12.0)*rhof*vx;
    f[14] = f[13] - (1.0/12.0)*rhof*vx;
    for (int k=0;k<15;k++) F[FLBM::FIDX(a[0].Ncells,a[0].Nneigh,0,k,n)] = f[k];
    real r=0.0; real3 v=make_real3(0.0,0.0,0.0);
    for (int k=0;k<15;k++) { r += f[k]; v = v + f[k]*a[0].C[k]; }
    Rho[n] = r; Vel[n] = (a[0].Cs/r)*v;
}

__global__ void t01_top(real * F, real * Rho, real3 * Vel, FLBM::lbm_aux const * a)
{
    int t = threadIdx.x + blockIdx.x*blockDim.x;
    int Nx=(int)a[0].Nx, Ny=(int)a[0].Ny, Nz=(int)a[0].Nz;
    if (t >= Nx*Ny) return;
    int i=t%Nx, j=t/Nx;
    size_t n = (size_t)i + (size_t)j*Nx + (size_t)(Nz-1)*Nx*Ny;
    real f[15];
    for (int k=0;k<15;k++) f[k] = F[FLBM::FIDX(a[0].Ncells,a[0].Nneigh,0,k,n)];
    f[6]=f[5]; f[8]=f[10]; f[9]=f[7]; f[12]=f[14]; f[13]=f[11];
    for (int k=0;k<15;k++) F[FLBM::FIDX(a[0].Ncells,a[0].Nneigh,0,k,n)] = f[k];
    real r=0.0; real3 v=make_real3(0.0,0.0,0.0);
    for (int k=0;k<15;k++) { r += f[k]; v = v + f[k]*a[0].C[k]; }
    Rho[n] = r; Vel[n] = (a[0].Cs/r)*v;
}
#endif

void Setup (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double scale = std::min(10.0*dom.Time/dat.Tf,1.0);
#ifdef USE_CUDA
    size_t B = dom.Nthread;
    int Ny=(int)dom.LBMDOM.Ndim(1), Nz=(int)dom.LBMDOM.Ndim(2);
    int Nx=(int)dom.LBMDOM.Ndim(0);
    t01_inlet <<<(Ny*Nz+B-1)/B,B>>>(dom.LBMDOM.pF,dom.LBMDOM.pRho,dom.LBMDOM.pVel,
                                    dat.dVel,(real)scale,dom.LBMDOM.plbmaux);
    t01_outlet<<<(Ny*Nz+B-1)/B,B>>>(dom.LBMDOM.pF,dom.LBMDOM.pRho,dom.LBMDOM.pVel,
                                    dat.rho,dom.LBMDOM.plbmaux);
    t01_top   <<<(Nx*Ny+B-1)/B,B>>>(dom.LBMDOM.pF,dom.LBMDOM.pRho,dom.LBMDOM.pVel,
                                    dom.LBMDOM.plbmaux);
#else
    // Cells with prescribed velocity (left boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<dom.LBMDOM.Ndim(1); ++i) //y-axis
    for (size_t j=0; j<dom.LBMDOM.Ndim(2); ++j) //z-axis
    {
        double * f = dom.LBMDOM.F[0][0][i][j];
        double rho = (f[0]+f[3]+f[4]+f[5]+f[6]+ 2.0*(f[2]+f[8]+f[10]+f[12]+f[14]))/(1.0-scale*dat.Vel[j]/dom.LBMDOM.Cs);
        f[1] = f[2] + (2.0/3.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        f[7] = f[8] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        f[9] = f[10] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        f[11]= f[12] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        f[13]= f[14] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        dom.LBMDOM.Vel[0][0][i][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][0][i][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][0][i][j] +=  dom.LBMDOM.F[0][0][i][j][k];
            dom.LBMDOM.Vel[0][0][i][j] +=  dom.LBMDOM.F[0][0][i][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][0][i][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][0][i][j];
    }

    // Cells with prescribed density (right boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<dom.LBMDOM.Ndim(1); ++i) //y-axis
    for (size_t j=0; j<dom.LBMDOM.Ndim(2); ++j) //z-axis
    {
        double * f = dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j];
        double vx = (f[0]+f[3]+f[4]+f[5]+f[6] + 2.0*(f[1]+f[7]+f[9]+f[11]+f[13]))/dat.rho - 1.0;
        f[2] = f[1] - (2.0/3.0)*dat.rho*vx;
        f[8] = f[7] - (1.0/12.0)*dat.rho*vx;
        f[10]= f[9] - (1.0/12.0)*dat.rho*vx;
        f[12]= f[11] - (1.0/12.0)*dat.rho*vx;
        f[14]= f[13] - (1.0/12.0)*dat.rho*vx;
        dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j] +=  dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j][k];
            dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] +=  dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j];
    }

    // free slip (top boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<dom.LBMDOM.Ndim(0); ++i) //x-axis
    for (size_t j=0; j<dom.LBMDOM.Ndim(1); ++j) //y-axis
    {
        double * f = dom.LBMDOM.F[0][i][j][dom.LBMDOM.Ndim(2)-1];
        f[6]  = f[5];
        f[8]  = f[10];
        f[9]  = f[7];
        f[12] = f[14];
        f[13] = f[11];
        dom.LBMDOM.Vel[0][i][j][dom.LBMDOM.Ndim(2)-1] = OrthoSys::O;
        dom.LBMDOM.Rho[0][i][j][dom.LBMDOM.Ndim(2)-1] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][i][j][dom.LBMDOM.Ndim(2)-1] +=  dom.LBMDOM.F[0][i][j][dom.LBMDOM.Ndim(2)-1][k];
            dom.LBMDOM.Vel[0][i][j][dom.LBMDOM.Ndim(2)-1] +=  dom.LBMDOM.F[0][i][j][dom.LBMDOM.Ndim(2)-1][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][i][j][dom.LBMDOM.Ndim(2)-1] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][i][j][dom.LBMDOM.Ndim(2)-1];
    }
#endif
}

void Report (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double x     = 0.0;
    double z     = 0.0;
    double xf    = 0.0; // Dx/b
    double zf    = 0.0; // Dz/b
    double x1    = 0.0; // position in x-axis _simulation
    double z1    = 0.0; // position in z-axis _simulation
    double vel2  = 0.0; // velocity of fluid
    double vel   = 0.0; 
    double force = 0.0;
    int    i     = 0;
    int    j     = 0;
    if (dom.idx_out==0)
    {
        String fs1;
        fs1.Printf("tlbmmpm_cu_01position.res");
        dat.oss_ss1.open(fs1.CStr());
        dat.oss_ss1 << Util::_10_6 << "X-position" << Util::_8s << "Z-position" << std::endl;

        String fs2;
        fs2.Printf("tlbmmpm_cu_01time.res");
        dat.oss_ss2.open(fs2.CStr());
        dat.oss_ss2 << Util::_10_6 << "Time" << Util::_8s << "Re" << Util::_8s << "CD" << Util::_8s << "deflection-x"<< Util::_8s << "deflection-z" << std::endl;
    }

    if (dom.idx_out==99)
    {
       //Position in the middle of this elastic beam
        for (size_t ip=0;ip<dat.middlePar.Size();ip++)
        {
            x1 = (dat.middlePar[ip]->x(0)-dat.xc)/dat.Lz;
            z1 = dat.middlePar[ip]->x(2)/dat.Lz;
            dat.oss_ss1 << Util::_10_6 << x1 << Util::_8s << z1 << std::endl;
        } 
    }

    else
    {
        // calculate Re 
        for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
        for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
        for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
        {
            vel += dom.LBMDOM.Vel[0][ix][iy][iz](0);
        }
        vel2 = vel/(dom.LBMDOM.Ndim(0)*dom.LBMDOM.Ndim(1)*dom.LBMDOM.Ndim(2));
        // Reynolds number based on the maximum of the parabolic inlet
        // profile and the flap width used by the MDPI paper (20 mm).
        double Re   = dat.umax*dat.Ly/dat.nu;

        // calculate CD
        for (size_t ip=0; ip < dom.MPMDOM.Corners.Size(); ip++)
        {
            force = force + dom.MPMDOM.Corners[ip]->h(0);
        }
        double CD = force/(0.5*dat.rho*vel2*vel2*dat.Ly*dat.Lz);

        // calculate deflection
        for (size_t ip=0;ip<dat.EndBeamPar.Size();ip++)
        {
            x = x+dat.EndBeamPar[ip]->x(0)-dat.x0[ip];
            z = z+dat.EndBeamPar[ip]->x(2)-dat.x2[ip];
        }
        x /= dat.EndBeamPar.Size();
        z /= dat.EndBeamPar.Size();
        xf = x/dat.Lz;
        zf = z/dat.Lz;
        dat.oss_ss2 << Util::_10_6 << dom.Time << Util::_8s << Re << Util::_8s << CD << Util::_8s << xf << Util::_8s << zf << std::endl;
    }           
}

int main(int argc, char **argv) try
{
    //Number of cores
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);

    //Properties of LBM
    // Cubic tank.  Optional argv[2] sets Ldom; default 0.3 m -> n=200 at
    // dx=1.5 mm.  n=133 is 0.2 m, n=334 is 0.5 m.
    double Ldom = 0.3;
    if (argc>2) Ldom = atof(argv[2]);
    double nu = 8.20e-4;
    double dx = 0.0015;
    size_t n  = size_t(std::round(Ldom/dx));
    size_t nx = n;
    size_t ny = n;
    size_t nz = n;
    double xc = 0.5*nx*dx;

    //Properties of MPM
    // Refined slab mesh: nominal edge Lz/ndiv = 1.5625 mm.  The generated
    // surface triangles have min edge 1.937 mm, still > the 1.5 mm LBM cell.
    size_t ndiv = 64; //number of divisions along the longest dimension
    // Optional argv[6] refines the body mesh; argv[7] is the MPM grid factor
    // beta (default 1.5).  Keep the defaults for the 5 mm half-thickness slab:
    // ndiv=80/beta=0.75 and beta=0.5 both produced huge corner forces and
    // aborted immediately, so the coarser default is the stable (but
    // thickness-under-resolved) choice.
    if (argc>6) ndiv = atoi(argv[6]);
    double Nu   = 0.3;
    double  E   = 1.23e6;
    double  K   = E/(3.0*(1.0-2.0*Nu));
    double  G   = E/(2.0*(1.0+Nu));
    double rho  = 1030.0;   
    double Lx   = 0.01;
    // Optional argv[5] sets the slab thickness in x (default 0.01 m).
    // Use 0.005 m for the half-thickness run with the same inlet velocity.
    if (argc>5) Lx = atof(argv[5]);
    double Ly   = 0.02;
    double Lz   = 0.1;
    double Cs   = sqrt(E/rho);
    double h    = Lz/ndiv;
    double dt   = 1.0*h/Cs;
    double Dx   = 2.0*Lz/ndiv;
    double Bc   = Dx;

    LBMMPM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;

    dat.Ly = Ly;
    dat.Lz = Lz;
    dat.xc = xc;
    dat.nu = nu;

    dom.MPMDOM.AddRectangularBeamMesh(-1, Vec3_t(0.5*nx*dx-0.5*Lx,0.5*ny*dx-0.5*Ly,-Bc), Vec3_t(0.5*nx*dx+0.5*Lx,0.5*ny*dx+0.5*Ly,Lz), rho, ndiv);
    double beta = 1.5;
    if (argc>7) beta = atof(argv[7]);
    // Keep the MPM box padding independent of the slab thickness: 10*Lx
    // gives only 5 cm for the half-thickness slab and truncates its transient.
    double PadX = Lz;
    dom.MPMDOM.ResizeDomainMesh(Vec3_t(0.5*nx*dx-PadX,0.5*ny*dx-2.0*Ly,-1.0*Lz), Vec3_t(0.5*nx*dx+PadX,0.5*ny*dx+2.0*Ly,2.0*Lz),beta);
    
    //Setting properties for the material points
    for (size_t ip=0; ip < dom.MPMDOM.Particles.Size(); ip++)
    {
        dom.MPMDOM.Particles[ip]->K = K;
        dom.MPMDOM.Particles[ip]->G = G;
        if (dom.MPMDOM.Particles[ip]->x(2)<0.0)
        {
            dom.MPMDOM.Particles[ip]->FixVeloc();
            dom.MPMDOM.Particles[ip]->Tag = -2;
        }
        if (dom.MPMDOM.Particles[ip]->x(2)>Lz*(ndiv-0.5)/ndiv)
        {
            dom.MPMDOM.Particles[ip]->Tag = -3;
            dat.EndBeamPar.Push(dom.MPMDOM.Particles[ip]);
            dat.x0.Push(dom.MPMDOM.Particles[ip]->x(0));
            dat.x2.Push(dom.MPMDOM.Particles[ip]->x(2));
        }
        //if (dom.MPMDOM.Particles[ip]->x(0) > (1/2*nx*dx-Lx/10) && dom.MPMDOM.Particles[ip]->x(0) < (1/2*nx*dx+Lx/10))
        if (std::abs(dom.MPMDOM.Particles[ip]->x(0)-xc) < 0.001)
        {
            dom.MPMDOM.Particles[ip]->Tag = -4;
            dat.middlePar.Push(dom.MPMDOM.Particles[ip]);
        }
    }

    //Fixing positions and forces over corner points 
    for (size_t ic=0; ic < dom.MPMDOM.Corners.Size(); ic++)
    {
        if ((*dom.MPMDOM.Corners[ic]->x)(2)<0.0)
        {
            dom.MPMDOM.Corners[ic]->FixVeloc();
        }
    }

    //Setting properties for nodes
    for (size_t in=0; in < dom.MPMDOM.Nnodes; in++)
    {
        Vec3_t xn;
        dom.MPMDOM.NodePosition(in,xn);
        if (xn(2)<0.0)
        {
            dom.MPMDOM.Nodes[in].FixVeloc();
        }
    }

    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,1220.0/*rho*/,v);
        if (iz==0)
        {
            dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
        }
    }  

    dat.Vel = new double[nz];
#ifdef USE_CUDA
    dat.dVel = NULL;
#endif
    // Uniform inlet, as in the MDPI paper: Re = v*w/nu with w=Ly.
    // Current target: v=1.815 m/s -> Re=44.27 and a ~39.6 mm centreline
    // tip deflection at Ldom=0.3 m.  Optional argv[3] overrides v (m/s).
    double v    = 1.815;
    if (argc>3) v = atof(argv[3]);
    dat.umax = v;
    for (size_t i=0;i<nz;i++)
    {
        dat.Vel[i] = v;
    }
    dat.rho = 1220.0;
#ifdef USE_CUDA
    cudaMalloc(&dat.dVel,nz*sizeof(double));
    cudaMemcpy(dat.dVel,dat.Vel,nz*sizeof(double),cudaMemcpyHostToDevice);
#endif
    
    double Tf = 1.5; // Final Time
    dat.Tf = Tf;
    dat.fluidforce = Vec3_t(5.0e-3*dx/dt/Tf,0.0,0.0);
    
    dom.MPMDOM.Gn = 1.0e0;
    bool RenderVideo = true;
    if (argc>4) RenderVideo = atoi(argv[4])!=0;
    dom.Solve(Tf,Tf/100,Setup,Report,"tlbmmpm_cu_01",RenderVideo,Nproc);
}
MECHSYS_CATCH

