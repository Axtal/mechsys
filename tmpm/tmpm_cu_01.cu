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
// Elastic beam static bending


// MechSys
#include <mechsys/mpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

struct UserData
{
    bool                      cpdi; //Use the CPDI body-mesh path
    thrust::device_vector<size_t> ForcePar;
    thrust::device_vector<size_t> ForceCor;
    Array<MPM::Particle *> PlaneBeamPar0; //Particles at each plane of the beam
    Array<double >         x0;
    Array<double >         x2; 
    double                       Tf; //Final simulation time
    size_t                     ndiv;
    size_t *                 pnFpar;
    size_t *                 pnFcor;
    double                       Lx;
    double                       Ly;
    double                       Lz;
    double                       E;
    real                         Fb; //Maximun applied body force
    real *                       pFb;
    std::ofstream                oss_ss; // The output file
};

__global__ void SetForce(size_t * ForcePar, size_t * nForcePar, MPM::ParticleCU * Par, real * Fb, MPM::mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=nForcePar[0]) return;
    size_t ip = ForcePar[ic];
    real F = fmin(10.0*Fb[0]*mpmaux[0].Time/mpmaux[0].Tf,Fb[0])/nForcePar[0];
    Par[ip].b = make_real3(0.0,0.0,-F);
}

__global__ void SetCornerForce(size_t * ForceCor, size_t * nForceCor, MPM::CornerCU * Cor, real * Fb, MPM::mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=nForceCor[0]) return;
    size_t iCor = ForceCor[ic];
    real F = fmin(10.0*Fb[0]*mpmaux[0].Time/mpmaux[0].Tf,Fb[0])/nForceCor[0];
    Cor[iCor].h = make_real3(0.0,0.0,-F);
}

void Setup (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.cpdi)
    {
        if (dat.ForceCor.size()>0)
            SetCornerForce<<<dat.ForceCor.size()/dom.Nthread+1,dom.Nthread>>>(thrust::raw_pointer_cast(dat.ForceCor.data()),dat.pnFcor,dom.pCornersCU,dat.pFb,dom.pmpmaux);
    }
    else
    {
        if (dat.ForcePar.size()>0)
            SetForce<<<dat.ForcePar.size()/dom.Nthread+1,dom.Nthread>>>(thrust::raw_pointer_cast(dat.ForcePar.data()),dat.pnFpar,dom.pParticlesCU,dat.pFb,dom.pmpmaux);
    }
}

void Report (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD)); 

    double Dx   = 0.5*dat.Lx/dat.ndiv;
    

    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("tmpm01_Comparison.res");
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "y(mpm)" << "\n";
    }


    if (dom.idx_out==99)
    {
        for (size_t i=0;i<dat.x0.Size();i++)
        {
            double x = dat.x0[i]-Dx;
            double z = dat.x2[i]-dat.PlaneBeamPar0[i]->x(2);
            double yt =  2.0*dat.Fb*pow(x,2)*(3.0*(dat.Lx-Dx)-x)/(dat.E*pow(dat.Lz,3)*dat.Ly);
            dat.oss_ss << Util::_8s << x << Util::_8s << yt <<Util::_8s << z<< "\n";

        }
    }
}

int main(int argc, char **argv) try
{

    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
   
    size_t Nproc = 0.75*omp_get_max_threads();
    bool   cpdi  = false; // standard MPM by default; pass cpdi to enable CPDI
    for (int i=1; i<argc; ++i)
    {
        std::string arg(argv[i]);
        if (arg=="cpdi" || arg=="--cpdi") cpdi = true;
        else if (arg=="mpm" || arg=="noncpdi" || arg=="nocpdi" || arg=="--mpm") cpdi = false;
        else Nproc = atoi(argv[i]);
    }
    dat.cpdi = cpdi;
    std::cout << "  MPM mode                          =  " << (cpdi ? "CPDI" : "non-CPDI") << std::endl;
    //size_t ndiv = 22; //number of divisions per x lenght
    size_t ndiv = cpdi ? 125 : 50; // CPDI uses the refined body mesh; keep the non-CPDI baseline
    double Nu	= 0.3;
    double	E	= 1.23e6;
    double  K   = E/(3.0*(1.0-2.0*Nu));
    double	G	= E/(2.0*(1.0+Nu));
    double rho  = 1030.0;
    double Cs	= sqrt(E/rho);
    double Lx   = 1.0e-1;   
    double Ly   = 0.2e-1;
    double Lz   = 0.2e-1;
    double h    = Lx/ndiv;
    double dt   = 1.0*h/Cs;
    double Dx   = 2.0*Lx/ndiv;
    double Bc   = Dx;

    dat.ndiv    = ndiv;
    dat.Lx      = Lx;
    dat.Ly      = Ly;
    dat.Lz      = Lz;
    dat.E       = E;

    if (cpdi)
    {
        dom.AddRectangularBeamMesh(-1, Vec3_t(0.0,0.0,0.0), Vec3_t(Lx,Ly,Lz), rho, ndiv);
        dom.ResizeDomainMesh(Vec3_t(-0.5*Lx,-0.5*Lx,-0.5*Lx),Vec3_t(1.5*Lx,1.5*Lx,1.5*Lx), 1.0);
    }
    else
    {
        dom.AddRectangularBeam(-1, Vec3_t(0.0,0.0,0.0), Vec3_t(Lx,Ly,Lz), rho, ndiv);
        dom.ResizeDomain(Vec3_t(-0.5*Lx,-0.5*Lx,-0.5*Lx),Vec3_t(1.5*Lx,1.5*Lx,1.5*Lx), Dx);
    }

    //Setting properties for the material points
    for (size_t ip=0; ip < dom.Particles.Size(); ip++)
    {
        dom.Particles[ip]->K = K;
        dom.Particles[ip]->G = G;
        if (dom.Particles[ip]->x(0)<Bc)
        {
            dom.Particles[ip]->FixVeloc();
            dom.Particles[ip]->Tag = -2;
        }
        if (!cpdi && dom.Particles[ip]->x(0)>Lx*(ndiv-1.0)/ndiv)
        {
            dom.Particles[ip]->Tag = -3;
            dat.ForcePar.push_back(ip);
        }
        if (dom.Particles[ip]->x(2)<0.9*Lz && dom.Particles[ip]->x(2)>0.1*Lz)
        {
            dat.PlaneBeamPar0.Push(dom.Particles[ip]);
            dat.x0.Push(dom.Particles[ip]->x(0));
            dat.x2.Push(dom.Particles[ip]->x(2));
        }
    }
    // CPDI: the last cross-section of the rectangular beam is an external
    // face, so every corner in the last mesh slice is already a boundary
    // corner.  Selecting by x alone is enough here; no boundary-facet mask
    // is needed.  Corners outside this slice are never loaded, and the base
    // corners below Bc are fixed.
    if (cpdi)
    {
        for (size_t ic=0; ic < dom.Corners.Size(); ic++)
        {
            if ((*dom.Corners[ic]->x)(0)<Bc)
            {
                dom.Corners[ic]->FixVeloc();
            }
            if ((*dom.Corners[ic]->x)(0)>Lx*(ndiv-0.1)/ndiv)
            {
                dat.ForceCor.push_back(ic);
            }
        }
    }

    size_t nforpar = dat.ForcePar.size();
    size_t nforcor = dat.ForceCor.size();
    std::cout << "  Load points: particles = " << nforpar
              << ", boundary corners = " << nforcor << std::endl;
    dat.Fb      = 5.0e-1;  // Force applied at the end (reduced from 1.0e-1)
    cudaMalloc(&dat.pFb, sizeof(real));
    cudaMemcpy(dat.pFb, &dat.Fb, sizeof(real), cudaMemcpyHostToDevice);
    cudaMalloc(&dat.pnFpar, sizeof(size_t));
    cudaMemcpy(dat.pnFpar, &nforpar, sizeof(size_t), cudaMemcpyHostToDevice);
    cudaMalloc(&dat.pnFcor, sizeof(size_t));
    cudaMemcpy(dat.pnFcor, &nforcor, sizeof(size_t), cudaMemcpyHostToDevice);

    //Setting properties for nodes
    for (size_t in=0; in < dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        if (xn(0)<Bc)
        {
            dom.Nodes[in].FixVeloc();
        }
    }

    dom.Gn = 1.0e0;
    double Tf = 5.0; // Final Time
    dat.Tf = Tf;     
    dom.Solve(Tf,dt,Tf/100,Setup,Report,"tmpm01",true,Nproc);

    return 0;
}
MECHSYS_CATCH

