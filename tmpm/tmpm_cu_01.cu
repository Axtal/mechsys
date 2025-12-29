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

struct UserData
{
    thrust::device_vector<size_t> ForcePar;
    Array<MPM::Particle *> PlaneBeamPar0; //Particles at each plane of the beam
    Array<double >         x0;
    Array<double >         x2; 
    double                       Tf; //Final simulation time
    size_t                     ndiv;
    size_t *                 pnFpar;
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
    real F = fmin(10.0*Fb[0]*mpmaux[0].Time/mpmaux[0].Tf,Fb[0]);
    Par[ip].b = make_real3(0.0,0.0,-F);
}

void Setup (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    SetForce<<<dat.ForcePar.size()/dom.Nthread+1,dom.Nthread>>>(thrust::raw_pointer_cast(dat.ForcePar.data()),dat.pnFpar,dom.pParticlesCU,dat.pFb,dom.pmpmaux);
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
            double yt =  2.0*dat.Fb*dat.ForcePar.size()*pow(x,2)*(3.0*(dat.Lx-Dx)-x)/(dat.E*pow(dat.Lz,3)*dat.Ly);
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
    if (argc>1) Nproc = atoi(argv[1]);
    //size_t ndiv = 22; //number of divisions per x lenght
    size_t ndiv = 50; //number of divisions per x lenght
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

    dom.AddRectangularBeam(-1, Vec3_t(0.0,0.0,0.0), Vec3_t(Lx,Ly,Lz), rho, ndiv);
    dom.ResizeDomain(Vec3_t(-0.5*Lx,-0.5*Lx,-0.5*Lx),Vec3_t(1.5*Lx,1.5*Lx,1.5*Lx), Dx);

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
        if (dom.Particles[ip]->x(0)>Lx*(ndiv-1.0)/ndiv)
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
    size_t nforpar = dat.ForcePar.size();
    dat.Fb      = 1.0e-1/dat.ForcePar.size();  // Force applied at the end
    cudaMalloc(&dat.pFb, sizeof(real));
    cudaMemcpy(dat.pFb, &dat.Fb, sizeof(real), cudaMemcpyHostToDevice);
    cudaMalloc(&dat.pnFpar, sizeof(size_t));
    cudaMemcpy(dat.pnFpar, &nforpar, sizeof(size_t), cudaMemcpyHostToDevice);

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

