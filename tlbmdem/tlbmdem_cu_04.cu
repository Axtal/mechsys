/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2021 Sergio Galindo                                    *
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
// falling body into SC fluid


// MechSys
#include <mechsys/lbmdem/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    thrust::device_vector<real3> g; 
    real3 *    pg;
};

__global__ void ApplyGravity(real3 const * pg, real3 * BForce, real * Rho, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;
    BForce[ic] = Rho[ic]*pg[0];
}

void Setup (LBMDEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    ApplyGravity<<<dom.LBMDOM.Ncells*dom.LBMDOM.Nl/dom.Nthread+1,dom.Nthread>>>(dat.pg,dom.LBMDOM.pBForce,dom.LBMDOM.pRho,dom.LBMDOM.plbmaux);
}

void Report (LBMDEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
}

int main(int argc, char **argv) try
{
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>=2) Nproc = atoi(argv[1]);

        
    size_t nx  = 100;
    size_t ny  = 100;
    size_t nz  = 100;
    double dx  = 1.0;
    double dt  = 1.0;
    double R   = 10.0;
    double g   = 0.001;
    Array<double> nu(2);
    nu[0] = 0.01;
    nu[1] = 0.01;
    //double nu  = 0.1;
    LBMDEM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.LBMDOM.G [0]= -200.0;
    dom.LBMDOM.Gs[0]=   -0.0;
    dom.LBMDOM.G [1]=    0.0;
    dom.LBMDOM.Gs[1]=    0.0;
    dom.LBMDOM.Gmix =  0.001;
    dom.UserData = &dat;
    dat.g.resize(1);
    dat.g[0]     = make_real3(0.0,0.0,-g);
    dat.pg       = thrust::raw_pointer_cast(dat.g.data());
    
    double rhos = 4000.0;

    dom.DEMDOM.AddSphere(1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.7*dx*nz),R,rhos);
    //dom.DEMDOM.AddFromOBJ(1,"damshape.obj",dx,rhos,15.0);
    //dom.DEMDOM.GetParticle( 1)->Position(Vec3_t(0.5*dx*nx,0.5*dx*ny,0.6*dx*nz));
    //Quaternion_t q;
    //NormalizeRotation (M_PI,OrthoSys::e0,q);
    //dom.DEMDOM.GetParticle( 1)->Rotate(q,dom.DEMDOM.GetParticle( 1)->x);
    dom.DEMDOM.GetParticle( 1)->Ff = dom.DEMDOM.GetParticle(1)->Props.m*Vec3_t(0.0,0.0,-g);
    dom.DEMDOM.AddPlane(-1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.0),dx,nx*dx,ny*dx,rhos);
    dom.DEMDOM.GetParticle(-1)->FixVeloc();
    dom.DEMDOM.GetParticle(-1)->Bdry = true;


    double rho1 = 1200.0;
    double rho2 =  100.0;

    
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        if (iz<=nz/2)
        {
            dom.LBMDOM.Initialize(0,idx,       rho1,OrthoSys::O);
            dom.LBMDOM.Initialize(1,idx,       0.01,OrthoSys::O);
        }
        else          
        {
            dom.LBMDOM.Initialize(0,idx,       0.01,OrthoSys::O);
            dom.LBMDOM.Initialize(1,idx,       rho2,OrthoSys::O);
        }

        if (iz==0||iz==nz-1) 
        {
            dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
            dom.LBMDOM.IsSolid[1][ix][iy][iz] = true;
        }
        //if (iy==0||iy==ny-1) 
        //{
            //dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
            //dom.LBMDOM.IsSolid[1][ix][iy][iz] = true;
        //}
        //if (ix==0||ix==nx-1) 
        //{
            //dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
            //dom.LBMDOM.IsSolid[1][ix][iy][iz] = true;
        //}
    }   

    dom.Alpha = 2.0*dx;
    dom.PeriodicX= true;
    dom.PeriodicY= true;


    double Tf = 4.0e4;
    dom.Solve(Tf,Tf/400,Setup,Report,"tlbmdem_cu_04",true,Nproc);
}
MECHSYS_CATCH

