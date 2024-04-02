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
// Drag coefficient of sphere.


// MechSys
#include <mechsys/lbmdem/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    real3              * pacc;
    double                 nu;
    double                  R;
};

__global__ void Setup(real3 * BForce, real const * Rho, real3 const * acc, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;
    BForce[ic] = Rho[ic]*acc[0];
}

void Setup (LBMDEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    Setup<<<dom.LBMDOM.Ncells/dom.Nthread+1,dom.Nthread>>>(dom.LBMDOM.pBForce,dom.LBMDOM.pRho,dat.pacc,dom.LBMDOM.plbmaux);
}

void Report (LBMDEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_force.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "Tx" << Util::_8s << "Ty" << Util::_8s << "Tz" << Util::_8s << "Vx" << Util::_8s << "Fx" << Util::_8s << "F" << Util::_8s << "Rho" << Util::_8s << "Re" << Util::_8s << "CD" << Util::_8s << "CDsphere \n";
    }
    if (!dom.Finished)
    {
        double M    = 0.0;
        double Vx   = 0.0;
        size_t nc   = 0;
        Vec3_t Flux = OrthoSys::O;
        for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
        for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
        for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
        {
            if (dom.LBMDOM.IsSolid[0][ix][iy][iz]||dom.LBMDOM.Gamma[ix][iy][iz]>1.0e-8) continue;
            Vx   += (1.0-dom.LBMDOM.Gamma[ix][iy][iz])*dom.LBMDOM.Vel[0][ix][iy][iz](0);
            Flux += dom.LBMDOM.Rho[0][ix][iy][iz]*dom.LBMDOM.Vel[0][ix][iy][iz];
            M    += dom.LBMDOM.Rho[0][ix][iy][iz];
            nc++;
        }
        Vx  /=dom.LBMDOM.Ncells;
        Flux/=M;
        M   /=nc;
        double CD  = 2.0*dom.DEMDOM.Particles[0]->F(0)/(Flux(0)*Flux(0)/M*M_PI*dat.R*dat.R);
        double Re  = 2.0*Flux(0)*dat.R/(M*dat.nu);
        double CDt = 24.0/Re + 6.0/(1.0+sqrt(Re)) + 0.4;
        if (Flux(0)<1.0e-12) 
        {
            Flux = OrthoSys::O;
            CD = 0.0;
            Re = 0.0;
        }

        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dom.DEMDOM.Particles[0]->F(0) << Util::_8s << dom.DEMDOM.Particles[0]->F(1) << Util::_8s << dom.DEMDOM.Particles[0]->F(2) << Util::_8s << dom.DEMDOM.Particles[0]->T(0) << Util::_8s << dom.DEMDOM.Particles[0]->T(1) << Util::_8s << dom.DEMDOM.Particles[0]->T(2) << Util::_8s << Vx << Util::_8s << Flux(0) << Util::_8s << norm(Flux) << Util::_8s << M << Util::_8s << Re << Util::_8s << CD << Util::_8s << CDt << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

int main(int argc, char **argv) try
{
    size_t Nproc = 0.75*omp_get_max_threads();

        
    double nu = 1.0;
    size_t nx = 241;
    size_t ny = 61;
    size_t nz = 61;
    double dx = 0.4;
    double dt = 1.6e-2;
    double R = 3.6;
    LBMDEM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;
    dat.R  = R;
    dat.nu = nu;
    
    dom.DEMDOM.AddSphere(-1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.5*dx*nz),R,1.0);
    dom.DEMDOM.GetParticle(-1)->FixVeloc();


    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,1.0/*rho*/,v);
    }   

    real3 acc = make_real3(1.5e-2,0.0,0.0);
    cudaMalloc(&dat.pacc, sizeof(real3));
    cudaMemcpy(dat.pacc, &acc, sizeof(real3), cudaMemcpyHostToDevice);

    dom.Alpha = 2.0*dx;

    double Tf = 1.0e4;
    dom.Solve(Tf,Tf/200,Setup,Report,"tlbmdem_cu_01",true,Nproc);
}
MECHSYS_CATCH

