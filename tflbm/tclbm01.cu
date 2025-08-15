/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

// FLow in a pipe with obstacle


// MechSys
#include <mechsys/flbm/Domain.h>

struct UserData
{
    thrust::device_vector<real> Vel; 
    real                     * pVel; 
    real                       vmax; 
    real                        rho; 
};

__global__ void Left_BC(real * VBC, bool * IsSolid, real * F, real3 * Vel, real * Rho, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ny) return;
    
    size_t ib  = ic*lbmaux[0].Nx;
    if (!IsSolid[ib])
    {
        //" Initialize(ib,F,Rho,Vel,1.0,(real3)(VBC[ic],0.0,0.0),lbmaux); \n"
        size_t iv = ib*lbmaux[0].Nneigh;
        real * f = F + iv;
	    real rho = (f[0]+f[2]+f[4] + 2.0*(f[3]+f[6]+f[7]))/(1.0-VBC[ic]);
	    f[1] = f[3] + (2.0/3.0)*rho*VBC[ic];
	    f[5] = f[7] + (1.0/6.0)*rho*VBC[ic] - 0.5*(f[2]-f[4]);
	    f[8] = f[6] + (1.0/6.0)*rho*VBC[ic] + 0.5*(f[2]-f[4]);
        Rho   [ib] = 0.0;
        Vel   [ib] = make_real3(0.0,0.0,0.0);
        for(size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ib] += F[iv + k];
            Vel[ib] = Vel[ib] + F[iv + k]*lbmaux[0].C[k];
        }
        Vel[ib] = (1.0/Rho[ib])*Vel[ib];
    }
}

__global__ void Right_BC(bool * IsSolid, real * F, real3 * Vel, real * Rho, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ny) return;

    size_t ib  = ic*lbmaux[0].Nx + lbmaux[0].Nx-1;
    if (!IsSolid[ib])
    {
        //" Initialize(ib,F,Rho,Vel,1.0,Vel[ib],lbmaux); \n"
        size_t iv  = ib*lbmaux[0].Nneigh;
        real * f = F + iv;
        real rho = 1.0;
	    real vx = -1.0 + (f[0]+f[2]+f[4] + 2.0f*(f[1]+f[5]+f[8]))/rho;
	    f[3] = f[1] - (2.0/3.0)*rho*vx; 
	    f[7] = f[5] - (1.0/6.0)*rho*vx + 0.5*(f[2]-f[4]);
	    f[6] = f[8] - (1.0/6.0)*rho*vx - 0.5*(f[2]-f[4]);
        Rho   [ib] = 0.0f;
        Vel   [ib] = make_real3(0.0,0.0,0.0);
        for(size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ib] += F[iv + k];
            Vel[ib] = Vel[ib] + F[iv + k]*lbmaux[0].C[k];
        }
        Vel[ib] = (1.0/Rho[ib])*Vel[ib];
    }
}

void Setup (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    Left_BC<<<dom.Ndim(1)/256+1,256>>>(dat.pVel,dom.pIsSolid,dom.pF,dom.pVel,dom.pRho,dom.plbmaux);
    //cudaDeviceSynchronize();
    Right_BC<<<dom.Ndim(1)/256+1,256>>>(dom.pIsSolid,dom.pF,dom.pVel,dom.pRho,dom.plbmaux);
    //cudaDeviceSynchronize();
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    double u_max  = 0.1;                // Poiseuille's maximum velocity
    double Re     = 40000.0;                  // Reynold's number
    if (argc==2) Re=atof(argv[1]);
    size_t nx = 2400;
    size_t ny = 300;
    int radius = ny/10 + 1;           // radius of inner circle (obstacle)
    double nu     = u_max*(2*radius)/Re; // viscocity
    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    
    UserData dat;
    Dom.UserData = &dat;

    dat.vmax = u_max;

    dat.Vel.resize(ny);
    for (size_t i=0;i<ny;i++)
    {
        // set parabolic profile
        double L  = ny - 2;                       // channel width in cell units
        double yp = i - 1.5;                      // ordinate of cell
        double vx = dat.vmax*4/(L*L)*(L*yp - yp*yp); // horizontal velocity
        dat.Vel[i] = vx;
    }

    dat.pVel = thrust::raw_pointer_cast(dat.Vel.data());
    
    dat.rho  = 1.0;

	// set inner obstacle
	int obsX   = ny/2;   // x position
	int obsY   = ny/2+3; // y position
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        if ((i-obsX)*(i-obsX)+(j-obsY)*(j-obsY)<radius*radius)
        {
            Dom.IsSolid[0][i][j][0] = true;
        }
    }

    //Assigning solid boundaries at top and bottom
    for (size_t i=0;i<nx;i++)
    {
        Dom.IsSolid[0][i][0][0]    = true;
        Dom.IsSolid[0][i][ny-1][0] = true;
    }


    double rho0 = 1.0;
    Vec3_t v0(0.08,0.0,0.0);

    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    {
        iVec3_t idx(ix,iy,0);
        Dom.Initialize(0,idx,rho0,v0);
    }  
    
    Dom.Solve(40000.0,400.0,Setup,NULL,"tclbm01",true,Nproc);
}
MECHSYS_CATCH
