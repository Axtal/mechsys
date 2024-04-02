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
    double      * Vel;
    double        vmax;
    double        rho;
    #ifdef USE_OCL
    cl::Buffer        bBCVel;
    cl::Program       UserProgram;
    #endif
};

void Setup (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    // Cells with prescribed velocity
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
	for (size_t i=0; i<dom.Ndim(1); ++i)
	{
        double * f = dom.F[0][0][i][0];
		double rho = (f[0]+f[2]+f[4] + 2.0*(f[3]+f[6]+f[7]))/(1.0-dat.Vel[i]);
		f[1] = f[3] + (2.0/3.0)*rho*dat.Vel[i];
		f[5] = f[7] + (1.0/6.0)*rho*dat.Vel[i] - 0.5*(f[2]-f[4]);
		f[8] = f[6] + (1.0/6.0)*rho*dat.Vel[i] + 0.5*(f[2]-f[4]);
        dom.Vel[0][0][i][0] = OrthoSys::O;
        dom.Rho[0][0][i][0] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][0][i][0] +=  dom.F[0][0][i][0][k];
            dom.Vel[0][0][i][0] +=  dom.F[0][0][i][0][k]*dom.C[k];
        }
        dom.Vel[0][0][i][0] /= dom.Rho[0][0][i][0];
	}

	// Cells with prescribed density
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
	for (size_t i=0; i<dom.Ndim(1); ++i)
	{
        double * f = dom.F[0][dom.Ndim(0)-1][i][0];
		double vx = -1.0 + (f[0]+f[2]+f[4] + 2.0*(f[1]+f[5]+f[8]))/dat.rho;
		f[3] = f[1] - (2.0/3.0)*dat.rho*vx; 
		f[7] = f[5] - (1.0/6.0)*dat.rho*vx + 0.5*(f[2]-f[4]);
		f[6] = f[8] - (1.0/6.0)*dat.rho*vx - 0.5*(f[2]-f[4]);
        dom.Vel[0][dom.Ndim(0)-1][i][0] = OrthoSys::O;
        dom.Rho[0][dom.Ndim(0)-1][i][0] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][dom.Ndim(0)-1][i][0] +=  dom.F[0][dom.Ndim(0)-1][i][0][k];
            dom.Vel[0][dom.Ndim(0)-1][i][0] +=  dom.F[0][dom.Ndim(0)-1][i][0][k]*dom.C[k];
        }
        dom.Vel[0][dom.Ndim(0)-1][i][0] /= dom.Rho[0][dom.Ndim(0)-1][i][0];
	}
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    if (argc==2) Nproc=atoi(argv[1]);
    double u_max  = 0.1;                // Poiseuille's maximum velocity
    double Re     = 400.0;                  // Reynold's number
    size_t nx = 2400;
    size_t ny = 300;
    int radius = ny/10 + 1;           // radius of inner circle (obstacle)
    double nu     = u_max*(2*radius)/Re; // viscocity
    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    
    UserData dat;
    Dom.UserData = &dat;

    dat.vmax = u_max;

    dat.Vel = new double[ny];
    for (size_t i=0;i<ny;i++)
    {
        // set parabolic profile
        double L  = ny - 2;                       // channel width in cell units
        double yp = i - 1.5;                      // ordinate of cell
        double vx = dat.vmax*4/(L*L)*(L*yp - yp*yp); // horizontal velocity
        dat.Vel[i] = vx;
    }
    
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
    
    Dom.Solve(40000.0,400.0,Setup,NULL,"tflbm01",true,Nproc);


}
MECHSYS_CATCH
