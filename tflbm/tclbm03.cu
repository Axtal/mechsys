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

// Multicomponent bubble

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/flbm/Domain.h>

using std::cout;
using std::endl;

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

void Setup (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    ApplyGravity<<<dom.Ncells*dom.Nl/256+1,256>>>(dat.pg,dom.pBForce,dom.pRho,dom.plbmaux);
    cudaDeviceSynchronize();
}

int main(int argc, char **argv) try
{
    Array<double> nu(2);
    nu[0] = 1.0/6.0;
    nu[1] = 1.0/6.0;
    //nu[1] = 1.0/30.0;

    size_t nx = 100, ny = 100;
    //size_t nx = 80, ny = 80;

    // Setting top and bottom wall as solid
    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.g.resize(1);
    dat.g[0]     = make_real3(0.0,-0.0001,0.0);
    dat.pg       = thrust::raw_pointer_cast(dat.g.data());
    for (size_t i=0;i<nx;i++)
    {
        Dom.IsSolid[0][i][0   ][0] = true;
        Dom.IsSolid[0][i][ny-1][0] = true;
        Dom.IsSolid[1][i][0   ][0] = true;
        Dom.IsSolid[1][i][ny-1][0] = true;
    }

    // Set inner drop
    int obsX = nx/2, obsY = ny/2;
    int radius =  nx/8.0;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
		Vec3_t V;  V = 0.0, 0.0, 0.0;
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
		{
            Dom.Initialize(0,iVec3_t(i,j,0),1300.0,V);
            Dom.Initialize(1,iVec3_t(i,j,0),   0.1,V);
            //Dom.Initialize(0,iVec3_t(i,j,0),   0.1,V);
            //Dom.Initialize(1,iVec3_t(i,j,0), 100.0,V);
		}
		else
		{
            //Dom.Initialize(0,iVec3_t(i,j,0),1300.0,V);
            //Dom.Initialize(1,iVec3_t(i,j,0),   0.1,V);
            Dom.Initialize(0,iVec3_t(i,j,0),   0.1,V);
            Dom.Initialize(1,iVec3_t(i,j,0), 100.0,V);
		}
    }

    // Set parameters
    Dom.G [0]    = -200.0;
    Dom.Gs[0]    = -400.0;
    Dom.G [1]    =  0.0;
    Dom.Gs[1]    =  400.0;
    Dom.Gmix     =  0.001;
    String filekey("tclbm03");
    String filename(filekey+".hdf5");
    // if tclbm03.hdf5 exists then load that file to continue
    if (Util::FileExists(filename)) Dom.Load(filekey.CStr());
    Dom.Solve(1.0e4,1.0e2,Setup,NULL,filekey.CStr(),true,1);
    Dom.Save(filekey.CStr());


    return 0;
}
MECHSYS_CATCH
