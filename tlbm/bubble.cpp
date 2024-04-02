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

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/lbm/Domain.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    if (argc==2) Nproc=atoi(argv[1]);
    size_t nx = 100;
    size_t ny = 100;
    size_t nz = 100;
    //size_t nz = 1;
    double nu = 1.0/6.0;
    double dx = 1.0;
    double dt = 1.0;
    double Tf = 10000.0;
    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    //LBM::Domain Dom(D3Q19, nu, iVec3_t(nx,ny,nz), dx, dt);
    //LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,nz), dx, dt);
    //Dom.Lat[0].Rhoref =  2.0;
    //Dom.Lat[0].G      = -4.0;
    //Dom.Lat[0].Gs     =  0.0;
    Dom.Lat[0].G      = -200.0;
    //Dom.Lat[0].Gs     = -100.0;
    
    //for (size_t i=0;i<nx;i++)
    //for (size_t j=0;j<ny;j++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(i,0   ,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,ny-1,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,j,0   ))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,j,ny-1))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(0   ,i,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(ny-1,i,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    //}

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
	for (size_t k=0; k<nz; ++k)
	{
		
		double rho0 = (200.0 +(1.0*rand())/RAND_MAX)*dx*dx;
		//double rho0 = (400.0 +(2.0*rand())/RAND_MAX)*dx*dx;
		//double rho0 = (2.0 +(0.02*rand())/RAND_MAX)*dx*dx;
		Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize (rho0, OrthoSys::O);
	}

	// Solve
    Dom.Solve(Tf,50.0,NULL,NULL,"bubble",true,Nproc);
}
MECHSYS_CATCH
