
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

// Simulation of single component bubble

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/flbm/Domain.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    if (argc==2) Nproc=atoi(argv[1]);
    size_t nx = 100;
    size_t ny = 100;
    size_t nz = 100;
    //size_t nx = 10;
    //size_t ny = 10;
    //size_t nz = 1;
    double nu = 1.0/6.0;
    //double nu = 0.001;
    double dx = 1.0;
    double dt = 1.0;
    double Tf = 10000.0;
    //double Tf = 200.0;
    FLBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    //FLBM::Domain Dom(D3Q19, nu, iVec3_t(nx,ny,nz), dx, dt);
    //FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,nz), dx, dt);
    //Dom.Rhoref[0] =  2.0;
    //Dom.G     [0] = -4.0;
    //Dom.Gs    [0] =  0.0;
    Dom.G[0]     = -200.0;
    //Dom.Sc       = 0.0;
    
	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
	for (size_t k=0; k<nz; ++k)
	{
		
		double rho0 = (200.0 +(1.0*rand())/RAND_MAX)*dx*dx;
		//double rho0 = (400.0 +(2.0*rand())/RAND_MAX)*dx*dx;
		//double rho0 = (2.0 +(0.02*rand())/RAND_MAX)*dx*dx;
		Dom.Initialize (0,iVec3_t(i,j,k),rho0, OrthoSys::O);
	}

    //Dom.SolidCube(Vec3_t (100,100,100),20);

	// Solve
    Dom.Solve(Tf,Tf/200,NULL,NULL,"tflbm02",true,Nproc);
}
MECHSYS_CATCH
