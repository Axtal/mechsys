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

// Advection diffusion equation

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/flbm/Domain.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    Array<double> nu(2);
    nu[0] = 0.1666; // This is the fluid viscosity
    nu[1] = 0.16666; // This is the diffusion coefficient


    size_t nx = 400, ny = 400;
    //size_t nx = 80, ny = 80;

    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);

    int obsX = nx/2, obsY = ny/2;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
		Vec3_t V;  V = 0.0, 0.0, 0.0;
        Dom.Initialize(0,iVec3_t(i,j,0),1.0,V);
        double r2  = pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0);
        double con = 1.0*exp(-r2/100.0);
        Dom.Initialize(1,iVec3_t(i,j,0),con,V);
    }

    Dom.Solve(1.0e5,1.0e3,NULL,NULL,"tclbm06",true,1,AdvectionDiffusion);


    return 0;
}
MECHSYS_CATCH
