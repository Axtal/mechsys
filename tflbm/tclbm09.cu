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

// Shallow Water Solver

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/flbm/Domain.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    size_t nx = 400, ny = 400;
    double nu = 0.16;
    double dx = 1.0;
    double dt = 1.0;

    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt, ShallowWater);

    double obsX = nx/2, obsY = ny/2;
    double Rext = nx/10;
    double h0   = 0.1;
    Dom.g       = 0.1;
    


	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
        double r     = sqrt(pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0));
        //double r     = sqrt(pow((int)(j)-obsY,2.0));
        //double r     = std::max(fabs((int)(i)-obsX),fabs((int)(j)-obsY));
        double h = 1.0 - h0*exp(-r*r/(Rext*Rext));
        //if (i==nx/10) h += 0.5;

        Dom.InitializeSW(iVec3_t(i,j,0),h,OrthoSys::O);
        Dom.BForce[0][i][j][0] = OrthoSys::O;
        for (size_t k=1;k<9;k++)
        {
            size_t ni = (size_t)((int)i + (int)Dom.C[k](0) + (int)nx)%nx;
            size_t nj = (size_t)((int)j + (int)Dom.C[k](1) + (int)ny)%ny;
            double nr = sqrt(pow((int)(ni)-obsX,2.0) + pow((int)(nj)-obsY,2.0));
            double nh = h0*exp(-nr*nr/(Rext*Rext));
            Dom.BForce[0][i][j][0] = Dom.BForce[0][i][j][0] + 3.0*Dom.W[k]*nh*Dom.C[k]/dx;
        }
    }

    Dom.Solve(1.0e5,1.0e3,NULL,NULL,"tclbm09",true,1);


    return 0;
}
MECHSYS_CATCH
