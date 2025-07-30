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

double ComputeZb(double x, double y)
{
    return 0.8 * exp(-5.0 * pow(x - 1.0, 2) - 50.0 * pow(y - 0.5, 2));
}

int main(int argc, char **argv) try
{
    size_t nx = 400, ny = 200;
    double nu = 0.16;
    double dx = 5.e-3;
    double dt = 2.e-5;

    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt, ShallowWater);

    Dom.g       = 9.8;
    Dom.Sc      = 0.3;
    

    double dh = 0.1;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
        double x = i * dx;
        double y = j * dx;
        double zb = ComputeZb(x, y);

        double h = 1.0 - zb + dh*exp(-(i*dx-0.4)*(i*dx-0.4)/(0.01*0.01));

        Dom.InitializeSW(iVec3_t(i,j,0),h,OrthoSys::O);
        Dom.BForce[0][i][j][0] = OrthoSys::O;
        for (size_t k=1;k<9;k++)
        {
            size_t ni = (size_t)((int)(i)+ (int)Dom.C[k](0) + (int)(nx))%nx;
            size_t nj = (size_t)((int)(j)+ (int)Dom.C[k](1) + (int)(ny))%ny;
            double xn = ni * dx;
            double yn = nj * dx;
            double zbn = ComputeZb(xn, yn);
            // double grad_zb = zbn - zb;
            Dom.BForce[0][i][j][0] = Dom.BForce[0][i][j][0] + 3.0 * Dom.W[k] * zbn * Dom.C[k] / dx;
        }
    }
    
    
    Dom.Solve(10.0,1.0e-2,NULL,NULL,"tclbm09",true,1);

    return 0;
}
MECHSYS_CATCH
