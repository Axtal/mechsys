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

// Phase Field for Bubble

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
    nu[1] = 0.1666; // This is the mobility


    size_t nx = 400, ny = 400;
    //size_t nx = 80, ny = 80;

    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0, PhaseField);

    Dom.rho[0] = 0.5;
    Dom.rho[1] = 1.0;
    Dom.thick  = 5.0;
    Dom.sigma  = 0.1;

    double obsX = nx/2, obsY = ny/2;
    double Rext = nx/10;


	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
        double r     = sqrt(pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0));
        //double r     = sqrt(pow((int)(j)-obsY,2.0));

        double smear = 0.5*(1.0-tanh(2.0*(r-Rext)/Dom.thick));
        //double pre   = 1.0;
        double pre   = 1.0*smear + (1.0-Dom.sigma/Rext)*(1.0-smear);
        double phi   = smear;
        Dom.Initialize(iVec3_t(i,j,0),pre,phi,OrthoSys::O);
    }

    //Dom.WriteXDMF("tclbm08");        

    Dom.Solve(1.0e5,1.0e3,NULL,NULL,"tclbm08",true,1);


    return 0;
}
MECHSYS_CATCH
