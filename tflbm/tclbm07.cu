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

// Phase Field for Ice equation single bubble freezing

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/flbm/Domain.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    Array<double> nu(3);
    nu[0] = 0.1666; // This is the fluid viscosity
    nu[1] = 0.1666; // This is the mobility
    nu[2] = 0.1666; // This is the thermal conductivity


    size_t nx = 400, ny = 400;
    //size_t nx = 80, ny = 80;

    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0, PhaseFieldIce);

    Dom.rho[0] = 0.5;
    Dom.rho[1] = 1.0;
    Dom.rho[2] = 1.0;
    Dom.cap[0] = 1.0;
    Dom.cap[1] = 1.0;
    Dom.cap[2] = 1.0;
    Dom.kap[0] = 0.1666;
    Dom.kap[1] = 0.1666;
    Dom.kap[2] = 0.1666;
    Dom.thick  = 5.0;
    Dom.sigma  = 0.1;
    Dom.Ts     = 0.5;
    Dom.Tl     = 1.5;
    Dom.L      = 1.0;

    double obsX = nx/2, obsY = ny/2;
    double Rext = nx/10;
    double Rint = nx/10;
    double Hi   = 0.3;
    //double Hi   = 3.0;
    double Hl   = 3.0;


	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
        double r     = sqrt(pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0));
        //double r     = sqrt(pow((int)(j)-obsY,2.0));

        double smear = 0.5*(1.0-tanh(2.0*(r-Rext)/Dom.thick));
        //double pre   = 1.0;
        double pre   = 1.0*smear + (1.0-Dom.sigma/Rext)*(1.0-smear);
        double phi   = smear;
        smear        = 0.5*(1.0-tanh(2.0*(r-Rint)/Dom.thick));
        double H     = smear*Hi + (1.0-smear)*Hl;
        Dom.Initialize(iVec3_t(i,j,0),pre,H,phi,OrthoSys::O);
    }

    //Dom.WriteXDMF("tclbm07");        

    Dom.Solve(1.0e5,1.0e3,NULL,NULL,"tclbm07",true,1);


    return 0;
}
MECHSYS_CATCH
