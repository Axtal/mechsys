/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Torres                                     *
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

//Multicomponent bubble interacting with solid surface

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/flbm/Domain.h>
#include <mechsys/util/numstreams.h>

using std::cout;
using std::endl;


int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    double dt    = 1.0;
    double nua   = 0.16;
    double Gmix  = 1.0;
    double G0    =-0.2;
    double G1    = 0.2;
    double R     = 25.0;
    double rho0  = 2.0;
    double rho1  = 2.0;
    double Tf    = 100000.0;
    double dtout = 1000.0;
    if (argc>=2) Nproc  =atoi(argv[ 1]);
    if (argc>=3)
    {
        dt     =atof(argv[ 2]);
        nua    =atof(argv[ 3]);
        Gmix   =atof(argv[ 4]);
        G0     =atof(argv[ 5]);
        G1     =atof(argv[ 6]);
        R      =atof(argv[ 7]);
        rho0   =atof(argv[ 8]);
        rho1   =atof(argv[ 9]);
        Tf     =atof(argv[10]);
        dtout  =atof(argv[11]);
    }
    Array<double> nu(2);
    nu[0] = nua;
    nu[1] = nua;

    size_t nx = 200, ny = 200;

    // Setting top and bottom wall as solid
    FLBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, dt);
    //Dom.Sc = 0.0;
    int obsX = nx/2, obsY = ny/8;
    int radius =  (int) R;

	for (size_t i=0; i<nx; ++i)
    {
        Dom.IsSolid[0][i][0   ][0] = true;
        Dom.IsSolid[1][i][0   ][0] = true;
        Dom.IsSolid[0][i][ny-1][0] = true;
        Dom.IsSolid[1][i][ny-1][0] = true;
	    for (size_t j=0; j<ny; ++j)
        {
	    	if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
	    	{
                Dom.Initialize(0,iVec3_t(i,j,0),0.999*rho0,OrthoSys::O);
                Dom.Initialize(1,iVec3_t(i,j,0),0.001*rho0,OrthoSys::O);
	    	}
	    	else
	    	{
                Dom.Initialize(0,iVec3_t(i,j,0),0.001*rho1,OrthoSys::O);
                Dom.Initialize(1,iVec3_t(i,j,0),0.999*rho1,OrthoSys::O);
	    	}
        }
    }

    // Set parameters
    Dom.Gmix  = Gmix;
    Dom.Gs[0] = G0;
    Dom.Gs[1] = G1;

    Dom.Solve(Tf,dtout,NULL,NULL,"tflbm04",true,Nproc);
    return 0;
}
MECHSYS_CATCH
