/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2014 Sergio Galindo                                    *
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

/////////////////////// Test 01 the Biot Savart law

// MechSys
#include <mechsys/emlbm2/Domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    int nx = 400;
    int ny = 400;
    int nz = 10;
    EMLBM::Domain Dom(iVec3_t(nx,ny,nz), 1.0, 1.0);
    Dom.Step = 1;


    double J0 = -1.0e-3;
    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        double dx=i-nx/2,dy=j-ny/2;
        Dom.Lat.GetCell(iVec3_t(i,j,k))->Jf = OrthoSys::e2*J0*exp(-0.75*(dx*dx+dy*dy));
    }
    Dom.Solve(2000.0,20.0,NULL,NULL,"temlbm01",true,nproc);
    

    return 0;
}
MECHSYS_CATCH

