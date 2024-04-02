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

// MechSys
#include <mechsys/lbm/Domain.h>

using std::cout;
using std::endl;

void Report(LBM::Domain & dom, void * UD)
{
    for (size_t i=0;i<dom.Disks.Size();i++)
    {
        std::cout << dom.Disks[i]->X << std::endl;
    }
}
int main(int argc, char **argv) try
{
    size_t Nproc  = 1; 
    size_t nx     = 200;
    size_t ny     = 200;
    double nu     = 1.0e-3;
    double dx     = 1.0;
    double dt     = 1.0;    
    double vb     = 0.2;
    double wb     = 0.0;
    double rc     = 10.0;
    double Kn     = 0.0e0;
    double Tf     = 1.0e3;
    if (argc>=2) Nproc = atoi(argv[1]);
    if (argc>=3) Tf    = atof(argv[2]);

    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt);
    Dom.Alpha    = 2.0*rc;
    //Dom.Sc       = 0.0;
    //Assigning solid boundaries at top and bottom
    //for (size_t i=0;i<nx;i++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    //}
    
    //Dom.AddDisk(0,Vec3_t(0.2*nx,0.2*ny,0.0),Vec3_t(vb,0.0,0.0),Vec3_t(0.0,0.0,wb),3.0,rc,1.0);
    Dom.AddDisk(0,Vec3_t(0.3*nx,0.5*ny,0.0),Vec3_t(vb ,0.0,0.0),Vec3_t(0.0,0.0,wb),3.0,rc,1.0);
    Dom.AddDisk(1,Vec3_t(0.6*nx,0.5*ny,0.0),Vec3_t(0.0,0.0,0.0),Vec3_t(0.0,0.0,wb),3.0,rc,1.0);

    for (size_t i=0;i<Dom.Disks.Size();i++)
    {
        Dom.Disks[i]->Kn = Kn;
        //std::cout << Dom.Disks[i]->M << std::endl;
    }

    for (size_t i=0;i<ny;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0,i,0))->IsSolid = true;
    }

    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
    }

    Dom.Solve(Tf,0.01*Tf,NULL,NULL,"magnus",true,Nproc);

}
MECHSYS_CATCH
