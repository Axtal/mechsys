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

struct UserData
{
    Vec3_t             g;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t j=0;j<dom.Lat.Size();j++)
    for (size_t i=0;i<dom.Lat[j].Ncells;i++)
    {
        Cell * c = dom.Lat[j].Cells[i];
        c->BForcef = c->Density()*dat.g;
    }
}

int main(int argc, char **argv) try
{

    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    //Array<double> nu(2);
    //nu[0] = 0.15;
    //nu[1] = 0.15;
    double  nu = 0.15;

    size_t nx = 400, ny = 200;
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.g           = 0.0,-1.0e-5,0.0;

    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }
    for (size_t i=0;i<ny;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    }
	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
		//if ((i>nx/2.0-10)&&(i<(nx/2.0+10))&&(j<ny/5.0))
        //{
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.1,OrthoSys::O);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.1,OrthoSys::O);
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->IsSolid = true;
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->IsSolid = true;
        //}
		//else if ((i<nx/8.0)&&(i>=1)&&(j<9.0*ny/10.0)&&(j>=1)) 
		if ((i<nx/8.0)&&(i>=1)&&(j<9.0*ny/10.0)&&(j>=1)) 
        {
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(3000,OrthoSys::O);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.001,OrthoSys::O);
        }
        else
        {
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.001  ,OrthoSys::O);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(1.0,OrthoSys::O);
        }
    }

    // Set parameters
    Dom.Lat[0].G  = -200.0;
    Dom.Lat[0].Gs = -200.0;
    //Dom.Lat[1].G  =  0.0;
    //Dom.Lat[1].Gs =  0.0;
    //Dom.Gmix      =  0.001;

    Dom.Solve(5000,50.0,Setup,NULL,"dam",true,nproc);
}
MECHSYS_CATCH
