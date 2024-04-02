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

/////////////////////// Test 03 the dipole antenna

// MechSys
#include <mechsys/emlbm/Domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;

struct UserData
{
    double w;
    double J0;
    int nx;
    int ny;
    int nz;
};
void Setup (EMLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (int i=0;i<dat.nx;i++)
    for (int j=0;j<dat.ny;j++)
    for (int k=0;k<dat.nz;k++)
    {
        double dx=i-dat.nx/2,dy=j-dat.ny/2;
        dom.Lat[0].GetCell(iVec3_t(i,j,k))->J[3] = 0.5*dat.J0*sin(dat.w*dom.Time)*exp(-0.75*(dx*dx+dy*dy))*(tanh(double(k)-2.0*dat.nz/5.0)+tanh(3.0*dat.nz/5.0-double(k)));
    }
}

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    int nx = 100;
    int ny = 100;
    int nz = 100;
    EMLBM::Domain Dom(D3Q7, 0.5, iVec3_t(nx,ny,nz), 0.25, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.w  = 2*M_PI*0.02;
    dat.J0 = 1.0e-3;
    Dom.Step = 1;
    dat.nx = nx;
    dat.ny = ny;
    dat.nz = nz;
    //for (int i=0;i<nx;i++)
    //for (int j=0;j<ny;j++)
    //for (int k=0;k<nz;k++)
    //{
        //double dx=i-nx/2,dy=j-ny/2;
        //Dom.Lat[0].GetCell(iVec3_t(i,j,k))->J[3] = J0*exp(-0.75*(dx*dx+dy*dy));
    //}
    //Dom.WriteXDMF("test");
    Dom.Solve(200.0,2.0,&Setup,NULL,"temlbm03",true,nproc);

    return 0;
}
MECHSYS_CATCH

