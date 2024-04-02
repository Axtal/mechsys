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
#include <mechsys/emlbm2/Domain.h>
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
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(dom.Nproc)
#endif
    for (int j=0;j<dat.ny;j++)
    for (int i=0;i<dat.nx;i++)
    {
        dom.Lat.GetCell(iVec3_t(i,j,dat.nz-1))->Initialize(0.0,OrthoSys::O,OrthoSys::O,OrthoSys::O);
        for (int k=0;k<dat.nz;k++)
        {
            //double dx=i-dat.nx/2,dy=j-dat.ny/2;
            //dom.Lat.GetCell(iVec3_t(i,j,k))->Jf = OrthoSys::e2*dat.J0*sin(dat.w*dom.Time)*exp(-0.75*(dx*dx+dy*dy))*(tanh(double(k)-2.0*dat.nz/5.0)+tanh(3.0*dat.nz/5.0-double(k)));
            double dz=k;
            dom.Lat.GetCell(iVec3_t(i,j,k))->Jf = OrthoSys::e0*dat.J0*sin(dat.w*dom.Time)*exp(-0.75*(dz*dz));
        }
    }
}
int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    int nx = 100;
    int ny = 100;
    int nz = 100;
    EMLBM::Domain Dom(iVec3_t(nx,ny,nz), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    
    double lambda = 20.0;//wavelenght in LBM cells.
    dat.w  = 2*M_PI/(2.0*lambda);
    dat.J0 = 1.0e-4;
    Dom.Step = 1;
    dat.nx = nx;
    dat.ny = ny;
    dat.nz = nz;
    //double Sig0 = 1.0e6;
    for (int i=0;i<dat.nx;i++)
    for (int j=0;j<dat.ny;j++)
    for (int k=0;k<dat.nz;k++)
    {
        Dom.Lat.GetCell(iVec3_t(i,j,k))->Eps = 1.0;
        //double dx=i-nx/2,dy=j-ny/2;
        //Dom.Lat.GetCell(iVec3_t(i,j,k))->Sig = Sig0*exp(-0.75*(dx*dx+dy*dy))*(tanh(double(k)-2.0*nz/5.0)+tanh(3.0*nz/5.0-double(k)));
    }
    Dom.Solve(200.0,2.0,&Setup,NULL,"temlbm03",true,nproc);
    
    return 0;
}
MECHSYS_CATCH

