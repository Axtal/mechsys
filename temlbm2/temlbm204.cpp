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

/////////////////////// Test 04 Dielectric Sphere in Uniform field

// MechSys
#include <mechsys/emlbm2/Domain.h>
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;


struct UserData
{
    Vec3_t E0;
    Vec3_t B0;
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
    for (int k=0;k<dat.nz;k++)
    {
        dom.Lat.GetCell(iVec3_t(0,j,k))->Initialize(0.0,OrthoSys::O,dat.E0,dat.B0);
    }
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(dom.Nproc)
#endif
    for (int i=0;i<dat.nx;i++)
    for (int k=0;k<dat.nz;k++)
    {
        dom.Lat.GetCell(iVec3_t(i,0,k))->Initialize(0.0,OrthoSys::O,dat.E0,dat.B0);
    }
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(dom.Nproc)
#endif
    for (int i=0;i<dat.nx;i++)
    for (int j=0;j<dat.ny;j++)
    {
        dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize(0.0,OrthoSys::O,dat.E0,dat.B0);
    }
}

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    int nx = 121;
    int ny = 121;
    int nz = 121;
    EMLBM::Domain Dom(iVec3_t(nx,ny,nz), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Step = 1;
    //double Sig0 = 1.0e6;
    double sigma = 0.1;
    Vec3_t    E0 = Vec3_t(0.0,0.0,sigma);
    //Vec3_t    B0 = Vec3_t(0.0,-sqrt(2.0)*E0(2),0.0);
    Vec3_t    B0 = OrthoSys::O;
    dat.E0 = E0;
    dat.B0 = B0;
    dat.nx = nx;
    dat.ny = ny;
    dat.nz = nz;
    int obsX = nx/2;
    int obsY = ny/2;
    int obsZ = nz/2;
    int    R = 10;


    //for (int i=0;i<nx;i++)
    //for (int j=0;j<ny;j++)
    //{
        //for (int k=0;k<nz;k++)
        //{
            //if (pow(i-obsX,2.0)+pow(j-obsY,2.0)+pow(k-obsZ,2.0)<R*R)
            //{
                //Dom.Lat.GetCell(iVec3_t(i,j,k))->Eps = 5.0;
            //}
        //}
    //}
    Vec3_t Xs(obsX,obsY,obsZ);
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(nproc)
#endif
    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Vec3_t Xc(i,j,k);
        double gamma = DEM::SphereCube(Xs,Xc,R,1.0)/12.0;
        Dom.Lat.GetCell(iVec3_t(i,j,k))->Eps = 4.0*gamma+1.0;
    }
    
    Dom.Solve(5.0e3,5.0e1,&Setup,NULL,"temlbm04",true,nproc);
    
    return 0;
}
MECHSYS_CATCH

