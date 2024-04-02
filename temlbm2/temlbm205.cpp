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

/////////////////////// Test 05 The conductive medium

// MechSys
#include <mechsys/emlbm2/Domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;

struct UserData
{
    double ome;
    double Tf;
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
    for (int k=0;k<dat.nz;k++)
    {
        dom.Lat.GetCell(iVec3_t(dat.nx-1,j,k))->Initialize(0.0,OrthoSys::O,OrthoSys::O,OrthoSys::O);
        for (int i=0;i<dat.nx;i++)
        {
            double dx=i;
            dom.Lat.GetCell(iVec3_t(i,j,k))->Jf = OrthoSys::e1*dat.J0*sin(dat.ome*dom.Time)*exp(-0.75*(dx*dx));
        }
    }
}

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    int nx = 2000;
    int ny = 1;
    int nz = 1;
    double Sigsol   = 10.0;
    double Tf       = 10000.0;
    double dtOut    = 50.0;
    double J0       = 1.0e-4;
    double lambda   = 200.0;
    EMLBM::Domain Dom(iVec3_t(nx,ny,nz), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Step = 1;
    dat.Tf   = Tf;
    dat.ome  = 2*M_PI/(sqrt(2.0)*lambda);
    dat.J0   = J0;
    dat.nx = nx;
    dat.ny = ny;
    dat.nz = nz;

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Dom.Lat.GetCell(iVec3_t(i,j,k))->Initialize(0.0,OrthoSys::O,OrthoSys::O,OrthoSys::O);
        Dom.Lat.GetCell(iVec3_t(i,j,k))->Sig = Sigsol*0.5*(tanh(i-nx/4)+1.0);
    }
    Dom.Solve(Tf,dtOut,&Setup,NULL,"temlbm05",true,nproc);

    return 0;
}
MECHSYS_CATCH

