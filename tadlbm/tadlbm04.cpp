/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2015 Sergio Galindo                                    *
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

// 3D diffusion

//STD
#include<iostream>

// MechSys
#include <mechsys/adlbm/Domain.h>

struct UserData
{
    Array<Cell *> xmin;
    Array<Cell *> xmax;
    Array<Cell *> Source;
    double        rhomin;
    double        rhomax;
};

void Setup (ADLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0;i<dat.xmin.Size();i++)
    {
        Cell * c = dat.xmin[i];
        if(c->IsSolid) continue;
        c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*dat.rhomin);
        c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2] +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*dat.rhomin);
        c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*dat.rhomin);
        c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*dat.rhomin);
        c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*dat.rhomin);
        c->CalcProp();
    }
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0;i<dat.xmax.Size();i++)
    {
        Cell * c = dat.xmax[i];
        if(c->IsSolid) continue;
        c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-dat.rhomax));
        c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*dat.rhomax);
        c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*dat.rhomax) ;
        c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*dat.rhomax);
        c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*dat.rhomax);
        c->CalcProp();
    }
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
	for (size_t i=0; i<dat.Source.Size(); ++i)
    {
		Cell * c = dat.Source[i];
        for (size_t k=0;k<c->Nneigh;k++)
        {
            c->G[k] += -c->W[k]*0.01;
        }
    }
}

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    size_t nx = 50;                       // Dimension in x
    size_t ny = 50;                       // Dimension in y
    size_t nz = 50;                       // Dimension in x
    double nu     = 0.16;                 // viscocity
    double dif    = 0.001;                 // Overall diffusion coefficient
    
    ADLBM::Domain Dom(D3Q15,nu, dif, iVec3_t(nx,ny,nz), 1.0, 1.0); // Declare the LBM domain
    UserData dat;
    Dom.UserData = &dat;

	int obsX   = nx/2;   // x position
	int obsY   = ny/2; // y position
	int obsZ   = nz/2; // y position
    int radius = nz/10;

    double rho0 = 1.0;
    for (size_t i=0;i<Dom.Lat.Ncells;i++)
    {
        Dom.Lat.Cells[i]->Initialize(rho0, rho0, OrthoSys::O);
    }
    

    dat.rhomin = 1.01;
    dat.rhomax = 1.0;

    for (size_t k=0;k<nz;k++)
    for (size_t j=0;j<ny;j++)
    {
        dat.xmin.Push(Dom.Lat.GetCell(iVec3_t(0   ,j,k)));
        dat.xmax.Push(Dom.Lat.GetCell(iVec3_t(nx-1,j,k)));
        for (size_t i=0;i<nx;i++)
        {
            if ((i-obsX)*(i-obsX)+(j-obsY)*(j-obsY)+(k-obsZ)*(k-obsZ)<radius*radius)
            {
                dat.Source.Push(Dom.Lat.GetCell(iVec3_t(i,j,k)));
            }
        }
    }

    //Solving
    Dom.Time = 0.0;
    Dom.Solve(40000.0,80.0,Setup,NULL,"tadlbm04",true,nproc);
 
}
MECHSYS_CATCH

