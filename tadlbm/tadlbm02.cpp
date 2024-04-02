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

// Heat Diffusion

//STD
#include<iostream>

// MechSys
#include <mechsys/adlbm/Domain.h>

struct UserData
{
    Array<Cell *> Left;
    Array<Cell *> Right;
    double        rho;
    double        templeft;
    double        tempright;
};

void Setup (ADLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
	for (size_t i=0; i<dat.Left.Size(); ++i)
	{
		Cell * c = dat.Left[i];
        c->Initialize(1.0,dat.templeft,c->Vel);
		//if (c->IsSolid) continue;
		//double rho = (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[3]+c->F[6]+c->F[7]))/(1.0-dat.Vel[i]);
		//c->F[1] = c->F[3] + (2.0/3.0)*rho*dat.Vel[i];
		//c->F[5] = c->F[7] + (1.0/6.0)*rho*dat.Vel[i] - 0.5*(c->F[2]-c->F[4]);
		//c->F[8] = c->F[6] + (1.0/6.0)*rho*dat.Vel[i] + 0.5*(c->F[2]-c->F[4]);
        //for (size_t k=0;k<c->Nneigh;k++)
        //{
            //c->G[k] = dom.Lat.Cells[c->Neighs[1]]->G[k];
        //}
        //c->CalcProp();
	}

	// Cells with prescribed density
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
	for (size_t i=0; i<dat.Right.Size(); ++i)
	{
		Cell * c = dat.Right[i];
        c->Initialize(1.0,dat.tempright,c->Vel);
		//if (c->IsSolid) continue;
		//double vx = -1.0 + (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[1]+c->F[5]+c->F[8]))/dat.rho;
		//c->F[3] = c->F[1] - (2.0/3.0)*dat.rho*vx; 
		//c->F[7] = c->F[5] - (1.0/6.0)*dat.rho*vx + 0.5*(c->F[2]-c->F[4]);
		//c->F[6] = c->F[8] - (1.0/6.0)*dat.rho*vx - 0.5*(c->F[2]-c->F[4]);
        //for (size_t k=0;k<c->Nneigh;k++)
        //{
            //c->G[k] = dom.Lat.Cells[c->Neighs[3]]->G[k];
        //}
        //c->CalcProp();
	}
}

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    double dif    = 10.0;
    size_t nx = 800;
    size_t ny = 100;
    ADLBM::Domain Dom(D2Q9,0.16, dif, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.templeft   = 2.0;
    dat.tempright  = 1.0;

    //std::cout << "1" << std::endl;

    //Assigning the left and right cells
    for (size_t i=0;i<ny;i++)
    {
        dat.Left .Push(Dom.Lat.GetCell(iVec3_t(0   ,i,0)));
        dat.Right.Push(Dom.Lat.GetCell(iVec3_t(nx-1,i,0)));
    }

    //Initializing values
    double rho0 = 1.0;
    for (size_t i=0;i<Dom.Lat.Ncells;i++)
    {
        Dom.Lat.Cells[i]->Initialize(rho0, dat.tempright, OrthoSys::O);
    }

    //Assigning all the domain as solid
    for (size_t i=0;i<Dom.Lat.Ncells;i++)
    {
        Dom.Lat.Cells[i]->IsSolid = true;
    }



    //Solving
    Dom.Time = 0.0;
    Dom.Solve(40000.0,80.0,Setup,NULL,"tadlbm02",true,nproc);
 
}
MECHSYS_CATCH

