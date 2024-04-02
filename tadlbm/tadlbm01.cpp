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

// FLow in a pipe with obstacle

//STD
#include<iostream>

// MechSys
#include <mechsys/adlbm/Domain.h>

struct UserData
{
    Array<Cell *> Left;
    Array<Cell *> Right;
    Array<Cell *> Source;
    Array<double> Vel;
    double        vmax;
    double        rho;
};

void Setup (ADLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //for (size_t i=0;i<dat.Left.Size();i++)
    //{
        //dat.Left [i]->Initialize(dat.Left [i]->RhoBC,dat.Left [i]->VelBC);
        //dat.Right[i]->Initialize(dat.Right[i]->RhoBC,dat.Right[i]->VelBC);
    //}

	// Cells with prescribed velocity
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
	for (size_t i=0; i<dat.Left.Size(); ++i)
	{
		Cell * c = dat.Left[i];
        c->Vel   = Vec3_t(dat.Vel[i],0.0,0.0);
        c->Initialize(1.0,c->Temp,c->Vel);
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
        c->Initialize(1.0,c->Temp,c->Vel);
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
    double u_max  = 0.1;                // Poiseuille's maximum velocity
    double Re     = 40000.0;                  // Reynold's number
    size_t nx = 800;
    size_t ny = 100;
    int radius = ny/10 + 1;           // radius of inner circle (obstacle)
    double nu     = u_max*(2*radius)/Re; // viscocity
    double dif    = 0.01;
    
    ADLBM::Domain Dom(D2Q9,nu, dif, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;

    //std::cout << "1" << std::endl;

    dat.vmax = u_max;
    //Assigning the left and right cells
    for (size_t i=0;i<ny;i++)
    {
        dat.Left .Push(Dom.Lat.GetCell(iVec3_t(0   ,i,0)));
        dat.Right.Push(Dom.Lat.GetCell(iVec3_t(nx-1,i,0)));
        
        // set parabolic profile
        double L  = ny - 2;                       // channel width in cell units
        double yp = i - 1.5;                      // ordinate of cell
        double vx = dat.vmax*4/(L*L)*(L*yp - yp*yp); // horizontal velocity
        dat.Vel.Push(vx);
    }
    dat.rho  = 1.0;

	// set inner obstacle
	int obsX   = nx/2;   // x position
	int obsY   = ny/2+3; // y position
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        if ((i-obsX)*(i-obsX)+(j-obsY)*(j-obsY)<radius*radius)
        {
            Dom.Lat.GetCell(iVec3_t(i,j,0))->IsSolid = true;
            Dom.Lat.GetCell(iVec3_t(i,j,0))->IsNodif = true;
        }
    }

    //Initializing values
    double rho0 = 1.0;
    Vec3_t v0(0.08,0.0,0.0);
    for (size_t i=0;i<Dom.Lat.Ncells;i++)
    {
        Dom.Lat.Cells[i]->Initialize(rho0, rho0, v0);
    }

    //set a drop of ink in the beginning
	obsX   = nx/6;   // x position
	obsY   = ny/2+3; // y position
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        if ((i-obsX)*(i-obsX)+(j-obsY)*(j-obsY)<radius*radius)
        {
            //Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize(1.0,2.0,v0);
            dat.Source.Push(Dom.Lat.GetCell(iVec3_t(i,j,0)));
        }
    }

    //Assigning solid boundaries at top and bottom
    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat.GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat.GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }



    //std::cout << "3" << std::endl;
    //Solving
    Dom.Time = 0.0;
    Dom.Solve(40000.0,80.0,Setup,NULL,"tadlbm01",true,nproc);
 
}
MECHSYS_CATCH

