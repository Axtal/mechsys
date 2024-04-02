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

// FLow in a pipe with temperature diffusion

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
    double        templeft;
    double        tempright;
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
        c->Initialize(1.0,c->Temp,c->Vel);
        for (size_t k=0;k<c->Nneigh;k++)
        {
            c->G[k] = dom.Lat.Cells[c->Neighs[3]]->G[k];
        }
        //c->Initialize(1.0,dat.tempright,c->Vel);
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
    double u_max  = 0.1;                // Poiseuille's maximum velocity
    double Re     = 100.0;                  // Reynold's number
    size_t nx = 400;                       // Dimension in x
    size_t ny = 200;                       // Dimension in y
    size_t W  = 50;                        // Width of the pipe
    double nu     = u_max*(2.0*W)/Re; // viscocity
    double dif    = 0.01;             // Overall diffusion coefficient
    
    ADLBM::Domain Dom(D2Q9,nu, dif, iVec3_t(nx,ny,1), 1.0, 1.0); // Declare the LBM domain
    UserData dat;
    Dom.UserData = &dat;
    dat.templeft   = 2.0;                             // Temperature at the Right hand side
    dat.tempright  = 1.0;                             // Initial temperature 

    //std::cout << "1" << std::endl;

    dat.vmax = u_max;
    //Assigning the left and right cells

    for (size_t i=0;i<ny;i++)
    {
        dat.Left .Push(Dom.Lat.GetCell(iVec3_t(0   ,i,0)));
        dat.Right.Push(Dom.Lat.GetCell(iVec3_t(nx-1,i,0)));
        
        // set parabolic profile
        double yp = i-ny/2+W/2;                      // ordinate of cell
        double vx = dat.vmax*4/(W*W)*(W*yp - yp*yp); // horizontal velocity
        if (vx<0.0) vx = 0.0;
        dat.Vel.Push(vx);
    }
    //Initializing values
    double rho0 = 1.0;
    Vec3_t v0(0.08,0.0,0.0);
    for (size_t i=0;i<Dom.Lat.Ncells;i++)
    {
        Dom.Lat.Cells[i]->Initialize(rho0, dat.tempright, OrthoSys::O);
    }

    //Assigning solid boundaries at top and bottom
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        if ((j<ny/2-W/2)||(j>ny/2+W/2))
        {
            Dom.Lat.GetCell(iVec3_t(i,j,0))->IsSolid = true;
            Dom.Lat.GetCell(iVec3_t(i,j,0))->Dif = 0.1;
        }
    }



    //std::cout << "3" << std::endl;
    //Solving
    Dom.Time = 0.0;
    Dom.Solve(40000.0,80.0,Setup,NULL,"tadlbm03",true,nproc);
 
}
MECHSYS_CATCH

