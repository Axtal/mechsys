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
// Sinking disks

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    Vec3_t             g;
    Array<Cell *> Bottom;
    Array<Cell *>   Gate;
    bool          Closed;
    double            Tf;
    double            Kn;
    Vec3_t          Xmin;
    Vec3_t          Xmax;
    double           rho;
    double         Xgate;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c = dom.Lat[0].Cells[i];
        c->BForcef = c->Density()*dat.g;
        c = dom.Lat[1].Cells[i];
        c->BForcef = c->Density()*dat.g;
    }
    //double alpha = 1.0;
    //double TC    = 3.0*dat.Tf/4.0;
    //double rho = 500.0*(1.0 + alpha*4.0*dom.Time*(TC-dom.Time)/(TC*TC));
    //if (dom.Time>TC) rho = 500.0;
    //#pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    //for (size_t i=0;i<dat.Bottom.Size();i++)
    //{
        //Cell * c = dat.Bottom[i];
        //c->RhoBC = dat.rho;
		//if (c->IsSolid) continue;
		//double vy = -1.0 + (c->F[0]+c->F[1]+c->F[3] + 2.0*(c->F[4]+c->F[7]+c->F[8]))/c->RhoBC;
		//c->F[2] = c->F[4] - (2.0/3.0)*c->RhoBC*vy; 
		//c->F[6] = c->F[8] - (1.0/6.0)*c->RhoBC*vy - 0.5*(c->F[3]-c->F[1]);
		//c->F[5] = c->F[7] - (1.0/6.0)*c->RhoBC*vy + 0.5*(c->F[3]-c->F[1]);
        //c->Rho = c->VelDen(c->Vel);
        //dat.Bottom[i]->Initialize(dat.rho,OrthoSys::O);
    //}

    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0;i<dom.Disks.Size();i++)
    {
        dom.Disks[i]->Ff = dom.Disks[i]->M*dat.g;
        double delta;
        delta =   dat.Xmin(0) + dom.Lat[0].dx - dom.Disks[i]->X(0) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(0) += dat.Kn*delta;
        delta = - dat.Xmax(0) - dom.Lat[0].dx + dom.Disks[i]->X(0) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(0) -= dat.Kn*delta;
        delta =   dat.Xmin(1) + dom.Lat[0].dx - dom.Disks[i]->X(1) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(1) += dat.Kn*delta;
        delta = - dat.Xmax(1) - dom.Lat[0].dx + dom.Disks[i]->X(1) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(1) -= dat.Kn*delta;
        if (dat.Closed)
        {
            delta = - dat.Xgate - dom.Lat[0].dx + dom.Disks[i]->X(0) + dom.Disks[i]->R;
            if (delta > 0.0)  dom.Disks[i]->Ff(0) -= dat.Kn*delta;
        }
    }
    
    return ;

    if (dom.Time>0.1*dat.Tf&&dat.Closed)
    {
        #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
        for (size_t i=0;i<dat.Gate.Size();i++)
        {
            Cell * c = dat.Gate[i];
            c->IsSolid = false;

            //iVec3_t ineigh = c->Index + iVec3_t(1,0,0);
            //Cell * nc= dom.Lat[0].GetCell(ineigh);
            //c->Initialize(nc->Rho,nc->Vel);

            //for (size_t j = 1;j<c->Nneigh;j++)
            //{
                //c->Ftemp[j] = c->F[j];
            //}
            //for (size_t j = 1;j<c->Nneigh;j++)
            //{
                //c->F[j]     = c->Ftemp[c->Op[j]];
            //}
           
            //c->F[1] = c->F[3];
            //c->F[5] = c->F[7];
            //c->F[8] = c->F[6];

            c->Rho = c->VelDen(c->Vel);
        }
        dat.Closed = false;
    }
}

void Report(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

}


int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    if (argc==2) Nproc=atoi(argv[1]);
    size_t nx = 1600;
    size_t ny =  800;
    double nu = 0.16;
    double dx = 1.0;
    double dt = 1.0;
    double Tf = 5.0e3;
    double rho= 1000.0;
    size_t Npx= 9;
    size_t Npy =8;
    double R   =15.0;
    Array<double> nuv(2);
    nuv[0] = nu;
    nuv[1] = nu;
    LBM::Domain Dom(D2Q9, nuv, iVec3_t(nx,ny,1), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    //Dom.Lat[0].Rhoref   =  2.0;
    //Dom.Lat[0].G        = -2.0;
    //Dom.Lat[0].Gs       = -0.0;
    Dom.Lat[0].G        = -100.0;
    Dom.Lat[0].Gs       =    0.0;
    Dom.Lat[1].G        =    0.0;
    Dom.Lat[1].Gs       =    0.0;
    Dom.Gmix            = 0.001;
    dat.g               = 0.0,-0.0001,0.0;
    dat.Tf              = Tf;
    dat.Xmin            = 0.0,0.0,0.0;
    dat.Xmax            = nx*dx,ny*dx,0.0;
    dat.Kn              = 1.0e2*rho/500.0;
    dat.rho             = rho;
    dat.Closed          = true;

    //Set solid boundaries
    for (size_t i=0;i<nx;i++)
    {
        //dat.Bottom.Push(Dom.Lat[0].GetCell(iVec3_t(i,0   ,0)));
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
    }
    for (size_t i=0;i<ny;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    }

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        if (j<0.5*ny&&i<0.3*nx)
        {
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize( 500.0 ,v0);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(   0.01,v0);
        }
        //if (j<0.5*ny) Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(1200.0,v0);
        else
        {          
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(   0.01,v0);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(   5.0,v0);
        }
        if (i>=0.3*nx&&i<=0.3*nx+1&&j>0&&j<ny-1) 
        {
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->IsSolid = true;
            dat.Gate.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0)));
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->IsSolid = true;
            dat.Gate.Push(Dom.Lat[1].GetCell(iVec3_t(i,j,0)));
            //std::cout << i << " " << j << std::endl;
        }
        //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(1.0e-6,v0);
    }

    dat.Xgate = 0.3*nx;

    //Dom.AddDisk(0,Vec3_t(0.50*nx,0.6*ny,0.0),Vec3_t(0.0,0.0,0.0),OrthoSys::O,rho,0.1*ny,dt);
    
    double fraction = 0.8; 
    for (size_t i=0;i<Npx;i++)
    for (size_t j=0;j<Npy;j++)
    {
        double ran = 0.9 + (0.1*rand())/RAND_MAX;
        Vec3_t x0((1.5*i+2.0)*2.0*R,(1.5*j+1.0)*2.0*R,0.0);
        if (rand()<fraction*RAND_MAX) 
        {
            Dom.AddDisk(0,x0,Vec3_t(0.0,0.0,0.0),OrthoSys::O,rho,ran*R,dt);
            Dom.Disks[Dom.Disks.Size()-1]->Kn = 1.0*dat.Kn;
            Dom.Disks[Dom.Disks.Size()-1]->Kn = 0.5*dat.Kn;
            Dom.Disks[Dom.Disks.Size()-1]->Gn = 0.016;
        }
    }
    
    //Solving
    Dom.Solve(Tf,Tf/200.0,Setup,Report,"test06",true,Nproc);
}
MECHSYS_CATCH

