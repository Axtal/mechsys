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
    double Kn;
    Vec3_t g;
    Vec3_t Xmin;
    Vec3_t Xmax;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c = dom.Lat[0].Cells[i];
        c->BForcef = c->Density()*dat.g;
    }
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0;i<dom.Disks.Size();i++)
    {
        dom.Disks[i]->Ff = dom.Disks[i]->M*dat.g;
        double delta;
        delta =   dat.Xmin(0) - dom.Disks[i]->X(0) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(0) += dat.Kn*delta;
        delta = - dat.Xmax(0) + dom.Disks[i]->X(0) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(0) -= dat.Kn*delta;
        delta =   dat.Xmin(1) - dom.Disks[i]->X(1) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(1) += dat.Kn*delta;
        delta = - dat.Xmax(1) + dom.Disks[i]->X(1) + dom.Disks[i]->R;
        if (delta > 0.0)  dom.Disks[i]->Ff(1) -= dat.Kn*delta;
    }
}


int main(int argc, char **argv) try
{
    size_t Nproc  = 1; 
    size_t nx = 200;
    size_t ny = 200;
    double nu = 0.1;
    double dx = 1.0;
    double dt = 1.0;
    double rho= 100.0;
    if (argc>=2) Nproc = atoi(argv[1]);
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Lat[0].G    = -200.0;
    Dom.Lat[0].Gs   = -200.0;
    Dom.Sc          = 0.0;
    dat.g        = 0.0,-0.001,0.0;
    dat.Xmin     = 0.0,0.0,0.0;
    dat.Xmax     = nx*dx,ny*dx,0.0;
    dat.Kn       = 1.0e5*rho/500.0;

    //Set solid boundaries
    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }
    for (size_t i=0;i<ny;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    }

    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        if (j<ny/2.0) Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(1300.0,v0);
        else          Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(  50.0,v0);
    }

    Dom.AddDisk(0,Vec3_t(0.30*nx,0.6*ny,0.0),Vec3_t(0.0,0.0,0.0),OrthoSys::O,rho,0.1*ny,dt);
    Dom.AddDisk(0,Vec3_t(0.40*nx,0.8*ny,0.0),Vec3_t(0.0,0.0,0.0),OrthoSys::O,rho,0.1*ny,dt);
    Dom.AddDisk(0,Vec3_t(0.50*nx,0.6*ny,0.0),Vec3_t(0.0,0.0,0.0),OrthoSys::O,rho,0.1*ny,dt);
    Dom.AddDisk(0,Vec3_t(0.60*nx,0.8*ny,0.0),Vec3_t(0.0,0.0,0.0),OrthoSys::O,rho,0.1*ny,dt);
    Dom.AddDisk(0,Vec3_t(0.70*nx,0.6*ny,0.0),Vec3_t(0.0,0.0,0.0),OrthoSys::O,rho,0.1*ny,dt);
    //Dom.Disks[Dom.Disks.Size()-1]->FixVelocity();

    //Solving
    Dom.Solve(10000.0,50.0,Setup,NULL,"test05","true",Nproc);
    //Dom.Solve(10000.0,20.0,NULL,NULL,"test05");
}
MECHSYS_CATCH

