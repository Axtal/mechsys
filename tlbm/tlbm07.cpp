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
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c   = dom.Lat[0].Cells[i];
        c->BForcef = c->Rho*dat.g;
        c          = dom.Lat[1].Cells[i];
        c->BForcef = c->Rho*dat.g;
    }
    //for (size_t i=0;i<dom.Particles.Size();i++)
    //{
        //dom.Particles[i]->Ff = dom.Particles[i]->Props.m*dat.g;
        //double delta;
        //delta =   dat.Xmin(0) - dom.Particles[i]->x(0) + dom.Particles[i]->Props.R;
        //if (delta > 0.0)  dom.Particles[i]->Ff(0) += dat.Kn*delta;
        //delta = - dat.Xmax(0) + dom.Particles[i]->x(0) + dom.Particles[i]->Props.R;
        //if (delta > 0.0)  dom.Particles[i]->Ff(0) -= dat.Kn*delta;
        //delta =   dat.Xmin(1) - dom.Particles[i]->x(1) + dom.Particles[i]->Props.R;
        //if (delta > 0.0)  dom.Particles[i]->Ff(1) += dat.Kn*delta;
        //delta = - dat.Xmax(1) + dom.Particles[i]->x(1) + dom.Particles[i]->Props.R;
        //if (delta > 0.0)  dom.Particles[i]->Ff(1) -= dat.Kn*delta;
        //delta =   dat.Xmin(2) - dom.Particles[i]->x(2) + dom.Particles[i]->Props.R;
        //if (delta > 0.0)  dom.Particles[i]->Ff(2) += dat.Kn*delta;
        //delta = - dat.Xmax(2) + dom.Particles[i]->x(2) + dom.Particles[i]->Props.R;
        //if (delta > 0.0)  dom.Particles[i]->Ff(2) -= dat.Kn*delta;
    //}
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    if (argc==2) Nproc=atoi(argv[1]);
    size_t nx = 100;
    size_t ny = 100;
    size_t nz = 100;
    double nu = 0.01;
    double dx = 1.0;
    double dt = 1.0;
    double rho= 3000.0;
    
    Array<double> Nu(2);
    Nu[0] = nu;
    Nu[1] = nu;

    LBM::Domain Dom(D3Q15, Nu, iVec3_t(nx,ny,nz), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Step = 1;
    Dom.Lat[0].G    = -200.0;
    Dom.Lat[0].Gs   = -0.0;
    Dom.Lat[1].G    = 0.0;
    Dom.Lat[1].Gs   = 0.0;
    Dom.Gmix        =  0.001;
    dat.g           = 0.0,-0.001,0.0;
    dat.Xmin        = 0.0,0.0,0.0;
    dat.Xmax        = nx*dx,ny*dx,nz*dx;
    dat.Kn          = 1.0e4*rho/500.0;

    //Set solid boundaries
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<nz;j++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,j))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,j))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,j))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,ny-1,j))->IsSolid = true;
    }
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,j,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,j,nz-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0   ,j,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,j,i))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,j,0   ))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,j,nz-1))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(0   ,j,i))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(nx-1,j,i))->IsSolid = true;
    }

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        //if (j<0.6*ny)
        if (j<ny/2.0)
        {
            Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(2300.0,v0);
            //Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(3000.0,v0);
            Dom.Lat[1].GetCell(iVec3_t(i,j,k))->Initialize(0.01  ,v0);
        }
        else
        {
            Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(0.01 ,v0);
            Dom.Lat[1].GetCell(iVec3_t(i,j,k))->Initialize(100.0,v0);
        }
    }

    Dom.GenBox(-6,1.02*nx,1.02*ny,1.02*nz,1.0,1.1,rho);
    Dom.Center(Vec3_t(0.5*nx*dx,0.5*ny*dx,0.5*nz*dx));
    Dom.GetParticle( -6)->FixVeloc();
    Dom.GetParticle( -7)->FixVeloc();
    Dom.GetParticle( -8)->FixVeloc();
    Dom.GetParticle( -9)->FixVeloc();
    Dom.GetParticle(-10)->FixVeloc();
    Dom.GetParticle(-11)->FixVeloc();
    Dom.GetParticle( -6)->Bdry = true;
    Dom.GetParticle( -7)->Bdry = true;
    Dom.GetParticle( -8)->Bdry = true;
    Dom.GetParticle( -9)->Bdry = true;
    Dom.GetParticle(-10)->Bdry = true;
    Dom.GetParticle(-11)->Bdry = true;

    Array<int> deltag;
    deltag.Push(-8);
    Dom.DelParticles(deltag);

    Dom.AddSphere(-1,Vec3_t(0.58*nx*dx,0.65*ny*dx,0.58*nz*dx),0.1*ny,0.4*rho);
    Dom.AddSphere(-2,Vec3_t(0.57*nx*dx,0.85*ny*dx,0.57*nz*dx),0.1*ny,rho    );
    Dom.AddSphere(-3,Vec3_t(0.43*nx*dx,0.65*ny*dx,0.42*nz*dx),0.1*ny,rho    );
    Dom.AddSphere(-4,Vec3_t(0.36*nx*dx,0.85*ny*dx,0.63*nz*dx),0.1*ny,0.3*rho);
    Dom.AddSphere(-5,Vec3_t(0.70*nx*dx,0.65*ny*dx,0.40*nz*dx),0.1*ny,0.6*rho);

    //Dom.AddRecBox(-1,OrthoSys::O,Vec3_t(0.2*nx*dx,0.2*ny*dx,0.1*nz*dx),1.0,0.5*rho);
    //Dom.GetParticle(-1)->Position(Vec3_t(0.5*nx*dx,0.8*ny*dx,0.5*nz*dx));
//
    //Dom.AddFromJson(-2, "test.msh", 1.0, 0.4*rho, 1.5);
    //Dom.GetParticle(-2)->Position(Vec3_t(0.35*nx*dx,0.65*ny*dx,0.35*nz*dx));


    for (size_t i=0;i<Dom.Particles.Size();i++)
    {
        Dom.Particles[i]->Ff = Dom.Particles[i]->Props.m*dat.g;
        Dom.Particles[i]->Props.Kn = 1.0*dat.Kn;
        Dom.Particles[i]->Props.Kt = 0.5*dat.Kn;
        Dom.Particles[i]->Props.Gn = 0.16;
    }

    //Dom.WriteXDMF("test07");
    //Solving
    Dom.Solve(4000.0,10.0,Setup,NULL,"test07",true,Nproc);
}
MECHSYS_CATCH

