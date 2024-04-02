/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang                                         *
 * Copyright (C) 2020 Siqi Sun                                          *
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
// Fish tail loaded from mesh


// MechSys
#include <mechsys/mpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    Array<MPM::Particle *> ForcePar; //Forced particles
    Array<MPM::Corner   *> ForceCor; //Forced particles
    size_t                     ndiv; //number of divisions per z lenght
    double                       Vb; //Maximun applied body velocity 
    double                       Nc; //Number of cycles for the vibration
    double                       Tf; //Final simulation time
    double                       Bc; //Boundary between fixed points and moving points
    double                       Lz; //Length of skeleton beam
    double                     xmin; //Minimum x-coordinate for skeleton beam  
    double                     xmax; //Maximum x-coordinate for skeleton beam
    double                     ymin; //Minimum y-coordinate for skeleton beam
    double                     ymax; //Maximum y-coordinate for skeleton beam
    double                     zmin; //Minimum z-coordinate for skeleton beam
};

void Setup (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    double w = 2.0*M_PI*dat.Nc/dat.Tf; 
    double vt = dat.Vb/dat.Lz*cos(w*dom.Time);
    Vec3_t Om(vt,0.0,0.0); //Velocity
    Vec3_t xf(dat.xmin+(dat.xmax-dat.xmin)/2,dat.ymin+(dat.ymax-dat.ymin)/2,dat.Bc); //Coordinate of the center point
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)   
    for (size_t ip=0; ip<dat.ForcePar.Size(); ip++)
    {
        Vec3_t xb = dat.ForcePar[ip]->x-xf;
        dat.ForcePar[ip]->vf = cross(Om,xb);
    } 
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)   
    for (size_t ic=0; ic<dat.ForceCor.Size(); ic++)
    {
        Vec3_t xb = *dat.ForceCor[ic]->x-xf;
        dat.ForceCor[ic]->vf = cross(Om,xb);
    }
}

int main(int argc, char **argv) try
{
    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
   
    size_t Nproc = 0.75*omp_get_max_threads();
    //size_t Nproc = 1;
    if (argc>1) Nproc = atoi(argv[1]);
    size_t ndiv = 20; //number of divisions per z lenght 
    //size_t ndiv = 5; //number of divisions per z lenght 
    double K    = 10.0e6;
    double Nu   = 0.48;
    double  E   = (3.0*(1.0-2.0*Nu))*K;
    double  G   = E/(2.0*(1.0+Nu));
    double rho  = 3000.0;
    double scale= 1.0e0;
    double Cs   = sqrt(scale*K/rho);
    double h    = 1.0/ndiv;
    double dt   = 1.0*h/Cs;
    double Lx   = 0.6;
    double Ly   = 0.6;
    dat.ndiv    = ndiv;

    //dom.AddFromJson(-1,"fish.msh",rho,1.0,ndiv); //Adding the fishtile
    dom.AddFromJsonMesh(-1,"fish.msh",OrthoSys::O,OrthoSys::O,OrthoSys::e0,0.0,rho,1.0,ndiv); //Adding the fishtile
    Vec3_t xmin,xmax;
    dom.BoundingBox(xmin,xmax);
    double Dz   = 2.0*(xmax(2)-xmin(2))/ndiv;
    double Bc   = Dz+xmin(2);
    double Lz   = (xmax(2)-xmin(2))*1.0/3.0;
    dat.Bc      = Bc;
    dat.Lz      = Lz;
    dat.xmin    = xmin(0);
    dat.xmax    = xmax(0);
    dat.ymin    = xmin(1);
    dat.ymax    = xmax(1);
    dat.zmin    = xmin(2);
    dom.ResizeDomainMesh(xmin-1.0*(xmax(2)-xmin(2)),xmax+1.0*(xmax(2)-xmin(2)), 0.5); //Mesh size
    
    //Setting properties for the material points
    for (size_t ip=0; ip < dom.Particles.Size(); ip++)
    {
        dom.Particles[ip]->K = K;
        dom.Particles[ip]->G = G;

        if (dom.Particles[ip]->x(2)<Bc)
        {
            dom.Particles[ip]->FixVeloc();
            dom.Particles[ip]->Tag = -2;
        }

        if (dom.Particles[ip]->x(0)>(xmin(0)+(xmax(0)-xmin(0))/2-Lx/2) && dom.Particles[ip]->x(0)<(xmin(0)+(xmax(0)-xmin(0))/2+Lx/2) && dom.Particles[ip]->x(1)>(xmin(1)+(xmax(1)-xmin(1))/2-Ly/2) && dom.Particles[ip]->x(1)<(xmin(1)+(xmax(1)-xmin(1))/2+Ly/2) && dom.Particles[ip]->x(2)>Bc && dom.Particles[ip]->x(2)<(xmin(2)+Lz))
        {
            dom.Particles[ip]->Tag = -3;
		    dom.Particles[ip]->FixVeloc();
            dat.ForcePar.Push(dom.Particles[ip]);
        }
    }

    //Fixing positions and forces over corner points 
    for (size_t ic=0; ic < dom.Corners.Size(); ic++)
    {
        if ((*dom.Corners[ic]->x)(2)<Bc)
        {
            dom.Corners[ic]->FixVeloc();
        }
        if ((*dom.Corners[ic]->x)(0)>(xmin(0)+(xmax(0)-xmin(0))/2-Lx/2) && (*dom.Corners[ic]->x)(0)<(xmin(0)+(xmax(0)-xmin(0))/2+Lx/2) && (*dom.Corners[ic]->x)(1)>(xmin(1)+(xmax(1)-xmin(1))/2-Ly/2) && (*dom.Corners[ic]->x)(1)<(xmin(1)+(xmax(1)-xmin(1))/2+Ly/2) && (*dom.Corners[ic]->x)(2)>Bc && (*dom.Corners[ic]->x)(2)<(xmin(2)+Lz))
        {
		    dom.Corners[ic]->FixVeloc();
            dat.ForceCor.Push(dom.Corners[ic]);
        }
    }

    //Setting properties for nodes
    for (size_t in=0; in < dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        if (xn(2)<Bc)
        {
            dom.Nodes[in].FixVeloc();
        }
    }

    //dom.ParticleToNode();
    dom.Gn = 0.0;
    double Tf = 100.0; // Final Time
    dat.Tf = Tf;  
    dat.Vb = 1.0e-2;  //velocity applied at particles   
    dat.Nc = 5.0; // number of cycles for vibrations
  

    dom.Solve(Tf,dt,Tf/100,Setup,NULL,"tmpm03",true,Nproc);

}
MECHSYS_CATCH

