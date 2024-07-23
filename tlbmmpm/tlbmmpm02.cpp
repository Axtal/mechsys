/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang
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
// Settled cube.


// MechSys
#include <mechsys/lbmmpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    double                         * Vel;
    double                           rho;
    double                          rhof;
    double                            Tf;
    double                            Ly;
    double                            Lz;
    double                            nu;
    double                            CD;
    double                            Re;
    Vec3_t                    fluidforce;
    Array<MPM::Particle *>    EndBeamPar;
    Array<double >                    x0;
    Array<double >                    x2; 
    std::ofstream                oss_ss1;

};

void Setup (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
    for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
    for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
    {
        dom.LBMDOM.BForce[0][ix][iy][iz] = dom.LBMDOM.Rho[0][ix][iy][iz]*dat.fluidforce;
    }
}

void Report (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs1;
        fs1.Printf("tlbmmpm02.res");
        dat.oss_ss1.open(fs1.CStr());
        dat.oss_ss1 << Util::_10_6 << "Time" << Util::_8s << "Re" << Util::_8s << "CD" << Util::_8s << "CD(formula)" << "\n";
    }
    else 
    {
        double x     = 0.0;
        double z     = 0.0;
        double xf    = 0.0; // Dx/b
        double zf    = 0.0; // Dz/b
        double vel   = 0.0; // velocity of fluid
        double force = 0.0;
        int    i     = 0;
        int    j     = 0;

        // calculate Re 
        for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
        for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
        for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
        {
            vel += dom.LBMDOM.Vel[0][ix][iy][iz](0);
            i = i+1;
        }
        vel /= i;
        double Re   = 1.2407*vel*dat.Ly/dat.nu;
        double CD1  = 128*pow(Re,-0.8)/5.5;

        // calculate CD
        for (size_t ip=0; ip < dom.MPMDOM.Corners.Size(); ip++)
        {
            force = force + dom.MPMDOM.Corners[ip]->h(0);
        }
        double CD = force/(0.5*dat.rhof*vel*vel*dat.Ly*dat.Lz*1.2090);
        dat.oss_ss1 << Util::_10_6 << dom.Time << Util::_8s << Re << Util::_8s << CD << Util::_8s << CD1 << std::endl;
    }
}

int main(int argc, char **argv) try
{
    //Number of cores
    //if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    size_t Nproc = 1; 
    if (argc>=2) Nproc=atoi(argv[1]);
    //String filekey  (argv[1]);
    //String filename (filekey+".inp");
    //if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    //ifstream infile(filename.CStr());
    //double bforce;
    //{
        //infile >> bforce;           infile.ignore(200,'\n');
    //}
    double bforce = 2.5e-1;
    size_t nx = 241;
    size_t ny = 61;
    size_t nz = 61;
    double dx = 0.1;
    double rhof = 1000.0; // density of fluid

    //Properties of MPM
    size_t ndiv = 1; //number of divisions per x lenght
    double K    = 10.0e4; //Bulk modulus
    double Nu   = 0.3; //Poisson ratio
    double  E   = (3.0*(1.0-2.0*Nu))*K; //Young modulus
    double  G   = E/(2.0*(1.0+Nu)); //Shear modulus
    double rho  = 694.95; //density of solid
    double Lx   = 14.5*dx; //length of the beam in x   
    double Ly   = 14.5*dx;
    double Lz   = 14.5*dx;
    double Cs   = sqrt(E/rho); //Speed of sound
    double h    = Ly/ndiv; 
    double h2   = Lz/ndiv; //length per mesh
    double dt   = 0.2;
    double Dx   = 2.0*Lz/ndiv;
    double Bc   = Dx;
    double nu   = dx*dx/(dt*24.0); //viscosity

    LBMMPM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;

    dat.Ly = Ly;
    dat.Lz = Lz;
    dat.nu = nu;
    dat.rhof = rhof;
    dom.MPMDOM.AddRectangularBeamMesh(-1, Vec3_t(0.5*nx*dx-0.5*Lx,0.5*ny*dx-0.5*Ly,0.5*nz*dx-0.5*Lz), Vec3_t(0.5*nx*dx+0.5*Lx,0.5*ny*dx+0.5*Ly,0.5*nz*dx+0.5*Lz), rho, ndiv);
    dom.MPMDOM.ResizeDomain(Vec3_t(0.5*nx*dx-5.0*Lx,0.5*ny*dx-5.0*Ly,0.5*nz*dx-5.0*Lz), Vec3_t(0.5*nx*dx+5.0*Lx,0.5*ny*dx+5.0*Ly,0.5*nz*dx+5.0*Lz),Dx);

    //Setting properties for the material points
    double v0 = 1.0;
    for (size_t ip=0; ip < dom.MPMDOM.Particles.Size(); ip++)
    {
        dom.MPMDOM.Particles[ip]->K = K;
        dom.MPMDOM.Particles[ip]->G = G;
        dom.MPMDOM.Particles[ip]->FixVeloc();
        dom.MPMDOM.Particles[ip]->Tag = -2;
    }
    dom.MPMDOM.Gn = 1.0e0;

    //Fixing positions and forces over corner points 
    for (size_t ic=0; ic < dom.MPMDOM.Corners.Size(); ic++)
    {
        if ((*dom.MPMDOM.Corners[ic]->x)(2)<0.0)
        {
            dom.MPMDOM.Corners[ic]->FixVeloc();
        }
    }

    //Setting properties for nodes
    for (size_t in=0; in < dom.MPMDOM.Nnodes; in++)
    {
        Vec3_t xn;
        dom.MPMDOM.NodePosition(in,xn);
        dom.MPMDOM.Nodes[in].FixVeloc();
    }

    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,1000.0/*rho*/,v);
    }  

    dat.Vel = new double[ny];
    double umax    = 1.0e-3*dx/dt;
    dat.rho = 1000.0;
    double Tf = 1.0e5*dt; // Final Time
    dat.Tf = Tf;
    dat.fluidforce = Vec3_t(bforce*dx/dt/Tf,0.0,0.0);
    dom.Solve(Tf,Tf/100,Setup,Report,"tlbmmpm02",true,Nproc);
}
MECHSYS_CATCH
