/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang
 * Copyright (C) 2020 Siqi Sun                                         *
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
// Elastic beam.


// MechSys
#include <mechsys/lbmmpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    double                         * Vel;
    double                           rho;
    double                            Tf;
    double                            Ly;
    double                            Lz;
    double                            nu;
    double                            CD;
    double                            Re;
    Vec3_t                    fluidforce;
    Array<MPM::Particle *>    EndBeamPar;
    Array<MPM::Particle *>     middlePar;
    Array<double >                    x0;
    Array<double >                    x2; 
    std::ofstream                oss_ss1; // The output of Cd & Re & Dx & Dy changed with Time
    std::ofstream                oss_ss2; // The output of x&y position 
};

void Setup (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double scale = std::min(10.0*dom.Time/dat.Tf,1.0);

    // Cells with prescribed velocity (left boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<dom.LBMDOM.Ndim(1); ++i) //y-axis
    for (size_t j=0; j<dom.LBMDOM.Ndim(2); ++j) //z-axis
    {
        //double c  = 1;
        double * f = dom.LBMDOM.F[0][0][i][j];
        double rho = (f[0]+f[3]+f[4]+f[5]+f[6]+ 2.0*(f[2]+f[8]+f[10]+f[12]+f[14]))/(1.0-scale*dat.Vel[j]/dom.LBMDOM.Cs);
        f[1] = f[2] + (2.0/3.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        f[7] = f[8] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        f[9] = f[10] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        f[11]= f[12] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        f[13]= f[14] + (1.0/12.0)*rho*scale*dat.Vel[j]/dom.LBMDOM.Cs;
        dom.LBMDOM.Vel[0][0][i][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][0][i][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][0][i][j] +=  dom.LBMDOM.F[0][0][i][j][k];
            dom.LBMDOM.Vel[0][0][i][j] +=  dom.LBMDOM.F[0][0][i][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][0][i][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][0][i][j];
    }

    // Cells with prescribed density (right boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<dom.LBMDOM.Ndim(1); ++i) //y-axis
    for (size_t j=0; j<dom.LBMDOM.Ndim(2); ++j) //z-axis
    {
        double * f = dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j];
        double vx = (f[0]+f[3]+f[4]+f[5]+f[6] + 2.0*(f[1]+f[7]+f[9]+f[11]+f[13]))/dat.rho - 1.0;
        f[2] = f[1] - (2.0/3.0)*dat.rho*vx; 
        f[8] = f[7] - (1.0/12.0)*dat.rho*vx;
        f[10]= f[9] - (1.0/12.0)*dat.rho*vx;
        f[12]= f[11] - (1.0/12.0)*dat.rho*vx; 
        f[14]= f[13] - (1.0/12.0)*dat.rho*vx;
        dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j] +=  dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j][k];
            dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] +=  dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j];
    }

    // free slip (top boundary)
    #pragma omp parallel for schedule(static) num_threads(dom.LBMDOM.Nproc)
    for (size_t i=0; i<dom.LBMDOM.Ndim(0); ++i) //x-axis
    for (size_t j=0; j<dom.LBMDOM.Ndim(1); ++j) //y-axis
    {
        double * f = dom.LBMDOM.F[0][i][j][dom.LBMDOM.Ndim(2)-1];
        f[6]  = f[5];
        f[8]  = f[10];
        f[9]  = f[7];
        f[12] = f[14];
        f[13] = f[11];
        dom.LBMDOM.Vel[0][i][j][dom.LBMDOM.Ndim(2)-1] = OrthoSys::O;
        dom.LBMDOM.Rho[0][i][j][dom.LBMDOM.Ndim(2)-1] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][i][j][dom.LBMDOM.Ndim(2)-1] +=  dom.LBMDOM.F[0][i][j][dom.LBMDOM.Ndim(2)-1][k];
            dom.LBMDOM.Vel[0][i][j][dom.LBMDOM.Ndim(2)-1] +=  dom.LBMDOM.F[0][i][j][dom.LBMDOM.Ndim(2)-1][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][i][j][dom.LBMDOM.Ndim(2)-1] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][i][j][dom.LBMDOM.Ndim(2)-1];

    }

}

void Report (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double x     = 0.0;
    double z     = 0.0;
    double xf    = 0.0; // Dx/b
    double zf    = 0.0; // Dz/b
    double x1    = 0.0; // position in x-axis _simulation
    double z1    = 0.0; // position in z-axis _simulation
    double vel2  = 0.0; // velocity of fluid
    double vel   = 0.0; 
    double force = 0.0;
    int    i     = 0;
    int    j     = 0;
    if (dom.idx_out==0)
    {
        String fs1;
        fs1.Printf("tlbmmpm01position.res");
        dat.oss_ss1.open(fs1.CStr());
        dat.oss_ss1 << Util::_10_6 << "X-position" << Util::_8s << "Z-position" << "\n";

        String fs2;
        fs2.Printf("tlbmmpm01time.res");
        dat.oss_ss2.open(fs2.CStr());
        dat.oss_ss2 << Util::_10_6 << "Time" << Util::_8s << "Re" << Util::_8s << "CD" << Util::_8s << "deflection-x"<< Util::_8s << "deflection-z"<< "\n";
    }

    if (dom.idx_out==99)
    {
       //Position in the middle of this elastic beam
        for (size_t ip=0;ip<dat.middlePar.Size();ip++)
        {
            x1 = (dat.middlePar[ip]->x(0)-0.15)/dat.Lz;
            z1 = dat.middlePar[ip]->x(2)/dat.Lz;
            dat.oss_ss1 << Util::_10_6 << x1 << Util::_8s << z1 << "\n";
        } 
    }

    else
    {
        // calculate Re 
        for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
        for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
        for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
        {
            vel += dom.LBMDOM.Vel[0][ix][iy][iz](0);
        }
        vel2 = vel/(dom.LBMDOM.Ndim(0)*dom.LBMDOM.Ndim(1)*dom.LBMDOM.Ndim(2));
        double Re   = vel2*dat.Ly/dat.nu;

        // calculate CD
        for (size_t ip=0; ip < dom.MPMDOM.Corners.Size(); ip++)
        {
            force = force + dom.MPMDOM.Corners[ip]->h(0);
        }
        double CD = force/(0.5*dat.rho*vel2*vel2*dat.Ly*dat.Lz);

        // calculate deflection
        for (size_t ip=0;ip<dat.EndBeamPar.Size();ip++)
        {
            x = x+dat.EndBeamPar[ip]->x(0)-dat.x0[ip];
            z = z+dat.EndBeamPar[ip]->x(2)-dat.x2[ip];
        }
        x /= dat.EndBeamPar.Size();
        z /= dat.EndBeamPar.Size();
        xf = x/dat.Lz;
        zf = z/dat.Lz;
        dat.oss_ss2 << Util::_10_6 << dom.Time << Util::_8s << Re << Util::_8s << CD << Util::_8s << xf << Util::_8s << zf << "\n";
    }           
}

int main(int argc, char **argv) try
{
    //Number of cores
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);

    //Properties of LBM
    size_t nx = 200;
    size_t ny = 200;
    size_t nz = 200;
    double nu = 8.20e-4;
    double dx = 0.0015;

    //Properties of MPM
    size_t ndiv = 40; //number of divisions per x lenght
    double Nu   = 0.3;
    double  E   = 1.23e6;
    double  K   = E/(3.0*(1.0-2.0*Nu));
    double  G   = E/(2.0*(1.0+Nu));
    double rho  = 1030.0;   
    double Lx   = 0.01;
    double Ly   = 0.02;
    double Lz   = 0.1;
    double Cs   = sqrt(E/rho);
    double h    = Lz/ndiv;
    double dt   = 1.0*h/Cs;
    double Dx   = 2.0*Lz/ndiv;
    double Bc   = Dx;

    LBMMPM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;

    dat.Ly = Ly;
    dat.Lz = Lz;
    dat.nu = nu;

    dom.MPMDOM.AddRectangularBeamMesh(-1, Vec3_t(0.5*nx*dx-0.5*Lx,0.5*ny*dx-0.5*Ly,-Bc), Vec3_t(0.5*nx*dx+0.5*Lx,0.5*ny*dx+0.5*Ly,Lz), rho, ndiv);
    dom.MPMDOM.ResizeDomainMesh(Vec3_t(0.5*nx*dx-10.0*Lx,0.5*ny*dx-2.0*Ly,-1.0*Lz), Vec3_t(0.5*nx*dx+10.0*Lx,0.5*ny*dx+2.0*Ly,2.0*Lz),1.5);
    
    //Setting properties for the material points
    for (size_t ip=0; ip < dom.MPMDOM.Particles.Size(); ip++)
    {
        dom.MPMDOM.Particles[ip]->K = K;
        dom.MPMDOM.Particles[ip]->G = G;
        if (dom.MPMDOM.Particles[ip]->x(2)<0.0)
        {
            dom.MPMDOM.Particles[ip]->FixVeloc();
            dom.MPMDOM.Particles[ip]->Tag = -2;
        }
        if (dom.MPMDOM.Particles[ip]->x(2)>Lz*(ndiv-0.5)/ndiv)
        {
            dom.MPMDOM.Particles[ip]->Tag = -3;
            dat.EndBeamPar.Push(dom.MPMDOM.Particles[ip]);
            dat.x0.Push(dom.MPMDOM.Particles[ip]->x(0));
            dat.x2.Push(dom.MPMDOM.Particles[ip]->x(2));
        }
        //if (dom.MPMDOM.Particles[ip]->x(0) > (1/2*nx*dx-Lx/10) && dom.MPMDOM.Particles[ip]->x(0) < (1/2*nx*dx+Lx/10))
        if (dom.MPMDOM.Particles[ip]->x(0) > (0.15-0.001) && dom.MPMDOM.Particles[ip]->x(0) < (0.15+0.001))
        {
            dom.MPMDOM.Particles[ip]->Tag = -4;
            dat.middlePar.Push(dom.MPMDOM.Particles[ip]);
        }
    }

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
        if (xn(2)<0.0)
        {
            dom.MPMDOM.Nodes[in].FixVeloc();
        }
    }

    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,1220.0/*rho*/,v);
        if (iz==0)
        {
            dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
        }
    }  

    dat.Vel = new double[nz];
    double Re = 9.0;
    double umax  = 1.5*Re*nu/Ly;
    for (size_t i=0;i<nz;i++)
    {
        dat.Vel[i] = (i-(nz-1))*(i-(nz-1))*(-umax/((nz-1.5)*(nz-1.5)))+umax;
    }
    dat.rho = 1220.0;
    
    double Tf = 1.5; // Final Time
    dat.Tf = Tf;
    dat.fluidforce = Vec3_t(5.0e-3*dx/dt/Tf,0.0,0.0);
    
    dom.MPMDOM.Gn = 1.0e0;
    dom.Solve(Tf,Tf/100,Setup,Report,"tlbmmpm01",true,Nproc);
}
MECHSYS_CATCH

