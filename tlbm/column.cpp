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

// Hydraulic column test for force scheme validation

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    std::ofstream      oss_ss;           ///< file for particle data
    Vec3_t G;                            ///< gravity
    double Cs2;                          ///< Speed of sound squared
    double L;                            ///< Domain height in length
    double rhof;                         ///< Fluid density

};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c   = dom.Lat[0].Cells[i];
        c->BForcef = c->Rho*dat.G;
    }
}

void Report (LBM::Domain & dom, void * UD)
{
    //A file will be writtn at several time steps with the density across the height of the column compared with the theoretical value
    UserData & dat = (*static_cast<UserData *>(UD));
    String fs;
    fs.Printf("%s_%04d_density.res",dom.FileKey.CStr(),dom.idx_out);
    dat.oss_ss.open(fs.CStr());
    dat.oss_ss << Util::_8s << "Density" << Util::_8s <<  "Theory" << Util::_8s << "Force" << Util::_8s << "Force Formula \n";
    double g    = fabs(dat.G(1));
    double dx   = dom.Lat[0].dx;
    double rho0 = dat.rhof*g*dat.L/(dat.Cs2*(exp(dat.L*g/dat.Cs2)-1)); //Value of the density at the top from the analytical model
    for (size_t i = 1; i<dom.Lat[0].Ndim(1); i++)
    {
        double theo = rho0*exp((dat.L- dx*(i+1))*g/dat.Cs2); //Theretical value at a given height
        dat.oss_ss << Util::_8s << dom.Lat[0].GetCell(iVec3_t(0,i,0))->Rho << Util::_8s << theo << Util::_8s << dom.Disks[0]->F(1) << Util::_8s << M_PI*dom.Disks[0]->R*dom.Disks[0]->R*dat.rhof*g << std::endl;
    }
    dat.oss_ss.close();
    std::cout << dom.Disks[0]->X(1) << " " << dom.Disks[0]->V(1) << " " << dom.Disks[0]->F(1) << std::endl;
}

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc>=2) nproc=atoi(argv[1]);
    size_t nx = 200;
    size_t ny = 200;
    double dx = 1.0;
    double dt = 1.0;
    double Tf = 800000.0;
    double nu = 0.16e0;

    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    dat.G   = 0.0,-0.00001,0.0;
    dat.Cs2 = dx*dx/(3.0*dt*dt);
    dat.L   = dx*ny;
    dat.rhof= 1.0;


    //Initializing values
    double g = fabs(dat.G(1));
    double rho0 = dat.rhof*g*dat.L/(dat.Cs2*(exp(dat.L*g/dat.Cs2)-1)); //Value of the density at the top from the analytical model
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        size_t ix = Dom.Lat[0].Cells[i]->Index(0);
        size_t iy = Dom.Lat[0].Cells[i]->Index(1);
        double rho = rho0*exp((dat.L- dx*(iy+1))*g/dat.Cs2);
        if (sqrt((ix-0.5*nx)*(ix-0.5*nx)+(iy-0.5*ny)*(iy-0.5*ny))<nx/20.0)
        {
            Dom.Lat[0].Cells[i]->Initialize(dat.rhof, OrthoSys::O);
        }
        else
        {
            Dom.Lat[0].Cells[i]->Initialize(rho, OrthoSys::O);
        }
    }

    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0,0))->IsSolid = true;
    }

    Dom.AddDisk(-1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.0),OrthoSys::O,OrthoSys::O,2.0,nx/20.0,dt);
    //Dom.Disks[0]->FixVeloc();

    //Solving
    Dom.Solve(Tf,Tf/200.0,Setup,Report,"column",true,nproc);
 
}
MECHSYS_CATCH


