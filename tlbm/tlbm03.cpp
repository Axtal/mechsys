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
// Magnus effect

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    Vec3_t                acc;
};

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_force.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << " x" << Util::_8s << " y" << Util::_8s << " z";
        dat.oss_ss                          << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz";
        dat.oss_ss                          << Util::_8s << "wx" << Util::_8s << "wy" << Util::_8s << "wz";
        dat.oss_ss                          << Util::_8s << "fx" << Util::_8s << "fy" << Util::_8s << "fz \n";
    }
    if (!dom.Finished) 
    {
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dom.Particles[0]->x(0) << Util::_8s << dom.Particles[0]->x(1) << Util::_8s << dom.Particles[0]->x(2);
        dat.oss_ss <<                            Util::_8s << dom.Particles[0]->v(0) << Util::_8s << dom.Particles[0]->v(1) << Util::_8s << dom.Particles[0]->v(2);
        dat.oss_ss <<                            Util::_8s << dom.Particles[0]->w(0) << Util::_8s << dom.Particles[0]->w(1) << Util::_8s << dom.Particles[0]->w(2);
        dat.oss_ss <<                            Util::_8s << dom.Particles[0]->F(0) << Util::_8s << dom.Particles[0]->F(1) << Util::_8s << dom.Particles[0]->F(2) << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}


int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    size_t nx = 100;
    size_t ny = 100;
    size_t nz = 100;
    double nu = 0.01;
    double Tf = 10000.0;
    double Dp = 0.001;
    if (argc>=2)
    {
        Nproc = atoi(argv[1]);
        Tf    = atof(argv[2]);
    }
    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), /*dx*/1.0, /*dt*/1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.acc      = Vec3_t(Dp,0.0,0.0);
    Dom.AddSphere(-1,Vec3_t(0.2*nx,0.5*ny,0.5*nz),0.1*nx,3.0);
    Dom.Particles[0]->v = Vec3_t(0.02,0.0,0.0);
    Dom.Particles[0]->w = Vec3_t(0.0,0.0,0.0);
    double rho0 = 0.12;
    Vec3_t v0(0.0,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
    }


    //Solving
    Dom.Solve(Tf,0.01*Tf,NULL,Report,"test03",true,Nproc);
}
MECHSYS_CATCH

