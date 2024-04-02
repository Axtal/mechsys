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

/////////////////////// Test 02 the turning table

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;
using Util::PI;
using Util::FmtErr;
using DEM::Domain;
struct UserData
{
    DEM::Particle *p;   // The sphere
    double omega;   // The angular velocity of the turntable
    double x0;      // The initial position in x
    double vx;      // velocity in the y axis
    std::ofstream      oss_ss; // The output file
};

void Report (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    double wr = (2.0/7.0)*dat.omega;  // frequency of the orbit
    // output triaxial test data
    // header
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_trayectory.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "xt" << Util::_8s << "yt" << "\n";

    }
    else 
    {
        double xt = dat.x0 + dat.vx/wr*sin(wr*dom.Time);
        double yt = dat.vx/wr -(dat.vx/wr)*cos(wr*dom.Time);
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dat.p->x(0) << Util::_8s << dat.p->x(1) << Util::_8s << xt << Util::_8s << yt << "\n";
    }
}

int main(int argc, char **argv) try
{
    UserData dat;
    DEM::Domain dom(&dat);
    double dt = 1.0e-5;      // time step
    double g = 980.0;        // gravity acceleration
    double angveltab = 10.0; // angular velocity of the sliding block
    double inipos = -4.0;    // Initial position in the x axis
    double inivel = 1.0;     // Initial velocity of the particle
    double rho = 7.8;        // Density of the ball
    double R = 0.635;        // Radius of the ball
    dat.omega = angveltab;
    dat.x0 = inipos;
    dat.vx = -inivel;
    // Adding the turning table
    dom.AddPlane(/*Tag*/-1,OrthoSys::O,0.1,20.0,20.0,1.0);

    // Get the position with tag -1 (the plane) and fix it
    DEM::Particle *p = dom.GetParticle(-1);
    p->FixVeloc();
    p->w = 0.0,0.0,angveltab;

    // Add the sphere above the plane with an spin that garantees no initial sliding and also gravity
    dom.AddSphere(/*Tag*/-2,/*initial pos*/ Vec3_t(inipos,0.0,0.1+R),/*radius*/R,/*rho*/rho);
    p = dom.GetParticle(-2);
    p->v = -inivel,0.0,0.0;
    p->w = angveltab*inipos/R,-inivel/R , 0.0;
    dom.Initialize(dt);
    p->Ff = 0.0,0.0,-p->Props.m*g;
    dat.p = p;

    //Set Gt = 0 the friction coefficient and the rolling coefficient equal to zero
    Dict B;
    B.Set(-1,"Gt Mu Eta Beta Kn Kt",0.0,0.15,0.0,0.,1.0e8,1.0e8);
    B.Set(-2,"Gt Mu Eta Beta Kn Kt",0.0,0.15,0.0,0.,1.0e8,1.0e8);
    dom.SetProps(B);



    dom.Solve(/*final time*/30.0,/*time step*/dt,/*Output step*/0.1,NULL,&Report,/*file key*/"test_02",/*Render video?*/2);

    return 0;
}
MECHSYS_CATCH

