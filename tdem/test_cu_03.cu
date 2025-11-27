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

/////////////////////// Test 03 column collapse for testing efficiency

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

int main(int argc, char **argv) try
{
    DEM::Domain dom;
    size_t Nproc  = 1;       // Number of cores
    size_t Nthread= 256;    // Number of CUDA threads
    double dt = 1.0e-5;     // time step
    double g = 9.8;         // gravity acceleration
    double mu = 0.4;        // friction coefficient
    double theta = 0.0*M_PI/180.0;
    if (argc>=2) Nproc      =atoi(argv[1]);
    if (argc>=3) dom.Alpha  =atof(argv[2]);
    if (argc>=4) dom.Nthread=atoi(argv[3]);
    //dom.Alpha = 100.0;
    //dom.Beta  = 4.0;
    dom.AddVoroPack (-1, 0.1, 20.0,20.0,240.0, 20,20,240, 1.0, false, true, 1000, 0.8);
    Vec3_t Xmin(-0.5*10.0,-0.5*10.0,-0.5*120.0);
    Vec3_t Xmax = -Xmin;
    //dom.GenSpheresBox (-1, Xmin, Xmax, 0.5, 1.0, "HCP", 1000, 0.8, 0.8);
    //Vec3_t Xmin,Xmax;
    dom.BoundingBoxTag(Xmin,Xmax,-1);
    dom.AddPlane(/*Tag*/-2,Vec3_t(0.0,0.0,Xmin(2)-0.1),0.1,400.0,400.0,1.0);
    DEM::Particle *p = dom.GetParticle(-2);
    p->FixVeloc();
    for (size_t ip=0;ip<dom.Particles.Size();ip++)
    {
        dom.Particles[ip]->Ff = Vec3_t(0.0,dom.Particles[ip]->Props.m*9.8*sin(theta),-dom.Particles[ip]->Props.m*9.8*cos(theta));
    }
    //dom.Xmax =  6.0;
    //dom.Xmin = -6.0;    
    //dom.Ymax =  6.0;
    //dom.Ymin = -6.0;    

    Dict B;
    B.Set(-1,"Gn Gt Kn Kt Mu",-0.2,0.0,1.0e7,5.0e6,mu);
    B.Set(-2,"Gn Gt Kn Kt Mu",-0.2,0.0,1.0e7,5.0e6,mu);
    dom.SetProps(B);

    double Tf = 10.0;
    dom.Solve(/*final time*/Tf,/*time step*/dt,/*Output step*/Tf/100,NULL,NULL,/*file key*/"test_cu_03",/*Render video?*/2,Nproc);

    return 0;
}
MECHSYS_CATCH

