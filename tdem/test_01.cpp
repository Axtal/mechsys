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

/////////////////////// Test 01 the sliding block

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;
using Util::PI;
using Util::FmtErr;
using DEM::Domain;

int main(int argc, char **argv) try
{
    DEM::Domain dom;
    double dt = 1.0e-5;     // time step
    double g = 9.8;         // gravity acceleration
    double velocity = 10.0; // velocity of the sliding block.
    double inipos = -5.0;   // Initial position in the x axis
    double mu = 0.1;        // friction coefficient
    dom.AddCube(/*Tag*/-1,/*initial position*/Vec3_t(inipos,0.0,0.7), /*spheroradius*/0.1,/*side length*/ 1.0,  /*rho*/1.0, /*orientation*/PI/2.0,/*axis of orientation*/&OrthoSys::e0);
    dom.AddPlane(/*Tag*/-2,OrthoSys::O,0.1,100.0,100.0,1.0);
    // Get the position with tag -2 (the plane) and fix it
    DEM::Particle *p = dom.GetParticle(-2);
    p->FixVeloc();
    // Assign gravity and a velocity to particle -1
    p = dom.GetParticle(-1);
    p->v  = velocity , 0.0, 0.0;
    dom.Initialize(dt);
    p->Ff = 0.0,0.0,-p->Props.m*g;
    // setting the tangential viscosity coeffcient Gt to 0 and the friction coefficient to the desired value
    Dict B;
    B.Set(-1,"Gt Kn Kt",0.0,1.0e8,5.0e7);
    B.Set(-2,"Gt Kn Kt",0.0,1.0e8,5.0e7);
    dom.SetProps(B);

    //Assign a different friciton coefficient for different classess of particles
    std::pair<int,int> pr (-1/*tag1*/,-2/*tag2*/);
    dom.FricCoeff[pr] = mu;

    dom.Solve(/*final time*/10.0,/*time step*/dt,/*Output step*/0.1,NULL,NULL,/*file key*/"test_01",/*Render video?*/false);

    //Get the particle final position
    double finpos = p->x(0);

    //Comparing with theoretical formula
    cout << "  Observed L = " << finpos - inipos<< " Theoretical L = " << velocity*velocity/(2*mu*g) << "\n";

    return 0;
}
MECHSYS_CATCH

