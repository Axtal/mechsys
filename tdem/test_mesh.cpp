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
    //String filename = "test.msh";
    String filename = "test.obj";
    if (argc>1) filename = argv[1];

    // test the mesh particle generation
    double R = 0.1;
    Domain d;
    //d.AddFromJson(-1, filename.CStr(), R, 1.0, 1.0);
    d.AddFromOBJ(-1, filename.CStr(), R, 1.0, 1.0, true);
    d.Particles[0]->Position(OrthoSys::O);
    double epsilon = R/pow(d.Particles[0]->Props.V,1.0/3.0);


    std::cout << "[1;33m\n--- Results for mass properties from surface integrals --------------------------------------[0m\n";
    cout << "  Volume            = " << d.Particles[0]->Props.V << "\n";
    cout << "  Center of mass    = " << d.Particles[0]->x(0)    << ", " << d.Particles[0]->x(1) << ", " << d.Particles[0]->x(2) << "\n";
    cout << "  Moment of inertia = " << d.Particles[0]->I(0)    << ", " << d.Particles[0]->I(1) << ", " << d.Particles[0]->I(2) << "\n";
    cout << "  Quaternion        = " << d.Particles[0]->Q(0)    << ", " << d.Particles[0]->Q(1) << ", " << d.Particles[0]->Q(2) << ", " << d.Particles[0]->Q(3) << "\n";
    cout << endl;



    d.Particles[0]->Props.V *= 1.0 + 3.0*epsilon;
    d.Particles[0]->I       *= 1.0 + 5.0*epsilon;

    std::cout << "[1;33m\n--- Results for mass properties from surface integrals corrected by spheroradius ------------[0m\n";
    cout << "  Volume            = " << d.Particles[0]->Props.V << "\n";
    cout << "  Center of mass    = " << d.Particles[0]->x(0)    << ", " << d.Particles[0]->x(1) << ", " << d.Particles[0]->x(2) << "\n";
    cout << "  Moment of inertia = " << d.Particles[0]->I(0)    << ", " << d.Particles[0]->I(1) << ", " << d.Particles[0]->I(2) << "\n";
    cout << "  Quaternion        = " << d.Particles[0]->Q(0)    << ", " << d.Particles[0]->Q(1) << ", " << d.Particles[0]->Q(2) << ", " << d.Particles[0]->Q(3) << "\n";
    cout << endl;

    d.Particles[0]->CalcProps(30000);
    // output
    std::cout << "[1;33m\n--- Results for mass properties from Monte Carlo integration -------------------------------[0m\n";
    cout << "  Volume            = " << d.Particles[0]->Props.V << "\n";
    cout << "  Center of mass    = " << d.Particles[0]->x(0)    << ", " << d.Particles[0]->x(1) << ", " << d.Particles[0]->x(2) << "\n";
    cout << "  Moment of inertia = " << d.Particles[0]->I(0)    << ", " << d.Particles[0]->I(1) << ", " << d.Particles[0]->I(2) << "\n";
    cout << "  Quaternion        = " << d.Particles[0]->Q(0)    << ", " << d.Particles[0]->Q(1) << ", " << d.Particles[0]->Q(2) << ", " << d.Particles[0]->Q(3) << "\n";
    cout << endl;

    //d.AddFromJson(-2, "test.msh", 0.01, 1.0, 1.0, true);
    //Vec3_t t(20.0,0.0,0.0);
    //d.Particles[1]->Translate(t);

    //Quaternion_t q;
    //NormalizeRotation(M_PI/3.0,OrthoSys::e2,q);
    //Mat3_t Rot;
    //RotationMatrix(q,Rot);
    //std::cout << Rot << std::endl;
    
    // draw
    d.WriteXDMF ("test_mesh");
    d.Save      ("test_mesh");
    return 0;
}
MECHSYS_CATCH
