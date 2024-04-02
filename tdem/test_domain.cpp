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
    // test the varius Domain and Particle constructors.
    Domain d;
    srand(1000);

    // add cube
    d.AddCube (-1, Vec3_t(20,0,0),1.0,30.0,3.0);
    Array<Vec3_t> Vresult;
    DEM::Particle * Pa = d.Particles[0];
    Array<Vec3_t> Vtemp(Pa->Verts.Size());
    Array<Vec3_t> Vres (Pa->Verts.Size());
    for (size_t j=0;j<Pa->Verts.Size();j++)
    {
        Vtemp[j] = *Pa->Verts[j];
        Vres [j] = *Pa->Verts[j];
    }
    DEM::Dilation(Vtemp,Pa->EdgeCon,Pa->FaceCon,Vres,Pa->Props.R);
    for (size_t j=0;j<Pa->Verts.Size();j++)
    {
        *Pa->Verts[j] = Vres[j];
    }
    Pa->poly_calc_props(Vres);
    double cube_vol = (4./3.)*PI + 3.*PI*30. + pow(30.0,3.0) + 6*pow(30.0,2.0);


    // add rectangular box
    Vec3_t L(10.0,40.0,90.0);  // An AC Clarke monolith
    d.AddRecBox (-1, Vec3_t(0,20,0),L,1.0,1.0);
    Pa = d.Particles[1];
    Vtemp.Resize(Pa->Verts.Size());
    Vres .Resize(Pa->Verts.Size());
    for (size_t j=0;j<Pa->Verts.Size();j++)
    {
        Vtemp[j] = *Pa->Verts[j];
        Vres [j] = *Pa->Verts[j];
    }
    DEM::Dilation(Vtemp,Pa->EdgeCon,Pa->FaceCon,Vres,Pa->Props.R);
    for (size_t j=0;j<Pa->Verts.Size();j++)
    {
        *Pa->Verts[j] = Vres[j];
    }
    Pa->poly_calc_props(Vres);
    double box_vol = (L(0)+2.0)*(L(1)+2.0)*(L(2)+2.0);
    Quaternion_t temp;
    Conjugate (Pa->Q,temp);
    Pa->Rotate    (temp,Pa->x);

    // add rice
    d.AddRice (-1, Vec3_t(40,0,0),1.,10.,1.,0.,NULL);
    double rice_vol = (4./3.)*PI + PI*10.0;
    Vec3_t rice_I((1./3.)*PI*100+(1./12.)*PI*1000+0.75*PI*10+(8./15.)*PI,
                  (1./3.)*PI*100+(1./12.)*PI*1000+0.75*PI*10+(8./15.)*PI,
                   0.5*PI*10+(8./15.)*PI);

    // add tetrahedron

    //vertices
    Array<Vec3_t> V(4);
    V[0] = -0.019294513013112358, 1.4, 2.989488614614216;
    V[1] = -0.3719223593595582,   1.4, 4.4;
    V[2] =  0.3333333333333335,   1.4, 4.4;
    V[3] = -0.3719223593595582,  -0.7157670780786753, 4.4;
    
     //edges
    Array<Array <int> > E(6);
    for (size_t i=0; i<6; ++i) E[i].Resize(2);
    E[0] = 0, 1;
    E[1] = 1, 2;
    E[2] = 2, 0;
    E[3] = 0, 3;
    E[4] = 1, 3;
    E[5] = 2, 3;

     //face
    Array<Array <int> > F;
    F.Resize(4);
    for (size_t i=0; i<4; ++i) F[i].Resize(3);
    F[0] = 0, 1, 2;
    F[1] = 0, 3, 1;
    F[2] = 0, 2, 3;
    F[3] = 1, 3, 2;
    d.Particles.Push (new DEM::Particle(-1,V,E,F,OrthoSys::O,OrthoSys::O,0.1,1.0));


    Mesh::Generic mesh(/*NDim*/3);
    mesh.SetSize(4,1);
    mesh.SetVert(1,-1, 4.0, 2.0, 4.1);
    mesh.SetVert(0,-1, 0.0,-2.0, 4.1);
    mesh.SetVert(2,-1, 4.0,-2.0, 0.1);
    mesh.SetVert(3,-1, 0.0, 2.0, 0.1);
    mesh.SetCell(0,-1,Array<int>(0,1,2,3));
    d.GenFromMesh (mesh,/*R*/0.1,/*rho*/1.0,true,false);

    d.Particles[0]->Index=0;
    d.Particles[1]->Index=1;
    d.Particles[2]->Index=2;
    d.Particles[3]->Index=3;
    d.Particles[4]->Index=4;

    // initialize
    //d.Particles[0]->Initialize(5000);
    //d.Particles[1]->Initialize(5000);
    //d.Particles[2]->Initialize(5000);
    //d.Particles[3]->Initialize(5000);

    //d.Initialize (/*dt*/0.0);

    // check
    double cube_tol_vol = 1.0;
    double  box_tol_vol = 1.0;
    double rice_tol_vol = 0.1;
    double rice_tol_I   = 2.5;
    double cube_err_vol = fabs(d.Particles[0]->Props.V - cube_vol);
    double box_err_vol  = fabs(d.Particles[1]->Props.V -  box_vol);
    double rice_err_vol = fabs(d.Particles[2]->Props.V - rice_vol);
    double rice_err_I   = fabs(d.Particles[2]->I(0) - rice_I(0)) +
                          fabs(d.Particles[2]->I(1) - rice_I(1)) +
                          fabs(d.Particles[2]->I(2) - rice_I(2));

    // output
    std::cout << "[1;33m\n--- Results ----------------------------------------------------[0m\n";
    cout << "  Cube Volume            = " << d.Particles[0]->Props.V << " (" << cube_vol << ") ==> Error = " << FmtErr(cube_err_vol,cube_tol_vol) << "\n";
    cout << "  Cube Center of mass    = " << d.Particles[0]->x(0)    << ", " << d.Particles[0]->x(1) << ", " << d.Particles[0]->x(2) << "\n";
    cout << "  Cube Moment of inertia = " << d.Particles[0]->I(0)    << ", " << d.Particles[0]->I(1) << ", " << d.Particles[0]->I(2) << "\n";
    cout << "  Cube Quaternion        = " << d.Particles[0]->Q(0)    << ", " << d.Particles[0]->Q(1) << ", " << d.Particles[0]->Q(2) << ", " << d.Particles[0]->Q(3) << "\n";
    cout << endl;
    cout << "  Box  Volume            = " << d.Particles[1]->Props.V << " (" << box_vol  << ") ==> Error = " << FmtErr(box_err_vol,box_tol_vol) << "\n";
    cout << "  Box  Center of mass    = " << d.Particles[1]->x(0)    << ", " << d.Particles[1]->x(1) << ", " << d.Particles[1]->x(2) << "\n";
    cout << "  Box  Moment of inertia = " << d.Particles[1]->I(0)    << ", " << d.Particles[1]->I(1) << ", " << d.Particles[1]->I(2) << "\n";
    cout << "  Box  Quaternion        = " << d.Particles[1]->Q(0)    << ", " << d.Particles[1]->Q(1) << ", " << d.Particles[1]->Q(2) << ", " << d.Particles[1]->Q(3) << "\n";
    cout << endl;
    cout << "  Rice Volume            = " << d.Particles[2]->Props.V << " (" << rice_vol << ") ==> Error = " << FmtErr(rice_err_vol,rice_tol_vol) << "\n";
    cout << "  Rice Center of mass    = " << d.Particles[2]->x(0)    << ", " << d.Particles[2]->x(1) << ", " << d.Particles[2]->x(2) << endl;
    cout << "  Rice Moment of inertia = " << d.Particles[2]->I(0)    << " (" << rice_I(0) << "), " << d.Particles[2]->I(1) << " (" << rice_I(1) << "), " << d.Particles[2]->I(2) << " (" << rice_I(2) << ") ==> Error = " << FmtErr(rice_err_I,rice_tol_I) << "\n";
    cout << "  Rice Quaternion        = " << d.Particles[2]->Q(0)    << ", " << d.Particles[2]->Q(1) << ", " << d.Particles[2]->Q(2) << ", " << d.Particles[2]->Q(3) << "\n";
    cout << endl;
    cout << "  Tetra Volume            = " << d.Particles[3]->Props.V<< " ( is inside? " << d.Particles[3]->IsInside(d.Particles[3]->x) << ")\n";
    cout << "  Tetra Center of mass    = " << d.Particles[3]->x(0)   << ", " << d.Particles[3]->x(1) << ", " << d.Particles[3]->x(2) << "\n";
    cout << "  Tetra Moment of inertia = " << d.Particles[3]->I(0)   << ", " << d.Particles[3]->I(1) << ", " << d.Particles[3]->I(2) << "\n";
    cout << "  Tetra Quaternion        = " << d.Particles[3]->Q(0)   << ", " << d.Particles[3]->Q(1) << ", " << d.Particles[3]->Q(2) << ", " << d.Particles[2]->Q(3) << "\n";
    cout << endl;
    cout << "  Tetra Volume            = " << d.Particles[4]->Props.V<< " ( is inside? " << d.Particles[4]->IsInside(d.Particles[4]->x) << ")\n";
    cout << "  Tetra Center of mass    = " << d.Particles[4]->x(0)   << ", " << d.Particles[4]->x(1) << ", " << d.Particles[4]->x(2) << "\n";
    cout << "  Tetra Moment of inertia = " << d.Particles[4]->I(0)   << ", " << d.Particles[4]->I(1) << ", " << d.Particles[4]->I(2) << "\n";
    cout << "  Tetra Quaternion        = " << d.Particles[4]->Q(0)   << ", " << d.Particles[4]->Q(1) << ", " << d.Particles[4]->Q(2) << ", " << d.Particles[4]->Q(3) << "\n";
    cout << endl;
    
    // draw
    d.WriteXDMF ("test_domain");

    // results
    if ((rice_err_vol>rice_tol_vol) || (rice_err_I>rice_tol_I) || (cube_err_vol>cube_tol_vol)) return 1;
    else return 0;
}
MECHSYS_CATCH
