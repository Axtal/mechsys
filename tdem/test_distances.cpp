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

// Std lib
#include <math.h>

// MechSys
#include <mechsys/dem/graph.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using namespace DEM;

int main(int argc, char **argv) try
{
    double error = 0.0;
    double tol   = 1.0e-15;
    //Periodic cell dimensions
    Vec3_t Per(10.0,10.0,10.0);

    // test vertex-edge distance
    {
        // vertex and edge
        Vec3_t V(21,1,0),V1(1,1,1),V2(0,0,0);
        Vec3_t *a = &V1,*b = &V2;
        Edge   E(a,b);


        // distance
        Vec3_t Xi,Xf,S;
        Distance (V,E,Xi,Xf,S,Per);

        // Report the points and branch vector
        cout << "Vertex-Edge Xi = " << Xi << std::endl;
        cout << "Vertex-Edge Xf = " << Xf << std::endl;
        cout << "Vertex-Edge S  = " << S  << std::endl;

        // check
        double d = Distance(V,E,Per);
        error += fabs(d - sqrt(2.0/3.0));
    }

    // test edge-edge distance
    {
        // two edges
        //Vec3_t a(-11.0,-1.0,11.0),b(-9.0,1.0,11.0),c(-1.0,1.0,0.0),d(1.0,-1.0,0.0);
        Vec3_t a(11.0,0.0,0.0),b(0.0,0.0,0.0),c(2.0,-2.0,0.0),d(2.0,2.0,0.0);
        Edge E1(a,b);
        Edge E2(c,d);

        // distance
        Vec3_t Xi,Xf,S;
        Distance (E1,E2,Xi,Xf,S,Per);

        // Report the points and branch vector
        cout << "Edge-Edge Xi = " << Xi << std::endl;
        cout << "Edge-Edge Xf = " << Xf << std::endl;
        cout << "Edge-Edge S  = " << S  << std::endl;

        error += fabs(Distance(E1,E2,Per)-1.0);
    }

    // test vertex-face distance
    {
        // vertex
        Vec3_t V(10.5,10.5,1.0);

        // face
        Array<Vec3_t> C(4); // connectivity
        C[0] = -1.0,-1.0,0.0;
        C[1] = -1.0, 1.0,0.0;
        C[2] =  1.0, 1.0,0.0;
        C[3] =  1.0,-1.0,0.0;

        Face F(C);

        // distance
        Vec3_t Xi,Xf,S;
        Distance (V,F,Xi,Xf,S,Per);

        // Report the points and branch vector
        cout << "Vertex-Face Xi = " << Xi << std::endl;
        cout << "Vertex-Face Xf = " << Xf << std::endl;
        cout << "Vertex-Face S  = " << S  << std::endl;

        error += fabs(Distance(V,F,Per)-1.0);
    }

    // test torus-vertex distance
    {
        // vertex
        Vec3_t V(10.0,  10.0, 0.0);

        // torus
        Vec3_t P0( 0.0, 0.0, 0.0);
        Vec3_t P1( 1.0, 0.0, 0.0);
        Vec3_t P2( 0.0, 0.0, 1.0);

        Torus T( P0, P1, P2);

        // distance
        Vec3_t Xi,Xf,s;
        Distance (V,T,Xi,Xf,s,OrthoSys::O);

        error += fabs(Distance(V,T,OrthoSys::O)-sqrt(100+81));
    }

    // test cylinder-vertex distance
    {
        // vertex
        Vec3_t V(-10.0, 0.5, -10.0);

        // torus
        Vec3_t P00( 0.0, 0.0, 0.0);
        Vec3_t P01( 1.0, 0.0, 0.0);
        Vec3_t P02( 0.0, 0.0, 1.0);
        Vec3_t P03( 0.0, 0.0,-1.0);

        Torus T0( P00, P01, P02);

        //Vec3_t P10( 0.0, 1.0, 0.0);
        //Vec3_t P11( 0.5, 1.0, 0.0);
        //Vec3_t P12( 0.0, 1.0, 0.5);
        //Vec3_t P13( 0.0, 1.0,-0.5);
        Vec3_t P10( 0.0, 1.0, 0.0);
        Vec3_t P11( 1.0, 1.0, 0.0);
        Vec3_t P12( 0.0, 1.0, 1.0);
        Vec3_t P13( 0.0, 1.0,-1.0);

        Torus T1( P10, P11, P12);

        Cylinder C(T0,T1,P03,P13);

        // distance
        Vec3_t Xi,Xf,s;
        Distance(V,C,Xi,Xf,s,OrthoSys::O);

        error += fabs(Distance(V,C,OrthoSys::O)-sqrt(200)+1);
    }

    cout << "error = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    return (error>tol ? 1 : 0);
}
MECHSYS_CATCH
