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

#include <math.h>
#include <random>

#include <gsl/gsl_linalg.h>
// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;


int main(int argc, char **argv) try
{   
    DEM::Domain dom;
    Mesh::Unstructured mesh(/*NDim*/3);
    mesh.Set(12/*number of points to define the goemetry*/,12/*number of faces to define it*/,1/*number of regions*/,0/*number of holes*/);

    double Vmax = 0.001; //maximun volume of cells, controls resolution
    mesh.SetReg(0,0,Vmax,0.01/*Lx*/,0.01/*Ly*/,0.01/*Lz*/); // The values of L must be apoint inside the region to be meshed

    //Defining the points of the geometry
    Vec3_t x00 = Vec3_t(0.0  , 0.0  , 0.0);
    Vec3_t x01 = Vec3_t(1.0  , 0.0  , 0.0);
    Vec3_t x02 = Vec3_t(0.5  , 0.866, 0.0);
    Vec3_t x03 = Vec3_t(0.071, 0.071, 0.0);
    Vec3_t x04 = Vec3_t(0.929, 0.071, 0.0);
    Vec3_t x05 = Vec3_t(0.5  , 0.795, 0.0);
    Vec3_t x06 = Vec3_t(0.0  , 0.0  , 0.5);
    Vec3_t x07 = Vec3_t(1.0  , 0.0  , 0.5);
    Vec3_t x08 = Vec3_t(0.5  , 0.866, 0.5);
    Vec3_t x09 = Vec3_t(0.071, 0.071, 0.5);
    Vec3_t x10 = Vec3_t(0.929, 0.071, 0.5);
    Vec3_t x11 = Vec3_t(0.5  , 0.795, 0.5);
    mesh.SetPnt( 0/*index of the point*/,-1/*tag of the point*/,x00(0),x00(1),x00(2)); 
    mesh.SetPnt( 1/*index of the point*/,-1/*tag of the point*/,x01(0),x01(1),x01(2)); 
    mesh.SetPnt( 2/*index of the point*/,-1/*tag of the point*/,x02(0),x02(1),x02(2)); 
    mesh.SetPnt( 3/*index of the point*/,-1/*tag of the point*/,x03(0),x03(1),x03(2)); 
    mesh.SetPnt( 4/*index of the point*/,-1/*tag of the point*/,x04(0),x04(1),x04(2)); 
    mesh.SetPnt( 5/*index of the point*/,-1/*tag of the point*/,x05(0),x05(1),x05(2)); 
    mesh.SetPnt( 6/*index of the point*/,-1/*tag of the point*/,x06(0),x06(1),x06(2)); 
    mesh.SetPnt( 7/*index of the point*/,-1/*tag of the point*/,x07(0),x07(1),x07(2)); 
    mesh.SetPnt( 8/*index of the point*/,-1/*tag of the point*/,x08(0),x08(1),x08(2)); 
    mesh.SetPnt( 9/*index of the point*/,-1/*tag of the point*/,x09(0),x09(1),x09(2)); 
    mesh.SetPnt(10/*index of the point*/,-1/*tag of the point*/,x10(0),x10(1),x10(2)); 
    mesh.SetPnt(11/*index of the point*/,-1/*tag of the point*/,x11(0),x11(1),x11(2)); 

    //Setting the faces by indexes of the points
    mesh.SetFac( 0, -2, Array<int>(0,1,7,6)/*array of indexes of the points defined before*/);
    mesh.SetFac( 1, -2, Array<int>(0,2,8,6)/*array of indexes of the points defined before*/);
    mesh.SetFac( 2, -2, Array<int>(1,7,8,2)/*array of indexes of the points defined before*/);
    mesh.SetFac( 3, -2, Array<int>(3,4,10,9)/*array of indexes of the points defined before*/);
    mesh.SetFac( 4, -2, Array<int>(4,10,11,5)/*array of indexes of the points defined before*/);
    mesh.SetFac( 5, -2, Array<int>(3,9,11,5)/*array of indexes of the points defined before*/);
    mesh.SetFac( 6, -2, Array<int>(0,3,5,2)/*array of indexes of the points defined before*/);
    mesh.SetFac( 7, -2, Array<int>(0,1,4,3)/*array of indexes of the points defined before*/);
    mesh.SetFac( 8, -2, Array<int>(1,2,5,4)/*array of indexes of the points defined before*/);
    mesh.SetFac( 9, -2, Array<int>(6,9,11,8)/*array of indexes of the points defined before*/);
    mesh.SetFac(10, -2, Array<int>(6,9,10,7)/*array of indexes of the points defined before*/);
    mesh.SetFac(11, -2, Array<int>(7,8,11,10)/*array of indexes of the points defined before*/);

    //Generate the mesh
    mesh.Generate();

    //Translate mesh into DEM particles
    double R = 0.001;          // spheroradius
    double rho = 3000.0;      // density of material
    dom.GenFromMesh (mesh,/*R*/R,/*rho*/rho,true/*Cohesion*/,false);

    //Plot it to test it
    dom.Dilate = true; //Will dilate the sphero-tetrehedra so it will resemble a solid body, this is just for visualization
    dom.WriteXDMF("brick"); 
}
MECHSYS_CATCH
