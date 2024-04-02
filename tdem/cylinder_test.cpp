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
#include <mechsys/linalg/matvec.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/mesh/structured.h>

using std::cout;
using std::endl;
using std::ofstream;
using DEM::Domain;

int main(int argc, char **argv) try
{
    if (argc!=2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".hdf5");
     //set the simulation domain ////////////////////////////////////////////////////////////////////////////
    //
    Domain d;
    Mesh::Unstructured mesh(/*NDim*/2);
    mesh.Set    (8, 8, 1, 1);            // 8 points, 8 segments, 1 region, 1 hole
    //mesh.SetReg (0, -1, -1.0, 0.2, 0.2); // id, tag, max{area}, x, y <<<<<<< regions
    //mesh.SetHol (0, 2.5, 2.5);           // id, x, y <<<<<<< holes
    //mesh.SetPnt (0, -1, 0.0, 0.0);       // id, vtag, x, y <<<<<< points
    //mesh.SetPnt (1, -2, 5.0, 0.0);       // id, vtag, x, y
    //mesh.SetPnt (2, -3, 5.0, 5.0);       // id, vtag, x, y
    //mesh.SetPnt (3, -4, 0.0, 5.0);       // id, vtag, x, y
    //mesh.SetPnt (4,  0, 2.0, 2.0);       // id, vtag, x, y
    //mesh.SetPnt (5,  0, 3.0, 2.0);       // id, vtag, x, y
    //mesh.SetPnt (6,  0, 3.0, 3.0);       // id, vtag, x, y
    //mesh.SetPnt (7,  0, 2.0, 3.0);       // id, vtag, x, y
    //mesh.SetSeg (0, -10,  0, 1);         // id, etag, L, R <<<<<<<<<<<< segments
    //mesh.SetSeg (1, -20,  1, 2);         // id, etag, L, R
    //mesh.SetSeg (2, -30,  2, 3);         // id, etag, L, R
    //mesh.SetSeg (3, -40,  3, 0);         // id, etag, L, R
    //mesh.SetSeg (4,   0,  4, 5);         // id, etag, L, R
    //mesh.SetSeg (5,   0,  5, 6);         // id, etag, L, R
    //mesh.SetSeg (6,   0,  6, 7);         // id, etag, L, R
    //mesh.SetSeg (7,   0,  7, 4);         // id, etag, L, R
    //mesh.Generate ();
//
    //d.GenFromMesh(mesh,/*spheroradius*/0.1,/*density*/1.0,/*iscohesive*/true,/*montecarlo mass properties*/false,/*thickness*/2.0);
    //d.Center(Vec3_t(0.0,8.0,0.0));
//
    //d.GenFromMesh(mesh,/*spheroradius*/0.1,/*density*/1.0,/*iscohesive*/true,/*montecarlo mass properties*/false,/*thickness*/2.0);

    
    d.AddVoroPack (-1, 0.1, 4,4,4, 4,4,4, 3.0, true, true, 1000, 1.0,0.0);
    d.Center(Vec3_t(0.0,4.0,0.0));
    d.AddPlane(/*tag*/-2,/*position*/Vec3_t(0.0,-0.2,0.0),/*spheroradius*/0.2,/*Lx*/100,/*Ly*/100,/*rho*/1.0,/*angle*/M_PI/2.0,/*axis*/&OrthoSys::e0);
    d.Initialize();
    d.Save(filekey.CStr());
    //Array<double> X,Y,D;
    //d.GetGSD(X,Y,D);

    double Kn          = 1.0e5;
    double Kt          = 3.3e4;
    double Bn          = 1.0e5;
    double Bt          = 3.3e4;
    double Bm          = 3.3e4;
    double Gn          = 16.0;
    double Gt          = 8.0;
    double eps         = 0.08;

    Dict B;
    B.Set(-1,"Bn Bt Bm Gn Gt eps Kn Kt",Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt);
    B.Set(-2,"Bn Bt Bm Gn Gt eps Kn Kt",Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt);
    d.SetProps(B);
    //ofstream fg("granulometry.txt");
    //for (size_t i = 0;i < X.Size() ; i++ )
    //{
        //fg << Util::_10_6 << X[i] << Util::_8s << Y[i] << std::endl;
    //}
    //fg.close();


    //Fix the plane
    DEM::Particle * p = d.GetParticle(-2,true);
    p->FixVeloc();

     //Initialize the gravity on the particles
    for (size_t i=0;i<d.Particles.Size();i++)
    {
        d.Particles[i]->Ff = d.Particles[i]->Props.m*Vec3_t(0.0,-9.8,0.0);
    }

    d.Solve     (/*tf*/10.0, /*dt*/0.00005, /*dtOut*/0.1, NULL, NULL, filekey.CStr());

    return 0;
}
MECHSYS_CATCH
