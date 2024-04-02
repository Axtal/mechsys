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

// GSL
#include <gsl/gsl_linalg.h>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using std::map;
using std::pair;
using DEM::Domain;

int main(int argc, char **argv) try
{
    // camera coordinates
    //double cx=10., cy=10., cz=2.;
    //if (argc>1) cx = atof(argv[1]);
    //if (argc>2) cy = atof(argv[2]);
    //if (argc>3) cz = atof(argv[3]);

    // domain
    Domain dom;
    //dom.CamPos = cx,cy,cz;
    DEM::Particle * p0 = new DEM::Particle(-1, "../tlbm/dolphin", 0.1, 3.0,10.0);
    dom.Particles.Push(p0);

    // nonconvex particle
    //Mesh::Generic mesh(/*NDim*/3);
    //mesh.ReadMesh("nonconvex", /*IsShell*/true); // read file nonconvex.mesh
    //DEM::Particle * p0 = new DEM::Particle(-1, mesh, 0.1, 3.0);
    //dom.Particles.Push (p0);

    // another particle
    //Array<Vec3_t> V;
    //Array< Array<int> > E,F;
    //V.Push (Vec3_t(0.,0.,1.1));
    //V.Push (Vec3_t(1.,0.,1.1));
    //V.Push (Vec3_t(1.,1.,1.1));
    //V.Push (Vec3_t(0.,1.,1.1));
    //V.Push (Vec3_t(0.,0.,1.5));
    //V.Push (Vec3_t(1.,0.,1.5));
    //V.Push (Vec3_t(1.,1.,1.5));
    //V.Push (Vec3_t(0.,1.,1.5));
    //E.Push (Array<int>(0,4));
    //E.Push (Array<int>(1,5));
    //E.Push (Array<int>(2,6));
    //E.Push (Array<int>(3,7));
    //E.Push (Array<int>(4,5));
    //E.Push (Array<int>(5,6));
    //E.Push (Array<int>(6,7));
    //E.Push (Array<int>(7,4));
    //F.Push (Array<int>(4,5,6,7));
    //DEM::Particle * p1 = new DEM::Particle(-1, V,E,F, Vec3_t(0.,0.,0.), Vec3_t(0.,0.,0.), 0.1);
    //dom.Particles.Push (p1);

    // standard particles
    //dom.AddSphere (-1, Vec3_t(2.,0.,0.), /*R*/1.0, /*rho*/1.0);
    //dom.AddCube   (-1, Vec3_t(3.,0.,0.), /*R*/0.1, /*L*/0.5, /*rho*/1.0);
    //dom.AddTetra  (-1, Vec3_t(4.,0.,0.), /*R*/0.1, /*L*/0.5, /*rho*/1.0);
    //dom.AddRice   (-1, Vec3_t(5.,0.,0.), /*R*/0.1, /*L*/0.5, /*rho*/1.0);
    //dom.AddPlane  (-1, Vec3_t(6.,0.,0.), /*R*/0.1, /*Lx*/1.0, /*Ly*/1.0, /*rho*/1.0);

    // domain
    //dom.WritePOV ("nonconvex");
    dom.Save     ("nonconvex");
    dom.WriteXDMF("nonconvex");
    return 0;
}
MECHSYS_CATCH
