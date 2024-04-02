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

// Std Lib
#include <iostream>

// MechSys
#include <mechsys/lbm/dem.h>

using std::cout;
using std::endl;

/* Flow past a cylinder with obstacle */

int main(int argc, char **argv) try
{
    // constants
    double u_max  = 0.1;                 // Poiseuille's maximum velocity
    double Re     = 100;                 // Reynold's number
    int    nx     = 400;                 // cell dimension
    int    ny     = 100;                 // cell dimension
    int    radius = ny/10 + 1;           // radius of inner circle (obstacle)
	double h      = 1;                   // grid space
	double dt     = 1;                   // timestep
    double nu     = u_max*(2*radius)/Re; // viscocity

	// allocate lattice
	LBM::Lattice l(/*filekey*/"cylinder", /*3d?*/false, nu, nx, ny, /*nz*/1, dt, h);             
					
	// set walls (top and bottom)
    l.SetTopSolid    ();
    l.SetBottomSolid ();

	// set inner obstacle
	int obsX = ny/2;   // x position
	int obsY = ny/2+3; // y position
    LBM::Disk Ball (Vec3_t(obsX, obsY, 0.0), Vec3_t(0.0,0.0,0.0), radius, 100.0, 1000.0, dt);
	Ball.DrawDisk  (l, dt);

	// define boundary conditions
	for (size_t j=0; j<l.Ny(); j++)
	{
        // set parabolic profile
        double L  = ny - 2;                       // channel width in cell units
        double yp = j - 1.5;                      // ordinate of cell
        double vx = u_max*4/(L*L)*(L*yp - yp*yp); // horizontal velocity
        double vy = 0.0;                          // vertical velocity
		Vec3_t v(vx, vy, 0.0);                    // velocity vector
		l.SetVelocityBC (0, j, v);                // set velocity BC to the left-side
		l.SetDensityBC  (nx-1, j, 1.0);           // set density BC to the right-side
	}

	// define initial conditions: velocity speed and density
	for (size_t i=0; i<l.Nx(); i++)
	for (size_t j=0; j<l.Ny(); j++)
	{
		double rho0 = 1.0;
		Vec3_t v0(0.08, 0.0, 0.0);
		l.GetCell(i,j)->Initialize (rho0, v0, l.Cs());
	}

	// solve
	l.Solve(/*tIni*/0.0, /*tFin*/5000.0,/*dtOut*/10.0);
}
MECHSYS_CATCH
