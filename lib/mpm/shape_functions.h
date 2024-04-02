/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang                                         *
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


#ifndef MECHSYS_MPM_SHAPE_H
#define MECHSYS_MPM_SHAPE_H

namespace MPM
{
//GIMP1D shape function, it takes the position of the particle x, the position of the node xc, the grid size lx, the particle size lpx and returns the
//value of the function in n and its gradient in gn
void GIMP(double x, double xc, double lx, double lpx, double & n, double & gn) 
{
	n = 0.;
	gn = 0.;
	double d = x-xc;
    double da = fabs(d);
    if (da<0.5*lpx)
    {   
        n  = 1.0-(4.0*d*d+lpx*lpx)/(4.0*lx*lpx);
        gn = -2.0*d/(lx*lpx);
    }
    else if (da<lx-0.5*lpx)
    {
        n  = 1.0 - da/lx;
        gn = -d/(da*lx);
    }
    else if (da<lx+0.5*lpx)
    {
        n  = (lx + 0.5*lpx - da)*(lx + 0.5*lpx - da)/(2.0*lx*lpx);
        gn = -d*(lx + 0.5*lpx - da)/(da*lx*lpx);
    }
}

//GIMP3D uses GIMP1D to produce the 3D vercion
void GIMP3D(Vec3_t const & x, Vec3_t const & xc, double lx, Vec3_t const & lp, double & n, Vec3_t & gn)
{
	double n0 =  0.0;
	double n1 =  0.0;
	double n2 =  0.0;
	double gn0 = 0.0;
	double gn1 = 0.0;
	double gn2 = 0.0;

	GIMP(x(0), xc(0), lx, lp(0), n0, gn0);
	GIMP(x(1), xc(1), lx, lp(1), n1, gn1);
	GIMP(x(2), xc(2), lx, lp(2), n2, gn2);

	n = n0*n1*n2;

	gn(0) = gn0* n1* n2;
	gn(1) =  n0*gn1* n2;
	gn(2) =  n0* n1*gn2;
}

// Traditional shape function
double Shape(double x, double xc, double lx)
{
    double sf = 0.0;
    double da = fabs(x - xc);
    if (da<lx) sf = 1.0 - da/lx;
    return sf;
}

double Shape3D(Vec3_t const & x, Vec3_t const & xc, double lx)
{
    return Shape(x(0),xc(0),lx)*Shape(x(1),xc(1),lx)*Shape(x(2),xc(2),lx);
}

}
#endif //MECHSYS_MPM_SHAPE_H
