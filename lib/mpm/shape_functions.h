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
enum ShapeFunctionType
{
    SHAPE_GIMP               = 0,
    SHAPE_BSPLINE_QUADRATIC  = 1,
    SHAPE_BSPLINE_CUBIC      = 2
};

// GIMP reference: Bardenhagen, S. G. and Kober, E. M. (2004), "The generalized interpolation material point method", CMES 5(6):477-495.
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

// Quadratic cardinal B-spline reference: de Boor, C. (1978), "A Practical Guide to Splines"; MPM use: Steffen, Kirby and Berzins (2008), "Analysis and reduction of quadrature errors in the material point method (MPM)".
//Quadratic B-spline shape function on a uniform grid.  It is C1-continuous and has support 1.5 cells wide in each direction.
void BSplineQuadratic(double x, double xc, double lx, double & n, double & gn)
{
    n = 0.0;
    gn = 0.0;
    double r  = (x-xc)/lx;
    double ar = fabs(r);

    if (ar<0.5)
    {
        n  = 0.75 - r*r;
        gn = -2.0*r/lx;
    }
    else if (ar<1.5)
    {
        double a = 1.5 - ar;
        n  = 0.5*a*a;
        gn = -(r/ar)*a/lx;
    }
}

// Cubic cardinal B-spline reference: de Boor, C. (1978), "A Practical Guide to Splines"; MPM use: Stomakhin, Schroeder, Chai, Teran and Selle (2013), "A material point method for snow simulation".
void BSplineCubic(double x, double xc, double lx, double & n, double & gn)
{
    n = 0.0;
    gn = 0.0;
    double r  = (x-xc)/lx;
    double ar = fabs(r);

    if (ar<1.0)
    {
        n  = 2.0/3.0 - r*r + 0.5*ar*ar*ar;
        gn = (-2.0*r + 1.5*r*ar)/lx;
    }
    else if (ar<2.0)
    {
        double a = 2.0 - ar;
        n  = a*a*a/6.0;
        gn = -0.5*(r/ar)*a*a/lx;
    }
}

//GIMP3D uses GIMP1D to produce the 3D version
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

//BSplineQuadratic3D uses BSplineQuadratic1D to produce the 3D version
void BSplineQuadratic3D(Vec3_t const & x, Vec3_t const & xc, double lx, double & n, Vec3_t & gn)
{
	double n0 =  0.0;
	double n1 =  0.0;
	double n2 =  0.0;
	double gn0 = 0.0;
	double gn1 = 0.0;
	double gn2 = 0.0;

	BSplineQuadratic(x(0), xc(0), lx, n0, gn0);
	BSplineQuadratic(x(1), xc(1), lx, n1, gn1);
    BSplineQuadratic(x(2), xc(2), lx, n2, gn2);

	n = n0*n1*n2;

	gn(0) = gn0* n1* n2;
	gn(1) =  n0*gn1* n2;
	gn(2) =  n0* n1*gn2;
}

//BSplineCubic3D uses BSplineCubic1D to produce the 3D version
void BSplineCubic3D(Vec3_t const & x, Vec3_t const & xc, double lx, double & n, Vec3_t & gn)
{
	double n0 =  0.0;
	double n1 =  0.0;
	double n2 =  0.0;
	double gn0 = 0.0;
	double gn1 = 0.0;
	double gn2 = 0.0;

	BSplineCubic(x(0), xc(0), lx, n0, gn0);
	BSplineCubic(x(1), xc(1), lx, n1, gn1);
    BSplineCubic(x(2), xc(2), lx, n2, gn2);

	n = n0*n1*n2;

	gn(0) = gn0* n1* n2;
	gn(1) =  n0*gn1* n2;
	gn(2) =  n0* n1*gn2;
}

double ShapeSupportRadius(int ShapeFunction, double lx, double lpx)
{
    if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC) return 1.5*lx;
    if (ShapeFunction==SHAPE_BSPLINE_CUBIC    ) return 2.0*lx;
    return lx + 0.5*lpx;
}

void ShapeNodeRange(int ShapeFunction, double x, double xlow, double lx, double lpx, size_t nmax, size_t & nlow, size_t & nhigh)
{
    double lo = 0.0;
    double hi = 0.0;
    if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC)
    {
        lo = floor((x-xlow)/lx - 0.5);
        hi = lo + 2.0;
    }
    else if (ShapeFunction==SHAPE_BSPLINE_CUBIC)
    {
        lo = floor((x-xlow)/lx) - 1.0;
        hi = lo + 3.0;
    }
    else
    {
        double r = ShapeSupportRadius(ShapeFunction,lx,lpx);
        lo = floor((x-r-xlow)/lx);
        hi = ceil ((x+r-xlow)/lx);
    }
    nlow  = (size_t)std::min(std::max(lo,0.0),(double)nmax-1.0);
    nhigh = (size_t)std::min(std::max(hi,0.0),(double)nmax-1.0);
}

void ShapeNodeRangeRaw(int ShapeFunction, double x, double xlow, double lx, double lpx, int & nlow, int & nhigh)
{
    double lo = 0.0;
    double hi = 0.0;
    if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC)
    {
        lo = floor((x-xlow)/lx - 0.5);
        hi = lo + 2.0;
    }
    else if (ShapeFunction==SHAPE_BSPLINE_CUBIC)
    {
        lo = floor((x-xlow)/lx) - 1.0;
        hi = lo + 3.0;
    }
    else
    {
        double r = ShapeSupportRadius(ShapeFunction,lx,lpx);
        lo = floor((x-r-xlow)/lx);
        hi = ceil ((x+r-xlow)/lx);
    }
    nlow  = (int)lo;
    nhigh = (int)hi;
}

inline size_t PeriodicNodeIndex(int n, size_t nmax)
{
    int m = n%(int)nmax;
    return (size_t)((m<0) ? m + (int)nmax : m);
}


// Precomputed separable 1-D shape data for one particle.  The 3-D values are
// formed with exactly the same multiplication order as GIMP3D/BSpline*3D, so
// the cached path is bitwise identical to the direct one.
struct ShapeCache1D
{
    double w[8];
    double g[8];
};

inline bool FillShapeCache3D(int ShapeFunction, Vec3_t const & x, Vec3_t const & xlow,
                             double lx, Vec3_t const & lp,
                             int const nlow[3], int const nhigh[3], ShapeCache1D c[3])
{
    for (int d=0; d<3; d++)
    {
        int span = nhigh[d]-nlow[d]+1;
        if (span<0 || span>8) return false;
        for (int k=0; k<span; k++)
        {
            double xc = (nlow[d]+k)*lx + xlow(d);
            if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC)
                BSplineQuadratic(x(d),xc,lx,c[d].w[k],c[d].g[k]);
            else if (ShapeFunction==SHAPE_BSPLINE_CUBIC)
                BSplineCubic(x(d),xc,lx,c[d].w[k],c[d].g[k]);
            else
                GIMP(x(d),xc,lx,lp(d),c[d].w[k],c[d].g[k]);
        }
    }
    return true;
}

void MPMShapeFunction3D(int ShapeFunction, Vec3_t const & x, Vec3_t const & xc, double lx, Vec3_t const & lp, double & n, Vec3_t & gn)
{
    if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC)
    {
        BSplineQuadratic3D(x,xc,lx,n,gn);
        return;
    }
    if (ShapeFunction==SHAPE_BSPLINE_CUBIC)
    {
        BSplineCubic3D(x,xc,lx,n,gn);
        return;
    }
    GIMP3D(x,xc,lx,lp,n,gn);
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

#ifdef USE_CUDA

__device__ __inline__ void GIMPCU(real x, real xc, real lx, real lpx, real & n, real & gn)
{
	n = 0.;
	gn = 0.;
	real d = x-xc;
    real da = fabs(d);
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

__device__ __inline__ void GIMP3DCU(real3 const & x, real3 const & xc, real lx, real3 const & lp, real & n, real3 & gn)
{
	real n0 =  0.0;
	real n1 =  0.0;
	real n2 =  0.0;
	real gn0 = 0.0;
	real gn1 = 0.0;
	real gn2 = 0.0;

	GIMPCU(x.x, xc.x, lx, lp.x, n0, gn0);
	GIMPCU(x.y, xc.y, lx, lp.y, n1, gn1);
	GIMPCU(x.z, xc.z, lx, lp.z, n2, gn2);

	n = n0*n1*n2;

	gn.x = gn0* n1* n2;
	gn.y =  n0*gn1* n2;
	gn.z =  n0* n1*gn2;
}

__device__ __inline__ void BSplineQuadraticCU(real x, real xc, real lx, real & n, real & gn)
{
    n = 0.0;
    gn = 0.0;
    real r  = (x-xc)/lx;
    real ar = fabs(r);

    if (ar<0.5)
    {
        n  = 0.75 - r*r;
        gn = -2.0*r/lx;
    }
    else if (ar<1.5)
    {
        real a = 1.5 - ar;
        n  = 0.5*a*a;
        gn = -(r/ar)*a/lx;
    }
}

__device__ __inline__ void BSplineCubicCU(real x, real xc, real lx, real & n, real & gn)
{
    n = 0.0;
    gn = 0.0;
    real r  = (x-xc)/lx;
    real ar = fabs(r);

    if (ar<1.0)
    {
        n  = 2.0/3.0 - r*r + 0.5*ar*ar*ar;
        gn = (-2.0*r + 1.5*r*ar)/lx;
    }
    else if (ar<2.0)
    {
        real a = 2.0 - ar;
        n  = a*a*a/6.0;
        gn = -0.5*(r/ar)*a*a/lx;
    }
}

__device__ __inline__ void BSplineQuadratic3DCU(real3 const & x, real3 const & xc, real lx, real & n, real3 & gn)
{
	real n0 =  0.0;
	real n1 =  0.0;
	real n2 =  0.0;
	real gn0 = 0.0;
	real gn1 = 0.0;
	real gn2 = 0.0;

	BSplineQuadraticCU(x.x, xc.x, lx, n0, gn0);
	BSplineQuadraticCU(x.y, xc.y, lx, n1, gn1);
    BSplineQuadraticCU(x.z, xc.z, lx, n2, gn2);

	n = n0*n1*n2;

	gn.x = gn0* n1* n2;
	gn.y =  n0*gn1* n2;
	gn.z =  n0* n1*gn2;
}

__device__ __inline__ void BSplineCubic3DCU(real3 const & x, real3 const & xc, real lx, real & n, real3 & gn)
{
	real n0 =  0.0;
	real n1 =  0.0;
	real n2 =  0.0;
	real gn0 = 0.0;
	real gn1 = 0.0;
	real gn2 = 0.0;

	BSplineCubicCU(x.x, xc.x, lx, n0, gn0);
	BSplineCubicCU(x.y, xc.y, lx, n1, gn1);
    BSplineCubicCU(x.z, xc.z, lx, n2, gn2);

	n = n0*n1*n2;

	gn.x = gn0* n1* n2;
	gn.y =  n0*gn1* n2;
	gn.z =  n0* n1*gn2;
}

__device__ __inline__ real ShapeSupportRadiusCU(int ShapeFunction, real lx, real lpx)
{
    if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC) return 1.5*lx;
    if (ShapeFunction==SHAPE_BSPLINE_CUBIC    ) return 2.0*lx;
    return lx + 0.5*lpx;
}

__device__ __inline__ void ShapeNodeRangeCU(int ShapeFunction, real x, real xlow, real lx, real lpx, size_t nmax, size_t & nlow, size_t & nhigh)
{
    real lo = 0.0;
    real hi = 0.0;
    if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC)
    {
        lo = floor((x-xlow)/lx - 0.5);
        hi = lo + 2.0;
    }
    else if (ShapeFunction==SHAPE_BSPLINE_CUBIC)
    {
        lo = floor((x-xlow)/lx) - 1.0;
        hi = lo + 3.0;
    }
    else
    {
        real r = ShapeSupportRadiusCU(ShapeFunction,lx,lpx);
        lo = floor((x-r-xlow)/lx);
        hi = ceil ((x+r-xlow)/lx);
    }
    nlow  = (size_t)fmin(fmax(lo,0.0),(real)nmax-1.0);
    nhigh = (size_t)fmin(fmax(hi,0.0),(real)nmax-1.0);
}

__device__ __inline__ void ShapeNodeRangeRawCU(int ShapeFunction, real x, real xlow, real lx, real lpx, int & nlow, int & nhigh)
{
    real lo = 0.0;
    real hi = 0.0;
    if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC)
    {
        lo = floor((x-xlow)/lx - 0.5);
        hi = lo + 2.0;
    }
    else if (ShapeFunction==SHAPE_BSPLINE_CUBIC)
    {
        lo = floor((x-xlow)/lx) - 1.0;
        hi = lo + 3.0;
    }
    else
    {
        real r = ShapeSupportRadiusCU(ShapeFunction,lx,lpx);
        lo = floor((x-r-xlow)/lx);
        hi = ceil ((x+r-xlow)/lx);
    }
    nlow  = (int)lo;
    nhigh = (int)hi;
}

__device__ __inline__ size_t PeriodicNodeIndexCU(int n, size_t nmax)
{
    int m = n%(int)nmax;
    return (size_t)((m<0) ? m + (int)nmax : m);
}

__device__ __inline__ void MPMShapeFunction3DCU(int ShapeFunction, real3 const & x, real3 const & xc, real lx, real3 const & lp, real & n, real3 & gn)
{
    if (ShapeFunction==SHAPE_BSPLINE_QUADRATIC)
    {
        BSplineQuadratic3DCU(x,xc,lx,n,gn);
        return;
    }
    if (ShapeFunction==SHAPE_BSPLINE_CUBIC)
    {
        BSplineCubic3DCU(x,xc,lx,n,gn);
        return;
    }
    GIMP3DCU(x,xc,lx,lp,n,gn);
}
#endif

}
#endif //MECHSYS_MPM_SHAPE_H
