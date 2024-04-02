/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
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

#ifndef MECHSYS_ROOT_H
#define MECHSYS_ROOT_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// MechSys
#include <mechsys/util/fatal.h>
#include <mechsys/util/string.h>

namespace Numerical
{

template<typename Instance>
class Root
{
public:
    // Typedefs
    typedef double (Instance::*pF)    (double x, void * UserData); ///< Callback function
    typedef double (Instance::*pdFdx) (double x, void * UserData); ///< Callback derivative of function w.r.t x

    /** Constructor. */
    Root (Instance * Inst, pF F, pdFdx dFdx=NULL);

    // Methods
    double Solve (double xa, double xb, double * x_guess=NULL, void * UserData=NULL); ///< Find root

    // Data
    double Tol;
    int    It;
    int    MaxIt;
    String Scheme;
    bool   Verbose;

private:
    Instance * _pInst; ///< Pointer to object
    pF         _pF;    ///< Pointer to function
    pdFdx      _pdFdx; ///< Pointer to derivative
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Instance>
inline Root<Instance>::Root (Instance * Inst, pF F, pdFdx dFdx)
    : Tol     (sqrt(DBL_EPSILON)),
      It      (0),
      MaxIt   (30),
      Scheme  ("Brent"),
      Verbose (false),
      _pInst  (Inst),
      _pF     (F),
      _pdFdx  (dFdx)
{
}

template<typename Instance>
inline double Root<Instance>::Solve (double xa, double xb, double * x_guess, void * UserData)
{
    if (Scheme=="Brent")
    {
        /* Based on ZEROIN C math library: http://www.netlib.org/c/

           Algorithm:
           G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
           computations. M., Mir, 1980, p.180 of the Russian edition

           The function makes use of the bissection procedure combined with
           the linear or quadric inverse interpolation.
           At every step program operates on three abscissae: a, b, and c:
            a: the last but one approximation
            b: the last and the best approximation to the root
            c: the last but one or even earlier approximation than a such that:
                1) |f(b)| <= |f(c)|
                2) f(b) and f(c) have opposite signs, i.e. b and c confine the root

           At every step this algorithm selects one of the two new approximations, the former being
           obtained by the bissection procedure and the latter resulting in the interpolation.
           If a,b and c are all different, the quadric interpolation is utilized, otherwise the linear one.
           If the latter (i.e. obtained by the interpolation) point is reasonable (i.e. lies within
           the current interval [b,c] not being too close to the boundaries) it is accepted.
           The bissection result is used in the other case. Therefore, the range of uncertainty is
           ensured to be reduced at least by the factor 1.6.

           Examples: <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/numerical/tst/tbrentroot.cpp?view=markup"> ttbrentroot.cpp Test Brent's method</a>
        */

        double a  = xa; // the last but one approximation
        double b  = xb; // the last and the best approximation to the root
        double c  = a;  // the last but one or even earlier approximation than a that
        double fa = (_pInst->*_pF) (a, UserData);
        double fb = (_pInst->*_pF) (b, UserData);
        double fc = fa;

        // Check input
        if ((fa>0.0 && fb>0.0) || (fa<0.0 && fb<0.0)) throw new Fatal("Root::Solve: Brent method: Root must be bracketed.");

        // Solve
        for (It=0; It<MaxIt; ++It)
        {
            // Distance from the last but one to the last approximation
            double prev_step = b-a;

            // Swap data for b to be the best approximation
            if (fabs(fc)<fabs(fb))
            {
                 a =  b;     b =  c;     c =  a;
                fa = fb;    fb = fc;    fc = fa;
            }
            double tol_act  = 2.0*DBL_EPSILON*fabs(b) + Tol/2.0;  // Actual tolerance
            double new_step = (c-b)/2.0;                          // Step at this iteration

            // Check for convergence
            if (Verbose) printf("It=%d, xnew=%g, f(xnew)=%g, err=%g\n", It, b, fb, fabs(new_step));
            if (fabs(new_step)<=tol_act || fb==0.0) return b; // Acceptable approx. is found

            // Decide if the interpolation can be tried
            if (fabs(prev_step)>=tol_act && fabs(fa)>fabs(fb)) // If prev_step was large enough and was in true direction
            {
                // Interpolation may be tried
                double p,q; // Interpolation step is calculated in the form p/q; division operations is delayed until the last moment
                double t1,cb,t2;
                cb = c-b;
                if(a==c) // If we have only two distinct points, linear interpolation can only be applied
                {
                    t1 = fb/fa;
                    p  = cb*t1;
                    q  = 1.0-t1;
                }
                else // Quadric inverse interpolation
                {
                    q = fa/fc;     t1 = fb/fc;     t2 = fb/fa;
                    p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
                    q = (q-1.0) * (t1-1.0) * (t2-1.0);
                }

                // p was calculated with the oposite sign; make p positive and assign possible minus to q
                if (p>0.0) q = -q;
                else       p = -p;

                // If b+p/q falls in [b,c] and isn't too large
                if (p<(0.75*cb*q-fabs(tol_act*q)/2.0) && p<fabs(prev_step*q/2.0))
                    new_step = p/q;// it is accepted

                // If p/q is too large then the bissection procedure can reduce [b,c] range to more extent
            }

            // Adjust the step to be not less than tolerance
            if (fabs(new_step)<tol_act)
            {
                if (new_step>0.0) new_step =  tol_act;
                else              new_step = -tol_act;
            }

            // Save the previous approx
             a =  b;
            fa = fb;

            // Do step to a new approxim.
            b += new_step;
            fb = (_pInst->*_pF) (b, UserData);

            // Adjust c for it to have a sign opposite to that of b
            if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0))
            {
                 c =  a;
                fc = fa;
            }
        }

        // Did not converge
        throw new Fatal("Root::Solve: Brent method did not converge after %d iterations",It);
    }
    else if (Scheme=="Newton")
    {
        if (_pdFdx==NULL) throw new Fatal("Root:Solve: Newton method: For this method, dFdx cannot be NULL");
        double fa = (_pInst->*_pF) (xa, UserData);
        double fb = (_pInst->*_pF) (xb, UserData);
        if (fabs(xb-xa)<Tol) throw new Fatal("Root:Solve: Newton method: |xb-xa| must be greater than %g",Tol);
        if (fabs(fb-fa)<Tol) throw new Fatal("Root:Solve: Newton method: |fb-fa| must be greater than %g",Tol);
        double alp  = (x_guess==NULL ? fa/(fa-fb) : ((*x_guess)-xa)/(xb-xa));
        double xsol = xa + alp*(xb-xa);
        for (It=0; It<MaxIt; ++It)
        {
            double F    = (_pInst->*_pF)    (xsol, UserData);
            double dFdx = (_pInst->*_pdFdx) (xsol, UserData);
            double xnew = xsol - F/dFdx;
            double err  = fabs((xnew-xsol)/(1.0+xnew));
            if (Verbose) printf("alp=%g, xsol=%g, F=%g, dFdx=%g, F/dFdx=%g, xnew=%g, err=%g\n", alp, xsol, F, dFdx, F/dFdx, xnew, err);
            if (err<Tol) return xnew;
            xsol = xnew;
        }

        // Did not converge
        throw new Fatal("Root::Solve: Newton-Rhapson method did not converge after %d iterations",It);
    }
    else throw new Fatal("Root::Solve: scheme==%s is not available",Scheme.CStr());
    return 0;
}

}; // namespace Numerical

#endif // MECHSYS_ROOT_H
