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

#ifndef MECHSYS_MONTECARLO_H
#define MECHSYS_MONTECARLO_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

namespace Numerical
{

/** Integration (quadrature) method. */
enum MCMethod
{
    PLAIN, ///< Plain Monte Carlo integration scheme
    MISER, ///< Miser method for strongly changing functions
    VEGAS  ///< Vegas for distribution adapative methods
};

template<typename Instance>
class MonteCarlo
{
public:
    // Constants
    static int WORKSPACE_SIZE; ///< Workspace size

    // Typedefs
    typedef double (Instance::*pFun) (double * x); ///< Callback function

    /** Constructor. */
    MonteCarlo (Instance * p2Inst, MCMethod method=VEGAS, size_t NCalls=500000);

    // Methods
    double Integrate (pFun p2Fun, double * ri, double * rs);           ///< Integrate F(x)
    double CallFun   (double  * r)   { return (_p2inst->*_p2fun)(r); } ///< Call F(x)

private:
    Instance * _p2inst; ///< Pointer to an instance
    pFun       _p2fun;  ///< Pointer to F(x) function
    MCMethod   _method; ///< Method of quadrature
    size_t     _NCALLS; ///< Number of calls for the MC method
};

template<typename Instance>
int MonteCarlo<Instance>::WORKSPACE_SIZE = 3;

/** Trick to pass pointers to member functions to GSL.
 * This will use (*params) to store the pointer to an instance of MonteCarlo, therefore (*params) will no longer be available. */
template<typename Instance>
double __monte_carlo_call_fun__(double * x,size_t dim, void *not_used_for_params)
{
    MonteCarlo<Instance> * p2inst = static_cast<MonteCarlo<Instance>*>(not_used_for_params);
    return p2inst->CallFun(x);
}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

template<typename Instance>
inline MonteCarlo<Instance>::MonteCarlo (Instance * p2Inst, MCMethod Method, size_t NCalls)
    : _p2inst(p2Inst), _method(Method), _NCALLS(NCalls)
{
}

template<typename Instance>
inline double MonteCarlo<Instance>::Integrate(pFun p2Fun, double * ri, double * rs)
{
    // set pointer to function
    _p2fun = p2Fun;

    // function for GSL
    gsl_monte_function f;
    f.f      = &__monte_carlo_call_fun__<Instance>;
    f.dim    = WORKSPACE_SIZE;
    f.params = this;
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);

    // solve
    double result,error;
    switch (_method)
    {
        case PLAIN:
        {
            gsl_monte_plain_state *wsp;
            wsp = gsl_monte_plain_alloc (WORKSPACE_SIZE);
            gsl_monte_plain_integrate(&f,ri,rs,3,_NCALLS,r,wsp,&result,&error);
            gsl_monte_plain_free (wsp);
            break;
        }
        case MISER:
        {
            gsl_monte_miser_state *wsm;
            wsm = gsl_monte_miser_alloc (WORKSPACE_SIZE);
            gsl_monte_miser_integrate(&f,ri,rs,3,_NCALLS,r,wsm,&result,&error);
            gsl_monte_miser_free (wsm);
            break;
        }
        case VEGAS:
        {
            gsl_monte_vegas_state *wsv;
            wsv = gsl_monte_vegas_alloc (WORKSPACE_SIZE);
            gsl_monte_vegas_integrate(&f,ri,rs,3,_NCALLS,r,wsv,&result,&error);
            //do
            //{
                //gsl_monte_vegas_integrate(&f,ri,rs,3,_NCALLS/5,r,wsv,&result,&error);
            //}
            //while(fabs(gsl_monte_vegas_chisq(wsv)-1.0)>0.01);
            gsl_monte_vegas_free (wsv);
            break;
        }
    };
    gsl_rng_free(r);

    // Results
    return result;
}

}; // namespace Numerical

#endif // MECHSYS_MONTECARLO_H
