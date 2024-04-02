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

#ifndef MECHSYS_NUMDIFF_H
#define MECHSYS_NUMDIFF_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

// MechSys
#include <mechsys/util/string.h>
#include <mechsys/util/fatal.h>

namespace Numerical
{

template<typename Instance>
class Diff // numerical differentiation
{
public:
    // Typedefs
    typedef double (Instance::*pFun) (double x); ///< Callback function: y=f(x)

    // Constructor
    Diff (Instance * p2Inst);

    // Methods
    double DyDx (pFun p2Fun, double AtX, double StepSize=1.0e-8); ///< Evaluate dy_dx at x

    // Internal methods
    double CallFun (double x) const { return (_p2inst->*_p2fun)(x); }

    // Data
    int    Method;     ///< -1=backward, 0=central, 1=forward
    double LastAbsErr; ///< Last Absolute Error

private:
    Instance * _p2inst;  ///< Pointer to an instance
    pFun       _p2fun;   ///< Pointer to instance function
};

/** Trick to pass pointers to member functions to GSL.
 * This will use (*params) to store the pointer to an instance of Diff, therefore (*params) will no longer be available. */
template<typename Instance>
double __numdiff_call_fun__ (double x, void * not_used_for_params)
{
    Diff<Instance> * p2inst = static_cast<Diff<Instance>*>(not_used_for_params);
    return p2inst->CallFun (x);
}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Instance>
inline Diff<Instance>::Diff (Instance * p2Inst)
    : Method     (0),
      LastAbsErr (0.),
      _p2inst    (p2Inst)
{}

template<typename Instance>
inline double Diff<Instance>::DyDx (pFun p2Fun, double AtX, double h)
{
    // set pointer to function
    _p2fun = p2Fun;

    // set GSL function
    gsl_function f;
    f.function = &__numdiff_call_fun__<Instance>;
    f.params   = this;

    // calculate derivative
    double result;
    if      (Method==-1) gsl_deriv_backward (&f, AtX, h, &result, &LastAbsErr);
    else if (Method==0)  gsl_deriv_central  (&f, AtX, h, &result, &LastAbsErr);
    else if (Method==1)  gsl_deriv_forward  (&f, AtX, h, &result, &LastAbsErr);
    else throw new Fatal("Numerical::Diff: Method==%d is invalid. Valid ones are -1=backward, 0=central, 1=forward",Method);
    return result;
}

}; // namespace Numerical

#endif // MECHSYS_NUMDIFF_H
