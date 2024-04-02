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

#ifndef MECHSYS_NUMERICAL_MIN_H
#define MECHSYS_NUMERICAL_MIN_H

// STL
#include <cmath>   // for sqrt
#include <cfloat>  // for DBL_EPSILON
#include <cstring> // for strcmp
#include <sstream> // for ostringstream

// GSL
#include <gsl/gsl_errno.h> // for GSL_SUCCESS
#include <gsl/gsl_multimin.h>

// MechSys
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>

namespace Numerical
{

template<typename Instance>
class Min
{
public:
    // Typedefs
    typedef double (Instance::*pF)  (double const X[]);                ///< calc F (returns F(X))
    typedef void   (Instance::*pdF) (double const X[], double dFdX[]); ///< calc dFdX

    // Constructor
    Min (Instance * p2Inst, pF p2F, pdF p2dF, size_t NDim, char const * Scheme="BFGS"); ///< p2dF can be NULL, in this case Scheme must be NMSimplex or NMSimplexRnd

    // Destructor
    ~Min ();

    // Methods
    void   SetScheme (char const * Scheme="BFGS"); ///< if p2dF is NULL, scheme must be NMSimplex or NMSimplexRnd
    double Find      (double X[], double StepSize=0.01, double Tol=1.0e-3, size_t MaxIt=1000); ///< Find min around X. Returns results in X

    // Data
    bool Debug;

    // Internal methods
    double _call_F  (double const X[])                { return (_p2inst->*_p2F) (X); }
    void   _call_dF (double const X[], double dFdX[]) { if (_p2dF==NULL) throw new Fatal("Numerical::_call_dF: function for the gradients (p2dF) is not available"); (_p2inst->*_p2dF) (X, dFdX); }

private:
    // Auxiliary variables
    Instance                  * _p2inst;
    pF                          _p2F;
    pdF                         _p2dF;
    size_t                      _ndim;
    gsl_multimin_function_fdf   _gsl_fdf;
    gsl_multimin_function       _gsl_f;
    gsl_multimin_fdfminimizer * _gsl_fdfmin;
    gsl_multimin_fminimizer   * _gsl_fmin;
    gsl_vector                * _x;
    gsl_vector                * _ss;
};

// Trick to pass pointers to member functions to GSL. This will use (*params) to store the pointer to an instance of Min, therefore (*params) will no longer be available.
template<typename Instance>
double __min_call_F__ (gsl_vector const * X, void * not_used_for_params)
{
    Min<Instance> * p2min = static_cast<Min<Instance>*>(not_used_for_params);
    return p2min->_call_F (X->data);
}

// Trick to pass pointers to member functions to GSL. This will use (*params) to store the pointer to an instance of Min, therefore (*params) will no longer be available.
template<typename Instance>
void __min_call_dF__ (gsl_vector const * X, void * not_used_for_params, gsl_vector * dFdX)
{
    Min<Instance> * p2min = static_cast<Min<Instance>*>(not_used_for_params);
    p2min->_call_dF (X->data, dFdX->data);
}

// Trick to pass pointers to member functions to GSL. This will use (*params) to store the pointer to an instance of Min, therefore (*params) will no longer be available.
template<typename Instance>
void __min_call_FdF__ (gsl_vector const * X, void * not_used_for_params, double * F, gsl_vector * dFdX)
{
    Min<Instance> * p2min = static_cast<Min<Instance>*>(not_used_for_params);
    (*F) = p2min->_call_F (X->data);
    p2min->_call_dF (X->data, dFdX->data);
}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Instance>
inline Min<Instance>::Min (Instance * p2Inst, pF p2F, pdF p2dF, size_t NDim, char const * Scheme)
    : Debug       (false),
      _p2inst     (p2Inst),
      _p2F        (p2F),
      _p2dF       (p2dF),
      _ndim       (NDim),
      _gsl_fdfmin (NULL),
      _gsl_fmin   (NULL),
      _ss         (NULL)
{
    // functions (with gradients)
    _gsl_fdf.n      = _ndim;
    _gsl_fdf.f      = &__min_call_F__<Instance>;
    _gsl_fdf.df     = &__min_call_dF__<Instance>;
    _gsl_fdf.fdf    = &__min_call_FdF__<Instance>;
    _gsl_fdf.params = this;

    // functions
    _gsl_f.n      = _ndim;
    _gsl_f.f      = &__min_call_F__<Instance>;
    _gsl_f.params = this;

    // allocate vector
    _x = gsl_vector_alloc (_ndim);

    // set scheme
    SetScheme (Scheme);
}

template<typename Instance>
inline Min<Instance>::~Min ()
{
    if (_gsl_fdfmin!=NULL) gsl_multimin_fdfminimizer_free (_gsl_fdfmin);
    if (_gsl_fmin  !=NULL) gsl_multimin_fminimizer_free   (_gsl_fmin);
    if (_x         !=NULL) gsl_vector_free                (_x);
    if (_ss        !=NULL) gsl_vector_free                (_ss);
}

template<typename Instance>
inline void Min<Instance>::SetScheme (char const * Scheme)
{
    // set scheme
    gsl_multimin_fdfminimizer_type const * fdftype = NULL;
    gsl_multimin_fminimizer_type   const * ftype   = NULL;
    if      (strcmp(Scheme,"ConjFR"      )==0) { fdftype = gsl_multimin_fdfminimizer_conjugate_fr;     }
    else if (strcmp(Scheme,"ConjPR"      )==0) { fdftype = gsl_multimin_fdfminimizer_conjugate_pr;     }
    else if (strcmp(Scheme,"BFGS"        )==0) { fdftype = gsl_multimin_fdfminimizer_vector_bfgs2;     }
    else if (strcmp(Scheme,"Steep"       )==0) { fdftype = gsl_multimin_fdfminimizer_steepest_descent; }
    else if (strcmp(Scheme,"NMSimplex"   )==0) { ftype   = gsl_multimin_fminimizer_nmsimplex2;         }
    else if (strcmp(Scheme,"NMSimplexRnd")==0) { ftype   = gsl_multimin_fminimizer_nmsimplex2rand;     }
    else throw new Fatal("Numerical::Min::SetScheme: Scheme %s is invalid. Valid ones are: ConjFR, ConjPR, BFGS, Steep, NMSimplex, NMSimplexRnd",Scheme);

    // clear previous minimizers
    if (_gsl_fdfmin!=NULL) { gsl_multimin_fdfminimizer_free (_gsl_fdfmin);  _gsl_fdfmin = NULL; }
    if (_gsl_fmin  !=NULL) { gsl_multimin_fminimizer_free   (_gsl_fmin);    _gsl_fmin   = NULL; }

    // set minimizer
    if (ftype==NULL) // with gradients
    {
        _gsl_fdfmin = gsl_multimin_fdfminimizer_alloc (fdftype, _ndim);
    }
    else
    {
        _gsl_fmin = gsl_multimin_fminimizer_alloc (ftype, _ndim);
        if (_ss==NULL) _ss = gsl_vector_alloc (_ndim);
    }
}

template<typename Instance>
inline double Min<Instance>::Find (double X[], double StepSize, double Tol, size_t MaxIt)
{
    // set initial X
    for (size_t i=0; i<_ndim; ++i) gsl_vector_set (_x, i, X[i]);

   // with gradients
    if (_gsl_fmin==NULL)
    {
        // set minimizer
        gsl_multimin_fdfminimizer_set (_gsl_fdfmin, &_gsl_fdf, _x, StepSize, Tol/10.0);

        // solve
        int    status=0;
        size_t it;
        for (it=0; it<MaxIt; ++it)
        {
            status = gsl_multimin_fdfminimizer_iterate (_gsl_fdfmin);
            if (status==GSL_ENOPROG) break;
            status = gsl_multimin_test_gradient (_gsl_fdfmin->gradient, Tol);
            if (status==GSL_SUCCESS)
            {
                for (size_t i=0; i<_ndim; ++i) X[i] = gsl_vector_get (_gsl_fdfmin->x, i);
                return _gsl_fdfmin->f;
            }
            if (Debug)
            {
                std::cout << Util::_6 << it;
                for (size_t i=0; i<_ndim; ++i) std::cout << Util::_8s << gsl_vector_get(_gsl_fdfmin->x, i);
                std::cout << TERM_CLR_MAGENTA_H << Util::_8s << _gsl_fdfmin->f << TERM_RST << std::endl; 
            }
        }

        // error
        throw new Fatal("Numerical::Min::Find: (With gradients) Could not find minimum (it=%zd, MaxIt=%zd, GSL_STATUS=%d)",it,MaxIt,status);
    }

    // no gradients
    else
    {
        // set initial step sizes
        gsl_vector_set_all (_ss, StepSize);

        // set minimizer
        gsl_multimin_fminimizer_set (_gsl_fmin, &_gsl_f, _x, _ss);

        // solve
        int    status=0;
        size_t it;
        double size;
        for (it=0; it<MaxIt; ++it)
        {
            status = gsl_multimin_fminimizer_iterate (_gsl_fmin);
            if (status==GSL_ENOPROG) break;
            size   = gsl_multimin_fminimizer_size (_gsl_fmin);
            status = gsl_multimin_test_size       (size, Tol);
            if (status==GSL_SUCCESS)
            {
                for (size_t i=0; i<_ndim; ++i) X[i] = gsl_vector_get (_gsl_fmin->x, i);
                return _gsl_fmin->fval;
            }
            if (Debug)
            {
                std::cout << Util::_6 << it;
                for (size_t i=0; i<_ndim; ++i) std::cout << Util::_8s << gsl_vector_get(_gsl_fmin->x, i);
                std::cout << TERM_CLR_MAGENTA_H << Util::_8s << _gsl_fmin->fval << TERM_RST << std::endl; 
            }
        }

        // error
        throw new Fatal("Numerical::Min::Find: (no gradients) Could not find minimum (it=%zd, MaxIt=%zd, GSL_STATUS=%d)",it,MaxIt,status);
    }
}

}; // namespace Numerical

#endif // MECHSYS_NUMERICAL_MIN_H
