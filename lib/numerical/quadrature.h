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

#ifndef MECHSYS_QUADRATURE_H
#define MECHSYS_QUADRATURE_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// GSL
#include <gsl/gsl_integration.h>

namespace Numerical
{

/** Integration (quadrature) method. */
enum QMethod
{
	QNG_T, ///< QNG  : non-adaptive Gauss-Kronrod integration
	QAGS_T ///< QAGS : adaptive integration with singularities
};

/** %Numerical integration.
Examples:
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/numerical/tst/tquadrature.cpp?view=markup"> ttquadrature.cpp Test numerical integration</a>
*/
template<typename Instance>
class Quadrature
{
public:
	// Constants
	static int WORKSPACE_SIZE; ///< Workspace size

	// Typedefs
	typedef double (Instance::*pFun) (double x) const; ///< Callback function

	/** Constructor. */
	 Quadrature (Instance const * p2Inst, pFun p2Fun, QMethod method=QAGS_T,double EPSREL=100.0*DBL_EPSILON); // p2Fun == &F(x)

	/** Destructor. */
	~Quadrature ();

	// Methods
	double Integrate (double a, double b);                          ///< Integrate F(x) in [a,b]
	double CallFun   (double x)   { return (_p2inst->*_p2fun)(x); } ///< Call F(x)
	double LastError () const     { return _last_error;           } ///< Return the last error after an Integrate call

private:
	Instance const            * _p2inst;     ///< Pointer to an instance
	pFun                        _p2fun;      ///< Pointer to F(x) function
	QMethod                     _method;     ///< Method of quadrature
	gsl_integration_workspace * _workspace;  ///< Workspace
	double                      _epsabs;     ///< Absolute error
	double                      _epsrel;     ///< Relative error
	double                      _last_error; ///< Last error after an Integrate call
	size_t                      _last_neval; ///< Last number of function evaluations used

}; // class Quadrature

template<typename Instance>
int Quadrature<Instance>::WORKSPACE_SIZE = 10000;

/** Trick to pass pointers to member functions to GSL.
 * This will use (*params) to store the pointer to an instance of Quadrature, 
 * therefore (*params) will no longer be available.
 */
template<typename Instance>
double __quadrature_call_fun__(double x, void *not_used_for_params)
{
	Quadrature<Instance> * p2inst = static_cast<Quadrature<Instance>*>(not_used_for_params);
	return p2inst->CallFun(x);

}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Instance>
inline Quadrature<Instance>::Quadrature(Instance const * p2Inst, pFun p2Fun, QMethod method, double EPSREL)
	: _p2inst     (p2Inst),
	  _p2fun      (p2Fun),
	  _method     (method),
	  _epsabs     (EPSREL),
	  _epsrel     (EPSREL),
	  _last_error (-1),
	  _last_neval (-1)
{
	/* Each algorithm computes an approximation to a definite integral of the form,
	 * 
	 *      I = \int_a^b f(x) w(x) dx
	 * 
	 * where w(x) is a weight function (for general integrands w(x)=1).
	 *
	 * The algorithms are built on pairs of quadrature rules, a higher order
	 * rule and a lower order rule.  The higher order rule is used to compute
	 * the best approximation to an integral over a small range.  The
	 * difference between the results of the higher order rule and the lower
	 * order rule gives an estimate of the error in the approximation.
	 *
	 * The user provides absolute and relative error bounds (epsabs, epsrel) which
	 * specify the following accuracy requirement
	 * 
	 *      |RESULT - I|  <= max(epsabs, epsrel |I|)
	 * 
	 * where RESULT is the numerical approximation obtained by the algorithm.
	 *
	 * The algorithms attempt to estimate the absolute error 
	 * ABSERR = |RESULT - I| in such a way that the following inequality holds,
	 * 
	 *      |RESULT - I| <= ABSERR <= max(epsabs, epsrel |I|)
	 * 
	 * The routines will fail to converge if the error bounds are too
	 * stringent, but always return the best approximation obtained up to that
	 * stage.
	 * 
	 *    The algorithms in QUADPACK use a naming convention based on the
	 * following letters,
	 * 
	 *      `Q' - quadrature routine
	 * 
	 *      `N' - non-adaptive integrator
	 *      `A' - adaptive integrator
	 * 
	 *      `G' - general integrand (user-defined)
	 *      `W' - weight function with integrand
	 * 
	 *      `S' - singularities can be more readily integrated
	 *      `P' - points of special difficulty can be supplied
	 *      `I' - infinite range of integration
	 *      `O' - oscillatory weight function, cos or sin
	 *      `F' - Fourier integral
	 *      `C' - Cauchy principal value
	 */

	_workspace = gsl_integration_workspace_alloc (WORKSPACE_SIZE);
}

template<typename Instance>
inline Quadrature<Instance>::~Quadrature()
{
	if (_workspace!=NULL) gsl_integration_workspace_free (_workspace);
}

template<typename Instance>
inline double Quadrature<Instance>::Integrate(double a, double b)
{
	// Function
	gsl_function f;
	f.function = &__quadrature_call_fun__<Instance>;
	f.params   = this;

	// Solve
	double result;
	switch (_method)
	{
		case QNG_T:
		{
			/* QNG non-adaptive Gauss-Kronrod integration
			 *
			 * The QNG algorithm is a non-adaptive procedure which uses fixed Gauss-Kronrod abscissae to sample the integrand at a
			 * maximum of 87 points. It is provided for fast integration of smooth functions.
			 *
			 * Function: int gsl_integration_qng (const gsl_function * F, double A, double B, double EPSABS, double EPSREL,
			 * double * RESULT, double * ABSERR, size_t * NEVAL)
			 *
			 * This function applies the Gauss-Kronrod 10-point, 21-point, 43-point and 87-point integration rules in succession
			 * until an estimate of the integral of f over (a,b) is achieved within the desired absolute and relative error limits,
			 * EPSABS and EPSREL. The function returns the final approximation, RESULT, an estimate of the absolute error, ABSERR
			 * and the number of function evaluations used, NEVAL.  The Gauss-Kronrod rules are designed in such a way that each
			 * rule uses all the results of its predecessors, in order to minimize the total number of function evaluations.
			 */

			gsl_integration_qng (&f, a, b, _epsabs, _epsrel, &result, &_last_error, &_last_neval);

			break;
		}
		case QAGS_T:
		{
			/* QAGS adaptive integration with singularities
			 *
			 * The presence of an integrable singularity in the integration region causes an adaptive routine to concentrate new
			 * subintervals around the singularity. As the subintervals decrease in size the successive approximations to the integral
			 * converge in a limiting fashion. This approach to the limit can be accelerated using an extrapolation procedure. The
			 * QAGS algorithm combines adaptive bisection with the Wynn epsilon-algorithm to speed up the integration of many types of
			 * integrable singularities.
			 *
			 * Function: int gsl_integration_qags (const gsl_function * f, double a, double b, double epsabs, double epsrel, size_t
			 * limit, gsl_integration_workspace * workspace, double * result, double * abserr)
			 *
			 * This function applies the Gauss-Kronrod 21-point integration rule adaptively until an estimate of the integral of f
			 * over (a,b) is achieved within the desired absolute and relative error limits, epsabs and epsrel. The results are
			 * extrapolated using the epsilon-algorithm, which accelerates the convergence of the integral in the presence of
			 * discontinuities and integrable singularities. The function returns the final approximation from the extrapolation,
			 * result, and an estimate of the absolute error, abserr. The subintervals and their results are stored in the memory
			 * provided by workspace. The maximum number of subintervals is given by limit, which may not exceed the allocated size of
			 * the workspace. 
			 */

			gsl_integration_qags (&f, a, b, _epsabs, _epsrel, WORKSPACE_SIZE, _workspace, &result, &_last_error); 

			break;
		}
	};

	// Results
	return result;
}

}; // namespace Numerical

#endif // MECHSYS_QUADRATURE_H
