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

#ifndef MECHSYS_MULTIPLEROOTS_H
#define MECHSYS_MULTIPLEROOTS_H

// STL
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// MechSys
#include <mechsys/numerical/brentroot.h>
#include <mechsys/numerical/quadrature.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>

namespace Numerical
{

/** Find all roots for F(x)=0 inside [A,B].
Examples:
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/numerical/tst/tmultipleroots.cpp?view=markup"> ttmultipleroots.cpp Test multipler roots solver</a>
*/
template<typename Instance>
class MultipleRoots
{
public:
	// Typedefs
	typedef double (Instance::*pFun) (double x) const; ///< Callback function

	/** Constructor. */
	MultipleRoots (Instance const * p2Inst, pFun F, pFun G, pFun H, double Tol=sqrt(DBL_EPSILON), double Gamma=0.001);

	// Methods
	void   Solve       (double A, double B, Array<double> & R, bool Single=false) const; ///< Solve for all roots (R) in [A,B]
	double NRoots      (double a, double b) const;                    ///< Computes the number of roots in [a,b]
	void   SetTol      (double Tol) { _tol = Tol;      }              ///< Error tolerance
	void   SetGam      (double Gam) { _gam = Gam;      }              ///< Small positive constant
	void   SetQuadError(double Err) { _qerr= Err;      }              ///< Quadrature Error tolerance
	int    Ncalls      () const     { return _N_calls; }              ///< Number of calls to NRoots function
	int    Mcalls      () const     { return _M_calls; }              ///< Number of calls to M(x) function
	void   ResetNcalls ()           { _N_calls = 0;    }              ///< Reset the number of calls to NRoots function
	void   ResetMcalls ()           { _M_calls = 0;    }              ///< Reset the number of calls to M(x) function

private:
	Instance const * _p2inst;  ///< Pointer to instance
	pFun             _F;       ///< Pointer to F(x) function
	pFun             _G;       ///< Pointer to G(x) = F'(x) function
	pFun             _H;       ///< Pointer to H(x) = F''(x) function
	double           _tol;     ///< Error tolerance
	double           _gam;     ///< Small positive constant
	double           _qerr;    ///< Quadrature Error tolerance
	mutable int      _N_calls; ///< Number of calls to NRoots function
	mutable int      _M_calls; ///< Number of calls to M function

   	// private functions
	double _M          (double x) const;                              ///< Auxiliar function used in order to compute the number of roots
	void   _find_roots (double a, double b, Array<double> & R) const; ///< Function for recursion during Solve

}; // class MultipleRoots


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

template<typename Instance>
inline MultipleRoots<Instance>::MultipleRoots(Instance const * p2Inst, pFun F, pFun G, pFun H, double Tol, double Gamma)
	: _p2inst  (p2Inst),
	  _F       (F),
	  _G       (G),
	  _H       (H),
	  _tol     (Tol),
	  _gam     (Gamma),
	  _qerr    (1.0e-3),
	  _N_calls (0),
	  _M_calls (0)
{
}

template<typename Instance>
inline void MultipleRoots<Instance>::Solve(double A, double B, Array<double> & R, bool Single) const
{
	// No roots
	R.Resize(0);

	// Single Root - using Newton's method
	if (Single)
	{
		double f0 = (_p2inst->*_F)(A);
		double f1 = (_p2inst->*_F)(B);
		if (f0*f1<0.0)
		{
			double ak = f0/(f0-f1);
			for (_N_calls=1; _N_calls<=30; ++_N_calls)
			{
				double fak = (_p2inst->*_F)(ak);
				if (fabs(fak)<_tol) { R.Push(ak); return; }
				ak += -fak/(_p2inst->*_G)(ak);
			}
			throw new Fatal(_("MultipleRoots::Solve: Did not converge for %d iterations"),30);
		}
	}
	else _find_roots(A,B,R); // Recursive method
}

template<typename Instance>
inline double MultipleRoots<Instance>::NRoots(double a, double b) const
{
	double f0 = (_p2inst->*_F)(a);
	double f1 = (_p2inst->*_F)(b);
	double g0 = (_p2inst->*_G)(a);
	double g1 = (_p2inst->*_G)(b);

	Quadrature<MultipleRoots> qu(this, &MultipleRoots::_M, /*method*/Numerical::QAGS_T, /*EPSREL*/0.001);
	double res1 = -_gam*qu.Integrate(a,b);
	double res2 = atan(_gam*(f0*g1-f1*g0)/(f0*f1+g0*g1*_gam*_gam+_tol));

	_N_calls++;

	return (res1+res2)/Util::Pi();
}


/* private */

template<typename Instance>
inline double MultipleRoots<Instance>::_M(double x) const
{
	double f = (_p2inst->*_F)(x);
	double g = (_p2inst->*_G)(x);
	double h = (_p2inst->*_H)(x);

	_M_calls++;

	double res = (f*h-g*g)/(f*f+_gam*_gam*g*g+_tol);
	//std::cout << "f=" << f << ", g=" << g << ", h=" << h << ", res=" << res << std::endl;
	if (res!=res) throw new Fatal (_("MultipleRoots::_M: integrand function is NaN"));
	return res;
}

template<typename Instance>
inline void MultipleRoots<Instance>::_find_roots(double a, double b, Array<double> & R) const
{
	double n  = NRoots(a,b);
	long   nr = lround(n);
	//std::cout << "a=" << a << ", b=" << b << ", nr=" << n << std::endl;
	if (nr==0) return;
	if (nr==1 && ((_p2inst->*_F)(a)*(_p2inst->*_F)(b)<0))
	{
		BrentRoot<Instance> br(_p2inst, _F);
		R.Push(br.Solve(a,b));
	}
	else
	{
		double mid = (a+b)/2.0;
		_find_roots (a        , mid , R);
		_find_roots (mid+_tol , b   , R);
	}
}


}; // namespace Numerical

#endif // MECHSYS_MULTIPLEROOTS_H
