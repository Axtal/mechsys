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

#ifndef MECHSYS_ODE_H
#define MECHSYS_ODE_H

// STL
#include <cmath>   // for sqrt
#include <cfloat>  // for DBL_EPSILON
#include <cstring> // for strcmp
#include <sstream> // for ostringstream

// GSL
#include <gsl/gsl_errno.h> // for GSL_SUCCESS
#include <gsl/gsl_odeiv.h>

// MechSys
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>

namespace Numerical
{

template<typename Instance>
class ODESolver
{
public:
    // Typedefs
    typedef int (Instance::*pFun) (double t, double const Y[], double dYdt[]); ///< Callback function
    typedef void (Instance::*pUpFun) (double t, double Y[]); ///< update function (for successful steps)

    // Constructor
    ODESolver (Instance * p2Inst, pFun p2Fun, size_t NEq, char const * Scheme="RKDP89",
               double stol=1.0e-6, double h=1.0, double EPSREL=DBL_EPSILON);

    // Destructor
    ~ODESolver ();

    // Methods
    void Evolve (double tf); ///< Evolve from t to tf
    void Evolve (double tf, double dtOut, char const * Filename);

    // Data
    double   t;      ///< Current time
    double * Y;      ///< Current vector of state variables
    pUpFun   UpFun;  ///< Update function
    String   Scheme; ///< Scheme name

    // ME data
    double STOL;
    size_t MaxSS;
    double mMin;
    double mMax;
    double T;
    double dTini;
    double dT;
    size_t SS;     ///< number of sub-steps
    size_t SSs;    ///< number of successful sub-steps

    // Internal methods
    int _call_fun (double time, double const y[], double dydt[]) { return (_p2inst->*_p2fun)(time, y, dydt); }

private:
    // Auxiliary variables
    size_t              _neq;
    Instance          * _p2inst;
    pFun                _p2fun;
    double              _h;
    gsl_odeiv_step    * _s;
    gsl_odeiv_control * _c;
    gsl_odeiv_evolve  * _e;
    gsl_odeiv_system    _sys;

    // for FE/ME
    bool     _FE; ///< Forward-Euler ?
    bool     _ME; ///< Modified-Euler (RK12) ?
    double * yFE;
    double * yME;
    double * dy1;
    double * dy2;
};

/** Trick to pass pointers to member functions to GSL.
 * This will use (*params) to store the pointer to an instance of ODESolver, therefore (*params) will no longer be available. */
template<typename Instance>
int __ode_call_fun__ (double time, double const y[], double dydt[], void * not_used_for_params)
{
    ODESolver<Instance> * p2inst = static_cast<ODESolver<Instance>*>(not_used_for_params);
    return p2inst->_call_fun (time, y, dydt);
}


#ifdef USE_BOOST_PYTHON

class PyODESolver
{
public:
    /** Constructor. */
    PyODESolver (BPy::object & TheInstance, char const * ClassName, char const * MethodName)
        : Solver(NULL), Instance(TheInstance)
    {
        Util::GetPyMethod (ClassName, MethodName, Method, "__main__");
    }

    /** Destructor. */
    ~PyODESolver () { if (Solver!=NULL) delete Solver; }

    /** Initialization method. */
    void Init (double t0, BPy::list const & Y0, char const * Scheme="RKDP89", double stol=1.0e-6, double h=1.0, double EPSREL=DBL_EPSILON)
    {
        NEq = BPy::len (Y0);
        if (Solver!=NULL) delete Solver;
        Solver = new ODESolver<PyODESolver> (this, &PyODESolver::Fun, NEq, Scheme, stol, h, EPSREL);
        Solver->t = t0;
        for (size_t i=0; i<NEq; ++i)
        {
            double val = BPy::extract<double>(Y0[i])();
            y   .append (val);
            dydt.append (0.0);
            Solver->Y[i] = val;
        }
    }

    /** Evolve to tf. */
    void Evolve (double tf)
    {
        if (Solver==NULL) throw new Fatal("PyODESolver::Evolve: Init must be called first");
        Solver->Evolve (tf);
        for (size_t i=0; i<NEq; ++i) y[i] = Solver->Y[i];
    }

    /** Evolve to tf and output each dtOut. */
    void EvolveOut (double tf, double dtOut, char const * Filename)
    {
        if (Solver==NULL) throw new Fatal("PyODESolver::EvolveOut: Init must be called first");
        Solver->Evolve (tf, dtOut, Filename);
        for (size_t i=0; i<NEq; ++i) y[i] = Solver->Y[i];
    }

    /** Internal callback function. */
    int Fun (double t, double const Y[], double dYdt[])
    {
        for (size_t i=0; i<NEq; ++i) y[i] = Y[i];
        Method (Instance, t,y,dydt);
        for (size_t i=0; i<NEq; ++i) dYdt[i] = BPy::extract<double>(dydt[i])();
        return GSL_SUCCESS;
    }

    // Data
    size_t                   NEq;      ///< Number of equations
    ODESolver<PyODESolver> * Solver;   ///< Pointer to ODE solver
    BPy::object              Instance; ///< Pointer to Python object
    BPy::object              Method;   ///< Python method to be called in callback function
    BPy::list                y;        ///< Array of unknowns
    BPy::list                dydt;     ///< Derivative of y w.r.t time
};

#endif


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


template<typename Instance>
inline ODESolver<Instance>::ODESolver (Instance * p2Inst, pFun p2Fun, size_t NEq, char const * sch, double stol, double h, double EPSREL)
    : t       (-1.0),
      Y       (NULL),
      UpFun   (NULL),
      Scheme  (sch),
      STOL    (stol),
      MaxSS   (2000),
      mMin    (0.2),
      mMax    (5.0),
      T       (0.0),
      dTini   (h),
      dT      (h),
      SS      (0),
      SSs     (0),
      _neq    (NEq),
      _p2inst (p2Inst),
      _p2fun  (p2Fun),
      _h      (h),
      _s      (NULL),
      _c      (NULL),
      _e      (NULL),
      _FE     (false),
      _ME     (false),
      yFE     (NULL),
      yME     (NULL),
      dy1     (NULL),
      dy2     (NULL)
{
    if (Scheme=="FE")
    {
        _FE = true;
        dy1 = new double [_neq];
    }
    else if (Scheme=="ME" || Scheme=="RK12")
    {
        _ME = true;
        yFE = new double [_neq];
        yME = new double [_neq];
        dy1 = new double [_neq];
        dy2 = new double [_neq];
    }
    else
    {
        // set scheme
        gsl_odeiv_step_type const * scheme = gsl_odeiv_step_rk8pd;
        if      (Scheme=="RK23"  ) { scheme = gsl_odeiv_step_rk2;    }
        else if (Scheme=="RK4"   ) { scheme = gsl_odeiv_step_rk4;    }
        else if (Scheme=="RKF45" ) { scheme = gsl_odeiv_step_rkf45;  }
        else if (Scheme=="RKCK45") { scheme = gsl_odeiv_step_rkck;   }
        else if (Scheme=="RKDP89") { scheme = gsl_odeiv_step_rk8pd;  }
        else if (Scheme=="RK2I"  ) { scheme = gsl_odeiv_step_rk2imp; }
        else if (Scheme=="RK4I"  ) { scheme = gsl_odeiv_step_rk4imp; }
        else if (Scheme=="G1"    ) { scheme = gsl_odeiv_step_gear1;  }
        else if (Scheme=="G2"    ) { scheme = gsl_odeiv_step_gear2;  }
        else throw new Fatal("ODESolver::ODESolver: Scheme %s is invalid. Valid ones are: FE, ME/RK12, RK23, RK4, RKF45, RKCK45, RKDP89, RK2I, RK4I, G1, G2",Scheme.CStr());

        // allocate auxiliary variables
        _s   = gsl_odeiv_step_alloc    (scheme, _neq);
        _c   = gsl_odeiv_control_y_new (stol, EPSREL);
        _e   = gsl_odeiv_evolve_alloc  (_neq);
        _sys.function  = __ode_call_fun__<Instance>;
        _sys.jacobian  = NULL;
        _sys.dimension = _neq;
        _sys.params    = this;
    }

    // allocate array
    Y = new double [_neq];
}

template<typename Instance>
inline ODESolver<Instance>::~ODESolver ()
{
    if (Y  !=NULL) delete [] Y;
    if (yFE!=NULL) delete [] yFE;
    if (yME!=NULL) delete [] yME;
    if (dy1!=NULL) delete [] dy1;
    if (dy2!=NULL) delete [] dy2;
    if (_e !=NULL) gsl_odeiv_evolve_free  (_e);
    if (_c !=NULL) gsl_odeiv_control_free (_c);
    if (_s !=NULL) gsl_odeiv_step_free    (_s);
}

template<typename Instance>
inline void ODESolver<Instance>::Evolve (double tf)
{
    // check
    if (t<0.0) throw new Fatal("ODESolver::Evolve: Initial values (t, Y) must be set first");

    // solve
    if (_FE)
    {
        while (t<tf)
        {
            double dt = (t+_h>tf ? tf-t : _h);
            (_p2inst->*_p2fun) (t, Y, dy1);
            for (size_t j=0; j<_neq; ++j) Y[j] += dy1[j]*dt;
            t += dt;
        }
    }
    else if (_ME)
    {
        // for each pseudo time T
        double Dtf = tf - t;
        T   = 0.0;
        dT  = dTini;
        SS  = 0;
        SSs = 0;
        for (SS=0; SS<MaxSS; ++SS)
        {
            // exit point
            if (T>=1.0) break;

            // FE increments
            double dt = dT*Dtf;
            (_p2inst->*_p2fun) (t, Y, dy1);
            for (size_t j=0; j<_neq; ++j) yFE[j] = Y[j] + dy1[j]*dt;

            // ME increments and local error
            (_p2inst->*_p2fun) (t+dt, yFE, dy2);
            double norm_yME = 0.0;
            double norm_err = 0.0;
            for (size_t j=0; j<_neq; ++j)
            {
                yME[j]    = Y[j] + 0.5*(dy1[j]+dy2[j])*dt;
                norm_yME += yME[j]*yME[j];
                norm_err += (yME[j]-yFE[j])*(yME[j]-yFE[j]);
            }

            // local error estimate
            double err = sqrt(norm_err)/(1.0+sqrt(norm_yME));

            // step multiplier
            double m = (err>0.0 ? 0.9*sqrt(STOL/err) : mMax);

            // update
            if (err<STOL)
            {
                // update state
                SSs++;
                T += dT;
                t += dt;
                for (size_t j=0; j<_neq; ++j) Y[j] = yME[j];

                // callback
                if (UpFun!=NULL) (_p2inst->*UpFun) (t, Y);

                // limit change on stepsize
                if (m>mMax) m = mMax;
            }
            else if (m<mMin) m = mMin;

            // change next step size
            dT = m * dT;

            // check for last increment
            if (dT>1.0-T) dT = 1.0-T;
        }
        if (SS>=MaxSS) throw new Fatal("ODESolver::Evolve: RK12 (Modified-Euler) did not converge after %d substeps",SS);
    }
    else
    {
        while (t<tf)
        {
           int status = gsl_odeiv_evolve_apply (_e, _c, _s, &_sys, &t, tf, &_h, Y);
           if (status!=GSL_SUCCESS) throw new Fatal("ODESolver::Evolve: gsl_odeiv_evolve_apply failed. status=%d",status);
           //if (h<hmin) { h = hmin; }
        }
    }
}

template<typename Instance>
inline void ODESolver<Instance>::Evolve (double tf, double dtOut, char const * Filename)
{
    // header and first output
    String buf;
    std::ostringstream oss;
    oss<<Util::_8s<<"t";  for (size_t i=0; i<_neq; ++i) { buf.Printf("y%zd",i);  oss<<Util::_8s<<buf;  }  oss<<std::endl;
    oss<<Util::_8s<<t;    for (size_t i=0; i<_neq; ++i) {                        oss<<Util::_8s<<Y[i]; }  oss<<std::endl;

    // solve
    double tout = t + dtOut;
    while (t<tf)
    {
        Evolve (tout);
        tout += dtOut;
        if (tout>tf) tout = tf;
        oss<<Util::_8s<<t;  for (size_t i=0; i<_neq; ++i) { oss<<Util::_8s<<Y[i]; }  oss<<std::endl;
    }

    // file
    std::ofstream of(Filename, std::ios::out);
    of << oss.str();
    of.close();
    printf("File <%s%s%s> written\n", TERM_CLR_BLUE_H, Filename, TERM_RST);
}

}; // namespace Numerical

#endif // MECHSYS_ODE_H
