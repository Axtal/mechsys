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

#ifndef MECHSYS_UMFPACK_H
#define MECHSYS_UMFPACK_H

// Std Lib
#include <cmath> // for pow

// UMFPACK
extern "C" {
#include <umfpack.h>
}

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/sparse_triplet.h>
#include <mechsys/linalg/sparse_matrix.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

#ifdef DO_DEBUG
  using std::cout;
  using std::endl;
#endif

/** \namespace UMFPACK Unsymmetric multifrontal sparse LU factorization package.
  %UMFPACK is a set of routines for solving unsymmetric sparse linear systems, Ax=b, using the Unsymmetric MultiFrontal method.
  See <a href="http://www.cise.ufl.edu/research/sparse/umfpack/">Tim Davis' web page.</a>
  
  Examples:
   - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tumfpack.cpp?view=markup">tumfpack.cpp Test UMFPACK solver</a>
 */
namespace UMFPACK
{

/** %UMFPACK internal error message.
 * \param info %UMFPACK internal error code
 */
String ErrorMsg(int info)
{
    switch (info)
    {
        case UMFPACK_ERROR_out_of_memory           : return String("UMFPACK_ERROR_out_of_memory (-1)");
        case UMFPACK_ERROR_invalid_Numeric_object  : return String("UMFPACK_ERROR_invalid_Numeric_object (-3)");
        case UMFPACK_ERROR_invalid_Symbolic_object : return String("UMFPACK_ERROR_invalid_Symbolic_object (-4)");
        case UMFPACK_ERROR_argument_missing        : return String("UMFPACK_ERROR_argument_missing (-5)");
        case UMFPACK_ERROR_n_nonpositive           : return String("UMFPACK_ERROR_n_nonpositive (-6)");
        case UMFPACK_ERROR_invalid_matrix          : return String("UMFPACK_ERROR_invalid_matrix (-8)");
        case UMFPACK_ERROR_different_pattern       : return String("UMFPACK_ERROR_different_pattern (-11)");
        case UMFPACK_ERROR_invalid_system          : return String("UMFPACK_ERROR_invalid_system (-13)");
        case UMFPACK_ERROR_invalid_permutation     : return String("UMFPACK_ERROR_invalid_permutation (-15)");
        case UMFPACK_ERROR_internal_error          : return String("UMFPACK_ERROR_internal_error (-911)");
        case UMFPACK_ERROR_file_IO                 : return String("UMFPACK_ERROR_file_IO (-17)");
        default                                    : { String r; r.Printf("Unknown UMFPACK_ERROR (%d)",info); return r; };
    }
}

/** Solves \f$ {X} \leftarrow [A]^{-1}{B} \f$.
 * \param A %Sparse matrix
 * \param B Right-hand side vector
 * \param X Result
 */
inline void Solve (Sparse::Matrix<double,int> const & A, Vec_t const & B, Vec_t & X)
{
    if (A.Rows()!=A.Cols())        throw new Fatal("UMFPACK::Solve: A (%d x %d) matrix must be squared.",A.Rows(),A.Cols());
    if (size(B)!=(size_t)A.Cols()) throw new Fatal("UMFPACK::Solve: B (%d x 1) vector must have a size equal to the number of columns of matrix A (%d)",size(B),A.Cols());
    if (size(X)!=size(B))          throw new Fatal("UMFPACK::Solve: X (%d x 1) vector must have the same size as B (%d)",size(X),size(B));
    double *null = (double *)NULL;
    void *symbolic, *numeric;
    int info = 0;
    info = umfpack_di_symbolic      (A.Rows(), A.Rows(), A.GetApPtr(), A.GetAiPtr(), A.GetAxPtr(), &symbolic, null, null);      if (info<0) throw new Fatal("UMFPACK::Solve: umfpack_di_symbolic failed. %s",ErrorMsg(info).CStr());
    info = umfpack_di_numeric       (A.GetApPtr(), A.GetAiPtr(), A.GetAxPtr(), symbolic, &numeric, null, null);                 if (info<0) throw new Fatal("UMFPACK::Solve: umfpack_di_numeric failed. %s",ErrorMsg(info).CStr());
           umfpack_di_free_symbolic (&symbolic);
    info = umfpack_di_solve         (UMFPACK_A, A.GetApPtr(), A.GetAiPtr(), A.GetAxPtr(), X.data, B.data, numeric, null, null); if (info<0) throw new Fatal("UMFPACK::Solve: umfpack_di_solve failed. %s",ErrorMsg(info).CStr());
           umfpack_di_free_numeric  (&numeric);
}
// TODO: Check: It seems that the following is not necessary
//       (A sparse matrix is automatically created from the Triplet)
inline void Solve (Sparse::Triplet<double,int> const & A, Vec_t const & B, Vec_t & X) { Solve(Sparse::Matrix<double,int>(A), B, X); }

inline double Det (Sparse::Matrix<double,int> const & A)
{
    double *null = (double *)NULL;
    void   *symbolic, *numeric;
    int    info = 0;
    double Mx, Ex;
    info = umfpack_di_symbolic        (A.Rows(), A.Rows(), A.GetApPtr(), A.GetAiPtr(), A.GetAxPtr(), &symbolic, null, null);  if (info<0) throw new Fatal("UMFPACK::Det: umfpack_di_symbolic failed. %s",ErrorMsg(info).CStr());
    info = umfpack_di_numeric         (A.GetApPtr(), A.GetAiPtr(), A.GetAxPtr(), symbolic, &numeric, null, null);             if (info<0) throw new Fatal("UMFPACK::Det: umfpack_di_numeric failed. %s",ErrorMsg(info).CStr());
           umfpack_di_free_symbolic   (&symbolic);
    info = umfpack_di_get_determinant (&Mx, &Ex, numeric, null);                                                              if (info<0) throw new Fatal("UMFPACK::Det: umfpack_di_numeric failed. %s",ErrorMsg(info).CStr());
    return Mx * pow (10.0, Ex);
}

class Sys
{
public:
    Sys (Sparse::Triplet<double,int> const & Atri)
    {
        if (Atri.Rows()!=Atri.Cols()) throw new Fatal("UMFPACK::Sys: A (%d x %d) matrix (triplet) must be square",Atri.Rows(),Atri.Cols());
        _A = new Sparse::Matrix<double,int> (Atri);
        double * null = (double*)NULL;
        int info = 0;
        info = umfpack_di_symbolic (_A->Rows(), _A->Rows(), _A->GetApPtr(), _A->GetAiPtr(), _A->GetAxPtr(), &_symbolic, null, null);  if (info<0) throw new Fatal("UMFPACK::Sys: umfpack_di_symbolic failed. %s",ErrorMsg(info).CStr());
        info = umfpack_di_numeric  (_A->GetApPtr(), _A->GetAiPtr(), _A->GetAxPtr(), _symbolic, &_numeric, null, null);                if (info<0) throw new Fatal("UMFPACK::Sys: umfpack_di_numeric failed. %s",ErrorMsg(info).CStr());
    }
    ~Sys () { delete _A; umfpack_di_free_symbolic(&_symbolic);  umfpack_di_free_numeric(&_numeric); }
    void Solve (Vec_t const & B, Vec_t & X)
    {
        if (size(B)!=(size_t)_A->Cols()) throw new Fatal("UMFPACK::Sys::Solve: B (%d x 1) vector must have a size equal to the number of columns of matrix A (%d)",size(B),_A->Cols());
        if (size(X)!=size(B))            throw new Fatal("UMFPACK::Sys::Solve: X (%d x 1) vector must have the same size as B (%d)",size(X),size(B));
        double * null = (double*)NULL;
        int info = umfpack_di_solve (UMFPACK_A, _A->GetApPtr(), _A->GetAiPtr(), _A->GetAxPtr(), X.data, B.data, _numeric, null, null); if (info<0) throw new Fatal("UMFPACK::Sys::Solve: umfpack_di_solve failed. %s",ErrorMsg(info).CStr());
    }
private:
    Sparse::Matrix<double,int> * _A;
    void                       * _symbolic;
    void                       * _numeric;
};

}; // namespace UMFPACK

#endif // MECHSYS_UMFPACK_H
