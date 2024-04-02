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

#ifndef MECHSYS_MUMPS_H
#define MECHSYS_MUMPS_H

// MUMPS
#ifdef HAS_MPI
 extern "C" {
   #include "dmumps_c.h"
 }
#endif

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/sparse_triplet.h>
#include <mechsys/util/fatal.h>

namespace MUMPS
{

/** Solves \f$ {X} \leftarrow [A]^{-1}{B} \f$.
 * \param A DISTRIBUTED %Sparse matrix
 * \param B DISTRIBUTED Right-hand side vector
 * \param X Complete result vector in all processors
 *
 * THIS METHOD ALTERS MATRIX A
 *
 */
inline void Solve (Sparse::Triplet<double,int> & A, Vec_t const & B, Vec_t & X, bool Prod=false)
{
#ifdef HAS_MPI
    // collect B from all processors into X==RHS in processor # 0
    int neq = A.Rows();
    if (Prod)
    {
        MPI::COMM_WORLD.Reduce (B.data, X.data, neq, MPI::DOUBLE, MPI::PROD, /*dest*/0);
    }
    else
    {
        MPI::COMM_WORLD.Reduce (B.data, X.data, neq, MPI::DOUBLE, MPI::SUM, /*dest*/0);
    }

    // initialize MUMPS
    DMUMPS_STRUC_C ms;
    ms.comm_fortran = -987654;
    ms.sym          =  0;    // 0=unsymmetric, 1=sym(pos-def), 2=symmetric(undef)
    ms.par          =  1;    // host also works
    ms.job          = -1;    // initialisation code
    dmumps_c (&ms);          // do initialize

    // set matrix and rhs
    A.IncIndices (1); // increment indices since MUMS is 1-based (FORTRAN)
    ms.n       = neq;
    ms.nz_loc  = A.Top();
    ms.irn_loc = A.GetAiPtr();
    ms.jcn_loc = A.GetAjPtr();
    ms.a_loc   = A.GetAxPtr();
    if (MPI::COMM_WORLD.Get_rank()==0) ms.rhs = X.data; // only proc # 0 needs the RHS

    // solve
    ms.icntl[1  -1] = -1;    // no output messages
    ms.icntl[2  -1] = -1;    // no warnings
    ms.icntl[3  -1] = -1;    // no global information
    ms.icntl[4  -1] = -1;    // message level
    ms.icntl[5  -1] =  0;    // assembled matrix (needed for distributed matrix)
    ms.icntl[7  -1] =  5;    // use metis for pivoting
    ms.icntl[8  -1] =  0;    // no scaling
    ms.icntl[18 -1] =  3;    // distributed matrix
    ms.job          =  6;    // reorder, factor and solve codes
    dmumps_c (&ms);          // do solve
    int info = ms.info[1-1]; // info code

    // finalize MUMPS
    ms.job = -2;    // finalize code
    dmumps_c (&ms); // do finalize

    // check for errors
    if (info<0)
    {
        switch (info)
        {
            case -6:  throw new Fatal("MUMPS::Solve: Error # -6: singular matrix");
            case -10: throw new Fatal("MUMPS::Solve: Error # -10: singular matrix");
            case -13: throw new Fatal("MUMPS::Solve: Error # -13: out of memory");
            default:  throw new Fatal("MUMPS::Solve: Error # %d: unknown error",info);
        }
    }

    // distribute solution
    MPI::COMM_WORLD.Bcast (X.data, neq, MPI::DOUBLE, /*from*/0);

#else
    throw new Fatal("MUMPS::Solve: This method requires the flag HAS_MPI");
#endif
}

}; // namespace MUMPS

#endif // MECHSYS_MUMPS_H
