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

#ifndef MECHSYS_SUPERLUD_H
#define MECHSYS_SUPERLUD_H

// SuperLU_DIST
#include <superlu_ddefs.h>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/sparse_crmatrix.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

#ifdef DO_DEBUG
  #include <iostream>
  #include <sstream>
  using std::cout;
  using std::endl;
#endif

namespace SuperLUd
{

/** \namespace SuperLUd Nonsymmetric sparse solver (distributed/parallel).
  %SuperLU is a general purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations on high performance machines.
  See <a href="http://crd.lbl.gov/~xiaoye/SuperLU/">Xiaoye's web page.</a>
 */

/** Solves \f$ {B} \leftarrow [A]^{-1}{B} \f$.
 * \param MyMinEq Minimum equation ID of this processor
 * \param MyMaxEq Maximum equation ID of this processor
 * \param MyNumEqs Number of equations of this processor
 * \param nProcs Number of processors
 * \param nDOFs Number of degrees of freedom (total size of A matrix)
 * \param AugLocalA Local sparse matrix (compressed-row format)
 * \param GlobalB Global right-hand side vector and solution vector
 */
inline void Solve(int MyMinEq, int MyMaxEq, int MyNumEqs, int nProcs, int nDOFs, Sparse::CRMatrix<double,int> & AugLocalA, Vec_t & GlobalB)
{

	// Sizes
	int m_loc = MyNumEqs; // local number of rows (equations)
	int m     = nDOFs;    // global number of rows
	int n     = nDOFs;    // global number of cols (must be equal to m for SuperLU_DIST)

	// 1) Set SuperLU distributed matrix ------------------------------------------------
	
	// Local Aq array
	int * Aq = new int [m_loc+1];
	for (int i=0; i<m_loc+1; ++i) Aq[i] = AugLocalA.Aq(MyMinEq+i)-AugLocalA.Aq(MyMinEq);
	int nz_loc = Aq[m_loc]; // number of non-zeros in the local matrix

    // Create matrix A in the format expected by SuperLU_DIST
	int         start = AugLocalA.Aq(MyMinEq);
	double const * Ay = &AugLocalA.GetAyPtr()[start];
	int    const * Aj = &AugLocalA.GetAjPtr()[start];
	NRformat_loc nrfl = { nz_loc, m_loc, MyMinEq, const_cast<double*>(Ay), Aq, const_cast<int*>(Aj) }; // Note: The initialization order of Aq and Aj is different of the order of Ai and Ap for the NCformat!
	SuperMatrix  a    = { /*Stype=CRMatrix_dist*/SLU_NR_loc, /*Dtype=double*/SLU_D, /*Mtype=general*/SLU_GE, m, n, &nrfl };

	// 2) Initialize the SuperLU process grid -------------------------------------------
	
    gridinfo_t           grid;  // grid info
	int                  nprow; // number of processors in each row of grid
	int                  npcol; // number of processors in each column of grid
	Util::FindBestSquare (nProcs, nprow,npcol);
	superlu_gridinit     (MPI::COMM_WORLD, nprow, npcol, &grid);

	// Bail out if I do not belong in the grid
	if (grid.iam>=nprow*npcol) throw new Fatal(_("SuperLUd::Gesv: SuperLU's grid was not created properly: grid.iam=%d > nprow*npcol=%d"),grid.iam,nprow*npcol);

	// 3) Solve the linear system -------------------------------------------------------

	// Set the input options
    superlu_options_t opts;
    opts.Fact              = DOFACT;
    opts.Equil             = YES;
    opts.ColPerm           = MMD_AT_PLUS_A; // NATURAL
    opts.Trans             = NOTRANS;
    opts.IterRefine        = NOREFINE; // DOUBLE
    opts.PrintStat         = NO;
    opts.RowPerm           = NOROWPERM; // LargeDiag
    opts.ReplaceTinyPivot  = YES;
    opts.SolveInitialized  = NO;
    opts.RefineInitialized = NO;

	// Initialize ScalePermstruct and LUstruct
    ScalePermstruct_t   scp;
    LUstruct_t          lus;
	ScalePermstructInit (m, n, &scp);
	LUstructInit        (m, n, &lus);

	// Initialize the statistics variables
    SuperLUStat_t stat;
	PStatInit     (&stat);

	// Call the linear equation solver
	int           info;
	double        berr; // size == nrhs
    SOLVEstruct_t svs;
	pdgssvx       (&opts, &a, &scp, &GlobalB.data[MyMinEq], /*ldb*/m, /*nrhs*/1, &grid, &lus, &svs, &berr, &stat, &info);

#ifdef DO_DEBUGx
	int rank = MPI::COMM_WORLD.Get_rank();
	//LinAlg::Matrix<double> D; AugLocalA.GetDense(D);
	std::ostringstream oss;
	//oss << "Processor #" << rank << " LocalA=\n" << D << endl;
	//oss << "Processor #" << rank << " LocalA=\n" << AugLocalA << endl;
	oss<<"Processor #" << rank << " ";
	oss<<"Aq={"; for (int i=0; i<m_loc+1; ++i) { oss<<Aq[i];                 if (i!=m_loc)    oss<<","; }; oss<<"} ";
	oss<<"Aj={"; for (int i=0; i<nz_loc;  ++i) { oss<<AugLocalA.Aj(start+i); if (i!=nz_loc-1) oss<<","; }; oss<<"} ";
	oss<<"Ay={"; for (int i=0; i<nz_loc;  ++i) { oss<<AugLocalA.Ay(start+i); if (i!=nz_loc-1) oss<<","; }; oss<<"} ";
	oss<<"B={";  for (int i=0; i<m_loc;   ++i) { oss<<GlobalB(MyMinEq+i);    if (i!=m_loc-1)  oss<<","; }; oss<<"}\n";
	cout<<oss.str();
	//dPrint_CompRowLoc_Matrix_dist(&a);
#endif

	// Print the statistics
	//PStatPrint (&opts, &stat, &grid);

	// 4) Deallocate storage ------------------------------------------------------------

	PStatFree                                 (&stat);
	ScalePermstructFree                       (&scp);
	Destroy_LU                                (n, &grid, &lus);
	LUstructFree                              (&lus);
	if (opts.SolveInitialized) dSolveFinalize (&opts, &svs);
	delete [] Aq;

	// 5) Release the SuperLU process grid ----------------------------------------------
	
	superlu_gridexit (&grid);

}

}; // namespace SuperLU

#endif // MECHSYS_SUPERLUD_H
