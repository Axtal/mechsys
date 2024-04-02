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

#ifndef MECHSYS_SPARSE_CRMATRIX_H
#define MECHSYS_SPARSE_CRMATRIX_H

// STL
#include <iostream>
#include <sstream>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/sparse_triplet.h>
#include <mechsys/util/array.h>
#include <mechsys/util/numstreams.h>

#ifdef DO_DEBUG
  using std::cout;
  using std::endl;
#endif

namespace Sparse
{

/** %Sparse matrix in Compressed-row (CR) format (sorted).

    Example:
	\verbatim

      rows:     _                         _   Ay:        Aj:     Aq:
       i=0 --> |   2    3    0    0    0   |   2  3      0 1      0
       i=1 --> |   3    0    4    0    6   |   3  4  6   0 2 4    2
  A =  i=2 --> |   0   -1   -3    2    0   |  -1 -3  2   1 2 3    5
       i=3 --> |   0    0    1    0    0   |   1         2        8
       i=4 --> |_  0    4    2    0    1  _|   4  2  1   1 2 4    9
                                                                 12
        n  = 5   (nrows==ncols)
        nz = 12  (number of non-zero values)

        Ay.size == nz      == 12 : Non-zero values
        Aj.size == nz      == 12 : Column indexes of the non-zero values
        Aq.size == nrows+1 == 6  : Pointers (indexes) to the columns indexes of non-zero values inside Aj
	\endverbatim

  Example:
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsparse.cpp?view=markup">tsparse.cpp Test sparse matrices</a>
*/
template<typename Value_T, typename Index_T>
class CRMatrix
{
public:
	/** Constructor (must call Set() method later). */
	 CRMatrix();

	/** Destructor. */
	~CRMatrix();

	/** Alternative constructor.
	 * Convert triplet (with duplicates) to compressed-row format
	 */
	CRMatrix(Triplet<Value_T,Index_T> const & T);

	// Methods
	void    Set  (Triplet<Value_T,Index_T> const & T);                                                                         ///< (Re)set this matrix for a given triplet (with duplicates)
	void    Set  (Index_T nRows, Index_T nCols, Index_T nZ, String const & strAy, String const & strAj, String const & strAq); ///< nRows: Number of rows, nCols: Number of columns, nZ: Number of non-zero values, Ay: Non-zero values, Aj: Columns indexes, Aq: Pointers to the column indexes
	Index_T Rows () const { return _nrows; } ///< Number of rows
	Index_T Cols () const { return _ncols; } ///< Number of columns
	Index_T nZ   () const { return _nz;    } ///< Number of non-zero values

	// Access methods
	Value_T         Ay       (Index_T iNz ) const; ///< Returns non-zero values
	Index_T         Aj       (Index_T iNz ) const; ///< Returns the columns indexes
	Index_T         Aq       (Index_T iRow) const; ///< Returns the pointers to columns indexes
	Value_T       * GetAyPtr ();                   ///< Access values (write)
	Index_T       * GetAjPtr ();                   ///< Access columns indexes (write)
	Index_T       * GetAqPtr ();                   ///< Access pointers to columns indexes (write)
	Value_T const * GetAyPtr () const;             ///< Access values (read)
	Index_T const * GetAjPtr () const;             ///< Access columns indexes (read)
	Index_T const * GetAqPtr () const;             ///< Access pointers to columns indexes (read)

	// Auxiliar methods
	void                    GetDense (Mat_t & D) const;                    ///< Convert this sparse structure to a Dense matrix
	void                    SetNS    (Util::NumStream & NS) { _ns = &NS; } ///< Set the NumStream, a structure to aid format output of numbers
	Util::NumStream const & NS       () const           { return (*_ns); } ///< Return the NumStream, a structure to aid format output of numbers

private:
	// Data
	Index_T           _nrows; ///< Number of rows
	Index_T           _ncols; ///< Number of columns
	Index_T           _nz;    ///< Number of non-zero values
	Value_T         * _Ay;    ///< Values (data).                   Size == nz
	Index_T         * _Aj;    ///< Columns indexes.                 Size == nz
	Index_T         * _Aq;    ///< Pointers to the columns indexes. Size == nrows + 1
	Util::NumStream * _ns;    ///< Structure to aid format output of numbers

}; // class CRMatrix


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor & Destructor

template<typename Value_T, typename Index_T>
inline CRMatrix<Value_T,Index_T>::CRMatrix()
	: _nrows (0),
	  _ncols (0),
	  _nz    (0),
	  _Ay    (NULL),
	  _Aj    (NULL),
	  _Aq    (NULL),
	  _ns    (&Util::_8s)
{}

template<typename Value_T, typename Index_T>
inline CRMatrix<Value_T,Index_T>::~CRMatrix()
{
	if (_Ay!=NULL) delete [] _Ay;
	if (_Aj!=NULL) delete [] _Aj;
	if (_Aq!=NULL) delete [] _Aq;
}

// Alternative constructor

template<typename Value_T, typename Index_T>
inline CRMatrix<Value_T,Index_T>::CRMatrix(Triplet<Value_T,Index_T> const & T)
	: _Ay (NULL),
	  _Aj (NULL),
	  _Aq (NULL),
	  _ns (&Util::_8s)
{
	Set(T);
}

// Methods

template<typename Value_T, typename Index_T>
inline void CRMatrix<Value_T,Index_T>::Set(Triplet<Value_T,Index_T> const & T)
{
	// Clean
	if (_Ay!=NULL) delete [] _Ay;
	if (_Aj!=NULL) delete [] _Aj;
	if (_Aq!=NULL) delete [] _Aq;

	// Allocate this matrix (With room for the duplicates. This space will be unused by the resulting CRMatrix)
	_nrows = T.Rows();
	_ncols = T.Cols();
	_nz    = T.Top();
	_Ay    = new Value_T [_nz];
	_Aj    = new Index_T [_nz];
	_Aq    = new Index_T [_nrows+1];

	// Temporary variables
	Index_T p, p1, p2, pdest; // pointers (indexes)
	Index_T nn = _nrows;  if (_ncols>nn) nn = _ncols;

	// Allocate workspaces
	Index_T * rq = new Index_T [_nrows+1]; // temporary row format
	Index_T * rj = new Index_T [_nz];      // temporary row format
	Value_T * ry = new Value_T [_nz];      // temporary row format
	Index_T * rc = new Index_T [_nrows];   // row count
	Index_T * w  = new Index_T [nn];       // workspace
	if (!rq || !rj || !ry || !w) throw new Fatal(_("Sparse::CRMatrix::CRMatrix(nRows,nCols,T): Out of memory."));

	// Count the entries in each row (also counting duplicates)
	for (Index_T i=0; i<_nrows; i++) w[i]=0; // use w as workspace for row counts (including duplicates)
	for (Index_T k=0; k<_nz; k++)
	{
		Index_T i = T.Ai(k);
		Index_T j = T.Aj(k);
		if (i<0 || i>=_nrows || j<0 || j>=_ncols) throw new Fatal(_("Sparse::CRMatrix::CRMatrix(nRows,nCols,T): (i<0 || i>=nRows || j<0 || j>=nCols) failed => triplet corresponds to an invalid matrix."));
		w[i]++;
	}

	// Compute the row pointers (with duplicates)
	rq[0] = 0;
	for (Index_T i=0; i<_nrows; i++)
	{
		rq[i+1] = rq[i] + w[i];
		w [i]   = rq[i];       // w is now equal to the row pointers
	}

	// Construct the row form (with duplicates)
	for (Index_T k=0 ; k<_nz; k++)
	{
		p     = w[T.Ai(k)]++; // rq stays the same, but w[i] is advanced to the start of row i+1
		rj[p] = T.Aj(k);
		ry[p] = T.Ax(k);
	}

	// Sum up duplicates (but do not change the size)
	for (Index_T j=0; j<_ncols; j++) w[j]=-1; // use w[j] to hold position in Aj and Ay of a_ij, for row i [
	for (Index_T i=0; i<_nrows; i++)
	{
		p1    = rq[i];
		p2    = rq[i+1];
		pdest = p1;
		for (p=p1; p<p2; p++) // At this point, w[j]<p1 holds true for all columns j, because Aj and Ay are stored in row oriented order
		{
			Index_T j  = rj[p];
			Index_T pj = w[j];
			if (pj>=p1) // this column index, j, is already in row i, at position pj
				ry[pj] += ry[p]; // sum the entry
			else // keep the entry
			{
				w[j] = pdest; // also keep track in w[j] of position of a_ij for case above
				if (pdest!=p) // no need to move the entry if pdest is equal to p
				{
					rj[pdest] = j;
					ry[pdest] = ry[p];
				}
				pdest++;
			}
		}
		rc[i] = pdest - p1; // will have the number of elemens without duplicates in each row (row count)
	} // done using w for position of a_ij ]

	// Construct the row form (without duplicates)
	_nz    = 0;
	_Aq[0] = 0;
	for (Index_T i=0; i<_nrows; i++)
	{
		_Aq[i+1] = _Aq[i] + rc[i];
		for (p=rq[i]; p<rq[i]+rc[i]; p++)
		{
			_Aj[_nz] = rj[p];
			_Ay[_nz] = ry[p];
			_nz++;
		}
	}

	// Cleanup
	delete [] rq;
	delete [] rj;
	delete [] ry;
	delete [] rc;
	delete [] w;

}

template<typename Value_T, typename Index_T>
inline void CRMatrix<Value_T,Index_T>::Set(Index_T nRows, Index_T nCols, Index_T nZ, String const & strAy, String const & strAj, String const & strAq)
{
	Array<Value_T> ay;
	Array<Index_T> aj;
	Array<Index_T> aq;

	// Values => ay
	std::istringstream iss_y(strAy);
	Value_T value;
	while (iss_y>>value) { ay.Push(value); }

	// Column indexes => aj
	std::istringstream iss_j(strAj);
	Index_T index;
	while (iss_j>>index) { aj.Push(index); }

	// Pointers to column indexes => aq
	std::istringstream iss_q(strAq);
	Index_T pointer;
	while (iss_q>>pointer) { aq.Push(pointer); }

#ifdef DO_DEBUG
	cout << "ay = " << strAy << " = "; for (size_t i=0; i<ay.Size(); i++) { cout << ay[i] << " "; } cout<<endl;
	cout << "aj = " << strAj << " = "; for (size_t i=0; i<aj.Size(); i++) { cout << aj[i] << " "; } cout<<endl;
	cout << "aq = " << strAq << " = "; for (size_t i=0; i<aq.Size(); i++) { cout << aq[i] << " "; } cout<<endl;
#endif

	if (static_cast<Index_T>(ay.Size())      !=nZ     ) throw new Fatal(_("Sparse::CRMatrix::Set: The size (%d) of strAy (values) must be equal to nz(=%d), where nz is the number of non-zero values."),ay.Size(),nZ);
	if (static_cast<Index_T>(aj.Size())      !=nZ     ) throw new Fatal(_("Sparse::CRMatrix::Set: The size (%d) of strAj (column indexes) must be equal to nz(=%d), where nz is the number of non-zero values."),aj.Size(),nZ);
	if (static_cast<Index_T>(aq.Size())      !=nRows+1) throw new Fatal(_("Sparse::CRMatrix::Set: The size (%d) of strAq (pointers to the column indexes) must be equal to nrows+1(=%d), where nrows is the number of rows."),aq.Size(),nRows+1);
	if (static_cast<Index_T>(aq[aq.Size()-1])!=nZ     ) throw new Fatal(_("Sparse::CRMatrix::Set: The last value of Aq(=%d) in strAq (pointers to the row indexes) must be equal to nz(=%d), where nz is the number of non-zero values."),aq[aq.Size()-1],nZ);

	_nrows = nRows;
	_ncols = nCols;
	_nz    = nZ;

	_Ay = new Value_T [nZ];
	_Aj = new Index_T [nZ];
	_Aq = new Index_T [nRows+1];

	for (Index_T k=0; k<nZ;      ++k) _Ay[k] = ay[k];
	for (Index_T k=0; k<nZ;      ++k) _Aj[k] = aj[k];
	for (Index_T k=0; k<nRows+1; ++k) _Aq[k] = aq[k];
	
}

// Access methods

template<typename Value_T, typename Index_T>
inline Value_T CRMatrix<Value_T,Index_T>::Ay(Index_T iNz) const
{
#ifndef DNDEBUG
	if ((iNz<0) || (iNz>=_nz)) throw new Fatal(_("Sparse::CRMatrix::Ay: Index for a non-zero value, iNz(=%d), must be greater than/or equal to zero and smaller than nz(=%d)"),iNz,_nz);
#endif
	return _Ay[iNz];
}

template<typename Value_T, typename Index_T>
inline Index_T CRMatrix<Value_T,Index_T>::Aj(Index_T iNz) const
{
#ifndef DNDEBUG
	if ((iNz<0) || (iNz>=_nz)) throw new Fatal(_("Sparse::CRMatrix::Aj: Index for an index to a non-zero value, iNz(=%d), must be greater than/or equal to zero and smaller than nz(=%d)"),iNz,_nz);
#endif
	return _Aj[iNz];
}

template<typename Value_T, typename Index_T>
inline Index_T CRMatrix<Value_T,Index_T>::Aq(Index_T iRow) const
{
#ifndef DNDEBUG
	if ((iRow<0) || (iRow>_nrows)) throw new Fatal(_("Sparse::CRMatrix::Aq: Index for a row, iRow(=%d), must be greater than/or equal to zero and smaller than nRows+1(=%d)"),iRow,_nrows+1);
#endif
	return _Aq[iRow];
}

template<typename Value_T, typename Index_T>
inline Value_T * CRMatrix<Value_T,Index_T>::GetAyPtr() 
{
#ifndef DNDEBUG
	if (_Ay==NULL) throw new Fatal(_("Sparse::CRMatrix::GetAyPtr: The matrix must be set prior to get the pointer to Ay (non-zero values)."));
#endif
	return _Ay;
}

template<typename Value_T, typename Index_T>
inline Index_T * CRMatrix<Value_T,Index_T>::GetAjPtr()
{
#ifndef DNDEBUG
	if (_Aj==NULL) throw new Fatal(_("Sparse::CRMatrix::GetAjPtr: The matrix must be set prior to get the pointer to Aj (column indexes of non-zero values)."));
#endif
	return _Aj;
}

template<typename Value_T, typename Index_T>
inline Index_T * CRMatrix<Value_T,Index_T>::GetAqPtr()
{
#ifndef DNDEBUG
	if (_Aq==NULL) throw new Fatal(_("Sparse::CRMatrix::GetAqPtr: The matrix must be set prior to get the pointer to Aq (pointers to the column indexes of non-zero values)."));
#endif
	return _Aq;
}

template<typename Value_T, typename Index_T>
inline Value_T const * CRMatrix<Value_T,Index_T>::GetAyPtr() const
{
#ifndef DNDEBUG
	if (_Ay==NULL) throw new Fatal(_("Sparse::CRMatrix::GetAyPtr: The matrix must be set prior to get the pointer to Ay (non-zero values)."));
#endif
	return _Ay;
}

template<typename Value_T, typename Index_T>
inline Index_T const * CRMatrix<Value_T,Index_T>::GetAjPtr() const
{
#ifndef DNDEBUG
	if (_Aj==NULL) throw new Fatal(_("Sparse::CRMatrix::GetAjPtr: The matrix must be set prior to get the pointer to Aj (column indexes of non-zero values)."));
#endif
	return _Aj;
}

template<typename Value_T, typename Index_T>
inline Index_T const * CRMatrix<Value_T,Index_T>::GetAqPtr() const
{
#ifndef DNDEBUG
	if (_Aq==NULL) throw new Fatal(_("Sparse::CRMatrix::GetAqPtr: The matrix must be set prior to get the pointer to Aq (pointers to the column indexes of non-zero values)."));
#endif
	return _Aq;
}

// Auxiliar methods

template<typename Value_T, typename Index_T>
inline void CRMatrix<Value_T,Index_T>::GetDense(Mat_t & D) const
{
	D.Resize(_nrows,_ncols);
	for (Index_T i=0; i<_nrows; ++i)
	{
		Index_T start = _Aq[i];
		Index_T endpo = _Aq[i+1]; // == end+1 (po: plus one)
		for (Index_T q=start; q<endpo; ++q)
			D(i,_Aj[q]) = _Ay[q];
	}
}


/** Outputs a sparse matrix in compressed-row format. */
template<typename Value_T, typename Index_T>
std::ostream & operator<< (std::ostream & os, const Sparse::CRMatrix<Value_T,Index_T> & M)
{
	os << "Ay =\n"; for (Index_T k=0; k<M.nZ();     ++k) os << M.NS()  <<M.Ay(k) << " "; os<<std::endl;
	os << "Aj =\n"; for (Index_T k=0; k<M.nZ();     ++k) os << Util::_3<<M.Aj(k) << " "; os<<std::endl;
	os << "Aq =\n"; for (Index_T i=0; i<M.Rows()+1; ++i) os << Util::_3<<M.Aq(i) << " "; os<<std::endl;
	os << std::endl;
#ifdef DO_DEBUG
	/*
	os << "Rows=" << M.Rows() << ", Cols=" << M.Cols() << ", nZ=" << M.nZ() << std::endl;
	for (Index_T i=0; i<M.Rows(); ++i)
	{
		os << "Row # " << i << std::endl;
		Index_T start = M.Aq(i);
		Index_T endpo = M.Aq(i+1); // end+1 (po: plus one)
		for (Index_T q=start; q<endpo; ++q)
		{
			Index_T j = M.Aj(q);
			os << "j=" << j << ", value=" << M.Ay(q) << std::endl;
		}
		os << std::endl;
	}
	*/
#endif
	return os;
}


}; // namespace Sparse

#endif // MECHSYS_SPARSE_CRMATRIX_H
