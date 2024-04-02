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

#ifndef MECHSYS_SPARSE_MATRIX_H
#define MECHSYS_SPARSE_MATRIX_H

// STL
#include <iostream>
#include <sstream>  // for istringstream, ostringstream

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/linalg/sparse_triplet.h>
#include <mechsys/util/array.h>
#include <mechsys/util/numstreams.h>

#ifdef DO_DEBUG
  using std::cout;
  using std::endl;
#endif

/** \namespace Sparse %Sparse matrices. */

namespace Sparse
{

/** %Sparse matrix in compressed-column (CC) format (sorted).

    Example:
    \verbatim

     columns:  0    1    2    3    4
               |    |    |    |    |
               V    V    V    V    V
            _                         _   rows:
           |   2    3    0    0    0   |   <--  i=0
           |   3    0    4    0    6   |   <--  i=1
       A = |   0   -1   -3    2    0   |   <--  i=2
           |   0    0    1    0    0   |   <--  i=3
           |_  0    4    2    0    1  _|   <--  i=4

     non-zero:
              Ax:  Ax:  Ax:  Ax:  Ax:
               2    3    4    2    6
               3   -1   -3         1
                    4    1
                         2                     =>  Ax = 2  3    3  -1  4    4 -3  1  2    2    6  1
              Ai:  Ai:  Ai:  Ai:  Ai:                   |       |           |             |    |
               0    0    1    2    1                    |       |           |             |    |
               1    2    2         4                    |       |           |             |    |
                    4    3                              |       |           |             |    |
                         4                     =>  Ai = 0  1    0   2  4    1  2  3  4    2    1  4
              Ap:  Ap:  Ap:  Ap:  Ap:                   |       |           |             |    |
               0    2    5    9    10    12    =>  Ap = 0       2           5             9   10     12

        n  = 5   (nrows==ncols)
        nz = 12  (number of non-zero values)

        Ax.size == nz      == 12 : Non-zero values
        Ai.size == nz      == 12 : Row indexes of the non-zero values
        Ap.size == ncols+1 == 6  : Pointers (indexes) to the row indexes of non-zero values inside Ai
    \endverbatim

  Example:
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsparse.cpp?view=markup">tsparse.cpp Test sparse matrices</a>
*/
template<typename Value_T, typename Index_T>
class Matrix
{
public:
    /** Constructor (must call Set() method later). */
     Matrix();

    /** Destructor. */
    ~Matrix();

    /** Alternative constructor.
     * Convert triplet (with duplicates) to compressed-column format
     */
    Matrix(Triplet<Value_T,Index_T> const & T);

    // Methods
    void    Set  (Triplet<Value_T,Index_T> const & T);                                                                         ///< (Re)set this matrix for a given triplet (with duplicates)
    void    Set  (Index_T nRows, Index_T nCols, Index_T nZ, String const & strAx, String const & strAi, String const & strAp); ///< nRows: Number of rows, nCols: Number of columns, nZ: Number of non-zero values, Ax: Non-zero values, Ai: Row indexes, Ap: Pointers to row indexes
    Index_T Rows () const { return _nrows; } ///< Number of rows
    Index_T Cols () const { return _ncols; } ///< Number of columns
    Index_T nZ   () const { return _nz;    } ///< Number of non-zero values

    // Access methods
    Value_T         Ax       (Index_T iNz ) const; ///< Returns non-zero values
    Index_T         Ai       (Index_T iNz ) const; ///< Returns the row indexes
    Index_T         Ap       (Index_T iCol) const; ///< Returns the pointers to row indexes
    Value_T       * GetAxPtr ();                   ///< Access values (write)
    Index_T       * GetAiPtr ();                   ///< Access row indexes (write)
    Index_T       * GetApPtr ();                   ///< Access pointers to row indexes (write)
    Value_T const * GetAxPtr () const;             ///< Access values (read)
    Index_T const * GetAiPtr () const;             ///< Access row indexes (read)
    Index_T const * GetApPtr () const;             ///< Access pointers to row indexes (read)

    // Auxiliar methods
    void                    GetDense  (Mat_t & D) const;                                ///< Convert this sparse structure to a Dense matrix
    void                    SetNS     (Util::NumStream & NS) { _ns = &NS; }             ///< Set the NumStream, a structure to aid format output of numbers
    Util::NumStream const & NS        () const           { return (*_ns); }             ///< Return the NumStream, a structure to aid format output of numbers
    void                    WriteSMAT (char const * FileKey, double Tol=1.0e-14) const; ///< Write .smat file for vismatrix

private:
    // Data
    Index_T           _nrows; ///< Number of rows
    Index_T           _ncols; ///< Number of columns
    Index_T           _nz;    ///< Number of non-zero values
    Value_T         * _Ax;    ///< Values (data).               Size == nz
    Index_T         * _Ai;    ///< Row indexes.                 Size == nz
    Index_T         * _Ap;    ///< Pointers to the row indexes. Size == ncols + 1
    Util::NumStream * _ns;    ///< Structure to aid format output of numbers

}; // class Matrix


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor & Destructor

template<typename Value_T, typename Index_T>
inline Matrix<Value_T,Index_T>::Matrix()
    : _nrows (0),
      _ncols (0),
      _nz    (0),
      _Ax    (NULL),
      _Ai    (NULL),
      _Ap    (NULL),
      _ns    (&Util::_8s)
{}

template<typename Value_T, typename Index_T>
inline Matrix<Value_T,Index_T>::~Matrix()
{
    if (_Ax!=NULL) delete [] _Ax;
    if (_Ai!=NULL) delete [] _Ai;
    if (_Ap!=NULL) delete [] _Ap;
}

// Alternative constructor

template<typename Value_T, typename Index_T>
inline Matrix<Value_T,Index_T>::Matrix(Triplet<Value_T,Index_T> const & T)
    : _Ax (NULL),
      _Ai (NULL),
      _Ap (NULL),
      _ns (&Util::_8s)
{
    Set(T);
}

// Methods

template<typename Value_T, typename Index_T>
inline void Matrix<Value_T,Index_T>::Set(Triplet<Value_T,Index_T> const & T)
{
    /* Based on Prof. Tim Davis' UMFPACK::UMF_triplet_map_x  inside umf_triplet.c */

    // Clean
    if (_Ax!=NULL) delete [] _Ax;
    if (_Ai!=NULL) delete [] _Ai;
    if (_Ap!=NULL) delete [] _Ap;

    // Allocate this matrix (with duplicates)
    _nrows = T.Rows();
    _ncols = T.Cols();
    _nz    = T.Top();
    _Ax    = new Value_T [_nz];
    _Ai    = new Index_T [_nz];
    _Ap    = new Index_T [_ncols+1];

    // Temporary variables
    Index_T p, p1, p2, pdest, cp; // pointers (indexes)
    Index_T nn = _nrows;  if (_ncols>nn) nn = _ncols;

    // Allocate workspaces
    Index_T * rp = new Index_T [_nrows+1];  // temporary row form
    Index_T * rj = new Index_T [_nz];       // temporary row form
    Value_T * rx = new Value_T [_nz];       // temporary row form
    Index_T * rc = new Index_T [_nrows];    // row count
    Index_T * w  = new Index_T [nn];        // workspace
    if (!rp || !rj || !rx || !rc || !w) throw new Fatal("Sparse::Matrix::Matrix(n,T): Out of memory.");

    // Count the entries in each row (also counting duplicates)
    for (Index_T i=0; i<_nrows; i++) w[i]=0; // use w as workspace for row counts (including duplicates)
    for (Index_T k=0; k<_nz; k++)
    {
        Index_T i = T.Ai(k);
        Index_T j = T.Aj(k);
        if (i<0 || i>=_nrows || j<0 || j>=_ncols) throw new Fatal("Sparse::Matrix::Matrix(n,T): (i<0 || i>=nRows || j<0 || j>=nCols) failed => triplet corresponds to an invalid matrix.");
        w[i]++;
    }

    // Compute the row pointers
    rp[0] = 0;
    for (Index_T i=0; i<_nrows; i++)
    {
        rp[i+1] = rp[i] + w[i];
        w [i]   = rp[i];       // w is now equal to the row pointers
    }

    // Construct the row form
    for (Index_T k=0 ; k<_nz; k++)
    {
        p     = w[T.Ai(k)]++; // rp stays the same, but w[i] is advanced to the start of row i+1
        rj[p] = T.Aj(k);
        rx[p] = T.Ax(k);
    }

    // Sum up duplicates
    for (Index_T j=0; j<_ncols; j++) w[j]=-1; // use w[j] to hold position in rj and rx of a_ij, for row i [
    for (Index_T i=0; i<_nrows; i++)
    {
        p1    = rp[i];
        p2    = rp[i+1];
        pdest = p1;
        for (p=p1; p<p2; p++) // At this point, w[j]<p1 holds true for all columns j, because rj and rx are stored in row oriented order
        {
            Index_T j  = rj[p];
            Index_T pj = w [j];
            if (pj>=p1) // this column index, j, is already in row i, at position pj
                rx[pj] += rx[p]; // sum the entry
            else // keep the entry
            {
                w[j] = pdest; // also keep track in w[j] of position of a_ij for case above
                if (pdest!=p) // no need to move the entry if pdest is equal to p
                {
                    rj[pdest] = j;
                    rx[pdest] = rx[p];
                }
                pdest++;
            }
        }
        rc[i] = pdest - p1 ;
    } // done using W for position of a_ij ]

    // Count the entries in each column
    for (Index_T j=0; j<_ncols; j++) w[j]=0; // [ use W as work space for column counts of A
    for (Index_T i=0; i<_nrows; i++)
        for (p=rp[i]; p<rp[i]+rc[i]; p++)
            w[rj[p]]++;

    // Create the column pointers
    _Ap[0] = 0;
    for (Index_T j=0; j<_ncols; j++) _Ap[j+1] = _Ap[j] + w[j]; // done using W as workspace for column counts of A ]
    for (Index_T j=0; j<_ncols; j++)   w[j]   = _Ap[j];

    // Construct the column form
    _nz = 0;
    for (Index_T i=0; i<_nrows; i++)
    {
        for (p=rp[i]; p<rp[i] + rc[i]; p++)
        {
            cp      = w[rj[p]]++;
            _Ai[cp] = i;
            _Ax[cp] = rx[p];
            _nz++;
        }
    }

    // Cleanup
    delete [] rp;
    delete [] rj;
    delete [] rx;
    delete [] rc;
    delete [] w;

}

template<typename Value_T, typename Index_T>
inline void Matrix<Value_T,Index_T>::Set(Index_T nRows, Index_T nCols, Index_T nZ, String const & strAx, String const & strAi, String const & strAp)
{
    Array<Value_T> ax;
    Array<Index_T> ai;
    Array<Index_T> ap;

    // Values => ax
    std::istringstream iss_x(strAx);
    Value_T value;
    while (iss_x>>value) { ax.Push(value); }

    // Row indexes => ai
    std::istringstream iss_i(strAi);
    Index_T index;
    while (iss_i>>index) { ai.Push(index); }

    // Pointers to row indexes => ap
    std::istringstream iss_p(strAp);
    Index_T pointer;
    while (iss_p>>pointer) { ap.Push(pointer); }

#ifdef DO_DEBUG
    cout << "ax = " << strAx << " = "; for (size_t i=0; i<ax.Size(); i++) { cout << ax[i] << " "; } cout<<endl;
    cout << "ai = " << strAi << " = "; for (size_t i=0; i<ai.Size(); i++) { cout << ai[i] << " "; } cout<<endl;
    cout << "ap = " << strAp << " = "; for (size_t i=0; i<ap.Size(); i++) { cout << ap[i] << " "; } cout<<endl;
#endif

    if (static_cast<Index_T>(ax.Size())      !=nZ     ) throw new Fatal("Sparse::Matrix::Set: The size (%d) of strAx (values) must be equal to nz(=%d), where nz is the number of non-zero values.",ax.Size(),nZ);
    if (static_cast<Index_T>(ai.Size())      !=nZ     ) throw new Fatal("Sparse::Matrix::Set: The size (%d) of strAi (row indexes) must be equal to nz(=%d), where nz is the number of non-zero values.",ai.Size(),nZ);
    if (static_cast<Index_T>(ap.Size())      !=nCols+1) throw new Fatal("Sparse::Matrix::Set: The size (%d) of strAp (pointers to the row indexes) must be equal to ncols+1(=%d), where ncols is the number of columns.",ap.Size(),nCols+1);
    if (static_cast<Index_T>(ap[ap.Size()-1])!=nZ     ) throw new Fatal("Sparse::Matrix::Set: The last value of Ap(=%d) in strAp (pointers to the row indexes) must be equal to nz(=%d), where nz is the number of non-zero values.",ap[ap.Size()-1],nZ);

    _nrows = nRows;
    _ncols = nCols;
    _nz    = nZ;

    _Ax = new Value_T [nZ];
    _Ai = new Index_T [nZ];
    _Ap = new Index_T [nCols+1];

    for (Index_T k=0; k<nZ;      ++k) _Ax[k] = ax[k];
    for (Index_T k=0; k<nZ;      ++k) _Ai[k] = ai[k];
    for (Index_T k=0; k<nCols+1; ++k) _Ap[k] = ap[k];
    
}

// Access methods

template<typename Value_T, typename Index_T>
inline Value_T Matrix<Value_T,Index_T>::Ax(Index_T iNz) const
{
#ifndef DNDEBUG
    if ((iNz<0) || (iNz>=_nz)) throw new Fatal("Sparse::Matrix::Ax: Index for a non-zero value, iNz(=%d), must be greater than/or equal to zero and smaller than nz(=%d)",iNz,_nz);
#endif
    return _Ax[iNz];
}

template<typename Value_T, typename Index_T>
inline Index_T Matrix<Value_T,Index_T>::Ai(Index_T iNz) const
{
#ifndef DNDEBUG
    if ((iNz<0) || (iNz>=_nz)) throw new Fatal("Sparse::Matrix::Ai: Index for a row index to a non-zero value, iNz(=%d), must be greater than/or equal to zero and smaller than nz(=%d)",iNz,_nz);
#endif
    return _Ai[iNz];
}

template<typename Value_T, typename Index_T>
inline Index_T Matrix<Value_T,Index_T>::Ap(Index_T iCol) const
{
#ifndef DNDEBUG
    if ((iCol<0) || (iCol>_ncols)) throw new Fatal("Sparse::Matrix::Ap: Index for a column, iCol(=%d), must be greater than/or equal to zero and smaller than nCols+1(=%d)",iCol,_ncols+1);
#endif
    return _Ap[iCol];
}

template<typename Value_T, typename Index_T>
inline Value_T * Matrix<Value_T,Index_T>::GetAxPtr() 
{
#ifndef DNDEBUG
    if (_Ax==NULL) throw new Fatal("Sparse::Matrix::GetAxPtr: The matrix must be set prior to get the pointer to Ax (non-zero values).");
#endif
    return _Ax;
}

template<typename Value_T, typename Index_T>
inline Index_T * Matrix<Value_T,Index_T>::GetAiPtr()
{
#ifndef DNDEBUG
    if (_Ai==NULL) throw new Fatal("Sparse::Matrix::GetAiPtr: The matrix must be set prior to get the pointer to Ai (row indexes of non-zero values).");
#endif
    return _Ai;
}

template<typename Value_T, typename Index_T>
inline Index_T * Matrix<Value_T,Index_T>::GetApPtr()
{
#ifndef DNDEBUG
    if (_Ap==NULL) throw new Fatal("Sparse::Matrix::GetApPtr: The matrix must be set prior to get the pointer to Ap (pointers to the row indexes of non-zero values).");
#endif
    return _Ap;
}

template<typename Value_T, typename Index_T>
inline Value_T const * Matrix<Value_T,Index_T>::GetAxPtr() const
{
#ifndef DNDEBUG
    if (_Ax==NULL) throw new Fatal("Sparse::Matrix::GetAxPtr: The matrix must be set prior to get the pointer to Ax (non-zero values).");
#endif
    return _Ax;
}

template<typename Value_T, typename Index_T>
inline Index_T const * Matrix<Value_T,Index_T>::GetAiPtr() const
{
#ifndef DNDEBUG
    if (_Ai==NULL) throw new Fatal("Sparse::Matrix::GetAiPtr: The matrix must be set prior to get the pointer to Ai (row indexes of non-zero values).");
#endif
    return _Ai;
}

template<typename Value_T, typename Index_T>
inline Index_T const * Matrix<Value_T,Index_T>::GetApPtr() const
{
#ifndef DNDEBUG
    if (_Ap==NULL) throw new Fatal("Sparse::Matrix::GetApPtr: The matrix must be set prior to get the pointer to Ap (pointers to the row indexes of non-zero values).");
#endif
    return _Ap;
}

// Auxiliar methods

template<typename Value_T, typename Index_T>
inline void Matrix<Value_T,Index_T>::GetDense(Mat_t & D) const
{
    D.change_dim(_nrows,_ncols);
    set_to_zero(D);
    for (Index_T j=0; j<_ncols; ++j)
    {
        Index_T start = _Ap[j];
        Index_T endpo = _Ap[j+1]; // == end+1
        for (Index_T p=start; p<endpo; ++p)
            D(_Ai[p],j) = _Ax[p];
    }
}

template<typename Value_T, typename Index_T>
inline void Matrix<Value_T,Index_T>::WriteSMAT (char const * FileKey, double Tol) const
{
    // find the number of really non-zero values
    size_t nz = 0;
    for (Index_T j=0; j<_ncols; ++j)
    {
        Index_T start = _Ap[j];
        Index_T endpo = _Ap[j+1]; // == end+1
        for (Index_T p=start; p<endpo; ++p)
            if (fabs(_Ax[p])>Tol) nz++;
    }

    // output
    std::ostringstream oss;
    char buf[256];
    sprintf(buf, "%d  %d  %zd\n", _nrows, _ncols, nz);
    oss << buf;
    for (Index_T j=0; j<_ncols; ++j)
    {
        Index_T start = _Ap[j];
        Index_T endpo = _Ap[j+1]; // == end+1
        for (Index_T p=start; p<endpo; ++p)
        {
            if (fabs(_Ax[p])>Tol)
            {
                sprintf(buf, "  %d  %d  %g\n", _Ai[p], j, _Ax[p]);
            }
        }
    }

    // write to file
    String fn(FileKey);  fn.append(".smat");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

/** Outputs a sparse matrix in compressed-column format. */
template<typename Value_T, typename Index_T>
std::ostream & operator<< (std::ostream & os, const Sparse::Matrix<Value_T,Index_T> & M)
{
    os << "Ax =\n"; for (Index_T k=0; k<M.nZ();     ++k) os << M.NS()  <<M.Ax(k) << " "; os<<std::endl;
    os << "Ai =\n"; for (Index_T k=0; k<M.nZ();     ++k) os << Util::_3<<M.Ai(k) << " "; os<<std::endl;
    os << "Ap =\n"; for (Index_T j=0; j<M.Cols()+1; ++j) os << Util::_3<<M.Ap(j) << " "; os<<std::endl;
    os << std::endl;
#ifdef DO_DEBUG
    /*
    os << "Rows=" << M.Rows() << ", Cols=" << M.Cols() << ", nZ=" << M.nZ() << std::endl;
    for (Index_T j=0; j<M.Cols(); ++j)
    {
        os << "Column # " << j << std::endl;
        Index_T start = M.Ap(j);
        Index_T endpo = M.Ap(j+1); // end+1
        for (Index_T p=start; p<endpo; ++p)
        {
            Index_T i = M.Ai(p);
            os << "i=" << i << ", value=" << M.Ax(p) << std::endl;
        }
        os << std::endl;
    }
    */
#endif
    return os;
}


}; // namespace Sparse

#endif // MECHSYS_SPARSE_MATRIX_H
