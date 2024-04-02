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

#ifndef MECHSYS_SPARSE_TRIPLET_H
#define MECHSYS_SPARSE_TRIPLET_H

// STL
#include <iostream>
#include <fstream>

// MechSys
#include <mechsys/util/numstreams.h>

namespace Sparse
{

/** %Sparse matrix in %Triplet format (Duplicates allowed; May be unsorted).

    Example:
    \verbatim
            _                         _   rows:
           |   2    3    0    0    0   |   <--  i=0
           |   3    0    4    0    6   |   <--  i=1
       A = |   0   -1   -3    2    0   |   <--  i=2
           |   0    0    1    0    0   |   <--  i=3
           |_  0    4    2    0    1  _|   <--  i=4

          ____ duplicates
         | |
    Ax = 1 1  3    3  -1  4    4 -3  1  2    2    6  1
    Ai = 0 0  1    0   2  4    1  2  3  4    2    1  4
    Aj = 0 0  0    1   1  1    2  2  2  2    3    4  4
                                                top _|
    size = 13
    top  = 12
    \endverbatim

  Example:
  - <a href="http://cvs.savannah.gnu.org/viewvc/mechsys/mechsys/lib/linalg/tst/tsparse.cpp?view=markup">tsparse.cpp Test sparse matrices</a>
*/
template<typename Value_T, typename Index_T>
class Triplet
{
public:
    /** Constructor (must call AllocSpace() and PushEntry() methods later). */
     Triplet();

    /** Constructing given other 4 triplets. */
     Triplet(Triplet<Value_T,Index_T> const & A, Triplet<Value_T,Index_T> const & B, Triplet<Value_T,Index_T> const & C, Triplet<Value_T,Index_T> const & D);

    /** Destructor. */
    ~Triplet();

    // Methods
    Index_T Rows        () const { return _nrows; }                   ///< Number of rows of the matrix from which this triplet was defined
    Index_T Cols        () const { return _ncols; }                   ///< Number of columns of the matrix from which this triplet was defined
    Index_T Size        () const { return _size;  }                   ///< Return the maximum number of components allowed (memory available)
    Index_T Top         () const { return _top;   }                   ///< Return current position of insertion of new components == current number of components
    void    AllocSpace  (Index_T nRows, Index_T nCols, Index_T Size); ///< Allocate memory for "Size" entries == number of non-zero values, including duplicates. The number of Rows and Columns are only saved to aid further format conversions. The "Size" must be any number, even bigger than Rows*Cols, in case there are duplicates.
    void    PushEntry   (Index_T i, Index_T j, Value_T x);            ///< Insert an entry into the arrays (increase top; top must be smaller than Size)
    void    ResetTop    (Index_T NewTop=0) { _top = NewTop; }         ///< Set Top to NewTop => clear the triplet after NewTop

    // Access methods
    Index_T         Ai       (Index_T k) const; ///< Row index
    Index_T         Aj       (Index_T k) const; ///< Column index
    Value_T         Ax       (Index_T k) const; ///< Non-zero value
    Index_T       * GetAiPtr ();                ///< Access row indexes (write)
    Index_T       * GetAjPtr ();                ///< Access column indexes (write)
    Value_T       * GetAxPtr ();                ///< Access non-zero values (write)
    Index_T const * GetAiPtr () const;          ///< Access row indexes (read)
    Index_T const * GetAjPtr () const;          ///< Access column indexes (read)
    Value_T const * GetAxPtr () const;          ///< Access non-zero values (read)

    // Auxiliar methods
    void                    SetNS      (Util::NumStream & NS) { _ns = &NS; }             ///< Set the NumStream, a structure to aid format output of numbers
    Util::NumStream const & NS         () const           { return (*_ns); }             ///< Return the NumStream, a structure to aid format output of numbers
    void                    WriteSMAT  (char const * FileKey, double Tol=1.0e-14) const; ///< Write .smat file for vismatrix
    void                    IncIndices (Index_T Delta);                                  ///< Increment i,j indices by Delta

private:
    // Data
    Index_T           _nrows; ///< Number of rows of the matrix from which this triplet was defined. This information is only saved to aid further format conversions.
    Index_T           _ncols; ///< Number of columns of the matrix from which this triplet was defined. This information is only saved to aid further format conversions.
    Index_T           _size;  ///< Maximum number of components allowed (memory available)
    Index_T           _top;   ///< Current position of insertion of new components == current number of components-1
    Index_T         * _Ai;    ///< Row indexes (can have duplicates)
    Index_T         * _Aj;    ///< Col indexes (can have duplicates)
    Value_T         * _Ax;    ///< Non zero values (can have duplicates)
    Util::NumStream * _ns;    ///< Structure to aid format output of numbers

}; // class Triplet


/** Add multiplication: Y += M*X */
template<typename Value_T, typename Index_T>
void AddMult (Sparse::Triplet<Value_T,Index_T> const & M, Vec_t const & X, Vec_t & Y)
{
    for (int k=0; k<M.Top(); ++k)
        Y(M.Ai(k)) += M.Ax(k) * X(M.Aj(k)); // Y += M*X
}

/** Subtract multiplication: Y -= M*X */
template<typename Value_T, typename Index_T>
void SubMult (Sparse::Triplet<Value_T,Index_T> const & M, Vec_t const & X, Vec_t & Y)
{
    for (int k=0; k<M.Top(); ++k)
        Y(M.Ai(k)) -= M.Ax(k) * X(M.Aj(k)); // Y -= M*X
}

/** Subtract multiplication: Y -= s*M*X */
template<typename Value_T, typename Index_T>
void SubMult (double s, Sparse::Triplet<Value_T,Index_T> const & M, Vec_t const & X, Vec_t & Y)
{
    for (int k=0; k<M.Top(); ++k)
        Y(M.Ai(k)) -= s * M.Ax(k) * X(M.Aj(k)); // Y -= M*X
}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor & Destructor

template<typename Value_T, typename Index_T>
inline Triplet<Value_T,Index_T>::Triplet()
    : _nrows (0),
      _ncols (0),
      _size  (0),
      _top   (0),
      _Ai    (NULL),
      _Aj    (NULL),
      _Ax    (NULL),
      _ns    (&Util::_8s)
{}

template<typename Value_T, typename Index_T>
inline Triplet<Value_T,Index_T>::~Triplet()
{
    if (_Ai!=NULL) delete [] _Ai;
    if (_Aj!=NULL) delete [] _Aj;
    if (_Ax!=NULL) delete [] _Ax;
}

template<typename Value_T, typename Index_T>
inline Triplet<Value_T,Index_T>::Triplet(Triplet<Value_T,Index_T> const & A, Triplet<Value_T,Index_T> const & B, Triplet<Value_T,Index_T> const & C, Triplet<Value_T,Index_T> const & D)
    : _nrows (0),
      _ncols (0),
      _size  (0),
      _top   (0),
      _Ai    (NULL),
      _Aj    (NULL),
      _Ax    (NULL),
      _ns    (&Util::_8s)
{
    /*
    std::cout << A.Rows() << " " << B.Rows() << " " << C.Rows() << " " << D.Rows() << " " << std::endl;
    std::cout << A.Cols() << " " << B.Cols() << " " << C.Cols() << " " << D.Cols() << " " << std::endl;
    */

    // check
    int m = A.Rows();
    int n = A.Cols();
    if (B.Rows()!=m || B.Cols()!=n ||
        C.Rows()!=m || C.Cols()!=n || 
        D.Rows()!=m || D.Cols()!=n) throw new Fatal("Sparse::Triplet: When construction with other 4 triplets, all triplets must have the same number of rows and columns");

    // set
    AllocSpace (m, n, A.Top()+B.Top()+C.Top()+D.Top());
    for (int k=0; k<A.Top(); ++k) PushEntry (A.Ai(k), A.Aj(k), A.Ax(k));
    for (int k=0; k<B.Top(); ++k) PushEntry (B.Ai(k), B.Aj(k), B.Ax(k));
    for (int k=0; k<C.Top(); ++k) PushEntry (C.Ai(k), C.Aj(k), C.Ax(k));
    for (int k=0; k<D.Top(); ++k) PushEntry (D.Ai(k), D.Aj(k), D.Ax(k));
}

// Methods

template<typename Value_T, typename Index_T>
inline void Triplet<Value_T,Index_T>::AllocSpace(Index_T nRows, Index_T nCols, Index_T Size)
{
    if (_Ai!=NULL) delete [] _Ai;
    if (_Aj!=NULL) delete [] _Aj;
    if (_Ax!=NULL) delete [] _Ax;
    _nrows = nRows;
    _ncols = nCols;
    _size  = Size;
    _top   = 0;
    _Ai    = new Index_T [_size];
    _Aj    = new Index_T [_size];
    _Ax    = new Value_T [_size];
}

template<typename Value_T, typename Index_T>
inline void Triplet<Value_T,Index_T>::PushEntry(Index_T i, Index_T j, Value_T x)
{
    if (_top>=_size) throw new Fatal("Sparse::Triplet::PushEntry: _top (%d) must be smaller than _size (%d)",_top,_size);
    _Ai[_top] = i;
    _Aj[_top] = j;
    _Ax[_top] = x;
    _top++;
}

// Access methods

template<typename Value_T, typename Index_T>
inline Index_T Triplet<Value_T,Index_T>::Ai(Index_T k) const
{
#ifndef DNDEBUG
    if ((k<0) || (k>=_top)) throw new Fatal("Sparse::Triplet::Ai: Index (k=%d) for a row index must be greater than/equal to zero and smaller than/ equal to top(=%d)",k,_top);
    //if ((k<0) || (k>_top)) throw new Fatal("Sparse::Triplet::Ai: Index (k=%d) for a row index must be greater than/equal to zero and smaller than top(=%d)",k,_top);
#endif
    return _Ai[k];
}

template<typename Value_T, typename Index_T>
inline Index_T Triplet<Value_T,Index_T>::Aj(Index_T k) const
{
#ifndef DNDEBUG
    if ((k<0) || (k>_top)) throw new Fatal("Sparse::Triplet::Aj: Index (k=%d) for a column index must be greater than/equal to zero and smaller than top(=%d)",k,_top);
#endif
    return _Aj[k];
}

template<typename Value_T, typename Index_T>
inline Value_T Triplet<Value_T,Index_T>::Ax(Index_T k) const
{
#ifndef DNDEBUG
    if ((k<0) || (k>=_top)) throw new Fatal("Sparse::Triplet::Ax: Index (k=%d) for a non-zero value must be greater than/or equal to zero and smaller than/equal to top(=%d)",k,_top);
    //if ((k<0) || (k>_top)) throw new Fatal("Sparse::Triplet::Ax: Index (k=%d) for a non-zero value must be greater than/or equal to zero and smaller than top(=%d)",k,_top);
#endif
    return _Ax[k];
}

template<typename Value_T, typename Index_T>
inline Index_T * Triplet<Value_T,Index_T>::GetAiPtr()
{
#ifndef DNDEBUG
    if (_Ai==NULL) throw new Fatal("Sparse::Triplet::GetAiPtr: The matrix must be allocated prior to get the pointer to Ai (row indexes).");
#endif
    return _Ai;
}

template<typename Value_T, typename Index_T>
inline Index_T * Triplet<Value_T,Index_T>::GetAjPtr()
{
#ifndef DNDEBUG
    if (_Aj==NULL) throw new Fatal("Sparse::Triplet::GetAjPtr: The matrix must be allocated prior to get the pointer to Aj (column indexes).");
#endif
    return _Aj;
}

template<typename Value_T, typename Index_T>
inline Value_T * Triplet<Value_T,Index_T>::GetAxPtr() 
{
#ifndef DNDEBUG
    if (_Ax==NULL) throw new Fatal("Sparse::Triplet::GetAxPtr: The matrix must be allocated prior to get the pointer to Ax (non-zero values, including duplicates).");
#endif
    return _Ax;
}

template<typename Value_T, typename Index_T>
inline Index_T const * Triplet<Value_T,Index_T>::GetAiPtr() const
{
#ifndef DNDEBUG
    if (_Ai==NULL) throw new Fatal("Sparse::Triplet::GetAiPtr: The matrix must be allocated prior to get the pointer to Ai (row indexes).");
#endif
    return _Ai;
}

template<typename Value_T, typename Index_T>
inline Index_T const * Triplet<Value_T,Index_T>::GetAjPtr() const
{
#ifndef DNDEBUG
    if (_Aj==NULL) throw new Fatal("Sparse::Triplet::GetAjPtr: The matrix must be allocated prior to get the pointer to Aj (column indexes).");
#endif
    return _Aj;
}

template<typename Value_T, typename Index_T>
inline Value_T const * Triplet<Value_T,Index_T>::GetAxPtr() const
{
#ifndef DNDEBUG
    if (_Ax==NULL) throw new Fatal("Sparse::Triplet::GetAxPtr: The matrix must be allocated prior to get the pointer to Ax (non-zero values, including duplicates).");
#endif
    return _Ax;
}

template<typename Value_T, typename Index_T>
inline void Triplet<Value_T,Index_T>::WriteSMAT (char const * FileKey, double Tol) const
{
    // find the number of really non-zero values
    size_t nz = 0;
    for (int k=0; k<Top(); ++k)
    {
        if (fabs(Ax(k))>Tol) nz++;
    }

    // output
    std::ostringstream oss;
    char buf[256];
    sprintf(buf, "%d  %d  %zd\n", Rows(), Cols(), nz);
    oss << buf;
    for (int k=0; k<Top(); ++k)
    {
        if (fabs(Ax(k))>Tol)
        {
            sprintf(buf, "  %d  %d  %g\n", Ai(k), Aj(k), Ax(k));
            oss << buf;
        }
    }

    // write to file
    String fn(FileKey);  fn.append(".smat");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

template<typename Value_T, typename Index_T>
inline void Triplet<Value_T,Index_T>::IncIndices (Index_T Delta)
{
    for (int k=0; k<Top(); ++k)
    {
        _Ai[k] += Delta;
        _Aj[k] += Delta;
    }
}

/** Ouptuts a sparse matrix in triplet format. */
template<typename Value_T, typename Index_T>
std::ostream & operator<< (std::ostream & os, const Sparse::Triplet<Value_T,Index_T> & T)
{
    os << "Rows=" << Util::_3<<T.Rows() << ", Cols=" << Util::_3<<T.Cols() << ", Top=" << Util::_3<<T.Top() << ", Size=" << Util::_3<<T.Size() << std::endl;
    for (Index_T k=0; k<T.Top(); ++k)
        os << Util::_3<< k << " #  i= " << Util::_3<<T.Ai(k) << ", j= " << Util::_3<<T.Aj(k) << ", value=" << T.NS()<< T.Ax(k) << std::endl;
    return os;
}


}; // namespace Sparse

#endif // MECHSYS_SPARSE_TRIPLET_H
