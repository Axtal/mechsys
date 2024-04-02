/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

/* Numerical evaluation of eigenvalues and eigenvectors for NxN SYMMETRIC
 *  matrices using the Jacobi-rotation Method.
 *
 * 2005-09-27 (C): Dorival M. Pedroso - Completely rewritten (3x3 => NxN); based on C numerical recipes
 * 2005-05-10 (C): Dorival M. Pedroso - C++
 * 2004-08-15 (C): Dorival M. Pedroso - Fortran 95
 * 2002-04-15 (C): Dorival M. Pedroso - Fortran 77
 */

#ifndef MECHSYS_LINALG_JACOBIROT_H
#define MECHSYS_LINALG_JACOBIROT_H

// STL
#include <cmath>

// MechSys
#include "mechsys/linalg/matvec.h"

/** Jacobi Transformation of a Symmetric Matrix<double>.
 * The Jacobi method consists of a sequence of orthogonal similarity transformations.
 * Each transformation (a Jacobi rotation) is just a plane rotation designed to annihilate one of the off-diagonal matrix elements.
 * Successive transformations undo previously set zeros, but the off-diagonal elements nevertheless get smaller and smaller.
 * Accumulating the product of the transformations as you go gives the matrix of eigenvectors (Q), while the elements of the final
 * diagonal matrix (A) are the eigenvalues.
 * The Jacobi method is absolutely foolproof for all real symmetric matrices.
 * For matrices of order greater than about 10, say, the algorithm is slower, by a significant constant factor, than the QR method.
 *
 * A = Q * L * Q.T
 *
 * \param A In/Out: is the matrix we seek for the eigenvalues (SYMMETRIC and square, i.e. Rows=Cols)
 * \param Q Out: is a matrix which columns are the eigenvectors
 * \param v Out: is a vector with the eigenvalues
 * \return The number of iterations
 */
template<typename MatrixType, typename VectorType>
inline int _jacobi_rot (int N, MatrixType & A, MatrixType & Q, VectorType & v, double errTol=DBL_EPSILON)
{
    const int    maxIt = 20;          // max number of iterations
    const double Zero  = DBL_EPSILON; // tolerance

    int    j,p,q;
    double theta,tau,t,sm,s,h,g,c;
    double * b = new double [N];
    double * z = new double [N];

    for (p=0; p<N; p++) // initialize Q to the identity matrix.
    {
        for (q=0; q<N; q++) Q(p,q) = 0.0;
        Q(p,p) = 1.0;
    }
    for (p=0; p<N; p++)
    {
        b[p] = v(p) = A(p,p); // initialize b and v to the diagonal of A
        z[p] = 0.0;           // this vector will accumulate terms of the form tapq as in equation (11.1.14).
    }

    for (int it=0; it<maxIt; it++)
    {
        sm = 0.0;
        for (p=0; p<N-1; p++) // sum off-diagonal elements.
        {
            for (q=p+1; q<N; q++)
            sm += fabs(A(p,q));
        }
        if (sm<errTol) // exit point
        {
            delete [] b;
            delete [] z;
            return it+1;
        }
        for (p=0; p<N-1; p++)
        {
            for (q=p+1; q<N; q++)
            {
                h = v(q)-v(p);
                if (fabs(h)<=Zero) t=1.0;
                else
                {
                    theta = 0.5*h/(A(p,q));
                    t     = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                    if (theta < 0.0) t=-t;
                }
                c      = 1.0/sqrt(1.0+t*t);
                s      = t*c;
                tau    = s/(1.0+c);
                h      = t*A(p,q);
                z[p]  -= h;
                z[q]  += h;
                v(p)  -= h;
                v(q)  += h;
                A(p,q) = 0.0;
                for (j=0; j<p; j++)  // case of rotations 0 <= j < p.
                {
                    g      = A(j,p);
                    h      = A(j,q);
                    A(j,p) = g - s*(h+g*tau);
                    A(j,q) = h + s*(g-h*tau);
                }
                for (j=p+1; j<q; j++) // case of rotations p < j < q.
                {
                    g      = A(p,j);
                    h      = A(j,q);
                    A(p,j) = g - s*(h+g*tau);
                    A(j,q) = h + s*(g-h*tau);
                }
                for (j=q+1; j<N; j++) //case of rotations q < j < N.
                {
                    g      = A(p,j);
                    h      = A(q,j);
                    A(p,j) = g - s*(h+g*tau);
                    A(q,j) = h + s*(g-h*tau);
                }
                for (j=0; j<N; j++) // Q matrix
                {
                    g      = Q(j,p);
                    h      = Q(j,q);
                    Q(j,p) = g - s*(h+g*tau);
                    Q(j,q) = h + s*(g-h*tau);
                }
            }
        }
        for (p=0; p<N; p++)
        {
            b[p] += z[p];
            z[p]  = 0.0;   // reinitialize z.
            v(p)  = b[p];  // update v with the sum of tapq,
        }
    }

    delete [] b;
    delete [] z;
    throw new Fatal("JacobiRot: Jacobi rotation dit not converge after %d iterations",maxIt+1);
    return maxIt+1;
}

/** Eigenvalues and Eigenvectors (columns of Q) of Matrix. */
inline int JacobiRot (Mat_t & A, Mat_t & Q, Vec_t & v, double errTol=DBL_EPSILON)
{
    int N = A.num_rows();
    if (N<2)                     throw new Fatal("JacobiRot: Matrix needs to be at least 2x2");
    if (A.num_cols()!=(size_t)N) throw new Fatal("JacobiRot: Matrix must be square");
    Q.change_dim (N,N);
    v.change_dim (N);
    return _jacobi_rot (N, A, Q, v, errTol);
}

/** Eigenvalues and Eigenvectors (columns of Q) of TinyMatrix. */
inline int JacobiRot (Mat3_t & A, Mat3_t & Q, Vec3_t & v, double errTol=DBL_EPSILON)
{
    return _jacobi_rot (/*N*/3, A, Q, v, errTol);
}

/** Eigenvalues and Eigenvectors (columns of Q) of 2nd order tensor in Mandel's basis. */
inline int JacobiRot (Vec_t const & TenA, Mat3_t & Q, Vec3_t & v, double errTol=DBL_EPSILON)
{
    Mat3_t A;
    Ten2Mat (TenA, A);
    return _jacobi_rot (/*N*/3, A, Q, v, errTol);
}

#ifdef USE_BOOST_PYTHON
inline int PyJacobiRot (BPy::list const & A, BPy::list & Q, BPy::list & V, double errTol=DBL_EPSILON)
{
    Mat_t a, q;
    Vec_t v;
    List2Mat (A, a);
    int it = JacobiRot (a, q, v, errTol);
    Mat2List (q, Q);
    Vec2List (v, V);
    return it;
}
#endif

#endif // MECHSYS_LINALG_JACOBIROT_H
