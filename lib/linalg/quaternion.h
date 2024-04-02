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
 * You should have received A copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_DEM_QUATERNION_H
#define MECHSYS_DEM_QUATERNION_H

// Std Lib
#include <cmath>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include <mechsys/linalg/matvec.h>


typedef blitz::TinyVector<double,4> Quaternion_t;

inline void NormalizeRotation (double Theta, Vec3_t const & Axis, Quaternion_t & C)
{
    if (norm(Axis)<1.0e-12) throw new Fatal("Quaternion: the norm of the axis is too small, please chose a different one");
    Vec3_t A = Axis/norm(Axis);
    C(0)     = cos(Theta/2.0);
    C(1)     = A(0)*sin(Theta/2.0);
    C(2)     = A(1)*sin(Theta/2.0);
    C(3)     = A(2)*sin(Theta/2.0);
}

inline void Conjugate (Quaternion_t const & A, Quaternion_t & C)
{
    C(0) =  A(0);
    C(1) = -A(1);
    C(2) = -A(2);
    C(3) = -A(3);
}

inline void GetVector (Quaternion_t const & A, Vec3_t & C)
{
    C(0) = A(1);
    C(1) = A(2);
    C(2) = A(3);
}

inline void SetQuaternion (double Scalar, Vec3_t const & A, Quaternion_t & C)
{
    C(0) = Scalar;
    C(1) = A(0);
    C(2) = A(1);
    C(3) = A(2);
}

inline void QuaternionProduct (Quaternion_t const & A, Quaternion_t const & B, Quaternion_t & C)
{
    Vec3_t t1,t2;
    GetVector (A,t1);
    GetVector (B,t2);
    double scalar = A(0)*B(0) - dot(t1,t2);
    Vec3_t vector = A(0)*t2 + B(0)*t1 + cross(t1,t2);
    SetQuaternion (scalar,vector,C);
}

inline void Rotation (Vec3_t const & A, Quaternion_t const & B, Vec3_t & C)
{
    Quaternion_t t1,t2,t3;
    SetQuaternion     (0.0,A,t1);
    QuaternionProduct (B,t1,t2);
    Conjugate         (B,t3);
    QuaternionProduct (t2,t3,t1);
    GetVector         (t1,C);
}

inline void RotationMatrix (Quaternion_t const & B, Mat3_t & C)
{
    double s = norm(B);
    C(0,0) = 1.0 - 2*s*(B(2)*B(2)+B(3)*B(3)); C(0,1) = 2*s*(B(1)*B(2) - B(3)*B(0)); C(0,2) = 2*s*(B(1)*B(3) + B(2)*B(0));
    C(1,1) = 1.0 - 2*s*(B(1)*B(1)+B(3)*B(3)); C(1,0) = 2*s*(B(1)*B(2) + B(3)*B(0)); C(1,2) = 2*s*(B(2)*B(3) - B(1)*B(0));
    C(2,2) = 1.0 - 2*s*(B(1)*B(1)+B(2)*B(2)); C(2,0) = 2*s*(B(1)*B(3) - B(2)*B(0)); C(2,1) = 2*s*(B(2)*B(3) + B(1)*B(0));
}

#ifdef USE_CUDA
//////////////////////////////// CUDA IMPLEMENTATION /////////////////////////////
__host__ __device__ void NormalizeRotation (real Theta, real3 const & Axis, real4 & C)
{
    //if (norm(Axis)<1.0e-12) throw new Fatal("Quaternion: the norm of the axis is too small, please chose a different one");
    real3 A = Axis/norm(Axis);
    C.x     = A.x*sin(Theta/2.0);
    C.y     = A.y*sin(Theta/2.0);
    C.z     = A.z*sin(Theta/2.0);
    C.w     = cos(Theta/2.0);
}

__host__ __device__ void Conjugate (real4 const & A, real4 & C)
{
    C.x = -A.x;
    C.y = -A.y;
    C.z = -A.z;
    C.w =  A.w;
}

__host__ __device__ void GetVector (real4 const & A, real3 & C)
{
    C.x = A.x;
    C.y = A.y;
    C.z = A.z;
}

__host__ __device__ void SetQuaternion (real Scalar, real3 const & A, real4 & C)
{
    C.x = A.x;
    C.y = A.y;
    C.z = A.z;
    C.w = Scalar;
}

__host__ __device__ void QuaternionProduct (real4 const & A, real4 const & B, real4 & C)
{
    real3 t1,t2;
    GetVector (A,t1);
    GetVector (B,t2);
    real scalar = A.w*B.w - dotreal3(t1,t2);
    real3 vector = A.w*t2 + B.w*t1 + cross(t1,t2);
    SetQuaternion (scalar,vector,C);
}

__host__ __device__ void Rotation (real3 const & A, real4 const & B, real3 & C)
{
    real4 t1,t2,t3;
    SetQuaternion     (0.0,A,t1);
    QuaternionProduct (B,t1,t2);
    Conjugate         (B,t3);
    QuaternionProduct (t2,t3,t1);
    GetVector         (t1,C);
}

__host__ __device__ real norm(real4 const & a)
{
    return sqrt(a.x*a.x+a.y*a.y+a.z*a.z+a.w*a.w);
}

__host__ __device__ real4 operator * (real const & a, real4 const & b)
{
    return make_real4(a*b.x,a*b.y,a*b.z,a*b.w);
}

__host__ __device__ real4 operator / (real4 const & a, real const & b)
{
    return make_real4(a.x/b,a.y/b,a.z/b,a.w/b);
}

__host__ __device__ inline real4 operator + (real4 const & a, real4 const & b)
{
    return make_real4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w);
}
#endif 

#endif // MECHSYS_DEM_QUATERNION_H
