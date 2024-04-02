/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
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

#ifndef MECHSYS_DEM_BASIC_H
#define MECHSYS_DEM_BASIC_H

// Stl
#include <algorithm>

// Some functions for the mass properties integration
inline double f1(size_t i, Vec3_t const & V0, Vec3_t const & V1, Vec3_t const & V2)
{
    return V0(i)+V1(i)+V2(i);
}

inline double f2(size_t i, Vec3_t const & V0, Vec3_t const & V1, Vec3_t const & V2)
{
    return V0(i)*V0(i)+V0(i)*V1(i)+V1(i)*V1(i)+V2(i)*f1(i,V0,V1,V2);
}

inline double f3(size_t i, Vec3_t const & V0, Vec3_t const & V1, Vec3_t const & V2)
{
    return V0(i)*V0(i)*V0(i)+V0(i)*V0(i)*V1(i)+V0(i)*V1(i)*V1(i)+V1(i)*V1(i)*V1(i)+V2(i)*f2(i,V0,V1,V2);
}

inline double g0(size_t i, Vec3_t const & V0, Vec3_t const & V1, Vec3_t const & V2)
{
    return f2(i,V0,V1,V2)+V0(i)*(f1(i,V0,V1,V2)+V0(i));
}

inline double g1(size_t i, Vec3_t const & V0, Vec3_t const & V1, Vec3_t const & V2)
{
    return f2(i,V0,V1,V2)+V1(i)*(f1(i,V0,V1,V2)+V1(i));
}

inline double g2(size_t i, Vec3_t const & V0, Vec3_t const & V1, Vec3_t const & V2)
{
    return f2(i,V0,V1,V2)+V2(i)*(f1(i,V0,V1,V2)+V2(i));
}

inline void CheckDestroGiro(Vec3_t & xp, Vec3_t & yp, Vec3_t & zp) // Check if the coordinate system is destrogiro and modify the system accordingly for the quaternion calculation
{
    bool destrogiro = dot(cross(xp,yp),zp)>0;
    if (!destrogiro)
    {
        xp = -xp;
        if (1+xp(0)+yp(1)+zp(2)>0) return;
        else 
        {
            xp = -xp;
            yp = -yp;
            if (1+xp(0)+yp(1)+zp(2)>0) return;
            else
            {
                yp = -yp;
                zp = -zp;
                if (1+xp(0)+yp(1)+zp(2)>0) return;
                else throw new Fatal("specialfunctions.h::CheckDestroGiro: The system cannot be transformed by a quaternion operation");
            }
        }
    }
}

inline double ReducedValue(double A, double B) //Reduced Value for two quantities
{
    double result = A*B;
    if (result>0.0) result/=(A+B);
    return result;
}

inline size_t HashFunction(size_t i, size_t j) // Hash function for fast pair lookup
{
    if ((i+j)%2==0) return (i+j+1)*((i+j)/2)+j;
    else            return (i+j)*((i+j+1)/2)+j;
}

#endif //MECHSYS_DEM_BASIC_H


