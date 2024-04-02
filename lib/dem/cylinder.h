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

#ifndef MECHSYS_DEM_CYLINDER_H
#define MECHSYS_DEM_CYLINDER_H

// Std lib
#include <cmath>

// MechSys
#include <mechsys/dem/torus.h>

namespace DEM
{

class Cylinder
{
public:
    // Constructor
    Cylinder (Torus * T0, Torus * T1, Vec3_t const * Y0, Vec3_t const * Y1);  ///< Initial T0 and final T1 tori
    Cylinder (Torus & T0, Torus & T1, Vec3_t const & Y0, Vec3_t const & Y1);  ///< Initial T0 and final T1 tori

    // Methods
    void UpdatedL ();                                                                                     ///< UdatedL for the tori and cylinder

    // Data
    Torus       * T0; ///< First torus
    Torus       * T1; ///< Second torus
    Vec3_t const *Y0; ///< Lower point of the first torus
    Vec3_t const *Y1; ///< Lower point of the second torus
    Vec3_t        X0; ///< Center of the first  torus
    Vec3_t        X1; ///< Center of the second torus
    double      Dmax;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Cylinder::Cylinder (Torus * TheT0, Torus * TheT1, Vec3_t const * TheY0, Vec3_t const * TheY1)
    : T0(TheT0), T1(TheT1)
{
    X0 = *T0->X0;
    X1 = *T1->X0;
    Y0 = TheY0;
    Y1 = TheY1;
    T0->X0 = &X0;
    T1->X0 = &X1;
    Dmax = sqrt(0.25*dot(X0 - X1,X0 - X1) + pow(std::max(T0->R,T1->R),2));
}

inline Cylinder::Cylinder (Torus & TheT0, Torus & TheT1, Vec3_t const & TheY0, Vec3_t const & TheY1)
    : T0(&TheT0), T1(&TheT1)
{
    X0 = *T0->X0;
    X1 = *T1->X0;
    Y0 = &TheY0;
    Y1 = &TheY1;
    T0->X0 = &X0;
    T1->X0 = &X1;
    Dmax = sqrt(0.25*dot(X0 - X1,X0 - X1) + pow(std::max(T0->R,T1->R),2));
}

inline void Cylinder::UpdatedL ()
{
    X0 = 0.5*(*T0->X1 + *Y0);
    X1 = 0.5*(*T1->X1 + *Y1);
}
}
#endif // MECHSYS_DEM_CYLINDER_H
