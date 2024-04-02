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

#ifndef MECHSYS_DEM_TORUS_H
#define MECHSYS_DEM_TORUS_H

// Std lib
#include <cmath>

// MechSys
#include <mechsys/dem/face.h>

namespace DEM
{

class Torus
{
public:
    // Constructor
    Torus (Vec3_t const * X0, Vec3_t const * X1, Vec3_t const * X2); ///< X0: Center of the torus Xi: endpoints of edge
    Torus (Vec3_t const & X0, Vec3_t const & X1, Vec3_t const & X2); ///< X0: Center of the torus Xi: endpoints of edge

    // Methods

    // Data
    Vec3_t const * X0; ///< Center of the Torus
    Vec3_t const * X1; ///< Right vector of the torus
    Vec3_t const * X2; ///< Left vector of the torus
    double          R; ///< Internal radius of the torus
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Torus::Torus (Vec3_t const * TheX0, Vec3_t const * TheX1, Vec3_t const * TheX2)
    : X0(TheX0), X1(TheX1), X2(TheX2)
{
    R = norm(*X0-*X1);
}

inline Torus::Torus (Vec3_t const & TheX0, Vec3_t const & TheX1, Vec3_t const & TheX2)
    : X0(&TheX0), X1(&TheX1), X2(&TheX2)
{
    R = norm(*X0-*X1);
}
}
#endif // MECHSYS_DEM_TORUS_H
