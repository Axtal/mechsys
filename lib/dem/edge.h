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

#ifndef MECHSYS_DEM_EDGE_H
#define MECHSYS_DEM_EDGE_H

// Std lib
#include <cmath>

// MechSys
#include <mechsys/linalg/quaternion.h>

namespace DEM
{

class Edge
{
public:
    // Constructor
    Edge (Vec3_t * X0, Vec3_t * X1); ///< Xi: endpoints of edge
    Edge (Vec3_t & X0, Vec3_t & X1); ///< Xi: endpoints of edge

    // Methods
    void UpdatedL  ();  

    // Data
    Vec3_t * X0; ///< Left endpoint
    Vec3_t * X1; ///< Right endpoint
    Vec3_t   dL; ///< Delta(X) = X1 - X0. difference Vector
    double Dmax; ///< Maximun length from edge centre
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Edge::Edge (Vec3_t * TheX0, Vec3_t * TheX1)
    : X0(TheX0), X1(TheX1), dL(*X1-*X0)
{
    Dmax = 0.5*norm(dL);
}

inline Edge::Edge (Vec3_t & TheX0, Vec3_t & TheX1)
    : X0(&TheX0), X1(&TheX1), dL(*X1-*X0)
{
    Dmax = 0.5*norm(dL);
}

inline void Edge::UpdatedL()
{
    dL = (*X1) - (*X0);
}

}
#endif // MECHSYS_DEM_EDGE_H
