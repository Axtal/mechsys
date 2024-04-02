/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2017 Maziar Gholami                                    *
 * Copyright (C) 2020 Mario Trujillo                                    *
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

#ifndef SPH_BOUNDARY_CONDITION_H
#define SPH_BOUNDARY_CONDITION_H

namespace SPH 
{

class Boundary
{
public:
    // Data
    double	inDensity;	///< Apply a certain density to inflow particles
    double	outDensity;	///< Apply a certain density to outflow particles
    double	allDensity;	///< Apply a certain density to outflow particles

    Vec3_t 	inv;		///< Apply a certain velocity to inflow particles
    Vec3_t 	outv;		///< Apply a certain velocity to outflow particles
    Vec3_t	allv;		///< Apply a certain velocity to all particle

    bool 	Periodic[3];	///< Considering periodic in all directions => 0=X, 1=Y, 2=Z

    int 	InOutFlow;	///< Considering inflow in all directions  by adding and deleting particles=> [0]=X, [1]=Y, [2]=Z and 0=none, 1=-
    double	InFlowLoc1;
    double	InFlowLoc2;
    double	InFlowLoc3;
    double	OutFlowLoc;
    double	cellfac;
    int		inoutcounter;
    bool	MassConservation;

    Array <int>	OutPart;
    Array <int>	InPart;

    Boundary();
};

inline Boundary::Boundary()
{
    allv	= 0.0;
    inv		= 0.0;
    outv	= 0.0;
    Periodic[0]=Periodic[1]=Periodic[2] = false;
    inDensity = 0.0;
    outDensity = 0.0;
    allDensity = 0.0;
    InOutFlow = 0;
    InFlowLoc1 = 0.0;
    InFlowLoc2 = 0.0;
    InFlowLoc3 = 0.0;
    OutFlowLoc = 0.0;
    cellfac = 3.9;
    inoutcounter = 0;
    MassConservation = false;
}

}; // namespace SPH

#endif // SPH_BOUNDARY_CONDITION_H
