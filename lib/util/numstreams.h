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

#ifndef MECHSYS_NUMSTREAMS_H
#define MECHSYS_NUMSTREAMS_H

// STL
#include <iostream>
#include <iomanip>

namespace Util
{

/** Number format via STL streams. */
struct NumStream
{
    bool BoolAlpha;  ///< Format output as a boolean?
    bool Integer;    ///< Format output as an integer?
    bool Scientific; ///< Format output as a scientific number?
    int  Width;      ///< Width of the output
    int  Precision;  ///< Precision for floating point numbers
    bool Reset;      ///< Reset ?
};

//                    bool  integ  scien   w   p  reset
NumStream _a     = {  true, false, false,  0,  0, false }; ///< Boolean
NumStream _2     = { false,  true, false,  2,  0, false }; ///< Integer
NumStream _3     = { false,  true, false,  3,  0, false }; ///< Integer
NumStream _4     = { false,  true, false,  4,  0, false }; ///< Integer
NumStream _6     = { false,  true, false,  6,  0, false }; ///< Integer
NumStream _8     = { false,  true, false,  8,  0, false }; ///< Integer
NumStream _3s    = { false, false,  true,  0,  3, false }; ///< Scientific
NumStream _6s    = { false, false,  true,  0,  6, false }; ///< Scientific
NumStream _8s    = { false, false,  true,  0,  8, false }; ///< Scientific
NumStream _15s   = { false, false,  true,  0, 15, false }; ///< Scientific
NumStream _3_2   = { false, false, false,  3,  2, false }; ///< General
NumStream _4_2   = { false, false, false,  4,  2, false }; ///< General
NumStream _6_2   = { false, false, false,  6,  3, false }; ///< General
NumStream _6_3   = { false, false, false,  6,  3, false }; ///< General
NumStream _6_4   = { false, false, false,  6,  4, false }; ///< General
NumStream _6_6   = { false, false, false,  6,  6, false }; ///< General
NumStream _8_0   = { false, false, false,  8,  0, false }; ///< General
NumStream _8_2   = { false, false, false,  8,  2, false }; ///< General
NumStream _8_3   = { false, false, false,  8,  3, false }; ///< General
NumStream _8_4   = { false, false, false,  8,  4, false }; ///< General
NumStream _8_6   = { false, false, false,  8,  6, false }; ///< General
NumStream _10_2  = { false, false, false, 10,  2, false }; ///< General
NumStream _10_3  = { false, false, false, 10,  3, false }; ///< General
NumStream _10_4  = { false, false, false, 10,  4, false }; ///< General
NumStream _10_6  = { false, false, false, 10,  6, false }; ///< General
NumStream _12_4  = { false, false, false, 12,  4, false }; ///< General
NumStream _12_6  = { false, false, false, 12,  6, false }; ///< General
NumStream _13_2  = { false, false, false, 13,  2, false }; ///< General
NumStream _13_6  = { false, false, false, 13,  6, false }; ///< General
NumStream _14_6  = { false, false, false, 14,  6, false }; ///< General
NumStream _15_2  = { false, false, false, 15,  2, false }; ///< General
NumStream _15_6  = { false, false, false, 15,  6, false }; ///< General
NumStream _20_15 = { false, false, false, 20, 15, false }; ///< General
NumStream _reset = { false, false, false,  0,  0, true  }; ///< Reset

/** Format the output. */
std::ostream & operator<< (std::ostream & os, NumStream const & NS)
{
    if (NS.Reset)
    {
        os << std::resetiosflags(std::ios_base::scientific) << std::resetiosflags(std::ios_base::fixed) << std::resetiosflags(std::ios_base::boolalpha);
    }
    else
    {
        if      (NS.BoolAlpha)  { os<<" "<<std::setw(6)<<std::boolalpha; }
        else if (NS.Integer)    { os<<" "<<std::setw(NS.Width); }
        else if (NS.Scientific) { os<<" "<<std::setw(NS.Precision+9)<<std::scientific<<std::setprecision(NS.Precision); } // add 9
        else                    { os<<" "<<std::setw(NS.Width+1)<<std::fixed<<std::setprecision(NS.Precision); } // add 1 (for sign)
    }
    return os;
}

}; // namespace Util

#endif // MECHSYS_NUMSTREAMS_H
