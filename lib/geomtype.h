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

#ifndef MECHSYS_GEOMTYPE_H
#define MECHSYS_GEOMTYPE_H

// Std Lib
#include <cstring> // for strcmp

// MechSys
#include <mechsys/util/string.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

enum GeomType { fra_t,   ///< Frame. Ex.: Truss, Beam, Spring
                psa_t,   ///< Plane-strain
                pse_t,   ///< Plane-stress
                axs_t,   ///< Axis-symmetric
                d2d_t,   ///< 2D
                d3d_t }; ///< 3D (solid)

String GTypeToStr (GeomType Type)
{
    String info;
         if (Type==fra_t) info = "fra";
    else if (Type==psa_t) info = "psa";
    else if (Type==pse_t) info = "pse";
    else if (Type==axs_t) info = "axs";
    else if (Type==d2d_t) info = "d2d";
    else if (Type==d3d_t) info = "d3d";
    return info;
}

GeomType StrToGType (const char * Str)
{
    GeomType gty;
         if (strcmp(Str,"fra")==0) gty = fra_t;
    else if (strcmp(Str,"psa")==0) gty = psa_t;
    else if (strcmp(Str,"pse")==0) gty = pse_t;
    else if (strcmp(Str,"axs")==0) gty = axs_t;
    else if (strcmp(Str,"d2d")==0) gty = d2d_t;
    else if (strcmp(Str,"d3d")==0) gty = d3d_t;
    else throw new Fatal("StrToGType: GeomType key %s is invalid",Str);
    return gty;
}

GeomType SDPairToGType (SDPair const & D, const char * Default)
{
    GeomType gty;
         if (D.HasKey("fra")) gty = fra_t;
    else if (D.HasKey("psa")) gty = psa_t;
    else if (D.HasKey("pse")) gty = pse_t;
    else if (D.HasKey("axs")) gty = axs_t;
    else if (D.HasKey("d2d")) gty = d2d_t;
    else if (D.HasKey("d3d")) gty = d3d_t;
    else                      gty = StrToGType(Default);
    return gty;
}

#endif // MECHSYS_GEOMTYPE_H
