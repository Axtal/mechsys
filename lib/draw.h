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

#ifndef MECHSYS_FEM_GRAPH_H
#define MECHSYS_FEM_GRAPH_H

// Std Lib
#include <iostream>

/** MatPlotLib. */
namespace MPL
{

inline void Header (std::ostream & os)
{
    // header
    os << "from numpy import array, sqrt\n";
    os << "from pylab import figure, text, show, axis, plot, grid, savefig\n";
    os << "from pylab import matplotlib as MPL\n\n";
    os << "PH = MPL.path.Path\n";
    os << "PP = MPL.patches\n";
    os << "PC = MPL.patches.PathPatch\n\n";

    // colors
    os << "# colors\n";
    os << "pink    = (250/255.0,204/255.0,228/255.0)\n";
    os << "dred    = (163/255.0,  0/255.0,  0/255.0)\n";
    os << "lblue   = (217/255.0,228/255.0,255/255.0)\n";
    os << "lgreen  = (100/255.0,241/255.0,193/255.0)\n";
    os << "dblue   = ( 45/255.0,  0/255.0,160/255.0)\n";
    os << "orange  = (241/255.0,125/255.0,  0/255.0)\n";
    os << "dpink   = (170/255.0, 94/255.0,137/255.0)\n";
    os << "lyellow = (234/255.0,228/255.0,179/255.0)\n\n";

    // figure, axis, and data
    os << "# figure, axis, and data\n";
    os << "fig = figure()\n";
    os << "ax  = fig.add_subplot(111)\n";
    os << "dat = []\n\n";
}

inline void AddPatch (std::ostream & os, char const * EdgeColor="dblue", char const * FaceColor="lblue")
{
    os << "if len(dat)>0:\n";
    os << "    cmd,vert = zip(*dat)\n";
    os << "    ph       = PH (vert, cmd)\n";
    os << "    pc       = PC (ph, edgecolor=" << EdgeColor << ", facecolor=" << FaceColor << ", linewidth=2)\n";
    os << "    ax.add_patch  (pc)\n\n";
}

inline void Show (std::ostream & os)
{
    os << "grid()\n";
    os << "axis('scaled')\n";
    os << "show()\n";
}

inline void SaveFig (char const * FileKey, std::ostream & os)
{
    String fn(FileKey);
    fn.append(".png");
    os << "grid()\n";
    os << "axis('scaled')\n";
    os << "savefig('"<< fn << "')\n";
}

}; // namespace MPL

#endif // MECHSYS_FEM_GRAPH_H
