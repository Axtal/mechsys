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

// Std Lib
#include <cstdlib> // for srand, rand and RAND_MAX
#include <ctime>   // for time

// VTK
#include <vtkPolygon.h>

// MechSys
#include <mechsys/dem/visualise.h>
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/axes.h>

using Util::PI;
using namespace DEM;

struct Data
{
     Data () : Vis(NULL) {}
    ~Data () { if (Vis!=NULL) delete Vis; }
    void Init (DEM::Domain const & Dom)
    {
        Vis = new Visualise(Dom, /*parts*/Array<int>(-10,true), /*walls*/Array<int>());
        VTK::Axes axe(1.0, /*hline*/false, /*reverse*/false, /*full*/true);
        Vis->AddTo  (Win);
        axe .AddTo  (Win);
        Win .Render ();
        Win .Render ();
        IdxOut = 0;
    }
    Visualise * Vis;
    VTK::Win    Win;
    int         IdxOut;
};

int main(int argc, char **argv) try
{
    // input
    String filekey("xpack-voro");
    if (argc>1) filekey = argv[1];
    std::cout << "Usage:\n  " << argv[0] << " filekey [.hdf5]" << std::endl;

    // user data and domain
    Data dat;
    DEM::Domain dom(&dat);

    // read particles
    dom.Load (filekey.CStr());

    // show
    dat.Init (dom);
    dat.Win.Show ();

    // end
    return 0;
}
MECHSYS_CATCH
