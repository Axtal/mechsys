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

// VTK
#include <vtkLightKit.h>

// C++ STL
#include <fstream>

// MechSys
#include <mechsys/dem/graph.h>
#include <mechsys/dem/domain.h>
#include <mechsys/dem/visualise.h>
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/axes.h>

int main(int argc, char **argv) try
{
    bool   readmsh = true;
    bool   showvtk = false;
    double R       = 0.5;
    if (argc>1) readmsh = atoi(argv[1]);
    if (argc>2) showvtk = atoi(argv[2]);
    if (argc>3) R       = atof(argv[3]);

    // DEM domain
    DEM::Domain d;

    // read mesh of particle
    if (readmsh)
    {
        DEM::Particle * p = new DEM::Particle;
        //p->ConstructFromJson(-2, "bucket.msh", R, 3.0, 10.0);
        p->ConstructFromJson(-2, "Rawbucket.msh", R, 3.0, 10.0);
        d.Particles.Push(p);
    }
    else
    {
        d.Particles.Push(new DEM::Particle(-2, "Rawbucket", R, 3.0, 10.0));
    }

    // debugging
    if (true)
    {
        std::cout << "EdgeCon =\n";
        for (size_t i=0; i<d.Particles[0]->EdgeCon.Size(); ++i)
        {
            std::cout << " " << d.Particles[0]->EdgeCon[i] << std::endl;
        }
        std::cout << "\nFaceCon =\n";
        for (size_t i=0; i<d.Particles[0]->FaceCon.Size(); ++i)
        {
            std::cout << " " << d.Particles[0]->FaceCon[i] << std::endl;
        }
    }

    // save output file for Blender
    std::ofstream of1("viewbucket.py", std::ios::out);
    DEM::BPYHeader(of1);
    d.Particles[0]->Draw(of1, "Red", true);
    of1.close();

    // save output file for POV
    std::ofstream of2("viewbucket.pov", std::ios::out);
    DEM::POVHeader(of2);
    d.Particles[0]->Draw(of2, "Red", false);
    of2.close();

    // visualisation
    if (showvtk) {
        bool showvert = true;
        bool showedge = true;
        DEM::Visualise vis(d, /*parts*/Array<int>(-2,true), /*walls*/Array<int>(-1,true), showvert, showedge);
        VTK::Axes axe(1.0, /*hline*/false, /*reverse*/false, /*full*/true);
        VTK::Win win;
        /*
        vis.PartFaces.SetMaterial (0.2, 0.1, 1.0, 200.0);
        vis.PartFaces.SetColor    ("cobalt", 1.0);
        double elevation = 60.0;
        double azimuth   = -45.0;./o./
        vtkLightKit *lgt = NULL;
        lgt->AddLightsToRenderer (win.GetRen());
        lgt->SetKeyLightAngle (elevation, azimuth);
        lgt->SetKeyToHeadRatio (0.01);
        lgt->SetKeyToFillRatio (1.2);
        */
		vis.AddTo(win);
		axe.AddTo(win);
		win.Render();
		win.Show();
	}

    return 0;
}
MECHSYS_CATCH
