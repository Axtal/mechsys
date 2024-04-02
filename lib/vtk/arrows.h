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

#ifndef MECHSYS_ARROWS_H
#define MECHSYS_ARROWS_H

// VTK
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkAppendPolyData.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkLODActor.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/sgrid.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/colors.h>
#include <mechsys/util/string.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class Arrows
{
public:
    // Constructor & Destructor
     Arrows () { _create();  SetGeometry();  SetScale(); }
    ~Arrows ();

    // Alternative constructor
    Arrows (VTK::SGrid & G, double ConPct=0.1, int Resolution=12);

    // Set methods
    Arrows & SetArrows   (VTK::SGrid & G);
    Arrows & SetGeometry (double ConPct=0.1, int Resolution=12);
    Arrows & SetColor    (char const * Name="blue", double Opacity=1.0);
    Arrows & SetScale    (double Factor=1.0) { _arrows->SetScaleFactor(Factor);  return (*this); }

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor (_arrows_actor); }

private:
    vtkConeSource     * _cone;
    vtkCylinderSource * _cylin;
    vtkGlyph3D        * _arrows;
    vtkPolyDataMapper * _arrows_mapper;
    vtkLODActor       * _arrows_actor;
    vtkLookupTable    * _ltable;
    void _create ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Arrows::Arrows (VTK::SGrid & G, double ConPct, int Resolution)
{
    _create     ();
    SetArrows   (G);
    SetGeometry (ConPct, Resolution);
    SetScale    ();
}

inline Arrows::~Arrows ()
{
    _cone          -> Delete();
    _cylin         -> Delete();
    _arrows        -> Delete();
    _arrows_mapper -> Delete();
    _arrows_actor  -> Delete();
    _ltable        -> Delete();
}

inline Arrows & Arrows::SetArrows (VTK::SGrid & G)
{
    vtkPolyData * polydata = vtkPolyData   ::New();
    polydata -> SetPoints                  (G.GetPoints());
    polydata -> GetPointData()->SetScalars (G.GetScalars());
    polydata -> GetPointData()->SetVectors (G.GetVectors());
    _arrows  -> SetInput                   (polydata);
    polydata -> Delete                     ();
    return (*this);
}

inline Arrows & Arrows::SetGeometry (double ConPct, int Resolution)
{
    double con_rad = 0.03;
    double cyl_rad = 0.015;
    double tot_len = 1.0;                   // length of arrow (glyph)
    double con_len = ConPct*tot_len;        // cone length/height
    double cyl_len = tot_len-con_len;       // cylinder length/height
    double con_cen = (tot_len-con_len)/2.0; // cone center
    double cyl_cen = -con_len/2.0;          // cylinder center

    _cone  -> SetCenter     (0.0, con_cen+tot_len/2.0, 0.0);
    _cone  -> SetHeight     (con_len);
    _cylin -> SetCenter     (0.0, cyl_cen+tot_len/2.0, 0.0);
    _cylin -> SetHeight     (cyl_len);
    _cone  -> SetRadius     (con_rad);
    _cylin -> SetRadius     (cyl_rad);
    _cone  -> SetResolution (Resolution);
    _cylin -> SetResolution (Resolution);

    return (*this);
}

inline Arrows & Arrows::SetColor (char const * Name, double Opacity)
{
    Vec3_t c = Colors::Get(Name);
    _ltable->SetNumberOfColors (2);
    _ltable->Build             ();
    _ltable->SetTableValue     (0, c(0), c(1), c(2));
    _ltable->SetTableValue     (1, c(0), c(1), c(2));
    _arrows_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline void Arrows::_create ()
{
    // create arrow
    _cone                     = vtkConeSource     ::New();
    _cylin                    = vtkCylinderSource ::New();
    vtkAppendPolyData * arrow = vtkAppendPolyData ::New();
	_cone -> SetDirection (0,1,0); // because cylinder direction is along y axis
	arrow -> AddInput     (_cone->GetOutput());
	arrow -> AddInput     (_cylin->GetOutput());

	// rotate around z axis because glyph3D needs the arrow along the x-axis
	vtkTransform       * rotate = vtkTransform       ::New();
	vtkTransformFilter * transf = vtkTransformFilter ::New();
	rotate -> RotateZ      (-90.0);
	transf -> SetTransform (rotate);
	transf -> SetInput     (arrow->GetOutput());

    // glyph
    _arrows        = vtkGlyph3D                   ::New();
    _arrows_mapper = vtkPolyDataMapper            ::New();
    _arrows_actor  = vtkLODActor                  ::New();
    _ltable        = vtkLookupTable               ::New();
    _arrows        -> SetSource                   (transf->GetPolyDataOutput());
    _arrows        -> SetScaleModeToScaleByVector ();
    _arrows        -> SetColorModeToColorByVector ();
    _arrows_mapper -> SetInputConnection          (_arrows->GetOutputPort());
    _arrows_mapper -> SetLookupTable              (_ltable);
    _arrows_actor  -> SetMapper                   (_arrows_mapper);

    // clean up
    arrow  -> Delete ();
    rotate -> Delete ();
    transf -> Delete ();
}

}; // namespace VTK

#endif // MECHSYS_ARROWS_H
