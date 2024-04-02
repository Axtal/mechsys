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

#ifndef MECHSYS_PLANE_H
#define MECHSYS_PLANE_H

// Std Lib
#include <string>

// VTK
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class Plane
{
public:
    // Constructor and destructor
     Plane () { _create(); }
    ~Plane ();

    // Alternative constructors
    Plane (Vec3_t const & Cen, Vec3_t const & n) { _create();  SetCen(Cen);  SetNormal(n); } // centre, and normal
    Plane (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2);                      // origin, point1, point2
    Plane (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2, Vec3_t const & n);    // origin, point1, point2, normal

    // Set methods
    Plane & SetCen       (Vec3_t const & Cen) { _plane->SetCenter (Cen(0), Cen(1), Cen(2));  return (*this); }
    Plane & SetNormal    (Vec3_t const & n)   { _plane->SetNormal (  n(0),   n(1),   n(2));  return (*this); }
    Plane & SetTuple     (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2);
    Plane & SetTuple     (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2, Vec3_t const & n);
    Plane & SetColor     (char const * Name="light_salmon", double Opacity=0.5);
    Plane & SetWireColor (char const * Name="black");
    Plane & SetWireWidth (int Width=1.0);

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor(_plane_actor);  win.AddActor(_wire_actor); }

private:
    vtkPlaneSource    * _plane;
    vtkPolyDataMapper * _plane_mapper;
    vtkActor          * _plane_actor;
    vtkPolyDataMapper * _wire_mapper;
    vtkActor          * _wire_actor;
    void _create ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Plane::Plane (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2)
{
    _create  ();
    SetTuple (Ori, Pt1, Pt2);
}

inline Plane::Plane (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2, Vec3_t const & n)
{
    _create  ();
    SetTuple (Ori, Pt1, Pt2, n);
}

inline Plane::~Plane ()
{
    _plane        -> Delete();
    _plane_mapper -> Delete();
    _plane_actor  -> Delete();
    _wire_mapper  -> Delete();
    _wire_actor   -> Delete();
}

inline Plane & Plane::SetTuple (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2)
{
    _plane -> SetOrigin (Ori(0), Ori(1), Ori(2));
    _plane -> SetPoint1 (Pt1(0), Pt1(1), Pt1(2));
    _plane -> SetPoint2 (Pt2(0), Pt2(1), Pt2(2));
    return (*this);
}

inline Plane & Plane::SetTuple (Vec3_t const & Ori, Vec3_t const & Pt1, Vec3_t const & Pt2, Vec3_t const & n)
{
    _plane -> SetOrigin (Ori(0), Ori(1), Ori(2));
    _plane -> SetPoint1 (Pt1(0), Pt1(1), Pt1(2));
    _plane -> SetPoint2 (Pt2(0), Pt2(1), Pt2(2));
    _plane -> SetNormal (  n(0),   n(1),   n(2));
    return (*this);
}

inline Plane & Plane::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _plane_actor->GetProperty()->SetColor   (c.data());
    _plane_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline Plane & Plane::SetWireColor (char const * Name) 
{
    Vec3_t c(Colors::Get(Name));
    _wire_actor->GetProperty()->SetColor (c.data());
    return (*this); 
}

inline Plane & Plane::SetWireWidth (int Width)
{
    _wire_actor->GetProperty()->SetLineWidth(Width);
    return (*this); 
}

inline void Plane::_create ()
{
    // plane
    _plane        = vtkPlaneSource    ::New();
    _plane_mapper = vtkPolyDataMapper ::New();
    _plane_actor  = vtkActor          ::New();
    _plane_mapper -> SetInputConnection (_plane->GetOutputPort());
    _plane_actor  -> SetMapper          (_plane_mapper);

    // borders
    _wire_mapper = vtkPolyDataMapper    ::New();
    _wire_actor  = vtkActor             ::New();
    _wire_mapper -> SetInput            (_plane->GetOutput());
    _wire_mapper -> ScalarVisibilityOff ();
    _wire_actor  -> SetMapper           (_wire_mapper);
    _wire_actor  -> GetProperty         ()->SetRepresentationToWireframe();

    // set mapper
    _plane_mapper -> SetResolveCoincidentTopologyPolygonOffsetParameters (0,1);
    _plane_mapper -> SetResolveCoincidentTopologyToPolygonOffset         ();
    _wire_mapper  -> SetResolveCoincidentTopologyPolygonOffsetParameters (1,1);
    _wire_mapper  -> SetResolveCoincidentTopologyToPolygonOffset         ();

    // same color for inside and outside edges
    _wire_mapper -> ScalarVisibilityOff          ();
    _wire_actor  -> GetProperty() -> SetAmbient  (1.0);
    _wire_actor  -> GetProperty() -> SetDiffuse  (0.0);
    _wire_actor  -> GetProperty() -> SetSpecular (0.0);

    // set colors and wire width
    SetColor     ();
    SetWireColor ();
    SetWireWidth ();
}

}; // namespace VTK

#endif // MECHSYS_PLANE_H
