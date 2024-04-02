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

#ifndef MECHSYS_CYLINDER_H
#define MECHSYS_CYLINDER_H

/** @file vtk/cylinder.h .*/

// Std Lib
#include <cmath>

// VTK
#include <vtkCylinderSource.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class Cylinder
{
public:
    // Constructor & Destructor
     Cylinder ();
    ~Cylinder ();

    // Alternative constructor
    Cylinder (Vec3_t const & X0, Vec3_t const & X1, double Radius, bool Capping=true, int Resolution=20);

    // Set methods
    Cylinder & SetRadius     (double Radius);
    Cylinder & SetPoints     (Vec3_t const & X0, Vec3_t const & X1);
    Cylinder & SetResolution (int Res=20);
    Cylinder & SetColor      (char const * Name="olive", double Opacity=1.0);
    Cylinder & SetWire       () { _cylin_actor->GetProperty()->SetRepresentationToWireframe();  return (*this); }

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor(_cylin_actor); }

private:
    vtkCylinderSource  * _cylin;
    vtkTransformFilter * _transform;
    vtkPolyDataMapper  * _cylin_mapper;
    vtkActor           * _cylin_actor;
    void _create (bool Cap=true);
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Cylinder::Cylinder ()
{
    _create       ();
    SetResolution ();
    SetColor      ();
}

inline Cylinder::Cylinder (Vec3_t const & X0, Vec3_t const & X1, double R, bool Cap, int Res)
{
    _create       (Cap);
    SetRadius     (R);
    SetPoints     (X0, X1);
    SetResolution (Res);
    SetColor      ();
}

inline Cylinder::~Cylinder ()
{
    _cylin        -> Delete();
    _transform    -> Delete();
    _cylin_mapper -> Delete();
    _cylin_actor  -> Delete();
}

inline Cylinder & Cylinder::SetRadius (double R)
{
    _cylin->SetRadius (R);
    return (*this);
}

inline Cylinder & Cylinder::SetPoints (Vec3_t const & X0, Vec3_t const & X1)
{
    // update length
    Vec3_t V(X1-X0);
    _cylin->SetHeight (Norm(V));

    // translate
    Vec3_t cen(X0+0.5*V); // center of cylin
    vtkTransform * affine = vtkTransform ::New();
    affine->Translate (cen(0), cen(1), cen(2));

    // rotate
    Vec3_t vy(0.0, 1.0, 0.0); // direction of cylinder source
    double angle = (180.0/Util::PI)*acos(dot(vy,V)/norm(V)); // angle of rotation
    if (angle>0.0)
    {
        Vec3_t axis = cross(vy, V); // axis of rotation
        if (norm(axis)>0.0)         // not parallel
        {
            affine->RotateWXYZ (angle, axis.data());
        }
        else // parallel and oposite (alpha=180)
        {
            affine->RotateWXYZ (angle, 0.0, 0.0, 1.0); // use z-direction for mirroring
        }
    }

    // tranform
    _transform->SetTransform (affine);

    // clean up
    affine->Delete ();
    return (*this);
}

inline Cylinder & Cylinder::SetResolution (int Res)
{
    _cylin->SetResolution (Res);
    return (*this);
}

inline Cylinder & Cylinder::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _cylin_actor->GetProperty()->SetColor   (c.data());
    _cylin_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline void Cylinder::_create (bool Cap)
{
    _cylin        = vtkCylinderSource  ::New();
    _transform    = vtkTransformFilter ::New();
    _cylin_mapper = vtkPolyDataMapper  ::New();
    _cylin_actor  = vtkActor           ::New();
    _cylin        -> SetCapping        (Cap);
    _transform    -> SetInput          (_cylin->GetOutput());
    _cylin_mapper -> SetInput          (_transform->GetPolyDataOutput());
    _cylin_actor  -> SetMapper         (_cylin_mapper);
}

}; // namespace VTK

#endif // MECHSYS_CYLINDER_H
