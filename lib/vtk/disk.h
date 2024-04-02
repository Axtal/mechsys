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

#ifndef MECHSYS_DISK_H
#define MECHSYS_DISK_H

// Std Lib
#include <cmath>

// VTK
#include <vtkDiskSource.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class Disk
{
public:
    // Constructor & Destructor
     Disk ();
    ~Disk ();

    // Alternative constructor
    Disk (Vec3_t const & X0, Vec3_t const & X1, double Rin, double Rout, bool Capping=true, int RRes=20, int CRes=30);

    // Set methods
    Disk & SetRadiusIn   (double Radius);
    Disk & SetRadiusOut  (double Radius);
    Disk & SetPoints     (Vec3_t const & X0, Vec3_t const & X1);
    Disk & SetResolution (int RRes=20, int CRes=30);
    Disk & SetColor      (char const * Name="rose_madder", double Opacity=1.0);
    Disk & SetWire       () { _disk_actor->GetProperty()->SetRepresentationToWireframe();  return (*this); }

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor(_disk_actor); }

private:
    vtkDiskSource      * _disk;
    vtkTransformFilter * _transform;
    vtkLinearExtrusionFilter * _extrusion;
    vtkPolyDataMapper  * _disk_mapper;
    vtkActor           * _disk_actor;
    void _create (bool Cap=true);
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Disk::Disk ()
{
    _create       ();
    SetResolution ();
    SetColor      ();
}

inline Disk::Disk (Vec3_t const & X0, Vec3_t const & X1, double Rin, double Rout, bool Cap, int RRes, int CRes)
{
    _create       (Cap);
    SetRadiusIn   (Rin);
    SetRadiusOut  (Rout);
    SetPoints     (X0, X1);
    SetResolution (RRes, CRes);
    SetColor      ();
}

inline Disk::~Disk ()
{
    _disk        -> Delete();
    _transform   -> Delete();
    _disk_mapper -> Delete();
    _disk_actor  -> Delete();
}

inline Disk & Disk::SetRadiusIn (double Rin)
{
    _disk->SetInnerRadius (Rin);
    return (*this);
}

inline Disk & Disk::SetRadiusOut (double Rout)
{
    _disk->SetOuterRadius (Rout);
    return (*this);
}

inline Disk & Disk::SetPoints (Vec3_t const & X0, Vec3_t const & X1)
{
    // update length
    Vec3_t V(X1-X0);

    // extrude along z
    double len = Norm(V);
    _extrusion->SetVector (0.0, 0.0, len);

    // translate
    Vec3_t cen(X0); // center of disk
    vtkTransform * affine = vtkTransform ::New();
    affine->Translate (cen(0), cen(1), cen(2));

    // rotate
    Vec3_t vy(0.0, 0.0, 1.0); // direction of extruded source
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

inline Disk & Disk::SetResolution (int RRes, int CRes)
{
    _disk->SetRadialResolution          (RRes);
    _disk->SetCircumferentialResolution (CRes);
    return (*this);
}

inline Disk & Disk::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _disk_actor->GetProperty()->SetColor   (c.data());
    _disk_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline void Disk::_create (bool Cap)
{
    _disk         = vtkDiskSource      ::New();
    _extrusion    = vtkLinearExtrusionFilter::New();
    _transform    = vtkTransformFilter ::New();
    _disk_mapper  = vtkPolyDataMapper  ::New();
    _disk_actor   = vtkActor           ::New();
    _extrusion   -> SetCapping        (Cap);
    _extrusion   -> SetInput          (_disk->GetOutput());
    _transform   -> SetInput          (_extrusion->GetOutput());
    _disk_mapper -> SetInput          (_transform->GetPolyDataOutput());
    _disk_actor  -> SetMapper         (_disk_mapper);
}

}; // namespace VTK

#endif // MECHSYS_DISK_H
