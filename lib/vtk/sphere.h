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

#ifndef MECHSYS_SPHERE_H
#define MECHSYS_SPHERE_H

// VTK
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class Sphere
{
public:
    // Constructor & Destructor
     Sphere ();
    ~Sphere ();

    // Alternative constructor
    Sphere (Vec3_t const & X, double R, int ThetaRes=20, int PhiRes=20);

    // Set methods
    Sphere & SetCenter     (Vec3_t const & X);
    Sphere & SetRadius     (double R);
    Sphere & SetResolution (int ThetaRes=20, int PhiRes=20);
    Sphere & SetColor      (char const * Name="brown", double Opacity=0.8);

    // Methods
    void AddTo (VTK::Win & win, bool RstCam=true) { win.AddActor(_sphere_actor, RstCam); }

private:
    vtkSphereSource   * _sphere;
    vtkPolyDataMapper * _sphere_mapper;
    vtkActor          * _sphere_actor;
    void _create ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Sphere::Sphere ()
{
    _create       ();
    SetResolution ();
    SetColor      ();
}

inline Sphere::Sphere (Vec3_t const & X, double R, int ThetaRes, int PhiRes)
{
    _create       ();
    SetCenter     (X);
    SetRadius     (R);
    SetResolution (ThetaRes, PhiRes);
    SetColor      ();
}

inline Sphere::~Sphere ()
{
    _sphere        -> Delete();
    _sphere_mapper -> Delete();
    _sphere_actor  -> Delete();
}

inline Sphere & Sphere::SetCenter (Vec3_t const & X)
{
    _sphere -> SetCenter (X(0), X(1), X(2));
    return (*this);
}

inline Sphere & Sphere::SetRadius (double R)
{
    _sphere -> SetRadius (R);
    return (*this);
}

inline Sphere & Sphere::SetResolution (int ThetaRes, int PhiRes)
{
    _sphere -> SetThetaResolution (ThetaRes);
    _sphere -> SetPhiResolution   (PhiRes);
    return (*this);
}

inline Sphere & Sphere::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _sphere_actor->GetProperty()->SetColor   (c.data());
    _sphere_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline void Sphere::_create ()
{
    _sphere        = vtkSphereSource     ::New();
    _sphere_mapper = vtkPolyDataMapper   ::New();
    _sphere_actor  = vtkActor            ::New();
    _sphere_mapper -> SetInputConnection (_sphere->GetOutputPort());
    _sphere_actor  -> SetMapper          (_sphere_mapper);
}

}; // namespace VTK

#endif // MECHSYS_SPHERE_H
