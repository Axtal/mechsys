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

#ifndef MECHSYS_SPHERES_H
#define MECHSYS_SPHERES_H

// VTK
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkLODActor.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkTextActor3D.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/colors.h>
#include <mechsys/util/string.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class Spheres
{
public:
    // Constructor & Destructor
     Spheres () : Ids(NULL) { _create();  SetResolution (); }
    ~Spheres ();

    // Alternative constructor
    Spheres (Array<Vec3_t> const & X, Array<double> const & R, int ThetaRes=20, int PhiRes=20);
    Spheres (Table const & XYZR, int ThetaRes=20, int PhiRes=20);
    Spheres (Array<double> const & X, Array<double> const & Y, Array<double> const & Z, Array<double> const & R, int ThetaRes=20, int PhiRes=20);

    // Set methods
    Spheres & SetSpheres    (Array<Vec3_t> const & X, Array<double> const * R=NULL);
    Spheres & SetSpheres    (Table const & XYZR);
    Spheres & SetSpheres    (Array<double> const & X, Array<double> const & Y, Array<double> const & Z, Array<double> const * R=NULL);
    Spheres & SetCenter     (int Id, Vec3_t const & X);
    Spheres & SetRadius     (int Id, double R);
    Spheres & SetResolution (int ThetaRes=20, int PhiRes=20);
    Spheres & SetColor      (char const * Name="red", double Opacity=0.8);

    // Alternative set methods
    void SetNumSpheres (int Num)                                        { _points->SetNumberOfPoints(Num);          _scalars->SetNumberOfTuples(Num); }
    void InsertSphere  (int Id, double x, double y, double z, double r) { _points->InsertPoint(Id,x,y,z);           _scalars->InsertTuple1(Id,2.0*r); }
    void InsertSphere  (int Id, Vec3_t const & X, double r)             { _points->InsertPoint(Id,X(0),X(1),X(2));  _scalars->InsertTuple1(Id,2.0*r); }
    void SetSphere     (int Id, Vec3_t const & X, double r)             { _points->SetPoint(Id,X(0),X(1),X(2));     _scalars->SetTuple1(Id,2.0*r);    }
    void Modified      ()                                               { _spheres->Modified(); }

    // Methods
    void ShowIds (double OriX=90, double OriY=90, double OriZ=45, double Scale=0.003, int SizePt=14, bool Shadow=true, char const * Color="black");
    void AddTo   (VTK::Win & win, bool RstCam=true);
    void DelFrom (VTK::Win & win);

    int const * Ids; // to be set externally: size = num points/spheres

private:
    vtkPoints              * _points;
    vtkDoubleArray         * _scalars;
    vtkSphereSource        * _sphere;
    vtkGlyph3D             * _spheres;
    vtkPolyDataMapper      * _spheres_mapper;
    vtkLODActor            * _spheres_actor;
    vtkLookupTable         * _ltable;
    Array<vtkTextActor3D*>   _text;
    void _create ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Spheres::Spheres (Array<Vec3_t> const & X, Array<double> const & R, int ThetaRes, int PhiRes)
    : Ids(NULL)
{
    _create       ();
    SetSpheres    (X, &R);
    SetResolution (ThetaRes, PhiRes);
}

inline Spheres::Spheres (Table const & XYZR, int ThetaRes, int PhiRes)
    : Ids(NULL)
{
    _create       ();
    SetSpheres    (XYZR);
    SetResolution (ThetaRes, PhiRes);
}

inline Spheres::Spheres (Array<double> const & X, Array<double> const & Y, Array<double> const & Z, Array<double> const & R, int ThetaRes, int PhiRes)
    : Ids(NULL)
{
    _create       ();
    SetSpheres    (X, Y, Z, &R);
    SetResolution (ThetaRes, PhiRes);
}

inline Spheres::~Spheres ()
{
    _points         -> Delete();
    _scalars        -> Delete();
    _sphere         -> Delete();
    _spheres        -> Delete();
    _spheres_mapper -> Delete();
    _spheres_actor  -> Delete();
    _ltable         -> Delete();
    for (size_t i=0; i<_text.Size(); ++i) _text[i] -> Delete();
}

inline Spheres & Spheres::SetSpheres (Array<Vec3_t> const & X, Array<double> const * R)
{
    _points  -> SetNumberOfPoints (X.Size());
    _scalars -> SetNumberOfTuples (X.Size());
    for (size_t i=0; i<X.Size(); ++i)
    {
        _points -> InsertPoint (i, X[i](0), X[i](1), X[i](2));
        if (R==NULL) _scalars -> InsertTuple1 (i, 1.0);
        else         _scalars -> InsertTuple1 (i, 2.0*(*R)[i]);
    }
    return (*this);
}

inline Spheres & Spheres::SetSpheres (Table const & XYZR)
{
    Array<double> const & X = XYZR("X");
    Array<double> const & Y = XYZR("Y");
    Array<double> const & Z = XYZR("Z");
    Array<double> const & R = XYZR("R");
    return SetSpheres (X, Y, Z, &R);
}

inline Spheres & Spheres::SetSpheres (Array<double> const & X, Array<double> const & Y, Array<double> const & Z, Array<double> const * R)
{
    if (Y.Size()!=X.Size()) throw new Fatal("Spheres::SetSpheres: X and Y arrays must have the same size");
    if (Z.Size()!=X.Size()) throw new Fatal("Spheres::SetSpheres: X and Z arrays must have the same size");
    _points  -> SetNumberOfPoints (X.Size());
    _scalars -> SetNumberOfTuples (X.Size());
    for (size_t i=0; i<X.Size(); ++i)
    {
        _points -> InsertPoint (i, X[i], Y[i], Z[i]);
        if (R==NULL) _scalars -> InsertTuple1 (i, 1.0);
        else         _scalars -> InsertTuple1 (i, 2.0*(*R)[i]);
    }
    return (*this);
}

inline Spheres & Spheres::SetCenter (int Id, Vec3_t const & X)
{
    _points -> SetPoint (Id, X(0), X(1), X(2));
    return (*this);
}

inline Spheres & Spheres::SetRadius (int Id, double R)
{
    _scalars -> SetTuple1 (Id, R);
    return (*this);
}

inline Spheres & Spheres::SetResolution (int ThetaRes, int PhiRes)
{
    _sphere -> SetThetaResolution (ThetaRes);
    _sphere -> SetPhiResolution   (PhiRes);
    return (*this);
}

inline Spheres & Spheres::SetColor (char const * Name, double Opacity)
{
    Vec3_t c = Colors::Get(Name);
    _ltable->SetNumberOfColors (2);
    _ltable->Build             ();
    _ltable->SetTableValue     (0, c(0), c(1), c(2));
    _ltable->SetTableValue     (1, c(0), c(1), c(2));
    _spheres_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline void Spheres::ShowIds (double OriX, double OriY, double OriZ, double Scale, int SizePt, bool Shadow, char const * Color)
{
    Vec3_t c(Colors::Get(Color));
    for (size_t i=0; i<_text.Size(); ++i) _text[i] -> Delete();
    String buf;
    _text.Resize (_points->GetNumberOfPoints());
    for (int i=0; i<_points->GetNumberOfPoints(); ++i)
    {
        if (Ids==NULL) buf.Printf ("%d",i);
        else           buf.Printf ("%d",Ids[i]);
        _text[i] = vtkTextActor3D                   ::New();
        _text[i] -> SetInput                        (buf.CStr());
        _text[i] -> SetPosition                     (_points->GetPoint(i));
        _text[i] -> SetOrientation                  (OriX, OriY, OriZ);
        _text[i] -> SetScale                        (Scale);
        _text[i] -> GetTextProperty()-> SetFontSize (SizePt);
        _text[i] -> GetTextProperty()-> SetShadow   (Shadow);
        _text[i] -> GetTextProperty()-> SetColor    (c.data());
    }
}

inline void Spheres::AddTo (VTK::Win & win, bool RstCam)
{
    win.AddActor (_spheres_actor, RstCam); 
    for (size_t i=0; i<_text.Size(); ++i) win.AddActor (reinterpret_cast<vtkActor*>(_text[i]), RstCam);
}

inline void Spheres::DelFrom (VTK::Win & win)
{
    win.DelActor (_spheres_actor); 
    for (size_t i=0; i<_text.Size(); ++i) win.DelActor (reinterpret_cast<vtkActor*>(_text[i]));
}

inline void Spheres::_create ()
{
    // points and scalars
    _points  = vtkPoints      ::New();
    _scalars = vtkDoubleArray ::New();
    _scalars -> SetNumberOfComponents (1);

    // polydata
    vtkPolyData * polydata = vtkPolyData   ::New();
    polydata -> SetPoints                  (_points);
    polydata -> GetPointData()->SetScalars (_scalars);

    // spheres
    _sphere         = vtkSphereSource              ::New();
    _spheres        = vtkGlyph3D                   ::New();
    _spheres_mapper = vtkPolyDataMapper            ::New();
    _spheres_actor  = vtkLODActor                  ::New();
    _ltable         = vtkLookupTable               ::New();
    _spheres        -> SetInput                    (polydata);
    _spheres        -> SetSource                   (_sphere->GetOutput());
    _spheres        -> SetScaleModeToScaleByScalar ();
    _spheres        -> SetColorModeToColorByScalar ();
    _spheres        -> SetScaleFactor              (1.0);
    _spheres_mapper -> SetInputConnection          (_spheres->GetOutputPort());
    _spheres_mapper -> SetLookupTable              (_ltable);
    _spheres_actor  -> SetMapper                   (_spheres_mapper);
    SetColor ();
}

}; // namespace VTK

#endif // MECHSYS_SPHERES_H
