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

#ifndef MECHSYS_ISOSURF_H
#define MECHSYS_ISOSURF_H

// VTK
#include <vtkMarchingContourFilter.h>
#include <vtkHedgeHog.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkLookupTable.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/sgrid.h>
#include <mechsys/vtk/win.h>
#include <mechsys/util/array.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class IsoSurf
{
public:
    // Constructor & Destructor
    IsoSurf (VTK::SGrid & G);

    // Destructor
    ~IsoSurf ();

    // Set methods
    void SetColors   (Array<String> const & Names, double Opacity=1.0);
    void SetColor    (char const * Name="blue",    double Opacity=1.0) { SetColors(Array<String>(Name,true), Opacity); }
    void SetValue    (double F=0.0)                                    { _isosurf->SetValue       (0,F);               }
    void GenValues   (int NSurfs, double fMin, double fMax)            { _isosurf->GenerateValues (NSurfs,fMin,fMax);  }
    void SetVecScale (double Factor=1.0)                               { _hedgehog->SetScaleFactor(Factor);            }
    void SetWire     ()                                                { _isosurf_actor->GetProperty()->SetRepresentationToWireframe(); }

    // Extra methods
    void SetMaterial (double Ambient, double Diffuse, double Specular, double SpecularPower);

    // Methods
    void AddTo (VTK::Win & win);

    // Data
    bool ShowIsoSurf;
    bool ShowVectors;

private:
    vtkMarchingContourFilter * _isosurf;
    vtkPolyDataMapper        * _isosurf_mapper;
    vtkActor                 * _isosurf_actor;
	vtkLookupTable           * _isosurf_lt;
    vtkHedgeHog              * _hedgehog;
    vtkPolyDataMapper        * _hedgehog_mapper;
    vtkActor                 * _hedgehog_actor;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline IsoSurf::IsoSurf (VTK::SGrid & G)
    : ShowIsoSurf (true),
      ShowVectors (true)
{
    // isosurf
    _isosurf        = vtkMarchingContourFilter ::New();
    _isosurf_mapper = vtkPolyDataMapper        ::New();
    _isosurf_actor  = vtkActor                 ::New();
	_isosurf_lt     = vtkLookupTable           ::New();
    _isosurf        -> SetInput                (G.GetGrid());
    _isosurf        -> ComputeNormalsOff       ();
    _isosurf        -> ComputeGradientsOff     ();
    _isosurf_mapper -> SetInputConnection      (_isosurf->GetOutputPort());
    _isosurf_mapper -> SetLookupTable          (_isosurf_lt);
    _isosurf_actor  -> SetMapper               (_isosurf_mapper);
    SetColor ();
    SetValue ();

    // hedgehog
    _hedgehog        = vtkHedgeHog         ::New();
    _hedgehog_mapper = vtkPolyDataMapper   ::New();
    _hedgehog_actor  = vtkActor            ::New();
    _hedgehog        -> SetInput           (G.GetGrid());
    _hedgehog_mapper -> SetInputConnection (_hedgehog->GetOutputPort());
    _hedgehog_actor  -> SetMapper          (_hedgehog_mapper);
    SetVecScale ();
}

inline IsoSurf::~IsoSurf ()
{
    _isosurf         -> Delete();
    _isosurf_mapper  -> Delete();
    _isosurf_actor   -> Delete();
    _isosurf_lt      -> Delete();
    _hedgehog        -> Delete();
    _hedgehog_mapper -> Delete();
    _hedgehog_actor  -> Delete();
}

inline void IsoSurf::SetColors (Array<String> const & Colors, double Opacity)
{
    _isosurf_lt->SetNumberOfColors (Colors.Size());
    _isosurf_lt->Build             ();
    for (size_t i=0; i<Colors.Size(); ++i)
    {
        Vec3_t c = Colors::Get(Colors[i]);
        _isosurf_lt->SetTableValue (i, c(0), c(1), c(2));
    }
	_isosurf_actor->GetProperty()->SetOpacity(Opacity);
}

inline void IsoSurf::AddTo (VTK::Win & win)
{
    if (ShowIsoSurf) win.AddActor (_isosurf_actor);
    if (ShowVectors) win.AddActor (_hedgehog_actor);
}

inline void IsoSurf::SetMaterial (double Ambient, double Diffuse, double Specular, double SpecularPower)
{
    _isosurf_actor->GetProperty()->SetAmbient       (Ambient);
    _isosurf_actor->GetProperty()->SetDiffuse       (Diffuse);
    _isosurf_actor->GetProperty()->SetSpecular      (Specular);
    _isosurf_actor->GetProperty()->SetSpecularPower (SpecularPower);
}

}; // namespace VTK

#endif // MECHSYS_ISOSURF_H
