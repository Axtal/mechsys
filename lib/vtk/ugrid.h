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

#ifndef MECHSYS_POLYGON_H
#define MECHSYS_POLYGON_H

// Std Lib
#include <string>

// VTK
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

class UGrid
{
public:
    // Constructors
    UGrid ();

    // Destructor
    ~UGrid ();

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor(_grid_actor);  if (_wire_actor!=NULL) win.AddActor(_wire_actor); }

    // Add methods
    int InsertNextPoint (double X, double Y, double Z) { return _pnts->InsertNextPoint (X,Y,Z);      }
    int InsertNextCell  (int Type, vtkIdList * PtIds)  { return _grid->InsertNextCell  (Type,PtIds); }

    // Set methods
    UGrid & SetColor     (char const * Name="peacock", double Opacity=0.8);
    void    SetPoint     (int Id, double X, double Y, double Z) { _pnts->SetPoint (Id, X,Y,Z); }
    void    Modified     () { _grid->Modified(); }
    void    SetWire      (bool WithWireframe=true);
    void    SetMaterial  (double Ambient, double Diffuse, double Specular, double SpecularPower);

private:
    vtkPoints           * _pnts;
    vtkUnstructuredGrid * _grid;
    vtkDataSetMapper    * _grid_mapper;
    vtkActor            * _grid_actor;
    vtkDataSetMapper    * _wire_mapper;
    vtkActor            * _wire_actor;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline UGrid::UGrid ()
    : _wire_actor (NULL)
{
    _pnts        = vtkPoints           ::New();
    _grid        = vtkUnstructuredGrid ::New();
    _grid_mapper = vtkDataSetMapper    ::New();
    _grid_actor  = vtkActor            ::New();
    _grid        -> SetPoints (_pnts);
    _grid_mapper -> SetInput  (_grid);
    _grid_actor  -> SetMapper (_grid_mapper);
    SetColor ();
}

inline UGrid::~UGrid ()
{
    _pnts        -> Delete();
    _grid        -> Delete();
    _grid_mapper -> Delete();
    _grid_actor  -> Delete();
    if (_wire_actor!=NULL)
    {
        _wire_mapper -> Delete();
        _wire_actor  -> Delete();
    }
}

inline UGrid & UGrid::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _grid_actor->GetProperty()->SetColor   (c(0), c(1), c(2));
    _grid_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline void UGrid::SetWire (bool WithWireframe)
{
    if (_wire_actor!=NULL)
    {
        _wire_mapper -> Delete();
        _wire_actor  -> Delete();
    }
    if (WithWireframe)
    {
        _wire_mapper = vtkDataSetMapper     ::New();
        _wire_actor  = vtkActor             ::New();
        _wire_mapper -> SetInput            (_grid);
        _wire_mapper -> ScalarVisibilityOff ();
        _wire_actor  -> SetMapper           (_wire_mapper);
        _wire_actor  -> GetProperty         ()->SetRepresentationToWireframe();
        Vec3_t c(Colors::Get("black"));
        _wire_actor->GetProperty()->SetColor (c.data());
    }
}

inline void UGrid::SetMaterial (double Ambient, double Diffuse, double Specular, double SpecularPower)
{
    _grid_actor->GetProperty()->SetAmbient       (Ambient);
    _grid_actor->GetProperty()->SetDiffuse       (Diffuse);
    _grid_actor->GetProperty()->SetSpecular      (Specular);
    _grid_actor->GetProperty()->SetSpecularPower (SpecularPower);
}

#endif // MECHSYS_POLYGON_H
