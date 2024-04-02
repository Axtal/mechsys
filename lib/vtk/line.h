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

#ifndef MECHSYS_LINE_H
#define MECHSYS_LINE_H

// VTK
#include <vtkLineSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>
#include <mechsys/util/string.h>
#include <mechsys/linalg/matvec.h>

namespace VTK
{

class Line
{
public:
    // Constructor & Destructor
     Line () { _create(); }
    ~Line ();

    // Alternative constructor
    Line (Vec3_t const & P, Vec3_t const & Q) { _create(P(0),P(1),P(2), Q(0),Q(1),Q(2)); }

    // Set methods
    //Line & SetLine      (Vec3_t const & P, Vec3_t const & Q);
    Line & SetWireColor (char const * Name="black") { Vec3_t c(Colors::Get(Name)); _line_actor->GetProperty()->SetColor(c.data()); return (*this); }
    Line & SetWireWidth (int Width=1.0)             { _line_actor->GetProperty()->SetLineWidth(Width); return (*this); }

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor(_line_actor); }

private:
    vtkLineSource     * _line_source;
    vtkPolyDataMapper * _line_mapper;
    vtkActor          * _line_actor;
    void _create (double P0=0, double P1=0, double P2=0, double Q0=1, double Q1=1, double Q2=1);
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Line::~Line ()
{
    _line_source -> Delete();
    _line_mapper -> Delete();
    _line_actor  -> Delete();
}

inline void Line::_create (double P0, double P1, double P2, double Q0, double Q1, double Q2)
{
    _line_source = vtkLineSource     ::New();
    _line_mapper = vtkPolyDataMapper ::New();
    _line_actor  = vtkActor          ::New();
    _line_source -> SetPoint1          (P0, P1, P2);
    _line_source -> SetPoint2          (Q0, Q1, Q2);
    _line_source -> Update             ();
    _line_mapper -> SetInputConnection (_line_source->GetOutputPort());
    _line_actor  -> SetMapper          (_line_mapper);
    SetWireColor();
}

}; // namespace VTK

#endif // MECHSYS_LINE_H
