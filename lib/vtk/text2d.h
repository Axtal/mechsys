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

#ifndef MECHSYS_TEXT2D_H
#define MECHSYS_TEXT2D_H

// VTK
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>

namespace VTK
{

class Text2D
{
public:
    // Constructor & Destructor
     Text2D () { _create(); }
    ~Text2D ();

    // Alternative constructor
    Text2D (int x, int y, char const * Txt) { _create();  SetPos(x,y);  SetText(Txt); }

    // Set methods
    Text2D & SetText  (char const * Txt)                 { _text_actor->SetInput    (Txt);  return (*this); }
    Text2D & SetPos   (int x, int y)                     { _text_actor->SetPosition2(x,y);  return (*this); }
    Text2D & SetProps (int SizePt=18, bool Shadow=false);
    Text2D & SetColor (char const * Name="black");

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor(reinterpret_cast<vtkActor*>(_text_actor)); }

private:
    vtkTextActor    * _text_actor;
    vtkTextProperty * _text_prop;
    void _create ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Text2D::~Text2D ()
{
    _text_actor -> Delete();
    _text_prop  -> Delete();
}

inline Text2D & Text2D::SetProps (int SizePt, bool Shadow)
{
    _text_prop -> SetFontSize (SizePt);
    _text_prop -> SetShadow   (Shadow);
    return (*this);
}

inline Text2D & Text2D::SetColor (char const * Name)
{
    Vec3_t c(Colors::Get(Name));
    _text_prop->SetColor (c.data());
    return (*this);
}

inline void Text2D::_create ()
{
    _text_prop  = vtkTextProperty ::New();
    _text_actor = vtkTextActor    ::New();
    _text_actor -> SetTextProperty (_text_prop);
    SetProps ();
    SetColor ();
}

}; // namespace VTK

#endif // MECHSYS_TEXT2D_H
