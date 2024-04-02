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

#ifndef MECHSYS_TEXT_H
#define MECHSYS_TEXT_H

// VTK
#include <vtkTextActor3D.h>
#include <vtkTextProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>

namespace VTK
{

class Text
{
public:
    // Constructor & Destructor
     Text () { _create(); }
    ~Text ();

    // Alternative constructor
    Text (Vec3_t const & X, char const * Txt) { _create();  SetPos(X);  SetText(Txt); }

    // Set methods
    Text & SetText  (char const * Txt)                                     { _text_actor->SetInput      (Txt);             return (*this); }
    Text & SetPos   (Vec3_t const & X)                                     { _text_actor->SetPosition   (X(0),X(1),X(2));  return (*this); }
    Text & SetOri   (double AlpX=90.0, double AlpY=90.0, double AlpZ=45.0) { _text_actor->SetOrientation(AlpX,AlpY,AlpZ);  return (*this); }
    Text & SetProps (double Scale=0.003, int SizePt=14, bool Shadow=true);
    Text & SetColor (char const * Name="blue");

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor(reinterpret_cast<vtkActor*>(_text_actor)); }

private:
    vtkTextActor3D  * _text_actor;
    vtkTextProperty * _text_prop;
    void _create ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Text::~Text ()
{
    _text_actor -> Delete();
    _text_prop  -> Delete();
}

inline Text & Text::SetProps (double Scale, int SizePt, bool Shadow)
{
    _text_actor -> SetScale    (Scale);
    _text_prop  -> SetFontSize (SizePt);
    _text_prop  -> SetShadow   (Shadow);
    return (*this);
}

inline Text & Text::SetColor (char const * Name)
{
    Vec3_t c(Colors::Get(Name));
    _text_prop->SetColor (c.data());
    return (*this);
}

inline void Text::_create ()
{
    _text_prop  = vtkTextProperty ::New();
    _text_actor = vtkTextActor3D  ::New();
    _text_actor -> SetTextProperty (_text_prop);
    SetOri   ();
    SetProps ();
    SetColor ();
}

}; // namespace VTK

#endif // MECHSYS_TEXT_H
