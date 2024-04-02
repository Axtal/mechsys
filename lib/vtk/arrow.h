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

#ifndef MECHSYS_ARROW3D_H
#define MECHSYS_ARROW3D_H

// Std Lib
#include <cmath>

// VTK
#include <vtkConeSource.h>
#include <vtkCylinderSource.h>
#include <vtkAppendPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>

// MechSys
#include <mechsys/vtk/win.h>
#include <mechsys/util/colors.h>
#include <mechsys/linalg/matvec.h>

/*                +-------------------------------+ 
 *                |            length             |
 *                +-----------------------+-------+
 *                |        bod_len        |tip_len|
 *                |                       |       |
 *                                        `.      ----+  
 *                                        | ``.       |
 *             +  +-----------------------|    ``.    |
 *     bod_rad |  |           +           |   +   >   | tip_rad   
 *             +  +-----------|-----------|   |_-'    |
 *                |           |           | _-|       |
 *                |           |           ''  |     --+  
 *                |           |               |
 *                +-----------+---------------+-------> y axis
 *                |           |               |    
 *                y0      y_bod_cen      y_tip_cen
 */

namespace VTK
{

class Arrow
{
public:
    // Constructor & Destructor
     Arrow ();
    ~Arrow ();

    // Alternative constructor
    Arrow (Vec3_t const & X0, Vec3_t const & V, double ConPct=0.1, double ConRad=0.03, double CylRad=0.015, int Res=20);

    // Set methods
    Arrow & SetGeometry   (double ConPct=0.1, double ConRad=0.03, double CylRad=0.015);
    Arrow & SetResolution (int Resolution=20);
    Arrow & SetColor      (char const * Name="yellow", double Opacity=1.0);
    Arrow & SetVector     (Vec3_t const & X0, Vec3_t const & V);
    Arrow & SetPoints     (Vec3_t const & X0, Vec3_t const & X1);

    // Methods
    void AddTo (VTK::Win & win) { win.AddActor(_arrow_actor); }

private:
    double               _con_pct;
    double               _tot_len;
    vtkConeSource      * _cone;
    vtkCylinderSource  * _cylin;
    vtkTransformFilter * _transform;
    vtkAppendPolyData  * _arrow;
    vtkPolyDataMapper  * _arrow_mapper;
    vtkActor           * _arrow_actor;
    void _create        ();
    void _update_length ();
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Arrow::Arrow ()
    : _tot_len (1.0)
{
    _create       ();
    SetGeometry   ();
    SetResolution ();
    SetColor      ();
}

inline Arrow::Arrow (Vec3_t const & X0, Vec3_t const & V, double ConPct, double ConRad, double CylRad, int Res)
    : _tot_len (norm(V))
{
    _create       ();
    SetGeometry   (ConPct, ConRad, CylRad);
    SetResolution (Res);
    SetColor      ();
    SetVector     (X0, V);
}

inline Arrow::~Arrow ()
{
    _cone         -> Delete();
    _cylin        -> Delete();
    _transform    -> Delete();
    _arrow        -> Delete();
    _arrow_mapper -> Delete();
    _arrow_actor  -> Delete();
}

inline Arrow & Arrow::SetGeometry (double ConPct, double ConRad, double CylRad)
{
    _con_pct = ConPct;
    _update_length ();
    _cone  -> SetRadius (ConRad);
    _cylin -> SetRadius (CylRad);
    return (*this);
}

inline Arrow & Arrow::SetResolution (int Res)
{
    _cone  -> SetResolution (Res);
    _cylin -> SetResolution (Res);
    return (*this);
}

inline Arrow & Arrow::SetColor (char const * Name, double Opacity)
{
    Vec3_t c(Colors::Get(Name));
    _arrow_actor->GetProperty()->SetColor   (c.data());
    _arrow_actor->GetProperty()->SetOpacity (Opacity);
    return (*this);
}

inline Arrow & Arrow::SetVector (Vec3_t const & X0, Vec3_t const & V)
{
    // update length
    _tot_len = Norm(V);
    _update_length ();

    // translate
    Vec3_t cen(X0+0.5*V); // center of arrow
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

inline Arrow & Arrow::SetPoints(Vec3_t const & X0, Vec3_t const & X1)
{
    SetVector (X0, X1-X0);
    return (*this);
}

inline void Arrow::_create ()
{
    _cone         = vtkConeSource      ::New();
    _cylin        = vtkCylinderSource  ::New();
    _arrow        = vtkAppendPolyData  ::New();
    _transform    = vtkTransformFilter ::New();
    _arrow_mapper = vtkPolyDataMapper  ::New();
    _arrow_actor  = vtkActor           ::New();
    _cone         -> SetDirection      (0.0, 1.0, 0.0);
    _arrow        -> AddInput          (_cone->GetOutput());
    _arrow        -> AddInput          (_cylin->GetOutput());
    _transform    -> SetInput          (_arrow->GetOutput());
    _arrow_mapper -> SetInput          (_transform->GetPolyDataOutput());
    _arrow_actor  -> SetMapper         (_arrow_mapper);
}

inline void Arrow::_update_length ()
{
    double con_len = _con_pct*_tot_len;      // cone length/height
    double cyl_len = _tot_len-con_len;       // cylinder length/height
    double con_cen = (_tot_len-con_len)/2.0; // cone center
    double cyl_cen = -con_len/2.0;           // cylinder center

    _cone  -> SetCenter (0.0, con_cen, 0.0);
    _cone  -> SetHeight (con_len);
    _cylin -> SetCenter (0.0, cyl_cen, 0.0);
    _cylin -> SetHeight (cyl_len);
}

}; // namespace VTK

#endif // MECHSYS_ARROW3D_H
