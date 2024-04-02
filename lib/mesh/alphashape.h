/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_MESH_ALPHASHAPE_H
#define MECHSYS_MESH_ALPHASHAPE_H

// STL
#include <list>
#include <map>

// CGAL
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/mesh.h>

namespace Mesh
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel                      CGAL_K;
typedef CGAL::Alpha_shape_vertex_base_2<CGAL_K>                                  CGAL_VB_2;
typedef CGAL::Triangulation_face_base_2<CGAL_K>                                  CGAL_TF_2;
typedef CGAL::Alpha_shape_face_base_2<CGAL_K,CGAL_TF_2>                          CGAL_FB_2;
typedef CGAL::Triangulation_default_data_structure_2<CGAL_K,CGAL_VB_2,CGAL_FB_2> CGAL_TDS_2;
typedef CGAL::Delaunay_triangulation_2<CGAL_K,CGAL_TDS_2>                        CGAL_DT_2;
typedef CGAL::Alpha_shape_2<CGAL_DT_2>                                           CGAL_ALPHA_SHAPE_2;
typedef CGAL::Alpha_shape_vertex_base_3<CGAL_K>                                  CGAL_VB_3;
typedef CGAL::Alpha_shape_cell_base_3<CGAL_K>                                    CGAL_FB_3;
typedef CGAL::Triangulation_data_structure_3<CGAL_VB_3,CGAL_FB_3>                CGAL_TDS_3;
typedef CGAL::Delaunay_triangulation_3<CGAL_K,CGAL_TDS_3>                        CGAL_DT_3;
typedef CGAL::Alpha_shape_3<CGAL_DT_3>                                           CGAL_ALPHA_SHAPE_3;

class AlphaShape : public virtual Mesh::Generic
{
public:
    // Constructor
    AlphaShape (int NDim) : Mesh::Generic(NDim) {}

    // Methods
    void ResetCloud    ();                                    ///< Reset cloud of points
    void AddCloudPoint (double X, double Y, double Z=0);      ///< Add new point to the list of points i of input PSLG. SetCloudSize MUST be called first.
    void Generate      (double Alpha=-1, bool Regular=false); ///< Generate

private:
    // Data
    std::list<CGAL_K::Point_2> _pts_2d; ///< List of input points (2D)
    std::list<CGAL_K::Point_3> _pts_3d; ///< List of input points (3D)
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void AlphaShape::ResetCloud ()
{
    _pts_2d.clear();
    _pts_3d.clear();
}

inline void AlphaShape::AddCloudPoint (double X, double Y, double Z)
{
    if (NDim==3) _pts_3d.push_back (CGAL_K::Point_3(X,Y,Z));
    else         _pts_2d.push_back (CGAL_K::Point_2(X,Y));
}

inline void AlphaShape::Generate (double Alpha, bool Regular)
{
    if (NDim==3)
    {
        // Alpha-shape structure
        CGAL_ALPHA_SHAPE_3 as(_pts_3d.begin(), _pts_3d.end(), /*alpha*/0, (Regular ? CGAL_ALPHA_SHAPE_3::REGULARIZED : CGAL_ALPHA_SHAPE_3::GENERAL));

        // Find optimal alpha
        double alp = Alpha;
        if (alp<0) alp = (*as.find_optimal_alpha(1));

        // Generate alpha shape
        as.set_alpha (alp);
        if (as.number_of_solid_components()!=1) throw new Fatal("AlphaShape::Generate: There is a problem with AlphaShape 3D");
    }
    else
    {
        // Alpha-shape structure
        CGAL_ALPHA_SHAPE_2 as(_pts_2d.begin(), _pts_2d.end(), /*alpha*/0, (Regular ? CGAL_ALPHA_SHAPE_2::REGULARIZED : CGAL_ALPHA_SHAPE_2::GENERAL));

        // Find optimal alpha
        double alp = Alpha;
        if (alp<0) alp = (*as.find_optimal_alpha(1));

        // Generate alpha shape
        as.set_alpha (alp);

        // Erase old mesh
        Erase ();

        // Set Vertices
        size_t id = 0;
        std::map<CGAL_ALPHA_SHAPE_2::Vertex_handle,size_t> vs;
        for (CGAL_ALPHA_SHAPE_2::Alpha_shape_vertices_iterator it=as.alpha_shape_vertices_begin(); it!=as.alpha_shape_vertices_end(); ++it)
        {
            PushVert (-100, (*it)->point().x(), (*it)->point().y());
            vs[(*it)] = id++;
        }

        // Set Elements
        for (CGAL_ALPHA_SHAPE_2::Alpha_shape_edges_iterator it=as.alpha_shape_edges_begin(); it!=as.alpha_shape_edges_end(); ++it)
        {
            PushCell (-1, Array<int>((int)vs[it->first->vertex(as.ccw(it->second))], (int)vs[it->first->vertex(as.cw(it->second))]));
        }
    }
}

}; // namespace Mesh


#endif // MECHSYS_MESH_ALPHASHAPE_H
