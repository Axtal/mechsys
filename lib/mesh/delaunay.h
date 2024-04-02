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

#ifndef MECHSYS_MESH_DELAUNAY_H
#define MECHSYS_MESH_DELAUNAY_H

// STL
#include <list>
#include <map>
#include <vector>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Unique_hash_map.h>

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/mesh.h>

namespace Mesh
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel CGAL_K;
typedef CGAL::Triangulation_3<CGAL_K>                       CGAL_Triangulation;
typedef CGAL_Triangulation::Cell_handle                     CGAL_Cell_handle;
typedef CGAL_Triangulation::Vertex_handle                   CGAL_Vertex_handle;
typedef CGAL_Triangulation::Locate_type                     CGAL_Locate_type;
typedef CGAL_Triangulation::Point                           CGAL_Point;

class Delaunay : public virtual Mesh::Generic
{
public:
    // Constructor
    Delaunay (int NDim) : Mesh::Generic(NDim) {}

    // Methods
    void Reset     ();                                                                          ///< Reset cloud of points
    void AddPoint  (double X, double Y, double Z=0);                                            ///< Add new point to the list of points
    void AddPoints (Array<double> const & X, Array<double> const & Y, Array<double> const & Z); ///< Set all points
    void Generate  ();                                                                          ///< Generate

private:
    // Data
    std::list<CGAL_Point> _pts; ///< List of input points (3D)
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void Delaunay::Reset ()
{
    _pts.clear();
}

inline void Delaunay::AddPoint (double X, double Y, double Z)
{
    _pts.push_back(CGAL_Point(X,Y,Z));
}

inline void Delaunay::AddPoints (Array<double> const & X, Array<double> const & Y, Array<double> const & Z)
{
    for (size_t i=0; i<X.Size(); ++i)
    {
        _pts.push_back(CGAL_Point(X[i], Y[i], Z[i]));
    }
}

inline void Delaunay::Generate ()
{
    // erase previous mesh
    Erase();

    // generate Delaunay triangulation using CGAL
    CGAL_Triangulation tr(_pts.begin(), _pts.end());
    //printf("number of vertices                = %zd\n", tr.number_of_vertices());
    //printf("number of cells (including ghost) = %zd\n", tr.number_of_cells());

    // outputs dimension and number of vertices
    size_t n = tr.number_of_vertices();
    if (n == 0) return;

    // set vertices handle
    std::vector<CGAL_Triangulation::Vertex_handle> vh(n+1);
    size_t i = 0;
    for (CGAL_Triangulation::Vertex_iterator it=tr.vertices_begin(), end=tr.vertices_end(); it!=end; ++it)
    {
        vh[i++] = it;
        printf("%g %g %g\n", vh[i-1]->point().x(), vh[i-1]->point().y(), vh[i-1]->point().z());
    }
    CGAL_triangulation_assertion(i==n+1);
    CGAL_triangulation_assertion(tr.is_infinite(vh[0]));

    // vertices
    std::map<CGAL_Triangulation::Vertex_handle, size_t> V;
    V[tr.infinite_vertex()] = 0;
    for (i=1; i<=n; i++)
    {
        PushVert(-1, vh[i]->point().x(), vh[i]->point().y(), vh[i]->point().z());
        V[vh[i]] = i;
        printf("vert %3zd: %g %g %g\n", Verts.Size()-1, vh[i]->point().x(), vh[i]->point().y(), vh[i]->point().z());
    }

    // cells handle
    std::map<CGAL_Triangulation::Cell_handle, size_t > C;
    i = 0;

    // write the cells
    Array<int> con(4);
    size_t m = tr.tds().number_of_cells();
    for (CGAL_Triangulation::Cell_iterator it=tr.tds().cells_begin(); it!=tr.tds().cells_end(); ++it)
    {
        C[it] = i++;
        bool has_inf_vert = false;
        for (int j=0; j<4; j++)
        {
            con[j] = V.find(it->vertex(j))->second - 1;
            // has infinite vertex?
            if (con[j] < 0) {
                has_inf_vert = true;
            }
        }
        if (has_inf_vert) {
            continue; // skip cell with infinite vertex
        }
        PushCell(-1, con);
        //printf("cell %3zd: %d %d %d %d\n", Cells.Size()-1, con[0], con[1], con[2], con[3]);
    }
    CGAL_triangulation_assertion(i==m);
}

}; // namespace Mesh


#endif // MECHSYS_MESH_DELAUNAY_H
