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

#ifndef MECHSYS_TRIANGULATE_H
#define MECHSYS_TRIANGULATE_H

// VTK
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkUnstructuredGridWriter.h>

// MechSys
#include <mechsys/vtk/sgrid.h>
#include <mechsys/util/array.h>

namespace VTK
{

void Triangulate (double const * X, double const * Y, double const * Z, size_t Size, char const * Filename, double Tol=1.0e-3)
{
    // check
    if (Size<1) throw new Fatal("VTK::Triangulate: Size of arrays must be greater than zero");

    // create VTK points
    vtkPoints * points = vtkPoints::New();
    points->Allocate (Size);
    for (size_t i=0; i<Size; ++i) points->InsertPoint (i, X[i], Y[i], Z[i]);

    // Create a 3D triangulation
    //   - The Tolerance is the distance that nearly coincident points are merged together.
    //   - Delaunay does better if points are well spaced.
    //   - The alpha value is the radius of circumcircles, circumspheres.
    //     Any mesh entity whose circumcircle is smaller than this value is output.
    vtkPolyData   * vertices = vtkPolyData   ::New();
    vtkDelaunay3D * delaunay = vtkDelaunay3D ::New();
    vertices -> SetPoints    (points);
    delaunay -> SetInput     (vertices);
    delaunay -> SetTolerance (Tol);
    delaunay -> SetAlpha     (0.2);
    delaunay -> BoundingTriangulationOff();

    // Write file
    vtkUnstructuredGridWriter * writer = vtkUnstructuredGridWriter ::New();
    writer -> SetInputConnection (delaunay->GetOutputPort());
    writer -> SetFileName        (Filename);
    writer -> Write              ();

    // Clean up
    points   -> Delete();
    vertices -> Delete();
    delaunay -> Delete();
    writer   -> Delete();
}

void Triangulate (Array<double> const & X, Array<double> const & Y, Array<double> const & Z, char const * Filename, double Tol=1.0e-3)
{
    if (X.Size()<1)         throw new Fatal("VTK::Triangulate: Size of arrays must be greater than zero");
    if (X.Size()!=Y.Size()) throw new Fatal("VTK::Triangulate: X, Y, and Z arrays must have the same size");
    if (X.Size()!=Z.Size()) throw new Fatal("VTK::Triangulate: X, Y, and Z arrays must have the same size");
    Triangulate (X.GetPtr(), Y.GetPtr(), Z.GetPtr(), X.Size(), Filename, Tol);
}

void Triangulate (VTK::SGrid & G, char const * Filename, double Tol=1.0e-3)
{
    // Create a 3D triangulation
    //   - The Tolerance is the distance that nearly coincident points are merged together.
    //   - Delaunay does better if points are well spaced.
    //   - The alpha value is the radius of circumcircles, circumspheres.
    //     Any mesh entity whose circumcircle is smaller than this value is output.
    vtkPolyData   * vertices = vtkPolyData   ::New();
    vtkDelaunay3D * delaunay = vtkDelaunay3D ::New();
    vertices -> SetPoints    (G.GetPoints());
    delaunay -> SetInput     (vertices);
    delaunay -> SetTolerance (Tol);
    delaunay -> SetAlpha     (0.2);
    delaunay -> BoundingTriangulationOff();

    // Write file
    vtkUnstructuredGridWriter * writer = vtkUnstructuredGridWriter ::New();
    writer -> SetInputConnection (delaunay->GetOutputPort());
    writer -> SetFileName        (Filename);
    writer -> Write              ();

    // Clean up
    vertices -> Delete();
    delaunay -> Delete();
    writer   -> Delete();
}

}; // namespace VTK

#endif // MECHSYS_TRIANGULATE_H
