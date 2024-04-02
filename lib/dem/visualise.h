/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
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

// VTK
#include <vtkPolygon.h>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/vtk/sphere.h>
#include <mechsys/vtk/cylinder.h>
#include <mechsys/vtk/ugrid.h>
#include <mechsys/vtk/win.h>

namespace DEM
{

class Visualise
{
public:
    // Constuctor & Destructor
     Visualise (DEM::Domain const & Dom, Array<int> const & TagParts, Array<int> const & TagWalls, bool ShowVert=false, bool ShowEdge=false, int ThRes=12, int PhiRes=12, int CylRes=12);
    ~Visualise ();

    // Methods
    void Update ();
    void AddTo  (VTK::Win & win);

    // Data
    DEM::Domain const & Dom;         // domain
    Array<int>          TagParts;    // tag of particles
    Array<int>          TagWalls;    // tag of walls
    bool                ShowVert;    // show all vertices of particles
    bool                ShowEdge;    // show all edges of particles
    String              PartColor;   // color of particles
    String              WallColor;   // color of walls
    double              PartOpacity; // opacity of particles
    double              WallOpacity; // opacity of walls

    // Spherical particles
    Array<VTK::Sphere*>   SphParts;  // spherical particles

    // Complex particles
    Array<VTK::Sphere*>   PartVerts; // particles: spheres
    Array<VTK::Cylinder*> PartEdges; // particles: cylinders
    UGrid                 PartFaces; // particles: all faces of particles (inner/outer)

    // Walls
    Array<VTK::Sphere*>   WallVerts; // walls: spheres
    Array<VTK::Cylinder*> WallEdges; // walls: cylinders
    UGrid                 WallFaces; // walls: faces
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Visualise::Visualise (DEM::Domain const & D, Array<int> const & TPs, Array<int> const & TWs, bool SV, bool SE, int ThRes, int PhiRes, int CylRes)
    : Dom         (D),
      ShowVert    (SV),
      ShowEdge    (SE),
      PartColor   ("brown"),
      WallColor   ("peacock"),
      PartOpacity (1.0),
      WallOpacity (0.1)
{
    TagParts = TPs;
    TagWalls = TWs;

    PartFaces.SetColor (PartColor.CStr(), PartOpacity);
    WallFaces.SetColor (WallColor.CStr(), WallOpacity);

    for (size_t i=0; i<Dom.Particles.Size(); ++i)
    {
        Particle & p = (*Dom.Particles[i]);

        // spherical particles
        if (p.Edges.Size()==0 && p.Faces.Size()==0)
        {
            SphParts.Push (new VTK::Sphere(p.x, p.Props.R, ThRes, PhiRes));
            SphParts[SphParts.Size()-1]->SetColor (PartColor.CStr(), PartOpacity);
        }

        else
        {
            // complex particles
            if (TagParts.Has(p.Tag))
            {
                // spheres
                if (ShowVert)
                {
                    for (size_t j=0; j<p.Verts.Size(); ++j)
                    {
                        PartVerts.Push (new VTK::Sphere((*p.Verts[j]), p.Props.R, ThRes, PhiRes));
                        PartVerts.Last () -> SetColor (PartColor.CStr(), PartOpacity);
                    }
                }

                // only edges
                if (p.Faces.Size()==0 || ShowEdge)
                {
                    for (size_t j=0; j<p.Edges.Size(); ++j)
                    {
                        PartEdges.Push (new VTK::Cylinder((*p.Edges[j]->X0), (*p.Edges[j]->X1), p.Props.R, /*cap*/true, CylRes));
                        PartEdges.Last () -> SetColor (PartColor.CStr(), PartOpacity);
                    }
                }

                // faces
                for (size_t j=0; j<p.Faces.Size(); ++j)
                {
                    Face & f = (*p.Faces[j]);
                    vtkPolygon * pl = vtkPolygon::New();
                    pl->GetPointIds()->SetNumberOfIds (f.Edges.Size());
                    for (size_t k=0; k<f.Edges.Size(); ++k)
                    {
                        Edge & e = (*f.Edges[k]);
                        int idx = PartFaces.InsertNextPoint ((*e.X0)(0), (*e.X0)(1), (*e.X0)(2));
                        pl->GetPointIds()->SetId (k, idx);
                    }
                    PartFaces.InsertNextCell (pl->GetCellType(), pl->GetPointIds());
                    pl->Delete();
                }
            }

            // walls
            else if (TagWalls.Has(p.Tag))
            {
                /*
                for (size_t j=0; j<p.Verts.Size(); ++j)
                {
                    WallVerts.Push (new VTK::Sphere((*p.Verts[j]), p.Props.R, ThRes, PhiRes));
                    WallVerts.Last () -> SetColor (WallColor.CStr(), WallOpacity);
                }

                for (size_t j=0; j<p.Edges.Size(); ++j)
                {
                    WallEdges.Push (new VTK::Cylinder((*p.Edges[j]->X0), (*p.Edges[j]->X1), p.Props.R, false, CylRes));
                    WallEdges.Last () -> SetColor (WallColor.CStr(), WallOpacity);
                }
                */

                for (size_t j=0; j<p.Faces.Size(); ++j)
                {
                    Face   & f = (*p.Faces[j]);
                    Vec3_t   n = cross(f.Edges[0]->dL, f.Edges[1]->dL);
                    n = n/norm(n);
                    vtkPolygon * pi = vtkPolygon::New();
                    vtkPolygon * po = vtkPolygon::New();
                    pi->GetPointIds()->SetNumberOfIds (f.Edges.Size());
                    po->GetPointIds()->SetNumberOfIds (f.Edges.Size());
                    for (size_t k=0; k<f.Edges.Size(); ++k)
                    {
                        Edge & e = (*f.Edges[k]);
                        Vec3_t vi((*e.X0) - p.Props.R*n);
                        Vec3_t vo((*e.X0) + p.Props.R*n);
                        int idx_i = WallFaces.InsertNextPoint (vi(0), vi(1), vi(2));
                        int idx_o = WallFaces.InsertNextPoint (vo(0), vo(1), vo(2));
                        pi->GetPointIds()->SetId (k, idx_i);
                        po->GetPointIds()->SetId (k, idx_o);
                    }
                    WallFaces.InsertNextCell (pi->GetCellType(), pi->GetPointIds());
                    WallFaces.InsertNextCell (po->GetCellType(), po->GetPointIds());
                    pi->Delete();
                    po->Delete();
                }
            }

            // unknown particle
            else throw new Fatal("Visualise::Visualise: Unknow particle: Tag = %d",p.Tag);
        }
    }
}

inline void Visualise::Update ()
{
    //Array<VTK::Sphere*>   * vts; // current vertices
    //Array<VTK::Cylinder*> * eds; // current edges
    //UGrid                 * fcs; // current faces

    // spherical particles
    size_t i_sph_part        = 0; // index of spherical particle

    // complex particles
    size_t i_part_vert       = 0; // index of particle vertex
    size_t i_part_edge       = 0; // index of particle edge
    size_t i_part_face_point = 0; // index of particle point at face

    // walls
    //size_t i_wall_vert       = 0; // index of wall vertex
    //size_t i_wall_edge       = 0; // index of wall edge
    size_t i_wall_face_point = 0; // index of wall point at face

    for (size_t i=0; i<Dom.Particles.Size(); ++i)
    {
        Particle & p = (*Dom.Particles[i]);

        // spherical particles
        if (p.Edges.Size()==0 && p.Faces.Size()==0)
        {
            SphParts[i_sph_part]->SetCenter (p.x);
            i_sph_part++;
        }

        else
        {
            // complex particles
            if (TagParts.Has(p.Tag))
            {
                // spheres
                if (ShowVert)
                {
                    for (size_t j=0; j<p.Verts.Size(); ++j)
                    {
                        PartVerts[i_part_vert]->SetCenter ((*p.Verts[j]));
                        i_part_vert++;
                    }
                }

                // only edges
                if (p.Faces.Size()==0 || ShowEdge)
                {
                    for (size_t j=0; j<p.Edges.Size(); ++j)
                    {
                        PartEdges[i_part_edge]->SetPoints ((*p.Edges[j]->X0), (*p.Edges[j]->X1));
                        i_part_edge++;
                    }
                }

                // faces
                for (size_t j=0; j<p.Faces.Size(); ++j)
                {
                    Face & f = (*p.Faces[j]);
                    for (size_t k=0; k<f.Edges.Size(); ++k)
                    {
                        Edge & e = (*f.Edges[k]);
                        PartFaces.SetPoint (i_part_face_point, (*e.X0)(0), (*e.X0)(1), (*e.X0)(2));
                        i_part_face_point++;
                    }
                }
            }

            // walls
            else if (TagWalls.Has(p.Tag))
            {
                /*
                for (size_t j=0; j<p.Verts.Size(); ++j)
                {
                    WallVerts[i_wall_vert]->SetCenter ((*p.Verts[j]));
                    i_wall_vert++;
                }

                for (size_t j=0; j<p.Edges.Size(); ++j)
                {
                    WallEdges[i_wall_edge]->SetPoints ((*p.Edges[j]->X0), (*p.Edges[j]->X1));
                    i_wall_edge++;
                }
                */

                for (size_t j=0; j<p.Faces.Size(); ++j)
                {
                    Face   & f = (*p.Faces[j]);
                    Vec3_t   n = cross(f.Edges[0]->dL, f.Edges[1]->dL);
                    n = n/norm(n);
                    for (size_t k=0; k<f.Edges.Size(); ++k)
                    {
                        Edge & e = (*f.Edges[k]);
                        Vec3_t vi((*e.X0) - p.Props.R*n);
                        Vec3_t vo((*e.X0) + p.Props.R*n);
                        WallFaces.SetPoint (i_wall_face_point++, vi(0), vi(1), vi(2));
                        WallFaces.SetPoint (i_wall_face_point++, vo(0), vo(1), vo(2));
                    }
                }
            }
        }
    }

    PartFaces.Modified();
    WallFaces.Modified();
}

inline Visualise::~Visualise ()
{
    for (size_t i=0; i<SphParts .Size(); ++i) delete SphParts [i];
    for (size_t i=0; i<PartVerts.Size(); ++i) delete PartVerts[i];
    for (size_t i=0; i<PartEdges.Size(); ++i) delete PartEdges[i];
    for (size_t i=0; i<WallVerts.Size(); ++i) delete WallVerts[i];
    for (size_t i=0; i<WallEdges.Size(); ++i) delete WallEdges[i];
}

inline void Visualise::AddTo (VTK::Win & win)
{

    PartFaces.AddTo (win);
    WallFaces.AddTo (win);
    for (size_t i=0; i<SphParts .Size(); ++i) SphParts [i]->AddTo (win);
    for (size_t i=0; i<PartVerts.Size(); ++i) PartVerts[i]->AddTo (win);
    for (size_t i=0; i<PartEdges.Size(); ++i) PartEdges[i]->AddTo (win);
    for (size_t i=0; i<WallVerts.Size(); ++i) WallVerts[i]->AddTo (win);
    for (size_t i=0; i<WallEdges.Size(); ++i) WallEdges[i]->AddTo (win);
}

}
