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

#ifndef MECHSYS_MESH_UNSTRUCTURED_H
#define MECHSYS_MESH_UNSTRUCTURED_H

/* LOCAL indexes of Vertices, Edges, and Faces

  2D:
             Nodes                 Edges

   y           2
   |           @                     @
   +--x       / \                   / \
           5 /   \ 4               /   \
            @     @             2 /     \ 1
           /       \             /       \
          /         \           /         \
         @-----@-----@         @-----------@
        0      3      1              0
*/

// STL
#include <iostream> // for cout, endl, ostream
#include <sstream>  // for ostringstream
#include <fstream>  // for ofstream
#include <cfloat>   // for DBL_EPSILON
#include <map>

// Jonathan R Shewchuk' Triangle
extern "C"
{
    #define REAL double
    #define ANSI_DECLARATORS
    #define VOID int
      #include "triangle.h"
    #undef REAL
    #undef ANSI_DECLARATORS
    #undef VOID
}

// Hang Si' Tetgen
#define TETLIBRARY
#include "tetgen.h"
#undef REAL
#undef TETLIBRARY

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/draw.h>

namespace Mesh
{


/////////////////////////////////////////////////////////////////////////////////////////// TriIO /////


/** JRS' Triangle Input/Output structure. */
typedef triangulateio TriIO;

/** HSI' Tetgen Input/Output structure. */
typedef tetgenio TetIO;

inline void TriAllocate (int NPoints, int NSegments, int NRegions, int NHoles, TriIO & Tio)
{
    // check
    if (NPoints<3)   throw new Fatal("Mesh::TriAllocate: At least 3 points are required. (%d is invalid)",NPoints);
    if (NSegments<3) throw new Fatal("Mesh::TriAllocate: At least 3 segments are required. (%d is invalid)",NSegments);

    // points
    Tio.pointlist       = (double*)malloc(NPoints*2*sizeof(double));
    Tio.numberofpoints  = NPoints;
    Tio.pointmarkerlist = (int*)malloc(NPoints*sizeof(int));

    // segments
    Tio.segmentlist       = (int*)malloc(NSegments*2*sizeof(int));
    Tio.segmentmarkerlist = (int*)malloc(NSegments * sizeof(int));
    Tio.numberofsegments  = NSegments;
    for (int i=0; i<NSegments; ++i) Tio.segmentmarkerlist[i]=0;

    // regions
    if (NRegions>0)
    {
        Tio.regionlist      = (double*)malloc(NRegions*4*sizeof(double));
        Tio.numberofregions = NRegions;
    }

    // holes
    if (NHoles>0)
    {
        Tio.holelist      = (double*)malloc(NHoles*2*sizeof(double));
        Tio.numberofholes = NHoles;
    }
}

inline void TriSetAllToNull (TriIO & Tio)
{
    // points
    Tio.pointlist               = NULL;
    Tio.pointattributelist      = NULL;
    Tio.pointmarkerlist         = NULL;
    Tio.numberofpoints          = 0;
    Tio.numberofpointattributes = 0;

    // triangles
    Tio.trianglelist               = NULL;
    Tio.triangleattributelist      = NULL;
    Tio.trianglearealist           = NULL;
    Tio.neighborlist               = NULL;
    Tio.numberoftriangles          = 0;
    Tio.numberofcorners            = 0;
    Tio.numberoftriangleattributes = 0;
    Tio.triedgemarks               = NULL;

    // segments
    Tio.segmentlist       = NULL;
    Tio.segmentmarkerlist = NULL;
    Tio.numberofsegments  = 0;

    // holes
    Tio.holelist      = NULL;
    Tio.numberofholes = 0;

    // regions
    Tio.regionlist      = NULL;
    Tio.numberofregions = 0;

    // edges
    Tio.edgelist       = NULL;
    Tio.edgemarkerlist = NULL;
    Tio.normlist       = NULL;
    Tio.numberofedges  = 0;
}

inline void TriDeallocateAll (TriIO & Tio)
{
    // Points
    if (Tio.pointlist          != NULL) free(Tio.pointlist);
    if (Tio.pointattributelist != NULL) free(Tio.pointattributelist);
    if (Tio.pointmarkerlist    != NULL) free(Tio.pointmarkerlist);

    // Triangles
    if (Tio.trianglelist          != NULL) free(Tio.trianglelist);
    if (Tio.triangleattributelist != NULL) free(Tio.triangleattributelist);
    if (Tio.trianglearealist      != NULL) free(Tio.trianglearealist);
    if (Tio.neighborlist          != NULL) free(Tio.neighborlist);
    if (Tio.triedgemarks          != NULL) free(Tio.triedgemarks);

    // Segments
    if (Tio.segmentlist       != NULL) free(Tio.segmentlist);
    if (Tio.segmentmarkerlist != NULL) free(Tio.segmentmarkerlist);

    // Holes
    if (Tio.holelist != NULL) free(Tio.holelist);

    // Regions
    if (Tio.regionlist != NULL) free(Tio.regionlist);

    // Edges
    if (Tio.edgelist       != NULL) free(Tio.edgelist);
    if (Tio.edgemarkerlist != NULL) free(Tio.edgemarkerlist);
    if (Tio.normlist       != NULL) free(Tio.normlist);

    // Clear all
    TriSetAllToNull (Tio);
}


/////////////////////////////////////////////////////////////////////////////////////////// Unstructured /////


class Unstructured : public virtual Mesh::Generic
{
public:
    // Constants
    static size_t FEM2TriPoint[]; ///< Map MechSys/FEM nodes to JRS-Triangle points
    static size_t FEM2TriEdge []; ///< Map MechSys/FEM nodes to JRS-Triangle edges
    static size_t FEM2TetPoint[]; ///< Map MechSys/FEM nodes to HSI-Tetgen points
    static size_t FEM2TetFace []; ///< Map MechSys/FEM nodes to HSI-Tetgen edges

    // Constructor
    Unstructured (int NDim);

    // Destructor
    ~Unstructured () { TriDeallocateAll(Tin); }

    /** 2D: Set Planar Straight Line Graph (PSLG)
     *  3D: Set Piecewise Linear Complex (PLC)
     * see tst/mesh01 for example */
    void Set    (size_t NPoints, size_t NSegmentsOrFacets, size_t NRegions, size_t NHoles);
    void SetReg (size_t iReg, int RTag, double MaxAreaOrVolume, double X, double Y, double Z=0.0);
    void SetHol (size_t iHol, double X, double Y, double Z=0.0);
    void SetPnt (size_t iPnt, int PTag, double X, double Y, double Z=0.0);
    void SetSeg (size_t iSeg, int ETag, int L, int R);
    void SetFac (size_t iFac, int FTag, Array<int> const & VertsOnFace);
    void SetFac (size_t iFac, int FTag, Array<int> const & Polygon1, Array<int> const & Polygon2);

    // Methods
    void Generate (bool O2=false, double GlobalMaxArea=-1, bool Quiet=true, double MinAngle=-1);       ///< Generate
    void WritePLY (char const * FileKey, bool Blender=true);                                           ///< (.ply)
    void GenBox   (bool O2=false, double MaxVolume=-1.0, double Lx=1.0, double Ly=1.0, double Lz=1.0); ///< Generate a cube with dimensions Lx,Ly,Lz and with tags on faces
    bool IsSet    () const;                                                                            ///< Check if points/edges/faces were already set

    // Alternative methods
    void Delaunay (Array<double> const & X, Array<double> const & Y, int Tag=-1); ///< Find Delaunay triangulation of a set of points
    void Delaunay (Array<double> const & X, Array<double> const & Y, Array<double> const & Z, int Tag=-1); ///< Find Delaunay tetrahedralization of a set of points in 3D

    // Data
    TriIO Tin; ///< Triangle structure: input PSLG
    TetIO Pin; ///< Tetgen structure: input PLC

#ifdef USE_BOOST_PYTHON
    void PySet (BPy::dict const & Dat);
#endif

private:
    // Data
    bool _lst_reg_set; ///< Was the last region (NRegions-1) set ?
    bool _lst_hol_set; ///< Was the last hole (NHoles-1) set ?
    bool _lst_pnt_set; ///< Was the last point (NPoints-1) set ?
    bool _lst_seg_set; ///< Was the last segment (NSegmentsOrFacets-1) set ?
    bool _lst_fac_set; ///< Was the last face (NSegmentsOrFacets-1) set ?
};

size_t Unstructured::FEM2TriPoint[] = {0,1,2,5,3,4};
size_t Unstructured::FEM2TriEdge [] = {0,1,2};
size_t Unstructured::FEM2TetPoint[] = {0,1,2,3,4,5,6,7,8,9};
size_t Unstructured::FEM2TetFace [] = {3,1,0,2};


/////////////////////////////////////////////////////////////////////////////////////////// PLC: Implementation /////

inline Unstructured::Unstructured (int NDim)
    : Mesh::Generic (NDim),
      _lst_reg_set  (false),
      _lst_hol_set  (false),
      _lst_pnt_set  (false),
      _lst_seg_set  (false),
      _lst_fac_set  (false)
{
    TriSetAllToNull  (Tin);
    Pin.deinitialize (); 
}

inline void Unstructured::Set (size_t NPoints, size_t NSegmentsOrFacets, size_t NRegions, size_t NHoles)
{

    //std::cout << "NRegions = " << NRegions << std::endl;
    //std::cout << "NPoints = " << NPoints << std::endl;
    //std::cout << "NFacets = " << NSegmentsOrFacets << std::endl;
    //std::cout << "NHoles = " << NHoles << std::endl;

    // check
    if (NPoints<3)           throw new Fatal("Mesh::Unstructured::Set: The number of points must be greater than 2. (%d is invalid)",NPoints);
    if (NSegmentsOrFacets<3) throw new Fatal("Mesh::Unstructured::Set: The number of segments or faces must be greater than 2. (%d is invalid)",NSegmentsOrFacets);
    if (NRegions<1)          throw new Fatal("Mesh::Unstructured::Set: The number of regions must be greater than 1. (%d is invalid)",NRegions);

    // flags
    _lst_reg_set = false;
    _lst_hol_set = false;
    _lst_pnt_set = false;
    _lst_seg_set = false;
    _lst_fac_set = false;

    if (NDim==2)
    {
        // erase previous PSLG
        TriDeallocateAll (Tin);

        // allocate PSLG
        TriAllocate (NPoints, NSegmentsOrFacets, NRegions, NHoles, Tin);
    }
    else if (NDim==3)
    {
        // erase previous PLC
        Pin.deinitialize ();

        // allocate PLC
        Pin.initialize ();
        
        // points
        Pin.firstnumber     = 0;
        Pin.numberofpoints  = NPoints;
        Pin.pointlist       = new double [NPoints*3];
        Pin.pointmarkerlist = new int [NPoints];

        // facets
        Pin.numberoffacets  = NSegmentsOrFacets;
        Pin.facetlist       = new TetIO::facet [NSegmentsOrFacets];
        Pin.facetmarkerlist = new int [NSegmentsOrFacets];

        // regions
        Pin.numberofregions = NRegions;
        Pin.regionlist      = new double [NRegions*5];

        // holes
        Pin.numberofholes = NHoles;
        Pin.holelist      = new double [NHoles*3];
    }
    else throw new Fatal("Unstructured::Set: NDim must be either 2 or 3. NDim==%d is invalid",NDim);
}

inline void Unstructured::SetReg (size_t iReg, int RTag, double MaxAreaOrVolume, double X, double Y, double Z)
{
    if (NDim==2)
    {
        Tin.regionlist[iReg*4  ] = X;
        Tin.regionlist[iReg*4+1] = Y;
        Tin.regionlist[iReg*4+2] = RTag;
        Tin.regionlist[iReg*4+3] = MaxAreaOrVolume;
        if ((int)iReg==Tin.numberofregions-1) _lst_reg_set = true;
    }
    else if (NDim==3)
    {
        Pin.regionlist[iReg*5  ] = X;
        Pin.regionlist[iReg*5+1] = Y;
        Pin.regionlist[iReg*5+2] = Z;
        Pin.regionlist[iReg*5+3] = RTag;
        Pin.regionlist[iReg*5+4] = MaxAreaOrVolume;
        if ((int)iReg==Pin.numberofregions-1) _lst_reg_set = true;
    }
}

inline void Unstructured::SetHol (size_t iHol, double X, double Y, double Z)
{
    if (NDim==2)
    {
        Tin.holelist[iHol*2  ] = X;
        Tin.holelist[iHol*2+1] = Y;
        if ((int)iHol==Tin.numberofholes-1) _lst_hol_set = true;
    }
    else if (NDim==3)
    {
        Pin.holelist[iHol*3  ] = X;
        Pin.holelist[iHol*3+1] = Y;
        Pin.holelist[iHol*3+2] = Z;
        if ((int)iHol==Pin.numberofholes-1) _lst_hol_set = true;
    }
}

inline void Unstructured::SetPnt (size_t iPnt, int PTag, double X, double Y, double Z)
{
    if (NDim==2)
    {
        Tin.pointlist[iPnt*2  ]   = X;
        Tin.pointlist[iPnt*2+1]   = Y;
        Tin.pointmarkerlist[iPnt] = PTag;
        if ((int)iPnt==Tin.numberofpoints-1) _lst_pnt_set = true;
    }
    else if (NDim==3)
    {
        Pin.pointlist[iPnt*3  ]   = X;
        Pin.pointlist[iPnt*3+1]   = Y;
        Pin.pointlist[iPnt*3+2]   = Z;
        Pin.pointmarkerlist[iPnt] = PTag;
        if ((int)iPnt==Pin.numberofpoints-1) _lst_pnt_set = true;
    }
}

inline void Unstructured::SetSeg (size_t iSeg, int ETag, int L, int R)
{
    if (NDim==3) throw new Fatal("Unstructured::SetSeg: This method must be called for 2D meshes only");
    Tin.segmentlist[iSeg*2  ]   = L;
    Tin.segmentlist[iSeg*2+1]   = R;
    Tin.segmentmarkerlist[iSeg] = ETag;
    if ((int)iSeg==Tin.numberofsegments-1) _lst_seg_set = true;
}

inline void Unstructured::SetFac (size_t iFac, int FTag, Array<int> const & VertsOnFace)
{
    if (NDim==2) throw new Fatal("Unstructured::SetSeg: This method must be called for 3D meshes only");

    Pin.facetmarkerlist[iFac] = FTag;
    TetIO::facet * f    = &Pin.facetlist[iFac];
    f->numberofpolygons = 1;
    f->polygonlist      = new TetIO::polygon [f->numberofpolygons];
    f->numberofholes    = 0;
    f->holelist         = NULL;

    // read vertices
    TetIO::polygon * p  = &f->polygonlist[0];
    int npoints         = static_cast<int>(VertsOnFace.Size());
    p->numberofvertices = npoints;
    p->vertexlist       = new int [npoints];
    for (int j=0; j<npoints; ++j) p->vertexlist[j] = VertsOnFace[j];

    if ((int)iFac==Pin.numberoffacets-1) _lst_fac_set = true;
}

inline void Unstructured::SetFac (size_t iFac, int FTag, Array<int> const & Polygon1, Array<int> const & Polygon2)
{
    if (NDim==2) throw new Fatal("Unstructured::SetSeg: This method must be called for 3D meshes only");

    Pin.facetmarkerlist[iFac] = FTag;
    TetIO::facet * f    = &Pin.facetlist[iFac];
    f->numberofpolygons = 2;
    f->polygonlist      = new TetIO::polygon [f->numberofpolygons];
    f->numberofholes    = 0;
    f->holelist         = NULL;

    // read vertices of polygon 1
    TetIO::polygon * p1  = &f->polygonlist[0];
    int np1              = static_cast<int>(Polygon1.Size());
    p1->numberofvertices = np1;
    p1->vertexlist       = new int [np1];
    for (int j=0; j<np1; ++j) p1->vertexlist[j] = Polygon1[j];

    // read vertices of polygon 2
    TetIO::polygon * p2  = &f->polygonlist[1];
    int np2              = static_cast<int>(Polygon2.Size());
    p2->numberofvertices = np2;
    p2->vertexlist       = new int [np2];
    for (int j=0; j<np2; ++j) p2->vertexlist[j] = Polygon2[j];

    if ((int)iFac==Pin.numberoffacets-1) _lst_fac_set = true;
}

inline void Unstructured::Generate (bool O2, double GlobalMaxArea, bool Quiet, double MinAngle)
{
    // check
    if (!IsSet()) throw new Fatal("Unstructured::Generate: Please, set the input data (regions,points,segments/facets) first.");

    // info
    Util::Stopwatch stopwatch(/*activated*/WithInfo);

    // parameters
    String prms("pzA"); // Q=quiet, p=poly, q=quality, z=zero
    if (Quiet)           prms.Printf("Q%s",   prms.CStr());
    if (GlobalMaxArea>0) prms.Printf("%sa%f", prms.CStr(), GlobalMaxArea);
    if (MinAngle>0)      prms.Printf("%sq%f", prms.CStr(), MinAngle);
    else                 prms.Printf("%sq",   prms.CStr());
    if (O2)              prms.Printf("%so2",  prms.CStr());
    prms.Printf("%sa", prms.CStr());

    if (NDim==2)
    {
        // generate
        TriIO tou;
        TriSetAllToNull (tou);
        triangulate (prms.CStr(), &Tin, &tou, NULL);

        // verts
        Verts.Resize (tou.numberofpoints);
        for (size_t i=0; i<Verts.Size(); ++i)
        {
            Verts[i]      = new Vertex;
            Verts[i]->ID  = i;
            Verts[i]->Tag = 0;
            Verts[i]->C   = tou.pointlist[i*2], tou.pointlist[i*2+1], 0.0;

            /* tou.pointmarkerlist[ipoint] will be equal to:
             * == edgeTag (<0) => on edge with tag <<<<<<<<<<<<<<<<<< REMOVED
             * == 0            => internal vertex (not on boundary)
             * == 1            => on boundary                   */
            int mark = tou.pointmarkerlist[i];
            if (mark<0)
            {
                Verts[i]->Tag = mark;
                TgdVerts.Push (Verts[i]);
            }
        }

        // cells
        Cells.Resize (tou.numberoftriangles);
        for (size_t i=0; i<Cells.Size(); ++i)
        {
            Cells[i]      = new Cell;
            Cells[i]->ID  = i;
            Cells[i]->Tag = tou.triangleattributelist[i*tou.numberoftriangleattributes];
            Cells[i]->PartID = 0; // partition id (domain decomposition)
            Cells[i]->V.Resize (tou.numberofcorners);
            for (size_t j=0; j<Cells[i]->V.Size(); ++j)
            {
                Share sha = {Cells[i],j};
                Cells[i]->V[j] = Verts[tou.trianglelist[i*tou.numberofcorners+FEM2TriPoint[j]]];
                Cells[i]->V[j]->Shares.Push (sha);
            }
            bool has_bry_tag = false;
            for (size_t j=0; j<3; ++j)
            {
                int edge_tag = tou.triedgemarks[i*3+FEM2TriEdge[j]];
                if (edge_tag<0)
                {
                    Cells[i]->BryTags[j] = edge_tag;
                    has_bry_tag          = true;
                }
            }
            if (has_bry_tag) TgdCells.Push (Cells[i]);
        }

        // clean up
        /* After triangulate (with -p switch), tou.regionlist gets the content of Tin.regionlist and
         * tou.holelist gets the content of Tin.holelist. Thus, these output variables must be set
         * to NULL in order to tell TriDeallocateAll to ignore them and do not double-free memory. */
        tou.regionlist      = NULL;
        tou.numberofregions = 0;
        tou.holelist        = NULL;
        tou.numberofholes   = 0;
        TriDeallocateAll (tou);
    }
    else
    {
        //for (size_t i=0; i<Pin.numberofpoints; ++i) std::cout << Util::_20_15 << Pin.pointlist[i*3] << " " << Util::_20_15 << Pin.pointlist[i*3+1] << " " << Util::_20_15 << Pin.pointlist[i*3+2] << std::endl;

        //prms = "CpJVVV"; // debug

        // generate
        prms.append("f");
        char sw[prms.size()+1];
        strcpy (sw, prms.CStr());
        TetIO pou;
        tetrahedralize (sw, &Pin, &pou);

        // verts
        Verts.Resize (pou.numberofpoints);
        for (size_t i=0; i<Verts.Size(); ++i)
        {
            Verts[i]      = new Vertex;
            Verts[i]->ID  = i;
            Verts[i]->Tag = 0;
            Verts[i]->C   = pou.pointlist[i*3], pou.pointlist[i*3+1], pou.pointlist[i*3+2];

            /* pou.pointmarkerlist[ipoint] will be equal to:
             * == faceTag (<0) => on face with tag <<<<<<<<<<<<<<<<<< REMOVED
             * == 0            => internal vertex (not on boundary)
             * == 1            => on boundary                   */
            int mark = pou.pointmarkerlist[i];
            if (mark<0)
            {
                Verts[i]->Tag = mark;
                TgdVerts.Push (Verts[i]);
            }
        }

        // cells
        Cells.Resize (pou.numberoftetrahedra);
        for (size_t i=0; i<Cells.Size(); ++i)
        {
            Cells[i]      = new Cell;
            Cells[i]->ID  = i;
            Cells[i]->Tag = pou.tetrahedronattributelist[i*pou.numberoftetrahedronattributes];
            Cells[i]->V.Resize (pou.numberofcorners);
            for (size_t j=0; j<Cells[i]->V.Size(); ++j)
            {
                Share sha = {Cells[i],j};
                Cells[i]->V[j] = Verts[pou.tetrahedronlist[i*pou.numberofcorners+FEM2TetPoint[j]]];
                Cells[i]->V[j]->Shares.Push (sha);
            }
        }

        // face tags
        for (std::map<int,tetgenio::facemarkers>::const_iterator p=pou.tetfacemarkers.begin(); p!=pou.tetfacemarkers.end(); ++p)
        {
            int  icell       = p->first;
            bool has_bry_tag = false;
            for (size_t j=0; j<4; ++j)
            {
                int face_tag = p->second.m[FEM2TetFace[j]];
                //std::cout << icell << " " << j << " " << face_tag << "\n";
                if (face_tag<0)
                {
                    Cells[icell]->BryTags[j] = face_tag;
                    has_bry_tag              = true;
                }
            }
            if (has_bry_tag) TgdCells.Push (Cells[icell]);
        }
    }

    // check
    if (Verts.Size()<1) throw new Fatal("Unstructured::Generate: Failed with %d vertices and %d cells", Verts.Size(), Cells.Size());
    if (Cells.Size()<1) throw new Fatal("Unstructured::Generate: Failed with %d vertices and %d cells", Verts.Size(), Cells.Size());

    // info
    if (WithInfo)
    {
        printf("\n%s--- Unstructured Mesh Generation --- %dD --- O%d -------------------------------------%s\n",TERM_CLR1,NDim,(O2?2:1),TERM_RST);
        printf("%s  %s    = %s%s\n",TERM_CLR4,(NDim==2?"Triangle command":"Tetgen command "),prms.CStr(),TERM_RST);
        printf("%s  Num of cells       = %zd%s\n", TERM_CLR2, Cells.Size(), TERM_RST);
        printf("%s  Num of vertices    = %zd%s\n", TERM_CLR2, Verts.Size(), TERM_RST);
    }
}

inline void Unstructured::WritePLY (char const * FileKey, bool Blender)
{
    // output string
    String fn(FileKey); fn.append(".ply");
    std::ostringstream oss;

    if (Blender)
    {
        // header
        oss << "import Blender\n";
        oss << "import bpy\n";

        // scene, mesh, and object
        oss << "scn = bpy.data.scenes.active\n";
        oss << "msh = bpy.data.meshes.new('unstruct_poly')\n";
        oss << "obj = scn.objects.new(msh,'unstruct_poly')\n";

        if (NDim==2)
        {
            // points
            oss << "pts = [";
            for (int i=0; i<Tin.numberofpoints; ++i)
            {
                oss << "[" << Tin.pointlist[i*2];
                oss << "," << Tin.pointlist[i*2+1] << ", 0.0]";
                if (i==Tin.numberofpoints-1) oss << "]\n";
                else                         oss << ",\n       ";
            }
            oss << "\n";

            // edges
            oss << "edg = [";
            for (int i=0; i<Tin.numberofsegments; ++i)
            {
                oss << "[" << Tin.segmentlist[i*2] << "," << Tin.segmentlist[i*2+1] << "]";
                if (i==Tin.numberofsegments-1) oss << "]\n";
                else                           oss << ",\n       ";
            }
            oss << "\n";
        }
        else
        {
            // points
            oss << "pts = [";
            for (int i=0; i<Pin.numberofpoints; ++i)
            {
                oss << "[" << Pin.pointlist[i*3];
                oss << "," << Pin.pointlist[i*3+1];
                oss << "," << Pin.pointlist[i*3+2] << "]";
                if (i==Pin.numberofpoints-1) oss << "]\n";
                else                         oss << ",\n       ";
            }
            oss << "\n";

            // edges
            oss << "edg = [";
            for (int i=0; i<Pin.numberoffacets; ++i)
            {
                TetIO::facet * f = &Pin.facetlist[i];
                for (int j=0; j<f->numberofpolygons; ++j)
                {
                    TetIO::polygon * p  = &f->polygonlist[j];
                    for (int k=1; k<p->numberofvertices; ++k)
                    {
                        oss << "[" << p->vertexlist[k-1] << "," << p->vertexlist[k] << "]";
                        if (k==p->numberofvertices-1) oss << ",[" << p->vertexlist[k] << "," << p->vertexlist[0] << "]";
                        if ((i==Pin.numberoffacets-1) && (j==f->numberofpolygons-1) && (k==p->numberofvertices-1)) oss << "]\n";
                        else oss << ",\n       ";
                    }
                }
            }
            oss << "\n";
        }

        // extend mesh
        oss << "msh.verts.extend(pts)\n";
        oss << "msh.edges.extend(edg)\n";
    }
    else // matplotlib
    {
        if (NDim==3) throw new Fatal("Unstructured::WritePLY: Method not available for 3D and MatPlotLib");

        // header
        MPL::Header (oss);

        // vertices and commands
        oss << "# vertices and commands\n";
        oss << "dat = []\n";
        for (int i=0; i<Tin.numberofsegments; ++i)
        {
            int I = Tin.segmentlist[i*2];
            int J = Tin.segmentlist[i*2+1];
            oss << "dat.append((PH.MOVETO, (" << Tin.pointlist[I*2] << "," << Tin.pointlist[I*2+1] << ")))\n";
            oss << "dat.append((PH.LINETO, (" << Tin.pointlist[J*2] << "," << Tin.pointlist[J*2+1] << ")))\n";
        }
        oss << "\n";

        // draw edges
        MPL::AddPatch (oss);

        // draw tags
        oss << "# draw tags\n";
        for (int i=0; i<Tin.numberofpoints; ++i)
        {
            int pt_tag = Tin.pointmarkerlist[i];
            if (pt_tag<0) oss << "ax.text(" << Tin.pointlist[i*2] << "," << Tin.pointlist[i*2+1] << ", " << pt_tag << ", ha='center', va='center', fontsize=14, backgroundcolor=lyellow)\n";
        }
        for (int i=0; i<Tin.numberofsegments; ++i)
        {
            int edge_tag = Tin.segmentmarkerlist[i];
            if (edge_tag<0)
            {
                int    I  = Tin.segmentlist[i*2];
                int    J  = Tin.segmentlist[i*2+1];
                double x0 = Tin.pointlist[I*2];
                double y0 = Tin.pointlist[I*2+1];
                double x1 = Tin.pointlist[J*2];
                double y1 = Tin.pointlist[J*2+1];
                double xm = (x0+x1)/2.0;
                double ym = (y0+y1)/2.0;
                oss << "ax.text(" << xm << "," << ym << ", " << edge_tag << ", ha='center', va='center', fontsize=14, backgroundcolor=pink)\n";
            }
        }
        oss << "\n";

        // show
        oss << "# show\n";
        oss << "axis ('scaled')\n";
        oss << "show ()\n";
    }

    // create file
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Unstructured::GenBox (bool O2, double MaxVolume, double Lx, double Ly, double Lz)
{
    Set    (8, 6, 1, 0);               // nverts, nfaces, nregs, nholes
    SetPnt (0, -1,   0.0,  0.0,  0.0); // id, vtag, x, y, z
    SetPnt (1, -2,    Lx,  0.0,  0.0);
    SetPnt (2, -3,    Lx,   Ly,  0.0);
    SetPnt (3, -4,   0.0,   Ly,  0.0);
    SetPnt (4, -5,   0.0,  0.0,   Lz);
    SetPnt (5, -6,    Lx,  0.0,   Lz);
    SetPnt (6, -7,    Lx,   Ly,   Lz);
    SetPnt (7, -8,   0.0,   Ly,   Lz);
    SetReg (0, -1, MaxVolume,  Lx/2., Ly/2., Lz/2.); // id, tag, max_vol, reg_x, reg_y, reg_z
    SetFac (0, -10, Array<int>(0,3,7,4));            // id, ftag, npolys, nverts, v0,v1,v2,v3
    SetFac (1, -20, Array<int>(1,2,6,5));
    SetFac (2, -30, Array<int>(0,1,5,4));
    SetFac (3, -40, Array<int>(2,3,7,6));
    SetFac (4, -50, Array<int>(0,1,2,3));
    SetFac (5, -60, Array<int>(4,5,6,7));
    Generate (O2);
}

inline bool Unstructured::IsSet () const
{
    if (NDim==2)
    {
        bool hol_ok = (Tin.numberofholes>0 ? _lst_hol_set : true);
        return (_lst_reg_set && hol_ok && _lst_pnt_set && _lst_seg_set);
    }
    if (NDim==3)
    {
        bool hol_ok = (Pin.numberofholes>0 ? _lst_hol_set : true);
        return (_lst_reg_set && hol_ok && _lst_pnt_set && _lst_fac_set);
    }
    return false;
}

inline void Unstructured::Delaunay (Array<double> const & X, Array<double> const & Y, int Tag)
{
    // check
    if (NDim==3)            throw new Fatal("Unstructured::Delaunay: This method is only available for 2D");
    if (X.Size()!=Y.Size()) throw new Fatal("Unstructured::Delaunay: Size of X and Y arrays must be equal (%d!=%d)",X.Size(),Y.Size());

    // erase previous PSLG
    TriDeallocateAll (Tin);

    // allocate only the pointlist array inside Triangle's IO structure
    size_t npoints = X.Size();
	Tin.pointlist      = (double*)malloc(npoints*2*sizeof(double));
	Tin.numberofpoints = npoints;

    // set points
    for (size_t i=0; i<npoints; ++i)
    {
        Tin.pointlist[0+i*2] = X[i];
        Tin.pointlist[1+i*2] = Y[i];
    }

    // triangulate
    TriIO tou;
    TriSetAllToNull (tou);
    triangulate ("Qz", &Tin, &tou, NULL); // Quiet, zero-based

    // verts
    Verts.Resize (tou.numberofpoints);
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        Verts[i]      = new Vertex;
        Verts[i]->ID  = i;
        Verts[i]->Tag = 0;
        Verts[i]->C   = tou.pointlist[i*2], tou.pointlist[i*2+1], 0.0;
    }

    // cells
    Cells.Resize (tou.numberoftriangles);
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        Cells[i]      = new Cell;
        Cells[i]->ID  = i;
        Cells[i]->Tag = Tag;
        Cells[i]->V.Resize (tou.numberofcorners);
        for (size_t j=0; j<Cells[i]->V.Size(); ++j)
        {
            Share sha = {Cells[i],j};
            Cells[i]->V[j] = Verts[tou.trianglelist[i*tou.numberofcorners+FEM2TriPoint[j]]];
            Cells[i]->V[j]->Shares.Push (sha);
        }
    }

    // clean up
    /* After triangulate (with -p switch), tou.regionlist gets the content of Tin.regionlist and
     * tou.holelist gets the content of Tin.holelist. Thus, these output variables must be set
     * to NULL in order to tell TriDeallocateAll to ignore them and do not double-free memory. */
    //tou.regionlist      = NULL;
    //tou.numberofregions = 0;
    //tou.holelist        = NULL;
    //tou.numberofholes   = 0;
    TriDeallocateAll (tou);
}

inline void Unstructured::Delaunay (Array<double> const & X, Array<double> const & Y, Array<double> const & Z, int Tag)
{
    // check
    if (NDim==2)            throw new Fatal("Unstructured::Delaunay: This method is only available for 3D");
    if (X.Size()!=Y.Size()) throw new Fatal("Unstructured::Delaunay: Size of X and Y arrays must be equal (%d!=%d)",X.Size(),Y.Size());
    if (Z.Size()!=Y.Size()) throw new Fatal("Unstructured::Delaunay: Size of Z and Y arrays must be equal (%d!=%d)",Z.Size(),Y.Size());

    // points
    Pin.deinitialize();
    Pin.initialize();
    size_t NPoints = X.Size();
    Pin.firstnumber     = 0;
    Pin.numberofpoints  = NPoints;
    Pin.pointlist       = new double [NPoints*3];
    Pin.pointmarkerlist = new int [NPoints];

    for (size_t i=0; i<X.Size(); ++i)
    {
        Pin.pointlist[3*i  ] = X[i];
        Pin.pointlist[3*i+1] = Y[i];
        Pin.pointlist[3*i+2] = Z[i];
    }

    TetIO pou;
    tetrahedralize ((char*)"Qz", &Pin, &pou);

    // verts
    Verts.Resize (pou.numberofpoints);
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        Verts[i]      = new Vertex;
        Verts[i]->ID  = i;
        Verts[i]->Tag = 0;
        Verts[i]->C   = pou.pointlist[i*3], pou.pointlist[i*3+1], pou.pointlist[i*3+2];

        /* pou.pointmarkerlist[ipoint] will be equal to:
         * == faceTag (<0) => on face with tag <<<<<<<<<<<<<<<<<< REMOVED
         * == 0            => internal vertex (not on boundary)
         * == 1            => on boundary                   */
        int mark = pou.pointmarkerlist[i];
        if (mark<0)
        {
            Verts[i]->Tag = mark;
            TgdVerts.Push (Verts[i]);
        }
    }

    // cells
    Cells.Resize (pou.numberoftetrahedra);
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        Cells[i]      = new Cell;
        Cells[i]->ID  = i;
        //Cells[i]->Tag = pou.tetrahedronattributelist[i*pou.numberoftetrahedronattributes];
        Cells[i]->Tag = Tag-i;
        Cells[i]->V.Resize (pou.numberofcorners);
        for (size_t j=0; j<Cells[i]->V.Size(); ++j)
        {
            Share sha = {Cells[i],j};
            Cells[i]->V[j] = Verts[pou.tetrahedronlist[i*pou.numberofcorners+FEM2TetPoint[j]]];
            Cells[i]->V[j]->Shares.Push (sha);
        }
    }

    // face tags
    for (std::map<int,tetgenio::facemarkers>::const_iterator p=pou.tetfacemarkers.begin(); p!=pou.tetfacemarkers.end(); ++p)
    {
        int  icell       = p->first;
        bool has_bry_tag = false;
        for (size_t j=0; j<4; ++j)
        {
            int face_tag = p->second.m[FEM2TetFace[j]];
            //std::cout << icell << " " << j << " " << face_tag << "\n";
            if (face_tag<0)
            {
                Cells[icell]->BryTags[j] = face_tag;
                has_bry_tag              = true;
            }
        }
        if (has_bry_tag) TgdCells.Push (Cells[icell]);
    }
}


#ifdef USE_BOOST_PYTHON

inline void Unstructured::PySet (BPy::dict const & Dat)
{
    BPy::list const & pts = BPy::extract<BPy::list>(Dat["pts"])(); // points
    BPy::list const & rgs = BPy::extract<BPy::list>(Dat["rgs"])(); // regions
    BPy::list const & hls = BPy::extract<BPy::list>(Dat["hls"])(); // holes
    BPy::list const & con = BPy::extract<BPy::list>(Dat["con"])(); /// segments/facets (connectivity)

    size_t NPoints           = BPy::len(pts);
    size_t NSegmentsOrFacets = BPy::len(con);
    size_t NRegions          = BPy::len(rgs);
    size_t NHoles            = BPy::len(hls);

    // allocate memory
    Set (NPoints, NSegmentsOrFacets, NRegions, NHoles);

    // read regions
    for (size_t i=0; i<NRegions; ++i) SetReg (i, BPy::extract<int   >(rgs[i][0])(),         // tag
                                                 BPy::extract<double>(rgs[i][1])(),         // MaxArea/MaxVolume
                                                 BPy::extract<double>(rgs[i][2])(),         // x
                                                 BPy::extract<double>(rgs[i][3])(),         // y
                                      (NDim==3 ? BPy::extract<double>(rgs[i][4])() : 0.0)); // z

    // read holes
    for (size_t i=0; i<NHoles; ++i) SetHol (i, BPy::extract<double>(hls[i][0])(),         // x
                                               BPy::extract<double>(hls[i][1])(),         // y
                                    (NDim==3 ? BPy::extract<double>(hls[i][2])() : 0.0)); // z

    // read points
    for (size_t i=0; i<NPoints; ++i) SetPnt (i, BPy::extract<int   >(pts[i][0])(),         // tag
                                                BPy::extract<double>(pts[i][1])(),         // x
                                                BPy::extract<double>(pts[i][2])(),         // y
                                     (NDim==3 ? BPy::extract<double>(pts[i][3])() : 0.0)); // y

    // read segments
    if (NDim==2)
    {
        for (size_t i=0; i<NSegmentsOrFacets; ++i) SetSeg (i, BPy::extract<int>(con[i][0])(),  // tag
                                                              BPy::extract<int>(con[i][1])(),  // L
                                                              BPy::extract<int>(con[i][2])()); // R
    }

    // read facets
    else
    {
        for (size_t i=0; i<NSegmentsOrFacets; ++i) SetFac (i,            BPy::extract<int      >(con[i][0])(),   // tag
                                                              Array<int>(BPy::extract<BPy::list>(con[i][1])())); // verts
    }
}

#endif

}; // namespace Mesh

#endif // MECHSYS_MESH_UNSTRUCTURED_H
