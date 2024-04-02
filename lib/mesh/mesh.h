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

#ifndef MECHSYS_MESH_H
#define MECHSYS_MESH_H

// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>   // for DBL_EPSILON
#include <cstdarg>  // for va_list, va_start, va_end
#include <cstdlib>  // for atoi
#include <map>

// Boost
#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>

// ParMETIS
#if defined(HAS_PARMETIS) && defined(HAS_MPI)
  #include <mpi.h>
extern "C" {
  #include <metis.h>
}
#endif

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/vtkcelltype.h>

namespace Mesh
{

#ifndef VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

#define NOEDGE    {-1,-1,-1}
#define NOFACE    {-1,-1,-1,-1}
#define NOEDGES   {NOEDGE,NOEDGE,NOEDGE,NOEDGE}
#define NOEDGES3D {NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE}
#define NOFACES   {NOFACE,NOFACE,NOFACE,NOFACE,NOFACE,NOFACE}

size_t MaxNVerts2D               = 8;
size_t NVertsToNEdges2D[]        = {0,0,0,3,4,0,3,0,4};
size_t NVertsToNVertsPerEdge2D[] = {0,0,0,2,2,0,3,0,3};
int    NVertsToEdge2D[][4/*edges at most*/][3/*verts per edge at most*/]=
{
    NOEDGES,NOEDGES,NOEDGES,               // 0,1,2 verts
    {{0,1,-1},{1,2,-1},{2,0,-1},NOEDGE},   // 3 verts => TRIANGLE
    {{0,1,-1},{1,2,-1},{2,3,-1},{3,0,-1}}, // 4 verts => QUAD
    NOEDGES,                               // 5 verts
    {{0,1,3},{1,2,4},{2,0,5},NOEDGE},      // 6 verts => O2 TRIANGLE
    NOEDGES,                               // 7 verts
    {{0,1,4},{1,2,5},{2,3,6},{3,0,7}}      // 8 verts => O2 QUAD
};

size_t MaxNVerts3D               = 20;
size_t NVertsToNFaces3D[]        = {0,0,0,0, 4, 0,0,0, 6, 0, 4, 0,0,0,0,0,0,0,0,0, 6};
size_t NVertsToNVertsPerFace3D[] = {0,0,0,0, 3, 0,0,0, 4, 0, 3, 0,0,0,0,0,0,0,0,0, 4}; // disregarding O2 nodes
int    NVertsToFace3D[][6/*faces at most*/][4/*verts per face at most*/]=
{
    NOFACES,NOFACES,NOFACES,NOFACES,                                         //  0,1,2,3 verts
    {{0,3,2,-1},{0,1,3,-1},{0,2,1,-1},{1,2,3,-1},NOFACE,NOFACE},             //  4 verts => TETRA
    NOFACES,NOFACES,NOFACES,                                                 //  5,6,7 verts
    {{0,4,7,3},{1,2,6,5},{0,1,5,4},{2,3,7,6},{0,3,2,1},{4,5,6,7}},           //  8 verts => HEX
    NOFACES,                                                                 //  9 verts
    {{0,3,2,-1},{0,1,3,-1},{0,2,1,-1},{1,2,3,-1},NOFACE,NOFACE},             // 10 verts => O2 TETRA
    NOFACES,NOFACES,NOFACES,NOFACES,NOFACES,NOFACES,NOFACES,NOFACES,NOFACES, // 11,12,13,14,15,16,17,18,19 verts
    {{0,4,7,3},{1,2,6,5},{0,1,5,4},{2,3,7,6},{0,3,2,1},{4,5,6,7}}            // 20 verts => O2 HEX
};

size_t NVertsToNEdges3D[] = {0,0,0,0, 6, 0,0,0, 12, 0, 6, 0,0,0,0,0,0,0,0,0, 12};
int    NVertsToEdge3D[][12/*edges at most*/][3/*verts per edge at most*/]=
{
    NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,                                                                       //  0,1,2,3 verts
    {{0,1,-1},{1,2,-1},{2,0,-1},{0,3,-1},{1,3,-1},{2,3,-1},NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE},             //  4 verts => TETRA
    NOEDGES3D,NOEDGES3D,NOEDGES3D,                                                                                 //  5,6,7 verts
    {{0,1,-1},{1,2,-1},{2,3,-1},{3,0,-1},{4,5,-1},{5,6,-1},{6,7,-1},{7,4,-1},{0,4,-1},{1,5,-1},{2,6,-1},{3,7,-1}}, //  8 verts => HEX
    NOEDGES3D,                                                                                                     //  9 verts
    {{0,1,4},{1,2,5},{2,0,6},{0,3,7},{1,3,8},{2,3,9},NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE},                   // 10 verts => O2 TETRA
    NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,                     // 11,12,13,14,15,16,17,18,19 verts
    {{0,1,8},{1,2,9},{2,3,10},{3,0,11},{4,5,12},{5,6,13},{6,7,14},{7,4,15},{0,4,16},{1,5,17},{2,6,18},{3,7,19}}    // 20 verts => O2 HEX
};

#define NONE 0
size_t NVertsToGeoNVerts2D[] = {NONE,              // 0 verts
                                NONE,              // 1 vert
                                 2,                // 2 verts
                                 3,                // 3
                                 4,                // 4
                                NONE,              // 5
                                 3,                // 6
                                NONE,              // 7
                                 4,                // 8
                                NONE,NONE,NONE,    // 9,10,11
                                NONE,NONE,NONE,    // 12,13,14
                                 3};               // 15

size_t NVertsToGeoNVerts3D[] = {NONE,                       //  0 verts
                                NONE,                       //  1 vert
                                 2,                         //  2 verts
                                NONE,                       //  3
                                 4,                         //  4
                                NONE,                       //  5
                                NONE,                       //  6
                                NONE,                       //  7
                                 8,                         //  8
                                NONE,                       //  9
                                 4,                         // 10
                                NONE,NONE,NONE,NONE,NONE,   // 11,12,13,14,15
                                NONE,NONE,NONE,NONE,        // 16,17,18,19
                                 8};                        // 20
#undef NONE

#define BRYKEY(num_verts,idx_cell,idx_bry)                                      \
    int vert_a, vert_b, vert_c=-1, vert_d=-1;                                   \
    if (NDim==2)                                                                \
    {                                                                           \
        vert_a = Cells[idx_cell]->V[NVertsToEdge2D[num_verts][idx_bry][0]]->ID; \
        vert_b = Cells[idx_cell]->V[NVertsToEdge2D[num_verts][idx_bry][1]]->ID; \
        Util::Sort (vert_a,vert_b);                                             \
    }                                                                           \
    else                                                                        \
    {                                                                           \
        vert_a = Cells[idx_cell]->V[NVertsToFace3D[num_verts][idx_bry][0]]->ID; \
        vert_b = Cells[idx_cell]->V[NVertsToFace3D[num_verts][idx_bry][1]]->ID; \
        vert_c = Cells[idx_cell]->V[NVertsToFace3D[num_verts][idx_bry][2]]->ID; \
        if (NVertsToNVertsPerFace3D[num_verts]>3)                               \
        vert_d = Cells[idx_cell]->V[NVertsToFace3D[num_verts][idx_bry][3]]->ID; \
        Util::Sort (vert_a,vert_b,vert_c,vert_d);                               \
    }                                                                           \
    BryKey_t brykey(vert_a,vert_b,vert_c,vert_d);

#undef NOEDGE
#undef NOFACE
#undef NOEDGES
#undef NOFACES

struct Cell; ///< Forward declaration due to the following definitions

struct Share
{
    Cell * C; ///< The cell
    size_t N; ///< Local node index. Example: 2D=>0,1,2,3, 3D=>0,1,2,3,4,5,6,7
};

struct Vertex
{
    size_t       ID;      ///< ID
    int          Tag;     ///< Tag
    Vec3_t       C;       ///< X, Y, and Z coordinates
    Array<Share> Shares;  ///< IDs of cells sharing this vertex
    Array<int>   PartIDs; ///< Partitions (domains) in which this vertex is located in
};

typedef std::map<int,int>                      BryTag_t;   ///< Map: edge/face ID => edge/face tag
typedef boost::tuple<int,int,int,int>          BryKey_t;   ///< Edge/face key = (left_node,right_node) for edges or (node0,node1,node2) for faces
typedef std::pair<int,Cell*>                   NeighDat_t; ///< Pair: local edge/face ID, neighbour cell
typedef std::map<BryKey_t, NeighDat_t>         Neighs_t;   ///< Map: edge/face key => neighbour data
typedef std::map<BryKey_t, Array<NeighDat_t> > BryCell_t;  ///< Map: edge/face key => cells sharing this boundary (edge/face)
typedef std::map<Vertex*,Array<Vertex*> >      Pin_t;      ///< Pin type

struct Cell
{
    size_t         ID;      ///< ID
    int            Tag;     ///< Cell tag. Required for setting up of attributes, for example.
    Array<Vertex*> V;       ///< Connectivity
    BryTag_t       BryTags; ///< Boundary (edge/face) tags: map iEdgeFace => Tag
    Neighs_t       Neighs;  ///< Neighbours information: map edge/face key => neighbour data (local edge/face ID, pointer to neighbour Cell)
    int            PartID;  ///< Partition (domain) ID
    Vec3_t         Xmin;    ///< Bounding box
    Vec3_t         Xmax;    ///< Bounding box
};

inline bool InsideCell(Vec3_t & V, Cell const * Ce)
{
    bool inside = false;
    double vtot = 0.0;
    Vec3_t v0 = Ce->V[0]->C;
    Vec3_t v1 = Ce->V[1]->C;
    Vec3_t v2 = Ce->V[2]->C;
    Vec3_t v3 = Ce->V[3]->C;
    Vec3_t a  = v1 - v0;
    Vec3_t b  = v2 - v0;
    Vec3_t c  = v3 - v0;
    double Vc = fabs(dot(a,cross(b,c)))/6.0;
    a = v1-V;
    b = v2-V;
    c = v3-V;
    double V0 = fabs(dot(a,cross(b,c)))/6.0;
    a = v0-V;
    b = v2-V;
    c = v3-V;
    double V1 = fabs(dot(a,cross(b,c)))/6.0;
    a = v1-V;
    b = v0-V;
    c = v3-V;
    double V2 = fabs(dot(a,cross(b,c)))/6.0;
    a = v1-V;
    b = v2-V;
    c = v0-V;
    double V3 = fabs(dot(a,cross(b,c)))/6.0;
    double dv = fabs(V0+V1+V2+V3-Vc)/Vc;
    if (dv<1.0e-8) inside = true;
    return inside;
}

inline bool InsideCell(Vec3_t & x, Cell const * Ce, Vec3_t & rst)
{
    Vec3_t x0 = Ce->V[0]->C;
    Vec3_t x1 = Ce->V[1]->C;
    Vec3_t x2 = Ce->V[2]->C;
    Vec3_t x3 = Ce->V[3]->C;

    Vec3_t B  = x0 - x;
    Mat3_t M;

    M(0,0) = x0(0) - x1(0); M(0,1) = x0(0) - x2(0); M(0,2) = x0(0) - x3(0);
    M(1,0) = x0(1) - x1(1); M(1,1) = x0(1) - x2(1); M(1,2) = x0(1) - x3(1);
    M(2,0) = x0(2) - x1(2); M(2,1) = x0(2) - x2(2); M(2,2) = x0(2) - x3(2);
    
    if (!SolAlt(M,B,rst)) return false;

    if (rst(0)>0.0&&rst(1)>0.0&&rst(2)>0.0&&rst(0)+rst(1)+rst(2)<1.0) return true;
}

inline int CellGetBryTag(int iSide, Cell const *c)
{
    BryTag_t::const_iterator p = c->BryTags.find(iSide);
    if (p==c->BryTags.end()) return 0;
    else                     return p->second;
}

inline int CellGetNeigh(int iSide, Cell const *c)
{
    for (Neighs_t::const_iterator p=c->Neighs.begin(); p!=c->Neighs.end(); ++p)
    {
        if (iSide==p->second.first) return p->second.second->ID;
    }
    return -1; // no neighbour
}

class Generic
{
public:
    // Constructor
    Generic (int TheNDim) : NDim(TheNDim), IsShell(false), WithInfo(true) {}

    // Destructor
    virtual ~Generic () { Erase(); }

    // Set methods
    void ReadMesh   (char const * FileKey, bool IsShell=false);           ///< (.mesh) Erase old mesh and read mesh from python file
    void SetSize    (size_t NumVerts, size_t NumCells);                   ///< Erase old mesh and set number of vertices
    void SetVert    (int iVert, int Tag, double X, double Y, double Z=0); ///< Set vertex
    void SetCell    (int iCell, int Tag, Array<int> const & Con);         ///< Set element ... => connectivity
    void SetBryTag  (int iCell, int iEdgeFace, int Tag);                  ///< Set element's edge or face tag
    void FindNeigh  ();                                                   ///< Find neighbours of each cell
    void GenO2Verts ();                                                   ///< Generate O2 (mid) vertices
    void Erase      ();                                                   ///< Erase current mesh (deallocate memory)

    // Tag methods (can be called after the mesh is created)
    void TagVertex    (int Tag, int Idx);                                                ///< Tag a vertex given its iVert
    int  TagVertex    (int Tag, double x, double y, double z=0.0, double Tol=1.0e-5);    ///< Tag a vertex. Returns the index of vertex in Verts array
    void TagVertsLine (int Tag, double x0, double y0, double AlpDeg, double Tol=1.0e-5); ///< Tag vertices along Line y = y0 + tan(AlpRad)*x
    void TagSegsLine  (int Tag, double x0, double y0, double AlpDeg, double Tol=1.0e-5); ///< Tag segments/edges along line y = y0 + tan(AlpRad)*x
    void TagHorzSegs  (int Tag, double y, double xMin, double xMax, double Tol=1.0e-5);  ///< Tag horizontal segment/edge inside xMin and xMax
    void GroundTags   (int LTag=-10, int RTag=-10, int BTag=-20, double Tol=1.0e-5);     ///< Set tags for the left and right (vertical) edges and bottom (horizontal) edge

    // Push methods
    int PushVert (int Tag, double X, double Y, double Z=0); ///< Push vertex
    int PushCell (int Tag, Array<int> const & Con);         ///< Push element

    // Method to extend mesh
    void AddLinCells (Array<int> const & IDsOrTags); ///< Set linear cells given edge tags (edg_tag,edg_tag,...) or pair of vertices (v0,v1,new_elem_tag, v0,v1,new_elem_tag, ...) or both
    void AddPin      (int VertexIdOrTag);

    // Methods
    void Check    (double Tol=1.0e-3) const;                         ///< Check if there are overlapping nodes
    void WriteVTU (char const * FileKey, int VolSurfOrBoth=0) const; ///< (.vtu) Write output file for ParaView. Vol=0, Surf=1, Both=2
    void WriteMPY (char const * FileKey, bool WithTags=true, bool WithIDs=true, bool WithShares=false,
                   char const * Extra=NULL) const;                   ///< (.mpy) Write Python script that calls mesh_drawing.py
    void WriteJSON(char const * FileKey) const;

    // Auxiliar methods
    int  FindVert    (int Tag)                                                 const; ///< Find ID of first vertex with Tag. Returns -1 if not found
    void BoundingBox (Vec3_t & Min, Vec3_t & Max)                              const; ///< Limits of mesh
    void ThrowError  (std::istringstream & iss, char const * Message)          const; ///< Used in ReadMesh
    void Adjacency   (Array<int> & Xadj, Array<int> & Adjncy, bool Full=false);       ///< Find list of adjacent elements
    void PartDomain  (int NParts, bool Full=false, int * Part=NULL);                  ///< Partition domain

    // Other methods
    void GenGroundSG (Array<double> const & X, Array<double> const & Y, double FootingLx=-1); ///< Generate ground square/box according to Smith and Griffiths' numbering
    void GenGroundSG (size_t Nx, size_t Ny, double Dx=1.0, double Dy=1.0);                    ///< Smith-Griffiths' ground
    void GenSector   (size_t Nr, size_t Nth, double r, double R, double ThetaRad);            ///< Generate a Circular sector
    void Quad4ToTri3 (bool BackDiagonal=true);              ///< Convert Quad4 mesh to Tri3 mesh
    void Quad8ToTri6 (bool BackDiagonal=true);              ///< Convert Quad8 mesh to Tri6 mesh
    void Tri6ToTri15 ();                                    ///< Convert Tri6 mesh to Tri15 mesh
    void Tri3ToTri6  ();                                    ///< Convert Tri3 mesh to Tri6 mesh
    void Refine      ();                                    ///< Refine frames: split Lin2 cells into two
    void Tri6Shape   (double r, double s, Vec_t & N) const; ///< Tri6 Shape functions (used in Tri6ToTri15). N(6) needs to be pre-resized

    // Data
    int            NDim;      ///< Space dimension
    bool           IsShell;   ///< Is shell mesh ? (only surface)
    Array<Vertex*> Verts;     ///< Vertices
    Array<Cell*>   Cells;     ///< Cells
    Array<Vertex*> TgdVerts;  ///< Tagged Vertices
    Array<Cell*>   TgdCells;  ///< Tagged Cells
    BryCell_t      Bry2Cells; ///< map: bry (edge/face ids) => neighbours cells
    Pin_t          Pins;      ///< Pins
    bool           WithInfo;  ///< Output information ?

#ifdef USE_BOOST_PYTHON
    void PySetCell       (int iCell, int Tag, BPy::list const Con) { SetCell (iCell, Tag, Array<int>(Con)); }
    void PyAddLinCells   (BPy::list const & IDsOrTags)             { AddLinCells (Array<int>(IDsOrTags)); }
    void PyGetVertsEdges (BPy::list & V, BPy::list & E) const;
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


std::ostream & operator<< (std::ostream & os, BryKey_t const & BK)
{
    os << "(" << boost::get<0>(BK) << "," << boost::get<1>(BK) << "," << boost::get<2>(BK) << "," << boost::get<3>(BK) << ")";
    return os;
}

std::ostream & operator<< (std::ostream & os, Generic const & M)
{
    // verts
    os << "V=[";
    for (size_t i=0; i<M.Verts.Size(); ++i)
    {
        if (M.Verts[i]==NULL) throw new Fatal("Mesh::operator<<: Vertex %d was not defined",i);
        os << "[" << Util::_4  << M.Verts[i]->ID;
        os << "," << Util::_4  << M.Verts[i]->Tag;
        os << "," << Util::_8s << M.Verts[i]->C(0);
        os << "," << Util::_8s << M.Verts[i]->C(1);  if (M.NDim==3)
        os << "," << Util::_8s << M.Verts[i]->C(2);
        if (i==M.Verts.Size()-1) os << "]]\n";
        else                     os << "],\n   ";
    }

    // cells
    os << "\nC=[";
    for (size_t i=0; i<M.Cells.Size(); ++i)
    {
        if (M.Cells[i]==NULL) throw new Fatal("Mesh::operator<<: Cell %d was not defined",i);
        os << "[" << Util::_4 << M.Cells[i]->ID;
        os << "," << Util::_4 << M.Cells[i]->Tag;
        os << ", [";
        for (size_t j=0; j<M.Cells[i]->V.Size(); ++j)
        {
            os << M.Cells[i]->V[j]->ID;
            if (j!=M.Cells[i]->V.Size()-1) os << ",";
        }
        os << "], {";
        size_t k = 0;
        for (BryTag_t::const_iterator p=M.Cells[i]->BryTags.begin(); p!=M.Cells[i]->BryTags.end(); ++p)
        {
            os << p->first << ":" << p->second;
            if (k!=M.Cells[i]->BryTags.size()-1) os << ",";
            k++;
        }
        os << "}";
        os << ", " << Util::_4 << M.Cells[i]->PartID;
        os << ", {";
        for (Neighs_t::const_iterator p=M.Cells[i]->Neighs.begin(); p!=M.Cells[i]->Neighs.end(); ++p)
        {
            if (p!=M.Cells[i]->Neighs.begin()) os << ",";
            os << "(" << boost::get<0>(p->first);
            os << "," << boost::get<1>(p->first);
            os << "," << boost::get<2>(p->first);
            os << "," << boost::get<3>(p->first);
            os << "):";
            os << "[" << p->second.first << "," << p->second.second->ID << "]";
        }
        os << "}";
        if (i==M.Cells.Size()-1) os << "]]\n";
        else                     os << "],\n   ";
    }
    return os;
}

inline int Generic::FindVert (int Tag) const
{
    for (size_t i=0; i<TgdVerts.Size(); ++i)
    {
        if (TgdVerts[i]->Tag==Tag) return TgdVerts[i]->ID;
    }
    return -1;
}

inline void Generic::BoundingBox (Vec3_t & Min, Vec3_t & Max) const
{
    Min = Verts[0]->C(0), Verts[0]->C(1), Verts[0]->C(2);
    Max = Verts[0]->C(0), Verts[0]->C(1), Verts[0]->C(2);
    for (size_t i=1; i<Verts.Size(); ++i)
    {
        if (Verts[i]->C(0)<Min(0)) Min(0) = Verts[i]->C(0);
        if (Verts[i]->C(1)<Min(1)) Min(1) = Verts[i]->C(1);
        if (Verts[i]->C(2)<Min(2)) Min(2) = Verts[i]->C(2);
        if (Verts[i]->C(0)>Max(0)) Max(0) = Verts[i]->C(0);
        if (Verts[i]->C(1)>Max(1)) Max(1) = Verts[i]->C(1);
        if (Verts[i]->C(2)>Max(2)) Max(2) = Verts[i]->C(2);
    }
}

inline void Generic::ThrowError (std::istringstream & iss, char const * Message) const
{
    size_t num = 1+iss.tellg();
    String str(iss.str());
    str.resize (num);
    size_t pos = str.find(" ");
    while (pos!=String::npos)
    {
        str.replace (pos,1,"");
        pos = str.find(" ",pos+1);
    }
    throw new Fatal("Generic::ReadMesh: Mesh file format invalid\n    %s\n    %s",Message,str.CStr());
}

inline void Generic::Adjacency (Array<int> & Xadj, Array<int> & Adjncy, bool Full)
{
    Xadj.Push (0);
    if (Full)
    {
        for (size_t i=0; i<Cells.Size(); ++i)
        {
            size_t nv = (NDim==2 ? NVertsToGeoNVerts2D[Cells[i]->V.Size()] : NVertsToGeoNVerts3D[Cells[i]->V.Size()]);
            Array<Cell*> neigh;
            for (size_t j=0; j<nv; ++j)
            {
                Array<Share> const & sha = Cells[i]->V[j]->Shares;
                for (size_t k=0; k<sha.Size(); ++k)
                {
                    if (sha[k].C!=Cells[i])
                    {
                        if (neigh.Find(sha[k].C)<0)
                        {
                            neigh .Push (sha[k].C);
                            Adjncy.Push (sha[k].C->ID);
                        }
                    }
                }
            }
            Xadj.Push (Xadj.Last()+neigh.Size());
        }
    }
    else
    {
        FindNeigh ();
        for (size_t i=0; i<Cells.Size(); ++i)
        {
            for (Neighs_t::const_iterator p=Cells[i]->Neighs.begin(); p!=Cells[i]->Neighs.end(); ++p)
            {
                Adjncy.Push (p->second.second->ID);
            }
            Xadj.Push (Xadj.Last()+Cells[i]->Neighs.size());
        }
    }
}

inline void Generic::PartDomain (int NParts, bool Full, int * Part)
{
    // allocate/set array with partition ids
    int  n        = Cells.Size();
    bool del_part = false;
    int  * part;
    if (Part==NULL)
    {
        part     = new int [n];
        del_part = true;
        if (NParts==1) for (int i=0; i<n; ++i) part[i] = 0;
    }
    else part = Part;

    // partition domain
    if (Part==NULL && NParts>1)
    {
#if defined(HAS_PARMETIS) && defined(HAS_MPI)
        Array<int> Xadj, Adjncy;
        Adjacency (Xadj, Adjncy, Full);
        int wgtflag    = 0; // No weights
        int numflag    = 0; // zero numbering
        int options[5] = {0,0,0,0,0};
        int edgecut;
        if (NParts<8) METIS_PartGraphRecursive (&n, Xadj.GetPtr(), Adjncy.GetPtr(), NULL, NULL, &wgtflag, &numflag, &NParts, options, &edgecut, part);
        else          METIS_PartGraphKway      (&n, Xadj.GetPtr(), Adjncy.GetPtr(), NULL, NULL, &wgtflag, &numflag, &NParts, options, &edgecut, part);
#else
        throw new Fatal("Generic::PartDomain: This method requires ParMETIS and MPI (if Part array is not provided)");
#endif
    }

    // find domains of elements
    for (int i=0; i<n; ++i) Cells[i]->PartID = part[i];
    if (del_part) delete [] part;

    // find domais of nodes
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        Array<Share> const & sha = Verts[i]->Shares;
        for (size_t k=0; k<sha.Size(); ++k)
        {
            if (Verts[i]->PartIDs.Find(sha[k].C->PartID)<0) Verts[i]->PartIDs.Push (sha[k].C->PartID);
        }
    }
}

inline void Generic::ReadMesh (char const * FileKey, bool Shell)
{
    // shell mesh ?
    IsShell = Shell;

    // open file
    String fn(FileKey); fn.append(".mesh");
    std::fstream fil(fn.CStr(), std::ios::in);
    if (!fil.is_open()) throw new Fatal("Generic::ReadMesh: Could not open file < %s >",fn.CStr());

    String line, str, comma, colon, buf;
    size_t vid,  cid,  bid;  // vertex id,  cell id,  boundary id
    int    vtag, ctag, btag; // vertex tag, cell tag, boundary tag
    double x,y,z=0.0;
    bool   reading_verts = false;
    bool   reading_cells = false;
    bool   verts_read    = false;
    while (!fil.eof())
    {
        // read line
        std::getline (fil,line);

        // add spaces around: = [ ] ,
        size_t pos;
        pos=line.find("="); while (pos!=String::npos) { line.replace(pos,1," = "); pos=line.find("=",pos+3); }
        pos=line.find("["); while (pos!=String::npos) { line.replace(pos,1," [ "); pos=line.find("[",pos+3); }
        pos=line.find("]"); while (pos!=String::npos) { line.replace(pos,1," ] "); pos=line.find("]",pos+3); }
        pos=line.find(","); while (pos!=String::npos) { line.replace(pos,1," , "); pos=line.find(",",pos+3); }
        pos=line.find("{"); while (pos!=String::npos) { line.replace(pos,1," { "); pos=line.find("{",pos+3); }
        pos=line.find("}"); while (pos!=String::npos) { line.replace(pos,1," } "); pos=line.find("}",pos+3); }
        pos=line.find(":"); while (pos!=String::npos) { line.replace(pos,1," : "); pos=line.find(":",pos+3); }

        // parse
        std::istringstream iss(line);
        while (iss>>str)
        {
            if (str=="V")
            {
                iss>>str>>str>>str; //  = [ [
                reading_cells = false;
                reading_verts = true;
            }
            else if (str=="C")
            {
                iss>>str>>str>>str; //  = [ [
                reading_verts = false;
                reading_cells = true;
            }
            if (reading_verts)
            {
                if (str!="[") ThrowError (iss, "Verts: opening bracket '[' lacking");
                while (iss>>vid)
                {
                    // read vertex data
                    iss >> comma >> vtag >> comma >> x >> comma >> y;
                    if (NDim==3) iss >> comma >> z;
                    iss >> str;
                    if (str!="]") ThrowError (iss, "Verts: closing brackets ']' lacking");
                    iss >> str;
                    if (str=="]") { reading_verts = false; verts_read = true; }
                    else if (str!=",") ThrowError (iss, "Verts: last vertex: ',' or ']' lacking");

                    // allocate vertex
                    int new_vid = PushVert (vtag, x, y, z);
                    if (new_vid!=(int)vid) throw new Fatal("Generic::ReadMesh: Verts IDs must be sequential. Vertex ID (%d) must be equal to (%d)",vid,new_vid);

                    break;
                }
            }
            if (reading_cells)
            {
                if (!verts_read) throw new Fatal("Generic::ReadMesh: List with vertices data must come before list with cells data");
                if (str!="[") ThrowError (iss, "Cells: opening brackets '[' lacking");
                while (iss>>cid)
                {
                    // read cell data
                    iss >> comma >> ctag >> comma >> str;
                    if (str!="[") ThrowError (iss, "Cells: opening brackets of connectivity '[' lacking");

                    // read connectivity
                    Array<int> con;
                    while (iss>>vid)
                    {
                        con.Push (Verts[vid]->ID);
                        iss >> str;
                        if (str=="]") break;
                        else if (str!=",") ThrowError (iss, "Cells: closing brackets of connectivity: ',' or ']' lacking");
                    }

                    // allocate cell
                    int new_cid = PushCell (ctag, con);
                    if (new_cid!=(int)cid) throw new Fatal("Generic::ReadMesh: Cells IDs must be sequential. Cell ID (%d) must be equal to (%d)",cid,new_cid);

                    iss >> comma;
                    iss >> str;
                    if (str!="{") ThrowError (iss, "Cells: opening braces of bry tags '{' lacking");

                    // read boundary tags
                    while (iss>>str)
                    {
                        if (str=="}") break;
                        else bid = atoi(str.CStr());
                        iss >> colon >> btag >> str;
                        Cells[cid]->BryTags[bid] = btag;
                        if (TgdCells.Find(Cells[cid])<0) TgdCells.Push (Cells[cid]);
                        if (str=="}") break;
                        else if (str!=",") ThrowError (iss, "Cells: closing braces of bry tags: ',' or '}' lacking");
                    }
                    iss >> str;
                    if (str!="]") ThrowError (iss, "Cells: closing brackets ']' lacking");
                    iss >> str;
                    if (str=="]") reading_cells = false;
                    else if (str!=",") ThrowError (iss, "Cells: last cell: ',' or ']' lacking");

                    break;
                }
            }
        }
    }
}

inline void Generic::SetSize (size_t NumVerts, size_t NumCells)
{
    // erase previous mesh
    Erase ();

    // nodes
    Verts.Resize    (NumVerts);
    Verts.SetValues (NULL);
    TgdVerts.Resize (0);

    // elements
    Cells.Resize    (NumCells);
    Cells.SetValues (NULL);
    TgdCells.Resize (0);
}

inline int Generic::PushVert (int Tag, double X, double Y, double Z)
{
    Verts.Push (NULL);
    int i = static_cast<int>(Verts.Size())-1;
    SetVert (i, Tag, X, Y, Z);
    return i;
}

inline void Generic::SetVert (int i, int Tag, double X, double Y, double Z)
{
    // set Verts
    if (Verts[i]==NULL) Verts[i] = new Vertex;
    Verts[i]->ID  = i;
    Verts[i]->Tag = Tag;
    Verts[i]->C   = X, Y, Z;

    // set TgdVerts
    if (Tag<0) TgdVerts.Push (Verts[i]);
}

inline int Generic::PushCell (int Tag, Array<int> const & Con)
{
    Cells.Push (NULL);
    int i = static_cast<int>(Cells.Size())-1;
    SetCell (i, Tag, Con);
    return i;
}

inline void Generic::SetCell (int i, int Tag, Array<int> const & Con)
{
    // number of vertices
    size_t NVerts = Con.Size();

    // check
    if (NDim==2) { if (NVertsToVTKCell2D[NVerts]<0) throw new Fatal("Generic::SetCell: In two-dimensions (NDim=2), Cell=%d with Tag=%d has invalid NVerts=%d",i,Tag,NVerts); }
    else         { if (NVertsToVTKCell3D[NVerts]<0) throw new Fatal("Generic::SetCell: In three-dimensions (NDim=3), Cell=%d with Tag=%d has invalid NVerts=%d",i,Tag,NVerts); }

    // set elems
    if (Cells[i]==NULL) Cells[i] = new Cell;
    Cells[i]->ID  = i;
    Cells[i]->Tag = Tag;
    Cells[i]->V.Resize (NVerts);

    // set connectivity
    double sum_z = 0.0;
    for (size_t j=0; j<NVerts; ++j)
    {
        int ivert = Con[j];
        if (Verts[ivert]==NULL) throw new Fatal("Generic::SetCell: Vert=%d of Cell %d, %d could not be found. Vertices must be set (with SetVert) before calling this method",ivert,i,Tag);
        Cells[i]->V[j] = Verts[ivert];
        Share sha = {Cells[i],j};
        Verts[ivert]->Shares.Push (sha);
        if (NDim==3) sum_z += Verts[ivert]->C(2);
    }

    // check connectivity (2D)
    if (NDim==2 && Cells[i]->V.Size()>2)
    {
        Vec3_t p0, p1, p2;
        for (size_t j=1; j<Cells[i]->V.Size()-1; ++j)
        {
            p0 = Cells[i]->V[j-1]->C - Cells[i]->V[j]->C;
            p1 = Cells[i]->V[j+1]->C - Cells[i]->V[j]->C;
            p2 = blitz::cross (p1,p0);
            if (p2(2)<0.0) throw new Fatal("Generic::SetCell: Order of vertices is incorrect (it must be counter-clockwise)");
            if (fabs(p2(0))>1.0e-10 || fabs(p2(1)>1.0e-10)) throw new Fatal("Generic::SetCell: In 2D, all vertices of cells must lie on the x-y plane");
        }
    }

    // check
    if (NDim==3 and fabs(sum_z)<1.0e-10) throw new Fatal("Generic::SetCell: In 3D, the sum of all z coordinates of a Cell (%d, %d) must not be zero",i,Tag);
}

inline void Generic::AddLinCells (Array<int> const & IDsOrTags)
{
    if (NDim==3) throw new Fatal("Generic::AddLinCells: Method not available for 3D yet");

    bool with_tags = (IDsOrTags[0]<0 ? true : false);
    size_t  nitems = (with_tags ? IDsOrTags.Size() : IDsOrTags.Size()/3);

    size_t idx = 0;
    for (size_t i=0; i<nitems; ++i)
    {
        if (with_tags) // tag given
        {
            int  etag  = IDsOrTags[idx];
            bool found = false;
            for (size_t j=0; j<TgdCells.Size(); ++j)
            {
                BryTag_t & eftags = TgdCells[j]->BryTags;
                for (BryTag_t::iterator p=eftags.begin(); p!=eftags.end(); ++p)
                {
                    if (etag==p->second)
                    {
                        int    eid    = p->first;
                        size_t nv     = TgdCells[j]->V.Size();
                        int    ivert0 = TgdCells[j]->V[NVertsToEdge2D[nv][eid][0]]->ID;
                        int    ivert1 = TgdCells[j]->V[NVertsToEdge2D[nv][eid][1]]->ID;
                        int    ivert2 = -1;
                        if (NVertsToNVertsPerEdge2D[nv]>2) ivert2 = TgdCells[j]->V[NVertsToEdge2D[nv][eid][2]]->ID;
                        if (ivert2<0)
                        {
                            Cells.Push (NULL);
                            SetCell (Cells.Size()-1, etag, Array<int>(ivert0,ivert1));
                        }
                        else
                        {
                            Cells.Push (NULL);  SetCell (Cells.Size()-1, etag, Array<int>(ivert0,ivert2));
                            Cells.Push (NULL);  SetCell (Cells.Size()-1, etag, Array<int>(ivert2,ivert1));
                        }
                        found  = true;
                        TgdCells[j]->BryTags.erase(p); // delete edge tag from element
                        break;
                    }
                }
            }
            if (!found) throw new Fatal("Generic::SetLinCells: Could not find cell with edge with tag = %d",etag);
            idx++;
        }
        else // vertices IDs given
        {
            int ivert0 = IDsOrTags[idx];  idx++;
            int ivert1 = IDsOrTags[idx];  idx++;
            int tag    = IDsOrTags[idx];  idx++;
            Cells.Push (NULL);
            SetCell (Cells.Size()-1, tag, Array<int>(ivert0,ivert1));
        }
    }
}

inline void Generic::AddPin (int VertIdOrTag)
{
    // find vertex
    Vertex * vert = NULL;
    if (VertIdOrTag<0) // vertex tag given
    {
        bool found = false;
        for (size_t i=0; i<TgdVerts.Size(); ++i)
        {
            if (VertIdOrTag==TgdVerts[i]->Tag)
            {
                vert  = Verts[TgdVerts[i]->ID];
                found = true;
                break;
            }
        }
        if (!found) throw new Fatal("Mesh::AddPin: Could not find vertex with tag = %d", VertIdOrTag);
    }
    else vert = Verts[VertIdOrTag];

    // break linear elements connected to this vert
    Array<size_t> lins_to_disconnect;
    size_t idx_lin = 0;
    for (size_t i=0; i<vert->Shares.Size(); ++i)
    {
        Cell * cell = vert->Shares[i].C;
        if (cell->V.Size()==2) // cell sharing this vertex is a lin2 (linear cell with 2 vertices)
        {
            if (idx_lin>0) // do this for the second lin2 sharing this vertex on
            {
                // add new vertex
                Verts.Push (NULL);
                //SetVert (Verts.Size()-1, vert->Tag, vert->C(0)-idx_lin*0.5, vert->C(1)-idx_lin*0.5, vert->C(2));
                SetVert (Verts.Size()-1, vert->Tag, vert->C(0), vert->C(1), vert->C(2));

                // set map of pins
                Pins[vert].Push (Verts[Verts.Size()-1]);

                // set connectivity of lin2
                int idx = vert->Shares[i].N; // local vertex ID
                cell->V[idx] = Verts[Verts.Size()-1];
                Share sha = {cell,(size_t)idx};
                Verts[Verts.Size()-1]->Shares.Push (sha);

                // set this lin2 to be disconnected from vert
                lins_to_disconnect.Push (i);
            }
            idx_lin++;
        }
    }
    if (idx_lin<2) throw new Fatal("Mesh::AddPin: Vertex %d (%d) must be connected to at least two linear elements",vert->ID,vert->Tag);

    // disconnect lins from vert
    Array<Share> old_shares(vert->Shares);
    vert->Shares.Resize (vert->Shares.Size()-lins_to_disconnect.Size());
    size_t k = 0;
    for (size_t i=0; i<old_shares.Size(); ++i)
    {
        if (lins_to_disconnect.Find(i)<0)
        {
            vert->Shares[k].C = old_shares[i].C;
            vert->Shares[k].N = old_shares[i].N;
            k++;
        }
    }
}

inline void Generic::SetBryTag (int i, int iEdgeFace, int Tag)
{
    if (Tag<0)
    {
        if (Cells[i]==NULL) throw new Fatal("Generic::SetBryTag: This method must be called after SetCell which allocates Cells");
        Cells[i]->BryTags[iEdgeFace] = Tag;
        TgdCells.XPush (Cells[i]);
    }
}

inline void Generic::FindNeigh ()
{
    // build Bry2Cells map
    Bry2Cells.clear();
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        size_t nverts = Cells[i]->V.Size();
        size_t nbrys  = (NDim==2 ? NVertsToNEdges2D[nverts] : NVertsToNFaces3D[nverts]);
        for (size_t j=0; j<nbrys; ++j)
        {
            BRYKEY(nverts,i,j)
            NeighDat_t neigh_dat(j,Cells[i]);
            Bry2Cells[brykey].Push (neigh_dat);
        }
    }

    // set neighbours information in cells
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        Cells[i]->Neighs.clear();
        size_t nverts = Cells[i]->V.Size();
        size_t nbrys  = (NDim==2 ? NVertsToNEdges2D[nverts] : NVertsToNFaces3D[nverts]);
        for (size_t j=0; j<nbrys; ++j)
        {
            BRYKEY(nverts,i,j)
            Array<NeighDat_t> const & neigh_dat = Bry2Cells[brykey];
            for (size_t k=0; k<neigh_dat.Size(); ++k)
            {
                if (neigh_dat[k].second!=Cells[i])
                {
                    NeighDat_t dat(j, neigh_dat[k].second);
                    Cells[i]->Neighs[brykey] = dat;
                    //std::cout << Cells[i]->ID << "  " << neigh_dat[k].second->ID << "  " << brykey << "\n";
                }
            }
        }
        //std::cout << "\n";
    }
}

inline void Generic::GenO2Verts ()
{
    if (IsShell) throw new Fatal("Generic::GenO2Verts: This mehtod is not ready for Shell meshes yet");

    // generate neighbours map
    if (Bry2Cells.empty()) FindNeigh ();

    //for (BryCell_t::const_iterator p=Bry2Cells.begin(); p!=Bry2Cells.end(); ++p) std::cout << p->first << std::endl;

    // create vertices
    if (NDim==2)
    {
        for (BryCell_t::const_iterator p=Bry2Cells.begin(); p!=Bry2Cells.end(); ++p) // loop over edges
        {
            // add vertex
            Vec3_t mid;
            mid = 0.5*(Verts[boost::get<0>(p->first)]->C + Verts[boost::get<1>(p->first)]->C);
            //std::cout << boost::get<0>(p->first) << ", " << boost::get<1>(p->first) << "   :   " << mid << "   ";
            Verts.Push (new Vertex);
            Verts[Verts.Size()-1]->ID  = Verts.Size()-1;
            Verts[Verts.Size()-1]->Tag = 0;
            Verts[Verts.Size()-1]->C   = mid;

            // set connectivity in cells
            for (size_t i=0; i<p->second.Size(); ++i) // loop over cells sharing this edge
            {
                int idx_edge = p->second[i].first;
                Cell *  cell = p->second[i].second;
                size_t    nv = cell->V.Size();
                int idx_vert = -1;
                if (NVertsToVTKCell2D[nv]==VTK_TRIANGLE || NVertsToVTKCell2D[nv]==VTK_QUAD)
                {
                    // allocate space for nv vertices
                    cell->V.PushN (NULL,nv);
                    idx_vert = NVertsToEdge2D[2*nv][idx_edge][2];
                }
                else if (NVertsToVTKCell2D[nv]==VTK_QUADRATIC_TRIANGLE || NVertsToVTKCell2D[nv]==VTK_QUADRATIC_QUAD)
                {
                    // OK. Space already allocated
                    idx_vert = NVertsToEdge2D[nv][idx_edge][2];
                }
                else throw new Fatal("GenO2Verts::GenO2Verts: NDim==2D. Cell must be of type VTK_TRIANGLE (tri3) or VTK_QUAD (quad4) in order to generate O2 vertices. Number of vertices = %d is invalid",nv);

                // set pointer to vertex
                cell->V[idx_vert] = Verts[Verts.Size()-1];

                // set shares information
                Share sha = {cell,(size_t)idx_vert};
                Verts[Verts.Size()-1]->Shares.Push (sha);
            }
        }
    }
    else throw new Fatal("Generic::GenO2Verts: Method not available for 3D yet");
}

inline void Generic::Erase ()
{
    for (size_t i=0; i<Verts.Size(); ++i) if (Verts[i]!=NULL) delete Verts[i]; // it is only necessary to delete nodes in Verts array
    for (size_t i=0; i<Cells.Size(); ++i) if (Cells[i]!=NULL) delete Cells[i]; // it is only necessary to delete elems in Cells array
    if (Verts   .Size()>0) Verts   .Resize(0);
    if (Cells   .Size()>0) Cells   .Resize(0);
    if (TgdVerts.Size()>0) TgdVerts.Resize(0);
    if (TgdCells.Size()>0) TgdCells.Resize(0);
}

inline void Generic::TagVertex (int Tag, int Idx)
{
    Verts[Idx]->Tag = Tag;
    TgdVerts.XPush (Verts[Idx]);
}

inline int Generic::TagVertex (int Tag, double x, double y, double z, double Tol)
{
    int idx = -1;
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        double d = sqrt(pow(Verts[i]->C(0)-x, 2.0) +
                        pow(Verts[i]->C(1)-y, 2.0) +
             (NDim==3 ? pow(Verts[i]->C(2)-z, 2.0) : 0.0));
        if (d<Tol) { idx=i; break; }
    }
    if (idx<0) throw new Fatal("Generic::TagVertex: Could not find Vertex at (%g,%g,%g)",x,y,z);
    Verts[idx]->Tag = Tag;
    TgdVerts.XPush (Verts[idx]);
    return idx;
}

inline void Generic::TagVertsLine (int Tag, double x0, double y0, double AlpDeg, double Tol)
{
    bool found = false;
    if (NDim==3) throw new Fatal("Generic::TagVertsLine: This method is only available for 2D meshes");
    if ((AlpDeg > -90.0) && (AlpDeg < 90.0))
    {
        double m = tan(AlpDeg*Util::PI/180.0);
        for (size_t i=0; i<Verts.Size(); ++i)
        {
            double x   = Verts[i]->C(0);
            double y   = Verts[i]->C(1);
            double err = y - (y0 + m*(x-x0));
            if (fabs(err)<Tol)
            {
                Verts[i]->Tag = Tag;
                TgdVerts.XPush (Verts[i]);
                found = true;
            }
        }
    }
    else
    {
        for (size_t i=0; i<Verts.Size(); ++i)
        {
            double x   = Verts[i]->C(0);
            double err = x - x0;
            if (fabs(err)<Tol)
            {
                Verts[i]->Tag = Tag;
                TgdVerts.XPush (Verts[i]);
                found = true;
            }
        }
    }
    if (!found) throw new Fatal("Generic::TagVertsLine: Could not find any vertex along line passing through (%g,%g)",x0,y0);
}

inline void Generic::TagSegsLine (int Tag, double x0, double y0, double AlpDeg, double Tol)
{
    bool found = false;
    if (NDim==3) throw new Fatal("Generic::TagSegsLine: This method is only available for 2D meshes");
    if ((AlpDeg > -90.0) && (AlpDeg < 90.0))
    {
        double m = tan(AlpDeg*Util::PI/180.0);
        for (size_t i=0; i<Cells.Size(); ++i)
        {
            size_t nverts = Cells[i]->V.Size();
            size_t nbrys  = NVertsToNEdges2D[nverts];
            for (size_t j=0; j<nbrys; ++j)
            {
                BRYKEY(nverts,i,j)
                double xa    = Verts[vert_a]->C(0);
                double ya    = Verts[vert_a]->C(1);
                double xb    = Verts[vert_b]->C(0);
                double yb    = Verts[vert_b]->C(1);
                double err_a = ya - (y0 + m*(xa-x0));
                double err_b = yb - (y0 + m*(xb-x0));
                if ((fabs(err_a)<Tol) && (fabs(err_b)<Tol))
                {
                    SetBryTag (i, j, Tag);
                    found = true;
                }
            }
        }
    }
    else
    {
        for (size_t i=0; i<Cells.Size(); ++i)
        {
            size_t nverts = Cells[i]->V.Size();
            size_t nbrys  = NVertsToNEdges2D[nverts];
            for (size_t j=0; j<nbrys; ++j)
            {
                BRYKEY(nverts,i,j)
                double xa    = Verts[vert_a]->C(0);
                double xb    = Verts[vert_b]->C(0);
                double err_a = xa-x0;
                double err_b = xb-x0;
                if ((fabs(err_a)<Tol) && (fabs(err_b)<Tol))
                {
                    SetBryTag (i, j, Tag);
                    found = true;
                }
            }
        }
    }
    if (!found) throw new Fatal("Generic::TagSegsLine: Could not find any edge on line passing through (%g,%g)",x0,y0);
}

inline void Generic::TagHorzSegs (int Tag, double y, double xMin, double xMax, double Tol)
{
    if (NDim==3) throw new Fatal("Generic::TagHSeg: This method is only available for 2D meshes");
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        size_t nverts = Cells[i]->V.Size();
        size_t nbrys  = NVertsToNEdges2D[nverts];
        for (size_t j=0; j<nbrys; ++j)
        {
            BRYKEY(nverts,i,j)
            double xa = Verts[vert_a]->C(0);
            double ya = Verts[vert_a]->C(1);
            double xb = Verts[vert_b]->C(0);
            double yb = Verts[vert_b]->C(1);
            if ((fabs(ya-y)<Tol) && (fabs(yb-y)<Tol))
            {
                if ((xa>=xMin) && (xa<=xMax) && (xb>=xMin) && (xb<=xMax))
                {
                    SetBryTag (i, j, Tag);
                }
            }
        }
    }
}

inline void Generic::GroundTags (int LTag, int RTag, int BTag, double Tol)
{
    if (NDim==3) throw new Fatal("Generic::GroundTags: This method is only available for 2D meshes");
    Vec3_t min, max;
    BoundingBox (min, max);
    double left_vert_x   = min(0); // left vertical line
    double right_vert_x  = max(0); // right vertical line
    double bottom_horz_y = min(1); // bottom horizontal line
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        size_t nverts = Cells[i]->V.Size();
        size_t nbrys  = NVertsToNEdges2D[nverts];
        for (size_t j=0; j<nbrys; ++j)
        {
            BRYKEY(nverts,i,j)
            double xa = Verts[vert_a]->C(0);
            double ya = Verts[vert_a]->C(1);
            double xb = Verts[vert_b]->C(0);
            double yb = Verts[vert_b]->C(1);

            // left vertical line
            if ((fabs(xa-left_vert_x)<Tol) && (fabs(xb-left_vert_x)<Tol))
            {
                SetBryTag (i, j, LTag);
            }
            
            // right vertical line
            if ((fabs(xa-right_vert_x)<Tol) && (fabs(xb-right_vert_x)<Tol))
            {
                SetBryTag (i, j, RTag);
            }

            // bottom horizontal line
            if ((fabs(ya-bottom_horz_y)<Tol) && (fabs(yb-bottom_horz_y)<Tol))
            {
                SetBryTag (i, j, BTag);
            }
        }
    }
}

inline void Generic::Check (double Tol) const
{
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        for (size_t j=i+1; j<Verts.Size(); ++j)
        {
            double dist = Norm(Verts[i]->C - Verts[j]->C);
            if (dist<Tol) throw new Fatal("Generic::Check: Found two vertices (%d and %d) with distance=%g smaller than %g",Verts[i]->ID,Verts[j]->ID,dist,Tol);
        }
    }
}

inline void Generic::WriteVTU (char const * FileKey, int VolSurfOrBoth) const
{
    if (IsShell) throw new Fatal("Generic::WriteVTU: This method is not ready for Shell meshes yet");

    // Vol=0, Surf=1, Both=2

    // boundary cells (for plotting bry tags)
    Array<Cell*> bcells;
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        int nverts = Cells[i]->V.Size();
        for (BryTag_t::const_iterator p=Cells[i]->BryTags.begin(); p!=Cells[i]->BryTags.end(); ++p)
        {
            int ibry = p->first;
            int btag = p->second;
            if (btag<0)
            {
                int ibcell  = bcells.Size();
                int nbverts = (NDim==3 ? NVertsToNVertsPerFace3D[nverts] : 2);
                bcells.Push (new Cell);
                bcells[ibcell]->ID  = Cells.Size()+ibcell;
                bcells[ibcell]->Tag = btag;
                for (int j=0; j<nbverts; ++j)
                {
                    bcells[ibcell]->V.Push (new Vertex);
                    bcells[ibcell]->V[j]->ID = (NDim==3 ? Cells[i]->V[NVertsToFace3D[nverts][ibry][j]]->ID : Cells[i]->V[NVertsToEdge2D[nverts][ibry][j]]->ID);
                }
            }
        }
    }

    // shares data
    size_t max_nshares = 0;
    size_t max_nparts  = 0;
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        if (Verts[i]->Shares .Size()>max_nshares) max_nshares = Verts[i]->Shares .Size();
        if (Verts[i]->PartIDs.Size()>max_nparts ) max_nparts  = Verts[i]->PartIDs.Size();
    }
    if (max_nshares<1) throw new Fatal("Mesh::WriteVTU: Max number of shares (%d) is wrong", max_nshares);

    // data
    String fn(FileKey); fn.append(".vtu");
    std::ostringstream oss;
    size_t nn = Verts.Size();                           // number of Nodes
    size_t nc = (VolSurfOrBoth!=1 ? Cells.Size() : 0);  // number of Cells
    size_t nb = bcells.Size();                          // number of boundaries (faces or edges)

    // constants
    size_t          nimax = 40;        // number of integers in a line
    size_t          nfmax =  6;        // number of floats in a line
    Util::NumStream nsflo = Util::_8s; // number format for floats

    // header
    oss << "<?xml version=\"1.0\"?>\n";
    oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    oss << "  <UnstructuredGrid>\n";
    oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << nc+nb << "\">\n";

    // nodes: coordinates
    oss << "      <Points>\n";
    oss << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    size_t k = 0; oss << "        ";
    for (size_t i=0; i<nn; ++i)
    {
        oss << "  " << nsflo <<          Verts[i]->C(0) << " ";
        oss <<         nsflo <<          Verts[i]->C(1) << " ";
        oss <<         nsflo << (NDim==3?Verts[i]->C(2):0.0);
        k++;
        VTU_NEWLINE (i,k,nn,nfmax/3-1,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Points>\n";

    // elements: connectivity, offsets, types
    oss << "      <Cells>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    if (VolSurfOrBoth!=1)
    for (size_t i=0; i<nc; ++i)
    {
        oss << "  ";
        for (size_t j=0; j<Cells[i]->V.Size(); ++j) oss << Cells[i]->V[j]->ID << " ";
        k++;
        VTU_NEWLINE (i,k,nc,nimax/Cells[i]->V.Size(),oss);
    }
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<nb; ++i)
    {
        oss << "  ";
        for (size_t j=0; j<bcells[i]->V.Size(); ++j) oss << bcells[i]->V[j]->ID << " ";
        k++;
        VTU_NEWLINE (i,k,nb,nimax/bcells[i]->V.Size(),oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    size_t offset = 0;
    if (VolSurfOrBoth!=1)
    for (size_t i=0; i<nc; ++i)
    {
        offset += Cells[i]->V.Size();
        oss << (k==0?"  ":" ") << offset;
        k++;
        VTU_NEWLINE (i,k,nc,nimax,oss);
    }
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<nb; ++i)
    {
        offset += bcells[i]->V.Size();
        oss << (k==0?"  ":" ") << offset;
        k++;
        VTU_NEWLINE (i,k,nb,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    if (VolSurfOrBoth!=1)
    for (size_t i=0; i<nc; ++i)
    {
        if (NDim==2) oss << (k==0?"  ":" ") << NVertsToVTKCell2D[Cells[i]->V.Size()];
        else         oss << (k==0?"  ":" ") << NVertsToVTKCell3D[Cells[i]->V.Size()];
        k++;
        VTU_NEWLINE (i,k,nc,nimax,oss);
    }
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<nb; ++i)
    {
        oss << (k==0?"  ":" ") << NVertsToVTKCell2D[bcells[i]->V.Size()];
        k++;
        VTU_NEWLINE (i,k,nb,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Cells>\n";

    // data -- nodes
    oss << "      <PointData Scalars=\"TheScalars\">\n";
    oss << "        <DataArray type=\"Int32\" Name=\"" << "tag" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<nn; ++i)
    {
        oss << (k==0?"  ":" ") << Verts[i]->Tag;
        k++;
        VTU_NEWLINE (i,k,nn,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"" << "shares" << "\" NumberOfComponents=\""<< max_nshares <<"\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<nn; ++i)
    {
        oss << "  ";
        for (size_t j=0; j<max_nshares; ++j)
        {
            if (j<Verts[i]->Shares.Size()) oss << Verts[i]->Shares[j].C->ID << " ";
            else                           oss << -1 << " ";
        }
        k++;
        VTU_NEWLINE (i,k,nn,nimax/max_nshares-1,oss);
    }
    oss << "        </DataArray>\n";
    if (max_nparts>0)
    {
        oss << "        <DataArray type=\"Int32\" Name=\"" << "part" << "\" NumberOfComponents=\""<< max_nparts <<"\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<nn; ++i)
        {
            oss << "  ";
            for (size_t j=0; j<max_nparts; ++j)
            {
                if (j<Verts[i]->PartIDs.Size()) oss << Verts[i]->PartIDs[j] << " ";
                else                            oss << -1 << " ";
            }
            k++;
            VTU_NEWLINE (i,k,nn,nimax/max_nparts-1,oss);
        }
        oss << "        </DataArray>\n";
    }
    oss << "      </PointData>\n";

    // data -- elements
    oss << "      <CellData Scalars=\"TheScalars\">\n";
    oss << "        <DataArray type=\"Int32\" Name=\"" << "tag" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    if (VolSurfOrBoth!=1)
    for (size_t i=0; i<nc; ++i)
    {
        oss << (k==0?"  ":" ") << Cells[i]->Tag;
        k++;
        VTU_NEWLINE (i,k,nc,nimax,oss);
    }
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<nb; ++i)
    {
        oss << (k==0?"  ":" ") << bcells[i]->Tag;
        k++;
        VTU_NEWLINE (i,k,nb,nimax,oss);
    }
    oss << "        </DataArray>\n";
    if (max_nparts>0)
    {
        oss << "        <DataArray type=\"Int32\" Name=\"" << "part" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        if (VolSurfOrBoth!=1)
        for (size_t i=0; i<nc; ++i)
        {
            oss << (k==0?"  ":" ") << Cells[i]->PartID;
            k++;
            VTU_NEWLINE (i,k,nc,nimax,oss);
        }
        oss << "        </DataArray>\n";
    }
    oss << "      </CellData>\n";

    // Bottom
    oss << "    </Piece>\n";
    oss << "  </UnstructuredGrid>\n";
    oss << "</VTKFile>" << std::endl;

    // Write to file
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Generic::WriteMPY (char const * FileKey, bool WithTags, bool WithIDs, bool WithShares, char const * Extra) const
{
    if (IsShell) throw new Fatal("Generic::WriteMPY: This mehtod is not ready for Shell meshes yet");

    // header
    String fn(FileKey); fn.append(".mpy");
    std::ostringstream oss;
    oss << "from msys_drawmesh import *\n\n";
    oss << (*this) << "\n";

    // pins
    oss << "pins = {";
    size_t npins = Pins.size();
    size_t k     = 0;
    for (Pin_t::const_iterator p=Pins.begin(); p!=Pins.end(); ++p)
    {
        oss << p->first->ID << ":[";
        for (size_t i=0; i<p->second.Size(); ++i)
        {
            if (i==p->second.Size()-1) oss << p->second[i]->ID << "]";
            else                       oss << p->second[i]->ID << ",";
        }
        if (k!=npins-1) oss << ",\n        ";
        k++;
    }
    oss << "}\n\n";

    // shares
    /*
    oss << "shares = {";
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        oss << Verts[i]->ID << ":[";
        for (size_t j=0; j<Verts[i]->Shares.Size(); ++j)
        {
            if (j==Verts[i]->Shares.Size()-1) oss << Verts[i]->Shares[j].C->ID << "]";
            else                              oss << Verts[i]->Shares[j].C->ID << ",";
        }
        if (Verts[i]->Shares.Size()>0)
        {
            if (i!=Verts.Size()-1) oss << ",\n          ";
        }
        else
        {
            if (i!=Verts.Size()-1) oss << "],\n          ";
            else                   oss << "]";
        }
    }
    oss << "}\n\n";
    */

    // drawing
    //oss << "d = DrawMesh(V,C,pins,shares)\n";
    oss << "d = DrawMesh(V,C,pins)\n";
    String prms;
    if (WithIDs)    prms.append("True,"); else prms.append("False,");
    if (WithTags)   prms.append("True");  else prms.append("False");
    //if (WithShares) prms.append("True,"); else prms.append("False");
    oss << "d.draw(" << prms << ")\n";
    if (Extra!=NULL) oss << Extra;
    oss << "d.show()\n";
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
    std::cout << "File <" << TERM_CLR_BLUE_H << fn << TERM_RST << "> written\n";
}

inline void Generic::WriteJSON (char const * FileKey) const
{
    // open stream
    String fn(FileKey); fn.append(".msh");
    std::ostringstream oss;
    oss << "{\n";

    // verts
    oss << "    \"verts\" : [\n";
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        if (Verts[i]==NULL) throw new Fatal("Generic::WriteJSON: Vertex %d is not defined",i);
        oss << "        { ";
        oss << "\"id\":"  << Util::_3 << Verts[i]->ID  << ", ";
        oss << "\"tag\":" << Util::_4 << Verts[i]->Tag << ", ";
        oss << "\"c\":[";
        if (NDim==3)
        {
            oss << Util::_20_15 << Verts[i]->C(0) << ", ";
            oss << Util::_20_15 << Verts[i]->C(1) << ", ";
            oss << Util::_20_15 << Verts[i]->C(2);
        }
        else
        {
            oss << Util::_20_15 << Verts[i]->C(0) << ", ";
            oss << Util::_20_15 << Verts[i]->C(1);
        }
        /*
        oss << "], \"shares\":[";
        for (size_t j=0; j<Verts[i]->Shares.Size(); ++j)
        {
            oss << Verts[i]->Shares[j].C->ID;
            if (j!=Verts[i]->Shares.Size()-1) oss << ", ";
        }
        oss << "], \"parts\":[";
        for (size_t j=0; j<Verts[i]->PartIDs.Size(); ++j)
        {
            oss << Verts[i]->PartIDs[j];
            if (j!=Verts[i]->PartIDs.Size()-1) oss << ", ";
        }
        */
        if (i==Verts.Size()-1) oss << "] }\n";
        else                   oss << "] },\n";
    }
    oss << "    ],\n";

    // cells
    oss << "    \"cells\" : [\n";
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        if (Cells[i]==NULL) throw new Fatal("Generic::WriteJSON: Cell %d is not defined",i);
        oss << "        { ";
        oss << "\"id\":"   << Util::_3 << Cells[i]->ID  << ", ";
        oss << "\"tag\":"  << Util::_3 << Cells[i]->Tag << ", ";
        oss << "\"gdim\":" << NDim << ", ";
        oss << "\"part\":" << Cells[i]->PartID << ", ";
        oss << "\"verts\":[";
        for (size_t j=0; j<Cells[i]->V.Size(); ++j)
        {
            oss << Util::_3 << Cells[i]->V[j]->ID;
            if (j!=Cells[i]->V.Size()-1) oss << ", ";
        }
        oss << "], \"ftags\":[";
        size_t nbrys = (NDim==2 ? NVertsToNEdges2D[Cells[i]->V.Size()] : NVertsToNFaces3D[Cells[i]->V.Size()]);
        for (size_t j=0; j<nbrys; ++j)
        {
            oss << Util::_3 << CellGetBryTag(j, Cells[i]);
            if (j!=nbrys-1) oss << ", ";
        }
        /*
        oss << "], \"neigh\":[";
        for (size_t j=0; j<nbrys; ++j)
        {
            oss << CellGetNeigh(j, Cells[i]);
            if (j!=nbrys-1) oss << ", ";
        }
        */
        if (i==Cells.Size()-1) oss << "] }\n";
        else                   oss << "] },\n";
    }
    oss << "    ]\n";

    // write file
    oss << "}\n";
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
    std::cout << "File <" << TERM_CLR_BLUE_H << fn << TERM_RST << "> written\n";
}

inline void Generic::GenGroundSG (Array<double> const & X, Array<double> const & Y, double FootingLx)
{
    // constants
    size_t nx = X.Size();                      // number of columns
    size_t ny = Y.Size();                      // number of rows
    size_t nc = (nx-1)*(ny-1);                 // number of cells
    size_t nv = nx*ny + nx*(ny-1) + (nx-1)*ny; // number of vertices

    // check
    if (NDim!=2) throw new Fatal("Mesh::Generic::GenGroundSG: This method is only available for 2D yet");
    for (size_t j=1; j<nx; ++j) if (X[j]<X[j-1]) throw new Fatal("Mesh::Generic::GenGroundSG: x-coordinates must be increasing");
    for (size_t j=1; j<ny; ++j) if (Y[j]>Y[j-1]) throw new Fatal("Mesh::Generic::GenGroundSG: y-coordinates must be decreasing");

    // set size
    SetSize (nv, nc);

    // generate mesh
    size_t idx_vert = 0;
    size_t idx_cell = 0;
    for (size_t i=0; i<nx-1; ++i)
    {
        // vertices: first column
        if (i==0)
        {
            for (size_t j=0; j<ny; ++j)
            {
                SetVert (idx_vert, 0, X[i], Y[j]);
                idx_vert++;
                if (j!=ny-1) // intermediate nodes
                {
                    double dy = Y[j+1] - Y[j];
                    SetVert (idx_vert, 0, X[i], Y[j]+dy/2.0);
                    idx_vert++;
                }
            }
        }

        // vertices: second column
        double dx = X[i+1] - X[i];
        for (size_t j=0; j<ny; ++j)
        {
            SetVert (idx_vert, 0, X[i]+dx/2.0, Y[j]);
            idx_vert++;
        }

        // vertices: third column
        for (size_t j=0; j<ny; ++j)
        {
            SetVert (idx_vert, 0, X[i]+dx, Y[j]);
            idx_vert++;
            if (j!=ny-1)
            {
                double dy = Y[j+1] - Y[j];
                SetVert (idx_vert, 0, X[i]+dx, Y[j]+dy/2.0);
                idx_vert++;
            }
        }

        // set cells
        for (size_t j=0; j<ny-1; ++j)
        {
            int a = (3*ny-1)*i + 2*j;
            int b = (3*ny-1)*i + (2*ny-1) + j;
            int c = a + 3*ny - 1;
            SetCell (idx_cell, -1, Array<int>(a+2, c+2, c, a,  b+1, c+1, b, a+1));
            if (i==0)    SetBryTag (idx_cell, 3, -10);
            if (i==nx-2) SetBryTag (idx_cell, 1, -20);
            if (j==ny-2) SetBryTag (idx_cell, 0, -30);
            if (j==0) // footing
            {
                if (X[i]+0.95*dx<=FootingLx)
                    SetBryTag (idx_cell, 2, -40);
            }
            idx_cell++;
        }
    }
}

inline void Generic::GenGroundSG (size_t Nx, size_t Ny, double Dx, double Dy)
{
    if (Nx<2) throw new Fatal("Generic::GenGroundSG: Number of columns along x must be at least two. Nx=%d is invalid",Nx);
    if (Ny<2) throw new Fatal("Generic::GenGroundSG: Number of columns along y must be at least two. Ny=%d is invalid",Ny);
    Array<double> X(Nx);
    Array<double> Y(Ny);
    for (size_t i=0; i<Nx; ++i) X[i] =  static_cast<double>(i)*Dx;
    for (size_t i=0; i<Ny; ++i) Y[i] = -static_cast<double>(i)*Dy;
    Y[0] = 0.0;
    GenGroundSG (X, Y);
}

inline void Generic::GenSector (size_t Nr, size_t Nth, double r, double R, double Theta)
{
    // check
    if (Nr<2)    throw new Fatal("Mesh::Generic::GenSector: Nr (number of lines/columns for each radii) must be greater than or equal to 2");
    if (Nth<2)   throw new Fatal("Mesh::Generic::GenSector: Nth (number of lines/rows for each theta) must be greater than or equal to 2");
    if (NDim!=2) throw new Fatal("Mesh::Generic::GenSector: This method is only available for 2D yet");

    // constants
    size_t nc  = (Nr-1)*(Nth-1);                   // number of cells
    size_t nv  = Nr*Nth + Nr*(Nth-1) + (Nr-1)*Nth; // number of vertices
    double dr  = (R-r)/(Nr-1);                     // delta radius
    double dth = Theta/(Nth-1);                    // delta theta

    // set size
    SetSize (nv, nc);

    // generate mesh
    size_t idx_vert = 0;
    size_t idx_cell = 0;
    for (size_t i=0; i<Nr-1; ++i)
    {
        double radius = r+i*dr;

        // vertices: first column
        if (i==0)
        {
            for (size_t j=0; j<Nth; ++j)
            {
                double th = Theta-j*dth;
                double x  = radius*cos(th);
                double y  = radius*sin(th);
                SetVert (idx_vert, (j==0 ? -100 : (j==Nth-1 ? -200 : 0)), x, y);
                idx_vert++;
                if (j!=Nth-1) // intermediate nodes
                {
                    th -= dth/2.0;
                    x   = radius*cos(th);
                    y   = radius*sin(th);
                    SetVert (idx_vert, 0, x, y);
                    idx_vert++;
                }
            }
        }

        // vertices: middle column
        radius += dr/2.0;
        for (size_t j=0; j<Nth; ++j)
        {
            double th = Theta-j*dth;
            double x  = radius*cos(th);
            double y  = radius*sin(th);
            SetVert (idx_vert, (j==0 ? -100 : (j==Nth-1 ? -200 : 0)), x, y);
            idx_vert++;
        }

        // vertices: last column
        radius += dr/2.0;
        for (size_t j=0; j<Nth; ++j)
        {
            double th = Theta-j*dth;
            double x  = radius*cos(th);
            double y  = radius*sin(th);
            SetVert (idx_vert, (j==0 ? -100 : (j==Nth-1 ? -200 : 0)), x, y);
            idx_vert++;
            if (j!=Nth-1)
            {
                th -= dth/2.0;
                x   = radius*cos(th);
                y   = radius*sin(th);
                SetVert (idx_vert, 0, x, y);
                idx_vert++;
            }
        }

        // set cells
        for (size_t j=0; j<Nth-1; ++j)
        {
            int a = (3*Nth-1)*i + 2*j;
            int b = (3*Nth-1)*i + (2*Nth-1) + j;
            int c = a + 3*Nth - 1;
            SetCell (idx_cell, -1, Array<int>(a+2, c+2, c, a,  b+1, c+1, b, a+1));
            if (i==0)     SetBryTag (idx_cell, 3, -10);
            if (i==Nr-2)  SetBryTag (idx_cell, 1, -20);
            if (j==Nth-2) SetBryTag (idx_cell, 0, -30);
            if (j==0)     SetBryTag (idx_cell, 2, -40);
            idx_cell++;
        }
    }
}

inline void Generic::Quad4ToTri3 (bool BackDiagonal)
{
    // copy pointers of old mesh
    Array<Cell*> old_cells(Cells.Size());
    for (size_t i=0; i<old_cells.Size(); ++i) old_cells[i] = Cells[i];
    Cells   .Resize (old_cells.Size()*2);
    TgdCells.Resize (0);
    Cells.SetValues (NULL);

    // erase shares information
    for (size_t i=0; i<Verts.Size(); ++i) Verts[i]->Shares.Resize(0);

    // new mesh
    int icell = 0;
    for (size_t i=0; i<old_cells.Size(); ++i)
    {
        // check
        if (old_cells[i]->V.Size()!=4) throw new Fatal("Mesh::Quad4ToTri3: This method only works for Quad4s (NVerts=%d is invalid)",old_cells[i]->V.Size());

        // indices of vertices of Quad4
        size_t i0 = old_cells[i]->V[0]->ID;
        size_t i1 = old_cells[i]->V[1]->ID;
        size_t i2 = old_cells[i]->V[2]->ID;
        size_t i3 = old_cells[i]->V[3]->ID;

        // new connectivities
        Array<int> con0(3), con1(3);
        if (BackDiagonal)
        {
            con0 = i0, i1, i3;
            con1 = i3, i1, i2;
        }
        else
        {
            con0 = i0, i1, i2;
            con1 = i0, i2, i3;
        }

        // new cells
        SetCell (icell++, old_cells[i]->Tag, con0);
        SetCell (icell++, old_cells[i]->Tag, con1);

        // boundary tags
        if (BackDiagonal)
        {
            for (BryTag_t::const_iterator p=old_cells[i]->BryTags.begin(); p!=old_cells[i]->BryTags.end(); ++p)
            {
                if (p->first==0) SetBryTag (icell-2, /*side*/0, /*tag*/p->second);
                if (p->first==1) SetBryTag (icell-1, /*side*/1, /*tag*/p->second);
                if (p->first==2) SetBryTag (icell-1, /*side*/2, /*tag*/p->second);
                if (p->first==3) SetBryTag (icell-2, /*side*/2, /*tag*/p->second);
            }
        }
        else
        {
            for (BryTag_t::const_iterator p=old_cells[i]->BryTags.begin(); p!=old_cells[i]->BryTags.end(); ++p)
            {
                if (p->first==0) SetBryTag (icell-2, /*side*/0, /*tag*/p->second);
                if (p->first==1) SetBryTag (icell-2, /*side*/1, /*tag*/p->second);
                if (p->first==2) SetBryTag (icell-1, /*side*/1, /*tag*/p->second);
                if (p->first==3) SetBryTag (icell-1, /*side*/2, /*tag*/p->second);
            }
        }
    }

    // delete old mesh
    for (size_t i=0; i<old_cells.Size(); ++i) delete old_cells[i];

    // info
    printf("%s  ------- After Quad4ToTri3 -------%s\n", TERM_GREEN, TERM_RST);
    printf("%s  Num of cells       = %zd%s\n", TERM_CLR2, Cells.Size(), TERM_RST);
    printf("%s  Num of vertices    = %zd%s\n", TERM_CLR2, Verts.Size(), TERM_RST);
}

inline void Generic::Quad8ToTri6 (bool BackDiagonal)
{
    // copy pointers of old mesh
    Array<Cell*> old_cells(Cells.Size());
    for (size_t i=0; i<old_cells.Size(); ++i) old_cells[i] = Cells[i];
    Cells   .Resize (old_cells.Size()*2);
    TgdCells.Resize (0);
    Cells.SetValues (NULL);

    // erase shares information
    for (size_t i=0; i<Verts.Size(); ++i) Verts[i]->Shares.Resize(0);

    // new mesh
    int icell = 0;
    for (size_t i=0; i<old_cells.Size(); ++i)
    {
        // check
        if (old_cells[i]->V.Size()!=8) throw new Fatal("Mesh::Quad8ToTri6: This method only works for Quad8s (NVerts=%d is invalid)",old_cells[i]->V.Size());

        // indices of vertices of Quad8
        size_t i0 = old_cells[i]->V[0]->ID;
        size_t i1 = old_cells[i]->V[1]->ID;
        size_t i2 = old_cells[i]->V[2]->ID;
        size_t i3 = old_cells[i]->V[3]->ID;
        size_t i4 = old_cells[i]->V[4]->ID;
        size_t i5 = old_cells[i]->V[5]->ID;
        size_t i6 = old_cells[i]->V[6]->ID;
        size_t i7 = old_cells[i]->V[7]->ID;

        // new vertex
        Vec3_t xnew = 0.5*(Verts[i0]->C + Verts[i2]->C);
        int iv = PushVert (/*tag*/0, xnew(0), xnew(1), xnew(2));

        // new connectivities
        Array<int> con0(6), con1(6);
        if (BackDiagonal)
        {
            con0 = i0, i1, i3,  i4, iv, i7;
            con1 = i3, i1, i2,  iv, i5, i6;
        }
        else
        {
            con0 = i0, i1, i2,  i4, i5, iv;
            con1 = i0, i2, i3,  iv, i6, i7;
        }

        // new cells
        SetCell (icell++, old_cells[i]->Tag, con0);
        SetCell (icell++, old_cells[i]->Tag, con1);

        // boundary tags
        if (BackDiagonal)
        {
            for (BryTag_t::const_iterator p=old_cells[i]->BryTags.begin(); p!=old_cells[i]->BryTags.end(); ++p)
            {
                if (p->first==0) SetBryTag (icell-2, /*side*/0, /*tag*/p->second);
                if (p->first==1) SetBryTag (icell-1, /*side*/1, /*tag*/p->second);
                if (p->first==2) SetBryTag (icell-1, /*side*/2, /*tag*/p->second);
                if (p->first==3) SetBryTag (icell-2, /*side*/2, /*tag*/p->second);
            }
        }
        else
        {
            for (BryTag_t::const_iterator p=old_cells[i]->BryTags.begin(); p!=old_cells[i]->BryTags.end(); ++p)
            {
                if (p->first==0) SetBryTag (icell-2, /*side*/0, /*tag*/p->second);
                if (p->first==1) SetBryTag (icell-2, /*side*/1, /*tag*/p->second);
                if (p->first==2) SetBryTag (icell-1, /*side*/1, /*tag*/p->second);
                if (p->first==3) SetBryTag (icell-1, /*side*/2, /*tag*/p->second);
            }
        }
    }

    // delete old mesh
    for (size_t i=0; i<old_cells.Size(); ++i) delete old_cells[i];

    // info
    printf("%s  ------- After Quad8ToTri6 -------%s\n", TERM_GREEN, TERM_RST);
    printf("%s  Num of cells       = %zd%s\n", TERM_CLR2, Cells.Size(), TERM_RST);
    printf("%s  Num of vertices    = %zd%s\n", TERM_CLR2, Verts.Size(), TERM_RST);
}

inline void Generic::Tri6ToTri15 ()
{
    // find neighbours
    FindNeigh ();

    // expand connectivity array
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        if (Cells[i]->V.Size()!=6) throw new Fatal("Mesh::Tri6ToTri15: This method only works for Tri6s (NVerts=%d is invalid)",Cells[i]->V.Size());
        for (size_t j=0; j<9; ++j) Cells[i]->V.Push (NULL);
    }

    // data
    Vec_t N(6);   // Tri6 shape functions
    Mat_t C(2,6); // matrix of coordinates
    Vec_t R(9);   // natural r coord of new nodes (6->14)
    Vec_t S(9);   // natural s coord of new nodes (6->14)
    R = 1.0/4.0 , 3.0/4.0 , 3.0/4.0 , 1.0/4.0 , 0.0     , 0.0     , 1.0/4.0 , 1.0/2.0 , 1.0/4.0;
    S = 0.0     , 0.0     , 1.0/4.0 , 3.0/4.0 , 3.0/4.0 , 1.0/4.0 , 1.0/4.0 , 1.0/4.0 , 1.0/2.0;

    // convert
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        // indices of vertices of Tri6
        size_t i0 = Cells[i]->V[0]->ID;
        size_t i1 = Cells[i]->V[1]->ID;
        size_t i2 = Cells[i]->V[2]->ID;
        size_t i3 = Cells[i]->V[3]->ID;
        size_t i4 = Cells[i]->V[4]->ID;
        size_t i5 = Cells[i]->V[5]->ID;

        // matrix of coordinates
        C = Verts[i0]->C(0), Verts[i1]->C(0), Verts[i2]->C(0), Verts[i3]->C(0), Verts[i4]->C(0), Verts[i5]->C(0),
            Verts[i0]->C(1), Verts[i1]->C(1), Verts[i2]->C(1), Verts[i3]->C(1), Verts[i4]->C(1), Verts[i5]->C(1);

        // border nodes
        for (size_t j=6; j<12; ++j)
        {
            if (Cells[i]->V[j]==NULL) // not set yet
            {
                // new vertex
                Tri6Shape (R(j-6), S(j-6), N);
                Vec_t x(C*N);
                int iv = PushVert (/*tag*/0, x(0), x(1));

                // set cell and shares
                Cells[i]->V[j] = Verts[iv];
                Share sha = {Cells[i],j};
                Verts[iv]->Shares.Push (sha);

                // set neighbours
                int idx  = (j-6)%2; // 0 or 1 (first or second node on side)
                int side = (j-6)/2; // side of new vertex
                for (Neighs_t::const_iterator p=Cells[i]->Neighs.begin(); p!=Cells[i]->Neighs.end(); ++p)
                {
                    if (side==p->second.first)
                    {
                        // find J: index of vertex in neighbour
                        Cell * neigh = p->second.second;
                        Neighs_t::const_iterator it = neigh->Neighs.find(p->first);
                        if (it==neigh->Neighs.end()) throw new Fatal("Mesh::Tri6ToTri15: __internal_error__");
                        int neigh_side = it->second.first;
                        int neigh_idx  = (1-idx) + neigh_side*2;
                        int J = 6 + neigh_idx;

                        // set cell and shares
                        neigh->V[J] = Verts[iv];
                        Share neigh_sha = {neigh,(size_t)J};
                        Verts[iv]->Shares.Push (neigh_sha);
                    }
                }
            }
        }

        // centre nodes
        for (size_t j=12; j<15; ++j)
        {
            // new vertex
            Tri6Shape (R(j-6), S(j-6), N);
            Vec_t x(C*N);
            int iv = PushVert (/*tag*/0, x(0), x(1));

            // set cell and shares
            Cells[i]->V[j] = Verts[iv];
            Share sha = {Cells[i],j};
            Verts[iv]->Shares.Push (sha);
        }
    }

    // info
    printf("%s  ------- After Tri6ToTri15 -------%s\n", TERM_GREEN, TERM_RST);
    printf("%s  Num of cells       = %zd%s\n", TERM_CLR2, Cells.Size(), TERM_RST);
    printf("%s  Num of vertices    = %zd%s\n", TERM_CLR2, Verts.Size(), TERM_RST);
}

inline void Generic::Tri3ToTri6 ()
{
    // find neighbours
    FindNeigh ();

    // expand connectivity array
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        if (Cells[i]->V.Size()!=3) throw new Fatal("Mesh::Tri3ToTri6: This method only works for Tri3s (NVerts=%d is invalid)",Cells[i]->V.Size());
        for (size_t j=0; j<3; ++j) Cells[i]->V.Push (NULL);
    }

    // convert
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        // border nodes
        for (size_t j=3; j<6; ++j)
        {
            if (Cells[i]->V[j]==NULL) // not set yet
            {
                // new vertex
                int side = j-3; // side of new vertex
                BRYKEY(/*nverts*/3,i,side)
                double xm = (Verts[vert_a]->C(0)+Verts[vert_b]->C(0))/2.0;
                double ym = (Verts[vert_a]->C(1)+Verts[vert_b]->C(1))/2.0;
                int    iv = PushVert (/*tag*/0, xm, ym);

                // set cell and shares
                Cells[i]->V[j] = Verts[iv];
                Share sha = {Cells[i],j};
                Verts[iv]->Shares.Push (sha);

                // set neighbours
                for (Neighs_t::const_iterator p=Cells[i]->Neighs.begin(); p!=Cells[i]->Neighs.end(); ++p)
                {
                    if (side==p->second.first)
                    {
                        // find J: index of vertex in neighbour
                        Cell * neigh = p->second.second;
                        Neighs_t::const_iterator it = neigh->Neighs.find(p->first);
                        if (it==neigh->Neighs.end()) throw new Fatal("Mesh::Tri3ToTri6: __internal_error__");
                        int neigh_side = it->second.first;
                        int J = 3 + neigh_side;

                        // set cell and shares
                        neigh->V[J] = Verts[iv];
                        Share neigh_sha = {neigh,(size_t)J};
                        Verts[iv]->Shares.Push (neigh_sha);
                    }
                }
            }
        }
    }

    // info
    printf("%s  ------- After Tri3ToTri6 --------%s\n", TERM_GREEN, TERM_RST);
    printf("%s  Num of cells       = %zd%s\n", TERM_CLR2, Cells.Size(), TERM_RST);
    printf("%s  Num of vertices    = %zd%s\n", TERM_CLR2, Verts.Size(), TERM_RST);
}

inline void Generic::Refine ()
{
    size_t old_sz = Cells.Size();
    for (size_t i=0; i<old_sz; ++i)
    {
        if (Cells[i]->V.Size()!=2) throw new Fatal("Mesh::Refine: This method only works for Lin2 cells (NVerts=%d is invalid)",Cells[i]->V.Size());

        // new vertex
        size_t va = Cells[i]->V[0]->ID;
        size_t vb = Cells[i]->V[1]->ID;
        double xm =            (Verts[va]->C(0)+Verts[vb]->C(0))/2.0;
        double ym =            (Verts[va]->C(1)+Verts[vb]->C(1))/2.0;
        double zm = (NDim==3 ? (Verts[va]->C(2)+Verts[vb]->C(2))/2.0 : 0.0);
        int    iv = PushVert (0, xm, ym, zm);

        // new cell
        PushCell (Cells[i]->Tag, Array<int>(va, iv));

        // set old cell
        Cells[i]->V[0] = Verts[iv];
    }

    // info
    printf("%s  ------- After Refine ------------%s\n", TERM_GREEN, TERM_RST);
    printf("%s  Num of cells       = %zd%s\n", TERM_CLR2, Cells.Size(), TERM_RST);
    printf("%s  Num of vertices    = %zd%s\n", TERM_CLR2, Verts.Size(), TERM_RST);
}

inline void Generic::Tri6Shape (double r, double s, Vec_t & N) const
{
    /*    s
     *    ^
     *    |
     *  2
     *    @,(0,1)
     *    | ',
     *    |   ',
     *    |     ',
     *    |       ',   4
     *  5 @          @
     *    |           ',
     *    |             ',
     *    |               ',
     *    |(0,0)            ', (1,0)
     *    @---------@---------@  --> r
     *  0           3          1
     */
    N(0) = 1.0-(r+s)*(3.0-2.0*(r+s));
    N(1) = r*(2.0*r-1.0);
    N(2) = s*(2.0*s-1.0);
    N(3) = 4.0*r*(1.0-(r+s));
    N(4) = 4.0*r*s;
    N(5) = 4.0*s*(1.0-(r+s));
}

#ifdef USE_BOOST_PYTHON
inline void Generic::PyGetVertsEdges (BPy::list & V, BPy::list & E) const
{
    if (NDim==3)
    {
        // vertices
        for (size_t i=0; i<Verts.Size(); ++i)
            V.append (BPy::make_tuple(Verts[i]->C(0), Verts[i]->C(1), Verts[i]->C(2)));

        // edges
        for (size_t i=0; i<Cells.Size(); ++i)
        {
            size_t nv = Cells[i]->V.Size();
            for (size_t j=0; j<NVertsToNEdges3D[nv]; ++j)
            {
                BPy::list pair;
                pair.append (Cells[i]->V[NVertsToEdge3D[nv][j][0]]->ID);
                pair.append (Cells[i]->V[NVertsToEdge3D[nv][j][1]]->ID);
                E.append    (pair);
            }
        }
    }
    else
    {
        // vertices
        for (size_t i=0; i<Verts.Size(); ++i)
            V.append (BPy::make_tuple(Verts[i]->C(0), Verts[i]->C(1), 0.0));

        // edges
        for (size_t i=0; i<Cells.Size(); ++i)
        {
            size_t nv = Cells[i]->V.Size();
            for (size_t j=0; j<NVertsToNEdges2D[nv]; ++j)
            {
                BPy::list pair;
                pair.append (Cells[i]->V[NVertsToEdge2D[nv][j][0]]->ID);
                pair.append (Cells[i]->V[NVertsToEdge2D[nv][j][1]]->ID);
                E.append    (pair);
            }
        }
    }
}
#endif

}; // namespace Mesh

#endif // MECHSYS_MESH_H
