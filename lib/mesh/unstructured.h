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

   This class is a wrapper around the Gmsh library (https://gmsh.info),
   built from source into <mechsys>/pkg/gmsh-4.15.2 by the install script.
   Both the 2D and the 3D generators are handled by Gmsh's built-in CAD
   kernel plus its mesh module.  The public API is kept identical to the
   former Triangle/Tetgen based implementation.
*/

// STL
#include <iostream> // for cout, endl, ostream
#include <sstream>  // for ostringstream
#include <fstream>  // for ofstream
#include <cfloat>   // for DBL_EPSILON
#include <cmath>    // for sqrt, cbrt, fabs
#include <cstdlib>  // for malloc, free
#include <vector>
#include <map>
#include <set>
#include <array>
#include <utility>
#include <functional>
#include <algorithm>

// Gmsh
#include <gmsh.h>

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/draw.h>

namespace Mesh
{


/////////////////////////////////////////////////////////////////////////////////////////// Gmsh helpers /////

// Number of nodes of a Gmsh element type (MSH numbering)
inline int GmshNVerts (int EType)
{
    switch (EType)
    {
        case  1: return  2; // 2-node line
        case  2: return  3; // 3-node triangle
        case  3: return  4; // 4-node quadrangle
        case  4: return  4; // 4-node tetrahedron
        case  5: return  8; // 8-node hexahedron
        case  8: return  3; // 3-node line (O2)
        case  9: return  6; // 6-node triangle (O2)
        case 10: return  9; // 9-node quadrangle (O2)
        case 11: return 10; // 10-node tetrahedron (O2)
        case 12: return 27; // 27-node hexahedron (O2)
        case 16: return  8; // 8-node quadrangle (O2)
        default: return -1;
    }
}

// Initialize the Gmsh API only once
inline void GmshEnsureInit ()
{
    if (!gmsh::isInitialized())
    {
        gmsh::initialize();
        gmsh::option::setNumber("General.Terminal",     1);
        gmsh::option::setNumber("General.Verbosity",    2);
        gmsh::option::setNumber("General.AbortOnError", 1);
    }
}

// Characteristic length equivalent to a maximum element area (2D)
inline double GmshSizeFromArea (double A)
{
    if (A<=0.0) return -1.0;
    return sqrt(4.0*A/sqrt(3.0));
}

// Characteristic length equivalent to a maximum element volume (3D)
inline double GmshSizeFromVolume (double V)
{
    if (V<=0.0) return -1.0;
    return cbrt(6.0*sqrt(2.0)*V);
}

// Signed area of a polygon given by a list of point indices
inline double PolyArea2D (std::vector<int> const & P, std::vector<double> const & Pts)
{
    double a = 0.0;
    size_t n = P.size();
    for (size_t i=0; i<n; ++i)
    {
        size_t j = (i+1)%n;
        a += Pts[P[i]*3]*Pts[P[j]*3+1] - Pts[P[j]*3]*Pts[P[i]*3+1];
    }
    return 0.5*a;
}

// Point in polygon test (ray casting)
inline bool PointInPoly2D (std::vector<int> const & P, std::vector<double> const & Pts, double X, double Y)
{
    bool inside = false;
    size_t n = P.size();
    for (size_t i=0, j=n-1; i<n; j=i++)
    {
        double xi = Pts[P[i]*3], yi = Pts[P[i]*3+1];
        double xj = Pts[P[j]*3], yj = Pts[P[j]*3+1];
        if (((yi>Y)!=(yj>Y)) && (X < (xj-xi)*(Y-yi)/(yj-yi)+xi)) inside = !inside;
    }
    return inside;
}

// Point in closed polyhedron test (ray casting). Orientation independent.
// Polys is a list of polygons (each a list of vertex indices); each polygon is
// fan-triangulated. A ray is cast in an arbitrary direction and the number of
// intersections with the surface is counted (odd => inside).
inline bool PointInPolyhedron (std::vector<std::vector<int> > const & Polys,
                               std::vector<double> const & Pts,
                               double X, double Y, double Z)
{
    Vec3_t p (X,Y,Z);
    Vec3_t d (1.0, 0.1234567, 0.7654321); // arbitrary, not aligned with the axes
    auto V = [&](int i)->Vec3_t { return Vec3_t(Pts[i*3],Pts[i*3+1],Pts[i*3+2]); };
    int count = 0;
    for (size_t q=0; q<Polys.size(); ++q)
    {
        std::vector<int> const & poly = Polys[q];
        if (poly.size()<3) continue;
        for (size_t t=1; t+1<poly.size(); ++t)
        {
            Vec3_t v0 = V(poly[0]);
            Vec3_t v1 = V(poly[t]);
            Vec3_t v2 = V(poly[t+1]);
            Vec3_t e1 = v1-v0;
            Vec3_t e2 = v2-v0;
            Vec3_t pv = cross(d,e2);
            double det = dot(e1,pv);
            if (fabs(det)<1.0e-14) continue;
            double inv = 1.0/det;
            Vec3_t tv = p-v0;
            double u = dot(tv,pv)*inv;
            if (u<0.0 || u>1.0) continue;
            Vec3_t qv = cross(tv,e1);
            double v = dot(d,qv)*inv;
            if (v<0.0 || u+v>1.0) continue;
            double tt = dot(e2,qv)*inv;
            if (tt>1.0e-12) count++;
        }
    }
    return (count%2)==1;
}

// Convex hull in 2D (monotone chain). Returns the hull in counter-clockwise order.
inline void ConvexHull2D (std::vector<double> const & Pts, std::vector<int> & Hull)
{
    size_t n = Pts.size()/3;
    std::vector<int> idx(n);
    for (size_t i=0; i<n; ++i) idx[i] = i;
    std::sort (idx.begin(), idx.end(), [&](int a, int b){
        if (Pts[a*3]!=Pts[b*3]) return Pts[a*3]<Pts[b*3];
        return Pts[a*3+1]<Pts[b*3+1];
    });
    auto cross = [&](int o, int a, int b)->double {
        return (Pts[a*3]-Pts[o*3])*(Pts[b*3+1]-Pts[o*3+1]) - (Pts[a*3+1]-Pts[o*3+1])*(Pts[b*3]-Pts[o*3]);
    };
    std::vector<int> H(2*n);
    size_t k = 0;
    for (size_t i=0; i<n; ++i) // lower hull
    {
        while (k>=2 && cross(H[k-2],H[k-1],idx[i])<=0) k--;
        H[k++] = idx[i];
    }
    for (size_t i=n-1, t=k+1; i>0; --i) // upper hull
    {
        while (k>=t && cross(H[k-2],H[k-1],idx[i-1])<=0) k--;
        H[k++] = idx[i-1];
    }
    H.resize(k-1); // last point == first
    Hull = H;
}

// Convex hull in 3D (incremental algorithm). Returns the outward-oriented triangular faces.
inline void ConvexHull3D (std::vector<double> const & Pts, std::vector<std::array<int,3> > & Faces)
{
    size_t n = Pts.size()/3;
    if (n<4) throw new Fatal("Mesh::ConvexHull3D: At least 4 points are required (%zd given)",n);

    auto P = [&](int i)->Vec3_t { return Vec3_t(Pts[i*3],Pts[i*3+1],Pts[i*3+2]); };
    auto D = [&](int a, int b)->Vec3_t { Vec3_t r = P(a)-P(b); return r; };

    // initial tetrahedron
    int i0 = 0;
    for (size_t i=1; i<n; ++i) if (P(i)(0)<P(i0)(0)) i0 = i;
    int i1 = -1; double dmax = -1.0;
    for (size_t i=0; i<n; ++i) if ((int)i!=i0)
    {
        Vec3_t d = D(i,i0);
        double dd = norm(d);
        if (dd>dmax) { dmax = dd; i1 = i; }
    }
    int i2 = -1; double amax = -1.0;
    for (size_t i=0; i<n; ++i) if ((int)i!=i0 && (int)i!=i1)
    {
        Vec3_t cr = cross(D(i1,i0),D(i,i0));
        double ar = norm(cr);
        if (ar>amax) { amax = ar; i2 = i; }
    }
    if (amax<1.0e-14) throw new Fatal("Mesh::ConvexHull3D: All points are collinear");
    int i3 = -1; double vmax = -1.0;
    Vec3_t nrm = cross(D(i1,i0),D(i2,i0));
    for (size_t i=0; i<n; ++i)
    {
        Vec3_t dd = D(i,i0);
        double v = fabs(dot(nrm,dd));
        if (v>vmax) { vmax = v; i3 = i; }
    }
    if (vmax<1.0e-14) throw new Fatal("Mesh::ConvexHull3D: All points are coplanar");

    // centroid (always inside the hull)
    Vec3_t center(0.0,0.0,0.0);
    for (size_t i=0; i<n; ++i) center += P(i);
    center = center/static_cast<double>(n);

    // append a face and make it outward (the centroid is always inside the hull)
    auto addFaceTo = [&](std::vector<std::array<int,3> > & Fv, int a, int b, int c)
    {
        Vec3_t ba = D(b,a);
        Vec3_t ca = D(c,a);
        Vec3_t nf = cross(ba,ca);
        Vec3_t cm = center-P(a);
        if (dot(nf,cm)>0.0) std::swap(b,c);
        std::array<int,3> f = {{a,b,c}};
        Fv.push_back(f);
    };

    std::vector<std::array<int,3> > F;
    addFaceTo(F,i0,i1,i2);
    addFaceTo(F,i0,i1,i3);
    addFaceTo(F,i0,i2,i3);
    addFaceTo(F,i1,i2,i3);

    for (size_t iq=0; iq<n; ++iq)
    {
        int q = iq;
        if (q==i0 || q==i1 || q==i2 || q==i3) continue;

        // visible faces
        std::vector<bool> vis(F.size(), false);
        bool any = false;
        for (size_t f=0; f<F.size(); ++f)
        {
            Vec3_t ba = D(F[f][1],F[f][0]);
            Vec3_t ca = D(F[f][2],F[f][0]);
            Vec3_t nf = cross(ba,ca);
            Vec3_t qa = D(q,F[f][0]);
            if (dot(nf,qa)>1.0e-12) { vis[f] = true; any = true; }
        }
        if (!any) continue;

        // horizon edges
        std::set<std::pair<int,int> > edges;
        for (size_t f=0; f<F.size(); ++f) if (vis[f])
        {
            edges.insert(std::make_pair(F[f][0],F[f][1]));
            edges.insert(std::make_pair(F[f][1],F[f][2]));
            edges.insert(std::make_pair(F[f][2],F[f][0]));
        }
        std::vector<std::pair<int,int> > horizon;
        for (std::set<std::pair<int,int> >::iterator p=edges.begin(); p!=edges.end(); ++p)
        {
            if (edges.find(std::make_pair(p->second,p->first))==edges.end()) horizon.push_back(*p);
        }

        // remove visible faces and add new ones
        std::vector<std::array<int,3> > NF;
        for (size_t f=0; f<F.size(); ++f) if (!vis[f]) NF.push_back(F[f]);
        for (size_t h=0; h<horizon.size(); ++h) addFaceTo (NF, horizon[h].first, horizon[h].second, q);
        F.swap(NF);
    }

    Faces = F;
}


/////////////////////////////////////////////////////////////////////////////////////////// Unstructured /////

class Unstructured : public virtual Mesh::Generic
{
public:
    // Constructor
    Unstructured (int NDim);

    // Destructor
    ~Unstructured () {}

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

#ifdef USE_BOOST_PYTHON
    void PySet (BPy::dict const & Dat);
#endif

private:
    // Auxiliar read methods
    void ReadGmsh     (std::map<int,int> const & EntReg, int DefTag);
    void ReadBryTags  ();

    // Data
    bool _lst_reg_set; ///< Was the last region (NRegions-1) set ?
    bool _lst_hol_set; ///< Was the last hole (NHoles-1) set ?
    bool _lst_pnt_set; ///< Was the last point (NPoints-1) set ?
    bool _lst_seg_set; ///< Was the last segment (NSegmentsOrFacets-1) set ?
    bool _lst_fac_set; ///< Was the last face (NSegmentsOrFacets-1) set ?

    // Input geometry
    std::vector<double>                             _pnts;    ///< Points: x,y,z concatenated
    std::vector<int>                                _ptag;    ///< Point tags
    std::vector<int>                                _segL;    ///< Segments: left node
    std::vector<int>                                _segR;    ///< Segments: right node
    std::vector<int>                                _segtag;  ///< Segment tags
    std::vector< std::vector< std::vector<int> > >  _facpoly; ///< Facets: list of polygons (vertex indices)
    std::vector<int>                                _factag;  ///< Facet tags
    struct Reg_t { int tag; double size; double x, y, z; };  ///< Region
    struct Hol_t { double x, y, z; };                        ///< Hole
    std::vector<Reg_t>                              _regs;    ///< Regions
    std::vector<Hol_t>                              _hols;    ///< Holes

    // Gmsh entity tags of the input boundary entities (for boundary (edge/face) tags)
    std::vector<int>                                _seg_curve; ///< Curve tag of each input segment (2D)
    std::vector<int>                                _fac_surf;  ///< Surface tag of each input facet (3D)

    // Node tag => vertex index
    std::map<size_t,int>                            _node2vert;

    // Points to be embedded in the interior of the domain (used by Delaunay)
    std::vector<int>                                _embed_pts;
};


/////////////////////////////////////////////////////////////////////////////////////////// PLC: Implementation /////

inline Unstructured::Unstructured (int NDim)
    : Mesh::Generic (NDim),
      _lst_reg_set  (false),
      _lst_hol_set  (false),
      _lst_pnt_set  (false),
      _lst_seg_set  (false),
      _lst_fac_set  (false)
{
    GmshEnsureInit ();
}

inline void Unstructured::Set (size_t NPoints, size_t NSegmentsOrFacets, size_t NRegions, size_t NHoles)
{
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

    // points
    _pnts.assign (NPoints*3, 0.0);
    _ptag.assign (NPoints, 0);
    _embed_pts.clear();

    // regions and holes
    _regs.assign (NRegions, Reg_t{0,-1.0,0.0,0.0,0.0});
    _hols.assign (NHoles,   Hol_t{0.0,0.0,0.0});

    if (NDim==2)
    {
        _segL.assign   (NSegmentsOrFacets, -1);
        _segR.assign   (NSegmentsOrFacets, -1);
        _segtag.assign (NSegmentsOrFacets,  0);
        _seg_curve.assign (NSegmentsOrFacets, -1);
    }
    else if (NDim==3)
    {
        _facpoly.assign (NSegmentsOrFacets, std::vector< std::vector<int> >());
        _factag.assign  (NSegmentsOrFacets, 0);
        _fac_surf.assign(NSegmentsOrFacets, -1);
    }
    else throw new Fatal("Unstructured::Set: NDim must be either 2 or 3. NDim==%d is invalid",NDim);
}

inline void Unstructured::SetReg (size_t iReg, int RTag, double MaxAreaOrVolume, double X, double Y, double Z)
{
    _regs[iReg].tag  = RTag;
    _regs[iReg].size = MaxAreaOrVolume;
    _regs[iReg].x    = X;
    _regs[iReg].y    = Y;
    _regs[iReg].z    = Z;
    if ((int)iReg==static_cast<int>(_regs.size())-1) _lst_reg_set = true;
}

inline void Unstructured::SetHol (size_t iHol, double X, double Y, double Z)
{
    _hols[iHol].x = X;
    _hols[iHol].y = Y;
    _hols[iHol].z = Z;
    if ((int)iHol==static_cast<int>(_hols.size())-1) _lst_hol_set = true;
}

inline void Unstructured::SetPnt (size_t iPnt, int PTag, double X, double Y, double Z)
{
    _pnts[iPnt*3  ] = X;
    _pnts[iPnt*3+1] = Y;
    _pnts[iPnt*3+2] = Z;
    _ptag[iPnt]     = PTag;
    if ((int)iPnt==static_cast<int>(_ptag.size())-1) _lst_pnt_set = true;
}

inline void Unstructured::SetSeg (size_t iSeg, int ETag, int L, int R)
{
    if (NDim==3) throw new Fatal("Unstructured::SetSeg: This method must be called for 2D meshes only");
    _segL[iSeg]   = L;
    _segR[iSeg]   = R;
    _segtag[iSeg] = ETag;
    if ((int)iSeg==static_cast<int>(_segL.size())-1) _lst_seg_set = true;
}

inline void Unstructured::SetFac (size_t iFac, int FTag, Array<int> const & VertsOnFace)
{
    if (NDim==2) throw new Fatal("Unstructured::SetSeg: This method must be called for 3D meshes only");
    std::vector<int> poly (VertsOnFace.Size());
    for (size_t j=0; j<VertsOnFace.Size(); ++j) poly[j] = VertsOnFace[j];
    _facpoly[iFac].assign (1, poly);
    _factag[iFac] = FTag;
    if ((int)iFac==static_cast<int>(_facpoly.size())-1) _lst_fac_set = true;
}

inline void Unstructured::SetFac (size_t iFac, int FTag, Array<int> const & Polygon1, Array<int> const & Polygon2)
{
    if (NDim==2) throw new Fatal("Unstructured::SetSeg: This method must be called for 3D meshes only");
    _facpoly[iFac].resize (2);
    _facpoly[iFac][0].resize (Polygon1.Size());
    for (size_t j=0; j<Polygon1.Size(); ++j) _facpoly[iFac][0][j] = Polygon1[j];
    _facpoly[iFac][1].resize (Polygon2.Size());
    for (size_t j=0; j<Polygon2.Size(); ++j) _facpoly[iFac][1][j] = Polygon2[j];
    _factag[iFac] = FTag;
    if ((int)iFac==static_cast<int>(_facpoly.size())-1) _lst_fac_set = true;
}

inline bool Unstructured::IsSet () const
{
    if (NDim==2)
    {
        bool hol_ok = (_hols.size()>0 ? _lst_hol_set : true);
        return (_lst_reg_set && hol_ok && _lst_pnt_set && _lst_seg_set);
    }
    if (NDim==3)
    {
        bool hol_ok = (_hols.size()>0 ? _lst_hol_set : true);
        return (_lst_reg_set && hol_ok && _lst_pnt_set && _lst_fac_set);
    }
    return false;
}

inline void Unstructured::Generate (bool O2, double GlobalMaxArea, bool Quiet, double MinAngle)
{
    // check
    if (!IsSet()) throw new Fatal("Unstructured::Generate: Please, set the input data (regions,points,segments/facets) first.");

    // info
    Util::Stopwatch stopwatch(/*activated*/WithInfo);

    // init gmsh
    GmshEnsureInit ();
    gmsh::clear();
    gmsh::model::add ("mechsys");
    gmsh::option::setNumber ("General.Terminal", Quiet ? 0 : 1);
    gmsh::option::setNumber ("Mesh.ElementOrder", O2 ? 2 : 1);
    gmsh::option::setNumber ("Mesh.MeshSizeExtendFromBoundary", 1);
    gmsh::option::setNumber ("Mesh.MeshSizeFromPoints",         1);
    gmsh::option::setNumber ("Mesh.MeshSizeFromCurvature",      0);
    if (MinAngle>0)
    {
        gmsh::option::setNumber ("Mesh.Algorithm",   6); // Frontal-Delaunay for quads/quality
        gmsh::option::setNumber ("Mesh.Algorithm3D", 4); // Frontal
    }

    // characteristic length
    double lc = -1.0;
    if (GlobalMaxArea>0) lc = (NDim==2 ? GmshSizeFromArea(GlobalMaxArea) : GmshSizeFromVolume(GlobalMaxArea));
    for (size_t i=0; i<_regs.size(); ++i)
    {
        if (_regs[i].size>0)
        {
            double rl = (NDim==2 ? GmshSizeFromArea(_regs[i].size) : GmshSizeFromVolume(_regs[i].size));
            if (lc<0.0 || rl<lc) lc = rl;
        }
    }
    if (lc<=0.0)
    {
        // No size constraint given: emulate the old (unrefined) Triangle/Tetgen behaviour
        // by asking Gmsh for a mesh without inserting additional points.
        double xmin=DBL_MAX, ymin=DBL_MAX, zmin=DBL_MAX;
        double xmax=-DBL_MAX, ymax=-DBL_MAX, zmax=-DBL_MAX;
        for (size_t i=0; i<_ptag.size(); ++i)
        {
            xmin=std::min(xmin,_pnts[i*3  ]); xmax=std::max(xmax,_pnts[i*3  ]);
            ymin=std::min(ymin,_pnts[i*3+1]); ymax=std::max(ymax,_pnts[i*3+1]);
            zmin=std::min(zmin,_pnts[i*3+2]); zmax=std::max(zmax,_pnts[i*3+2]);
        }
        double diag = sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));
        lc = 1.0e3*(diag>0.0 ? diag : 1.0);
    }
    gmsh::option::setNumber ("Mesh.MeshSizeMax", lc);

    // default region tag
    int def_tag = (_regs.size()>0 ? _regs[0].tag : 0);

    // add points
    size_t np = _ptag.size();
    for (size_t i=0; i<np; ++i)
        gmsh::model::geo::addPoint (_pnts[i*3], _pnts[i*3+1], _pnts[i*3+2], lc, static_cast<int>(i)+1);

    // entity tag => region tag
    std::map<int,int> ent_reg;
    std::vector<int>  vol_tags; // 3D volume entity tags

    if (NDim==2)
    {
        size_t ns = _segL.size();

        // add curves (lines)
        for (size_t k=0; k<ns; ++k)
            _seg_curve[k] = gmsh::model::geo::addLine (_segL[k]+1, _segR[k]+1, static_cast<int>(k)+1);

        // trace closed loops of the PSLG
        std::map<int, std::vector<int> > inc;
        for (size_t k=0; k<ns; ++k)
        {
            inc[_segL[k]].push_back (static_cast<int>(k));
            inc[_segR[k]].push_back (static_cast<int>(k));
        }
        struct Loop_t { std::vector<int> pts; std::vector<int> tags; double area; };
        std::vector<Loop_t> loops;
        std::vector<bool> used (ns, false);
        for (size_t k0=0; k0<ns; ++k0)
        {
            if (used[k0]) continue;
            int start = _segL[k0];
            int node  = _segR[k0];
            Loop_t L;
            L.pts.push_back (start);
            L.pts.push_back (node);
            L.tags.push_back (static_cast<int>(k0)+1); // forward (line was created L->R)
            used[k0] = true;
            bool closed = (node==start);
            for (size_t guard=0; guard<=ns && !closed; ++guard)
            {
                int nxt = -1; bool fwd = true;
                std::vector<int> const & v = inc[node];
                for (size_t a=0; a<v.size(); ++a)
                {
                    int kk = v[a];
                    if (used[kk]) continue;
                    nxt = kk;
                    fwd = (_segL[kk]==node);
                    break;
                }
                if (nxt<0) break;
                used[nxt] = true;
                node = (fwd ? _segR[nxt] : _segL[nxt]);
                L.pts.push_back (node);
                L.tags.push_back (fwd ? (nxt+1) : -(nxt+1));
                if (node==start) closed = true;
            }
            if (closed)
            {
                L.area = PolyArea2D (L.pts, _pnts);
                if (L.area<0.0) // make it counter-clockwise
                {
                    std::reverse (L.tags.begin(), L.tags.end());
                    for (size_t t=0; t<L.tags.size(); ++t) L.tags[t] = -L.tags[t];
                    L.area = -L.area;
                }
                loops.push_back (L);
            }
        }
        if (loops.size()<1) throw new Fatal("Unstructured::Generate: Could not trace any closed loop from the given segments");

        // classify holes: a hole loop contains a hole point and no region point
        std::vector<bool> is_hole (loops.size(), false);
        for (size_t i=0; i<loops.size(); ++i)
        {
            bool has_hol = false;
            for (size_t h=0; h<_hols.size(); ++h)
                if (PointInPoly2D (loops[i].pts, _pnts, _hols[h].x, _hols[h].y)) { has_hol = true; break; }
            bool has_reg = false;
            for (size_t r=0; r<_regs.size(); ++r)
                if (PointInPoly2D (loops[i].pts, _pnts, _regs[r].x, _regs[r].y)) { has_reg = true; break; }
            is_hole[i] = (has_hol && !has_reg);
        }

        // create plane surfaces (exterior + holes)
        std::vector<int> surf_tags;
        for (size_t i=0; i<loops.size(); ++i)
        {
            if (is_hole[i]) continue;
            std::vector<int> wires;
            wires.push_back (gmsh::model::geo::addCurveLoop (loops[i].tags));
            for (size_t j=0; j<loops.size(); ++j)
            {
                if (!is_hole[j]) continue;
                if (loops[j].pts.size()<1) continue;
                if (PointInPoly2D (loops[i].pts, _pnts, _pnts[loops[j].pts[0]*3], _pnts[loops[j].pts[0]*3+1]))
                    wires.push_back (gmsh::model::geo::addCurveLoop (loops[j].tags));
            }
            int s = gmsh::model::geo::addPlaneSurface (wires);
            surf_tags.push_back (s);

            // region tag: first region point inside this loop
            int rt = def_tag;
            for (size_t r=0; r<_regs.size(); ++r)
                if (PointInPoly2D (loops[i].pts, _pnts, _regs[r].x, _regs[r].y)) { rt = _regs[r].tag; break; }
            ent_reg[s] = rt;
        }
        if (surf_tags.size()<1) throw new Fatal("Unstructured::Generate: Could not create any surface");
    }
    else // NDim==3
    {
        size_t nf = _facpoly.size();

        // unique edges shared by facets
        std::map<std::pair<int,int>,int> edge_curve;
        std::map<int,std::vector<int> >   curve_facets;

        for (size_t i=0; i<nf; ++i)
        {
            std::vector<int> wires;
            for (size_t p=0; p<_facpoly[i].size(); ++p)
            {
                std::vector<int> const & poly = _facpoly[i][p];
                std::vector<int> wire;
                size_t nvp = poly.size();
                for (size_t j=0; j<nvp; ++j)
                {
                    int a = poly[j];
                    int b = poly[(j+1)%nvp];
                    int lo = std::min(a,b);
                    int hi = std::max(a,b);
                    std::pair<int,int> key (lo,hi);
                    int ct;
                    std::map<std::pair<int,int>,int>::iterator it = edge_curve.find (key);
                    if (it==edge_curve.end())
                    {
                        ct = gmsh::model::geo::addLine (lo+1, hi+1);
                        edge_curve[key] = ct;
                        curve_facets[ct] = std::vector<int>();
                    }
                    else ct = it->second;
                    curve_facets[ct].push_back (static_cast<int>(i));
                    wire.push_back (a< b ? ct : -ct); // line goes from lo to hi
                }
                if (wire.size()>0) wires.push_back (gmsh::model::geo::addCurveLoop (wire));
            }
            if (wires.size()<1) throw new Fatal("Unstructured::Generate: Facet %zd has no vertices",i);
            _fac_surf[i] = gmsh::model::geo::addPlaneSurface (wires);
        }

        // connected components of facets (surfaces sharing a curve)
        std::vector<int> par (nf);
        for (size_t i=0; i<nf; ++i) par[i] = i;
        std::function<int(int)> find = [&](int x)->int { while (par[x]!=x) { par[x]=par[par[x]]; x=par[x]; } return x; };
        for (std::map<int,std::vector<int> >::iterator it=curve_facets.begin(); it!=curve_facets.end(); ++it)
        {
            std::vector<int> const & fs = it->second;
            for (size_t a=1; a<fs.size(); ++a) par[find(fs[a])] = find(fs[0]);
        }
        std::map<int,std::vector<int> > comps;
        for (size_t i=0; i<nf; ++i) comps[find(static_cast<int>(i))].push_back (static_cast<int>(i));

        // create one surface loop per component, and keep its polygons/root
        std::map<int,int>                                root2shell;
        std::map<int,std::vector<std::vector<int> > >    comp_polys;
        for (std::map<int,std::vector<int> >::iterator it=comps.begin(); it!=comps.end(); ++it)
        {
            std::vector<int> surfaces;
            for (size_t a=0; a<it->second.size(); ++a)
            {
                surfaces.push_back (_fac_surf[it->second[a]]);
                for (size_t p=0; p<_facpoly[it->second[a]].size(); ++p)
                    comp_polys[it->first].push_back (_facpoly[it->second[a]][p]);
            }
            root2shell[it->first] = gmsh::model::geo::addSurfaceLoop (surfaces);
        }

        // classify each component: a hole contains a hole point and no region point
        std::map<int,bool> is_hole;
        for (std::map<int,std::vector<int> >::iterator it=comps.begin(); it!=comps.end(); ++it)
        {
            std::vector<std::vector<int> > const & polys = comp_polys[it->first];
            bool has_reg = false;
            for (size_t r=0; r<_regs.size() && !has_reg; ++r)
                if (PointInPolyhedron (polys,_pnts,_regs[r].x,_regs[r].y,_regs[r].z)) has_reg = true;
            bool has_hol = false;
            for (size_t h=0; h<_hols.size() && !has_hol; ++h)
                if (PointInPolyhedron (polys,_pnts,_hols[h].x,_hols[h].y,_hols[h].z)) has_hol = true;
            is_hole[it->first] = (has_hol && !has_reg);
        }

        // create a volume per region component, punching the holes it contains
        for (std::map<int,std::vector<int> >::iterator it=comps.begin(); it!=comps.end(); ++it)
        {
            if (is_hole[it->first]) continue;
            std::vector<int> shells;
            shells.push_back (root2shell[it->first]);
            for (std::map<int,std::vector<int> >::iterator jt=comps.begin(); jt!=comps.end(); ++jt)
            {
                if (!is_hole[jt->first]) continue;
                bool inside = false;
                for (size_t h=0; h<_hols.size() && !inside; ++h)
                {
                    if (PointInPolyhedron (comp_polys[jt->first],_pnts,_hols[h].x,_hols[h].y,_hols[h].z) &&
                        PointInPolyhedron (comp_polys[it->first],_pnts,_hols[h].x,_hols[h].y,_hols[h].z)) inside = true;
                }
                if (inside) shells.push_back (root2shell[jt->first]);
            }
            int vol = gmsh::model::geo::addVolume (shells);
            vol_tags.push_back (vol);
            ent_reg[vol] = def_tag;
        }
    }

    // synchronize the geometry
    gmsh::model::geo::synchronize ();

    // 3D: assign region tags by checking which volume contains each region point
    if (NDim==3 && _regs.size()>1)
    {
        for (size_t i=0; i<vol_tags.size(); ++i)
        {
            int vol = vol_tags[i];
            for (size_t r=0; r<_regs.size(); ++r)
            {
                std::vector<double> c (3);
                c[0]=_regs[r].x; c[1]=_regs[r].y; c[2]=_regs[r].z;
                if (gmsh::model::isInside (3, vol, c)>0) { ent_reg[vol] = _regs[r].tag; break; }
            }
        }
    }

    // embed interior points (Delaunay)
    if (_embed_pts.size()>0)
    {
        // do not let short boundary edges force a refined mesh
        gmsh::option::setNumber ("Mesh.MeshSizeExtendFromBoundary", 0);
        gmsh::option::setNumber ("Mesh.MeshSizeFromPoints",         0);
        gmsh::option::setNumber ("Mesh.MeshSizeFromCurvature",      0);

        std::vector<int> ptags (_embed_pts.size());
        for (size_t i=0; i<_embed_pts.size(); ++i) ptags[i] = _embed_pts[i]+1;
        gmsh::vectorpair ents;
        gmsh::model::getEntities (ents, NDim);
        for (size_t ie=0; ie<ents.size(); ++ie)
            gmsh::model::mesh::embed (0, ptags, NDim, ents[ie].second);
    }

    // generate mesh
    gmsh::model::mesh::generate (NDim);

    // read mesh
    ReadGmsh (ent_reg, def_tag);

    // tag input points
    for (size_t i=0; i<_ptag.size(); ++i)
    {
        if (_ptag[i]==0) continue;
        std::vector<std::size_t> ntags;
        std::vector<double> coord, param;
        gmsh::model::mesh::getNodes (ntags, coord, param, 0, static_cast<int>(i)+1, false, false);
        for (size_t j=0; j<ntags.size(); ++j)
        {
            std::map<size_t,int>::iterator it = _node2vert.find (ntags[j]);
            if (it==_node2vert.end()) continue;
            if (Verts[it->second]->Tag==0)
            {
                Verts[it->second]->Tag = _ptag[i];
                TgdVerts.Push (Verts[it->second]);
            }
        }
    }

    // boundary (edge/face) tags
    ReadBryTags ();

    // check
    if (Verts.Size()<1) throw new Fatal("Unstructured::Generate: Failed with %d vertices and %d cells", Verts.Size(), Cells.Size());
    if (Cells.Size()<1) throw new Fatal("Unstructured::Generate: Failed with %d vertices and %d cells", Verts.Size(), Cells.Size());

    // info
    if (WithInfo)
    {
        printf("\n%s--- Unstructured Mesh Generation (Gmsh) --- %dD --- O%d ---------------------%s\n",TERM_CLR1,NDim,(O2?2:1),TERM_RST);
        printf("%s  Element order      = %s%d%s\n", TERM_CLR2, TERM_CLR4, (O2?2:1), TERM_RST);
        printf("%s  Num of cells       = %zd%s\n", TERM_CLR2, Cells.Size(), TERM_RST);
        printf("%s  Num of vertices    = %zd%s\n", TERM_CLR2, Verts.Size(), TERM_RST);
    }
}

inline void Unstructured::ReadGmsh (std::map<int,int> const & EntReg, int DefTag)
{
    // clear previous mesh
    Erase ();

    // nodes
    std::vector<std::size_t> ntags;
    std::vector<double> coord, param;
    gmsh::model::mesh::getNodes (ntags, coord, param, -1, -1, false, false);
    Verts.Resize (ntags.size());
    _node2vert.clear();
    for (size_t i=0; i<ntags.size(); ++i)
    {
        Verts[i]      = new Vertex;
        Verts[i]->ID  = i;
        Verts[i]->Tag = 0;
        Verts[i]->C   = coord[i*3], coord[i*3+1], coord[i*3+2];
        _node2vert[ntags[i]] = static_cast<int>(i);
    }
    TgdVerts.Resize (0);

    // cells, grouped by entity (to recover the region tag)
    gmsh::vectorpair ents;
    gmsh::model::getEntities (ents, NDim);
    for (size_t ie=0; ie<ents.size(); ++ie)
    {
        int etag = ents[ie].second;
        int reg  = DefTag;
        std::map<int,int>::const_iterator it = EntReg.find (etag);
        if (it!=EntReg.end()) reg = it->second;

        std::vector<int> types;
        std::vector<std::vector<std::size_t> > etags, nodetags;
        gmsh::model::mesh::getElements (types, etags, nodetags, NDim, etag);
        for (size_t t=0; t<types.size(); ++t)
        {
            int nv = GmshNVerts (types[t]);
            if (nv<0) continue;
            size_t ne = etags[t].size();
            for (size_t e=0; e<ne; ++e)
            {
                Cells.Push (static_cast<Cell*>(0));
                size_t ic = Cells.Size()-1;
                Cells[ic]        = new Cell;
                Cells[ic]->ID    = ic;
                Cells[ic]->Tag   = reg;
                Cells[ic]->PartID= 0;
                Cells[ic]->V.Resize (nv);
                for (int j=0; j<nv; ++j)
                {
                    size_t nt = nodetags[t][e*nv+j];
                    int    iv = _node2vert[nt];
                    Share sha = {Cells[ic],static_cast<size_t>(j)};
                    Cells[ic]->V[j] = Verts[iv];
                    Verts[iv]->Shares.Push (sha);
                }
            }
        }
    }
}

inline void Unstructured::ReadBryTags ()
{
    // clear previous tagged cells
    TgdCells.Resize (0);

    if (NDim==2)
    {
        // map: sorted pair of vertex indices => (cell, local edge)
        std::map<std::pair<int,int>, std::pair<int,int> > emap;
        for (size_t ic=0; ic<Cells.Size(); ++ic)
        {
            size_t nv = Cells[ic]->V.Size();
            size_t ne = NVertsToNEdges2D[nv];
            for (size_t f=0; f<ne; ++f)
            {
                int a = Cells[ic]->V[NVertsToEdge2D[nv][f][0]]->ID;
                int b = Cells[ic]->V[NVertsToEdge2D[nv][f][1]]->ID;
                if (a>b) std::swap(a,b);
                emap[std::make_pair(a,b)] = std::make_pair(static_cast<int>(ic),static_cast<int>(f));
            }
        }

        for (size_t k=0; k<_seg_curve.size(); ++k)
        {
            if (_segtag[k]==0 || _seg_curve[k]<0) continue;
            std::vector<int> types;
            std::vector<std::vector<std::size_t> > etags, nodetags;
            gmsh::model::mesh::getElements (types, etags, nodetags, 1, _seg_curve[k]);
            for (size_t t=0; t<types.size(); ++t)
            {
                int nv = GmshNVerts (types[t]);
                if (nv<2) continue;
                size_t ne = etags[t].size();
                for (size_t e=0; e<ne; ++e)
                {
                    int a = _node2vert[nodetags[t][e*nv  ]];
                    int b = _node2vert[nodetags[t][e*nv+1]];
                    if (a>b) std::swap(a,b);
                    std::map<std::pair<int,int>, std::pair<int,int> >::iterator it = emap.find (std::make_pair(a,b));
                    if (it==emap.end()) continue;
                    Cells[it->second.first]->BryTags[it->second.second] = _segtag[k];
                }
            }
        }
    }
    else // NDim==3
    {
        // map: sorted triple of vertex indices => (cell, local face)
        std::map<std::array<int,3>, std::pair<int,int> > fmap;
        for (size_t ic=0; ic<Cells.Size(); ++ic)
        {
            size_t nv = Cells[ic]->V.Size();
            size_t nf = NVertsToNFaces3D[nv];
            size_t npf= NVertsToNVertsPerFace3D[nv];
            for (size_t f=0; f<nf; ++f)
            {
                std::array<int,3> key = {{-1,-1,-1}};
                for (size_t q=0; q<npf; ++q)
                    key[q] = Cells[ic]->V[NVertsToFace3D[nv][f][q]]->ID;
                std::sort (key.begin(), key.end());
                fmap[key] = std::make_pair(static_cast<int>(ic),static_cast<int>(f));
            }
        }

        for (size_t k=0; k<_fac_surf.size(); ++k)
        {
            if (_factag[k]==0 || _fac_surf[k]<0) continue;
            std::vector<int> types;
            std::vector<std::vector<std::size_t> > etags, nodetags;
            gmsh::model::mesh::getElements (types, etags, nodetags, 2, _fac_surf[k]);
            for (size_t t=0; t<types.size(); ++t)
            {
                int nv = GmshNVerts (types[t]);
                if (nv<3) continue;
                size_t ne = etags[t].size();
                for (size_t e=0; e<ne; ++e)
                {
                    std::array<int,3> key = {{ _node2vert[nodetags[t][e*nv  ]],
                                               _node2vert[nodetags[t][e*nv+1]],
                                               _node2vert[nodetags[t][e*nv+2]] }};
                    std::sort (key.begin(), key.end());
                    std::map<std::array<int,3>, std::pair<int,int> >::iterator it = fmap.find (key);
                    if (it==fmap.end()) continue;
                    Cells[it->second.first]->BryTags[it->second.second] = _factag[k];
                }
            }
        }
    }

    // collect tagged cells
    for (size_t ic=0; ic<Cells.Size(); ++ic)
        if (Cells[ic]->BryTags.size()>0) TgdCells.Push (Cells[ic]);
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

        // points
        oss << "pts = [";
        for (size_t i=0; i<_ptag.size(); ++i)
        {
            oss << "[" << _pnts[i*3] << "," << _pnts[i*3+1] << "," << _pnts[i*3+2] << "]";
            if (i==_ptag.size()-1) oss << "]\n";
            else                   oss << ",\n       ";
        }
        oss << "\n";

        // edges
        oss << "edg = [";
        if (NDim==2)
        {
            for (size_t i=0; i<_segL.size(); ++i)
            {
                oss << "[" << _segL[i] << "," << _segR[i] << "]";
                if (i==_segL.size()-1) oss << "]\n";
                else                   oss << ",\n       ";
            }
        }
        else
        {
            bool first = true;
            for (size_t i=0; i<_facpoly.size(); ++i)
            {
                for (size_t p=0; p<_facpoly[i].size(); ++p)
                {
                    std::vector<int> const & poly = _facpoly[i][p];
                    for (size_t k=1; k<poly.size(); ++k)
                    {
                        if (!first) oss << ",\n       ";
                        oss << "[" << poly[k-1] << "," << poly[k] << "]";
                        first = false;
                    }
                    if (poly.size()>1)
                    {
                        if (!first) oss << ",\n       ";
                        oss << "[" << poly[poly.size()-1] << "," << poly[0] << "]";
                        first = false;
                    }
                }
            }
            oss << "]\n";
        }
        oss << "\n";

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
        for (size_t i=0; i<_segL.size(); ++i)
        {
            int I = _segL[i];
            int J = _segR[i];
            oss << "dat.append((PH.MOVETO, (" << _pnts[I*3] << "," << _pnts[I*3+1] << ")))\n";
            oss << "dat.append((PH.LINETO, (" << _pnts[J*3] << "," << _pnts[J*3+1] << ")))\n";
        }
        oss << "\n";

        // draw edges
        MPL::AddPatch (oss);

        // draw tags
        oss << "# draw tags\n";
        for (size_t i=0; i<_ptag.size(); ++i)
        {
            if (_ptag[i]<0) oss << "ax.text(" << _pnts[i*3] << "," << _pnts[i*3+1] << ", " << _ptag[i] << ", ha='center', va='center', fontsize=14, backgroundcolor=lyellow)\n";
        }
        for (size_t i=0; i<_segL.size(); ++i)
        {
            if (_segtag[i]>=0) continue;
            int    I  = _segL[i];
            int    J  = _segR[i];
            double x0 = _pnts[I*3];
            double y0 = _pnts[I*3+1];
            double x1 = _pnts[J*3];
            double y1 = _pnts[J*3+1];
            double xm = (x0+x1)/2.0;
            double ym = (y0+y1)/2.0;
            oss << "ax.text(" << xm << "," << ym << ", " << _segtag[i] << ", ha='center', va='center', fontsize=14, backgroundcolor=pink)\n";
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

inline void Unstructured::Delaunay (Array<double> const & X, Array<double> const & Y, int Tag)
{
    // check
    if (NDim==3)            throw new Fatal("Unstructured::Delaunay: This method is only available for 2D");
    if (X.Size()!=Y.Size()) throw new Fatal("Unstructured::Delaunay: Size of X and Y arrays must be equal (%d!=%d)",X.Size(),Y.Size());

    // points
    size_t n = X.Size();
    std::vector<double> pts (n*3, 0.0);
    for (size_t i=0; i<n; ++i) { pts[i*3]=X[i]; pts[i*3+1]=Y[i]; }

    // convex hull
    std::vector<int> hull;
    ConvexHull2D (pts, hull);
    if (hull.size()<3) throw new Fatal("Unstructured::Delaunay: Could not build a convex hull");

    // set up the geometry and generate
    Set (n, hull.size(), 1, 0);
    std::vector<bool> on_hull (n, false);
    for (size_t i=0; i<hull.size(); ++i) on_hull[hull[i]] = true;
    for (size_t i=0; i<n; ++i)
    {
        SetPnt (i, 0, X[i], Y[i], 0.0);
        if (!on_hull[i]) _embed_pts.push_back (static_cast<int>(i));
    }
    for (size_t i=0; i<hull.size(); ++i) SetSeg (i, 0, hull[i], hull[(i+1)%hull.size()]);
    double cx = 0.0, cy = 0.0;
    for (size_t i=0; i<hull.size(); ++i) { cx += pts[hull[i]*3]; cy += pts[hull[i]*3+1]; }
    cx /= hull.size(); cy /= hull.size();
    SetReg (0, Tag, -1.0, cx, cy, 0.0);
    Generate ();
}

inline void Unstructured::Delaunay (Array<double> const & X, Array<double> const & Y, Array<double> const & Z, int Tag)
{
    // check
    if (NDim==2)            throw new Fatal("Unstructured::Delaunay: This method is only available for 3D");
    if (X.Size()!=Y.Size()) throw new Fatal("Unstructured::Delaunay: Size of X and Y arrays must be equal (%d!=%d)",X.Size(),Y.Size());
    if (Z.Size()!=Y.Size()) throw new Fatal("Unstructured::Delaunay: Size of Z and Y arrays must be equal (%d!=%d)",Z.Size(),Y.Size());

    // points
    size_t n = X.Size();
    std::vector<double> pts (n*3, 0.0);
    for (size_t i=0; i<n; ++i) { pts[i*3]=X[i]; pts[i*3+1]=Y[i]; pts[i*3+2]=Z[i]; }

    // convex hull
    std::vector<std::array<int,3> > faces;
    ConvexHull3D (pts, faces);
    if (faces.size()<4) throw new Fatal("Unstructured::Delaunay: Could not build a convex hull");

    // set up the geometry and generate
    Set (n, faces.size(), 1, 0);
    std::vector<bool> on_hull (n, false);
    for (size_t i=0; i<faces.size(); ++i)
        on_hull[faces[i][0]] = on_hull[faces[i][1]] = on_hull[faces[i][2]] = true;
    for (size_t i=0; i<n; ++i)
    {
        SetPnt (i, 0, X[i], Y[i], Z[i]);
        if (!on_hull[i]) _embed_pts.push_back (static_cast<int>(i));
    }
    for (size_t i=0; i<faces.size(); ++i)
    {
        SetFac (i, 0, Array<int>(faces[i][0], faces[i][1], faces[i][2]));
    }
    double cx=0.0, cy=0.0, cz=0.0;
    for (size_t i=0; i<n; ++i) { cx += X[i]; cy += Y[i]; cz += Z[i]; }
    cx /= n; cy /= n; cz /= n;
    SetReg (0, Tag, -1.0, cx, cy, cz);
    Generate ();
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
