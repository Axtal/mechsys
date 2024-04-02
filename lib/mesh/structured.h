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

#ifndef MECHSYS_MESH_STRUCTURED_H
#define MECHSYS_MESH_STRUCTURED_H

/* LOCAL indexes of VERTICES, EDGES, and FACES

  2D:
                Vertices                              Edges
   y
   |        3      6      2                             2
   +--x      @-----@-----@                        +-----------+
             |           |                        |           |
             |           |                        |           |
           7 @           @ 5                     3|           |1
             |           |                        |           |
             |           |                        |           |
             @-----@-----@                        +-----------+
            0      4      1                             0

  3D:
                  Vertices                            Faces
    z
    |           4        15        7
   ,+--y         @-------@--------@                 +----------------+
 x'            ,'|              ,'|               ,'|              ,'|
          12 @'  |         14 ,'  |             ,'  |  ___       ,'  |
           ,'    |16        ,@    |19         ,'    |,'5,'  [0],'    |
     5   ,'      @      6 ,'      @         ,'      |~~~     ,'      |
       @'=======@=======@'        |       +'===============+'  ,'|   |
       |      13 |      |         |       |   ,'|   |      |   |3|   |
       |         |      |  11     |       |   |2|   |      |   |,'   |
    17 |       0 @- - - | @- - - -@       |   |,'   +- - - | +- - - -+
       @       ,'       @       ,' 3      |       ,'       |       ,'
       |   8 @'      18 |     ,'          |     ,' [1]  ___|     ,'
       |   ,'           |   ,@ 10         |   ,'      ,'4,'|   ,'
       | ,'             | ,'              | ,'        ~~~  | ,'
       @-------@--------@'                +----------------+'
     1         9         2
*/

// Std Lib
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>   // for DBL_EPSILON
#include <cstdarg>  // for va_list, va_start, va_end
#include <map>      // for std::map

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/maps.h>

namespace Mesh
{

    /* TODO:
     *        1) Add additional checks, especially for 3D meshes
     *        2) Check second-order (o2) elements
     *        3) Add quality check and improvement
     *        4) Create list of neighbours blocks/edges
     */


/////////////////////////////////////////////////////////////////////////////////////////// Block /////


typedef std::map<int,SDPair> Side2Ctr_t; ///< Side to constraint map

class Block
{
public:
    // Constructor
    Block () : Nx(1), Ny(1), Nz(1) {}

    /** Ex:
     *      Set (2, -1, 4,                 // 2D, Tag, 4 vertices
     *           -1.0,  0.0, 0.0,          // vtag, x, y, [z,]
     *           -2.0,  1.0, 0.0,          // vtag, x, y, [z,]
     *           -3.0,  1.0, 1.0,          // vtag, x, y, [z,]
     *           -4.0,  0.0, 1.0,          // vtag, x, y, [z,]
     *           -10.0,-20.0,-30.0,-40.0); // boundary (edge/face) tags
     *  Note:
     *     After NVerts, all data must be (double)   */
    void Set (int NDim, int Tag, size_t NVerts, ...); ///< Set block

    // Methods
    void   SetNx       (size_t Nx, double Ax=0.0, bool NonLin=false);                ///< Set number of x divisions with all x weights equal to 1.0
    void   SetNy       (size_t Ny, double Ay=0.0, bool NonLin=false);                ///< Set number of y divisions with all y weights equal to 1.0
    void   SetNz       (size_t Nz, double Az=0.0, bool NonLin=false);                ///< Set number of z divisions with all z weights equal to 1.0
    void   SetNx       (Array<double> const & TheWx);                                ///< Set divisions along x given weights
    double Diagonal    () const;                                                     ///< Find diagonal of block
    void   GenMidNodes ();                                                           ///< Generate mid nodes of block
    bool   BryIDs      (size_t i, size_t j, size_t k, Array<int> & BryIDs) const;    ///< ID of boundary sides (edges or faces) of the block where i,j,k is located on
    int    GetVTag     (size_t i, size_t j, size_t k) const;                         ///< Get vertex tag
    void   AddArcCtr   (int Side, double x, double y, double R, double Tol=1.0e-10); ///< Add Arc constraint to nodes along side # Side
    void   AddCylXCtr  (int Side, double y, double z, double R, double Tol=1.0e-10); ///< Constarint: cylinder parallel to x
    void   ApplyCtr    (Array<int> const & Sides, Vertex & V) const;                 ///< Apply constraints to Vertex V which is belong to sides Sides

    // Data
    int           NDim;              ///< Space dimension
    int           Tag;               ///< Tag
    Mat_t         C;                 ///< Coordinates
    Array<int>    VTags;             ///< Vertex tags
    Array<int>    BryTags;           ///< Boundary (edge or Face) tags
    size_t        Nx,Ny,Nz;          ///< Number of divisons
    Array<double> Wx,Wy,Wz;          ///< Weights
    double        SumWx,SumWy,SumWz; ///< Sum of weights
    Array<bool>   SideHasCtr;        ///< Side has constraint ?
    Side2Ctr_t    Side2Ctr;          ///< Side to constraint map

    /* 2D:      _             _
     *     C = |  x0 x1 x2 x3  |
     *         |_ y0 y1 y2 y3 _|
     *     or
     *          _                          _
     *     C = |  x0 x1 x2 x3  x4 x5 x6 x7  |
     *         |_ y0 y1 y2 y3  y4 y5 y6 y7 _|
     *
     * 3D:      _                          _
     *         |  x0 x1 x2 x3 x4 x5 x6 x7   |
     *     C = |  y0 y1 y2 y3 y4 y5 y6 y7   |
     *         |_ z0 z1 z2 z3 z4 z5 z6 z7  _|
     *
     *     or   _                                         _
     *         |  x0 x1 x2 x3 x4 x5 x6 x7 ... x17 x18 x19  |
     *     C = |  y0 y1 y2 y3 y4 y5 y6 y7 ... y17 y18 y19  |
     *         |_ z0 z1 z2 z3 z4 z5 z6 z7 ... z17 z18 z19 _|
     */

#ifdef USE_BOOST_PYTHON
    void PySet (BPy::dict const & Dat);
#endif

private:
    void _initialize (size_t NVerts);
};


/////////////////////////////////////////////////////////////////////////////////////////// Structured /////


class Structured : public virtual Mesh::Generic
{
public:
    // Structures
    struct VertexInfo
    {
        bool         OnBry;     ///< On boundary ?
        int          BlkNum;    ///< Number of the block which created this vertex
        bool         Dupl;      ///< Duplicated vertex ?
        Array<int>   BlkBryIDs; ///< IDs of boudaries (edges/faces) where this vertex is sitting on
    };

    // Constructor
    Structured (int NDim) : Mesh::Generic(NDim) {}

    // Destructor
    ~Structured () { Erase(); }

    // Methods
    void Generate      (Array<Block> const & Blks, bool O2=false);                                                ///< Boundary marks are set first for Faces, then Edges, then Vertices (if any)
    void GenBox        (bool O2=true, int Nx=2, int Ny=2, int Nz=2, double Lx=1.0, double Ly=1.0, double Lz=1.0); ///< Generate a cube with dimensions Lx,Ly,Lz and with tags on faces
    void GenQRing      (bool O2=true, int Nx=2, int Ny=2, double r=1.0, double R=2.0,
                        size_t Nb=2, double Ax=1.0, bool NonLin=false, char const * WeightsX=NULL); ///< Generate a quarter of a ring
    void GenQRing      (bool O2, int Nx, int Ny, int Nz, double r, double R, double t);             ///< Generate a quarter of a 3D ring (t=thickness)
    void GenQDisk      (bool O2, int Nx1, int Nx2, int Ny, double r, double R);                     ///< Generate a quarter of a disk
    void GenQPlateHole (double r, double R, double L, int Nx1, int Nx2, int Ny1, int Ny2, bool O2); ///< Generate a quarter of a square plate with inner hole
    void ShapeFunc     (double r, double s, double t);                                              ///< Calculate shape function. Return results on N

    // Data
    Vec_t N; ///< Shape functions

#ifdef USE_BOOST_PYTHON
    void PyGenerate (BPy::list const & Blks, bool O2=false);
#endif
};


//////////////////////////////////////////////////////////////////////////////////// Block: Implementation /////


inline void Block::Set (int TheNDim, int TheTag, size_t NVerts, ...)
{
    // allocate data
    NDim = TheNDim;
    Tag  = TheTag;
    _initialize (NVerts);

    // read data
    va_list   arg_list;
    va_start (arg_list, NVerts);
    for (size_t i=0; i<NVerts; ++i)
    {
        VTags[i] = static_cast<int>(va_arg(arg_list,double));
        C(0,i)   = va_arg(arg_list,double);
        C(1,i)   = va_arg(arg_list,double);  if (NDim==3)
        C(2,i)   = va_arg(arg_list,double);
    }
    for (size_t i=0; i<BryTags.Size(); ++i) BryTags[i] = static_cast<int>(va_arg(arg_list,double));
    va_end (arg_list);

    // check connectivity (2D)
    if (NDim==2)
    {
        Vec3_t p0, p1, p2;
        for (size_t i=1; i<NVerts-1; ++i)
        {
            p0 = C(0,i-1)-C(0,i),  C(1,i-1)-C(1,i),  0.0;
            p1 = C(0,i+1)-C(0,i),  C(1,i+1)-C(1,i),  0.0;
            p2 = blitz::cross (p1,p0);
            if (p2(2)<0.0) throw new Fatal("Mesh::Block::Set: Order of vertices is incorrect (it must be counter-clockwise)");
        }
    }

    // generate mid nodes
    if (NDim==2 && NVerts==4) GenMidNodes();
    if (NDim==3 && NVerts==8) GenMidNodes();

    // constraints
    if (NDim==2)
    {
        SideHasCtr.Resize    (4);
        SideHasCtr.SetValues (false);
    }
    else
    {
        SideHasCtr.Resize    (6);
        SideHasCtr.SetValues (false);
    }
}

inline void Block::SetNx (size_t TheNx, double Ax, bool NonLin)
{
    if (TheNx<1) throw new Fatal("Block::SetNx: Number of divisions of block for structured mesh must be greater than or equal to 1. Nx=%zd is invalid",TheNx);
    Nx = TheNx;
    Wx.Resize(Nx);
    if (NonLin) for (size_t i=0; i<Nx; ++i) Wx[i] = pow(i+1.0,Ax);
    else        for (size_t i=0; i<Nx; ++i) Wx[i] = 1.0+Ax*i;

    SumWx = 0.0;
    for (size_t i=0; i<Nx; ++i) SumWx += Wx[i];
    Wx.Push(0.0); // extra value just to help loop over weights
}

inline void Block::SetNx (Array<double> const & TheWx)
{
    Nx = TheWx.Size();
    Wx = TheWx;

    SumWx = 0.0;
    for (size_t i=0; i<Nx; ++i) SumWx += Wx[i];
    Wx.Push(0.0); // extra value just to help loop over weights
}

inline void Block::SetNy (size_t TheNy, double Ay, bool NonLin)
{
    if (TheNy<1) throw new Fatal("Block::SetNy: Number of divisions of block for structured mesh must be greater than or equal to 1. Ny=%zd is invalid",TheNy);
    Ny = TheNy;
    Wy.Resize(Ny);
    if (NonLin) for (size_t i=0; i<Ny; ++i) Wy[i] = pow(i+1.0,Ay);
    else        for (size_t i=0; i<Ny; ++i) Wy[i] = 1.0+Ay*i;

    SumWy = 0.0;
    for (size_t i=0; i<Ny; ++i) SumWy += Wy[i];
    Wy.Push(0.0); // extra value just to help loop over weights
}

inline void Block::SetNz (size_t TheNz, double Az, bool NonLin)
{
    if (TheNz<1) throw new Fatal("Block::SetNz: Number of divisions of block for structured mesh must be greater than or equal to 1. Nz=%zd is invalid",TheNz);
    Nz = TheNz;
    Wz.Resize(Nz);
    if (NonLin) for (size_t i=0; i<Nz; ++i) Wz[i] = pow(i+1.0,Az);
    else        for (size_t i=0; i<Nz; ++i) Wz[i] = 1.0+Az*i;

    SumWz = 0.0;
    for (size_t i=0; i<Nz; ++i) SumWz += Wz[i];
    Wz.Push(0.0); // extra value just to help loop over weights
}

inline double Block::Diagonal () const
{
    if (NDim==3) return sqrt(pow(C(0,6)-C(0,0),2.0)+pow(C(1,6)-C(1,0),2.0)+pow(C(2,6)-C(2,0),2.0));
    else         return sqrt(pow(C(0,2)-C(0,0),2.0)+pow(C(1,2)-C(1,0),2.0));
}

inline void Block::GenMidNodes ()
{
    if (NDim==3)
    {
        for (size_t i=0; i<3; ++i)
        {
            C(i, 8) = (C(i,0) + C(i,1))/2.0;
            C(i, 9) = (C(i,1) + C(i,2))/2.0;
            C(i,10) = (C(i,2) + C(i,3))/2.0;
            C(i,11) = (C(i,3) + C(i,0))/2.0;

            C(i,12) = (C(i,4) + C(i,5))/2.0;
            C(i,13) = (C(i,5) + C(i,6))/2.0;
            C(i,14) = (C(i,6) + C(i,7))/2.0;
            C(i,15) = (C(i,7) + C(i,4))/2.0;

            C(i,16) = (C(i,0) + C(i,4))/2.0;
            C(i,17) = (C(i,1) + C(i,5))/2.0;
            C(i,18) = (C(i,2) + C(i,6))/2.0;
            C(i,19) = (C(i,3) + C(i,7))/2.0;
        }
    }
    else
    {
        for (size_t i=0; i<2; ++i)
        {
            C(i,4) = (C(i,0) + C(i,1))/2.0;
            C(i,5) = (C(i,1) + C(i,2))/2.0;
            C(i,6) = (C(i,2) + C(i,3))/2.0;
            C(i,7) = (C(i,3) + C(i,0))/2.0;
        }
    }
}

inline bool Block::BryIDs (size_t i, size_t j, size_t k, Array<int> & BryIDs) const
{
    bool onbry = false;
    if (NDim==2)
    {
        if (j==0 ) { if (BryIDs.Find(0)<0) BryIDs.Push(0); onbry = true; } // bottom
        if (i==Nx) { if (BryIDs.Find(1)<0) BryIDs.Push(1); onbry = true; } // right
        if (j==Ny) { if (BryIDs.Find(2)<0) BryIDs.Push(2); onbry = true; } // top
        if (i==0 ) { if (BryIDs.Find(3)<0) BryIDs.Push(3); onbry = true; } // left
    }
    else
    {
        if (i==0 ) { if (BryIDs.Find(0)<0) BryIDs.Push(0); onbry = true; } // behind
        if (i==Nx) { if (BryIDs.Find(1)<0) BryIDs.Push(1); onbry = true; } // front
        if (j==0 ) { if (BryIDs.Find(2)<0) BryIDs.Push(2); onbry = true; } // left
        if (j==Ny) { if (BryIDs.Find(3)<0) BryIDs.Push(3); onbry = true; } // right
        if (k==0 ) { if (BryIDs.Find(4)<0) BryIDs.Push(4); onbry = true; } // bottom
        if (k==Nz) { if (BryIDs.Find(5)<0) BryIDs.Push(5); onbry = true; } // top
    }
    return onbry;
}

inline int Block::GetVTag (size_t i, size_t j, size_t k) const
{
    if (NDim==2)
    {
        if (i==0  && j==0 ) return VTags[0];
        if (i==Nx && j==0 ) return VTags[1];
        if (i==Nx && j==Ny) return VTags[2];
        if (i==0  && j==Ny) return VTags[3];
    }
    else
    {
        if (i==0  && j==0  && k==0 ) return VTags[0];
        if (i==Nx && j==0  && k==0 ) return VTags[1];
        if (i==Nx && j==Ny && k==0 ) return VTags[2];
        if (i==0  && j==Ny && k==0 ) return VTags[3];
        if (i==0  && j==0  && k==Nz) return VTags[4];
        if (i==Nx && j==0  && k==Nz) return VTags[5];
        if (i==Nx && j==Ny && k==Nz) return VTags[6];
        if (i==0  && j==Ny && k==Nz) return VTags[7];
    }
    return 0;
}

inline void Block::AddArcCtr (int Side, double xc, double yc, double R, double Tol)
{
    SDPair data;
    data.Set ("type",1.0); // 1.0 = Arc
    data.Set ("xc", xc);
    data.Set ("yc", yc);
    data.Set ("R",  R);
    data.Set ("Tol",Tol);
    Side2Ctr[Side]   = data;
    SideHasCtr[Side] = true;
}

inline void Block::AddCylXCtr (int Side, double yc, double zc, double R, double Tol)
{
    SDPair data;
    data.Set ("type",2.0); // 2.0 = Cylinder parallel to x
    data.Set ("xc", yc);
    data.Set ("yc", zc);
    data.Set ("R",  R);
    data.Set ("Tol",Tol);
    Side2Ctr[Side]   = data;
    SideHasCtr[Side] = true;
}

inline void Block::ApplyCtr (Array<int> const & Sides, Vertex & V) const
{
    for (size_t m=0; m<Sides.Size(); ++m)
    {
        int side = Sides[m];
        if (SideHasCtr[side])
        {
            Side2Ctr_t::const_iterator it = Side2Ctr.find(side);
            double *X, *Y;
            if (static_cast<int>(it->second("type"))==1) // arc
            {
                X = &V.C(0);
                Y = &V.C(1);
            }
            else if (static_cast<int>(it->second("type"))==2) // cylinder parallel to x
            {
                X = &V.C(1);
                Y = &V.C(2);
            }
            else throw new Fatal("Mesh::Block: __internal_error__: Unknown constraint type (%d)",it->second("type"));
            double x   = (*X);
            double y   = (*Y);
            double xc  = it->second("xc");
            double yc  = it->second("yc");
            double R   = it->second("R");
            double tol = it->second("Tol");
            double F   = pow(x-xc,2.0) + pow(y-yc,2.0) - R*R;
            if (fabs(F)>tol)
            {
                //std::cout << "F_{before} = " << F;
                double nx = 2.0*(x-xc);
                double ny = 2.0*(y-yc);
                double A  = nx*nx + ny*ny;
                double B  = 2.0*nx*(x-xc) + 2.0*ny*(y-yc);
                double D  = B*B - 4.0*A*F;
                if (D<0.0) throw new Fatal("Mesh::Block::ApplyCtr: Constraint failed with D=%g",D);
                double a = (F<0.0 ? (-B+sqrt(D))/(2.0*A) : (-B-sqrt(D))/(2.0*A));
                (*X) = x + a*nx;
                (*Y) = y + a*ny;
                F = pow((*X)-xc,2.0) + pow((*Y)-yc,2.0) - R*R;
                if (fabs(F)>tol) throw new Fatal("Mesh::Block::ApplyCtr: Constraint failed with F=%g",F);
                //std::cout << "   F_{after} = " << F << std::endl;
            }
        }
    }
}

inline void Block::_initialize (size_t NVerts)
{
    if (NDim==2)
    {
        if (NVerts==4 || NVerts==8)
        {
            C.change_dim      (NDim,8);
            VTags.Resize      (8); // 8 vertices
            BryTags.Resize    (4); // 4 edges
            VTags.SetValues   (0);
            BryTags.SetValues (0);
        }
        else throw new Fatal("Block::_initialize: With NDim=2 Number of vertices must be either 4 or 8 (%d is invalid)",NVerts);
    }
    else if (NDim==3)
    {
        if (NVerts==8 || NVerts==20)
        {
            C.change_dim      (NDim,20);
            VTags.Resize      (20); // 20 vertices
            BryTags.Resize    (6);  // 6 faces
            VTags.SetValues   (0);
            BryTags.SetValues (0);
        }
        else throw new Fatal("Block::_initialize: With NDim=3 Number of vertices must be either 8 or 20 (%d is invalid)",NVerts);
    }
    else throw new Fatal("Block::_initialize: NDim=%d is invalid",NDim);
    SetNx (Nx);
    SetNy (Ny);
    SetNz (Nz);
}

std::ostream & operator<< (std::ostream & os, Block const & B)
{
    os << "B={'ndim':" << B.NDim << ", 'tag':" << B.Tag;
    os << ", 'nx':"   << B.Nx;
    os << ", 'ny':"   << B.Ny;  if (B.NDim==3)
    os << ", 'nz':"   << B.Nz;
    os << ", 'brytags':[";
    for (size_t i=0; i<B.BryTags.Size(); ++i)
    {
        os << B.BryTags[i];
        if (i!=B.BryTags.Size()-1) os << ",";
    }
    os << "],\n   'V':[";
    for (size_t i=0; i<B.C.num_cols(); ++i)
    {
        os << "[" << Util::_4  << B.VTags[i];
        os << "," << Util::_8s << B.C(0,i);
        os << "," << Util::_8s << B.C(1,i);  if (B.NDim==3)
        os << "," << Util::_8s << B.C(2,i);
        if (i==B.C.num_cols()-1) os << "]]}";
        else                     os << "],\n        ";
    }
    return os;
}

#ifdef USE_BOOST_PYTHON

inline void Block::PySet (BPy::dict const & Dat)
{
    // allocate data
    NDim = BPy::extract<int>(Dat["ndim"])();
    Tag  = BPy::extract<int>(Dat["tag"])();
    Nx   = BPy::extract<int>(Dat["nx"])();
    Ny   = BPy::extract<int>(Dat["ny"])();
    Nz   = (Dat.has_key("nz") ? BPy::extract<int>(Dat["nz"])() : 0);
    BPy::list const & verts = BPy::extract<BPy::list>(Dat["V"])();
    size_t nverts = BPy::len(verts);
    _initialize (nverts);

    // read data
    for (size_t i=0; i<nverts; ++i)
    {
        BPy::list const & line = BPy::extract<BPy::list>(verts[i])();
        VTags[i] = BPy::extract<int   >(line[0])();
        C(0,i)   = BPy::extract<double>(line[1])();
        C(1,i)   = BPy::extract<double>(line[2])();  if (NDim==3)
        C(2,i)   = BPy::extract<double>(line[3])();
    }
    if (Dat.has_key("brytags"))
    {
        BPy::list const & brytags = BPy::extract<BPy::list>(Dat["brytags"])();
        for (size_t i=0; i<BryTags.Size(); ++i) BryTags[i] = BPy::extract<int>(brytags[i])();
    }

    // check connectivity (2D)
    if (NDim==2)
    {
        Vec3_t p0, p1, p2;
        for (size_t i=1; i<nverts-1; ++i)
        {
            p0 = C(0,i-1)-C(0,i),  C(1,i-1)-C(1,i),  0.0;
            p1 = C(0,i+1)-C(0,i),  C(1,i+1)-C(1,i),  0.0;
            p2 = blitz::cross (p1,p0);
            if (p2(2)<0.0) throw new Fatal("Mesh::Block::PySet: Order of vertices is incorrect (it must be counter-clockwise)");
        }
    }

    // generate mid nodes
    if (NDim==2 && nverts==4) GenMidNodes();
    if (NDim==3 && nverts==8) GenMidNodes();
}

#endif


//////////////////////////////////////////////////////////////////////////////////// Structured: Implementation /////////


inline void Structured::Generate (Array<Block> const & Blks, bool O2)
{
    // data
    Array<Vertex*> verts;     // Vertices (with duplicates)
    Array<Vertex*> verts_bry; // Vertices on boundary (with duplicates)
    Array<Vertex*> verts_m1;  // X O2 Vertices (with duplicates)
    Array<Vertex*> verts_m2;  // Y O2 Vertices (with duplicates)
    Array<Vertex*> verts_m3;  // Z O2 Vertices (with duplicates)

    // info
    Util::Stopwatch stopwatch(/*activated*/WithInfo);

    // check
    if (Blks.Size()<1) throw new Fatal("Structured::Generate: Number of blocks must be greater than 0 (%d is invalid)",Blks.Size());

    // erase previous mesh
    Erase();

    // check if the first block is 3D
    if (Blks[0].NDim!=NDim) throw new Fatal("Structured::Generate: NDim=%d of blocks must be equal to NDim=%d of Mesh::Structured",Blks[0].NDim,NDim);
    N.change_dim((NDim==3 ? 20 : 8)); // resize the shape values vector

    // vertices information
    std::map<Vertex*,VertexInfo> vinfo; // map: pointer2vertex => info

    // generate vertices and elements (with duplicates)
    double min_diagonal = Blks[0].Diagonal(); // minimum diagonal among all elements
    for (size_t b=0; b<Blks.Size(); ++b)
    {
        // check if all blocks have the same space dimension
        if (Blks[b].NDim!=NDim) throw new Fatal("Structured::Generate: All blocks must have the same space dimension");

        // check Nx Ny Nz
        if (Blks[b].Nx<1) throw new Fatal("Structured::Generate: All blocks must have Nx>0");
        if (Blks[b].Ny<1) throw new Fatal("Structured::Generate: All blocks must have Ny>0");  if (NDim==3)
        if (Blks[b].Nz<1) throw new Fatal("Structured::Generate: All blocks must have Nz>0");

        // generate
        double t_prev = -2.0;
        double t      = -1.0; // initial Z natural coordinate
        for (size_t k=0; k<(NDim==3 ? Blks[b].Nz+1 : 1); ++k)
        {
            double s_prev = -2.0;
            double s      = -1.0; // initial Y natural coordinate
            for (size_t j=0; j<Blks[b].Ny+1; ++j)
            {
                double r_prev = -2.0;
                double r      = -1.0; // initial X natural coordinate
                for (size_t i=0; i<Blks[b].Nx+1; ++i)
                {
                    // new vertex
                    ShapeFunc (r,s,t);
                    Vec_t c(Blks[b].C * N);
                    Vertex * v      = new Vertex;
                    v->Tag          = 0;
                    v->C            = c(0), c(1), (NDim==3?c(2):0.0);
                    vinfo[v].BlkNum = b;
                    vinfo[v].Dupl   = false;
                    vinfo[v].OnBry  = Blks[b].BryIDs (i,j,k,vinfo[v].BlkBryIDs);
                    if (vinfo[v].OnBry)
                    {
                        v->Tag = Blks[b].GetVTag (i,j,k);
                        verts_bry.Push (v);
                        Blks[b].ApplyCtr (vinfo[v].BlkBryIDs, (*v)); // constraints
                    }
                    verts.Push (v);

                    // new O2 vertices
                    Vertex * v1 = NULL;
                    Vertex * v2 = NULL;
                    Vertex * v3 = NULL;
                    if (O2)
                    {
                        if (i!=0)
                        {
                            ShapeFunc ((r_prev+r)/2.0,s,t);
                            c                = Blks[b].C * N;
                            v1               = new Vertex;
                            v1->Tag          = 0;
                            v1->C            = c(0), c(1), (NDim==3?c(2):0.0);
                            vinfo[v1].BlkNum = b;
                            vinfo[v1].Dupl   = false;
                            size_t m = (Blks[b].Nx==1 ? 2 : i-1);
                            if (i==Blks[b].Nx) vinfo[v1].OnBry = Blks[b].BryIDs(m,j,k, vinfo[v1].BlkBryIDs);
                            else               vinfo[v1].OnBry = Blks[b].BryIDs(i,j,k, vinfo[v1].BlkBryIDs);
                            if (vinfo[v1].OnBry)
                            {
                                verts_bry.Push (v1);
                                Blks[b].ApplyCtr (vinfo[v1].BlkBryIDs, (*v1)); // constraints
                            }
                            verts_m1.Push (v1);
                        }
                        if (j!=0)
                        {
                            ShapeFunc (r,(s_prev+s)/2.0,t);
                            c                = Blks[b].C * N;
                            v2               = new Vertex;
                            v2->Tag          = 0;
                            v2->C            = c(0), c(1), (NDim==3?c(2):0.0);
                            vinfo[v2].BlkNum = b;
                            vinfo[v2].Dupl   = false;
                            size_t m = (Blks[b].Ny==1 ? 2 : j-1);
                            if (j==Blks[b].Ny) vinfo[v2].OnBry = Blks[b].BryIDs(i,m,k, vinfo[v2].BlkBryIDs);
                            else               vinfo[v2].OnBry = Blks[b].BryIDs(i,j,k, vinfo[v2].BlkBryIDs);
                            if (vinfo[v2].OnBry)
                            {
                                verts_bry.Push (v2);
                                Blks[b].ApplyCtr (vinfo[v2].BlkBryIDs, (*v2)); // constraints
                            }
                            verts_m2.Push (v2);
                        }
                        if (k!=0)
                        {
                            ShapeFunc (r,s,(t_prev+t)/2.0);
                            c                = Blks[b].C * N;
                            v3               = new Vertex;
                            v3->Tag          = 0;
                            v3->C            = c(0), c(1), (NDim==3?c(2):0.0);
                            vinfo[v3].BlkNum = b;
                            vinfo[v3].Dupl   = false;
                            size_t m = (Blks[b].Nz==1 ? 2 : k-1);
                            if (k==Blks[b].Nz) vinfo[v3].OnBry = Blks[b].BryIDs(i,j,m, vinfo[v3].BlkBryIDs);
                            else               vinfo[v3].OnBry = Blks[b].BryIDs(i,j,k, vinfo[v3].BlkBryIDs);
                            if (vinfo[v3].OnBry)
                            {
                                verts_bry.Push (v3);
                                Blks[b].ApplyCtr (vinfo[v3].BlkBryIDs, (*v3)); // constraints
                            }
                            verts_m3.Push (v3);
                        }
                    }

                    // new element
                    if (i!=0 && j!=0 && (NDim==3 ? k!=0 : true))
                    {
                        Cell * e    = new Cell;
                        e->ID       = Cells.Size();
                        e->Tag      = Blks[b].Tag;
                        e->PartID   = 0; // partition id (domain decomposition)
                        int nvonbry = 0; // number of vertices of this element on boundary
                        size_t  pnv = verts.Size()-1; // previous number of vertices
                        Cells.Push (e);
                        if (NDim==2)
                        {
                            // connectivity
                            e->V.Resize((O2?8:4));
                            e->V[0] = verts[pnv - 1 - (Blks[b].Nx+1)];
                            e->V[1] = verts[pnv     - (Blks[b].Nx+1)];
                            e->V[2] = verts[pnv];
                            e->V[3] = verts[pnv - 1];
                            if (O2)
                            {
                                size_t pnv1 = verts_m1.Size()-1;
                                size_t pnv2 = verts_m2.Size()-1;
                                e->V[4] = verts_m1[pnv1 - Blks[b].Nx];
                                e->V[5] = verts_m2[pnv2];
                                e->V[6] = verts_m1[pnv1];
                                e->V[7] = verts_m2[pnv2 - 1];
                            }
                            // shares and onbry information
                            for (size_t m=0; m<e->V.Size(); ++m)
                            {
                                Share s = {e,m};
                                if (vinfo[e->V[m]].OnBry) nvonbry++;
                                e->V[m]->Shares.Push (s);
                            }
                            // diagonal
                            double d = sqrt(pow(e->V[2]->C(0) - e->V[0]->C(0),2.0)+
                                            pow(e->V[2]->C(1) - e->V[0]->C(1),2.0));
                            if (d<min_diagonal) min_diagonal = d;
                        }
                        else
                        {
                            // connectivity
                            e->V.Resize((O2?20:8));
                            e->V[0] = verts[pnv - 1 - (Blks[b].Nx+1) - (Blks[b].Nx+1)*(Blks[b].Ny+1)];
                            e->V[1] = verts[pnv     - (Blks[b].Nx+1) - (Blks[b].Nx+1)*(Blks[b].Ny+1)];
                            e->V[2] = verts[pnv                      - (Blks[b].Nx+1)*(Blks[b].Ny+1)];
                            e->V[3] = verts[pnv - 1                  - (Blks[b].Nx+1)*(Blks[b].Ny+1)];
                            e->V[4] = verts[pnv - 1 - (Blks[b].Nx+1)];
                            e->V[5] = verts[pnv -     (Blks[b].Nx+1)];
                            e->V[6] = verts[pnv];
                            e->V[7] = verts[pnv - 1];
                            if (O2)
                            {
                                size_t pnv1 = verts_m1.Size()-1;
                                size_t pnv2 = verts_m2.Size()-1;
                                size_t pnv3 = verts_m3.Size()-1;
                                e->V[ 8] = verts_m1[pnv1 - Blks[b].Nx - Blks[b].Nx*(Blks[b].Ny+1)];
                                e->V[ 9] = verts_m2[pnv2              - Blks[b].Ny*(Blks[b].Nx+1)];
                                e->V[10] = verts_m1[pnv1              - Blks[b].Nx*(Blks[b].Ny+1)];
                                e->V[11] = verts_m2[pnv2 -1           - Blks[b].Ny*(Blks[b].Nx+1)];

                                e->V[12] = verts_m1[pnv1 - Blks[b].Nx];
                                e->V[13] = verts_m2[pnv2];
                                e->V[14] = verts_m1[pnv1];
                                e->V[15] = verts_m2[pnv2 - 1];

                                e->V[16] = verts_m3[pnv3 - 1 - (Blks[b].Nx+1)];
                                e->V[17] = verts_m3[pnv3 -     (Blks[b].Nx+1)];
                                e->V[18] = verts_m3[pnv3];
                                e->V[19] = verts_m3[pnv3 - 1];
                            }
                            // shares and onbry information
                            for (size_t m=0; m<e->V.Size(); ++m)
                            {
                                Share s = {e,m};
                                if (vinfo[e->V[m]].OnBry) nvonbry++;
                                e->V[m]->Shares.Push (s);
                            }
                            // diagonal
                            double d = sqrt(pow(e->V[6]->C(0) - e->V[0]->C(0),2.0)+
                                            pow(e->V[6]->C(1) - e->V[0]->C(1),2.0)+
                                            pow(e->V[6]->C(2) - e->V[0]->C(2),2.0));
                            if (d<min_diagonal) min_diagonal = d;
                        }
                        // apply boundary tags
                        if (nvonbry>0)
                        {
                            bool has_bry_tag = false;
                            for (size_t m=0; m<(NDim==2?4:8); ++m)
                            {
                                for (size_t n=0; n<vinfo[e->V[m]].BlkBryIDs.Size(); ++n)
                                {
                                    int side = vinfo[e->V[m]].BlkBryIDs[n];
                                    int tag  = Blks[b].BryTags[side];
                                    if (tag<0)
                                    {
                                        e->BryTags[side] = tag;
                                        has_bry_tag      = true;
                                    }
                                }
                            }
                            if (has_bry_tag) TgdCells.Push (e);
                        }
                    }
                    // next r
                    r_prev = r;
                    r += (2.0/Blks[b].SumWx) * Blks[b].Wx[i];
                }
                // next s
                s_prev = s;
                s += (2.0/Blks[b].SumWy) * Blks[b].Wy[j];
            }
            // next t
            t_prev = t;
            t += (NDim==3 ? (2.0/Blks[b].SumWz) * Blks[b].Wz[k] : 0.0);
        }
    }

    // remove duplicates
    long ncomp = 0;                  // number of comparisons
    long ndupl = 0;                  // number of duplicates
    double tol = min_diagonal*0.001; // tolerance to decide whether two vertices are coincident or not
    if (Blks.Size()>1)
    {
        for (size_t i=0; i<verts_bry.Size(); ++i)
        {
            for (size_t j=i+1; j<verts_bry.Size(); ++j)
            {
                if (vinfo[verts_bry[i]].BlkNum!=vinfo[verts_bry[j]].BlkNum) // vertices are located on different blocks
                {
                    // check distance
                    double dist = sqrt(          pow(verts_bry[i]->C(0)-verts_bry[j]->C(0),2.0)+
                                                 pow(verts_bry[i]->C(1)-verts_bry[j]->C(1),2.0)+
                                      (NDim==3 ? pow(verts_bry[i]->C(2)-verts_bry[j]->C(2),2.0) : 0.0));
                    if (dist<tol)
                    {
                        // mark duplicated
                        if (vinfo[verts_bry[j]].Dupl==false) // vertex not tagged as duplicated yet
                        {
                            vinfo[verts_bry[j]].Dupl = true;
                            // change elements' connectivities
                            for (size_t k=0; k<verts_bry[j]->Shares.Size(); ++k)
                            {
                                Cell * e = verts_bry[j]->Shares[k].C;
                                int    n = verts_bry[j]->Shares[k].N;
                                e->V[n]  = verts_bry[i];
                                Share sha = {e,(size_t)n};
                                verts_bry[i]->Shares.Push (sha);
                            }
                            ndupl++;
                        }
                    }
                    ncomp++;
                }
            }
        }
    }

    // set new array with non-duplicated vertices
    size_t k = 0;
    Verts.Resize (verts.Size() + verts_m1.Size() + verts_m2.Size() + verts_m3.Size() - ndupl);
    for (size_t i=0; i<verts.Size(); ++i)
    {
        if (vinfo[verts[i]].Dupl==false)
        {
            Verts[k]     = verts[i];
            Verts[k]->ID = k;
            if (vinfo[verts[i]].OnBry)
            {
                if (Verts[k]->Tag<0) TgdVerts.Push (Verts[k]);
            }
            k++;
        }
        else delete verts[i];
    }
    for (size_t i=0; i<verts_m1.Size(); ++i)
    {
        if (vinfo[verts_m1[i]].Dupl==false)
        {
            Verts[k]     = verts_m1[i];
            Verts[k]->ID = k;
            if (vinfo[verts_m1[i]].OnBry)
            {
                if (Verts[k]->Tag<0) TgdVerts.Push (Verts[k]);
            }
            k++;
        }
        else delete verts_m1[i];
    }
    for (size_t i=0; i<verts_m2.Size(); ++i)
    {
        if (vinfo[verts_m2[i]].Dupl==false)
        {
            Verts[k]     = verts_m2[i];
            Verts[k]->ID = k;
            if (vinfo[verts_m2[i]].OnBry)
            {
                if (Verts[k]->Tag<0) TgdVerts.Push (Verts[k]);
            }
            k++;
        }
        else delete verts_m2[i];
    }
    for (size_t i=0; i<verts_m3.Size(); ++i)
    {
        if (vinfo[verts_m3[i]].Dupl==false)
        {
            Verts[k]     = verts_m3[i];
            Verts[k]->ID = k;
            if (vinfo[verts_m3[i]].OnBry)
            {
                if (Verts[k]->Tag<0) TgdVerts.Push (Verts[k]);
            }
            k++;
        }
        else delete verts_m3[i];
    }

    // info
    if (WithInfo)
    {
        printf("\n%s--- Structured Mesh Generation --- %dD --- O%d ---------------------------------------%s\n",TERM_CLR1,NDim,(O2?2:1),TERM_RST);
        printf("%s  Num of comparisons = %ld%s\n", TERM_CLR5, ncomp,        TERM_RST);
        printf("%s  Num of duplicates  = %ld%s\n", TERM_CLR5, ndupl,        TERM_RST);
        printf("%s  Min diagonal       = %g%s\n",  TERM_CLR4, min_diagonal, TERM_RST);
        printf("%s  Num of cells       = %zd%s\n", TERM_CLR2, Cells.Size(), TERM_RST);
        printf("%s  Num of vertices    = %zd%s\n", TERM_CLR2, Verts.Size(), TERM_RST);
    }
}

inline void Structured::GenBox (bool O2, int Nx, int Ny, int Nz, double Lx, double Ly, double Lz)
{
    /*
                      4----------------7
                    ,'|              ,'|
                  ,'  |            ,'  |
                ,'    | -31   -10,'    |
              ,'      |        ,'      |
            5'===============6'        |
            |         |      |    -21  |
            |    -20  |      |         |
            |         0- - - | -  - - -3
            |       ,'       |       ,'
            |     ,' -11     |     ,'
            |   ,'        -30|   ,'
            | ,'             | ,'
            1----------------2'
    */
    Array<Block> blks(1);
    blks[0].Set (/*NDim*/3, /*Tag*/-1, /*NVert*/8,
                 -1.,  0.0, 0.0, 0.0,  // tag, x, y, z
                 -2.,   Lx, 0.0, 0.0, 
                 -3.,   Lx,  Ly, 0.0, 
                 -4.,  0.0,  Ly, 0.0,
                 -5.,  0.0, 0.0,  Lz,  // tag, x, y, z
                 -6.,   Lx, 0.0,  Lz, 
                 -7.,   Lx,  Ly,  Lz, 
                 -8.,  0.0,  Ly,  Lz,
                 -10.,-11.,-20.,-21.,-30.,-31.); // face tags
    blks[0].SetNx (Nx);
    blks[0].SetNy (Ny);
    blks[0].SetNz (Nz);
    NDim = 3;
    Generate (blks,O2);
}

inline void Structured::GenQRing (bool O2, int Nx, int Ny, double r, double R, size_t Nb, double Ax, bool NonLin, char const * WeightsX)
{
    Array<double> Wx;
    if (WeightsX!=NULL)
    {
        double wx;
        std::istringstream iss(WeightsX);
        while (iss>>wx) { Wx.Push(wx); }
    }

    Array<Block> blks(Nb);
    double alp = (Util::PI/2.0)/Nb;
    double bet = alp/2.0;
    double m   = (r+R)/2.0;
    for (size_t i=0; i<Nb; ++i)
    {
        double tag0 = (i==0    ? -10. : 0.);
        double tag1 =            -20.;
        double tag2 = (i==Nb-1 ? -30. : 0.);
        double tag3 =            -40.;
        blks[i].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/8,
                      0.,  r*cos(i*alp)      ,  r*sin(i*alp)     ,
                      0.,  R*cos(i*alp)      ,  R*sin(i*alp)     ,
                      0.,  R*cos((i+1)*alp)  ,  R*sin((i+1)*alp) ,
                      0.,  r*cos((i+1)*alp)  ,  r*sin((i+1)*alp) ,
                      0.,  m*cos(i*alp)      ,  m*sin(i*alp)     ,
                      0.,  R*cos(i*alp+bet)  ,  R*sin(i*alp+bet) ,
                      0.,  m*cos((i+1)*alp)  ,  m*sin((i+1)*alp) ,
                      0.,  r*cos(i*alp+bet)  ,  r*sin(i*alp+bet) ,
                     tag0,tag1,tag2,tag3);
        if (WeightsX!=NULL) blks[i].SetNx (Wx);
        else                blks[i].SetNx (Nx, Ax, NonLin);
        blks[i].SetNy (Ny);
        blks[i].AddArcCtr (/*side*/1, /*x*/0.0, /*y*/0.0, R);
        blks[i].AddArcCtr (/*side*/3, /*x*/0.0, /*y*/0.0, r);
    }
    NDim = 2;
    Generate (blks,O2);
}

inline void Structured::GenQRing (bool O2, int Nx, int Ny, int Nz, double r, double R, double t)
{
    double m = r + (R-r)/2.0;
    double a = r*cos(Util::PI/4.0);
    double b = R*cos(Util::PI/4.0);
    Array<Mesh::Block> blks(1);
    blks[0].Set (3, -1, 20,
                 -1. , 0.0  , r   , 0.0 , //  0
                 -2. , t    , r   , 0.0 , //  1
                 -3. , t    , R   , 0.0 , //  2
                 -4. , 0.0  , R   , 0.0 , //  3

                 -5. , 0.0  , 0.0 , r   , //  4
                 -6. , t    , 0.0 , r   , //  5
                 -7. , t    , 0.0 , R   , //  6
                 -8. , 0.0  , 0.0 , R   , //  7

                  0. , t/2. , r   , 0.0 , //  8
                  0. , t    , m   , 0.0 , //  9
                  0. , t/2. , R   , 0.0 , // 10
                  0. , 0.0  , m   , 0.0 , // 11

                  0. , t/2. , 0.0 , r   , // 12
                  0. , t    , 0.0 , m   , // 13
                  0. , t/2. , 0.0 , R   , // 14
                  0. , 0.0  , 0.0 , m   , // 15

                  0. , 0.0  , a   , a   , // 16
                  0. , t    , a   , a   , // 17
                  0. , t    , b   , b   , // 18
                  0. , 0.0  , b   , b   , // 19

                 -10.,-20.,-30.,-40.,-50.,-60.); // face tags

    blks[0].SetNx (Nx);
    blks[0].SetNy (Ny);
    blks[0].SetNz (Nz);

    blks[0].AddCylXCtr (2, 0.0, 0.0, r);
    blks[0].AddCylXCtr (3, 0.0, 0.0, R);

    NDim = 3;
    Generate (blks, O2);
}

inline void Structured::GenQDisk (bool O2, int Nx1, int Nx2, int Ny, double r, double R)
{
    Array<Block> blks(3);
    double m = r + (R-r)/2.0;
    double n = r/2.0;
    double a = r*cos(Util::PI/4.0);
    double b = r*cos(Util::PI/8.0);
    double c = r*sin(Util::PI/8.0);
    double A = R*cos(Util::PI/4.0);
    double B = R*cos(Util::PI/8.0);
    double C = R*sin(Util::PI/8.0);
    double d = m*cos(Util::PI/4.0);

    blks[0].Set (2, -1, 8,
                 0. , 0.0 , 0.0 , 
                 0. , r   , 0.0 , 
                 0. , a   , a   , 
                 0. , 0.0 , r   , 
                 0. , n   , 0.0 , 
                 0. , b   , c   , 
                 0. , c   , b   , 
                 0. , 0.0 , n   , 
                 -1., -11., -22., -3.);

    blks[1].Set (2, -1, 8,
                 0. , r   , 0.0 , 
                 0. , R   , 0.0 , 
                 0. , A   , A   , 
                 0. , a   , a   , 
                 0. , m   , 0.0 , 
                 0. , B   , C   , 
                 0. , d   , d   , 
                 0. , b   , c   , 
                 -1., -2., -33., -11.);

    blks[2].Set (2, -1, 8,
                 0. , a   , a   , 
                 0. , A   , A   , 
                 0. , 0.0 , R   , 
                 0. , 0.0 , r   , 
                 0. , d   , d   , 
                 0. , C   , B   , 
                 0. , 0.0 , m   , 
                 0. , c   , b   , 
                 -33., -2., -3., -22.);

    blks[0].SetNx (Nx1);
    blks[0].SetNy (Ny );
    blks[1].SetNx (Nx2);
    blks[1].SetNy (Ny );
    blks[2].SetNx (Nx2);
    blks[2].SetNy (Nx1);

    blks[0].AddArcCtr (1, 0.0, 0.0, r);
    blks[0].AddArcCtr (2, 0.0, 0.0, r);
    blks[1].AddArcCtr (1, 0.0, 0.0, R);
    blks[1].AddArcCtr (3, 0.0, 0.0, r);
    blks[2].AddArcCtr (1, 0.0, 0.0, R);
    blks[2].AddArcCtr (3, 0.0, 0.0, r);

    NDim = 2;
    Generate (blks,O2);
}

inline void Structured::GenQPlateHole (double r, double R, double L, int Nx1, int Nx2, int Ny1, int Ny2, bool O2)
{
    /*         ____________________________
     *        |    Ny1   |                 |
     *        |          |                 |
     *        |     3    |       4         |Ny2
     *        |__        |                 |
     *        |   `--,__ |                 |
     *        |    1    `,A ______ q_______|
     *        |__      d/  `,              |
     *            `--,/      \B,C          |
     *                `-   0  \    2       |Ny1
     *   _________   |  \      |           |
     *  c |          |   |     |           |
     *   _|____      |   |_Nx1_|____Nx2____|
     *               | |             p
     *        |---a--| ||   |  |           |
     *                 ||   |  |           |
     *        |---b----||   |  |           |
     *                  |   |  |           |
     *        |----r----|   |  |           |
     *                      |  |           |
     *        |------m------|  |           |
     *                         |           |
     *        |-------R--------|           |
     *                                     |
     *        |-------------L--------------|     */

    if (R-r<0.0) throw new Fatal("Mesh::Structured:GenQPlateHole: R=%g must be greater than r=%g",R,r);
    if (L-R<0.0) throw new Fatal("Mesh::Structured:GenQPlateHole: L=%g must be greater than R=%g",L,R);

    Array<Block> blks(5);
    double m = r + (R-r)/2.0;
    double a = r*cos(Util::PI/4.0);
    double b = r*cos(Util::PI/8.0);
    double c = r*sin(Util::PI/8.0);
    double A = R*cos(Util::PI/4.0);
    double B = R*cos(Util::PI/8.0);
    double C = R*sin(Util::PI/8.0);
    double d = m*cos(Util::PI/4.0);
    double p = R + (L-R)/2.0;
    double q = A + (L-A)/2.0;

    blks[0].Set (2, -1, 8,
                 0. , r   , 0.0 , 
                 0. , R   , 0.0 , 
                 0. , A   , A   , 
                 0. , a   , a   , 
                 0. , m   , 0.0 , 
                 0. , B   , C   , 
                 0. , d   , d   , 
                 0. , b   , c   , 
                 -10., 0., 0., -14.);

    blks[1].Set (2, -1, 8,
                 0. , a   , a   , 
                 0. , A   , A   , 
                 0. , 0.0 , R   , 
                 0. , 0.0 , r   , 
                 0. , d   , d   , 
                 0. , C   , B   , 
                 0. , 0.0 , m   , 
                 0. , c   , b   , 
                 0., 0., -13., -14.);

    blks[2].Set (2, -1, 8,
                 0. , R    , 0.0 , 
                 0. , L    , 0.0 , 
                 0. , L    , A   , 
                 0. , A    , A   , 
                 0. , p    , 0.0 , 
                 0. , L    , A/2., 
                 0. , q    , A   , 
                 0. , B    , C   , 
                 -10., -11., 0., 0.);

    blks[3].Set (2, -1, 8,
                 0. , A    , A   , 
                 0. , A    , L   , 
                 0. , 0.0  , L   , 
                 0. , 0.0  , R   , 
                 0. , A    , q   , 
                 0. , A/2. , L   , 
                 0. , 0.0  , p   , 
                 0. , C    , B   , 
                 0., -12., -13., 0.);

    blks[4].Set (2, -1, 8,
                 0. , A    , A   , 
                 0. , L    , A   , 
                 0. , L    , L   , 
                 0. , A    , L   , 
                 0. , q    , A   , 
                 0. , L    , q   , 
                 0. , q    , L   , 
                 0. , A    , q   , 
                 0., -11., -12., 0.);

    blks[0].SetNx(Nx1);  blks[0].SetNy(Ny1);
    blks[1].SetNx(Nx1);  blks[1].SetNy(Ny1);
    //blks[2].SetNx(Nx2);  blks[2].SetNy(Ny1);
    //blks[3].SetNx(Ny2);  blks[3].SetNy(Ny1);
    //blks[4].SetNx(Nx2);  blks[4].SetNy(Ny2);
    blks[2].SetNx(Nx2, 1.5, true);  blks[2].SetNy(Ny1           );
    blks[3].SetNx(Ny2, 1.5, true);  blks[3].SetNy(Ny1           );
    blks[4].SetNx(Nx2, 1.5, true);  blks[4].SetNy(Ny2, 1.5, true);

    blks[0].AddArcCtr(1, 0.0, 0.0, R);  blks[0].AddArcCtr(3, 0.0, 0.0, r);
    blks[1].AddArcCtr(1, 0.0, 0.0, R);  blks[1].AddArcCtr(3, 0.0, 0.0, r);
    blks[2].AddArcCtr(3, 0.0, 0.0, R);
    blks[3].AddArcCtr(3, 0.0, 0.0, R);

    NDim = 2;
    Generate (blks,O2);
}

inline void Structured::ShapeFunc (double r, double s, double t)
{
    if (NDim==2)
    {
        /*      3           6            2
         *        @---------@----------@
         *        |               (1,1)|
         *        |       s ^          |
         *        |         |          |
         *        |         |          |
         *      7 @         +----> r   @ 5
         *        |       (0,0)        |
         *        |                    |
         *        |                    |
         *        |(-1,-1)             |
         *        @---------@----------@
         *      0           4            1
         */
        double rp1=1.0+r; double rm1=1.0-r;
        double sp1=1.0+s; double sm1=1.0-s;

        N(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
        N(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
        N(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
        N(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
        N(4) = 0.50*sm1*(1.0-r*r);
        N(5) = 0.50*rp1*(1.0-s*s);
        N(6) = 0.50*sp1*(1.0-r*r);
        N(7) = 0.50*rm1*(1.0-s*s);
    }
    else
    {
        /*     t
         *     |
         *     +---s       4         15         7
         *   ,'              @_______@________@
         *  r              ,'|              ,'|     The origin is in the
         *            12 @'  |         14 ,'  |     centre of the cube:
         *             ,'    |16        ,@    |19      0 => (-1,-1)
         *       5   ,'      @      6 ,'      @        6 => ( 1, 1)
         *         @'=======@=======@'        |
         *         |      13 |      |         |
         *         |         |      |  11     |
         *      17 |       0 @______|_@_______@
         *         @       ,'       @       ,' 3
         *         |   8 @'      18 |     ,'
         *         |   ,'           |   ,@ 10
         *         | ,'             | ,'
         *         @_______@________@'
         *        1        9         2
         */
        N( 0) = 0.125*(1.0-r)  *(1.0-s)  *(1.0-t)  *(-r-s-t-2.0);
        N( 1) = 0.125*(1.0+r)  *(1.0-s)  *(1.0-t)  *( r-s-t-2.0);
        N( 2) = 0.125*(1.0+r)  *(1.0+s)  *(1.0-t)  *( r+s-t-2.0);
        N( 3) = 0.125*(1.0-r)  *(1.0+s)  *(1.0-t)  *(-r+s-t-2.0);
        N( 4) = 0.125*(1.0-r)  *(1.0-s)  *(1.0+t)  *(-r-s+t-2.0);
        N( 5) = 0.125*(1.0+r)  *(1.0-s)  *(1.0+t)  *( r-s+t-2.0);
        N( 6) = 0.125*(1.0+r)  *(1.0+s)  *(1.0+t)  *( r+s+t-2.0);
        N( 7) = 0.125*(1.0-r)  *(1.0+s)  *(1.0+t)  *(-r+s+t-2.0);

        N( 8) = 0.25 *(1.0-r*r)*(1.0-s)  *(1.0-t);
        N( 9) = 0.25 *(1.0+r)  *(1.0-s*s)*(1.0-t);
        N(10) = 0.25 *(1.0-r*r)*(1.0+s)  *(1.0-t);
        N(11) = 0.25 *(1.0-r)  *(1.0-s*s)*(1.0-t);

        N(12) = 0.25 *(1.0-r*r)*(1.0-s)  *(1.0+t);
        N(13) = 0.25 *(1.0+r)  *(1.0-s*s)*(1.0+t);
        N(14) = 0.25 *(1.0-r*r)*(1.0+s)  *(1.0+t);
        N(15) = 0.25 *(1.0-r)  *(1.0-s*s)*(1.0+t);

        N(16) = 0.25 *(1.0-r)  *(1.0-s)  *(1.0-t*t);
        N(17) = 0.25 *(1.0+r)  *(1.0-s)  *(1.0-t*t);
        N(18) = 0.25 *(1.0+r)  *(1.0+s)  *(1.0-t*t);
        N(19) = 0.25 *(1.0-r)  *(1.0+s)  *(1.0-t*t);
    }
}

#ifdef USE_BOOST_PYTHON

inline void Structured::PyGenerate (BPy::list const & Blks, bool O2)
{
    size_t nblks = BPy::len(Blks);
    Array<Block> blks(nblks);
    for (size_t i=0; i<nblks; ++i) blks[i].PySet (BPy::extract<BPy::dict>(Blks[i])());
    Generate (blks, O2);
}

#endif

}; // namespace Mesh

#endif // MECHSYS_MESH_STRUCTURED_H
