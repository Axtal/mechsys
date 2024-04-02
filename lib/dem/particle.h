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

#ifndef MECHSYS_DEM_PARTICLE_H
#define MECHSYS_DEM_PARTICLE_H

// Std lib
#include <iostream>
#include <fstream>
#ifdef USE_OMP
    #include <omp.h>
#endif

// boost => to read json files
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

// MechSys
#include <mechsys/dem/face.h>
#include <mechsys/dem/special_functions.h>
#include <mechsys/util/array.h>
#include <mechsys/numerical/montecarlo.h>
#include <mechsys/mesh/mesh.h>

using std::ifstream;

namespace DEM
{

struct ParticleProps
{
    double  Kn;   ///< Normal stiffness
    double  Kt;   ///< Tengential stiffness
    double  Bn;   ///< Spring constant for normal bonding
    double  Bt;   ///< Spring constant for tangential bonding
    double  Bm;   ///< Spring constant for torque bonding
    double  Gn;   ///< Normal viscous coefficient
    double  Gt;   ///< Tangential viscous coefficient
    double  Gv;   ///< Linear vleocity viscous coefficient
    double  Gm;   ///< Viscous coefficient for the torque
    double  Mu;   ///< Microscopic coefficient of friction
    double  eps;  ///< Maximun strain supported before breaking
    double  Beta; ///< Rolling stiffness coeffcient
    double  Eta;  ///< Plastic moment coefficient
    double  R;    ///< Spheroradius
    double  rho;  ///< Density
    double  V;    ///< Volume
    double  m;    ///< Mass
};

class Particle
{
public:
    // Additional auxiliary methods
    void init_default_values (int tag=-1, double r=0.1, double rho=1.0); ///< Initialises default values
    void poly_calc_props     (Array<Vec3_t> & V);                        ///< Calculates properties for polyhedra

    // Constructor
    Particle() {}
    Particle(int                         Tag,      ///< Tag of the particle
             Array<Vec3_t>       const & V,        ///< List of vertices
             Array<Array <int> > const & E,        ///< List of edges with connectivity
             Array<Array <int> > const & F,        ///< List of faces with connectivity
             Vec3_t              const & v0,       ///< Initial velocity
             Vec3_t              const & w0,       ///< Initial angular velocity
             double                      R,        ///< Spheroradius
             double                      rho=1.0); ///< Density of the material

    // Alternative constructors
    Particle (int Tag, Mesh::Generic const & M, double R, double rho=1.0);                       ///< Generats from mesh file
    Particle (int Tag, char const * TheFileKey, double R, double rho=1.0, double scale = 1.0);   ///< Generates from verts-face pair files

    // Destructor
    ~Particle ();

    // Extra constructor methods
    void ConstructFromJson (int Tag, char const * Filename, double R, double rho=1.0, double scale=1.0); ///< Reads from .msh (json) file format

    // Methods
    void   Initialize         (size_t i=0, size_t NCalls=5000);                               ///< Initialize this particle
    void   InitializeVelocity (double dt = 1.0);                                              ///< Initialize this particle
    void   Rotate             (double dt);                                                    ///< Apply rotation on the particle once the total torque is found
    void   Rotate             (Quaternion_t & Q, Vec3_t & V);                                 ///< Apply rotation given by Quaternion Q at point v
    void   Translate          (double dt);                                                    ///< Apply translation once the total force is found
    void   Translate          (Vec3_t const & t);                                                   ///< Apply translation by vector t
    void   Shrink             (double factor);                                                ///< Shrink the particle around the center of mass 
    void   Erode              (double R);                                                     ///< Erode the particle by a quantity R
    void   Position           (Vec3_t   V);                                                   ///< Position the particle at point V
    void   ResetDisplacements ();                                                             ///< Reset the displacements for the verlet algorithm
    double MaxDisplacement    ();                                                             ///< Maximun displacement for the verlet algorithm
    void   FixVeloc           (double vx=0.0, double vy=0.0, double vz=0.0);                  ///< Fix all velocities
    //bool   IsFree             () {return !vxf&&!vyf&&!vzf&&!wxf&&!wyf&&!wzf;};              ///< Ask if the particle has any constrain in its movement
    bool   IsFree             () {return (!vxf&&!vyf&&!vzf&&!wxf&&!wyf&&!wzf)||FixFree;};                       ///< Ask if the particle has any constrain in its movement
#ifdef USE_OMP
    omp_lock_t      lck;             ///< to protect variables in multithreading
#endif

    int             Tag;             ///< Tag of the particle
    size_t          Index;           ///< index of the particle in the domain
    int             Cluster;         ///< The number of the cohesive cluster the particle belongs to.
    bool            PropsReady;      ///< Are the properties calculated ready ?
    bool            Eroded;          ///< True if the particle has been eroded
    bool            Bdry;            ///< True if the particle is in contact with one of the boundary containers
    bool            Closed;          ///< True if the particle is a closed polyhedron
    bool            vxf, vyf, vzf;   ///< Fixed components of velocity
    bool            wxf, wyf, wzf;   ///< Fixed components of angular velocity
    bool            FixFree;         ///< Fixed to be free, even if there are some constrains 
    Vec3_t          x;               ///< Position of the center of mass
    Vec3_t          xb;              ///< Former position for the Verlet algorithm
    Vec3_t          v;               ///< Velocity
    Vec3_t          w;               ///< Angular velocity
    Vec3_t          wa;              ///< Angular acceleration
    Vec3_t          wb;              ///< Former angular velocity for the leap frog algorithm
    Vec3_t          F;               ///< Force over the particle
    Vec3_t          Flbm;            ///< Force over the particle by the lbm fluid (hydraulic force)
    Vec3_t          Ff;              ///< Fixed Force over the particle
    Vec3_t          T;               ///< Torque over the particle
    Vec3_t          Tf;              ///< Fixed Torque over the particle
    Vec3_t          I;               ///< Vector containing the principal components of the inertia tensor
    Quaternion_t    Q;               ///< The quaternion representing the rotation
    double          Erot;            ///< Rotational energy of the particle
    double          Ekin;            ///< Kinetical energy of the particle
    double          Dmax;            ///< Maximal distance from the center of mass to the surface of the body
    double          Diam;            ///< Diameter of the parallelogram containing the particle
#ifdef USE_CUDA
    size_t          Nvi;             ///< indexes for GPU arrays of vertices, edges and faces
    size_t          Nei;
    size_t          Nfi;
#endif

    ParticleProps       Props;       ///< Properties
    Array<Vec3_t*>      Verts;       ///< Vertices
    Array<Vec3_t*>      Vertso;      ///< Original postion of the Vertices
    Array<Array <int> > EdgeCon;     ///< Conectivity of Edges 
    Array<Array <int> > FaceCon;     ///< Conectivity of Faces 
    Array<Edge*>        Edges;       ///< Edges
    Array<Face*>        Faces;       ///< Faces
    Array<Torus*>       Tori;        ///< Toroidal features
    Array<Cylinder*>    Cylinders;   ///< Cylindrical features

    // Auxiliar methods
    void   CalcProps            (size_t NCalls=5000);                                   ///< Calculate properties: mass, center of mass, and moment of inertia
    bool   IsInside             (Vec3_t & V);                                           ///< Find whether the point V is inside the particle or not
    bool   IsInsideAlt          (Vec3_t & V);                                           ///< Find whether the point V is inside the particle or not
    bool   IsInsideFaceOnly     (Vec3_t & V, Vec3_t const & Per = OrthoSys::O);         ///< Find whether the point V is inside the polyhedron or not
    double IsInside             (double * V);                                           ///< Find whether the point V is inside the particle or not
    double MaxX                 ();                                                     ///< Find Maximun X coordinate
    double MaxY                 ();                                                     ///< Find Maximun Y coordinate
    double MaxZ                 ();                                                     ///< Find Maximun Y coordinate
    double MinX                 ();                                                     ///< Find Minimun X coordinate
    double MinY                 ();                                                     ///< Find Minimun Y coordinate
    double MinZ                 ();                                                     ///< Find Minimun Y coordinate

    // Integrants for the calc of properties
    double Vol (double * X); ///< Calculate the volume of the sample at X
    double Xc  (double * X); ///< Calculate the coordinates of the center of mass at X
    double Yc  (double * X); ///< Calculate the coordinates of the center of mass at X
    double Zc  (double * X); ///< Calculate the coordinates of the center of mass at X
    double Ixx (double * X); ///< Calculate the inertia tensor at X
    double Iyy (double * X); ///< Calculate the inertia tensor at X
    double Izz (double * X); ///< Calculate the inertia tensor at X
    double Ixy (double * X); ///< Calculate the inertia tensor at X
    double Ixz (double * X); ///< Calculate the inertia tensor at X
    double Iyz (double * X); ///< Calculate the inertia tensor at X

};


std::ostream & operator<< (std::ostream & os, Particle const & P)
{
    os << "Tag           = "  << P.Tag        << std::endl;
    os << "Index         = "  << P.Index      << std::endl;
    os << "PropsReady    = "  << P.PropsReady << std::endl;
    os << "vxf, vyf, vzf = "  << P.vxf << ", " << P.vyf << ", " << P.vzf << std::endl;
    os << "wxf, wyf, wzf = "  << P.wxf << ", " << P.wyf << ", " << P.wzf << std::endl;
    os << "x             = "  << PrintVector(P.x );
    os << "xb            = "  << PrintVector(P.xb);
    os << "v             = "  << PrintVector(P.v );
    os << "w             = "  << PrintVector(P.w );
    os << "wb            = "  << PrintVector(P.wb);
    os << "F             = "  << PrintVector(P.F );
    os << "Ff            = "  << PrintVector(P.Ff);
    os << "T             = "  << PrintVector(P.T );
    os << "Tf            = "  << PrintVector(P.Tf);
    os << "I             = "  << PrintVector(P.I );
    os << "Q             = "  << P.Q << std::endl;
    os << "Erot          = "  << P.Erot << std::endl;
    os << "Ekin          = "  << P.Ekin << std::endl;
    os << "Dmax          = "  << P.Dmax << std::endl;
    os << "Diam          = "  << P.Diam << std::endl;
    os << "Kn            = "  << P.Props.Kn   << std::endl; 
    os << "Kt            = "  << P.Props.Kt   << std::endl; 
    os << "Bn            = "  << P.Props.Bn   << std::endl; 
    os << "Bt            = "  << P.Props.Bt   << std::endl; 
    os << "Bm            = "  << P.Props.Bm   << std::endl; 
    os << "Gn            = "  << P.Props.Gn   << std::endl; 
    os << "Gt            = "  << P.Props.Gt   << std::endl; 
    os << "Gv            = "  << P.Props.Gv   << std::endl; 
    os << "Gm            = "  << P.Props.Gm   << std::endl; 
    os << "Mu            = "  << P.Props.Mu   << std::endl; 
    os << "eps           = "  << P.Props.eps  << std::endl; 
    os << "Beta          = "  << P.Props.Beta << std::endl;
    os << "Eta           = "  << P.Props.Eta  << std::endl; 
    os << "R             = "  << P.Props.R    << std::endl; 
    os << "rho           = "  << P.Props.rho  << std::endl; 
    os << "V             = "  << P.Props.V    << std::endl; 
    os << "m             = "  << P.Props.m    << std::endl; 
    return os;
}


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor and destructor

inline void Particle::init_default_values(int tag, double r, double rho)
{
    Tag        = tag;
    Cluster    = 0;
    PropsReady = false;
    Eroded     = false;
    Bdry       = false;
    Closed     = true;
    v          = 0.0,0.0,0.0;
    w          = 0.0,0.0,0.0;

    Props.Kn   = 1.0e4;
    Props.Kt   = 5.0e3;
    Props.Bn   = 1.0e4;
    Props.Bt   = 5.0e3;
    Props.Bm   = 5.0e3;
    Props.Gn   =-0.2;
    Props.Gt   = 0.0;
    Props.Gv   = 0.0;
    Props.Gm   = 0.0;
    Props.Mu   = 0.4;
    Props.eps  = 0.01;
    Props.Beta = 0.12;
    Props.Eta  = 1.0;
    Props.R    = r;
    Props.rho  = rho;

    vxf     = false;
    vyf     = false;
    vzf     = false;
    wxf     = false;
    wyf     = false;
    wzf     = false;
    FixFree = false;

    F      = 0.0,0.0,0.0;
    Flbm   = 0.0,0.0,0.0;
    T      = 0.0,0.0,0.0;
    Ff     = 0.0,0.0,0.0;
    Tf     = 0.0,0.0,0.0;
}

inline void Particle::poly_calc_props(Array<Vec3_t> & V)
{
    double vol; // volume of the polyhedron
    Vec3_t CM;  // Center of mass of the polyhedron
    Mat3_t It;  // Inertia tensor of the polyhedron
    PolyhedraMP(V,FaceCon,vol,CM,It); // Calculate the mass properties of the polyhedron
    x       = CM;
    Props.V = vol;
    Props.m = vol*Props.rho;

    Vec3_t xp,yp,zp;
    Eig(It,I,xp,yp,zp);
    I *= Props.rho;
    CheckDestroGiro(xp,yp,zp);
    Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    Q(1) = (yp(2)-zp(1))/(4*Q(0));
    Q(2) = (zp(0)-xp(2))/(4*Q(0));
    Q(3) = (xp(1)-yp(0))/(4*Q(0));

    Q = Q/norm(Q); // TODO: the other piece of code normalises Q, so I've added it here as well <<<<<<<<<<<<<<<<<

    Dmax = Distance(CM,V[0])+Props.R;
    for (size_t i=1; i<Verts.Size(); ++i)
    {
        if (Distance(CM,*Verts[i])+Props.R > Dmax) Dmax = Distance(CM,*Verts[i])+Props.R;
    }
    Ekin       = 0.0;
    Erot       = 0.0;
    PropsReady = true;
}

inline Particle::Particle (int TheTag, Array<Vec3_t> const & V, Array<Array <int> > const & E, Array<Array <int> > const & Fa, Vec3_t const & v0, Vec3_t const & w0, double TheR, double TheRho)
    : Tag(TheTag), Cluster(0), PropsReady(false), Eroded(false), v(v0), w(w0)
{
    // default values
    init_default_values(TheTag, TheR, TheRho);

    EdgeCon = E;
    FaceCon = Fa;

    for (size_t i=0; i<V.Size(); i++)
    {
        Verts.Push (new Vec3_t(V[i]));
        Vertso.Push (new Vec3_t(V[i]));
    }
    for (size_t i=0; i<Fa.Size(); i++)
    {
        Array<Vec3_t*> verts(Fa[i].Size());
        for (size_t j=0; j<Fa[i].Size(); ++j) verts[j] = Verts[Fa[i][j]];
        Faces.Push (new Face(verts));
    }
    for (size_t i=0; i<E.Size(); i++) Edges.Push (new Edge((*Verts[E[i][0]]), (*Verts[E[i][1]])));
#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
    
}

inline Particle::Particle (int TheTag, Mesh::Generic const & M, double TheR, double TheRho)
    : Tag(TheTag), Cluster(0), PropsReady(false), Eroded(false), v(Vec3_t(0.0,0.0,0.0)), w(Vec3_t(0.0,0.0,0.0))
{
    // default values
    init_default_values(TheTag, TheR, TheRho);

    // check if mesh is Shell
    if (!M.IsShell) throw new Fatal("Particle::Particle: Mesh must be of Shell type");

    // vertices
    size_t nv = M.Verts.Size();
    Array<Vec3_t> V;
    for (size_t i=0; i<nv; ++i)
    {
        Verts .Push (new Vec3_t(M.Verts[i]->C(0), M.Verts[i]->C(1), M.Verts[i]->C(2)));
        Vertso.Push (new Vec3_t(M.Verts[i]->C(0), M.Verts[i]->C(1), M.Verts[i]->C(2)));
        V     .Push (    Vec3_t(M.Verts[i]->C(0), M.Verts[i]->C(1), M.Verts[i]->C(2)));
    }

    // edges and faces
    typedef std::map<std::pair<int,int>,Edge*> Key2Edge_t;
    Key2Edge_t key2edge;        // map edge pair (v0,v1) to Edge* in Edges
    size_t nf = M.Cells.Size(); // number of faces: each cell is one face
    for (size_t i=0; i<nf; ++i)
    {
        // number of vertices per face
        size_t nvf = M.Cells[i]->V.Size();

        // edges
        size_t v0, v1;
        std::pair<int,int> keya, keyb;
        for (size_t j=0; j<nvf; ++j)
        {
            v0   = M.Cells[i]->V[j]->ID;
            v1   = M.Cells[i]->V[(j+1)%nvf]->ID;
            keya = std::make_pair(v0,v1);
            keyb = std::make_pair(v1,v0);
            Key2Edge_t::const_iterator ita = key2edge.find(keya);
            Key2Edge_t::const_iterator itb = key2edge.find(keyb);
            if (ita==key2edge.end() && itb==key2edge.end()) // new edge
            {
                Edges.Push (new Edge((*Verts[v0]), (*Verts[v1])));
                key2edge[keya] = Edges[Edges.Size()-1];
                EdgeCon.Push (Array<int>((int)v0,(int)v1)); // TODO: we may remove this
            }
        }

        // faces
        Array<Vec3_t*> verts(nvf);
        FaceCon.Push (Array<int>()); // TODO: we may remove this
        for (size_t j=0; j<nvf; ++j)
        {
            v0 = M.Cells[i]->V[j]->ID;
            verts[j] = Verts[v0];
            FaceCon[FaceCon.Size()-1].Push (v0); // TODO: we may remove this
        }
        Faces.Push (new Face(verts));
    }

    // calculate properties
    poly_calc_props(V);

#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

inline Particle::Particle(int TheTag, char const * TheFileKey, double TheR, double TheRho, double scale)
    : Tag(TheTag), Cluster(0), PropsReady(false), Eroded(false), v(Vec3_t(0.0,0.0,0.0)), w(Vec3_t(0.0,0.0,0.0))
{
    // default values
    init_default_values(TheTag, TheR, TheRho);

    String fnv(TheFileKey); fnv.append("_verts.mesh");
    String fnf(TheFileKey); fnf.append("_faces.mesh");
    ifstream fnvf(fnv.CStr());
    ifstream fnff(fnf.CStr());
    if (!Util::FileExists(fnv)) throw new Fatal("File <%s> not found",fnv.CStr());
    if (!Util::FileExists(fnf)) throw new Fatal("File <%s> not found",fnf.CStr());

    Array<Vec3_t>            V;
    Array<Array<int> >      Fa;
    Array<Array<int> >       E;
    Array<Array<int> >  VFlist;
    
    size_t ncol=0;
    Vec3_t Vtemp;
    while (!fnvf.eof())
    {
        fnvf >> Vtemp(ncol);
        ncol++;
        ncol = ncol%3;
        if (ncol==0)
        {
            V.Push(scale*Vtemp);
            //std::cout << V[V.Size()-1] << std::endl;
        }
    }
    Array<int> Ftemp;
    ncol = 0;
    while (!fnff.eof())
    {
        size_t tmp;
        fnff >> tmp;
        ncol++;
        if (ncol!=4)
        {
            Ftemp.Push(tmp-1);
        }
        else 
        {
            Fa.Push(Ftemp);
            Ftemp.Resize(0);
            ncol = 0;
        }
    }
    //for (size_t i=0;i<F.Size();i++)
    //{
        //for (size_t j=0;j<F[i].Size();j++)
        //{
            //std::cout << F[i][j] << " ";
        //}
        //std::cout << std::endl;
    //}

    VFlist.Resize(V.Size());

    for (size_t i=0; i<V.Size(); i++)
    {
        Verts .Push (new Vec3_t(V[i]));
        Vertso.Push (new Vec3_t(V[i]));
    }

    for (size_t i=0; i<Fa.Size(); i++)
    {
        Array<Vec3_t*> verts(Fa[i].Size());
        for (size_t j=0; j<Fa[i].Size(); ++j)
        {
            verts[j] = Verts[Fa[i][j]];
            VFlist[Fa[i][j]].Push(i);
        }
        Faces.Push (new Face(verts));
    }

    E.Resize(0);

    for (size_t i = 0; i < V.Size()-1; i++)
    {
        for (size_t j = i+1; j < V.Size(); j++)
        {
            bool first = true;
            for (size_t k = 0; k < VFlist[i].Size(); k++)
            {
            	// Checking if vertex i and j share face k
                if (VFlist[j].Find(VFlist[i][k])!=-1)
                {
                    if (!first)
                    {
                        Array<int> Eaux(2);
                        Eaux[0] = i;
                        Eaux[1] = j;
                        E.Push(Eaux);
                    }
                    first = false;
                }
            }
        }
    }

    for (size_t i=0; i<E.Size(); i++) Edges.Push (new Edge((*Verts[E[i][0]]), (*Verts[E[i][1]])));

    //Array<Array<int> > FaMP;

    for (size_t i=0;i<Faces.Size();i++)
    {
        Vec3_t X0,X1,X2;
        Faces[i]->Centroid(X0);
        Faces[i]->Normal(X1);
        X2 = X0 - X1;
        X1 = X0 + X1;
        size_t ntimesp = 0;
        //size_t ntimesn = 0;
        for (size_t j=0;j<Faces.Size();j++)
        {
            if (i==j) continue;
            if (Faces[j]->RayIntersect(X0,X1)) ntimesp++;
            //if (Faces[j]->RayIntersect(X0,X2)) ntimesn++;
        }
        //if (ntimesp%2==0&&ntimesn%2==0) continue;
        if (ntimesp%2!=0) 
        {           
            Array<Edge *> Etemp(Faces[i]->Edges.Size());
            Array<size_t> Fatemp(Fa[i].Size());
            for (size_t j=0;j<Faces[i]->Edges.Size();j++)
            {
                Vec3_t * P;
                Etemp[j] = Faces[i]->Edges[Faces[i]->Edges.Size()-1-j];
                P = Etemp[j]->X0;
                Etemp[j]->X0 = Etemp[j]->X1;
                Etemp[j]->X1 = P;
                Etemp[j]->UpdatedL();
                Fatemp[j] = Fa[i][Fa[i].Size()-1-j];
            }
            for (size_t j=0;j<Faces[i]->Edges.Size();j++)
            {
                Faces[i]->Edges[j] = Etemp[j];
                Fa[i][j] = Fatemp[j];
            }
        }
        //FaMP.Push(Fa[i]);
    }

    EdgeCon = E;
    FaceCon = Fa;

    // calculate properties
    poly_calc_props(V);

    //for (size_t i=0;i<Faces.Size();i++)
    //{
        //Vec3_t X0,X1;
        //Faces[i]->Centroid(X0);
        //Faces[i]->Normal(X1);
        //X1 = X0 + X1;
        //size_t ntimes = 0;
        //for (size_t j=0;j<Faces.Size();j++)
        //{
            //if (i==j) continue;
            //if (Faces[j]->RayIntersect(X0,X1)) ntimes++;
        //}
        //if (ntimes%2!=0) 
        //{           
            //std::cout << i << std::endl;
        //}
    //}

#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

inline Particle::~Particle()
{
    for (size_t i=0; i<Verts .Size(); ++i) delete Verts[i];
    for (size_t i=0; i<Vertso.Size(); ++i) delete Vertso[i];
    for (size_t i=0; i<Edges .Size(); ++i) delete Edges[i];
    for (size_t i=0; i<Faces .Size(); ++i) delete Faces[i];
}

inline void Particle::ConstructFromJson (int Tag, char const * Filename, double R, double Rho, double scale)
{
    // default values
    init_default_values(Tag, R, Rho);

    // read json file
    Array<Vec3_t> V;

    try {
        boost::property_tree::ptree pt;
        boost::property_tree::read_json(Filename, pt);
        BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("verts")) {
            Vec3_t coords;
            int i = 0;
            BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("c")) {
                coords[i] = scale * boost::lexical_cast<double>(b.second.data());
                i++;
            }
            Verts .Push(new Vec3_t(coords));
            Vertso.Push(new Vec3_t(coords));
            V.Push(coords);
        }
        BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("edges")) {
            Array<int> vids(2);
            int i = 0;
            BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("verts")) {
                vids[i] = boost::lexical_cast<int>(b.second.data());
                i++;
            }
            Edges.Push(new Edge((*Verts[vids[0]]), (*Verts[vids[1]])));
            EdgeCon.Push(vids);

        }
        BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("faces")) {
            Array<int>     vids;
            Array<Vec3_t*> verts;
            BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("verts")) {
                int vid = boost::lexical_cast<int>(b.second.data());
                vids .Push(vid);
                verts.Push(Verts[vid]);
            }
            Faces.Push(new Face(verts));
            FaceCon.Push(vids);
        }
        printf("[1;32mparticle.h: ConstructFromJson: finished[0m\n");
    } catch (std::exception & e) {
        throw new Fatal("particle.h: ConstructFromJson failed:\n\t%s", e.what());
    }

    // calculate properties
    poly_calc_props(V);

#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

// Methods

inline void Particle::Initialize (size_t i, size_t NCalls)
{
    // calc properties
    if (!PropsReady)
    {
        Index = i;
        CalcProps (NCalls);
    }
}

inline void Particle::InitializeVelocity (double dt)
{
    // initialize the particle for the Verlet algorithm
    xb = x-v*dt;
    wb = w;
    Ekin = 0.5*Props.m*dot(v,v);
    Erot = 0.5*(I(0)*w(0)*w(0)+I(1)*w(1)*w(1)+I(2)*w(2)*w(2));
}

inline void Particle::Rotate (double dt)
{
    double q0,q1,q2,q3,wx,wy,wz;
    q0 = 0.5*Q(0);
    q1 = 0.5*Q(1);
    q2 = 0.5*Q(2);
    q3 = 0.5*Q(3);

    Vec3_t Tt = T;

    if (wxf) Tt(0) = 0.0;
    if (wyf) Tt(1) = 0.0;
    if (wzf) Tt(2) = 0.0;

    if (norm(w)>1.0e-12)
    {
        //Vec3_t wn = w/norm(w);
        //Tt -= Props.Gm*w*Vec3_t(I(0)*wn(0)*wn(0),I(1)*wn(1)*wn(1),I(2)*wn(2)*wn(2));
        Tt -= Props.Gm*Vec3_t(I(0)*w(0),I(1)*w(1),I(2)*w(2));
    }

    wa(0)=(Tt(0)+(I(1)-I(2))*wb(1)*wb(2))/I(0);
    wa(1)=(Tt(1)+(I(2)-I(0))*wb(0)*wb(2))/I(1);
    wa(2)=(Tt(2)+(I(0)-I(1))*wb(1)*wb(0))/I(2);
    w = wb+0.5*dt*wa;
    wx = w(0);
    wy = w(1);
    wz = w(2);
    Quaternion_t dq(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz),qm;

    wb  = wb+wa*dt;
    qm  = Q+dq*(0.5*dt);
    q0  = 0.5*qm(0);
    q1  = 0.5*qm(1);
    q2  = 0.5*qm(2);
    q3  = 0.5*qm(3);
    wx  = wb(0);
    wy  = wb(1);
    wz  = wb(2);
    dq  = Quaternion_t(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz);
    Quaternion_t Qd = (qm+dq*0.5*dt),temp;
    Conjugate (Q,temp);
    Rotate    (temp,x);
    Q  = Qd/norm(Qd);
    Rotate (Q,x);
    Erot=0.5*(I(0)*wx*wx+I(1)*wy*wy+I(2)*wz*wz);

}

inline void Particle::Rotate (Quaternion_t & Q,Vec3_t & V)
{
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size(),nc = Cylinders.Size();
    for (size_t i = 0; i < nv; i++)
    {
        Vec3_t xt = *Verts[i]-V;
        Rotation(xt,Q,*Verts[i]);
        *Verts[i] += V;
    }
    for (size_t i = 0; i < ne; i++)
    {
        Edges[i]->UpdatedL();
    }
    for (size_t i = 0; i < nf; i++)
    {
        Faces[i]->UpdatedL();
    }
    for (size_t i = 0; i < nc; i++)
    {
        Cylinders[i]->UpdatedL();
    }
}

inline void Particle::Translate (double dt)
{
    Vec3_t Ft = F;
    if (vxf) Ft(0) = 0.0;
    if (vyf) Ft(1) = 0.0;
    if (vzf) Ft(2) = 0.0;

    //std::cout << "1" << std::endl;
    Ft -= Props.Gv*Props.m*v;

    if(Util::IsNan(norm(Ft))) 
    {
        std::cout << "Index    = " << Index   << std::endl;
        std::cout << "Position = " << x       << std::endl;
        std::cout << "Force    = " << Ft      << std::endl;
        std::cout << "Mass     = " << Props.m << std::endl;
        std::cout << "Inertia  = " << I       << std::endl;
        printf("Particle::Translate: The force is not a number %zd(%d), try reducing the time step \n",Index,Tag);
        throw new Fatal("Particle::Translate: The force is not a number %zd(%d)",Index,Tag);
    }
    Vec3_t temp,xa;
    xa    = 2.0*x - xb + Ft*(dt*dt/Props.m);
    temp  = xa - x;
    v    = 0.5*(xa - xb)/dt;
    xb   = x;
    x    = xa;
    Ekin = 0.5*Props.m*dot(v,v);

    //std::cout << "2" << std::endl;
    size_t nv = Verts.Size();
    for (size_t i = 0; i < nv; i++)
    {
        *Verts[i] += temp;
    }
    //std::cout << "3" << std::endl;
}

inline void Particle::Translate (Vec3_t const & V)
{
    size_t nv = Verts.Size();
    for (size_t i = 0; i < nv; i++)
    {
        *Verts [i] += V;
        *Vertso[i] += V;
    }
    x += V;
    xb += V;
}

inline void Particle::Position  (Vec3_t V)
{
    Vec3_t DV = V - x;
    Translate(DV);
}

inline void Particle::Shrink (double factor)
{
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size(),nc = Cylinders.Size();
    for (size_t i = 0; i < nv; i++)
    {
        *Verts[i] = factor*(*Verts[i] - x) + x;
    }
    for (size_t i = 0; i < ne; i++)
    {
        Edges[i]->UpdatedL();
    }
    for (size_t i = 0; i < nf; i++)
    {
        Faces[i]->UpdatedL();
    }
    for (size_t i = 0; i < nc; i++)
    {
        Cylinders[i]->UpdatedL();
    }
    Props.R *= factor;
    Props.V *= pow(factor,3.0);
    Props.m *= pow(factor,3.0);
    Dmax    *= factor;
    I       *= pow(factor,5.0);
}

inline void Particle::Erode (double R)
{
    Array<Vec3_t> V(Verts.Size());
    for (size_t i=0; i<Verts.Size(); i++)
    {
        V[i] = *Verts[i];
    }
    for (size_t i=0; i<Verts .Size(); ++i) delete Verts[i];
    for (size_t i=0; i<Vertso.Size(); ++i) delete Vertso[i];
    for (size_t i=0; i<Edges .Size(); ++i) delete Edges[i];
    for (size_t i=0; i<Faces .Size(); ++i) delete Faces[i];
    Erosion(V,EdgeCon,FaceCon,R);
    Verts.Resize(0);
    Vertso.Resize(0);
    Faces.Resize(0);
    Edges.Resize(0);
    for (size_t i=0; i<V.Size(); i++)
    {
        Verts.Push (new Vec3_t(V[i]));
        Vertso.Push (new Vec3_t(V[i]));
    }
    for (size_t i=0; i<FaceCon.Size(); i++)
    {
        Array<Vec3_t*> verts(FaceCon[i].Size());
        for (size_t j=0; j<FaceCon[i].Size(); ++j) verts[j] = Verts[FaceCon[i][j]];
        Faces.Push (new Face(verts));
    }
    
    for (size_t i=0; i<EdgeCon.Size(); i++) Edges.Push (new Edge((*Verts[EdgeCon[i][0]]), (*Verts[EdgeCon[i][1]])));
}

inline void Particle::ResetDisplacements ()
{
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        (*Vertso[i]) = (*Verts[i]);
    }
}

inline double Particle::MaxDisplacement ()
{
    double md = 0.0;
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        double mpd = Distance((*Vertso[i]),(*Verts[i]));
        if (mpd>md) md = mpd;
    }
    return md;
}

inline void Particle::FixVeloc (double vx, double vy, double vz)
{
    w   = OrthoSys::O;
    v   = vx, vy, vz;
    vxf = true; vyf = true; vzf = true; 
    wxf = true; wyf = true; wzf = true;
}

// Auxiliar methods

inline void Particle::CalcProps (size_t NCalls)
{
    if (Verts.Size()==1 && Edges.Size()==0 && Faces.Size()==0)
    {
        Props.V = (4./3.)*M_PI*Props.R*Props.R*Props.R;
        I = Vec3_t((8./15.)*M_PI*pow(Props.R,5.),(8./15.)*M_PI*pow(Props.R,5.),(8./15.)*M_PI*pow(Props.R,5.));
        x = *Verts[0];
        Q = 1,0,0,0;
        Props.m = Props.rho*Props.V;
        I*= Props.rho;
        Ekin = 0.5*Props.m*dot(v,v);
        Erot = 0.5*(I(0)*w(0)*w(0)+I(1)*w(1)*w(1)+I(2)*w(2)*w(2));
        Dmax = Props.R;
    }
    else 
    {
        Mat3_t It;
        double Xi[3] = { MinX() , MinY() , MinZ() };
        double Xs[3] = { MaxX() , MaxY() , MaxZ() };
        Numerical::MonteCarlo<Particle> MC(this, Numerical::VEGAS, NCalls);
        Props.V = MC.Integrate(&Particle::Vol, Xi,Xs);
        x(0)    = MC.Integrate(&Particle::Xc,  Xi,Xs)/Props.V;
        x(1)    = MC.Integrate(&Particle::Yc,  Xi,Xs)/Props.V;
        x(2)    = MC.Integrate(&Particle::Zc,  Xi,Xs)/Props.V;
        It(0,0) = MC.Integrate(&Particle::Ixx, Xi,Xs);
        It(1,1) = MC.Integrate(&Particle::Iyy, Xi,Xs);
        It(2,2) = MC.Integrate(&Particle::Izz, Xi,Xs);
        It(1,0) = MC.Integrate(&Particle::Ixy, Xi,Xs);
        It(2,0) = MC.Integrate(&Particle::Ixz, Xi,Xs);
        It(2,1) = MC.Integrate(&Particle::Iyz, Xi,Xs);
        It(0,1) = It(1,0);
        It(0,2) = It(2,0);
        It(1,2) = It(2,1);

        Vec3_t xp,yp,zp;
        Eig(It,I,xp,yp,zp);
        I *= Props.rho;
        CheckDestroGiro(xp,yp,zp);
        Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
        Q(1) = (yp(2)-zp(1))/(4*Q(0));
        Q(2) = (zp(0)-xp(2))/(4*Q(0));
        Q(3) = (xp(1)-yp(0))/(4*Q(0));
        Q = Q/norm(Q);
        Rotation(w,Q,wb);
        w = wb;
        Props.m = Props.rho*Props.V;
        Ekin = 0.5*Props.m*dot(v,v);
        Erot = 0.5*(I(0)*w(0)*w(0)+I(1)*w(1)*w(1)+I(2)*w(2)*w(2));
        Dmax = Distance(x,(*Verts[0]))+Props.R;
        for (size_t i=1; i<Verts.Size(); ++i)
        {
            if (Distance(x,(*Verts[i]))+Props.R > Dmax) Dmax = Distance(x,(*Verts[i]))+Props.R;
        }
    }
    PropsReady = true;
}

inline bool Particle::IsInside (Vec3_t & V)
{
    bool inside = false;
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size();
    if (Distance(x,V)>Dmax) return inside;
    for (size_t i = 0; i < nv; i++)
    {
        if (Distance(V,*Verts[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }

    for (size_t i = 0; i < ne; i++)
    {
        if (Distance(V,*Edges[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }
    for (size_t i = 0; i < nf; i++)
    {
        if (Distance(V,*Faces[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }
    if (nf>3)
    {
        size_t k = 0;
        double Mindistance = Distance(V,*Faces[k]);
        for (size_t i = 1; i < nf; i++)
        {
            if (Distance(V,*Faces[i])<Mindistance) 
            {
                k = i;
                Mindistance = Distance(V,*Faces[k]);
            }
        }
        Vec3_t ct(0,0,0);
        Faces[k]->Centroid(ct);
        Vec3_t pro = V - ct;
        Vec3_t nor;
        Faces[k]->Normal(nor);
        if (dot(pro,nor)<0) inside =true;
    }


    return inside;
}

inline bool Particle::IsInsideAlt(Vec3_t & V)
{
    bool inside = false;
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size();
    if (Distance(x,V)>Dmax) return inside;
    for (size_t i = 0; i < nv; i++)
    {
        if (Distance(V,*Verts[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }

    for (size_t i = 0; i < ne; i++)
    {
        if (Distance(V,*Edges[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }
    for (size_t i = 0; i < nf; i++)
    {
        if (Distance(V,*Faces[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }
    if (nf>3)
    {
        inside = true;
        Vec3_t D = V - x;
        for (size_t i = 0; i < nf; i++)
        {
            Vec3_t B = x - *Faces[i]->Edges[0]->X0;
            Mat3_t M;
            M(0,0) = -D(0);
            M(0,1) = (Faces[i]->Edges[0]->dL)(0);
            M(0,2) = (Faces[i]->Edges[1]->dL)(0);
            M(1,0) = -D(1);
            M(1,1) = (Faces[i]->Edges[0]->dL)(1);
            M(1,2) = (Faces[i]->Edges[1]->dL)(1);
            M(2,0) = -D(2);
            M(2,1) = (Faces[i]->Edges[0]->dL)(2);
            M(2,2) = (Faces[i]->Edges[1]->dL)(2);
            Vec3_t X;
            if (!SolAlt(M,B,X)) continue;
            if (X(0)>1.0||X(0)<0.0) continue;
            B = *Faces[i]->Edges[0]->X0 + Faces[i]->Edges[0]->dL*X(1)+Faces[i]->Edges[1]->dL*X(2);
            Vec3_t nor;
            Faces[i]->Normal(nor);
            bool test = true;
            for (size_t j=0; j<Faces[i]->Edges.Size(); j++) 
            {
                Vec3_t tmp = B-*Faces[i]->Edges[j]->X0;
                if (dot(cross(Faces[i]->Edges[j]->dL,tmp),nor)<0)
                {
                    test = false;
                }
            }
            if (test) inside = false;
        }
    }
    return inside;
}

inline bool Particle::IsInsideFaceOnly(Vec3_t & V, Vec3_t const & Per)
{
    size_t nf = Faces.Size();
    size_t ni = 0;             //Number of intersections
    Vec3_t ref((1.0*rand())/RAND_MAX,(1.0*rand())/RAND_MAX,(1.0*rand())/RAND_MAX);
    ref = Dmax*ref/norm(ref);
    Vec3_t B;
    BranchVec(x,V,B,Per);
    Vec3_t Vt = x + B;
    if (nf>3)
    {
        for (size_t i = 0; i < nf; i++)
        {
            if (Faces[i]->RayIntersect(Vt,ref)) ni++;
        }
    }
    if (ni%2==0) return false;
    else         return true;
}

inline double Particle::IsInside (double * V)
{
    Vec3_t p(V);
    return static_cast<double>(IsInside(p));
}

inline double Particle::MaxX ()
{
    double result = (*Verts[0])(0)+Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(0)+Props.R > result) result = (*Verts[i])(0)+Props.R; 
    }
    return result;
}

inline double Particle::MaxY ()
{
    double result = (*Verts[0])(1)+Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(1)+Props.R > result) result = (*Verts[i])(1)+Props.R; 
    }
    return result;
}

inline double Particle::MaxZ ()
{
    double result = (*Verts[0])(2)+Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(2)+Props.R > result) result = (*Verts[i])(2)+Props.R; 
    }
    return result;
}

inline double Particle::MinX ()
{
    double result = (*Verts[0])(0)-Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(0)-Props.R < result) result = (*Verts[i])(0)-Props.R; 
    }
    return result;
}

inline double Particle::MinY ()
{
    double result = (*Verts[0])(1)-Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(1)-Props.R < result) result = (*Verts[i])(1)-Props.R; 
    }
    return result;
}

inline double Particle::MinZ ()
{
    double result = (*Verts[0])(2)-Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(2)-Props.R < result) result = (*Verts[i])(2)-Props.R; 
    }
    return result;
}

// Integrants for the calculation of properties

inline double Particle::Vol (double * X)
{
    return IsInside(X);
}

inline double Particle::Xc (double * X)
{
    return X[0]*IsInside(X);
}

inline double Particle::Yc (double * X)
{
    return X[1]*IsInside(X);
}

inline double Particle::Zc (double * X)
{
    return X[2]*IsInside(X);
}

inline double Particle::Ixx (double * X)
{
    return ((X[1]-x(1))*(X[1]-x(1))+(X[2]-x(2))*(X[2]-x(2)))*IsInside(X);
}

inline double Particle::Iyy (double * X)
{
    return ((X[0]-x(0))*(X[0]-x(0))+(X[2]-x(2))*(X[2]-x(2)))*IsInside(X);
}

inline double Particle::Izz (double * X)
{
    return ((X[0]-x(0))*(X[0]-x(0))+(X[1]-x(1))*(X[1]-x(1)))*IsInside(X);
}

inline double Particle::Ixy (double * X)
{
    return -(X[0]-x(0))*(X[1]-x(1))*IsInside(X);
}

inline double Particle::Ixz (double * X)
{
    return -(X[0]-x(0))*(X[2]-x(2))*IsInside(X);
}

inline double Particle::Iyz (double * X)
{
    return -(X[1]-x(1))*(X[2]-x(2))*IsInside(X);
}

//////////////////////////////// CUDA IMPLEMENTATION /////////////////////////////

#ifdef USE_CUDA
struct ParticleCU
{
    //Data of particle class
    int            Tag;                                             ///< Tag of the particle
    size_t         Index;                                           ///< index of the particle in the domain
    bool           vxf, vyf, vzf;                                   ///< Fixed components of velocity
    bool           wxf, wyf, wzf;                                   ///< Fixed components of angular velocity
    bool           FixFree;                                         ///< Fixed to be free, even if there are some constrains 
    bool           Closed;                                          ///< Flag to say if the polyhedra is closed or not
    real           R;                                               ///< Spheroradious of particle
    real           m;                                               ///< Mass of particle
    real           Dmax;                                            ///< Maximun Diameter
    real3          Ff;                                              ///< Fixed Force over the particle
    real3          T;                                               ///< Torque over the particle
    real3          Tf;                                              ///< Fixed Torque over the particle
    real3          I;                                               ///< Vector containing the principal components of the inertia tensor
    size_t         Nvi;                                             ///< indexes for GPU arrays of vertices, edges and faces
    size_t         Nvf;                                               
    size_t         Nei;
    size_t         Nef;
    size_t         Nfi;
    size_t         Nff;
};

struct DynParticleCU
{
    //Data for the dynamics of the particle class
    real3          x;                                               ///< Position of the center of mass
    real3          xb;                                              ///< Former position for the Verlet algorithm
    real3          v;                                               ///< Velocity
    real3          w;                                               ///< Angular velocity
    real3          wa;                                              ///< Angular acceleration
    real3          wb;                                              ///< Former angular velocity for the leap frog algorithm
    real3          F;                                               ///< Force over the particle
    real3          Flbm;                                            ///< Force over the particle by the lbm fluid (hydraulic force)
    real4          Q;                                               ///< The quaternion representing the rotation
};

__host__ void UploadParticle(DEM::DynParticleCU & DPc, DEM::ParticleCU & Pcu,DEM::Particle & Par)
{
    {
    Pcu.Tag             = Par.Tag;
    Pcu.Index           = Par.Index;
    Pcu.vxf             = Par.vxf; 
    Pcu.vyf             = Par.vyf;
    Pcu.vzf             = Par.vzf;
    Pcu.wxf             = Par.wxf; 
    Pcu.wyf             = Par.wyf;
    Pcu.wzf             = Par.wzf;
    Pcu.FixFree         = Par.FixFree; 
    Pcu.Closed          = Par.Closed;
    Pcu.R               = Par.Props.R;
    Pcu.m               = Par.Props.m;
    Pcu.Dmax            = Par.Dmax;
    Pcu.Ff.x            = Par.Ff(0);
    Pcu.Ff.y            = Par.Ff(1);
    Pcu.Ff.z            = Par.Ff(2);
    Pcu.T.x             = Par.T(0);
    Pcu.T.y             = Par.T(1);
    Pcu.T.z             = Par.T(2);
    Pcu.Tf.x            = Par.Tf(0);
    Pcu.Tf.y            = Par.Tf(1);
    Pcu.Tf.z            = Par.Tf(2);
    Pcu.I.x             = Par.I(0);
    Pcu.I.y             = Par.I(1);
    Pcu.I.z             = Par.I(2);
    }
    {
    DPc.x.x             = Par.x(0);
    DPc.x.y             = Par.x(1);
    DPc.x.z             = Par.x(2);
    DPc.xb.x            = Par.xb(0);
    DPc.xb.y            = Par.xb(1);
    DPc.xb.z            = Par.xb(2); 
    DPc.v.x             = Par.v(0);
    DPc.v.y             = Par.v(1);
    DPc.v.z             = Par.v(2);
    DPc.w.x             = Par.w(0);
    DPc.w.y             = Par.w(1);
    DPc.w.z             = Par.w(2);
    DPc.wb.x            = Par.wb(0);
    DPc.wb.y            = Par.wb(1);
    DPc.wb.z            = Par.wb(2);
    DPc.wa.x            = Par.wa(0);
    DPc.wa.y            = Par.wa(1);
    DPc.wa.z            = Par.wa(2);
    DPc.F.x             = Par.F(0);
    DPc.F.y             = Par.F(1);
    DPc.F.z             = Par.F(2);
    DPc.Flbm.x          = Par.Flbm(0);
    DPc.Flbm.y          = Par.Flbm(1);
    DPc.Flbm.z          = Par.Flbm(2);
    DPc.Q.x             = Par.Q(1);
    DPc.Q.y             = Par.Q(2);
    DPc.Q.z             = Par.Q(3);
    DPc.Q.w             = Par.Q(0);
    }
}

__host__ void DnloadParticle(DEM::DynParticleCU & DPc,DEM::Particle & Par)
{
    {
    Par.x(0)            = DPc.x.x   ;
    Par.x(1)            = DPc.x.y   ;
    Par.x(2)            = DPc.x.z   ;
    Par.xb(0)           = DPc.xb.x  ;
    Par.xb(1)           = DPc.xb.y  ;
    Par.xb(2)           = DPc.xb.z  ;
    Par.v(0)            = DPc.v.x   ;
    Par.v(1)            = DPc.v.y   ;
    Par.v(2)            = DPc.v.z   ;
    Par.w(0)            = DPc.w.x   ;
    Par.w(1)            = DPc.w.y   ;
    Par.w(2)            = DPc.w.z   ;
    Par.wb(0)           = DPc.wb.x  ;
    Par.wb(1)           = DPc.wb.y  ;
    Par.wb(2)           = DPc.wb.z  ;
    Par.wa(0)           = DPc.wa.x  ;
    Par.wa(1)           = DPc.wa.y  ;
    Par.wa(2)           = DPc.wa.z  ;
    Par.F(0)            = DPc.F.x   ;
    Par.F(1)            = DPc.F.y   ;
    Par.F(2)            = DPc.F.z   ;
    Par.Flbm(0)         = DPc.Flbm.x;
    Par.Flbm(1)         = DPc.Flbm.y;
    Par.Flbm(2)         = DPc.Flbm.z;
    Par.Q(1)            = DPc.Q.x   ;
    Par.Q(2)            = DPc.Q.y   ;
    Par.Q(3)            = DPc.Q.z   ;
    Par.Q(0)            = DPc.Q.w   ;
    }
    Par.Ekin = 0.5*Par.Props.m*dot(Par.v,Par.v);
    Par.Erot = 0.5*(Par.I(0)*Par.w(0)*Par.w(0)+Par.I(1)*Par.w(1)*Par.w(1)+Par.I(2)*Par.w(2)*Par.w(2));
}

__host__ __device__ void Translate(real3 * Verts, real3 * Vertso, DEM::ParticleCU & Pcu, DEM::DynParticleCU & DPc, real3 const & trans)
{
    DPc.xb = DPc.xb + trans;
    DPc.x  = DPc.x  + trans;
    for (size_t iv=Pcu.Nvi;iv<Pcu.Nvf;iv++)
    {
        Verts [iv] = Verts [iv] + trans;
        Vertso[iv] = Vertso[iv] + trans;
    }

}
#endif 

}
#endif // MECHSYS_DEM_PARTICLE_H
