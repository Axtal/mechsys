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

#ifndef MECHSYS_DEM_FACE_H
#define MECHSYS_DEM_FACE_H

// MechSys
#include <mechsys/dem/edge.h>
#include <mechsys/dem/graph.h>
#include <mechsys/dem/basic_functions.h>
#include <mechsys/util/array.h>

namespace DEM
{

class Face
{
public:
    // Constructor
    Face (Array<Edge *> E); ///< E: Edges of the face
    Face (Array<Vec3_t> & V);
    Face (Array<Vec3_t*> & V);
    Face () {};

    // Destructor
    ~Face ();


    // Methods
    void UpdatedL                   ();                             ///< UdatedL for each edge
    void Normal                     (Vec3_t & N);                   ///< Calculates the normal vecto rof the face
    void Centroid                   (Vec3_t & C);                   ///< Calculates the centroid of the polygonal face
    void GetVec                     (size_t n,Vec3_t & V);          ///< Returns the nth vector to the vector V
    double Area                     ();                             ///< Calculates the area of the face
    bool RayIntersect               (Vec3_t & XO, Vec3_t & X1);     ///< Check if a ray starting at X0 going trough X1 intersects the face

    // Data
    Array<Edge*> Edges;    ///< Edges
    bool         Allocate; ///< It allocates memory or not
    double       Dmax;     ///< Maximun length from face centre
    Vec3_t       Nor;      ///< Normal vector pointing out kept in memory
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Face::Face (Array<Edge *> E)
{
    Edges = E;
    Allocate = false;
    Vec3_t C;
    Centroid (C);
    Dmax = 0.0;
    for (size_t i=0; i<Edges.Size(); i++)
    {
        if (norm(C - *Edges[i]->X0)>Dmax) Dmax = norm(C - *Edges[i]->X0);
    }
    UpdatedL();
    //std::cout << Nor << " " << C << std::endl;
}

inline Face::Face (Array<Vec3_t> & V)
{
    for (size_t i = 0; i < V.Size() ; i++)
    {
        Edges.Push(new Edge(&V[i],&V[(i+1)%V.Size()]));
    }
    Allocate = true;
    Vec3_t C;
    Centroid (C);
    Dmax = 0.0;
    for (size_t i=0; i<Edges.Size(); i++)
    {
        if (norm(C - *Edges[i]->X0)>Dmax) Dmax = norm(C - *Edges[i]->X0);
    }
    UpdatedL();
    //std::cout << Nor << " " << C << std::endl;
}

inline Face::Face(Array<Vec3_t*> & V)
{
    for (size_t i = 0; i < V.Size() ; i++)
    {
        Edges.Push(new Edge(V[i],V[(i+1)%V.Size()]));
    }
    Allocate = true;
    Vec3_t C;
    Centroid (C);
    Dmax = 0.0;
    for (size_t i=0; i<Edges.Size(); i++)
    {
        if (norm(C - *Edges[i]->X0)>Dmax) Dmax = norm(C - *Edges[i]->X0);
    }
    UpdatedL();
    //std::cout << Nor << " " << C << std::endl;
}

inline Face::~Face ()
{
    if (Allocate) 
    {
        for (size_t i = 0; i<Edges.Size();i++)
        {
            delete Edges[i];
        }
    }
}

inline void Face::UpdatedL()
{
    for (size_t i = 0; i<Edges.Size();i++)
    {
        Edges[i]->UpdatedL();
    }
    Normal(Nor);
}

inline void Face::Normal(Vec3_t & N)
{
    N = cross(Edges[0]->dL, Edges[1]->dL);
    N = N/norm(N);
}

inline void Face::Centroid(Vec3_t & N)
{
    N = Vec3_t(0.0,0.0,0.0);
    for (size_t i=0; i<Edges.Size(); i++)
    {
        N += *Edges[i]->X0;
    }
    N/=Edges.Size();
}

inline void Face::GetVec(size_t n, Vec3_t & V)
{
    if (n>=Edges.Size()) throw new Fatal("DEM::Face::GetVec: The selected number is greater than the number of existing vertices");
    V = *Edges[n]->X0;
}

inline double Face::Area()
{
    Vec3_t N;
    Normal(N);
    double area=0;
    for (size_t i=0; i<Edges.Size(); i++)
    {
        area += 0.5*dot(N,cross(*Edges[i]->X0,*Edges[(i+1)%Edges.Size()]->X0));
    }
    return area;
}

inline bool Face::RayIntersect(Vec3_t & X0,Vec3_t & X1)
{
    Vec3_t D = X1 - X0;
    Vec3_t B = X0 - *Edges[0]->X0;
    Mat3_t M;
    M(0,0)   = -D(0);
    M(0,1)   = (Edges[0]->dL)(0);
    M(0,2)   = (Edges[1]->dL)(0);
    M(1,0)   = -D(1);
    M(1,1)   = (Edges[0]->dL)(1);
    M(1,2)   = (Edges[1]->dL)(1);
    M(2,0)   = -D(2);
    M(2,1)   = (Edges[0]->dL)(2);
    M(2,2)   = (Edges[1]->dL)(2);
    Vec3_t X;
    if (!SolAlt(M,B,X)) return false;
    if (X(0)<0.0)       return false;
    B = *Edges[0]->X0 + Edges[0]->dL*X(1) + Edges[1]->dL*X(2);
    Vec3_t nor;
    Normal(nor);
    for (size_t j=0; j<Edges.Size(); j++) 
    {
        Vec3_t tmp = B-*Edges[j]->X0;
        if (dot(cross(Edges[j]->dL,tmp),nor)<0)
        {
            return false;
        }
    }
    return true;
}

}
#endif // MECHSYS_DEM_FACE_H
