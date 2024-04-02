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

#ifndef MECHSYS_DEM_DISTANCE_H
#define MECHSYS_DEM_DISTANCE_H

// MechSys
#include <mechsys/dem/edge.h>
#include <mechsys/dem/face.h>
#include <mechsys/dem/cylinder.h>

namespace DEM
{

/// The Branch functions produce branch vector shifted by the periodic boundaries.

inline void BranchVec (Vec3_t const & V0, Vec3_t const & V1, Vec3_t & Branch, Vec3_t const & Per = OrthoSys::O)
{
    Branch = V1-V0;
    if (fabs(Per(0))>0.0) Branch(0) -= round(Branch(0)/Per(0))*Per(0);
    if (fabs(Per(1))>0.0) Branch(1) -= round(Branch(1)/Per(1))*Per(1);
    if (fabs(Per(2))>0.0) Branch(2) -= round(Branch(2)/Per(2))*Per(2);
}

inline void BranchVecDis(Vec3_t const & V0, Vec3_t const & V1, Vec3_t & Branch, Vec3_t & S, Vec3_t const & Per = OrthoSys::O)
{
    Branch = V1-V0;
    S      = OrthoSys::O;
    if (fabs(Per(0))>0.0) S(0) -= round(Branch(0)/Per(0))*Per(0);
    if (fabs(Per(1))>0.0) S(1) -= round(Branch(1)/Per(1))*Per(1);
    if (fabs(Per(2))>0.0) S(2) -= round(Branch(2)/Per(2))*Per(2);
    Branch += S;
}

/// The Distance functions evaluate the distance between different goemetric features. They give the points Xi and Xf os the points of minimun
/// distance between the geometric features

inline void Distance (Vec3_t const & V, Edge const & E, Vec3_t & Xi, Vec3_t & Xf)
{
    double t = (dot(V,E.dL)-dot((*E.X0),E.dL))/(dot(E.dL,E.dL));
    Xi = V;
    if      (t<0) Xf = (*E.X0);
    else if (t>1) Xf = (*E.X1);
    else          Xf = (*E.X0) + E.dL*t;
}

inline void Distance (Vec3_t const & V, Edge const & E, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    Vec3_t B;
    BranchVec(*E.X0,V,B,Per);
    double t = dot(B,E.dL)/(dot(E.dL,E.dL));
    Xi = V;
    if      (t<0) Xf = (*E.X0);
    else if (t>1) Xf = (*E.X1);
    else          Xf = (*E.X0) + E.dL*t;
    BranchVec(Xi,Xf,S,Per);
}

inline void Distance (Edge const & E, Vec3_t const & V, Vec3_t & Xi, Vec3_t & Xf)
{
    Distance (V,E,Xf,Xi);
}

inline void Distance (Edge const & E, Vec3_t const & V, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    Distance (V,E,Xf,Xi,S,Per);
    S*=-1.0;
}

inline void Distance (Edge const & E0, Edge const & E1, Vec3_t & Xi, Vec3_t & Xf)
{
    double a = dot(E0.dL, (*E0.X0)-(*E1.X0));
    double b = dot(E1.dL, (*E0.X0)-(*E1.X0));
    double c = dot(E0.dL, E0.dL);
    double d = dot(E1.dL, E1.dL);
    double e = dot(E0.dL, E1.dL);
    double t = (c*b-e*a)/(c*d-e*e);
    double s = (e*b-a*d)/(c*d-e*e);
    
    if ((s>0) && (s<1) && (t>0) && (t<1)) 
    {
        Xi = (*E0.X0)+E0.dL*s;
        Xf = (*E1.X0)+E1.dL*t;
    }
    else 
    {
        Vec3_t xi1,xf1;
        Vec3_t xi2,xf2;
        Vec3_t xi3,xf3;
        Vec3_t xi4,xf4;
        Distance ((*E0.X0),E1,xi1,xf1);
        Distance ((*E0.X1),E1,xi2,xf2);
        Distance ((*E1.X0),E0,xf3,xi3);
        Distance ((*E1.X1),E0,xf4,xi4);
        double l1 = norm(xf1-xi1);
        double l2 = norm(xf2-xi2);
        double l3 = norm(xf3-xi3);
        double l4 = norm(xf4-xi4);
        if ((l1<=l2) && (l1<=l3) && (l1<=l4))
        {   
            Xi = xi1;
            Xf = xf1;
        }
        if ((l2<=l1) && (l2<=l3) && (l2<=l4)) 
        {   
            Xi = xi2;
            Xf = xf2;
        }
        if ((l3<=l1) && (l3<=l2) && (l3<=l4)) 
        {   
            Xi = xi3;
            Xf = xf3;
        }
        if ((l4<=l1) && (l4<=l2) && (l4<=l3)) 
        {   
            Xi = xi4;
            Xf = xf4;
        }
    }
}

inline void Distance (Edge const & E0, Edge const & E1, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    Vec3_t B;
    BranchVec(*E1.X0,*E0.X0,B,Per);
    double a = dot(E0.dL, B);
    double b = dot(E1.dL, B);
    double c = dot(E0.dL, E0.dL);
    double d = dot(E1.dL, E1.dL);
    double e = dot(E0.dL, E1.dL);
    double t = (c*b-e*a)/(c*d-e*e);
    double s = (e*b-a*d)/(c*d-e*e);
    
    if ((s>0) && (s<1) && (t>0) && (t<1)) 
    {
        Xi = (*E0.X0)+E0.dL*s;
        Xf = (*E1.X0)+E1.dL*t;
        BranchVec(Xi,Xf,S,Per);
    }
    else 
    {
        Vec3_t xi1,xf1,s1;
        Vec3_t xi2,xf2,s2;
        Vec3_t xi3,xf3,s3;
        Vec3_t xi4,xf4,s4;
        Distance ((*E0.X0),E1,xi1,xf1,s1,Per);
        Distance ((*E0.X1),E1,xi2,xf2,s2,Per);
        Distance (E0,(*E1.X0),xi3,xf3,s3,Per);
        Distance (E0,(*E1.X1),xi4,xf4,s4,Per);
        double l1 = norm(s1);
        double l2 = norm(s2);
        double l3 = norm(s3);
        double l4 = norm(s4);
        if ((l1<=l2) && (l1<=l3) && (l1<=l4))
        {   
            Xi = xi1;
            Xf = xf1;
            S  = s1;
        }
        if ((l2<=l1) && (l2<=l3) && (l2<=l4)) 
        {   
            Xi = xi2;
            Xf = xf2;
            S  = s2;
        }
        if ((l3<=l1) && (l3<=l2) && (l3<=l4)) 
        {   
            Xi = xi3;
            Xf = xf3;
            S  = s3;
        }
        if ((l4<=l1) && (l4<=l2) && (l4<=l3)) 
        {   
            Xi = xi4;
            Xf = xf4;
            S  = s4;
        }
    }
}

inline void Distance (Vec3_t const & V, Face const & F, Vec3_t & Xi, Vec3_t & Xf)
{
    // find the normal to face
    Vec3_t nor = cross(F.Edges[0]->dL, F.Edges[1]->dL);
    nor = nor/norm(nor);

    // find the projection
    double a   = dot(*F.Edges[0]->X0-V, F.Edges[0]->dL);
    double b   = dot(F.Edges[0]->dL   , F.Edges[0]->dL);
    double c   = dot(F.Edges[0]->dL   , F.Edges[1]->dL);
    double d   = dot(*F.Edges[0]->X0-V, F.Edges[1]->dL);
    double f   = dot(F.Edges[1]->dL   , F.Edges[1]->dL);
    double s   = (c*d-a*f)/(b*f-c*c);
    double t   = (a*c-b*d)/(b*f-c*c);
    Vec3_t pro = *F.Edges[0]->X0 + s*F.Edges[0]->dL + t*F.Edges[1]->dL;

    // check if vertex is inside
    bool inside = true;
    for (size_t i=0; i<F.Edges.Size(); i++) 
    {
        Vec3_t tmp = pro-*F.Edges[i]->X0;
        if (dot(cross(F.Edges[i]->dL,tmp),nor)<0) inside = false;
    }

    // get Xi, Xf
    if (inside)
    {
        Xi = V;
        Xf = pro;
    }
    else // compare V against each edge
    {
        Distance (V, (*F.Edges[0]), pro, nor); // pro=Xi, nor=Xf
        double lt = norm(nor-pro);
        double ld = lt;
        Xi = pro;
        Xf = nor;
        for (size_t i=1; i<F.Edges.Size(); i++) 
        {
            Distance (V, (*F.Edges[i]), pro, nor);
            lt = norm(nor-pro);
            if (lt<ld)
            {
                Xi = pro;
                Xf = nor;
                ld = lt;
            }
        }
    }
}

inline void Distance (Vec3_t const & V, Face const & F, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    // find the normal to face
    Vec3_t nor = cross(F.Edges[0]->dL, F.Edges[1]->dL);
    nor = nor/norm(nor);

    // find the projection
    Vec3_t B;
    BranchVec(V,(*F.Edges[0]->X0),B,Per);
    double a   = dot(B, F.Edges[0]->dL);
    double b   = dot(F.Edges[0]->dL   , F.Edges[0]->dL);
    double c   = dot(F.Edges[0]->dL   , F.Edges[1]->dL);
    double d   = dot(B, F.Edges[1]->dL);
    double f   = dot(F.Edges[1]->dL   , F.Edges[1]->dL);
    double s   = (c*d-a*f)/(b*f-c*c);
    double t   = (a*c-b*d)/(b*f-c*c);
    Vec3_t pro = *F.Edges[0]->X0 + s*F.Edges[0]->dL + t*F.Edges[1]->dL;

    // check if vertex is inside
    bool inside = true;
    for (size_t i=0; i<F.Edges.Size(); i++) 
    {
        Vec3_t tmp = pro-*F.Edges[i]->X0;
        if (dot(cross(F.Edges[i]->dL,tmp),nor)<0) inside = false;
    }

    // get Xi, Xf
    if (inside)
    {
        Xi = V;
        Xf = pro;
        BranchVec(Xi,Xf,S,Per);
    }
    else // compare V against each edge
    {
        Vec3_t st;
        Distance (V, (*F.Edges[0]), pro, nor, st, Per); // pro=Xi, nor=Xf
        double lt = norm(st);
        double ld = lt;
        Xi = pro;
        Xf = nor;
        S  = st;
        for (size_t i=1; i<F.Edges.Size(); i++) 
        {
            Distance (V, (*F.Edges[i]), pro, nor, st, Per);
            lt = norm(st);
            if (lt<ld)
            {
                Xi = pro;
                Xf = nor;
                S  = st;
                ld = lt;
            }
        }
    }
}

inline void Distance (Face const & F, Vec3_t const & V, Vec3_t & Xi, Vec3_t & Xf)
{
    Distance (V,F,Xf,Xi);
}

inline void Distance (Face const & F, Vec3_t const & V, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{    
    Distance (V,F,Xf,Xi,S,Per);
    S*=-1.0;
}

inline void Distance (Vec3_t const & V0, Vec3_t const & V1, Vec3_t & Xi, Vec3_t & Xf)
{
    Xi = V0;
    Xf = V1;
}

inline void Distance (Vec3_t const & V0, Vec3_t const & V1, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    Xi = V0;
    Xf = V1;
    BranchVec(Xi,Xf,S,Per);
}

inline void Distance (Vec3_t const & V0, Torus const & T1, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    Xi = V0;
    double theta1 = atan(dot(Xi-*T1.X0,*T1.X2-*T1.X0)/dot(Xi-*T1.X0,*T1.X1-*T1.X0));
    double theta2 = theta1 + M_PI;
    Vec3_t P1,P2;
    P1 = *T1.X0 + cos(theta1)*(*T1.X1-*T1.X0) + sin(theta1)*(*T1.X2-*T1.X0);
    P2 = *T1.X0 + cos(theta2)*(*T1.X1-*T1.X0) + sin(theta2)*(*T1.X2-*T1.X0);
    double dist1 = norm(Xi-P1);
    double dist2 = norm(Xi-P2);
    if (dist1<dist2) Xf = P1;
    else             Xf = P2;
}

inline void Distance (Torus const & T1, Vec3_t const & V0, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    Distance (V0,T1,Xf,Xi,S,Per);
}

inline void Distance (Vec3_t const & V0, Cylinder const & C1, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    Xi = V0;
    Vec3_t xn,yn;
    xn = *C1.T1->X1-*C1.T1->X0;
    yn = *C1.T1->X2-*C1.T1->X0;
    xn/= norm(xn);
    yn/= norm(yn);
    double theta1 = atan(dot(Xi-*C1.T1->X0,yn)/dot(Xi-*C1.T1->X0,xn));
    double theta2 = theta1 + M_PI;
    double DR  = C1.T1->R - C1.T0->R;
    Vec3_t DX  = *C1.T1->X0-*C1.T0->X0;
    double DX2 = dot(DX,DX);
    Vec3_t DV  = V0-*C1.T0->X0;
    double DVDX = dot(DX,DV);
    double s1,s2;
    s1 = (DVDX + DR*dot(DV,cos(theta1)*xn+sin(theta1)*yn) - C1.T0->R*DR)/(DX2+DR*DR);
    s2 = (DVDX + DR*dot(DV,cos(theta2)*xn+sin(theta2)*yn) - C1.T0->R*DR)/(DX2+DR*DR);
    if (s1>1.0) s1 = 1.0;
    if (s1<0.0) s1 = 0.0;
    if (s2>1.0) s2 = 1.0;
    if (s2<0.0) s2 = 0.0;
    Vec3_t P1,P2;
    P1 = C1.X0 + s1*DX + (C1.T0->R + s1*DR)*(cos(theta1)*xn+sin(theta1)*yn);
    P2 = C1.X0 + s2*DX + (C1.T0->R + s2*DR)*(cos(theta2)*xn+sin(theta2)*yn);
    double dist1 = norm(Xi-P1);
    double dist2 = norm(Xi-P2);
    if (dist1<dist2) Xf = P1;
    else             Xf = P2;
}

inline void Distance (Cylinder const & C1, Vec3_t const & V0, Vec3_t & Xi, Vec3_t & Xf, Vec3_t & S, Vec3_t const & Per)
{
    Distance (V0,C1,Xf,Xi,S,Per);
}

/// The following Distance functions return only the distance as a number without the positions of the vectors

inline double Distance (Vec3_t const & V0, Vec3_t const & V1)
{
    return norm(V1-V0);
}

inline double Distance (Vec3_t const & V0, Vec3_t const & V1, Vec3_t const & Per)
{
    Vec3_t dx = V1-V0;
    if (fabs(Per(0))>0.0) dx(0) -= round(dx(0)/Per(0))*Per(0);
    if (fabs(Per(1))>0.0) dx(1) -= round(dx(1)/Per(1))*Per(1);
    if (fabs(Per(2))>0.0) dx(2) -= round(dx(2)/Per(2))*Per(2);
    return norm(dx);
}

inline double Distance (Edge const & E, Vec3_t const & V)
{
    Vec3_t Xi,Xf;
    Distance (E,V,Xi,Xf);
    return norm(Xf-Xi);
}

inline double Distance (Edge const & E, Vec3_t const & V, Vec3_t const & Per)
{
    Vec3_t Xi,Xf,S;
    Distance (E,V,Xi,Xf,S,Per);
    return norm(S);
}

inline double Distance (Vec3_t const & V, Edge const & E)
{
    return Distance (E,V);
}

inline double Distance (Vec3_t const & V, Edge const & E, Vec3_t const & Per)
{
    return Distance (E,V,Per);
}

inline double Distance (Edge const & E0, Edge const & E1)
{
    Vec3_t Xi,Xf;
    Distance (E0,E1,Xi,Xf);
    return norm(Xf-Xi);
}

inline double Distance (Edge const & E0, Edge const & E1, Vec3_t const & Per)
{
    Vec3_t Xi,Xf,S;
    Distance (E0,E1,Xi,Xf,S,Per);
    return norm(S);
}

inline double Distance (Face const & F, Vec3_t const & V)
{
    Vec3_t Xi,Xf;
    Distance (F,V,Xi,Xf);
    return norm(Xf-Xi);
}

inline double Distance (Face const & F, Vec3_t const & V, Vec3_t const & Per)
{
    Vec3_t Xi,Xf,S;
    Distance (F,V,Xi,Xf,S,Per);
    return norm(S);
}

inline double Distance (Vec3_t const & V, Face const & F)
{
    return Distance (F,V);
}

inline double Distance (Vec3_t const & V, Face const & F,Vec3_t const & Per)
{
    return Distance (F,V,Per);
}

inline double Distance (Vec3_t const & V0, Torus const & T1, Vec3_t const & Per)
{
    Vec3_t Xi,Xf,S;
    Distance (V0,T1,Xi,Xf,S,Per);
    return norm(Xf-Xi);
}

inline double Distance (Torus const & T1, Vec3_t const & V0, Vec3_t const & Per)
{
    Vec3_t Xi,Xf,S;
    Distance (T1,V0,Xi,Xf,S,Per);
    return norm(Xf-Xi);
}

inline double Distance (Vec3_t const & V0, Cylinder const & C1, Vec3_t const & Per)
{
    Vec3_t Xi,Xf,S;
    Distance (C1,V0,Xi,Xf,S,Per);
    return norm(Xf-Xi);
}

inline double Distance (Cylinder const & C1, Vec3_t const & V0, Vec3_t const & Per)
{
    Vec3_t Xi,Xf,S;
    Distance (C1,V0,Xi,Xf,S,Per);
    return norm(Xf-Xi);
}

/// The Overlap functions evaluate if two geometric features are close enough to overlap

inline bool Overlap (Vec3_t & V0, Vec3_t & V1, double R0, double R1)
{
    return true;
}

inline bool Overlap (Vec3_t & V0, Edge   & E1, double R0, double R1)
{
    double dist = norm(V0 - 0.5*(*E1.X0 + *E1.X1));
    return (dist < R0 + R1 + E1.Dmax);
}

inline bool Overlap (Vec3_t & V0, Edge   & E1, double R0, double R1, Vec3_t const & Per)
{
    Vec3_t B;
    BranchVec(V0,0.5*(*E1.X0 + *E1.X1),B,Per);
    double dist = norm(B);
    return (dist < R0 + R1 + E1.Dmax);
}

inline bool Overlap (Edge   & E0, Vec3_t & V1, double R0, double R1)
{
    return (Overlap(V1,E0,R1,R0));
}

inline bool Overlap (Edge   & E0, Vec3_t & V1, double R0, double R1, Vec3_t const & Per)
{
    return (Overlap(V1,E0,R1,R0,Per));
}

inline bool Overlap (Edge   & E0, Edge   & E1, double R0, double R1)
{
    double dist = norm(0.5*(*E0.X0 + *E0.X1) - 0.5*(*E1.X0 + *E1.X1));
    return (dist < R0 + R1 + E0.Dmax + E1.Dmax);
}

inline bool Overlap (Edge   & E0, Edge   & E1, double R0, double R1, Vec3_t const & Per)
{
    Vec3_t B;
    BranchVec(0.5*(*E0.X0 + *E0.X1),0.5*(*E1.X0 + *E1.X1),B,Per);
    double dist = norm(B);
    return (dist < R0 + R1 + E0.Dmax + E1.Dmax);
}

inline bool Overlap (Vec3_t & V0, Face   & F1, double R0, double R1)
{
    Vec3_t C = OrthoSys::O;
    for (size_t i=0; i<F1.Edges.Size(); i++)
    {
        C += *F1.Edges[i]->X0;
    }
    C/=F1.Edges.Size();
    double dist = norm(C - V0);
    return (dist < R0 + R1 + F1.Dmax);
}

inline bool Overlap (Vec3_t & V0, Face   & F1, double R0, double R1, Vec3_t const & Per)
{
    Vec3_t C = OrthoSys::O;
    for (size_t i=0; i<F1.Edges.Size(); i++)
    {
        C += *F1.Edges[i]->X0;
    }
    C/=F1.Edges.Size();
    Vec3_t B;
    BranchVec(V0,C,B,Per);
    double dist = norm(B);
    return (dist < R0 + R1 + F1.Dmax);
}

inline bool Overlap (Face   & F0, Vec3_t & V1, double R0, double R1)
{
    return (Overlap(V1,F0,R1,R0));
}

inline bool Overlap (Face   & F0, Vec3_t & V1, double R0, double R1, Vec3_t const & Per)
{
    return (Overlap(V1,F0,R1,R0,Per));
}

inline bool Overlap (Vec3_t & V0, Torus  & T1, double R0, double R1, Vec3_t const & Per)
{
    double dist = norm(V0 - *T1.X0);
    return (dist < R0 + R1 + T1.R);
}

inline bool Overlap (Torus  & T0, Vec3_t & V1, double R0, double R1, Vec3_t const & Per)
{
    return (Overlap(V1,T0,R1,R0,Per));
}

inline bool Overlap (Vec3_t & V0, Cylinder & C1, double R0, double R1, Vec3_t const & Per)
{
    double dist = norm(V0 - 0.5*(*C1.T0->X0 + *C1.T1->X0));
    return (dist < R0 + R1 + C1.Dmax);
}

inline bool Overlap (Cylinder & C0, Vec3_t & V1, double R0, double R1, Vec3_t const & Per)
{
    return (Overlap(V1,C0,R1,R0,Per));
}

#ifdef USE_CUDA
/////////////////////////////DEM CUDA implementation////////////////////

__host__ __device__ void BranchVec (real3 const & V0, real3 const & V1, real3 & Branch, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    Branch = V1-V0;
    if (fabs(Per.x)>0.0) Branch.x -= round(Branch.x/Per.x)*Per.x;
    if (fabs(Per.y)>0.0) Branch.y -= round(Branch.y/Per.y)*Per.y;
    if (fabs(Per.z)>0.0) Branch.z -= round(Branch.z/Per.z)*Per.z;
}

__host__ __device__ void BranchVecDis (real3 const & V0, real3 const & V1, real3 & Branch, real3 & S, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    Branch = V1-V0;
    S      = make_real3(0.0,0.0,0.0);
    if (fabs(Per.x)>0.0) S.x -= round(Branch.x/Per.x)*Per.x;
    if (fabs(Per.y)>0.0) S.y -= round(Branch.y/Per.y)*Per.y;
    if (fabs(Per.z)>0.0) S.z -= round(Branch.z/Per.z)*Per.z;
    Branch = Branch + S;
}

__host__ __device__ bool OverlapEE(size_t const * Edges, real3 const * Verts, size_t e0, size_t e1, real const & Dmax0, real const & Dmax1, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    real3 x0   = 0.5*(Verts[Edges[2*e0+1]]+Verts[Edges[2*e0  ]]);
    real3 x1   = 0.5*(Verts[Edges[2*e1+1]]+Verts[Edges[2*e1  ]]);
    real3 B;
    BranchVec(x0,x1,B,Per);
    return (norm(B) < Dmax0 + Dmax1);
}

__host__ __device__ bool OverlapVF(size_t const * Faces, size_t const * Facid, real3 const * Verts, real3 const & V0, size_t f1,real const & Dmax0, real const & Dmax1, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    real3 x1 = make_real3(0.0,0.0,0.0);
    size_t nvertex = Faces[2*f1+1] - Faces[2*f1];
    for (size_t i=0;i<nvertex;i++)
    {
        size_t ig  = i             + Faces[2*f1];
        x1 = x1 + Verts[Facid[ig]];
    }
    x1 = x1/nvertex;
    real3 B;
    BranchVec(V0,x1,B,Per);
    return (norm(B) < Dmax0 + Dmax1);
}

__host__ __device__ bool OverlapFV(size_t const * Faces, size_t const * Facid, real3 const * Verts, size_t f0, real3 const & V1, real const & Dmax0, real const & Dmax1, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    real3 x0 = make_real3(0.0,0.0,0.0);
    size_t nvertex = Faces[2*f0+1] - Faces[2*f0];
    for (size_t i=0;i<nvertex;i++)
    {
        size_t ig  = i             + Faces[2*f0];
        x0 = x0 + Verts[Facid[ig]];
    }
    x0 = x0/nvertex;
    real3 B;
    BranchVec(V1,x0,B,Per);
    return (norm(B) < Dmax0 + Dmax1);
}

__host__ __device__ void DistanceVE(real3 const & e0, real3 const & e1, real3 const & v0, real3 & S, real3 & Xf, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    real3 B;
    real3 dL = e1 - e0;
    BranchVec(e0,v0,B,Per);
    real t = dotreal3(B,dL)/(dotreal3(dL,dL));
    if      (t<0) Xf = e0;
    else if (t>1) Xf = e1;
    else          Xf = e0 + t*dL;
    BranchVec(v0,Xf,S,Per);
}

__host__ __device__ void DistanceEE(size_t const * Edges, real3 const * Verts, size_t e0, size_t e1, real3 & Xi, real3 & Xf, real3 & S, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    real3 B;
    BranchVec(Verts[Edges[2*e1  ]],Verts[Edges[2*e0  ]],B,Per);
    real3 dL0 = Verts[Edges[2*e0+1]] - Verts[Edges[2*e0  ]];
    real3 dL1 = Verts[Edges[2*e1+1]] - Verts[Edges[2*e1  ]];
    real a = dotreal3(dL0, B);
    real b = dotreal3(dL1, B);
    real c = dotreal3(dL0, dL0);
    real d = dotreal3(dL1, dL1);
    real e = dotreal3(dL0, dL1);
    real t = (c*b-e*a)/(c*d-e*e);
    real s = (e*b-a*d)/(c*d-e*e);
    
    if ((s>0) && (s<1) && (t>0) && (t<1)) 
    {
        Xi = (Verts[Edges[2*e0  ]])+s*dL0;
        Xf = (Verts[Edges[2*e1  ]])+t*dL1;
        BranchVec(Xi,Xf,S,Per);
    }
    else 
    {
        real3 xi1,xf1,s1;
        real3 xi2,xf2,s2;
        real3 xi3,xf3,s3;
        real3 xi4,xf4,s4;
        xi1   = Verts[Edges[2*e0  ]];
        xi2   = Verts[Edges[2*e0+1]];
        xi3   = Verts[Edges[2*e1  ]];
        xi4   = Verts[Edges[2*e1+1]];
        DistanceVE(Verts[Edges[2*e1  ]],Verts[Edges[2*e1+1]],xi1,s1,xf1,Per);
        DistanceVE(Verts[Edges[2*e1  ]],Verts[Edges[2*e1+1]],xi2,s2,xf2,Per);
        DistanceVE(Verts[Edges[2*e0  ]],Verts[Edges[2*e0+1]],xi3,s3,xf3,Per);
        DistanceVE(Verts[Edges[2*e0  ]],Verts[Edges[2*e0+1]],xi4,s4,xf4,Per);
        real l1 = norm(s1);
        real l2 = norm(s2);
        real l3 = norm(s3);
        real l4 = norm(s4);
        if ((l1<=l2) && (l1<=l3) && (l1<=l4))
        {   
            Xi = xi1;
            Xf = xf1;
            S  = s1;
        }
        if ((l2<=l1) && (l2<=l3) && (l2<=l4)) 
        {   
            Xi = xi2;
            Xf = xf2;
            S  = s2;
        }
        if ((l3<=l1) && (l3<=l2) && (l3<=l4)) 
        {   
            Xi = xf3;
            Xf = xi3;
            S  = -1.0*s3;
        }
        if ((l4<=l1) && (l4<=l2) && (l4<=l3)) 
        {   
            Xi = xf4;
            Xf = xi4;
            S  = -1.0*s4;
        }
    }
}

__host__ __device__ void DistanceVF(size_t const * Faces, size_t const * Facid, real3 const * Verts, real3 const & V0, size_t f1, real3 & Xf, real3 & S, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    // find the normal to face
    real3  dL0 = Verts[Facid[Faces[2*f1]+1]]-Verts[Facid[Faces[2*f1]  ]];
    real3  dL1 = Verts[Facid[Faces[2*f1]+2]]-Verts[Facid[Faces[2*f1]+1]];
    real3  nor = cross(dL0, dL1);
    nor = nor/norm(nor);

    // find the projection
    real3 B;
    BranchVec(V0,Verts[Facid[Faces[2*f1]  ]],B,Per);
    real  a   = dotreal3(B  , dL0);
    real  b   = dotreal3(dL0, dL0);
    real  c   = dotreal3(dL0, dL1);
    real  d   = dotreal3(B  , dL1);
    real  f   = dotreal3(dL1, dL1);
    real  s   = (c*d-a*f)/(b*f-c*c);
    real  t   = (a*c-b*d)/(b*f-c*c);
    real3 pro = Verts[Facid[Faces[2*f1]  ]] + s*dL0 + t*dL1;

    // check if vertex is inside
    bool inside = true;
    size_t nvertex = Faces[2*f1+1] - Faces[2*f1];
    for (size_t i=0;i<nvertex;i++)
    {
        size_t ig  = i             + Faces[2*f1];
        size_t ip  = (i+1)%nvertex + Faces[2*f1];
        real3  tmp = pro-Verts[Facid[ig]];
        real3  dL  = Verts[Facid[ip]]-Verts[Facid[ig]];
        if (dotreal3(cross(dL,tmp),nor)<0) inside = false;
    }

    // get Xi, Xf
    if (inside)
    {
        Xf = pro;
        BranchVec(V0,Xf,S,Per);
    }
    else // compare V against each edge
    {
        real3 st;
        real3 E0 = Verts[Facid[Faces[2*f1]  ]];
        real3 E1 = Verts[Facid[Faces[2*f1]+1]];
        DistanceVE(E0,E1,V0,st,nor,Per);
        double lt = norm(st);
        double ld = lt;
        Xf = nor;
        S  = st;
        for (size_t i=0;i<nvertex;i++)
        {
            size_t ig  = i             + Faces[2*f1];
            size_t ip  = (i+1)%nvertex + Faces[2*f1];
            E0 = Verts[Facid[ig]];
            E1 = Verts[Facid[ip]];
            DistanceVE(E0,E1,V0,st,nor,Per);
            lt = norm(st);
            if (lt<ld)
            {
                Xf = nor;
                S  = st;
                ld = lt;
            }
        }
    }
}

__host__ __device__ void DistanceFV(size_t const * Faces, size_t const * Facid, real3 const * Verts, size_t f0, real3 const & V1, real3 & Xi, real3 & S, real3 const & Per = make_real3(0.0,0.0,0.0))
{
    DistanceVF(Faces, Facid, Verts, V1, f0, Xi, S, Per);
    S = -1.0*S;
}
#endif

}
#endif // MECHSYS_DEM_DISTANCE_H
