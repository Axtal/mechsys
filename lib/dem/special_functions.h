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

#ifndef MECHSYS_DEM_SPECIAL_H
#define MECHSYS_DEM_SPECIAL_H

// MechSys
#include <mechsys/dem/distance.h>

// Stl
#include <algorithm>

namespace DEM
{

inline void PolyhedraMP(Array<Vec3_t> & V, Array<Array <int> > & F, double & vol, Vec3_t & CM, Mat3_t & It) // Calculate mass properties of general polyhedra
{
    vol = 0.0;
    CM = 0.0,0.0,0.0;
    It(0,0) = 0.0;
    It(1,1) = 0.0;
    It(2,2) = 0.0;
    It(1,0) = 0.0;
    It(2,0) = 0.0;
    It(2,1) = 0.0;
    It(0,1) = 0.0;
    It(0,2) = 0.0;
    It(1,2) = 0.0;
    Array<Face*> Faces;
    for (size_t i=0; i<F.Size(); i++)
    {
        Array<Vec3_t*> verts(F[i].Size());
        for (size_t j=0; j<F[i].Size(); ++j) verts[j] = (&V[F[i][j]]);
        Faces.Push (new Face(verts));
    }
    for (size_t i=0;i<F.Size();i++)
    {
        Vec3_t V0;
        Faces[i]->Centroid(V0);
        for (size_t j=0;j<Faces[i]->Edges.Size();j++)
        {
            Vec3_t V1 = *Faces[i]->Edges[j]->X0;
            Vec3_t V2 = *Faces[i]->Edges[(j+1)%Faces[i]->Edges.Size()]->X0;
            Vec3_t d = cross(Vec3_t(V1-V0),Vec3_t(V2-V0));
            vol += d(2)*(f1(2,V0,V1,V2))/6;
            CM += Vec3_t(d(0)*f2(0,V0,V1,V2),d(1)*f2(1,V0,V1,V2),d(2)*f2(2,V0,V1,V2))/24.0;
            It(0,0) += (d(1)*f3(1,V0,V1,V2)+d(2)*f3(2,V0,V1,V2))/60.0;
            It(1,1) += (d(0)*f3(0,V0,V1,V2)+d(2)*f3(2,V0,V1,V2))/60.0;
            It(2,2) += (d(1)*f3(1,V0,V1,V2)+d(0)*f3(0,V0,V1,V2))/60.0;
            It(1,0) -= d(0)*(V0(1)*g0(0,V0,V1,V2)+V1(1)*g1(0,V0,V1,V2)+V2(1)*g2(0,V0,V1,V2))/120.0;
            It(2,1) -= d(1)*(V0(2)*g0(1,V0,V1,V2)+V1(2)*g1(1,V0,V1,V2)+V2(2)*g2(1,V0,V1,V2))/120.0;
            It(0,2) -= d(2)*(V0(0)*g0(2,V0,V1,V2)+V1(0)*g1(2,V0,V1,V2)+V2(0)*g2(2,V0,V1,V2))/120.0;
        }
    }
    CM/=vol;
    It(0,0) -= (CM(1)*CM(1)+CM(2)*CM(2))*vol;
    It(1,1) -= (CM(0)*CM(0)+CM(2)*CM(2))*vol;
    It(2,2) -= (CM(0)*CM(0)+CM(1)*CM(1))*vol;
    It(1,0) += CM(0)*CM(1)*vol;
    It(0,2) += CM(0)*CM(2)*vol;
    It(2,1) += CM(2)*CM(1)*vol;
    It(0,1)  = It(1,0);
    It(2,0)  = It(0,2);
    It(1,2)  = It(2,1);

}

inline void Erosion(Array<Vec3_t> & V, Array<Array<int> > & E, Array<Array <int> > & F, double R) // Mathematical morphology erosion
{
    if (V.Size()<=3) throw new Fatal("special_functions.h::Erosion There are no enough vertices to work with");
    Array<Face*> Faces;
    for (size_t i=0; i<F.Size(); i++)
    {
        Array<Vec3_t*> verts(F[i].Size());
        for (size_t j=0; j<F[i].Size(); ++j) verts[j] = (&V[F[i][j]]);
        Faces.Push (new Face(verts));
    }
    Array<Vec3_t> Normal(3);
    Array<Array <int> > VFlist;
    Array<Vec3_t> Vtemp;
    //First each trio of faces is moved inward and the intersection is checked
    for (size_t i=0; i<F.Size()-2; i++)
    {
        Normal[0] = cross(Faces[i]->Edges[0]->dL,Faces[i]->Edges[1]->dL);
        Normal[0] = Normal[0]/norm(Normal[0]);
        Normal[0] = *Faces[i]->Edges[0]->X0 - R*Normal[0];
        for (size_t j=i+1; j<F.Size()-1; j++) 
        {
            Normal[1] = cross(Faces[j]->Edges[0]->dL,Faces[j]->Edges[1]->dL);
            Normal[1] = Normal[1]/norm(Normal[1]);
            Normal[1] = *Faces[j]->Edges[0]->X0 - R*Normal[1];
            for (size_t k=j+1; k<F.Size(); k++)
            {
                Normal[2] = cross(Faces[k]->Edges[0]->dL,Faces[k]->Edges[1]->dL);
                Normal[2] = Normal[2]/norm(Normal[2]);
                Normal[2] = *Faces[k]->Edges[0]->X0 - R*Normal[2];
                Mat_t M(9,9);
                Vec_t X(9),B(9);
                for (size_t l = 0;l<9;l++)
                {
                    for (size_t m = 0;m<9;m++)
                    {
                        M(l,m) = 0;
                    }
                    X(l) = 0;
                    B(l) = 0;
                }
                M(0,0) = Faces[i]->Edges[0]->dL(0);
                M(0,1) = Faces[i]->Edges[1]->dL(0);
                M(0,6) = -1;
                M(1,0) = Faces[i]->Edges[0]->dL(1);
                M(1,1) = Faces[i]->Edges[1]->dL(1);
                M(1,7) = -1;
                M(2,0) = Faces[i]->Edges[0]->dL(2);
                M(2,1) = Faces[i]->Edges[1]->dL(2);
                M(2,8) = -1;
                M(3,2) = Faces[j]->Edges[0]->dL(0);
                M(3,3) = Faces[j]->Edges[1]->dL(0);
                M(3,6) = -1;
                M(4,2) = Faces[j]->Edges[0]->dL(1);
                M(4,3) = Faces[j]->Edges[1]->dL(1);
                M(4,7) = -1;
                M(5,2) = Faces[j]->Edges[0]->dL(2);
                M(5,3) = Faces[j]->Edges[1]->dL(2);
                M(5,8) = -1;
                M(6,4) = Faces[k]->Edges[0]->dL(0);
                M(6,5) = Faces[k]->Edges[1]->dL(0);
                M(6,6) = -1;
                M(7,4) = Faces[k]->Edges[0]->dL(1);
                M(7,5) = Faces[k]->Edges[1]->dL(1);
                M(7,7) = -1;
                M(8,4) = Faces[k]->Edges[0]->dL(2);
                M(8,5) = Faces[k]->Edges[1]->dL(2);
                M(8,8) = -1;
                B(0) = -Normal[0](0);
                B(1) = -Normal[0](1);
                B(2) = -Normal[0](2);
                B(3) = -Normal[1](0);
                B(4) = -Normal[1](1);
                B(5) = -Normal[1](2);
                B(6) = -Normal[2](0);
                B(7) = -Normal[2](1);
                B(8) = -Normal[2](2);
                try { Sol(M,B,X); }
                catch (Fatal * fatal) { continue; }
                Vec3_t Inter;
                Inter(0) = X(6);
                Inter(1) = X(7);
                Inter(2) = X(8);
                bool inside = true;
                for (size_t l = 0;l<F.Size();l++)
                {
                    Vec3_t ct(0,0,0);
                    for (size_t m = 0;m<Faces[l]->Edges.Size();m++)
                    {
                        ct += *Faces[l]->Edges[m]->X0;
                    }
                    ct /= Faces[l]->Edges.Size();
                    Vec3_t temp = Inter - ct;
                    Vec3_t N = cross(Faces[l]->Edges[0]->dL,Faces[l]->Edges[1]->dL);
                    if (dot(temp,N)>0) inside = false;
                }
                if (inside) 
                {
                    bool belong = true;
                    for (size_t l = 0;l<F.Size();l++)
                    {
                        if ((Distance(Inter,*Faces[l])<R)&&(l!=i)&&(l!=j)&&(l!=k))
                        {
                            belong = false;
                        }
                    }
                    if (belong)
                    {
                        Vtemp.Push(Inter);
                        Array<int> VFlistaux(0);
                        VFlistaux.Push(i);
                        VFlistaux.Push(j);
                        VFlistaux.Push(k);
                        VFlist.Push(VFlistaux);
                    }
                }
            }
        }
    }
    V = Vtemp;
    if (V.Size()<=3) throw new Fatal("The erosion gave too few vertices to build a convex hull, try a smaller erosion parameter");
    // cm will be a point inside the polyhedron
    Vec3_t cm(0,0,0);
    for (size_t i = 0; i < V.Size(); i++)
    {
        cm += V[i];
    }
    cm /=V.Size();
    // Here the edges are constructed cheking if the vertices share a face
    E.Resize(0);
    for (size_t i = 0; i < V.Size()-1; i++)
    {
        for (size_t j = i+1; j < V.Size(); j++)
        {
        	// Since a pair of vertices share two faces then the edge is created only for the first one found
            bool first = true;
            for (size_t k = 0; k < 3; k++)
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
    // Creating faces, in order to do this we assume that the number of eroded faces is less or equal to the number 
    // of original faces. If a vertex belongs to the original, an eroded face is pushed.
    Array<Array <int> > Ftemp;
    Array <int> Faux;
    for (size_t i = 0;i < F.Size();i++)
    {
        Faux.Resize(0);
        for (size_t j = 0;j<V.Size();j++)
        {
            if (VFlist[j].Find(i)!=-1)
            {
                Faux.Push(j);
            }
        }
        if(Faux.Size()!=0) Ftemp.Push(Faux);
    }
    
    //Now the vertices are ordered to ensure the the normal vector point outwars
    for (size_t i=0;i<Ftemp.Size();i++)
    {
        Array<double> Angles(Ftemp[i].Size());
        Vec3_t ct(0,0,0);
        for (size_t j = 0;j<Ftemp[i].Size();j++)
        {
            ct += V[Ftemp[i][j]];
        }
        ct /= Ftemp[i].Size();
        Vec3_t axis = V[Ftemp[i][0]] - ct;
        Vec3_t inward = cm - ct;
        Angles[0]=0.0;
        for (size_t j = 1;j<Ftemp[i].Size();j++)
        {
            Vec3_t t1 = V[Ftemp[i][j]] - ct;
            double arg = dot(axis,t1)/(norm(axis)*norm(t1));
            if (arg> 1.0) arg= 1.0;
            if (arg<-1.0) arg=-1.0;
            Angles[j] = acos(arg);
            Vec3_t t2 = cross(axis,t1);
            if (dot(t2,inward)>0) Angles[j] = 2*M_PI - Angles[j];
        }
        bool swap = false;
        while (!swap)
        {
            swap = true;
            for (size_t j=0;j<Ftemp[i].Size()-1;j++)
            {
                if (Angles[j]>Angles[j+1])
                {
                    double temp1 = Angles[j];
                    Angles[j] = Angles[j+1];
                    Angles[j+1] = temp1;
                    int temp2 = Ftemp[i][j];
                    Ftemp[i][j] = Ftemp[i][j+1];
                    Ftemp[i][j+1] = temp2;
                    swap = false;
                }
            }
        }
    }
    F = Ftemp;
    for (size_t i=0; i<F.Size(); i++)
    {
        delete Faces[i];
    }
}

inline void Dilation(Array<Vec3_t> & V, Array<Array<int> > & E, Array<Array <int> > & F, Array<Vec3_t> & Vresult, double R) // Mathematical morphology dilation 
{
    if (F.Size()<=3) throw new Fatal("special_functions.h::Dilation There are no enough faces to work with");

    //Detect the set of faces that a vertex belong to
    Array<Array <size_t> > FaceSet(V.Size());
    for (size_t i=0;i<V.Size();i++)
    {
        for (size_t j=0;j<F.Size();j++)
        {
            if (F[j].Has(i)) FaceSet[i].Push(j);
        }
    }

    //Build new vertices, the number of vertices is exactly the same as the original set
    Vresult.Resize(V.Size());
    Array<Vec3_t> Normal(3);

    Array<Face*> Faces;
    for (size_t i=0; i<F.Size(); i++)
    {
        Array<Vec3_t*> verts(F[i].Size());
        for (size_t j=0; j<F[i].Size(); ++j) verts[j] = (&V[F[i][j]]);
        Faces.Push (new Face(verts));
    }
    for (size_t n=0;n<V.Size();n++)
    {
        size_t i = FaceSet[n][0];
        size_t j = FaceSet[n][1];
        size_t k = FaceSet[n][2];
        Normal[0] = cross(Faces[i]->Edges[0]->dL,Faces[i]->Edges[1]->dL);
        Normal[0] = Normal[0]/norm(Normal[0]);
        Normal[0] = *Faces[i]->Edges[0]->X0 + R*Normal[0];
        Normal[1] = cross(Faces[j]->Edges[0]->dL,Faces[j]->Edges[1]->dL);
        Normal[1] = Normal[1]/norm(Normal[1]);
        Normal[1] = *Faces[j]->Edges[0]->X0 + R*Normal[1];
        Normal[2] = cross(Faces[k]->Edges[0]->dL,Faces[k]->Edges[1]->dL);
        Normal[2] = Normal[2]/norm(Normal[2]);
        Normal[2] = *Faces[k]->Edges[0]->X0 + R*Normal[2];
        Mat_t M(9,9);
        Vec_t X(9),B(9);
        M(0,0) = Faces[i]->Edges[0]->dL(0);
        M(0,1) = Faces[i]->Edges[1]->dL(0);
        M(0,6) = -1;
        M(1,0) = Faces[i]->Edges[0]->dL(1);
        M(1,1) = Faces[i]->Edges[1]->dL(1);
        M(1,7) = -1;
        M(2,0) = Faces[i]->Edges[0]->dL(2);
        M(2,1) = Faces[i]->Edges[1]->dL(2);
        M(2,8) = -1;
        M(3,2) = Faces[j]->Edges[0]->dL(0);
        M(3,3) = Faces[j]->Edges[1]->dL(0);
        M(3,6) = -1;
        M(4,2) = Faces[j]->Edges[0]->dL(1);
        M(4,3) = Faces[j]->Edges[1]->dL(1);
        M(4,7) = -1;
        M(5,2) = Faces[j]->Edges[0]->dL(2);
        M(5,3) = Faces[j]->Edges[1]->dL(2);
        M(5,8) = -1;
        M(6,4) = Faces[k]->Edges[0]->dL(0);
        M(6,5) = Faces[k]->Edges[1]->dL(0);
        M(6,6) = -1;
        M(7,4) = Faces[k]->Edges[0]->dL(1);
        M(7,5) = Faces[k]->Edges[1]->dL(1);
        M(7,7) = -1;
        M(8,4) = Faces[k]->Edges[0]->dL(2);
        M(8,5) = Faces[k]->Edges[1]->dL(2);
        M(8,8) = -1;
        B(0) = -Normal[0](0);
        B(1) = -Normal[0](1);
        B(2) = -Normal[0](2);
        B(3) = -Normal[1](0);
        B(4) = -Normal[1](1);
        B(5) = -Normal[1](2);
        B(6) = -Normal[2](0);
        B(7) = -Normal[2](1);
        B(8) = -Normal[2](2);
        Sol(M,B,X);
        Vec3_t Inter;
        Inter(0) = X(6);
        Inter(1) = X(7);
        Inter(2) = X(8);
        Vresult[n] = Inter;
        //std::cout << Inter << std::endl;
    }
    for (size_t i=0; i<F.Size(); i++)
    {
        delete Faces[i];
    }
}

inline double SphereCube (Vec3_t & Xs, Vec3_t & Xc, double R, double dx) // Calculates the volume of intersection between a Cube and a Sphere
{
    Array<Vec3_t> P(8);
    P[0] = Xc - 0.5*dx*OrthoSys::e0 - 0.5*dx*OrthoSys::e1 + 0.5*dx*OrthoSys::e2; 
    P[1] = Xc + 0.5*dx*OrthoSys::e0 - 0.5*dx*OrthoSys::e1 + 0.5*dx*OrthoSys::e2;
    P[2] = Xc + 0.5*dx*OrthoSys::e0 + 0.5*dx*OrthoSys::e1 + 0.5*dx*OrthoSys::e2;
    P[3] = Xc - 0.5*dx*OrthoSys::e0 + 0.5*dx*OrthoSys::e1 + 0.5*dx*OrthoSys::e2;
    P[4] = Xc - 0.5*dx*OrthoSys::e0 - 0.5*dx*OrthoSys::e1 - 0.5*dx*OrthoSys::e2; 
    P[5] = Xc + 0.5*dx*OrthoSys::e0 - 0.5*dx*OrthoSys::e1 - 0.5*dx*OrthoSys::e2;
    P[6] = Xc + 0.5*dx*OrthoSys::e0 + 0.5*dx*OrthoSys::e1 - 0.5*dx*OrthoSys::e2;
    P[7] = Xc - 0.5*dx*OrthoSys::e0 + 0.5*dx*OrthoSys::e1 - 0.5*dx*OrthoSys::e2;
    
    double dmin = 2*R;
    double dmax = 0.0;
    for (size_t j=0;j<P.Size();j++)
    {
        double dist = norm(P[j] - Xs);
        if (dmin>dist) dmin = dist;
        if (dmax<dist) dmax = dist;
    }
    if (dmin > R + dx) return 0.0;
    
    if (dmax < R)
    {
        return 12.0*dx;
    }
    
    double len = 0.0;
    for (size_t j=0;j<4;j++)
    {
        Vec3_t D;
        double a; 
        double b; 
        double c; 
        D = P[(j+1)%4] - P[j];
        a = dot(D,D);
        b = 2*dot(P[j]-Xs,D);
        c = dot(P[j]-Xs,P[j]-Xs) - R*R;
        if (b*b-4*a*c>0.0)
        {
            double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
            double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
            if (ta>1.0&&tb>1.0) continue;
            if (ta<0.0&&tb<0.0) continue;
            if (ta<0.0) ta = 0.0;
            if (tb>1.0) tb = 1.0;
            len += norm((tb-ta)*D);
        }
        D = P[(j+1)%4 + 4] - P[j + 4];
        a = dot(D,D);
        b = 2*dot(P[j + 4]-Xs,D);
        c = dot(P[j + 4]-Xs,P[j + 4]-Xs) - R*R;
        if (b*b-4*a*c>0.0)
        {
            double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
            double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
            if (ta>1.0&&tb>1.0) continue;
            if (ta<0.0&&tb<0.0) continue;
            if (ta<0.0) ta = 0.0;
            if (tb>1.0) tb = 1.0;
            len += norm((tb-ta)*D);
        }
        D = P[j+4] - P[j];
        a = dot(D,D);
        b = 2*dot(P[j]-Xs,D);
        c = dot(P[j]-Xs,P[j]-Xs) - R*R;
        if (b*b-4*a*c>0.0)
        {
            double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
            double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
            if (ta>1.0&&tb>1.0) continue;
            if (ta<0.0&&tb<0.0) continue;
            if (ta<0.0) ta = 0.0;
            if (tb>1.0) tb = 1.0;
            len += norm((tb-ta)*D);
        }
    }
    return len;
}

inline double DiskSquare (Vec3_t & Xs, Vec3_t & Xc, double R, double dx) // Calculates the area of intersection between a Square an a Disk
{
    Array<Vec3_t> P(4);
    P[0] = Xc - 0.5*dx*OrthoSys::e0 - 0.5*dx*OrthoSys::e1 + 0.5*dx*OrthoSys::e2; 
    P[1] = Xc + 0.5*dx*OrthoSys::e0 - 0.5*dx*OrthoSys::e1 + 0.5*dx*OrthoSys::e2;
    P[2] = Xc + 0.5*dx*OrthoSys::e0 + 0.5*dx*OrthoSys::e1 + 0.5*dx*OrthoSys::e2;
    P[3] = Xc - 0.5*dx*OrthoSys::e0 + 0.5*dx*OrthoSys::e1 + 0.5*dx*OrthoSys::e2;
    
    double dmin = 2*R;
    double dmax = 0.0;
    for (size_t j=0;j<P.Size();j++)
    {
        double dist = norm(P[j] - Xs);
        if (dmin>dist) dmin = dist;
        if (dmax<dist) dmax = dist;
    }
    if (dmin > R + dx) return 0.0;
    
    if (dmax < R)
    {
        return 4.0*dx;
    }
    double len = 0.0;
    for (size_t j=0;j<4;j++)
    {
        Vec3_t D = P[(j+1)%4] - P[j];
        double a = dot(D,D);
        double b = 2*dot(P[j]-Xs,D);
        double c = dot(P[j]-Xs,P[j]-Xs) - R*R;
        if (b*b-4*a*c>0.0)
        {
            double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
            double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
            if (ta>1.0&&tb>1.0) continue;
            if (ta<0.0&&tb<0.0) continue;
            if (ta<0.0) ta = 0.0;
            if (tb>1.0) tb = 1.0;
            len += norm((tb-ta)*D);
        }
    }
    return len;
}

inline size_t Pt2idx(iVec3_t const & iv, iVec3_t & Dim) // Calculates the index of the cell at coordinates iv for a cubic lattice of dimensions Dim
{
    return iv(0) + iv(1)*Dim(0) + iv(2)*Dim(0)*Dim(1);
}

inline void   idx2Pt(size_t n, iVec3_t & iv, iVec3_t & Dim) // Calculates the coordinates from the index
{
    iv(0) = n%Dim(0);
    iv(1) = (n/Dim(0))%(Dim(1));
    iv(2) = n/(Dim(0)*Dim(1));
}

}
#endif //MECHSYS_DEM_SPECIAL_H

