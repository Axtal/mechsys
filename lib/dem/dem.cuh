/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2023 Sergio Galindo                                    *
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

/////////////////////////////DEM CUDA implementation////////////////////

#ifndef MECHSYS_DEM_CUH
#define MECHSYS_DEM_CUH

//Mechsys
#include <mechsys/dem/interacton.h>
namespace DEM
{

struct dem_aux
{
    size_t ncoint = 0;   ///< number of common interactons
    size_t nverts = 0;   ///< number of total vertices
    size_t nfacid = 0;   ///< number of vertices list for faces
    size_t nfaces = 0;   ///< number of total faces
    size_t nedges = 0;   ///< number of total edges
    size_t nfvint = 0;   ///< number of total face vertex interactions
    size_t nvfint = 0;   ///< number of total vertex face interactions
    size_t neeint = 0;   ///< number of total edge edge interactions
    size_t nvvint = 0;   ///< number of total vertex vertex interactions
    size_t nparts = 0;   ///< number of particles
    real   dt;           ///< Time step
    real   Time;         ///< Time clock
    size_t iter   = 0;   ///< Iteration clock;
    real3  Per;          ///< Vector with the periodic boundary condition information
    real   Xmin   = 0.0;   
    real   Ymin   = 0.0;
    real   Zmin   = 0.0;
    real   Xmax   = 0.0;
    real   Ymax   = 0.0;
    real   Zmax   = 0.0;
    bool   px     = false;
    bool   py     = false;
    bool   pz     = false;

};

typedef void (*ForceVV_ptr_t)(InteractonCU const *, ComInteractonCU *, DynInteractonCU *, ParticleCU *
        , DynParticleCU *, dem_aux const *, void *);

typedef void (*ForceEE_ptr_t)(size_t const *, real3 const *, InteractonCU const *, ComInteractonCU *, DynInteractonCU *, ParticleCU *
        , DynParticleCU *, dem_aux const *, void *);

typedef void (*ForceVF_ptr_t)(size_t const *, size_t const *, real3 const *, InteractonCU const *, ComInteractonCU *, DynInteractonCU *, ParticleCU *
        , DynParticleCU *, dem_aux const *, void *);

typedef void (*ForceFV_ptr_t)(size_t const *, size_t const *, real3 const *, InteractonCU const *, ComInteractonCU *, DynInteractonCU *, ParticleCU *
        , DynParticleCU *, dem_aux const *, void *);

typedef void (*Translate_ptr_t)(real3 *, ParticleCU const *, DynParticleCU *, dem_aux const *, void *);

typedef void (*Rotate_ptr_t)   (real3 *, ParticleCU const *, DynParticleCU *, dem_aux const *, void *);

typedef void (*Reset_ptr_t)    (ParticleCU *, DynParticleCU *, InteractonCU const * , ComInteractonCU *, dem_aux const *, void *);

__global__ void CalcForceVV(InteractonCU const * __restrict__ Int, ComInteractonCU * __restrict__ CInt, DynInteractonCU * __restrict__ DIntVV, ParticleCU * __restrict__ Par,
        DynParticleCU * __restrict__ DPar, dem_aux const * __restrict__ demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nvvint) return;
    size_t id = DIntVV[ic].Idx;
    size_t i1 = CInt [id].I1;
    size_t i2 = CInt [id].I2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    real3  Branch;
    if (Int[id].BothFree) BranchVec(xf,xi,Branch,demaux[0].Per);
    else Branch = xi-xf;

    real dist  = norm(Branch);
    real delta = r1 + r2 - dist;

    DIntVV[ic].Fn = make_real3(0.0,0.0,0.0);

    if (delta>0.0)
    {
        real3  n   = -1.0*Branch/dist;
        real   d   = (r1*r1-r2*r2+dist*dist)/(2.0*dist);
        real3  x1c = xi+d*n;
        real3  x2c = xf-(dist-d)*n;

        real3 t1,t2,x1,x2;
        Rotation(DPar[i1].w,DPar[i1].Q,t1);
        Rotation(DPar[i2].w,DPar[i2].Q,t2);
        x1 = x1c- xi;
        x2 = x2c- xf;
        real3 vrel = (DPar[i1].v+cross(t1,x1))-(DPar[i2].v+cross(t2,x2));
        real3 vt   = vrel - dotreal3(n,vrel)*n;
        
        DIntVV[ic].Fn  = Int[id].Kn*delta*n;
        DIntVV[ic].Ft  = DIntVV[ic].Ft + (Int[id].Kt*demaux[0].dt)*vt;
        DIntVV[ic].Ft  = DIntVV[ic].Ft - dotreal3(DIntVV[ic].Ft,n)*n;

        real  Ftn = norm(DIntVV[ic].Ft);
        real  Fnn = norm(DIntVV[ic].Fn);
        real3 tan = (Ftn>0.0) ? DIntVV[ic].Ft/Ftn : DIntVV[ic].Ft;
        if (Ftn>Int[id].Mu*Fnn)
        {
            DIntVV[ic].Ft = Int[id].Mu*Fnn*tan;
        }

        real3 vr = r1*r2*cross((t1 - t2),n)/(r1+r2);
        DIntVV[ic].Fr  = DIntVV[ic].Fr + (Int[id].Beta*Int[id].Kt*demaux[0].dt)*vr;
        DIntVV[ic].Fr  = DIntVV[ic].Fr - dotreal3(DIntVV[ic].Fr,n)*n;

        real  Frn = norm(DIntVV[ic].Fr);
        tan = (Frn>0.0) ? DIntVV[ic].Fr/Frn : DIntVV[ic].Fr;
        if (Frn>Int[id].Eta*Int[id].Mu*Fnn)
        {
            DIntVV[ic].Fr = Int[id].Eta*Int[id].Mu*Fnn*tan;
        }
        
        DIntVV[ic].F = DIntVV[ic].Fn + DIntVV[ic].Ft + Int[id].Gn*dotreal3(n,vrel)*n + Int[id].Gt*vt;

        real3 T1,T2,T, Tt;
        Tt = cross (x1,DIntVV[ic].F) + r1*cross(n,DIntVV[ic].Fr);
        real4 q;
        Conjugate (DPar[i1].Q,q);
        Rotation  (Tt,q,T);
        T1 = -1.0*T;
        Tt = cross (x2,DIntVV[ic].F) + r2*cross(n,DIntVV[ic].Fr);
        Conjugate (DPar[i2].Q,q);
        Rotation  (Tt,q,T);
        T2 =      T;

        atomicAdd(&CInt[id].Fnnet.x, DIntVV[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntVV[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntVV[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntVV[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntVV[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntVV[ic].Ft.z);

        atomicAdd(&DPar[i1].F.x,-DIntVV[ic].F .x);
        atomicAdd(&DPar[i1].F.y,-DIntVV[ic].F .y);
        atomicAdd(&DPar[i1].F.z,-DIntVV[ic].F .z);
        atomicAdd(&DPar[i2].F.x, DIntVV[ic].F .x);
        atomicAdd(&DPar[i2].F.y, DIntVV[ic].F .y);
        atomicAdd(&DPar[i2].F.z, DIntVV[ic].F .z);
        atomicAdd(& Par[i1].T.x,            T1.x);
        atomicAdd(& Par[i1].T.y,            T1.y);
        atomicAdd(& Par[i1].T.z,            T1.z);
        atomicAdd(& Par[i2].T.x,            T2.x);
        atomicAdd(& Par[i2].T.y,            T2.y);
        atomicAdd(& Par[i2].T.z,            T2.z);
    }
}

__global__ void CalcForceVV_Hertz(InteractonCU const * __restrict__ Int, ComInteractonCU * __restrict__ CInt, DynInteractonCU * __restrict__ DIntVV, ParticleCU * __restrict__ Par,
        DynParticleCU * __restrict__ DPar, dem_aux const * __restrict__ demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nvvint) return;

    size_t id = DIntVV[ic].Idx;
    size_t i1 = CInt [id].I1;
    size_t i2 = CInt [id].I2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    real3  Branch;

    if (Int[id].BothFree) BranchVec(xf,xi,Branch,demaux[0].Per);
    else Branch = xi-xf;

    real dist  = norm(Branch);
    real delta = r1 + r2 - dist;

    DIntVV[ic].Fn = make_real3(0.0,0.0,0.0);

    if (delta>0.0)
    {
        real3  n   = -1.0*Branch/dist;
        real   d   = (r1*r1-r2*r2+dist*dist)/(2.0*dist);
        real3  x1c = xi+d*n;
        real3  x2c = xf-(dist-d)*n;

        real3 t1,t2,x1,x2;
        Rotation(DPar[i1].w,DPar[i1].Q,t1);
        Rotation(DPar[i2].w,DPar[i2].Q,t2);
        x1 = x1c- xi;
        x2 = x2c- xf;
        real3 vrel = (DPar[i1].v+cross(t1,x1))-(DPar[i2].v+cross(t2,x2));
        real3 vt   = vrel - dotreal3(n,vrel)*n;

        real sqrtdelta = sqrt(delta);
        DIntVV[ic].Fn  = Int[id].Kn*sqrtdelta*delta*n;
        DIntVV[ic].Ft  = DIntVV[ic].Ft + (Int[id].Kt*sqrtdelta*demaux[0].dt)*vt;
        DIntVV[ic].Ft  = DIntVV[ic].Ft - dotreal3(DIntVV[ic].Ft,n)*n;

        real  Ftn = norm(DIntVV[ic].Ft);
        real  Fnn = norm(DIntVV[ic].Fn);
        real3 tan = (Ftn>0.0) ? DIntVV[ic].Ft/Ftn : DIntVV[ic].Ft;
        if (Ftn>Int[id].Mu*Fnn)
        {
            DIntVV[ic].Ft = Int[id].Mu*Fnn*tan;
        }

        real3 vr = r1*r2*cross((t1 - t2),n)/(r1+r2);
        DIntVV[ic].Fr  = DIntVV[ic].Fr + (Int[id].Beta*Int[id].Kt*sqrtdelta*demaux[0].dt)*vr;
        DIntVV[ic].Fr  = DIntVV[ic].Fr - dotreal3(DIntVV[ic].Fr,n)*n;

        real  Frn = norm(DIntVV[ic].Fr);
        tan = (Frn>0.0) ? DIntVV[ic].Fr/Frn : DIntVV[ic].Fr;
        if (Frn>Int[id].Eta*Int[id].Mu*Fnn)
        {
            DIntVV[ic].Fr = Int[id].Eta*Int[id].Mu*Fnn*tan;
        }
        
        DIntVV[ic].F = DIntVV[ic].Fn + DIntVV[ic].Ft + Int[id].Gn*sqrt(sqrtdelta)*dotreal3(n,vrel)*n + Int[id].Gt*sqrt(sqrtdelta)*vt;

        real3 T1,T2,T, Tt;
        Tt = cross (x1,DIntVV[ic].F) + r1*cross(n,DIntVV[ic].Fr);
        real4 q;
        Conjugate (DPar[i1].Q,q);
        Rotation  (Tt,q,T);
        T1 = -1.0*T;
        Tt = cross (x2,DIntVV[ic].F) + r2*cross(n,DIntVV[ic].Fr);
        Conjugate (DPar[i2].Q,q);
        Rotation  (Tt,q,T);
        T2 =      T;

        atomicAdd(&CInt[id].Fnnet.x, DIntVV[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntVV[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntVV[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntVV[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntVV[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntVV[ic].Ft.z);

        atomicAdd(&DPar[i1].F.x,-DIntVV[ic].F .x);
        atomicAdd(&DPar[i1].F.y,-DIntVV[ic].F .y);
        atomicAdd(&DPar[i1].F.z,-DIntVV[ic].F .z);
        atomicAdd(&DPar[i2].F.x, DIntVV[ic].F .x);
        atomicAdd(&DPar[i2].F.y, DIntVV[ic].F .y);
        atomicAdd(&DPar[i2].F.z, DIntVV[ic].F .z);
        atomicAdd(& Par[i1].T.x,            T1.x);
        atomicAdd(& Par[i1].T.y,            T1.y);
        atomicAdd(& Par[i1].T.z,            T1.z);
        atomicAdd(& Par[i2].T.x,            T2.x);
        atomicAdd(& Par[i2].T.y,            T2.y);
        atomicAdd(& Par[i2].T.z,            T2.z);
    }
}

__global__ void CalcForceEE(size_t const * __restrict__ Edges, real3 const * __restrict__ Verts, InteractonCU const * __restrict__ Int, ComInteractonCU * __restrict__ CInt, DynInteractonCU * __restrict__ DIntEE,
        ParticleCU * __restrict__ Par, DynParticleCU * __restrict__ DPar, dem_aux const * __restrict__ demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].neeint) return;
    size_t id = DIntEE[ic].Idx;
    size_t i1 = CInt  [id].I1;
    size_t i2 = CInt  [id].I2;
    size_t f1 = DIntEE[ic].IF1;
    size_t f2 = DIntEE[ic].IF2;
    real  dm1 = DIntEE[ic].Dmax1;
    real  dm2 = DIntEE[ic].Dmax2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    
    DIntEE[ic].Fn = make_real3(0.0,0.0,0.0);
    
    real3 s;
    real3 Pert = make_real3(0.0,0.0,0.0);
    if (Int[id].BothFree) Pert = demaux[0].Per;
    if (!OverlapEE(Edges,Verts,f1,f2,dm1,dm2,Pert)) return;
    DistanceEE(Edges,Verts,f1,f2,xi,xf,s,Pert);
    real dist  = norm(s);
    real delta = r1 + r2 - dist;
    if (delta>0)
    {
        real3  n = s/dist;
        real   d = (r1*r1-r2*r2+dist*dist)/(2*dist);
        real3  x1c = xi+d*n;
        real3  x2c = xf-(dist-d)*n;
        real3 t1,t2,x1,x2;
        Rotation(DPar[i1].w,DPar[i1].Q,t1);
        Rotation(DPar[i2].w,DPar[i2].Q,t2);
        x1 = x1c - DPar [i1].x;
        x2 = x2c - DPar [i2].x;
        real3 vrel = (DPar[i1].v+cross(t1,x1))-(DPar[i2].v+cross(t2,x2));
        real3 vt   = vrel - dotreal3(n,vrel)*n;

        DIntEE[ic].Fn  = Int[id].Kn*delta*n;
        DIntEE[ic].Ft  = DIntEE[ic].Ft + (Int[id].Kt*demaux[0].dt)*vt;
        DIntEE[ic].Ft  = DIntEE[ic].Ft - dotreal3(DIntEE[ic].Ft,n)*n;

        real  Ftn = norm(DIntEE[ic].Ft);
        real  Fnn = norm(DIntEE[ic].Fn);
        real3 tan = (Ftn>0.0) ? DIntEE[ic].Ft/Ftn : DIntEE[ic].Ft;
        if (Ftn>Int[id].Mu*Fnn)
        {
            DIntEE[ic].Ft = Int[id].Mu*Fnn*tan;
        }

        DIntEE[ic].F = DIntEE[ic].Fn + DIntEE[ic].Ft + Int[id].Gn*dotreal3(n,vrel)*n + Int[id].Gt*vt;

        real3 T1,T2,T, Tt;
        Tt = cross (x1,DIntEE[ic].F);
        real4 q;
        Conjugate (DPar[i1].Q,q);
        Rotation  (Tt,q,T);
        T1 = -1.0*T;
        Tt = cross (x2,DIntEE[ic].F);
        Conjugate (DPar[i2].Q,q);
        Rotation  (Tt,q,T);
        T2 =      T;

        atomicAdd(&CInt[id].Fnnet.x, DIntEE[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntEE[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntEE[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntEE[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntEE[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntEE[ic].Ft.z);

        atomicAdd(&DPar[i1].F.x,-DIntEE[ic].F .x);
        atomicAdd(&DPar[i1].F.y,-DIntEE[ic].F .y);
        atomicAdd(&DPar[i1].F.z,-DIntEE[ic].F .z);
        atomicAdd(&DPar[i2].F.x, DIntEE[ic].F .x);
        atomicAdd(&DPar[i2].F.y, DIntEE[ic].F .y);
        atomicAdd(&DPar[i2].F.z, DIntEE[ic].F .z);
        atomicAdd(& Par[i1].T.x,            T1.x);
        atomicAdd(& Par[i1].T.y,            T1.y);
        atomicAdd(& Par[i1].T.z,            T1.z);
        atomicAdd(& Par[i2].T.x,            T2.x);
        atomicAdd(& Par[i2].T.y,            T2.y);
        atomicAdd(& Par[i2].T.z,            T2.z);
    }
}

__global__ void CalcForceVF(size_t const * __restrict__ Faces, size_t const * __restrict__ Facid, real3 const * __restrict__ Verts, InteractonCU const * __restrict__ Int, ComInteractonCU * __restrict__ CInt,
        DynInteractonCU * __restrict__ DIntVF, ParticleCU * __restrict__ Par, DynParticleCU * __restrict__ DPar, dem_aux const * __restrict__ demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nvfint) return;
    size_t id = DIntVF[ic].Idx;
    size_t i1 = CInt  [id].I1;
    size_t i2 = CInt  [id].I2;
    size_t f1 = DIntVF[ic].IF1;
    size_t f2 = DIntVF[ic].IF2;
    real  dm1 = DIntVF[ic].Dmax1;
    real  dm2 = DIntVF[ic].Dmax2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    
    DIntVF[ic].Fn = make_real3(0.0,0.0,0.0);

    real3 s;
    xi = Verts[f1];
    real3 Pert = make_real3(0.0,0.0,0.0);
    if (Int[id].BothFree) Pert = demaux[0].Per;
    if (!OverlapVF(Faces,Facid,Verts,xi,f2,dm1,dm2,Pert)) return;
    DistanceVF(Faces,Facid,Verts,xi,f2,xf,s,Pert);
    real dist  = norm(s);
    real delta = r1 + r2 - dist;
    if (delta>0)
    {
        real3  n = s/dist;
        real   d = (r1*r1-r2*r2+dist*dist)/(2*dist);
        real3  x1c = xi+d*n;
        real3  x2c = xf-(dist-d)*n;
        real3 t1,t2,x1,x2;
        Rotation(DPar[i1].w,DPar[i1].Q,t1);
        Rotation(DPar[i2].w,DPar[i2].Q,t2);
        x1 = x1c - DPar [i1].x;
        x2 = x2c - DPar [i2].x;
        real3 vrel = (DPar[i1].v+cross(t1,x1))-(DPar[i2].v+cross(t2,x2));
        real3 vt   = vrel - dotreal3(n,vrel)*n;

        DIntVF[ic].Fn  = Int[id].Kn*delta*n;
        DIntVF[ic].Ft  = DIntVF[ic].Ft + (Int[id].Kt*demaux[0].dt)*vt;
        DIntVF[ic].Ft  = DIntVF[ic].Ft - dotreal3(DIntVF[ic].Ft,n)*n;

        real  Ftn = norm(DIntVF[ic].Ft);
        real  Fnn = norm(DIntVF[ic].Fn);
        real3 tan = (Ftn>0.0) ? DIntVF[ic].Ft/Ftn : DIntVF[ic].Ft;
        if (Ftn>Int[id].Mu*Fnn)
        {
            DIntVF[ic].Ft = Int[id].Mu*Fnn*tan;
        }

        DIntVF[ic].F = DIntVF[ic].Fn + DIntVF[ic].Ft + Int[id].Gn*dotreal3(n,vrel)*n + Int[id].Gt*vt;

        real3 T1,T2,T, Tt;
        Tt = cross (x1,DIntVF[ic].F);
        real4 q;
        Conjugate (DPar[i1].Q,q);
        Rotation  (Tt,q,T);
        T1 = -1.0*T;
        Tt = cross (x2,DIntVF[ic].F);
        Conjugate (DPar[i2].Q,q);
        Rotation  (Tt,q,T);
        T2 =      T;

        atomicAdd(&CInt[id].Fnnet.x, DIntVF[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntVF[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntVF[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntVF[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntVF[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntVF[ic].Ft.z);

        atomicAdd(&DPar[i1].F.x,-DIntVF[ic].F .x);
        atomicAdd(&DPar[i1].F.y,-DIntVF[ic].F .y);
        atomicAdd(&DPar[i1].F.z,-DIntVF[ic].F .z);
        atomicAdd(&DPar[i2].F.x, DIntVF[ic].F .x);
        atomicAdd(&DPar[i2].F.y, DIntVF[ic].F .y);
        atomicAdd(&DPar[i2].F.z, DIntVF[ic].F .z);
        atomicAdd(& Par[i1].T.x,            T1.x);
        atomicAdd(& Par[i1].T.y,            T1.y);
        atomicAdd(& Par[i1].T.z,            T1.z);
        atomicAdd(& Par[i2].T.x,            T2.x);
        atomicAdd(& Par[i2].T.y,            T2.y);
        atomicAdd(& Par[i2].T.z,            T2.z);
    }
}

__global__ void CalcForceFV(size_t const * __restrict__ Faces, size_t const * __restrict__ Facid, real3 const * __restrict__ Verts, InteractonCU const * __restrict__ Int, ComInteractonCU * __restrict__ CInt,
        DynInteractonCU * __restrict__ DIntFV, ParticleCU * __restrict__ Par, DynParticleCU * __restrict__ DPar, dem_aux const * __restrict__ demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nfvint) return;
    size_t id = DIntFV[ic].Idx;
    size_t i1 = CInt  [id].I1;
    size_t i2 = CInt  [id].I2;
    size_t f1 = DIntFV[ic].IF1;
    size_t f2 = DIntFV[ic].IF2;
    real  dm1 = DIntFV[ic].Dmax1;
    real  dm2 = DIntFV[ic].Dmax2;
    real   r1 = Par  [i1].R;
    real   r2 = Par  [i2].R;
    real3  xi = DPar [i1].x;
    real3  xf = DPar [i2].x;
    
    DIntFV[ic].Fn = make_real3(0.0,0.0,0.0);

    real3 s;
    xf = Verts[f2];
    real3 Pert = make_real3(0.0,0.0,0.0);
    if (Int[id].BothFree) Pert = demaux[0].Per;
    if (!OverlapFV(Faces,Facid,Verts,f1,xf,dm1,dm2,Pert)) return;
    DistanceFV(Faces,Facid,Verts,f1,xf,xi,s,Pert);
    real dist  = norm(s);
    real delta = r1 + r2 - dist;
    if (delta>0)
    {
        real3  n = s/dist;
        real   d = (r1*r1-r2*r2+dist*dist)/(2*dist);
        real3  x1c = xi+d*n;
        real3  x2c = xf-(dist-d)*n;
        real3 t1,t2,x1,x2;
        Rotation(DPar[i1].w,DPar[i1].Q,t1);
        Rotation(DPar[i2].w,DPar[i2].Q,t2);
        x1 = x1c - DPar [i1].x;
        x2 = x2c - DPar [i2].x;
        real3 vrel = (DPar[i1].v+cross(t1,x1))-(DPar[i2].v+cross(t2,x2));
        real3 vt   = vrel - dotreal3(n,vrel)*n;

        DIntFV[ic].Fn  = Int[id].Kn*delta*n;
        DIntFV[ic].Ft  = DIntFV[ic].Ft + (Int[id].Kt*demaux[0].dt)*vt;
        DIntFV[ic].Ft  = DIntFV[ic].Ft - dotreal3(DIntFV[ic].Ft,n)*n;

        real  Ftn = norm(DIntFV[ic].Ft);
        real  Fnn = norm(DIntFV[ic].Fn);
        real3 tan = (Ftn>0.0) ? DIntFV[ic].Ft/Ftn : DIntFV[ic].Ft;
        if (Ftn>Int[id].Mu*Fnn)
        {
            DIntFV[ic].Ft = Int[id].Mu*Fnn*tan;
        }

        DIntFV[ic].F = DIntFV[ic].Fn + DIntFV[ic].Ft + Int[id].Gn*dotreal3(n,vrel)*n + Int[id].Gt*vt;
        
        real3 T1,T2,T, Tt;
        Tt = cross (x1,DIntFV[ic].F);
        real4 q;
        Conjugate (DPar[i1].Q,q);
        Rotation  (Tt,q,T);
        T1 = -1.0*T;
        Tt = cross (x2,DIntFV[ic].F);
        Conjugate (DPar[i2].Q,q);
        Rotation  (Tt,q,T);
        T2 =      T;

        atomicAdd(&CInt[id].Fnnet.x, DIntFV[ic].Fn.x);
        atomicAdd(&CInt[id].Fnnet.y, DIntFV[ic].Fn.y);
        atomicAdd(&CInt[id].Fnnet.z, DIntFV[ic].Fn.z);
        atomicAdd(&CInt[id].Ftnet.x, DIntFV[ic].Ft.x);
        atomicAdd(&CInt[id].Ftnet.y, DIntFV[ic].Ft.y);
        atomicAdd(&CInt[id].Ftnet.z, DIntFV[ic].Ft.z);

        atomicAdd(&DPar[i1].F.x,-DIntFV[ic].F .x);
        atomicAdd(&DPar[i1].F.y,-DIntFV[ic].F .y);
        atomicAdd(&DPar[i1].F.z,-DIntFV[ic].F .z);
        atomicAdd(&DPar[i2].F.x, DIntFV[ic].F .x);
        atomicAdd(&DPar[i2].F.y, DIntFV[ic].F .y);
        atomicAdd(&DPar[i2].F.z, DIntFV[ic].F .z);
        atomicAdd(& Par[i1].T.x,            T1.x);
        atomicAdd(& Par[i1].T.y,            T1.y);
        atomicAdd(& Par[i1].T.z,            T1.z);
        atomicAdd(& Par[i2].T.x,            T2.x);
        atomicAdd(& Par[i2].T.y,            T2.y);
        atomicAdd(& Par[i2].T.z,            T2.z);
    }
}

__global__ void Translate(real3 * Verts, ParticleCU const * Par, DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nparts) return;
    real3 Ft = DPar[ic].F;
    if (Par[ic].vxf) Ft.x = 0.0;
    if (Par[ic].vyf) Ft.y = 0.0;
    if (Par[ic].vzf) Ft.z = 0.0;
    //if (ic==0) printf("Fz = %g \n",DPar[ic].Flbm.z);
    real3 temp,xa;
    xa    = 2.0*DPar[ic].x - DPar[ic].xb + (demaux[0].dt*demaux[0].dt/Par[ic].m)*Ft;
    temp  = xa - DPar[ic].x;
    DPar[ic].v    = 0.5*(xa - DPar[ic].xb)/demaux[0].dt;
    DPar [ic].xb  = DPar[ic].x;
    DPar[ic].x    = xa;


    //if (isnan(norm(DPar[ic].x))) printf("ic %lu it %lu \n",ic,demaux[0].iter);

    for (size_t iv=Par[ic].Nvi;iv<Par[ic].Nvf;iv++)
    {
        Verts [iv] = Verts [iv] + temp;
    }
}

__global__ void Rotate(real3 * Verts, ParticleCU const * Par, DynParticleCU * DPar, dem_aux const * demaux, void * extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nparts) return;
    real q0,q1,q2,q3,wx,wy,wz;

    q0 = 0.5*DPar[ic].Q.w;
    q1 = 0.5*DPar[ic].Q.x;
    q2 = 0.5*DPar[ic].Q.y;
    q3 = 0.5*DPar[ic].Q.z;

    real3 Tt = Par[ic].T;

    if (Par[ic].wxf) Tt.x = 0.0;
    if (Par[ic].wyf) Tt.y = 0.0;
    if (Par[ic].wzf) Tt.z = 0.0;

    //if (isnan(norm(Tt)))
    //{
        //printf("ic %d \n",ic);
        //printf("it %d \n",demaux[0].iter);
    //}

    DPar[ic].wa.x=(Tt.x+(Par[ic].I.y-Par[ic].I.z)*DPar[ic].wb.y*DPar[ic].wb.z)/Par[ic].I.x;
    DPar[ic].wa.y=(Tt.y+(Par[ic].I.z-Par[ic].I.x)*DPar[ic].wb.x*DPar[ic].wb.z)/Par[ic].I.y;
    DPar[ic].wa.z=(Tt.z+(Par[ic].I.x-Par[ic].I.y)*DPar[ic].wb.y*DPar[ic].wb.x)/Par[ic].I.z;
    DPar[ic].w = DPar[ic].wb+0.5*demaux[0].dt*DPar[ic].wa;

    wx = DPar[ic].w.x;
    wy = DPar[ic].w.y;
    wz = DPar[ic].w.z;
    real4 dq,qm;
    dq.w = -(q1*wx+q2*wy+q3*wz);
    dq.x = q0*wx-q3*wy+q2*wz;
    dq.y = q3*wx+q0*wy-q1*wz;
    dq.z = -q2*wx+q1*wy+q0*wz;

    DPar[ic].wb  = DPar[ic].wb+demaux[0].dt*DPar[ic].wa;
    qm  = DPar[ic].Q+(0.5*demaux[0].dt)*dq;

    q0 = 0.5*qm.w;
    q1 = 0.5*qm.x;
    q2 = 0.5*qm.y;
    q3 = 0.5*qm.z;

    wx  = DPar[ic].wb.x;
    wy  = DPar[ic].wb.y;
    wz  = DPar[ic].wb.z;
    
    dq.w = -(q1*wx+q2*wy+q3*wz);
    dq.x = q0*wx-q3*wy+q2*wz;
    dq.y = q3*wx+q0*wy-q1*wz;
    dq.z = -q2*wx+q1*wy+q0*wz;

    real4 Qd = qm+0.5*demaux[0].dt*dq,temp;
    Conjugate(DPar[ic].Q,temp);

    // Fused: un-rotate by the old orientation and re-rotate by the new one in a
    // single pass over this particle's vertices.  This was two passes, each
    // reading and writing every vertex.  The arithmetic is kept in exactly the
    // original order -- in particular xt2 is recomputed as (y + x) - x rather
    // than reusing the pre-addition value, because (a+x)-x is not a in general
    // and reusing it would change the last bits.
    real4 Qnew = Qd/norm(Qd);
    for (size_t iv=Par[ic].Nvi;iv<Par[ic].Nvf;iv++)
    {
        real3 xt = Verts[iv] - DPar[ic].x;
        real3 y;
        Rotation(xt,temp,y);
        y = y + DPar[ic].x;
        real3 xt2 = y - DPar[ic].x;
        real3 z;
        Rotation(xt2,Qnew,z);
        Verts[iv] = z + DPar[ic].x;
    }

    DPar[ic].Q = Qnew;
}

__global__ void Reset (ParticleCU * Par, DynParticleCU * DPar, InteractonCU const * Int, ComInteractonCU * CInt, dem_aux const * demaux, void *
        extraparams)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic<demaux[0].nparts)
    {
        DPar[ic].F    = Par[ic].Ff;
        DPar[ic].Flbm = Par[ic].Flbmf;
        Par [ic].T    = Par[ic].Tf;
    }
    else if (ic<demaux[0].nparts+demaux[0].ncoint)
    {
        size_t id = ic-demaux[0].nparts;
        CInt[id].Fnnet = Int[id].Fnf;
        CInt[id].Ftnet = Int[id].Ftf;
    }
    else return;
}

__global__ void MaxD(real3 const * Verts, real3 const * Vertso, real * maxd, dem_aux * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    // Reduction note: this used to write one value per VERTEX (nverts entries),
    // which the host then scanned with thrust::max_element every step.  The max
    // of a set is exact and independent of reduction order, so reducing inside
    // the block and writing one value per BLOCK gives the identical answer while
    // cutting the traffic by ~blockDim.  The comparison keeps thrust's
    // "update only when strictly less" rule so ties (and NaNs) behave as before.
    real m = 0.0;
    if (ic<demaux[0].nverts)
    {
        if (ic==0) 
        {
            demaux[0].Time += demaux[0].dt;
            demaux[0].iter++;
        }
        m = norm(Vertso[ic]-Verts[ic]);
    }
    #pragma unroll
    for (int off=16;off>0;off>>=1)
    {
        real o = __shfl_down_sync(0xffffffffu,m,off);
        if (m<o) m = o;
    }
    __shared__ real smax[32];
    size_t lane = threadIdx.x & 31u;
    size_t wid  = threadIdx.x >> 5u;
    if (lane==0) smax[wid] = m;
    __syncthreads();
    if (wid==0)
    {
        size_t nw = (blockDim.x+31u)>>5u;
        m = (lane<nw) ? smax[lane] : 0.0;
        #pragma unroll
        for (int off=16;off>0;off>>=1)
        {
            real o = __shfl_down_sync(0xffffffffu,m,off);
            if (m<o) m = o;
        }
        if (lane==0) maxd[blockIdx.x] = m;
    }
    //if (maxd[ic]>0.0)
    //{
        //printf("ic %d \n",ic);
        //printf("md %g \n",maxd[ic]);
    //}
}

__global__ void ResetMaxD(real3 * Verts, real3 * Vertso, real * maxd, ParticleCU const * Par, DynParticleCU * DPar, dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=demaux[0].nparts) return;

    bool isfree = ((!Par[ic].vxf&&!Par[ic].vyf&&!Par[ic].vzf&&!Par[ic].wxf&&!Par[ic].wyf&&!Par[ic].wzf)||Par[ic].FixFree);

    real3 dis = make_real3(0.0,0.0,0.0);

    if (isfree)
    {
        if (demaux[0].px)
        {
            if (DPar[ic].x.x< demaux[0].Xmin) dis.x = demaux[0].Xmax - demaux[0].Xmin;
            if (DPar[ic].x.x>=demaux[0].Xmax) dis.x = demaux[0].Xmin - demaux[0].Xmax;
        }
        if (demaux[0].py)
        {
            if (DPar[ic].x.y< demaux[0].Ymin) dis.y = demaux[0].Ymax - demaux[0].Ymin;
            if (DPar[ic].x.y>=demaux[0].Ymax) dis.y = demaux[0].Ymin - demaux[0].Ymax;
        }
        if (demaux[0].pz)
        {
            if (DPar[ic].x.z< demaux[0].Zmin) dis.z = demaux[0].Zmax - demaux[0].Zmin;
            if (DPar[ic].x.z>=demaux[0].Zmax) dis.z = demaux[0].Zmin - demaux[0].Zmax;
        }
    }

    DPar[ic].x  = DPar[ic].x  + dis;
    DPar[ic].xb = DPar[ic].xb + dis;

    for (size_t iv=Par[ic].Nvi;iv<Par[ic].Nvf;iv++)
    {
        Verts [iv] = Verts [iv] + dis;
        Vertso[iv] = Verts [iv];
        maxd  [iv] = 0.0;
    }

}

}
#endif //MECHSYS_DEM_CUH
