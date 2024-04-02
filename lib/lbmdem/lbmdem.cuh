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

/////////////////////////////LBMDEM CUDA implementation////////////////////

#ifndef MECHSYS_LBMDEM_CUH
#define MECHSYS_LBMDEM_CUH

//Mechsys
#include <mechsys/linalg/matvec.h>
#include <mechsys/dem/dem.cuh>
#include <mechsys/flbm/lbm.cuh>

namespace LBMDEM
{

struct lbmdem_aux
{
    size_t nvc;  ///< number of vertex cell pairs
    size_t nfc;  ///< number of face cell pairs
};

struct ParCellPairCU
{
    size_t Ic; ///< cell index
    size_t Ip; ///< particle index
    size_t Nfi;///< initial feature
    size_t Nff;///< final feature
};

__device__ real cudaSphereCube(real3 & Xs, real3 & Xc, real R, real dx)
{
    real3 P[8];
    P[0] = Xc - 0.5*dx*make_real3(1.0,0.0,0.0) - 0.5*dx*make_real3(0.0,1.0,0.0) + 0.5*dx*make_real3(0.0,0.0,1.0); 
    P[1] = Xc + 0.5*dx*make_real3(1.0,0.0,0.0) - 0.5*dx*make_real3(0.0,1.0,0.0) + 0.5*dx*make_real3(0.0,0.0,1.0);
    P[2] = Xc + 0.5*dx*make_real3(1.0,0.0,0.0) + 0.5*dx*make_real3(0.0,1.0,0.0) + 0.5*dx*make_real3(0.0,0.0,1.0);
    P[3] = Xc - 0.5*dx*make_real3(1.0,0.0,0.0) + 0.5*dx*make_real3(0.0,1.0,0.0) + 0.5*dx*make_real3(0.0,0.0,1.0);
    P[4] = Xc - 0.5*dx*make_real3(1.0,0.0,0.0) - 0.5*dx*make_real3(0.0,1.0,0.0) - 0.5*dx*make_real3(0.0,0.0,1.0); 
    P[5] = Xc + 0.5*dx*make_real3(1.0,0.0,0.0) - 0.5*dx*make_real3(0.0,1.0,0.0) - 0.5*dx*make_real3(0.0,0.0,1.0);
    P[6] = Xc + 0.5*dx*make_real3(1.0,0.0,0.0) + 0.5*dx*make_real3(0.0,1.0,0.0) - 0.5*dx*make_real3(0.0,0.0,1.0);
    P[7] = Xc - 0.5*dx*make_real3(1.0,0.0,0.0) + 0.5*dx*make_real3(0.0,1.0,0.0) - 0.5*dx*make_real3(0.0,0.0,1.0);
    
    real dmin = 2*R;
    real dmax = 0.0;
    for (size_t j=0;j<8;j++)
    {
        real dist = norm(P[j] - Xs);
        if (dmin>dist) dmin = dist;
        if (dmax<dist) dmax = dist;
    }
    if (dmin > R + dx) return 0.0;
    
    if (dmax < R)
    {
        return 12.0*dx;
    }
    
    real len = 0.0;
    for (size_t j=0;j<4;j++)
    {
        real3 D;
        real a; 
        real b; 
        real c; 
        D = P[(j+1)%4] - P[j];
        a = dotreal3(D,D);
        b = 2*dotreal3(P[j]-Xs,D);
        c = dotreal3(P[j]-Xs,P[j]-Xs) - R*R;
        if (b*b-4*a*c>0.0)
        {
            real ta = (-b - sqrt(b*b-4*a*c))/(2*a);
            real tb = (-b + sqrt(b*b-4*a*c))/(2*a);
            if (ta>1.0&&tb>1.0) continue;
            if (ta<0.0&&tb<0.0) continue;
            if (ta<0.0) ta = 0.0;
            if (tb>1.0) tb = 1.0;
            len += norm((tb-ta)*D);
        }
        D = P[(j+1)%4 + 4] - P[j + 4];
        a = dotreal3(D,D);
        b = 2*dotreal3(P[j + 4]-Xs,D);
        c = dotreal3(P[j + 4]-Xs,P[j + 4]-Xs) - R*R;
        if (b*b-4*a*c>0.0)
        {
            real ta = (-b - sqrt(b*b-4*a*c))/(2*a);
            real tb = (-b + sqrt(b*b-4*a*c))/(2*a);
            if (ta>1.0&&tb>1.0) continue;
            if (ta<0.0&&tb<0.0) continue;
            if (ta<0.0) ta = 0.0;
            if (tb>1.0) tb = 1.0;
            len += norm((tb-ta)*D);
        }
        D = P[j+4] - P[j];
        a = dotreal3(D,D);
        b = 2*dotreal3(P[j]-Xs,D);
        c = dotreal3(P[j]-Xs,P[j]-Xs) - R*R;
        if (b*b-4*a*c>0.0)
        {
            real ta = (-b - sqrt(b*b-4*a*c))/(2*a);
            real tb = (-b + sqrt(b*b-4*a*c))/(2*a);
            if (ta>1.0&&tb>1.0) continue;
            if (ta<0.0&&tb<0.0) continue;
            if (ta<0.0) ta = 0.0;
            if (tb>1.0) tb = 1.0;
            len += norm((tb-ta)*D);
        }
    }
    return len;
}

__global__ void cudaReset(bool const * IsSolid, real * Gammaf, real * Gamma, real * Omeis,FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;

    if (IsSolid[ic])  Gamma[ic] = 1.0;
    else              Gamma[ic] = Gammaf[ic];
    //for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    //{
        //Omeis[ic*lbmaux[0].Nneigh + k] = 0.0;
    //}
}

__global__ void cudaImprintLatticeVC(size_t const * PaCeV, DEM::ParticleCU * Par, DEM::DynParticleCU * DPar, real const * Rho, real * Gamma, real * Omeis, real const *
        F, DEM::dem_aux const * demaux, FLBM::lbm_aux const * lbmaux, lbmdem_aux const * lbmdemaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmdemaux[0].nvc) return;

    size_t ice = PaCeV[2*ic];
    size_t icx =  ice%lbmaux[0].Nx;
    size_t icy = (ice/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ice/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;

    size_t ip  = PaCeV[2*ic+1];

    real3  C = lbmaux[0].dx*make_real3(real(icx),real(icy),real(icz));
    real3  Xs,B;
    real   len = 12.0*lbmaux[0].dx;
    real3  Pert = demaux[0].Per;
    bool isfree = ((!Par[ip].vxf&&!Par[ip].vyf&&!Par[ip].vzf&&!Par[ip].wxf&&!Par[ip].wyf&&!Par[ip].wzf)||Par[ip].FixFree);
    if (!isfree) Pert = make_real3(0.0,0.0,0.0);
    DEM::BranchVec(DPar[ip].x,C,B,Pert);
    if (norm(B)>Par[ip].Dmax) return;
    Xs  = C-B;
    len = cudaSphereCube(Xs,C,Par[ip].Dmax,lbmaux[0].dx);
    if (fabs(len)<1.0e-12) return;
    real Tau = lbmaux[0].Tau[0];
    real gamma  = len/(12.0*lbmaux[0].dx);
    if (gamma<Gamma[ice]) return;
    Gamma[ice] = gamma;
    real3 tmp;
    Rotation(DPar[ip].w,DPar[ip].Q,tmp);
    real3 VelP   = DPar[ip].v + cross(tmp,B);
    real rho = Rho[ice];
    real Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
    //real Bn  = gamma;
    size_t ncells = lbmaux[0].Nneigh;
    real3 Flbm = make_real3(0.0,0.0,0.0);
    for (size_t k=0;k<ncells;k++)
    {
        real Fvpp     = FLBM::FeqFluid(lbmaux[0].Op[k],rho,VelP,lbmaux);
        real Fvp      = FLBM::FeqFluid(k              ,rho,VelP,lbmaux);
        real Omega    = F[ice*lbmaux[0].Nneigh + lbmaux[0].Op[k]] - Fvpp - (F[ice*lbmaux[0].Nneigh + k] - Fvp);
        Omeis[ice*lbmaux[0].Nneigh + k] = Omega;
        Flbm = Flbm - Bn*Omega*lbmaux[0].Cs*lbmaux[0].Cs*lbmaux[0].dx*lbmaux[0].dx*lbmaux[0].C[k];
    }
    real3 Tlbm,Tt;
    Tt =           cross(B,Flbm);
    real4 q;
    Conjugate    (DPar[ip].Q,q);
    Rotation     (Tt,q,Tlbm);

    atomicAdd(&DPar[ip].F   .x,Flbm.x);
    atomicAdd(&DPar[ip].F   .y,Flbm.y);
    atomicAdd(&DPar[ip].F   .z,Flbm.z);
    atomicAdd(&DPar[ip].Flbm.x,Flbm.x);
    atomicAdd(&DPar[ip].Flbm.y,Flbm.y);
    atomicAdd(&DPar[ip].Flbm.z,Flbm.z);
    atomicAdd(& Par[ip].T.x   ,Tlbm.x);
    atomicAdd(& Par[ip].T.y   ,Tlbm.y);
    atomicAdd(& Par[ip].T.z   ,Tlbm.z);
}

__global__ void cudaImprintLatticeFC(ParCellPairCU const * PaCe, size_t const * PaCeF, size_t const * Faces, size_t const * Facid, real3 const * Verts, DEM::ParticleCU * Par, DEM::DynParticleCU * DPar
        , real const * Rho, real * Gamma, real * Omeis, real const * F, DEM::dem_aux const * demaux, FLBM::lbm_aux const * lbmaux, lbmdem_aux const * lbmdemaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmdemaux[0].nfc) return;

    size_t ice = PaCe[ic].Ic;
    size_t icx =  ice%lbmaux[0].Nx;
    size_t icy = (ice/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ice/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;

    size_t ip  = PaCe[ic].Ip;

    real3  C = lbmaux[0].dx*make_real3(real(icx),real(icy),real(icz));
    real3  Xs,B,S;
    real   len = 12.0*lbmaux[0].dx,minl=Par[ip].Dmax;
    real3  Pert = demaux[0].Per;
    bool isfree = ((!Par[ip].vxf&&!Par[ip].vyf&&!Par[ip].vzf&&!Par[ip].wxf&&!Par[ip].wyf&&!Par[ip].wzf)||Par[ip].FixFree);
    if (!isfree) Pert = make_real3(0.0,0.0,0.0);
    DEM::BranchVec(DPar[ip].x,C,B,Pert);
    if (norm(B)>Par[ip].Dmax) return;
    size_t igeo = PaCe[ic].Nff - PaCe[ic].Nfi;
    if (igeo>0)
    {
        real3 xi;
        size_t f1 = PaCeF[PaCe[ic].Nfi];
        DEM::DistanceFV(Faces,Facid,Verts,f1,C,xi,S,Pert);
        minl = norm(S);
        Xs   = xi;
        real3  dL0 = Verts[Facid[Faces[2*f1]+1]]-Verts[Facid[Faces[2*f1]  ]];
        real3  dL1 = Verts[Facid[Faces[2*f1]+2]]-Verts[Facid[Faces[2*f1]+1]];
        real3  Nor = cross(dL0, dL1);
        Nor        = Nor/norm(Nor);
        for (size_t iff=PaCe[ic].Nfi+1;iff<PaCe[ic].Nff;iff++)
        {
            f1 = PaCeF[iff];
            real3 St;
            DEM::DistanceFV(Faces,Facid,Verts,f1,C,xi,St,Pert);
            if (norm(St) < minl)
            {
                S  = St;
                minl = norm(S);
                Xs = xi;
                dL0 = Verts[Facid[Faces[2*f1]+1]]-Verts[Facid[Faces[2*f1]  ]];
                dL1 = Verts[Facid[Faces[2*f1]+2]]-Verts[Facid[Faces[2*f1]+1]];
                Nor = cross(dL0, dL1);
                Nor = Nor/norm(Nor);
            }
        }
        real dotpro = dotreal3(S,Nor);
        if (dotpro>0.0||fabs(dotpro)<0.95*minl||(Par[ip].Nff-Par[ip].Nfi<4)||!Par[ip].Closed)
        {
            Xs = C-S;
            len = cudaSphereCube(Xs,C,Par[ip].R,lbmaux[0].dx);
        }
    }

    if (fabs(len)<1.0e-12) return;
    real Tau = lbmaux[0].Tau[0];
    real gamma  = len/(12.0*lbmaux[0].dx);
    if (gamma<Gamma[ice]) return;
    Gamma[ice] = gamma;
    real3 tmp;
    Rotation(DPar[ip].w,DPar[ip].Q,tmp);
    real3 VelP   = DPar[ip].v + cross(tmp,B);
    real rho = Rho[ice];
    real Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
    //real Bn  = gamma;
    size_t ncells = lbmaux[0].Nneigh;
    real3 Flbm = make_real3(0.0,0.0,0.0);
    for (size_t k=0;k<ncells;k++)
    {
        real Fvpp     = FLBM::FeqFluid(lbmaux[0].Op[k],rho,VelP,lbmaux);
        real Fvp      = FLBM::FeqFluid(k              ,rho,VelP,lbmaux);
        real Omega    = F[ice*lbmaux[0].Nneigh + lbmaux[0].Op[k]] - Fvpp - (F[ice*lbmaux[0].Nneigh + k] - Fvp);
        Omeis[ice*lbmaux[0].Nneigh + k] = Omega;
        Flbm = Flbm - Bn*Omega*lbmaux[0].Cs*lbmaux[0].Cs*lbmaux[0].dx*lbmaux[0].dx*lbmaux[0].C[k];
    }
    real3 Tlbm,Tt;
    Tt =           cross(B,Flbm);
    real4 q;
    Conjugate    (DPar[ip].Q,q);
    Rotation     (Tt,q,Tlbm);

    atomicAdd(&DPar[ip].F   .x,Flbm.x);
    atomicAdd(&DPar[ip].F   .y,Flbm.y);
    atomicAdd(&DPar[ip].F   .z,Flbm.z);
    atomicAdd(&DPar[ip].Flbm.x,Flbm.x);
    atomicAdd(&DPar[ip].Flbm.y,Flbm.y);
    atomicAdd(&DPar[ip].Flbm.z,Flbm.z);
    atomicAdd(& Par[ip].T.x   ,Tlbm.x);
    atomicAdd(& Par[ip].T.y   ,Tlbm.y);
    atomicAdd(& Par[ip].T.z   ,Tlbm.z);
}

__global__ void cudaStream2(bool const * IsSolid, real * Gamma, real * Omeis, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;
    if (ic==0)
    {
        lbmaux[0].Time += lbmaux[0].dt;
    }
    //for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    //{
        //F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k];
    //}
    Rho   [ic] = 0.0;
    Vel   [ic] = make_real3(0.0,0.0,0.0);
    BForce[ic] = make_real3(0.0,0.0,0.0);
    //Gamma[ic]  = (real) IsSolid[ic];
    if (!IsSolid[ic])
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ic] += F[ic*lbmaux[0].Nneigh + k];
            Vel[ic] = Vel[ic] + F[ic*lbmaux[0].Nneigh + k]*lbmaux[0].C[k];
            Omeis[ic*lbmaux[0].Nneigh + k] = 0.0;
            //if (ic==0) printf("k: %lu %f %f %f %f \n",k,Rho[ic],Vel[ic].x,Vel[ic].y,Vel[ic].z);
        }
        Vel[ic] = lbmaux[0].Cs/Rho[ic]*Vel[ic];
    }
}
}
#endif //MECHSYS_LBMDEM_CUH
