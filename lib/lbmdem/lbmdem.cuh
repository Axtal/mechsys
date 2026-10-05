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

// The distribution buffers here share FLBM's direction-major layout (see the
// comment on FLBM::FIDX in flbm/lbm.cuh); bring the index helper into scope so
// the kernels below read the same way they did when that layout was implicit.
using FLBM::FIDX;

struct lbmdem_aux
{
    size_t nvc;  ///< number of vertex cell pairs
    size_t nfc;  ///< number of face cell pairs
    real  Fconv; ///< Force conversion factor
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

__global__ void cudaReset(bool const * __restrict__ IsSolid, real const * __restrict__ Gammaf, real * __restrict__ Gamma, real * __restrict__ Omeis, FLBM::lbm_aux * __restrict__ lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;
    if (IsSolid[ic])  Gamma[ic] = 1.0;
    else              Gamma[ic] = Gammaf[ic];
}

/*
//#ifdef USE_IBB
__global__ void cudaCheckOutsideVC(size_t const * PaCeV, DEM::ParticleCU * Par, DEM::DynParticleCU * DPar, real * Gamma, 
        int * Inside, DEM::dem_aux const * demaux, FLBM::lbm_aux const * lbmaux, lbmdem_aux const * lbmdemaux)
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
    real3  Pert = demaux[0].Per;
    bool isfree = ((!Par[ip].vxf&&!Par[ip].vyf&&!Par[ip].vzf&&!Par[ip].wxf&&!Par[ip].wyf&&!Par[ip].wzf)||Par[ip].FixFree);
    if (!isfree) Pert = make_real3(0.0,0.0,0.0);
    DEM::BranchVec(DPar[ip].x,C,B,Pert);
    if (norm(B)>Par[ip].Dmax&&Inside[ice]==ip)
    {
        Inside[ice] = -ip-2;
        //if (icx<100&&icy==100&&icz==100) printf("ice %lu iter %lu \n",ice,lbmaux[0].iter);
        //printf("ice %lu iter %lu \n",ice,lbmaux[0].iter);
    }
}

__global__ void cudaCheckInsideVC(size_t const * PaCeV, DEM::ParticleCU * Par, DEM::DynParticleCU * DPar, real * Gamma, 
        int * Inside, DEM::dem_aux const * demaux, FLBM::lbm_aux const * lbmaux, lbmdem_aux const * lbmdemaux)
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
    real3  Pert = demaux[0].Per;
    bool isfree = ((!Par[ip].vxf&&!Par[ip].vyf&&!Par[ip].vzf&&!Par[ip].wxf&&!Par[ip].wyf&&!Par[ip].wzf)||Par[ip].FixFree);
    if (!isfree) Pert = make_real3(0.0,0.0,0.0);
    DEM::BranchVec(DPar[ip].x,C,B,Pert);
    if (norm(B)<Par[ip].Dmax)
    {
        Gamma [ice] = 1.0;
        Inside[ice] = ip;
        //printf("ice %lu \n",ice);
    }
}

__global__ void cudaRefill(DEM::ParticleCU * Par, DEM::DynParticleCU * DPar,real * F, real * Rho, real3 * Vel, int * Inside, DEM::dem_aux const * demaux, FLBM::lbm_aux const * lbmaux, lbmdem_aux const * lbmdemaux)
{
    size_t ice = threadIdx.x + blockIdx.x * blockDim.x;
    if (ice>=lbmaux[0].Nl*lbmaux[0].Ncells) return;

    const size_t Ncells = lbmaux[0].Ncells;
    const size_t Nn     = lbmaux[0].Nneigh;
    const size_t icl    = ice/Ncells;
    const size_t ice0   = ice - icl*Ncells;

    size_t icx =  ice%lbmaux[0].Nx;
    size_t icy = (ice/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ice/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;
    

    if (Inside[ice]<=-2)
    {
        for (size_t k = 0; k < lbmaux[0].Nneigh ; k++)
        {
            F[FIDX(Ncells,Nn,icl,k,ice0)] = 0.0;
        }
        size_t naem = 0;
        for (size_t k = 1; k < lbmaux[0].Nneigh ; k++)
        {
            size_t inx = (size_t)((int)icx +   (int)lbmaux[0].C[k ].x + (int)lbmaux[0].Nx)%lbmaux[0].Nx;
            size_t iny = (size_t)((int)icy +   (int)lbmaux[0].C[k ].y + (int)lbmaux[0].Ny)%lbmaux[0].Ny;
            size_t inz = (size_t)((int)icz +   (int)lbmaux[0].C[k ].z + (int)lbmaux[0].Nz)%lbmaux[0].Nz;
            size_t imx = (size_t)((int)icx + 2*(int)lbmaux[0].C[k ].x + (int)lbmaux[0].Nx)%lbmaux[0].Nx;
            size_t imy = (size_t)((int)icy + 2*(int)lbmaux[0].C[k ].y + (int)lbmaux[0].Ny)%lbmaux[0].Ny;
            size_t imz = (size_t)((int)icz + 2*(int)lbmaux[0].C[k ].z + (int)lbmaux[0].Nz)%lbmaux[0].Nz;
            size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;
            size_t im  = imx + imy*lbmaux[0].Nx + imz*lbmaux[0].Nx*lbmaux[0].Ny;

            if (Inside[in]>=0||Inside[in]<=-2) continue;
            if (Inside[im]>=0||Inside[im]<=-2) im = in;
            real rhon = 0.0;
            real rhom = 0.0;
            for (size_t kt = 0; kt < lbmaux[0].Nneigh ; kt++)
            {
                F[FIDX(Ncells,Nn,icl,kt,ice0)] += 2.0*F[FIDX(Ncells,Nn,0,kt,in)] - F[FIDX(Ncells,Nn,0,kt,im)];
                rhon += F[FIDX(Ncells,Nn,0,kt,in)];
                rhom += F[FIDX(Ncells,Nn,0,kt,im)];
            }
            naem++;
            //if (ice==3740896&&lbmaux[0].iter==1)
            //{
                //printf("Rho %g %g in %lu im %lu k %lu \n",rhon,rhom,in,im,k);
            //}
        }
        if (naem==0)
        {
            int ip = -2-Inside[ice];
            real3  C = lbmaux[0].dx*make_real3(real(icx),real(icy),real(icz));
            real3  B;
            real3  Pert = demaux[0].Per;
            bool isfree = ((!Par[ip].vxf&&!Par[ip].vyf&&!Par[ip].vzf&&!Par[ip].wxf&&!Par[ip].wyf&&!Par[ip].wzf)||Par[ip].FixFree);
            if (!isfree) Pert = make_real3(0.0,0.0,0.0);
            DEM::BranchVec(DPar[ip].x,C,B,Pert);
            real rho = Rho[ice];
            real3 tmp;
            Rotation(DPar[ip].w,DPar[ip].Q,tmp);
            real3 VelP   = DPar[ip].v + cross(tmp,B);
            for (size_t k = 1; k < lbmaux[0].Nneigh ; k++)
            {
                F[FIDX(Ncells,Nn,icl,k,ice0)] = FeqFluid(k,rho,VelP,lbmaux);
            }
            naem = 1;
        }

        //if (ice==3740896&&lbmaux[0].iter==1)
        //{
            //printf("naem %lu \n",naem);
        //}

        Rho[ice] = 0.0;
        Vel[ice] = make_real3(0.0,0.0,0.0);
        for (size_t k = 0; k < lbmaux[0].Nneigh ; k++)
        {
            F[FIDX(Ncells,Nn,icl,k,ice0)] = fabs(F[FIDX(Ncells,Nn,icl,k,ice0)])/naem;
            Rho[ice] += F[FIDX(Ncells,Nn,icl,k,ice0)];
            Vel[ice] = Vel[ice] + F[FIDX(Ncells,Nn,icl,k,ice0)]*lbmaux[0].C[k];
        }
        Vel[ice] = lbmaux[0].Cs/Rho[ice]*Vel[ice];
        //if (ice==3740896&&lbmaux[0].iter==1)
        //{
            //printf("Rho %g \n",Rho[ice]);
        //}

    }
}

__global__ void cudaImprintLatticeVC(size_t const * PaCeV, DEM::ParticleCU * Par, DEM::DynParticleCU * DPar, real const * Rho, real * Gamma, real * Omeis, real *
        F, int * Inside, DEM::dem_aux const * demaux, FLBM::lbm_aux const * lbmaux, lbmdem_aux const * lbmdemaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmdemaux[0].nvc) return;

    size_t ice = PaCeV[2*ic];
    size_t icx =  ice%lbmaux[0].Nx;
    size_t icy = (ice/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ice/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;
    const size_t Ncells_ = lbmaux[0].Ncells;
    const size_t Nn_     = lbmaux[0].Nneigh;
    const size_t icl_    = ice/Ncells_;
    const size_t ice0_   = ice - icl_*Ncells_;

    size_t ip  = PaCeV[2*ic+1];

    real3  C = lbmaux[0].dx*make_real3(real(icx),real(icy),real(icz));
    real3  Xs,B;
    real3  Pert = demaux[0].Per;
    bool isfree = ((!Par[ip].vxf&&!Par[ip].vyf&&!Par[ip].vzf&&!Par[ip].wxf&&!Par[ip].wyf&&!Par[ip].wzf)||Par[ip].FixFree);
    if (!isfree) Pert = make_real3(0.0,0.0,0.0);
    DEM::BranchVec(DPar[ip].x,C,B,Pert);
    Xs  = C-B;
    if (norm(B)>Par[ip].Dmax+2.0*lbmaux[0].dx||Inside[ice]==ip) return;
    if (Inside[ice]<=-2) Inside[ice] = -1;
    //if (Inside[ice]==ip)
    //{
        //real3 tmp;
        //Rotation(DPar[ip].w,DPar[ip].Q,tmp);
        //real rho = Rho[ice];
        //real3 VelP   = DPar[ip].v + cross(tmp,B);
        //for (size_t k = 1; k < lbmaux[0].Nneigh ; k++)
        //{
            //real Fvpp    = FLBM::FeqFluid(lbmaux[0].Op[k],rho,VelP,lbmaux);
            //real Fvp     = FLBM::FeqFluid(k              ,rho,VelP,lbmaux);
            //real Omega   = F[FIDX(Ncells_,Nn_,icl_,lbmaux[0].Op[k],ice0_)] - Fvpp - (F[FIDX(Ncells_,Nn_,icl_,k,ice0_)] - Fvp);
            //Omeis[FIDX(Ncells_,Nn_,icl_,k,ice0_)] = Omega;
        //}
        //return;
    //}
   
    real ld = lbmaux[0].dx;
    real Cs = lbmaux[0].Cs;

    real3 Flbm = make_real3(0.0,0.0,0.0);
    real rho = Rho[ice];
    real3 tmp;
    Rotation(DPar[ip].w,DPar[ip].Q,tmp);
    for (size_t k = 1; k < lbmaux[0].Nneigh ; k++)
    {
        size_t inx = (size_t)((int)icx + (int)lbmaux[0].C[k ].x + (int)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((int)icy + (int)lbmaux[0].C[k ].y + (int)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((int)icz + (int)lbmaux[0].C[k ].z + (int)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t ine = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;
        if(Inside[ine]!=ip) continue;
        real3   Cn = ld*make_real3(real(inx),real(iny),real(inz));
        real3   dC;
        //DEM::BranchVec(DPar[ip].x,Cn,dC,Pert);

        //Solving the ray intersecting sphere problem
        dC     = ld*lbmaux[0].C[k];
        real a =     dotreal3(dC,dC);
        real b = 2.0*dotreal3(dC,B );
        real c =     dotreal3(B ,B ) - Par[ip].Dmax*Par[ip].Dmax;

        real r = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);

        if (b*b<4.0*a*c) continue;
        //real r1 = (-b - sqrt(b*b - 4.0*a*c))/(2.0*a);
        //real r2 = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a);

        //real r;

        //if      (r1>0.0&&r1<=1.0) r = r1;
        //else if (r2>0.0&&r2<=1.0) r = r2;
        //else    continue;
        
        size_t ko  = lbmaux[0].Op[k];
        size_t ifx = (size_t)((int)icx + (int)lbmaux[0].C[ko].x + (int)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t ify = (size_t)((int)icy + (int)lbmaux[0].C[ko].y + (int)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t ifz = (size_t)((int)icz + (int)lbmaux[0].C[ko].z + (int)lbmaux[0].Nz)%lbmaux[0].Nz;

        size_t ife = ifx + ify*lbmaux[0].Nx + ifz*lbmaux[0].Nx*lbmaux[0].Ny;
        real3 Xw     = C + r*(dC) - Xs;
        real3 VelP   = DPar[ip].v + cross(tmp,Xw);
        
        //MPM IBB
        F[FIDX(Ncells_,Nn_,icl_,ko,ice0_)] = (r*F[FIDX(Ncells_,Nn_,0,ko,ife)] + (1.0-r)*F[FIDX(Ncells_,Nn_,icl_,k,ice0_)] 
                + r*F[FIDX(Ncells_,Nn_,0,k,ine)] + 6.0*rho*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[ko],VelP)/Cs)/(1.0+r);


        //F[FIDX(Ncells_,Nn_,0,k,ine)]  = F[FIDX(Ncells_,Nn_,icl_,k,ice0_)];

        Flbm = Flbm + ld*ld*Cs*(F[FIDX(Ncells_,Nn_,0,k,ine)]*(Cs*lbmaux[0].C[k]-VelP) - F[FIDX(Ncells_,Nn_,icl_,ko,ice0_)]*(Cs*lbmaux[0].C[ko]-VelP));
        //Flbm = Flbm + ld*ld*Cs*(F[FIDX(Ncells_,Nn_,0,k,ine)]*(Cs*lbmaux[0].C[k]) - F[FIDX(Ncells_,Nn_,icl_,ko,ice0_)]*(Cs*lbmaux[0].C[ko]));
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
//#else
*/
//#endif
//Functions for the smoothing function of LBM DEM interaction
typedef real (*FuncBn_ptr)(real,real);
 
__device__ __inline__ real BnSmooth(real gamma, real tau)
{
    return (gamma*(tau-0.5))/((1.0-gamma)+(tau-0.5));
}

__device__ __inline__ real BnLadd(real gamma, real tau)
{
    return floor(gamma);
}

__device__ FuncBn_ptr d_fBnSmooth = BnSmooth;

__device__ FuncBn_ptr d_fBnLadd   = BnLadd;

// ---------------------------------------------------------------------------
// Fused streaming + macroscopic update.
//
// The original per-step sequence was
//     FLBM::cudaStream1(...)   // scatter: Ftemp[ic + C[k]][k] = F[ic][k]
//     swap(F,Ftemp)
//     LBMDEM::cudaStream2(...) // Rho, Vel, reset BForce/Omeis
// which streams the whole distribution array in two extra full passes.
//
// Streaming is a pure lattice translation, so the scatter done by cudaStream1
// is exactly equivalent to the gather done here:
//     F[ic][k] = Fpost[ic - C[k]][k]
// Doing the gather allows the streaming, the density/velocity accumulation and
// the per-step resets to share a single pass over memory, removing one full
// read plus one full write of the distribution array per time step.
//
// The Omeis reset is additionally restricted to cells with Gamma>0. Only
// ImprintLattice* ever writes a non-zero Omeis, and it always leaves Gamma>0
// on those same cells, so "Omeis==0 wherever Gamma==0" is preserved by
// induction and the bulk of the grid no longer needs to be zeroed.
// ---------------------------------------------------------------------------
template <size_t NN>
__global__ void cudaFusedStreamMacro (bool const * __restrict__ IsSolid, real const * __restrict__ Gamma, real * __restrict__ Omeis,
        real const * __restrict__ Fpost, real * __restrict__ F, real3 * __restrict__ BForce, real3 * __restrict__ Vel,
        real * __restrict__ Rho, FLBM::lbm_aux const * __restrict__ lbmaux)
{
    const size_t Nneigh = (NN>0) ? NN : lbmaux[0].Nneigh;
    const size_t Nx     = lbmaux[0].Nx;
    const size_t Ny     = lbmaux[0].Ny;
    const size_t Nz     = lbmaux[0].Nz;
    const size_t Ncells = lbmaux[0].Ncells;
    const size_t Nxy    = Nx*Ny;
    const real   Cs     = lbmaux[0].Cs;

    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*Ncells) return;

    size_t icx =  ic%Nx;
    size_t icy = (ic/Nx)%Ny;
    size_t icz = (ic/Nxy)%Nz;
    size_t icl =  ic/Ncells;

    BForce[ic] = make_real3(0.0,0.0,0.0);
    Rho    [ic] = 0.0;
    Vel    [ic] = make_real3(0.0,0.0,0.0);

    const size_t ic0 = ic - icl*Ncells;
    real  rho = 0.0;
    real3 vel = make_real3(0.0,0.0,0.0);

    #pragma unroll
    for (size_t k=0;k<Nneigh;k++)
    {
        // Periodic wrap.  (icx - Cx) lies in (-Nx, 2Nx) for every supported
        // lattice, so at most one adjustment per axis is needed.  An integer
        // modulo by a run-time divisor costs ~20 instructions and this loop runs
        // three of them per discrete velocity, i.e. 45 per cell; the conditional
        // form is exact and replaces each with a predictably-predicated add.
        int inx = (int)icx - (int)lbmaux[0].C[k].x;
        int iny = (int)icy - (int)lbmaux[0].C[k].y;
        int inz = (int)icz - (int)lbmaux[0].C[k].z;
        if (inx<0) inx += (int)Nx; else if (inx>=(int)Nx) inx -= (int)Nx;
        if (iny<0) iny += (int)Ny; else if (iny>=(int)Ny) iny -= (int)Ny;
        if (inz<0) inz += (int)Nz; else if (inz>=(int)Nz) inz -= (int)Nz;
        size_t in0 = (size_t)inx + (size_t)iny*Nx + (size_t)inz*Nxy;
        real   f   = Fpost[FIDX(Ncells,Nneigh,icl,k,in0)];
        F[FIDX(Ncells,Nneigh,icl,k,ic0)] = f;
        rho        += f;
        vel         = vel + f*lbmaux[0].C[k];
    }

    if (!IsSolid[ic])
    {
        Rho[ic] = rho;
        Vel[ic] = Cs/rho*vel;
        if (Gamma[ic]>1.0e-12)
        {
            #pragma unroll
            for (size_t k=0;k<Nneigh;k++) Omeis[FIDX(Ncells,Nneigh,icl,k,ic0)] = 0.0;
        }
    }
}

// cudaCollideSCDEM<NN,Smooth>
//
// NN     : number of discrete velocities, known at compile time so that the
//          population loops are fully unrolled and NonEq[] is kept in
//          registers instead of being spilled to local memory.
//          NN==0 keeps the previous run-time bounded behaviour.
// Smooth : selects BnSmooth (true) or BnLadd (false). This used to be a
//          __device__ function pointer passed as a kernel argument, which the
//          compiler had to turn into an indirect call (no inlining, forced
//          register spilling around the call) for every single cell.
template <size_t NN, bool Smooth>
__global__ void cudaCollideSCDEM(bool const * __restrict__ IsSolid, real * __restrict__ F, real * __restrict__ Ftemp,
        real3 const * __restrict__ BForce, real3 const * __restrict__ Vel, real const * __restrict__ Rho,
        real const * __restrict__ Gamma, real const * __restrict__ Omeis, FLBM::lbm_aux const * __restrict__ lbmaux)
{
    const size_t Nneigh = (NN>0) ? NN : lbmaux[0].Nneigh;
    const size_t Ncells = lbmaux[0].Ncells;
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;

    if (!IsSolid[ic])
    {
        const real  dt     = lbmaux[0].dt;
        const real  Cs     = lbmaux[0].Cs;
        const real  Cs2    = Cs*Cs;
        const real  tau0   = lbmaux[0].Tau[0];
        const real  Sc     = lbmaux[0].Sc;
        const real  rho    = Rho[ic];
        const real  gamma  = Gamma[ic];

        real3 vel   = Vel[ic]+dt*(tau0/rho)*BForce[ic];
        real  tau   = tau0;
        real  Bn    = Smooth ? BnSmooth(gamma,tau) : BnLadd(gamma,tau);

        real  NonEq[(NN>0)?NN:27];
        real  Q = 0.0;

        // Loop invariants of the equilibrium distribution. Fc is exactly the
        // sub-expression the compiler used to hoist out of FLBM::FeqFluid, so
        // keeping it explicit costs nothing and stays bit-for-bit identical.
        const real VdotV = dotreal3(vel,vel);
        const real Fc    = 1.5*VdotV/Cs2;
        // The two per-direction divisions are replaced by products with the
        // loop-invariant reciprocals.  Not bit-for-bit identical with the
        // literal division (each quotient may move by <=1 ulp), but it removes
        // 2 double divisions per discrete velocity per cell, which is a real
        // cost on GPUs that run FP64 at a reduced rate.
        const real iCs   = 1.0/Cs;
        const real iCs2  = 1.0/Cs2;
        #pragma unroll
        for (size_t k=0;k<Nneigh;k++)
        {
            // Same as FLBM::FeqFluid(k,rho,vel,lbmaux) with the loop invariants
            // hoisted out (identical arithmetic ordering).
            real VdotC = dotreal3(vel,lbmaux[0].C[k]);
            real Feq   = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC*iCs + 4.5*VdotC*VdotC*iCs2 - Fc);
            NonEq[k]   = F[FIDX(Ncells,Nneigh,0,k,ic)] - Feq;
            Q         += NonEq[k]*NonEq[k]*lbmaux[0].EEk[k];
        }
        Q = sqrt(2.0*Q);
        tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));
        // tau is fixed for this cell, so a single reciprocal replaces the Nneigh
        // divisions below (`x/tau` -> `x*itau`; each quotient may move <=1 ulp).
        const real itau = 1.0/tau;

        // Fold the (constant during the fixed point iteration) non-equilibrium
        // term into NonEq[] once. The original recomputed
        //     noneq = (1-Bn)*NonEq[k]/tau - Bn*Omeis[k]
        // from NonEq[] and Omeis[] on every pass of the loop below; since Bn,
        // tau, NonEq[] and Omeis[] do not change inside the loop the value is
        // identical, so the recomputation and the extra Omeis reads are free to
        // drop.
        #pragma unroll
        for (size_t k=0;k<Nneigh;k++)
        {
            // The bulk of the grid has gamma==0 hence Bn==0, for which the
            // Omeis term vanishes exactly. Skipping the load keeps a 15x wider
            // read of Omeis out of the whole interior of the domain.
            real Ome = 0.0;
            if (Bn!=0.0) Ome = Omeis[FIDX(Ncells,Nneigh,0,k,ic)];
            NonEq[k] = (1.0 - Bn)*NonEq[k]*itau - Bn*Ome;
        }

        bool valid = true;
        real alpha = 1.0;
        size_t numit = 0;
        while (valid&&numit<2)
        {
            valid = false;
            #pragma unroll
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[FIDX(Ncells,Nneigh,0,k,ic)] = F[FIDX(Ncells,Nneigh,0,k,ic)] - alpha*(NonEq[k]);
                if (Ftemp[FIDX(Ncells,Nneigh,0,k,ic)]<0.0)
                {
                    real temp = F[FIDX(Ncells,Nneigh,0,k,ic)]/NonEq[k];
                    if (temp<alpha) alpha = temp;
                    valid = true;
                }
            }
            if (valid) numit++;
        }
    }
    else
    {
        #pragma unroll
        for (size_t k=0;k<Nneigh;k++)
        {
            Ftemp[FIDX(Ncells,Nneigh,0,k,ic)] = F[FIDX(Ncells,Nneigh,0,lbmaux[0].Op[k],ic)]; 
        }
    }
}

template <size_t NN, bool Smooth>
__global__ void cudaCollideMPDEM(bool const * __restrict__ IsSolid, real * __restrict__ F, real * __restrict__ Ftemp,
        real3 const * __restrict__ BForce, real3 const * __restrict__ Vel, real const * __restrict__ Rho,
        real const * __restrict__ Gamma, real const * __restrict__ Omeis, FLBM::lbm_aux const * __restrict__ lbmaux)
{
    const size_t Nneigh = (NN>0) ? NN : lbmaux[0].Nneigh;
    const size_t Ncells = lbmaux[0].Ncells;
    const real   Cs0    = lbmaux[0].Cs;
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;

    real3 Vmix = make_real3(0.0,0.0,0.0);
    real  den  = 0.0;
    for(size_t il=0;il<lbmaux[0].Nl;il++)
    {
        Vmix = Vmix + (Rho[ic+il*lbmaux[0].Ncells]/lbmaux[0].Tau[il])*Vel[ic+il*lbmaux[0].Ncells];
        den  = den  + Rho[ic+il*lbmaux[0].Ncells]/lbmaux[0].Tau[il];
    }
    Vmix = Vmix/den;

    real  gamma = Gamma[ic];
    for(size_t il=0;il<lbmaux[0].Nl;il++)
    {
        if (!IsSolid[ic+il*lbmaux[0].Ncells])
        {
            real  rho   = Rho[ic+il*lbmaux[0].Ncells];
            real3 vel   = Vmix + (lbmaux[0].dt*lbmaux[0].Tau[il]/rho)*BForce[ic+il*lbmaux[0].Ncells];
            real  VdotV = dotreal3(vel,vel);
            real  tau   = lbmaux[0].Tau[il];
            // tau is fixed for this cell, so a single reciprocal replaces the Nneigh
            // divisions below (`x/tau` -> `x*itau`; each quotient may move <=1 ulp).
            const real itau = 1.0/tau;
            real  Bn    = Smooth ? BnSmooth(gamma,tau) : BnLadd(gamma,tau);
            const real Fc = 1.5*VdotV/(Cs0*Cs0);
            bool valid = true;
            real alphal = 1.0;
            real alphat = 1.0;
            size_t numit = 0;
            while (valid)
            {
                numit++;
                valid = false;
                alphal = alphat;
                #pragma unroll
                for (size_t k=0;k<Nneigh;k++)
                {
                    real VdotC = dotreal3(vel,lbmaux[0].C[k]);
                    const real iCs   = 1.0/Cs0;
                    const real iCs2  = 1.0/(Cs0*Cs0);
                    real Feq   = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC*iCs + 4.5*VdotC*VdotC*iCs2 - Fc);
                    size_t idx = FIDX(Ncells,Nneigh,il,k,ic);
                    real Ome   = Omeis[idx];
                    real NonEq = (1.0-Bn)*(F[idx]-Feq)*itau - Bn*Ome;
                    Ftemp[idx] = F[idx] - alphal*NonEq;
                    if (Ftemp[idx]<0.0&&numit<2)
                    {
                        real temp = F[idx]/NonEq;
                        if (temp<alphat) alphat = temp;
                        valid = true;
                    }
                    if (Ftemp[idx]<0.0&&numit>=2)
                    {
                        Ftemp[idx] = 0.0;
                    }
                }
            }
        }
        else
        {
            #pragma unroll
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[FIDX(Ncells,Nneigh,il,k,ic)] = F[FIDX(Ncells,Nneigh,il,lbmaux[0].Op[k],ic)]; 
            }
        }
    } 
}

template <typename FuncBn>
__global__ void cudaImprintLatticeVC(FuncBn fBn, size_t const * PaCeV, DEM::ParticleCU * Par, DEM::DynParticleCU * DPar, real const * Rho, real * Gamma, real * Omeis, real const *
        F, DEM::dem_aux const * demaux, FLBM::lbm_aux const * lbmaux, lbmdem_aux const * lbmdemaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmdemaux[0].nvc) return;

    size_t ice = PaCeV[2*ic];
    size_t icx =  ice%lbmaux[0].Nx;
    size_t icy = (ice/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ice/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;
    const size_t Ncells_ = lbmaux[0].Ncells;
    const size_t Nn_     = lbmaux[0].Nneigh;
    const size_t icl_    = ice/Ncells_;
    const size_t ice0_   = ice - icl_*Ncells_;

    size_t ip  = PaCeV[2*ic+1];

    real3  C = lbmaux[0].dx*make_real3(real(icx),real(icy),real(icz));
    real3  Xs,B;
    real   len = 12.0*lbmaux[0].dx;
    real3  Pert = demaux[0].Per;
    bool isfree = ((!Par[ip].vxf&&!Par[ip].vyf&&!Par[ip].vzf&&!Par[ip].wxf&&!Par[ip].wyf&&!Par[ip].wzf)||Par[ip].FixFree);
    if (!isfree) Pert = make_real3(0.0,0.0,0.0);
    DEM::BranchVec(DPar[ip].x,C,B,Pert);
    if (norm(B)>Par[ip].Dmax+1.74*lbmaux[0].dx) return;
    Xs  = C-B;
    len = cudaSphereCube(Xs,C,Par[ip].Dmax,lbmaux[0].dx);
    if (fabs(len)<1.0e-12) return;
    real gamma  = len/(12.0*lbmaux[0].dx);
    if (gamma<Gamma[ice]) return;
    Gamma[ice] = gamma;
    real3 tmp;
    Rotation(DPar[ip].w,DPar[ip].Q,tmp);
    real3 VelP   = DPar[ip].v + cross(tmp,B);

    size_t ncells = lbmaux[0].Nneigh;
    real3 Flbm = make_real3(0.0,0.0,0.0);
    for (size_t il=0;il<lbmaux[0].Nl;il++)
    {
        real Bn  = fBn(gamma,lbmaux[0].Tau[il]);
        real rho = Rho[ice + il*lbmaux[0].Ncells];
        for (size_t k=0;k<ncells;k++)
        {
            real Fvpp     = FLBM::FeqFluid(lbmaux[0].Op[k],rho,VelP,lbmaux);
            real Fvp      = FLBM::FeqFluid(k              ,rho,VelP,lbmaux);
            real Omega    = F[FIDX(Ncells_,Nn_,icl_+il,lbmaux[0].Op[k],ice0_)] - Fvpp - (F[FIDX(Ncells_,Nn_,icl_+il,k,ice0_)] - Fvp);
            Omeis[FIDX(Ncells_,Nn_,icl_+il,k,ice0_)] = Omega;
            Flbm = Flbm - lbmdemaux[0].Fconv*Bn*Omega*lbmaux[0].Cs*lbmaux[0].Cs*lbmaux[0].dx*lbmaux[0].dx*lbmaux[0].C[k];
        }
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
    const size_t Ncells_ = lbmaux[0].Ncells;
    const size_t Nn_     = lbmaux[0].Nneigh;
    const size_t icl_    = ice/Ncells_;
    const size_t ice0_   = ice - icl_*Ncells_;

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

    //if (icx==50&&icy==50&&icz==6) printf("gamma: %g dist = %g igeo = %d ip = %d \n",len/(12.0*lbmaux[0].dx),minl,igeo,ip);
    if (fabs(len)<1.0e-12) return;
    real gamma  = len/(12.0*lbmaux[0].dx);
    if (gamma<Gamma[ice]) return;
    Gamma[ice] = gamma;
    real3 tmp;
    Rotation(DPar[ip].w,DPar[ip].Q,tmp);
    real3 VelP   = DPar[ip].v + cross(tmp,B);
    //real Bn  = gamma;
    size_t ncells = lbmaux[0].Nneigh;
    real3 Flbm = make_real3(0.0,0.0,0.0);
    for (size_t il=0;il<lbmaux[0].Nl;il++)
    {
        real Bn  = (gamma*(lbmaux[0].Tau[il]-0.5))/((1.0-gamma)+(lbmaux[0].Tau[il]-0.5));
        real rho = Rho[ice + il*lbmaux[0].Ncells];
        for (size_t k=0;k<ncells;k++)
        {
            real Fvpp     = FLBM::FeqFluid(lbmaux[0].Op[k],rho,VelP,lbmaux);
            real Fvp      = FLBM::FeqFluid(k              ,rho,VelP,lbmaux);
            real Omega    = F[FIDX(Ncells_,Nn_,icl_,lbmaux[0].Op[k],ice0_)] - Fvpp - (F[FIDX(Ncells_,Nn_,icl_,k,ice0_)] - Fvp);
            Omeis[FIDX(Ncells_,Nn_,icl_,k,ice0_)] = Omega;
            Flbm = Flbm - Bn*Omega*lbmaux[0].Cs*lbmaux[0].Cs*lbmaux[0].dx*lbmaux[0].dx*lbmaux[0].C[k];
        }
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

__global__ void cudaApplyForcesSC(uint3 * pCellPairs, bool const * IsSolid, real * Gamma, real3 * BForce, real const * Rho, FLBM::lbm_aux const * lbmaux)
{
    size_t icp = threadIdx.x + blockIdx.x * blockDim.x;
    if (icp>=lbmaux[0].NCPairs) return;
    size_t ic = pCellPairs[icp].x;
    size_t in = pCellPairs[icp].y;
    size_t k  = pCellPairs[icp].z;

    for (size_t il=0;il<lbmaux[0].Nl;il++)
    {
        real psic = 1.0;
        real psin = 1.0;
        real G    = lbmaux[0].G[il];
        bool chksolidc = false;
        bool chksolidn = false;
        if (fabs(G)<1.0e-12) continue;
        if (!IsSolid[ic]&&fabs(1.0-Gamma[ic])>1.0e-12) psic = lbmaux[0].Psi[il]*exp(-lbmaux[0].Rhoref[il]/Rho[ic+il*lbmaux[0].Ncells]);
        else
        {            
            G    = lbmaux[0].Gs[il];
            chksolidc = true;
        }
        if (!IsSolid[in]&&fabs(1.0-Gamma[in])>1.0e-12) psin = lbmaux[0].Psi[il]*exp(-lbmaux[0].Rhoref[il]/Rho[in+il*lbmaux[0].Ncells]);
        else
        {
            G    = lbmaux[0].Gs[il];
            chksolidn = true;
        }

        if (chksolidc&&chksolidn) G = 0.0;
        
        real3 bforce = (-G*lbmaux[0].W[k]*psic*psin)*lbmaux[0].C[k];

        atomicAdd(&BForce[ic+il*lbmaux[0].Ncells].x, bforce.x);
        atomicAdd(&BForce[ic+il*lbmaux[0].Ncells].y, bforce.y);
        atomicAdd(&BForce[ic+il*lbmaux[0].Ncells].z, bforce.z);
        atomicAdd(&BForce[in+il*lbmaux[0].Ncells].x,-bforce.x);
        atomicAdd(&BForce[in+il*lbmaux[0].Ncells].y,-bforce.y);
        atomicAdd(&BForce[in+il*lbmaux[0].Ncells].z,-bforce.z);
    }

    //printf("BForce %d %f %f %f \n",ic,BForce[ic].x,BForce[ic].y,BForce[ic].z);
}

__global__ void cudaApplyForcesSCMP(uint3 * pCellPairs, bool const * IsSolid, real * Gamma, real3 * BForce, real const * Rho, FLBM::lbm_aux const * lbmaux)
{
    size_t icp = threadIdx.x + blockIdx.x * blockDim.x;
    if (icp>=lbmaux[0].NCPairs) return;
    size_t ic = pCellPairs[icp].x;
    size_t in = pCellPairs[icp].y;
    size_t k  = pCellPairs[icp].z;

    bool chksolidc = false;
    bool chksolidn = false;

    for (size_t il=0;il<lbmaux[0].Nl;il++)
    {
        real psic = 0.0;
        real psin = 0.0;
        real G    = lbmaux[0].G[il];
        if (fabs(G)<1.0e-12) continue;
        if (!IsSolid[ic+il*lbmaux[0].Ncells]&&fabs(1.0-Gamma[ic])>1.0e-12) psic = lbmaux[0].Psi[il]*exp(-lbmaux[0].Rhoref[il]/Rho[ic+il*lbmaux[0].Ncells]);
        if (!IsSolid[in+il*lbmaux[0].Ncells]&&fabs(1.0-Gamma[in])>1.0e-12) psin = lbmaux[0].Psi[il]*exp(-lbmaux[0].Rhoref[il]/Rho[in+il*lbmaux[0].Ncells]);

        real3 bforce = (-G*lbmaux[0].W[k]*psic*psin)*lbmaux[0].C[k];

        atomicAdd(&BForce[ic+il*lbmaux[0].Ncells].x, bforce.x);
        atomicAdd(&BForce[ic+il*lbmaux[0].Ncells].y, bforce.y);
        atomicAdd(&BForce[ic+il*lbmaux[0].Ncells].z, bforce.z);
        atomicAdd(&BForce[in+il*lbmaux[0].Ncells].x,-bforce.x);
        atomicAdd(&BForce[in+il*lbmaux[0].Ncells].y,-bforce.y);
        atomicAdd(&BForce[in+il*lbmaux[0].Ncells].z,-bforce.z);
    }
    
    for (size_t il1=0    ;il1<lbmaux[0].Nl-1;il1++)
    for (size_t il2=il1+1;il2<lbmaux[0].Nl  ;il2++)
    {
        chksolidc = false;
        chksolidn = false;
        real psic = 1.0;
        real psin = 1.0;
        real G    = lbmaux[0].Gmix;
        if (!IsSolid[ic+il1*lbmaux[0].Ncells]&&fabs(1.0-Gamma[ic])>1.0e-12) psic = Rho[ic+il1*lbmaux[0].Ncells];
        else
        {
            G    = lbmaux[0].Gs[il2];
            chksolidc = true;
        }
        if (!IsSolid[in+il2*lbmaux[0].Ncells]&&fabs(1.0-Gamma[in])>1.0e-12) psin = Rho[in+il2*lbmaux[0].Ncells];
        else
        {
            G    = lbmaux[0].Gs[il1];
            chksolidn = true;
        }

        if (chksolidc&&chksolidn) G = 0.0;

        real3 bforce = (-G*lbmaux[0].W[k]*psic*psin)*lbmaux[0].C[k];

        atomicAdd(&BForce[ic+il1*lbmaux[0].Ncells].x, bforce.x);
        atomicAdd(&BForce[ic+il1*lbmaux[0].Ncells].y, bforce.y);
        atomicAdd(&BForce[ic+il1*lbmaux[0].Ncells].z, bforce.z);
        atomicAdd(&BForce[in+il2*lbmaux[0].Ncells].x,-bforce.x);
        atomicAdd(&BForce[in+il2*lbmaux[0].Ncells].y,-bforce.y);
        atomicAdd(&BForce[in+il2*lbmaux[0].Ncells].z,-bforce.z);

        chksolidc = false;
        chksolidn = false;
        psic = 1.0;
        psin = 1.0;
        G    = lbmaux[0].Gmix;
        if (!IsSolid[ic+il2*lbmaux[0].Ncells]&&fabs(1.0-Gamma[ic])>1.0e-12) psic = Rho[ic+il2*lbmaux[0].Ncells];
        else
        {            
            G    = lbmaux[0].Gs[il1];
            chksolidc = true;
        }
        if (!IsSolid[in+il1*lbmaux[0].Ncells]&&fabs(1.0-Gamma[in])>1.0e-12) psin = Rho[in+il1*lbmaux[0].Ncells];
        else
        {
            G    = lbmaux[0].Gs[il2];
            chksolidn = true;
        }

        if (chksolidc&&chksolidn) G = 0.0;

        bforce = (-G*lbmaux[0].W[k]*psic*psin)*lbmaux[0].C[k];

        atomicAdd(&BForce[ic+il2*lbmaux[0].Ncells].x, bforce.x);
        atomicAdd(&BForce[ic+il2*lbmaux[0].Ncells].y, bforce.y);
        atomicAdd(&BForce[ic+il2*lbmaux[0].Ncells].z, bforce.z);
        atomicAdd(&BForce[in+il1*lbmaux[0].Ncells].x,-bforce.x);
        atomicAdd(&BForce[in+il1*lbmaux[0].Ncells].y,-bforce.y);
        atomicAdd(&BForce[in+il1*lbmaux[0].Ncells].z,-bforce.z);
    }
    
}

__global__ void cudaStream2(bool const * IsSolid, real * Gamma, real * Omeis, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;
    if (ic==0)
    {
        lbmaux[0].Time += lbmaux[0].dt;
        lbmaux[0].iter++;
    }
    const size_t NcellsD = lbmaux[0].Ncells;
    const size_t NnD     = lbmaux[0].Nneigh;
    const size_t iclD    = ic/NcellsD;
    const size_t ic0D    = ic - iclD*NcellsD;
    BForce[ic] = make_real3(0.0,0.0,0.0);
    Rho   [ic] = 0.0;
    Vel   [ic] = make_real3(0.0,0.0,0.0);
    if (!IsSolid[ic])
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Omeis[FIDX(NcellsD,NnD,iclD,k,ic0D)] = 0.0;
            Rho[ic] += F[FIDX(NcellsD,NnD,iclD,k,ic0D)];
            Vel[ic] = Vel[ic] + F[FIDX(NcellsD,NnD,iclD,k,ic0D)]*lbmaux[0].C[k];
        }
        Vel[ic] = lbmaux[0].Cs/Rho[ic]*Vel[ic];
    }
}
}
#endif //MECHSYS_LBMDEM_CUH
