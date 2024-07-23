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

__global__ void cudaReset(bool const * IsSolid, real * Gammaf, real * Gamma, real * Omeis,FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;
    //if (Gamma[ic]<0) printf("ice %lu iter %lu \n",ic,lbmaux[0].iter);
    //if (ic==484) printf("Gamma %g Gammaf %g iter %lu \n",Gamma[ic],Gammaf[ic], lbmaux[0].iter);
    if (IsSolid[ic])  Gamma[ic] = 1.0;
    else              Gamma[ic] = Gammaf[ic];
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        Omeis[ic*lbmaux[0].Nneigh + k] = 0.0;
    }
}

#ifdef USE_IBB
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

    size_t icx =  ice%lbmaux[0].Nx;
    size_t icy = (ice/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ice/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;
    

    if (Inside[ice]<=-2)
    {
        for (size_t k = 0; k < lbmaux[0].Nneigh ; k++)
        {
            F[ice*lbmaux[0].Nneigh + k] = 0.0;
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
                F[ice*lbmaux[0].Nneigh + kt] += 2.0*F[in*lbmaux[0].Nneigh + kt] - F[im*lbmaux[0].Nneigh + kt];
                rhon += F[in*lbmaux[0].Nneigh + kt];
                rhom += F[im*lbmaux[0].Nneigh + kt];
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
                F[ice*lbmaux[0].Nneigh + k] = FeqFluid(k,rho,VelP,lbmaux);
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
            F[ice*lbmaux[0].Nneigh + k] = fabs(F[ice*lbmaux[0].Nneigh + k])/naem;
            Rho[ice] += F[ice*lbmaux[0].Nneigh + k];
            Vel[ice] = Vel[ice] + F[ice*lbmaux[0].Nneigh + k]*lbmaux[0].C[k];
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
            //real Omega   = F[ice*lbmaux[0].Nneigh + lbmaux[0].Op[k]] - Fvpp - (F[ice*lbmaux[0].Nneigh + k] - Fvp);
            //Omeis[ice*lbmaux[0].Nneigh + k] = Omega;
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
        F[ice*lbmaux[0].Nneigh + ko] = (r*F[ife*lbmaux[0].Nneigh + ko] + (1.0-r)*F[ice*lbmaux[0].Nneigh + k] 
                + r*F[ine*lbmaux[0].Nneigh + k] + 6.0*rho*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[ko],VelP)/Cs)/(1.0+r);


        //F[ine*lbmaux[0].Nneigh + k]  = F[ice*lbmaux[0].Nneigh + k];

        Flbm = Flbm + ld*ld*Cs*(F[ine*lbmaux[0].Nneigh + k]*(Cs*lbmaux[0].C[k]-VelP) - F[ice*lbmaux[0].Nneigh + ko]*(Cs*lbmaux[0].C[ko]-VelP));
        //Flbm = Flbm + ld*ld*Cs*(F[ine*lbmaux[0].Nneigh + k]*(Cs*lbmaux[0].C[k]) - F[ice*lbmaux[0].Nneigh + ko]*(Cs*lbmaux[0].C[ko]));
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
#else
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
    if (norm(B)>Par[ip].Dmax+1.74*lbmaux[0].dx) return;
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
#ifndef USE_LADD
    real Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
#else
    real Bn  = floor(gamma);
#endif
    size_t ncells = lbmaux[0].Nneigh;
    real3 Flbm = make_real3(0.0,0.0,0.0);
    for (size_t k=0;k<ncells;k++)
    {
        real Fvpp     = FLBM::FeqFluid(lbmaux[0].Op[k],rho,VelP,lbmaux);
        real Fvp      = FLBM::FeqFluid(k              ,rho,VelP,lbmaux);
        real Omega    = F[ice*lbmaux[0].Nneigh + lbmaux[0].Op[k]] - Fvpp - (F[ice*lbmaux[0].Nneigh + k] - Fvp);
        Omeis[ice*lbmaux[0].Nneigh + k] = Omega;
        Flbm = Flbm - lbmdemaux[0].Fconv*Bn*Omega*lbmaux[0].Cs*lbmaux[0].Cs*lbmaux[0].dx*lbmaux[0].dx*lbmaux[0].C[k];
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
#endif

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
    #ifndef USE_IBB
    if (ic==0)
    {
        lbmaux[0].Time += lbmaux[0].dt;
        lbmaux[0].iter++;
    }
    BForce[ic] = make_real3(0.0,0.0,0.0);
    #endif
    //for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    //{
        //F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k];
    //}
    Rho   [ic] = 0.0;
    Vel   [ic] = make_real3(0.0,0.0,0.0);
    //Gamma[ic]  = (real) IsSolid[ic];
    if (!IsSolid[ic])
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ic] += F[ic*lbmaux[0].Nneigh + k];
            Vel[ic] = Vel[ic] + F[ic*lbmaux[0].Nneigh + k]*lbmaux[0].C[k];
        }
        Vel[ic] = lbmaux[0].Cs/Rho[ic]*Vel[ic];
    }
}
}
#endif //MECHSYS_LBMDEM_CUH
