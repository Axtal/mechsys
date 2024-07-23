/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2016 Sergio Galindo                                    *
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

/////////////////////////////LBM CUDA implementation////////////////////

#ifndef MECHSYS_LBM_CUH
#define MECHSYS_LBM_CUH

//Mechsys
#include <mechsys/linalg/matvec.h>

namespace FLBM
{
struct lbm_aux
{
    size_t     Nl;          ///< Number of Lattices
    size_t     Nneigh;      ///< Number of Neighbors
    size_t     NCPairs;     ///< Number of cell pairs
    size_t     Nx;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Ny;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Nz;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Ncells;      ///< Integer vector with the dimensions of the LBM domain
    size_t     Op[27];      ///< Array with the opposite directions for bounce back calculation
    size_t     iter;        ///< counter for the number of iterations
    real3      C[27];       ///< Collection of discrete velocity vectors
    real       EEk[27];     ///< Dyadic product of discrete velocities for LES calculation
    real       W[27];       ///< Collection of discrete weights
    real       M[729];      ///< Matrix for MRT
    real       Mi[729];     ///< Inverse matrix for MRT
    real       S[27];       ///< Vector with relaxation times for MRT
    real       Tau[3];      ///< Collection of characteristic collision times
    real       Cs;          ///< Lattice speed
    real       Sc;          ///< Smagorinsky constant
    real       dx;          ///< grid size
    real       dt;          ///< time step
    real       Time;        ///< Time clock
    //These parameters are for Shan Chen type of simulations
    real       G[2];        ///< Collection of cohesive constants for multiphase simulation
    real       Gs[2];       ///< Collection of cohesive constants for multiphase simulation
    real       Rhoref[2];   ///< Collection of cohesive constants for multiphase simulation
    real       Psi[2];      ///< Collection of cohesive constants for multiphase simulation
    real       Gmix;        ///< Repulsion constant for multicomponent simulation
    //These parameters are for the Phase Field Ice model 0 solid 1 liquid 2 gas
    real       rho[3];      ///< Density of the phases
    real       cap[3];      ///< Heat capcity for each phase
    real       kap[3];      ///< heat conductivity for each phase    
    real       thick;       ///< thickness of the phase field interfase;
    real       sigma;       ///< surface tension of the interfase
    real       Ts;          ///< Solidus temperature
    real       Tl;          ///< Liquidus temperature
    real       L;           ///< Latent heat
};

__device__ real FeqFluid(size_t const & k, real const & rho, real3 const & vel, lbm_aux const * lbmaux)
{
    real VdotV = dotreal3(vel,vel);
    real VdotC = dotreal3(vel,lbmaux[0].C[k]);
    real Cs    = lbmaux[0].Cs;
    return lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

__global__ void cudaCheckUpLoad (lbm_aux const * lbmaux)
{
    /*
    printf("Nl          %lu \n",  lbmaux[0].Nl     );
    printf("Nneigh      %lu \n",  lbmaux[0].Nneigh );
    printf("NCP         %lu \n",  lbmaux[0].NCPairs);
    printf("Dim      %d %lu \n",0,lbmaux[0].Nx );
    printf("Dim      %d %lu \n",1,lbmaux[0].Ny );
    printf("Dim      %d %lu \n",2,lbmaux[0].Nz );
    printf("Ncells      %lu \n"  ,lbmaux[0].Ncells );
    printf("Sc          %f \n"   ,lbmaux[0].Sc );
    printf("Cs          %f \n"   ,lbmaux[0].Cs );
    printf("dx          %f \n"   ,lbmaux[0].dx );
    printf("dt          %f \n"   ,lbmaux[0].dt );
    
    for (size_t i=0;i < lbmaux[0].Nl;i++)
    {
        printf("Tau     %d %f \n",  i, lbmaux[0].Tau[i]   );
        printf("G       %d %f \n",  i, lbmaux[0].G[i]     );
        printf("Gs      %d %f \n",  i, lbmaux[0].Gs[i]    );
    }

    for (size_t i=0;i < lbmaux[0].Nneigh;i++)
    {
        printf("C      %d %f %f %f \n",i,lbmaux[0].C[i].x,lbmaux[0].C[i].y,lbmaux[0].C[i].z);
    }
    for (size_t i=0;i < lbmaux[0].Nneigh;i++)
    {
        printf("Wk     %d %f       \n",i,lbmaux[0].W[i]);
        printf("EEk    %d %f       \n",i,lbmaux[0].EEk[i]);
        printf("Op     %d %lu      \n",i,lbmaux[0].Op[i]);
        printf("Sk     %d %f       \n",i,lbmaux[0].S[i]);
    }
    for (size_t i=0;i < lbmaux[0].Nneigh;i++)
    {
        for (size_t j=0;j < lbmaux[0].Nneigh;j++)
        {
            printf("M      %d %lu %f    \n",i,j,lbmaux[0].M [j+i*lbmaux[0].Nneigh]);
        }
    }
    for (size_t i=0;i < lbmaux[0].Nneigh;i++)
    {
        for (size_t j=0;j < lbmaux[0].Nneigh;j++)
        {
            printf("Mi     %d %lu %f    \n",i,j,lbmaux[0].Mi[j+i*lbmaux[0].Nneigh]);
        }
    }

    //real Feq = FeqFluid(3,1.0,(real3)(0.2,0.0,0.0),lbmaux);
    //printf(" %f \n",Feq);
    */
}

__global__ void cudaCollideSC(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux const * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;

    if (!IsSolid[ic])
    {
        real3 vel   = Vel[ic]+lbmaux[0].dt*(lbmaux[0].Tau[0]/Rho[ic])*BForce[ic];
        real  rho   = Rho[ic];
        real  tau   = lbmaux[0].Tau[0];
        
        real  NonEq[27];
        //real  Feq  [27];
        real  Q = 0.0;
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            //real VdotC = dot(vel,lbmaux[0].C[k]);
            //Feq  [k]     = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
            //Feq[k]       = FeqFluid(k,rho,vel,lbmaux);
            //NonEq[k]     = F[ic*lbmaux[0].Nneigh + k] - Feq[k];
            NonEq[k]     = F[ic*lbmaux[0].Nneigh + k] - FeqFluid(k,rho,vel,lbmaux);
            Q           += NonEq[k]*NonEq[k]*lbmaux[0].EEk[k];
        }
        Q = sqrt(2.0*Q);
        tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*lbmaux[0].Sc/rho));

        bool valid = true;
        real alpha = 1.0;
        size_t numit = 0;
        while (valid&&numit<2)
        {
            valid = false;
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k] - alpha*(NonEq[k]/tau);
                if (Ftemp[ic*lbmaux[0].Nneigh + k]<0.0)
                {
                    //real temp = F[ic*lbmaux[0].Nneigh + k]/(NonEq[k]/tau - Fk);
                    real temp = tau*F[ic*lbmaux[0].Nneigh + k]/(NonEq[k]);
                    if (temp<alpha) alpha = temp;
                    valid = true;
                }
            }
            if (valid) numit++;
        }
    }
    else
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
        }
    }
    //for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    //{
        //F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k]; 
    //}
}

__global__ void cudaCollideSCDEM(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, real * Gamma, real * Omeis, lbm_aux const * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;

    if (!IsSolid[ic])
    {
        real3 vel   = Vel[ic]+lbmaux[0].dt*(lbmaux[0].Tau[0]/Rho[ic])*BForce[ic];
        real  rho   = Rho[ic];
        real  tau   = lbmaux[0].Tau[0];
        real  gamma = Gamma[ic];
#ifndef USE_LADD
        real  Bn    = (gamma*(tau-0.5))/((1.0-gamma)+(tau-0.5));
        //real Bn = gamma;
#else
        real Bn = floor(gamma);
#endif

        real  NonEq[27];
        //real  Feq  [27];
        real  Q = 0.0;

        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            //real VdotC = dot(vel,lbmaux[0].C[k]);
            //Feq  [k]     = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
            //Feq[k]       = FeqFluid(k,rho,vel,lbmaux);
            //NonEq[k]     = F[ic*lbmaux[0].Nneigh + k] - Feq[k];
            NonEq[k]     = F[ic*lbmaux[0].Nneigh + k] - FeqFluid(k,rho,vel,lbmaux);
            Q           += NonEq[k]*NonEq[k]*lbmaux[0].EEk[k];
        }
        Q = sqrt(2.0*Q);
        tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*lbmaux[0].Sc/rho));

        bool valid = true;
        real alpha = 1.0;
        size_t numit = 0;
        while (valid&&numit<2)
        {
            valid = false;
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                real Ome   = Omeis[ic*lbmaux[0].Nneigh + k];
                #ifndef USE_IBB
                real noneq = (1.0 - Bn)*NonEq[k]/tau-Bn*Ome;
                #else
                real noneq = NonEq[k]/tau;
                #endif
                Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k] - alpha*(noneq);
                //if((ic==555||ic==556||ic==557||ic==558)&&(k==1)) printf("ic %lu Ftk %g Fk %g iter %lu \n",ic,Ftemp[ic*lbmaux[0].Nneigh + k],F[ic*lbmaux[0].Nneigh + k],lbmaux[0].iter); 
                if (Ftemp[ic*lbmaux[0].Nneigh + k]<0.0)
                {
                    //real temp = F[ic*lbmaux[0].Nneigh + k]/(NonEq[k]/tau - Fk);
                    real temp = F[ic*lbmaux[0].Nneigh + k]/noneq;
                    if (temp<alpha) alpha = temp;
                    valid = true;
                }
            }
            if (valid) numit++;
        }
    }
    else
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
        }
    }
    //for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    //{
        //F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k]; 
    //}
}

__global__ void cudaCollideMP(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux const * lbmaux)
{
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

    for(size_t il=0;il<lbmaux[0].Nl;il++)
    {
        if (!IsSolid[ic+il*lbmaux[0].Ncells])
        {
            real  rho   = Rho[ic+il*lbmaux[0].Ncells];
            real3 vel   = Vmix + (lbmaux[0].dt*lbmaux[0].Tau[il]/rho)*BForce[ic+il*lbmaux[0].Ncells];
            real  VdotV = dotreal3(vel,vel);
            real  tau   = lbmaux[0].Tau[il];
            bool valid = true;
            real alphal = 1.0;
            real alphat = 1.0;
            size_t numit = 0;
            while (valid)
            {
                numit++;
                valid = false;
                alphal = alphat;
                for (size_t k=0;k<lbmaux[0].Nneigh;k++)
                {
                    real VdotC = dotreal3(vel,lbmaux[0].C[k]);
                    real Cs    = lbmaux[0].Cs;
                    real Feq   = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                    size_t idx = ic*lbmaux[0].Nneigh + il*lbmaux[0].Ncells*lbmaux[0].Nneigh + k;
                    Ftemp[idx] = F[idx] - alphal*(F[idx]-Feq)/tau;
                    if (Ftemp[idx]<0.0&&numit<2)
                    {
                        real temp = tau*fabs(F[idx]/(F[idx]-Feq));
                        if (temp<alphat) alphat = temp;
                        valid = true;
                    }
                    if (Ftemp[idx]<0.0&&numit>=2)
                    {
                        Ftemp[idx] = 0.0;
                        //size_t icx = ic%lbmaux[0].Nx;
                        //size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
                        //size_t icz = (ic/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;
                        //printf("%lu %lu %lu %g %lu %g \n",icx,icy,icz,Ftemp[idx],lbmaux[0].iter,alphal);
                    }
                }
            }
        }
        else
        {
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                Ftemp[ic*lbmaux[0].Nneigh + il*lbmaux[0].Ncells*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + il*lbmaux[0].Ncells*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
            }
        }
    } 
}

__global__ void cudaCollideSC_MRT(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux const * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;
    if (!IsSolid[ic])
    {
        real  Cs  = lbmaux[0].Cs;
        real  dx  = lbmaux[0].dx;
        real  dt  = lbmaux[0].dt;
        real  rho = Rho[ic];
        real3 vel = Vel[ic];

        real * f  = F     + ic*lbmaux[0].Nneigh;
        real * ft = Ftemp + ic*lbmaux[0].Nneigh;
        real fneq[27];

        MtVecMul(lbmaux[0].M,f,ft,lbmaux[0].Nneigh,lbmaux[0].Nneigh);

        ft[ 0] = 0.0; 
        ft[ 1] = lbmaux[0].S[ 1]*(ft[ 1] + rho - rho*dotreal3(vel,vel)/(Cs*Cs));
        ft[ 2] = lbmaux[0].S[ 2]*(ft[ 2] - rho);
        ft[ 3] = 0.0;
        ft[ 4] = lbmaux[0].S[ 4]*(ft[ 4] + 7.0/3.0*rho*vel.x/Cs); 
        ft[ 5] = 0.0;
        ft[ 6] = lbmaux[0].S[ 6]*(ft[ 6] + 7.0/3.0*rho*vel.y/Cs); 
        ft[ 7] = 0.0;
        ft[ 8] = lbmaux[0].S[ 8]*(ft[ 8] + 7.0/3.0*rho*vel.z/Cs); 
        ft[ 9] = lbmaux[0].S[ 9]*(ft[ 9] - rho*(2.0*vel.x*vel.x-vel.y*vel.y-vel.z*vel.z)/(Cs*Cs));
        ft[10] = lbmaux[0].S[10]*(ft[10] - rho*(vel.y*vel.y-vel.z*vel.z)/(Cs*Cs));
        ft[11] = lbmaux[0].S[11]*(ft[11] - rho*(vel.x*vel.y)/(Cs*Cs));
        ft[12] = lbmaux[0].S[12]*(ft[12] - rho*(vel.y*vel.z)/(Cs*Cs));
        ft[13] = lbmaux[0].S[13]*(ft[13] - rho*(vel.x*vel.z)/(Cs*Cs));
        ft[14] = lbmaux[0].S[14]* ft[14];

        MtVecMul(lbmaux[0].Mi,ft,fneq,lbmaux[0].Nneigh,lbmaux[0].Nneigh);

        bool valid = true;
        real alpha = 1.0;
        size_t numit = 0;
        while (valid&&numit<2)
        {
            valid = false;
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k] - alpha*(fneq[k] - 3.0*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[k],BForce[ic])*dt*dt/dx);
                if (Ftemp[ic*lbmaux[0].Nneigh + k]<0.0)
                {
                    real temp = F[ic*lbmaux[0].Nneigh + k]/(fneq[k] - 3.0*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[k],BForce[ic])*dt*dt/dx);
                    if (temp<alpha) alpha = temp;
                    valid = true;
                }
            }
            if (valid) numit++;
        }
    }
    else
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
        }
    }
}

__global__ void cudaCollideAD(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux const * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;

    real Cs    = lbmaux[0].Cs;
    for(size_t il=0;il<lbmaux[0].Nl;il++)
    {
        if (!IsSolid[ic+il*lbmaux[0].Ncells])
        {
            real  rho   = Rho[ic+il*lbmaux[0].Ncells];
            real3 vel   = Vel[ic] + (lbmaux[0].dt*lbmaux[0].Tau[il]/Rho[ic])*BForce[ic+il*lbmaux[0].Ncells]; //ic for the velocity because for the AD
                                                                                                         //solver only velocity of phase 0 is
                                                                                                         //important
            real  VdotV = dotreal3(vel,vel);
            real  tau   = lbmaux[0].Tau[il];
            bool valid = true;
            real alpha = 1.0;
            size_t numit = 0;
            while (valid&&numit<2)
            {
                valid = false;
                for (size_t k=0;k<lbmaux[0].Nneigh;k++)
                {
                    real VdotC = dotreal3(vel,lbmaux[0].C[k]);
                    real Feq   = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                    size_t idx = ic*lbmaux[0].Nneigh + il*lbmaux[0].Ncells*lbmaux[0].Nneigh + k;
                    Ftemp[idx] = F[idx] - alpha*(F[idx]-Feq)/tau;
                    if (Ftemp[idx]<0.0)
                    {
                        real temp = tau*F[idx]/(F[idx]-Feq);
                        if (temp<alpha) alpha = temp;
                        valid = true;
                    }
                }
                if (valid) numit++;
            }
        }
        else
        {
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                Ftemp[ic*lbmaux[0].Nneigh + il*lbmaux[0].Ncells*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + il*lbmaux[0].Ncells*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
            }
        }
    } 
}

__global__ void cudaCollidePFI(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux const * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;

    real Cs    = lbmaux[0].Cs;
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ic/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;
        //Calculation of gradients and temporal derivates of key variables


    //Phase field variables
    real3 vel   = Vel[ic]; // it is ic since only the velocity of the first layer is relavant
    real  phi   = Rho[ic+1*lbmaux[0].Ncells];
    real3 dFdt  = (phi*vel - Vel[ic+1*lbmaux[0].Ncells])/lbmaux[0].dt; //Here we will use the Vel of lattice one as a temporal variable for the
                                                                       //phase field flux
    real lambda = 4.0*phi*(1.0-phi)/lbmaux[0].thick;
    real a      = 1.5*lbmaux[0].sigma*lbmaux[0].thick;
    real b      = 12.0*lbmaux[0].sigma/lbmaux[0].thick;

    // Enthalpy variables
    real H    = Rho[ic + 2*lbmaux[0].Ncells];
    real Hs   = lbmaux[0].Ts*lbmaux[0].cap[0];
    real Hl   = lbmaux[0].Tl*lbmaux[0].cap[1]+lbmaux[0].L;
    real fl   = Vel[ic+2*lbmaux[0].Ncells].x; //For now the x component of this vector will hold the liquid fraction;
    real fs   = 1.0-fl;
    real fla  = Vel[ic+2*lbmaux[0].Ncells].y; //For now the x component of this vector will hold the liquid fraction;
    real dfdt = -(fl-fla)/lbmaux[0].dt;

    real Cp   = phi*(fs*lbmaux[0].cap[0] + fl*lbmaux[0].cap[1])+fl*(1.0-phi)*lbmaux[0].cap[2];
    real Temp = H/Cp;
    if      (H>=Hs&&H<=Hl) Temp = lbmaux[0].Ts + (H-Hs)/(Hl-Hs)*(lbmaux[0].Tl-lbmaux[0].Ts);
    else if (H>Hl)         Temp = lbmaux[0].Tl + (H-Hl)/Cp;

    //Density variables
    real rho  = phi*(fs*lbmaux[0].rho[0] + fl*lbmaux[0].rho[1])+fl*(1.0-phi)*lbmaux[0].rho[2];
    real pre  = Rho[ic];
    real S    = Vel[ic+2*lbmaux[0].Ncells].z;
    
    //Gradients
    real3 gradphi = make_real3(0.0,0.0,0.0);
    real3 gradrho = make_real3(0.0,0.0,0.0);
    real3 gradpre = make_real3(0.0,0.0,0.0);
    real3 gradS   = make_real3(0.0,0.0,0.0);
    real  delphi  = 0.0;

    for (size_t k=1;k<lbmaux[0].Nneigh;k++)
    {
        size_t inx = (size_t)((int)icx + (int)lbmaux[0].C[k].x + (int)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((int)icy + (int)lbmaux[0].C[k].y + (int)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((int)icz + (int)lbmaux[0].C[k].z + (int)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;

        real pren = Rho[in];
        real phin = Rho[in+1*lbmaux[0].Ncells];
        real fln  = Vel[in+2*lbmaux[0].Ncells].x;
        real fsn  = 1.0-fln;
        real rhon = phin*(fsn*lbmaux[0].rho[0] + fln*lbmaux[0].rho[1])+fln*(1.0-phin)*lbmaux[0].rho[2];
        real Sn   = Vel[in+2*lbmaux[0].Ncells].z;
        
        gradpre  = gradpre + lbmaux[0].W[k]*pren*lbmaux[0].C[k];
        gradphi  = gradphi + lbmaux[0].W[k]*phin*lbmaux[0].C[k];
        gradrho  = gradrho + lbmaux[0].W[k]*rhon*lbmaux[0].C[k];
        gradS    = gradS   + lbmaux[0].W[k]*Sn  *lbmaux[0].C[k];
        delphi  +=       2.0*lbmaux[0].W[k]*(phin-phi);
    }
    gradphi = (3.0/lbmaux[0].dx)*gradphi;
    gradrho = (3.0/lbmaux[0].dx)*gradrho;
    gradpre = (3.0/lbmaux[0].dx)*gradpre;
    gradS   = (3.0/lbmaux[0].dx)*gradS  ;
    delphi *= 3.0/(lbmaux[0].dx*lbmaux[0].dx);

    real3 n    = (1.0/norm(gradphi))*gradphi;
    if (norm(gradphi)<1.0e-12) n = make_real3(0.0,0.0,0.0);
    real  divu = (1.0-lbmaux[0].rho[0]/lbmaux[0].rho[1])*dfdt;
    
    real  tau   = lbmaux[0].Tau[0];
    real  taup  = lbmaux[0].Tau[1];
    real  kappa = phi*(fs*lbmaux[0].kap[0] + fl*lbmaux[0].kap[1])+fl*(1.0-phi)*lbmaux[0].kap[2];
    real  tauh  = 3.0*kappa*lbmaux[0].dt/(lbmaux[0].dx*lbmaux[0].dx) + 0.5;
    real  VdotV                      = dotreal3(vel,vel);
    real  mu                         = 4.0*b*(phi)*(phi-1.0)*(phi-0.5) - a*delphi;
    BForce[ic]                       = mu*gradphi;
    real3 Fm                         = BForce[ic] + Cs*Cs*gradrho/3.0 + Cs*Cs*gradS/3.0 - gradpre;
    BForce[ic + 1*lbmaux[0].Ncells]  = -rho*fs*(BForce[ic + 1*lbmaux[0].Ncells])/lbmaux[0].dt;
    real  VdotF = dotreal3(vel,Fm);
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        real VdotC = dotreal3(vel,lbmaux[0].C[k]);
        real FdotC = dotreal3(Fm ,lbmaux[0].C[k]);
        //Navier Stokes equation
        real sk    = 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs);
        real Feq   = lbmaux[0].W[k]*(3.0*pre/(Cs*Cs) + rho*sk);
        real Fk    = lbmaux[0].W[k]*(S + 3.0*dotreal3(lbmaux[0].C[k],BForce[ic] + BForce[ic + 1*lbmaux[0].Ncells])/(Cs*Cs) + 4.5*VdotC*FdotC/(Cs*Cs) - 1.5*VdotF/(Cs*Cs));
        if (k==0)  
        {
            Feq    += -3.0*pre/(Cs*Cs);
            BForce[ic + 2*lbmaux[0].Ncells].x = 0.5*lbmaux[0].dt*S + tau*lbmaux[0].dt*Fk + rho*lbmaux[0].W[k]*sk;
        }
        size_t idx = ic*lbmaux[0].Nneigh + k;
        Ftemp[idx] = F[idx] - (F[idx]-Feq)/tau + (1.0-0.5/tau)*lbmaux[0].dt*Fk;
        if (IsSolid[ic]) Ftemp[idx] = F[ic*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 

        //Phase field equation
        idx  = ic*lbmaux[0].Nneigh + 1*lbmaux[0].Ncells*lbmaux[0].Nneigh + k;
        Feq   = lbmaux[0].W[k]*phi*(1.0 + 3.0*VdotC/Cs);
        real G     = 3.0*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[k],dFdt+Cs*Cs/3.0*lambda*n)/(Cs*Cs)
                     + lbmaux[0].W[k]*phi*divu;
        if (k==0) BForce[ic + 2*lbmaux[0].Ncells].y = 0.5*lbmaux[0].dt*phi*divu;
        Ftemp[idx] = F[idx] - (F[idx]-Feq)/taup + (1.0-0.5/taup)*lbmaux[0].dt*G;

        // Enthalpy equation
        idx   = ic*lbmaux[0].Nneigh + 2*lbmaux[0].Ncells*lbmaux[0].Nneigh + k;
        Feq   = lbmaux[0].W[k]*Cp*Temp*(1.0 + sk);
        if (k==0)  Feq += H - Cp*Temp;
        Ftemp[idx] = F[idx] - (F[idx]-Feq)/tauh;
    }


    //update for next time step
    Vel   [ic+1*lbmaux[0].Ncells]   = phi*vel;
    BForce[ic+2*lbmaux[0].Ncells].z = rho*divu + dotreal3(vel,gradrho);
}

__global__ void cudaApplyForcesSC(uint3 * pCellPairs, bool const * IsSolid, real3 * BForce, real const * Rho, lbm_aux const * lbmaux)
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
        if (fabs(G)<1.0e-12) continue;
        if (!IsSolid[ic]) psic = lbmaux[0].Psi[il]*exp(-lbmaux[0].Rhoref[il]/Rho[ic+il*lbmaux[0].Ncells]);
        if (!IsSolid[in]) psin = lbmaux[0].Psi[il]*exp(-lbmaux[0].Rhoref[il]/Rho[in+il*lbmaux[0].Ncells]);
        else              G    = lbmaux[0].Gs[il];
        
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

__global__ void cudaApplyForcesSCMP(uint3 * pCellPairs, bool const * IsSolid, real3 * BForce, real const * Rho, lbm_aux const * lbmaux)
{
    size_t icp = threadIdx.x + blockIdx.x * blockDim.x;
    if (icp>=lbmaux[0].NCPairs) return;
    size_t ic = pCellPairs[icp].x;
    size_t in = pCellPairs[icp].y;
    size_t k  = pCellPairs[icp].z;

    for (size_t il=0;il<lbmaux[0].Nl;il++)
    {
        real psic = 0.0;
        real psin = 0.0;
        real G    = lbmaux[0].G[il];
        if (fabs(G)<1.0e-12) continue;
        if (!IsSolid[ic+il*lbmaux[0].Ncells]) psic = lbmaux[0].Psi[il]*exp(-lbmaux[0].Rhoref[il]/Rho[ic+il*lbmaux[0].Ncells]);
        if (!IsSolid[in+il*lbmaux[0].Ncells]) psin = lbmaux[0].Psi[il]*exp(-lbmaux[0].Rhoref[il]/Rho[in+il*lbmaux[0].Ncells]);
        else              G    = lbmaux[0].Gs[il];

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
        real psic = 1.0;
        real psin = 1.0;
        real G    = lbmaux[0].Gmix;
        if (!IsSolid[ic+il1*lbmaux[0].Ncells]) psic = Rho[ic+il1*lbmaux[0].Ncells];
        else              G    = lbmaux[0].Gs[il2];
        if (!IsSolid[in+il2*lbmaux[0].Ncells]) psin = Rho[in+il2*lbmaux[0].Ncells];
        else              G    = lbmaux[0].Gs[il1];

        real3 bforce = (-G*lbmaux[0].W[k]*psic*psin)*lbmaux[0].C[k];

        atomicAdd(&BForce[ic+il1*lbmaux[0].Ncells].x, bforce.x);
        atomicAdd(&BForce[ic+il1*lbmaux[0].Ncells].y, bforce.y);
        atomicAdd(&BForce[ic+il1*lbmaux[0].Ncells].z, bforce.z);
        atomicAdd(&BForce[in+il2*lbmaux[0].Ncells].x,-bforce.x);
        atomicAdd(&BForce[in+il2*lbmaux[0].Ncells].y,-bforce.y);
        atomicAdd(&BForce[in+il2*lbmaux[0].Ncells].z,-bforce.z);

        psic = 1.0;
        psin = 1.0;
        G    = lbmaux[0].Gmix;
        if (!IsSolid[ic+il2*lbmaux[0].Ncells]) psic = Rho[ic+il2*lbmaux[0].Ncells];
        else              G    = lbmaux[0].Gs[il1];
        if (!IsSolid[in+il1*lbmaux[0].Ncells]) psin = Rho[in+il1*lbmaux[0].Ncells];
        else              G    = lbmaux[0].Gs[il2];

        bforce = (-G*lbmaux[0].W[k]*psic*psin)*lbmaux[0].C[k];

        atomicAdd(&BForce[ic+il2*lbmaux[0].Ncells].x, bforce.x);
        atomicAdd(&BForce[ic+il2*lbmaux[0].Ncells].y, bforce.y);
        atomicAdd(&BForce[ic+il2*lbmaux[0].Ncells].z, bforce.z);
        atomicAdd(&BForce[in+il1*lbmaux[0].Ncells].x,-bforce.x);
        atomicAdd(&BForce[in+il1*lbmaux[0].Ncells].y,-bforce.y);
        atomicAdd(&BForce[in+il1*lbmaux[0].Ncells].z,-bforce.z);
    }
    
}

__global__ void cudaStream1(real * F, real * Ftemp, real3 * BForce, lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;
#ifdef USE_IBB
    if (ic==0)
    {
        lbmaux[0].Time += lbmaux[0].dt;
        lbmaux[0].iter++;
    }
    BForce[ic] = make_real3(0.0,0.0,0.0);
#endif
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ic/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;
    size_t icl = ic/lbmaux[0].Ncells;

    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        size_t inx = (size_t)((int)icx + (int)lbmaux[0].C[k].x + (int)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((int)icy + (int)lbmaux[0].C[k].y + (int)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((int)icz + (int)lbmaux[0].C[k].z + (int)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny + icl*lbmaux[0].Ncells;
        Ftemp[in*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k];
        //if (ic==lbmaux[0].Ncells/2+lbmaux[0].Nx/2-1) printf("%g %lu %lu \n",F[ic*lbmaux[0].Nneigh + 0*lbmaux[0].Ncells + k],in,k);
    }
}

__global__ void cudaStream2(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;
    if (ic==0)
    {
        lbmaux[0].Time += lbmaux[0].dt;
        lbmaux[0].iter++;
    }
    BForce[ic] = make_real3(0.0,0.0,0.0);
    //for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    //{
        //F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k];
    //}
    Rho   [ic] = 0.0;
    Vel   [ic] = make_real3(0.0,0.0,0.0);
    if (!IsSolid[ic])
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ic] += F[ic*lbmaux[0].Nneigh + k];
            Vel[ic] = Vel[ic] + F[ic*lbmaux[0].Nneigh + k]*lbmaux[0].C[k];
            //if (ic==0) printf("k: %lu %f %f %f %f \n",k,Rho[ic],Vel[ic].x,Vel[ic].y,Vel[ic].z);
        }
        Vel[ic] = lbmaux[0].Cs/Rho[ic]*Vel[ic];
    }
}

__global__ void cudaStreamPFI2(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;
    if (ic==0)
    {
        lbmaux[0].Time += lbmaux[0].dt;
        lbmaux[0].iter++;
    }
    Rho   [ic                   ] = BForce[ic+2*lbmaux[0].Ncells].x;
    Rho   [ic+1*lbmaux[0].Ncells] = BForce[ic+2*lbmaux[0].Ncells].y;
    Rho   [ic+2*lbmaux[0].Ncells] = 0.0;
    Vel   [ic] = make_real3(0.0,0.0,0.0);
    //if (ic==lbmaux[0].Ncells/2+lbmaux[0].Nx/2) printf("%g %g %g %lu \n",Rho[ic],Rho[ic+1*lbmaux[0].Ncells],Rho[ic+2*lbmaux[0].Ncells],lbmaux[0].iter);
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        if (k!=0) Rho[ic                     ] += F[ic*lbmaux[0].Nneigh                                       + k];
                  Rho[ic + 1*lbmaux[0].Ncells] += F[ic*lbmaux[0].Nneigh + 1*lbmaux[0].Nneigh*lbmaux[0].Ncells + k];
                  Rho[ic + 2*lbmaux[0].Ncells] += F[ic*lbmaux[0].Nneigh + 2*lbmaux[0].Nneigh*lbmaux[0].Ncells + k];
        
        //if (ic==217*lbmaux[0].Nx+0) printf("%g %g %g %lu \n",F[ic*lbmaux[0].Nneigh + 0*lbmaux[0].Ncells + k],Ftemp[ic*lbmaux[0].Nneigh + 0*lbmaux[0].Ncells + k],Rho[ic+0*lbmaux[0].Ncells],k);
        Vel[ic] = Vel[ic] + F[ic*lbmaux[0].Nneigh + k]*lbmaux[0].C[k];
    }

    Rho[ic] *= lbmaux[0].Cs*lbmaux[0].Cs/3.0/(1-lbmaux[0].W[0]);

    //if (ic==216*lbmaux[0].Nx+199) printf("%g %g %g %lu \n",Rho[ic],Rho[ic+1*lbmaux[0].Ncells],Rho[ic+2*lbmaux[0].Ncells],lbmaux[0].iter);

    real H    = Rho[ic + 2*lbmaux[0].Ncells];
    real Hs   = lbmaux[0].Ts*lbmaux[0].cap[0];
    real Hl   = lbmaux[0].Tl*lbmaux[0].cap[1]+lbmaux[0].L;
    real fl   = 0.0;
    if      (H>=Hs&&H<=Hl) fl = (H-Hs)/(Hl-Hs);
    else if (H>Hl)         fl = 1.0;
    Vel[ic + 2*lbmaux[0].Ncells].y = Vel[ic + 2*lbmaux[0].Ncells].x;
    Vel[ic + 2*lbmaux[0].Ncells].x = fl;
    real fs = 1.0 - fl;

    real phi  = Rho[ic + 1*lbmaux[0].Ncells];
    real rho  = phi*(fs*lbmaux[0].rho[0] + fl*lbmaux[0].rho[1])+fl*(1.0-phi)*lbmaux[0].rho[2];
    BForce[ic + 1*lbmaux[0].Ncells] = (lbmaux[0].Cs/rho)*(Vel[ic] + 0.5*lbmaux[0].dt*BForce[ic]);
    Vel[ic] = (1.0-0.5*fs)*BForce[ic + 1*lbmaux[0].Ncells];

    Vel[ic+2*lbmaux[0].Ncells].z = BForce[ic+2*lbmaux[0].Ncells].z;
}
}
#endif //MECHSYS_LBM_CUH
