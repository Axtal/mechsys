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
    //These parameters are for the SW solver;
    real       g;           ///< Gravity acceleration
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

__device__ __inline__ real FeqFluid(size_t const & k, real const & rho, real3 const & vel, lbm_aux const * lbmaux)
{
    real VdotV = dotreal3(vel,vel);
    real VdotC = dotreal3(vel,lbmaux[0].C[k]);
    real Cs    = lbmaux[0].Cs;
    return lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

__device__ __inline__ real FeqSW   (size_t const & k, real const & h  , real3 const & vel, lbm_aux const * lbmaux)
{    
    real Cs    = lbmaux[0].Cs;
    real VdotV = dotreal3(vel,vel)/(Cs*Cs);
    real VdotC = dotreal3(vel,lbmaux[0].C[k])/Cs;
    if (k==0)
    {
        return h - 5.0/6.0*lbmaux[0].g*h*h/(Cs*Cs) - 2.0/3.0*h*VdotV;
    }
    else
    {
        return lbmaux[0].W[k]*h*(1.5*lbmaux[0].g*h/(Cs*Cs) + 3.0*VdotC + 4.5*VdotC*VdotC - 1.5*VdotV);
    }
}

// ---------------------------------------------------------------------------
// Layout of the distribution-function buffers (F and Ftemp).
//
// Each buffer holds Nl*Nneigh planes of Ncells values.  Population k of cell ic
// (a *local* cell index, 0 <= ic < Ncells) of lattice il lives at
//
//     F[(il*Nneigh + k)*Ncells + ic]
//
// that is, the direction index varies slowest: a structure of arrays.
//
// This used to be the other way round, F[ic*Nneigh + k], which keeps the
// populations of one cell contiguous.  That is pleasant on a CPU but bad on a
// GPU: for a fixed direction k, consecutive threads then read addresses
// Nneigh*8 bytes apart, so a 32-thread warp touches 32 different 128 byte
// sectors (D3Q15: 32*120 = 3840 bytes) in order to use 8 bytes out of each.
// With the direction-major layout a warp reading a fixed k is perfectly
// contiguous, and one sector serves 16 consecutive threads instead of one.
//
// The DRAM traffic is identical either way -- every value is still read exactly
// once -- so this is purely about transaction/sector efficiency.  That is why
// the win is large even though these kernels are nowhere near bandwidth bound.
// Measured on this machine's RTX A6000 for the fused streaming kernel at 192^3
// D3Q15 in double precision: 25.45 -> 4.19 ms/step, 6.07x, with bit-for-bit
// identical results.  Padding the per-cell block to a power of two (16 for
// D3Q15) changed nothing, which confirms the cause is coalescing rather than
// address arithmetic.
//
// On-disk files (hdf5/xmf) still store the populations cell-contiguously: the
// conversion happens in Domain::InitDevice and Domain::Save, and host code that
// uses the public F[il][nx][ny][nz][k] arrays is unaffected.
// ---------------------------------------------------------------------------
__host__ __device__ __forceinline__ size_t FIDX (size_t Ncells, size_t Nneigh, size_t il, size_t k, size_t ic)
{
    return (il*Nneigh + k)*Ncells + ic;
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

// Keeps the device-side time/iteration counters of lbm_aux up to date.  These
// used to live in a branch of cudaStream2 / cudaStreamPF2 / ...; splitting them
// out lets the hot fused kernel take lbm_aux as const __restrict__.
__global__ void cudaTick (lbm_aux * lbmaux)
{
    lbmaux[0].Time += lbmaux[0].dt;
    lbmaux[0].iter++;
}

// ---------------------------------------------------------------------------
// Fused streaming + macroscopic update.
//
// The original per-step sequence for the Navier-Stokes / Advection-Diffusion
// solvers was
//     cudaStream1(...)   // scatter: Ftemp[in][k] = F[ic][k], in = ic + C[k]
//     swap(F,Ftemp)
//     cudaStream2(...)   // re-reads F: Rho, Vel, reset BForce
// i.e. the whole distribution array was streamed twice, and the second pass
// re-read it purely to sum the populations.
//
// Streaming is a pure lattice translation, so the scatter is exactly the gather
//     F[ic][k] = Fpost[ic - C[k]][k]
// Writing it as a gather lets one thread own one destination cell, so the
// translation, the Rho/Vel accumulation and the BForce reset all share a single
// pass over memory.  This removes one full read plus one full write of the
// distribution array per time step.
//
// The F write is unconditional (solid cells take part in the stream, carrying
// the bounce-back populations written by the collide kernel); only the Rho/Vel
// store is restricted to fluid cells, exactly as cudaStream2 did.
// ---------------------------------------------------------------------------
template <size_t NN>
__global__ void cudaFusedStreamMacro (bool const * __restrict__ IsSolid, real const * __restrict__ Fpost, real * __restrict__ F,
        real3 * __restrict__ BForce, real3 * __restrict__ Vel, real * __restrict__ Rho, lbm_aux const * __restrict__ lbmaux)
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

    const size_t ic0 = ic - icl*Ncells;   // local cell index within the lattice
    real  rho = 0.0;
    real3 vel = make_real3(0.0,0.0,0.0);

    #pragma unroll
    for (size_t k=0;k<Nneigh;k++)
    {
        // Periodic wrap.  (icx - Cx) lies in (-Nx, 2Nx) for every supported
        // lattice, so at most one adjustment is needed per axis.  An integer
        // modulo by a run-time divisor costs ~20 instructions, and this loop
        // runs three of them per discrete velocity, i.e. 45 per cell.
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
    }
}

// cudaCollideSC<NN>
//
// NN : number of discrete velocities, known at compile time so that the
//      population loops are fully unrolled and NonEq[] is kept in registers
//      instead of being spilled to local memory (the loop bound used to be the
//      run-time value lbmaux[0].Nneigh, which forced a 216 byte local frame).
//      NN==0 keeps the previous run-time bounded behaviour.
template <size_t NN>
__global__ void cudaCollideSC(bool const * __restrict__ IsSolid, real * __restrict__ F, real * __restrict__ Ftemp,
        real3 const * __restrict__ BForce, real3 const * __restrict__ Vel, real const * __restrict__ Rho,
        lbm_aux const * __restrict__ lbmaux)
{
    const size_t Nneigh = (NN>0) ? NN : lbmaux[0].Nneigh;
    const size_t Ncells = lbmaux[0].Ncells;
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;

    if (!IsSolid[ic])
    {
        const real  dt   = lbmaux[0].dt;
        const real  Cs   = lbmaux[0].Cs;
        const real  Cs2  = Cs*Cs;
        const real  tau0 = lbmaux[0].Tau[0];
        const real  Sc   = lbmaux[0].Sc;
        const real  rho  = Rho[ic];

        real3 vel = Vel[ic]+dt*(tau0/rho)*BForce[ic];
        real  tau = tau0;

        real  NonEq[(NN>0)?NN:27];
        real  Q = 0.0;

        // Loop invariants of the equilibrium distribution.  Fc is exactly the
        // sub-expression the compiler used to common out of FeqFluid, so naming
        // it costs nothing and keeps the arithmetic identical.
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
            // Same as FeqFluid(k,rho,vel,lbmaux) with the invariants hoisted.
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

        bool valid = true;
        real alpha = 1.0;
        size_t numit = 0;
        while (valid&&numit<2)
        {
            valid = false;
            #pragma unroll
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[FIDX(Ncells,Nneigh,0,k,ic)] = F[FIDX(Ncells,Nneigh,0,k,ic)] - alpha*(NonEq[k]*itau);
                if (Ftemp[FIDX(Ncells,Nneigh,0,k,ic)]<0.0)
                {
                    real temp = tau*F[FIDX(Ncells,Nneigh,0,k,ic)]/(NonEq[k]);
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

// cudaCollideMP<NN> -- multi-component / multiphase collision.
// Same treatment as cudaCollideSC: compile-time Nneigh, fully unrolled loops,
// __restrict__ pointers and hoisted loop invariants.
template <size_t NN>
__global__ void cudaCollideMP(bool const * __restrict__ IsSolid, real * __restrict__ F, real * __restrict__ Ftemp,
        real3 const * __restrict__ BForce, real3 const * __restrict__ Vel, real const * __restrict__ Rho,
        lbm_aux const * __restrict__ lbmaux)
{
    const size_t Nneigh = (NN>0) ? NN : lbmaux[0].Nneigh;
    const real   Cs0    = lbmaux[0].Cs;
    const size_t Ncells = lbmaux[0].Ncells;
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=Ncells) return;

    real3 Vmix = make_real3(0.0,0.0,0.0);
    real  den  = 0.0;
    for(size_t il=0;il<lbmaux[0].Nl;il++)
    {
        Vmix = Vmix + (Rho[ic+il*Ncells]/lbmaux[0].Tau[il])*Vel[ic+il*Ncells];
        den  = den  + Rho[ic+il*Ncells]/lbmaux[0].Tau[il];
    }
    Vmix = Vmix/den;

    for(size_t il=0;il<lbmaux[0].Nl;il++)
    {
        if (!IsSolid[ic+il*Ncells])
        {
            real  rho   = Rho[ic+il*Ncells];
            real3 vel   = Vmix + (lbmaux[0].dt*lbmaux[0].Tau[il]/rho)*BForce[ic+il*Ncells];
            real  VdotV = dotreal3(vel,vel);
            real  tau   = lbmaux[0].Tau[il];
            // tau is fixed for this cell, so a single reciprocal replaces the Nneigh
            // divisions below (`x/tau` -> `x*itau`; each quotient may move <=1 ulp).
            const real itau = 1.0/tau;
            const real Fc    = 1.5*VdotV/(Cs0*Cs0);
            const real iCs   = 1.0/Cs0;
            const real iCs2  = 1.0/(Cs0*Cs0);
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
                    real Feq   = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC*iCs + 4.5*VdotC*VdotC*iCs2 - Fc);
                    size_t idx = FIDX(Ncells,Nneigh,il,k,ic);
                    Ftemp[idx] = F[idx] - alphal*(F[idx]-Feq)*itau;
                    if (Ftemp[idx]<0.0&&numit<2)
                    {
                        real temp = tau*fabs(F[idx]/(F[idx]-Feq));
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

__global__ void cudaCollideSC_MRT(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux const * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;

    // In the direction-major layout the populations of a cell are no longer
    // contiguous, so gather them into a local frame for the matrix products.
    const size_t NcellsM = lbmaux[0].Ncells;
    const size_t NnM     = lbmaux[0].Nneigh;
    const size_t iclM    = ic/NcellsM;
    const size_t ic0M    = ic - iclM*NcellsM;

    if (!IsSolid[ic])
    {
        real  Cs  = lbmaux[0].Cs;
        real  dx  = lbmaux[0].dx;
        real  dt  = lbmaux[0].dt;
        real  rho = Rho[ic];
        real3 vel = Vel[ic];

        real f   [27];
        real ft  [27];
        real fneq[27];
        for (size_t k=0;k<NnM;k++) f[k] = F[FIDX(NcellsM,NnM,iclM,k,ic0M)];

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
                Ftemp[FIDX(NcellsM,NnM,iclM,k,ic0M)] = F[FIDX(NcellsM,NnM,iclM,k,ic0M)] - alpha*(fneq[k] - 3.0*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[k],BForce[ic])*dt*dt/dx);
                if (Ftemp[FIDX(NcellsM,NnM,iclM,k,ic0M)]<0.0)
                {
                    real temp = F[FIDX(NcellsM,NnM,iclM,k,ic0M)]/(fneq[k] - 3.0*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[k],BForce[ic])*dt*dt/dx);
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
            Ftemp[FIDX(NcellsM,NnM,iclM,k,ic0M)] = F[FIDX(NcellsM,NnM,iclM,lbmaux[0].Op[k],ic0M)]; 
        }
    }
}

// cudaCollideAD<NN> -- advection-diffusion collision.
template <size_t NN>
__global__ void cudaCollideAD(bool const * __restrict__ IsSolid, real * __restrict__ F, real * __restrict__ Ftemp,
        real3 const * __restrict__ BForce, real3 const * __restrict__ Vel, real const * __restrict__ Rho,
        lbm_aux const * __restrict__ lbmaux)
{
    const size_t Nneigh = (NN>0) ? NN : lbmaux[0].Nneigh;
    const real   Cs0    = lbmaux[0].Cs;
    const real   Cs2    = Cs0*Cs0;
    const size_t Ncells = lbmaux[0].Ncells;
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=Ncells) return;

    for(size_t il=0;il<lbmaux[0].Nl;il++)
    {
        if (!IsSolid[ic+il*Ncells])
        {
            real  rho   = Rho[ic+il*Ncells];
            real3 vel   = Vel[ic] + (lbmaux[0].dt*lbmaux[0].Tau[il]/Rho[ic])*BForce[ic+il*Ncells]; // ic for the velocity because for the AD
                                                                                                   // solver only velocity of phase 0 is
                                                                                                   // important
            real  VdotV = dotreal3(vel,vel);
            real  tau   = lbmaux[0].Tau[il];
            // tau is fixed for this cell, so a single reciprocal replaces the Nneigh
            // divisions below (`x/tau` -> `x*itau`; each quotient may move <=1 ulp).
            const real itau = 1.0/tau;
            const real Fc = 1.5*VdotV/Cs2;
            const real iCs   = 1.0/Cs0;
            const real iCs2  = 1.0/(Cs0*Cs0);
            bool valid = true;
            real alpha = 1.0;
            size_t numit = 0;
            while (valid&&numit<2)
            {
                valid = false;
                #pragma unroll
                for (size_t k=0;k<Nneigh;k++)
                {
                    real VdotC = dotreal3(vel,lbmaux[0].C[k]);
                    real Feq   = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC*iCs + 4.5*VdotC*VdotC*iCs2 - Fc);
                    size_t idx = FIDX(Ncells,Nneigh,il,k,ic);
                    Ftemp[idx] = F[idx] - alpha*(F[idx]-Feq)*itau;
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
            #pragma unroll
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[FIDX(Ncells,Nneigh,il,k,ic)] = F[FIDX(Ncells,Nneigh,il,lbmaux[0].Op[k],ic)];
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
        // Periodic wrap, single conditional adjustment per axis instead of an
        // integer modulo by a run-time divisor (see cudaFusedStreamMacro).
        int inx = (int)icx + (int)lbmaux[0].C[k].x;
        int iny = (int)icy + (int)lbmaux[0].C[k].y;
        int inz = (int)icz + (int)lbmaux[0].C[k].z;
        if (inx<0) inx += (int)lbmaux[0].Nx; else if (inx>=(int)lbmaux[0].Nx) inx -= (int)lbmaux[0].Nx;
        if (iny<0) iny += (int)lbmaux[0].Ny; else if (iny>=(int)lbmaux[0].Ny) iny -= (int)lbmaux[0].Ny;
        if (inz<0) inz += (int)lbmaux[0].Nz; else if (inz>=(int)lbmaux[0].Nz) inz -= (int)lbmaux[0].Nz;
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
    // tau is fixed for this cell, so a single reciprocal replaces the Nneigh
    // divisions below (`x/tau` -> `x*itau`; each quotient may move <=1 ulp).
    const real itau = 1.0/tau;
    real  taup  = lbmaux[0].Tau[1];
    real  kappa = phi*(fs*lbmaux[0].kap[0] + fl*lbmaux[0].kap[1])+fl*(1.0-phi)*lbmaux[0].kap[2];
    real  tauh  = 3.0*kappa*lbmaux[0].dt/(lbmaux[0].dx*lbmaux[0].dx) + 0.5;
    real  VdotV                      = dotreal3(vel,vel);
    real  mu                         = 4.0*b*(phi)*(phi-1.0)*(phi-0.5) - a*delphi;
    BForce[ic]                       = BForce[ic] + mu*gradphi;
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
        real Fk    = lbmaux[0].W[k]*(S + 3.0*dotreal3(lbmaux[0].C[k],BForce[ic] + BForce[ic + 1*lbmaux[0].Ncells])/Cs + 4.5*VdotC*FdotC/(Cs*Cs) - 1.5*VdotF/(Cs*Cs));
        if (k==0)  
        {
            Feq    += -3.0*pre/(Cs*Cs);
            BForce[ic + 2*lbmaux[0].Ncells].x = 0.5*lbmaux[0].dt*S + tau*lbmaux[0].dt*Fk + rho*lbmaux[0].W[k]*sk;
        }
        size_t idx = FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic);
        Ftemp[idx] = F[idx] - (F[idx]-Feq)*itau + (1.0-0.5*itau)*lbmaux[0].dt*Fk;
        if (IsSolid[ic]) Ftemp[idx] = F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,lbmaux[0].Op[k],ic)]; 

        //Phase field equation
        idx  = FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,1,k,ic);
        Feq   = lbmaux[0].W[k]*phi*(1.0 + 3.0*VdotC/Cs);
        real G     = 3.0*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[k],dFdt+Cs*Cs/3.0*lambda*n)/(Cs*Cs)
                     + lbmaux[0].W[k]*phi*divu;
        if (k==0) BForce[ic + 2*lbmaux[0].Ncells].y = 0.5*lbmaux[0].dt*phi*divu;
        Ftemp[idx] = F[idx] - (F[idx]-Feq)/taup + (1.0-0.5/taup)*lbmaux[0].dt*G;

        // Enthalpy equation
        idx   = FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,2,k,ic);
        Feq   = lbmaux[0].W[k]*Cp*Temp*(1.0 + sk);
        if (k==0)  Feq += H - Cp*Temp;
        Ftemp[idx] = F[idx] - (F[idx]-Feq)/tauh;
    }


    //update for next time step. Care must be taken if vel is not zero at the beginning
    Vel   [ic+1*lbmaux[0].Ncells]   = phi*vel;
    BForce[ic+2*lbmaux[0].Ncells].z = rho*divu + dotreal3(vel,gradrho);
}

__global__ void cudaCollidePF(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux const * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;

    real Cs    = lbmaux[0].Cs;
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ic/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;

    //Phase field variables
    real3 vel   = Vel[ic]; // it is ic since only the velocity of the first layer is relavant
    real  phi   = Rho[ic+1*lbmaux[0].Ncells];
    real3 dFdt  = (phi*vel - Vel[ic+1*lbmaux[0].Ncells])/lbmaux[0].dt; //Here we will use the Vel of lattice one as a temporal variable for the
                                                                       //phase field flux
    real lambda = 4.0*phi*(1.0-phi)/lbmaux[0].thick;
    real a      = 1.5*lbmaux[0].sigma*lbmaux[0].thick;
    real b      = 12.0*lbmaux[0].sigma/lbmaux[0].thick;


    //Density variables
    real rho  = phi*lbmaux[0].rho[0] + (1.0-phi)*lbmaux[0].rho[1];
    real pre  = Rho[ic];
    real S    = BForce[ic+1*lbmaux[0].Ncells].x;
    
    //Gradients
    real3 gradphi = make_real3(0.0,0.0,0.0);
    real3 gradrho = make_real3(0.0,0.0,0.0);
    real3 gradpre = make_real3(0.0,0.0,0.0);
    real3 gradS   = make_real3(0.0,0.0,0.0);
    real  delphi  = 0.0;

    for (size_t k=1;k<lbmaux[0].Nneigh;k++)
    {
        // Periodic wrap, single conditional adjustment per axis instead of an
        // integer modulo by a run-time divisor (see cudaFusedStreamMacro).
        int inx = (int)icx + (int)lbmaux[0].C[k].x;
        int iny = (int)icy + (int)lbmaux[0].C[k].y;
        int inz = (int)icz + (int)lbmaux[0].C[k].z;
        if (inx<0) inx += (int)lbmaux[0].Nx; else if (inx>=(int)lbmaux[0].Nx) inx -= (int)lbmaux[0].Nx;
        if (iny<0) iny += (int)lbmaux[0].Ny; else if (iny>=(int)lbmaux[0].Ny) iny -= (int)lbmaux[0].Ny;
        if (inz<0) inz += (int)lbmaux[0].Nz; else if (inz>=(int)lbmaux[0].Nz) inz -= (int)lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;

        real pren = Rho[in];
        real phin = Rho[in+1*lbmaux[0].Ncells];
        real rhon = phin*lbmaux[0].rho[0]+(1.0-phin)*lbmaux[0].rho[1];
        real Sn   = BForce[in+1*lbmaux[0].Ncells].x;
        
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
    
    real  tau   = lbmaux[0].Tau[0];
    // tau is fixed for this cell, so a single reciprocal replaces the Nneigh
    // divisions below (`x/tau` -> `x*itau`; each quotient may move <=1 ulp).
    const real itau = 1.0/tau;
    real  taup  = lbmaux[0].Tau[1];
    real  VdotV                      = dotreal3(vel,vel);
    real  mu                         = 4.0*b*(phi)*(phi-1.0)*(phi-0.5) - a*delphi;
    BForce[ic]                       = BForce[ic] + mu*gradphi;
    real3 Fm                         = BForce[ic] + Cs*Cs*gradrho/3.0 + Cs*Cs*gradS/3.0 - gradpre;
    real  VdotF = dotreal3(vel,Fm);
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        real VdotC = dotreal3(vel,lbmaux[0].C[k]);
        real FdotC = dotreal3(Fm ,lbmaux[0].C[k]);
        //Navier Stokes equation
        real sk    = 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs);
        real Feq   = lbmaux[0].W[k]*(3.0*pre/(Cs*Cs) + rho*sk);
        real Fk    = lbmaux[0].W[k]*(S + 3.0*dotreal3(lbmaux[0].C[k],BForce[ic])/Cs + 4.5*VdotC*FdotC/(Cs*Cs) - 1.5*VdotF/(Cs*Cs));
        if (k==0)  
        {
            Feq    += -3.0*pre/(Cs*Cs);
            BForce[ic + 1*lbmaux[0].Ncells].z = 0.5*lbmaux[0].dt*S + tau*lbmaux[0].dt*Fk + rho*lbmaux[0].W[k]*sk;
        }
        size_t idx = FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic);
        Ftemp[idx] = F[idx] - (F[idx]-Feq)*itau + (1.0-0.5*itau)*lbmaux[0].dt*Fk;
        if (IsSolid[ic]) Ftemp[idx] = F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,lbmaux[0].Op[k],ic)]; 

        //Phase field equation
        idx  = FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,1,k,ic);
        Feq   = lbmaux[0].W[k]*phi*(1.0 + 3.0*VdotC/Cs);
        real G     = 3.0*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[k],dFdt+Cs*Cs/3.0*lambda*n)/(Cs*Cs);
        Ftemp[idx] = F[idx] - (F[idx]-Feq)/taup + (1.0-0.5/taup)*lbmaux[0].dt*G;
    }


    //update for next time step
    Vel   [ic+1*lbmaux[0].Ncells]   = phi*vel;
    BForce[ic+1*lbmaux[0].Ncells].y = dotreal3(vel,gradrho);
}

__global__ void cudaCollideSW(real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux const * lbmaux)
{   
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;

    const size_t Ncells = lbmaux[0].Ncells;
    const size_t Nn     = lbmaux[0].Nneigh;
    const size_t icl    = ic/Ncells;
    const size_t ic0    = ic - icl*Ncells;

    real3 vel   = Vel[ic];
    real  h     = Rho[ic];
    real  tau   = lbmaux[0].Tau[0];
    
    real  NonEq[27];
    real  Q = 0.0;
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        NonEq[k]     = F[FIDX(Ncells,Nn,icl,k,ic0)] - FeqSW(k,h,vel,lbmaux);
        Q           += NonEq[k]*NonEq[k]*lbmaux[0].EEk[k];
    }
    Q = sqrt(Q);
    tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*lbmaux[0].Sc*lbmaux[0].Sc/h));
    // tau is fixed for this cell, so a single reciprocal replaces the
    // Nneigh divisions below.
    const real itau = 1.0/tau;

    real3 Force = -lbmaux[0].g*h*BForce[ic] - lbmaux[0].cap[0]*norm(vel)*vel;

    bool valid = true;
    real alpha = 1.0;
    size_t numit = 0;
    while (valid&&numit<2)
    {
        valid = false;
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Ftemp[FIDX(Ncells,Nn,icl,k,ic0)] = F[FIDX(Ncells,Nn,icl,k,ic0)] - alpha*(NonEq[k]*itau -
                    lbmaux[0].dt*dotreal3(lbmaux[0].C[k],Force)/(6.0*lbmaux[0].Cs));
            if (Ftemp[FIDX(Ncells,Nn,icl,k,ic0)]<0.0)
            {
                real temp = tau*F[FIDX(Ncells,Nn,icl,k,ic0)]/(NonEq[k] - lbmaux[0].dt*dotreal3(lbmaux[0].C[k],Force)/(6.0*lbmaux[0].Cs));
                if (temp<alpha) alpha = temp;
                valid = true;
            }
        }
        if (valid) numit++;
    }
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
        else              G    = lbmaux[0].Gs[il];
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
        else              G    = lbmaux[0].Gs[il];
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

    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = (ic/(lbmaux[0].Nx*lbmaux[0].Ny))%lbmaux[0].Nz;
    size_t icl = ic/lbmaux[0].Ncells;
    size_t ic0 = ic - icl*lbmaux[0].Ncells;

    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        // Periodic wrap, single conditional adjustment per axis instead of an
        // integer modulo by a run-time divisor (see cudaFusedStreamMacro).
        int inx = (int)icx + (int)lbmaux[0].C[k].x;
        int iny = (int)icy + (int)lbmaux[0].C[k].y;
        int inz = (int)icz + (int)lbmaux[0].C[k].z;
        if (inx<0) inx += (int)lbmaux[0].Nx; else if (inx>=(int)lbmaux[0].Nx) inx -= (int)lbmaux[0].Nx;
        if (iny<0) iny += (int)lbmaux[0].Ny; else if (iny>=(int)lbmaux[0].Ny) iny -= (int)lbmaux[0].Ny;
        if (inz<0) inz += (int)lbmaux[0].Nz; else if (inz>=(int)lbmaux[0].Nz) inz -= (int)lbmaux[0].Nz;
        size_t in0 = (size_t)inx + (size_t)iny*lbmaux[0].Nx + (size_t)inz*lbmaux[0].Nx*lbmaux[0].Ny;
        Ftemp[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,in0)] = F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)];
        //if (ic==lbmaux[0].Ncells/2+lbmaux[0].Nx/2-1) printf("%g %lu %lu \n",F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,in0)],in0,k);
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
    const size_t icl = ic/lbmaux[0].Ncells;
    const size_t ic0 = ic - icl*lbmaux[0].Ncells;
    BForce[ic] = make_real3(0.0,0.0,0.0);
    //for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    //{
        //F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)] = Ftemp[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)];
    //}
    Rho   [ic] = 0.0;
    Vel   [ic] = make_real3(0.0,0.0,0.0);
    if (!IsSolid[ic])
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ic] += F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)];
            Vel[ic] = Vel[ic] + F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)]*lbmaux[0].C[k];
            //if (ic==0) printf("k: %lu %f %f %f %f \n",k,Rho[ic],Vel[ic].x,Vel[ic].y,Vel[ic].z);
        }
        Vel[ic] = lbmaux[0].Cs/Rho[ic]*Vel[ic];
    }
}

__global__ void cudaStreamPF2(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ncells) return;
    if (ic==0)
    {
        lbmaux[0].Time += lbmaux[0].dt;
        lbmaux[0].iter++;
    }
    Rho   [ic                   ] = BForce[ic+1*lbmaux[0].Ncells].z;
    Rho   [ic+1*lbmaux[0].Ncells] = 0.0;
    Vel   [ic] = make_real3(0.0,0.0,0.0);
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        if (k!=0) Rho[ic                     ] += F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic)];
                  Rho[ic + 1*lbmaux[0].Ncells] += F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,1,k,ic)];
        
        //if (ic==217*lbmaux[0].Nx+0) printf("%g %g %g %lu \n",F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic)],Ftemp[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic)],Rho[ic+0*lbmaux[0].Ncells],k);
        Vel[ic] = Vel[ic] + F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic)]*lbmaux[0].C[k];
    }
    Rho[ic] *= lbmaux[0].Cs*lbmaux[0].Cs/3.0/(1.0-lbmaux[0].W[0]);


    real phi  = Rho[ic + 1*lbmaux[0].Ncells];
    real rho  = phi*lbmaux[0].rho[0] + (1.0-phi)*lbmaux[0].rho[1];
    Vel[ic]   = (lbmaux[0].Cs/rho)*(Vel[ic]+0.5*lbmaux[0].dt*BForce[ic]);

    BForce[ic+1*lbmaux[0].Ncells].x = BForce[ic+1*lbmaux[0].Ncells].y;
    BForce[ic] = make_real3(0.0,0.0,0.0);
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
        if (k!=0) Rho[ic                     ] += F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic)];
                  Rho[ic + 1*lbmaux[0].Ncells] += F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,1,k,ic)];
                  Rho[ic + 2*lbmaux[0].Ncells] += F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,2,k,ic)];
        
        //if (ic==217*lbmaux[0].Nx+0) printf("%g %g %g %lu \n",F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic)],Ftemp[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic)],Rho[ic+0*lbmaux[0].Ncells],k);
        Vel[ic] = Vel[ic] + F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,0,k,ic)]*lbmaux[0].C[k];
    }

    Rho[ic] *= lbmaux[0].Cs*lbmaux[0].Cs/3.0/(1.0-lbmaux[0].W[0]);

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
    BForce[ic] = make_real3(0.0,0.0,0.0);
}

__global__ void cudaStreamSW2(bool const * IsSolid, real * F, real * Ftemp, real3 * BForce, real3 * Vel, real * Rho, lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;
    if (ic==0)
    {
        lbmaux[0].Time += lbmaux[0].dt;
        lbmaux[0].iter++;
    }
    const size_t icl = ic/lbmaux[0].Ncells;
    const size_t ic0 = ic - icl*lbmaux[0].Ncells;
    //BForce[ic] = make_real3(0.0,0.0,0.0);
    //for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    //{
        //F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)] = Ftemp[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)];
    //}
    Rho   [ic] = 0.0;
    Vel   [ic] = make_real3(0.0,0.0,0.0);
    if (!IsSolid[ic])
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ic] += F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)];
            Vel[ic] = Vel[ic] + F[FIDX(lbmaux[0].Ncells,lbmaux[0].Nneigh,icl,k,ic0)]*lbmaux[0].C[k];
            //if (ic==0) printf("k: %lu %f %f %f %f \n",k,Rho[ic],Vel[ic].x,Vel[ic].y,Vel[ic].z);
        }
        Vel[ic] = lbmaux[0].Cs/Rho[ic]*Vel[ic];
    }
}
}
#endif //MECHSYS_LBM_CUH
