/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang                                         *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful, but  *
 * WITHOUT ANY WARRANTY; without even the implied warranty of           *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     *
 * General Public License for more details.                             *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

/////////////////////////////LBMMPM CUDA implementation////////////////

#ifndef MECHSYS_LBMMPM_CUH
#define MECHSYS_LBMMPM_CUH

// Mechsys
#include <mechsys/flbm/lbm.cuh>
#include <mechsys/linalg/matvec.h>

namespace LBMMPM
{

using FLBM::FIDX;

struct FacPair
{
    unsigned int facet;
    unsigned int cell;
};

struct BodyPair
{
    unsigned int body;
    unsigned int cell;
};

// The LBM distribution buffers use FLBM's direction-major layout:
//   F[(k)*Ncells + ic]
// The MPM body mesh is described by corner indices only; current positions and
// velocities are uploaded every coupling step, so the same kernels work while
// the structure moves.

__device__ __inline__ real3 lbmmpm_cross(real3 const & a, real3 const & b)
{
    return make_real3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

__device__ __inline__ bool lbmmpm_inside_tetra(real3 const & x,
        real3 const & v0, real3 const & v1, real3 const & v2, real3 const & v3)
{
    real Vc = fabs(dotreal3(v1-v0,lbmmpm_cross(v2-v0,v3-v0)))/6.0;
    real V0 = fabs(dotreal3(v1-x,lbmmpm_cross(v2-x,v3-x)))/6.0;
    real V1 = fabs(dotreal3(v0-x,lbmmpm_cross(v2-x,v3-x)))/6.0;
    real V2 = fabs(dotreal3(v1-x,lbmmpm_cross(v0-x,v3-x)))/6.0;
    real V3 = fabs(dotreal3(v1-x,lbmmpm_cross(v2-x,v0-x)))/6.0;
    return fabs(V0+V1+V2+V3-Vc)/Vc < 1.0e-8;
}

__global__ void cudaResetGamma(real * __restrict__ Gamma, size_t Ncells)
{
    size_t ic = threadIdx.x + blockIdx.x*blockDim.x;
    if (ic < Ncells) Gamma[ic] = 0.0;
}

// Mark lattice cells whose centres are inside a MPM body cell.  The host
// supplies a compact (body cell, candidate lattice cell) list, so one thread
// performs exactly one tetrahedron test instead of one thread scanning a whole
// AABB serially.
__global__ void cudaImprintCells(real3 const * __restrict__ CornerPos,
        size_t const * __restrict__ CellVerts, BodyPair const * __restrict__ Pairs,
        size_t numPairs, real * __restrict__ Gamma, bool * __restrict__ Inside,
        FLBM::lbm_aux const * __restrict__ lbmaux)
{
    size_t ip = threadIdx.x + blockIdx.x*blockDim.x;
    if (ip >= numPairs) return;
    size_t ib = Pairs[ip].body;
    size_t ic = Pairs[ip].cell;
    size_t b = 4*ib;
    real3 v0 = CornerPos[CellVerts[b+0]];
    real3 v1 = CornerPos[CellVerts[b+1]];
    real3 v2 = CornerPos[CellVerts[b+2]];
    real3 v3 = CornerPos[CellVerts[b+3]];
    int Nx = (int)lbmaux[0].Nx, Ny = (int)lbmaux[0].Ny;
    int nx = (int)(ic % Nx);
    int ny = (int)((ic/Nx) % Ny);
    int nz = (int)(ic/((size_t)Nx*Ny));
    real ld = lbmaux[0].dx;
    if (lbmmpm_inside_tetra(make_real3(nx*ld,ny*ld,nz*ld),v0,v1,v2,v3))
    {
        Gamma[ic]  = 1.0;
        Inside[ic] = true;
    }
}

// Refill the cells flagged Inside from their surrounding fluid cells.
// This mirrors FLBM::Domain::ImprintLattice's second loop.
__global__ void cudaRefill(real * __restrict__ F, real * __restrict__ Rho,
        real3 * __restrict__ Vel, real * __restrict__ Gamma,
        bool * __restrict__ Inside, FLBM::lbm_aux const * __restrict__ lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x*blockDim.x;
    if (ic >= lbmaux[0].Ncells) return;
    if (!(Gamma[ic] < 0.5 && Inside[ic])) return;

    int Nx = (int)lbmaux[0].Nx, Ny = (int)lbmaux[0].Ny, Nz = (int)lbmaux[0].Nz;
    int nx = (int)(ic % Nx);
    int ny = (int)((ic/Nx) % Ny);
    int nz = (int)(ic/((size_t)Nx*Ny));
    int NN = (int)lbmaux[0].Nneigh;
    size_t Nc = lbmaux[0].Ncells;
    real Ft[27];
    for (int k=0; k<NN; k++) Ft[k] = 0.0;
    int naem = 0;
    for (int k=0; k<NN; k++)
    {
        int cx = (int)lbmaux[0].C[k].x;
        int cy = (int)lbmaux[0].C[k].y;
        int cz = (int)lbmaux[0].C[k].z;
        int nxf = nx+cx, nyf = ny+cy, nzf = nz+cz;
        int nxff = nx+2*cx, nyff = ny+2*cy, nzff = nz+2*cz;
        if (nxf<0 || nyf<0 || nzf<0 || nxf>=Nx || nyf>=Ny || nzf>=Nz) continue;
        if (nxff<0 || nyff<0 || nzff<0 || nxff>=Nx || nyff>=Ny || nzff>=Nz)
        {
            nxff=nxf; nyff=nyf; nzff=nzf;
        }
        size_t nf  = nxf + (size_t)nyf*Nx + (size_t)nzf*Nx*Ny;
        size_t nff = nxff + (size_t)nyff*Nx + (size_t)nzff*Nx*Ny;
        if (Gamma[nf] < 0.5 && !Inside[nf])
        {
            for (int kt=0; kt<NN; kt++)
                Ft[kt] += 2.0*F[(size_t)kt*Nc+nf] - F[(size_t)kt*Nc+nff];
            naem++;
        }
    }
    Rho[ic] = 0.0;
    Vel[ic] = make_real3(0.0,0.0,0.0);
    if (naem == 0)
    {
        for (int k=0;k<NN;k++) F[(size_t)k*Nc+ic] = 0.0;
        Rho[ic] = 1.0;
        Inside[ic] = false;
        return;
    }
    for (int k=0;k<NN;k++)
    {
        Ft[k] = fabs(Ft[k])/naem;
        F[(size_t)k*Nc+ic] = Ft[k];
        Rho[ic] += Ft[k];
        Vel[ic]   = Vel[ic] + Ft[k]*lbmaux[0].C[k];
    }
    Vel[ic] = (lbmaux[0].Cs/Rho[ic])*Vel[ic];
    Inside[ic] = false;
}

__device__ __inline__ bool lbmmpm_facet_hit(real3 const & x0, real3 const & x1,
        real3 const & x2, real3 const & xn, real3 const & nor,
        FLBM::lbm_aux const * lbmaux, int & kout, real & rout, real & sout, real & tout)
{
    int NN = (int)lbmaux[0].Nneigh;
    real ld = lbmaux[0].dx;
    real3 B = x0 - xn;
    for (int k=1; k<NN; k++)
    {
        Mat3 M(ld*lbmaux[0].C[k].x, ld*lbmaux[0].C[k].y, ld*lbmaux[0].C[k].z,
               x0.x-x1.x,           x0.y-x1.y,           x0.z-x1.z,
               x0.x-x2.x,           x0.y-x2.y,           x0.z-x2.z);
        real3 rst;
        if (!M.SolAlt(B,rst,1.0e-20)) continue;
        real r=rst.x, s=rst.y, t=rst.z;
        // Match the strict CPU test exactly.  LBMMPM CUDA is compiled with
        // --fmad=false so Mat3::SolAlt and the CPU SolAlt round the same way
        // on facet edges; with default FMA contraction the ~1e-16 edge noise
        // changes which surface cells are marked.
        if (!(r>0.0 && r<=1.0 && s>0.0 && t>0.0 && s+t<1.0)) continue;
        if (!(dotreal3(xn-x0,nor) > 0.0)) continue;
        kout=k; rout=r; sout=s; tout=t;
        return true;
    }
    return false;
}

// Mark the facet-overlapped cells first.  The forcing kernel below reads
// Gamma of the neighbouring cell; doing both jobs in one launch made that read
// race with the marking writes and silently dropped surface forces.
__global__ void cudaMarkFacets(real3 const * __restrict__ CornerPos,
        size_t const * __restrict__ FacVerts, FacPair const * __restrict__ Pairs,
        size_t numPairs, real * __restrict__ Gamma,
        FLBM::lbm_aux const * __restrict__ lbmaux)
{
    size_t ip = threadIdx.x + blockIdx.x*blockDim.x;
    if (ip >= numPairs) return;
    size_t ibf = Pairs[ip].facet;
    size_t ic  = Pairs[ip].cell;
    size_t i0 = FacVerts[3*ibf+0];
    size_t i1 = FacVerts[3*ibf+1];
    size_t i2 = FacVerts[3*ibf+2];
    real3 x0 = CornerPos[i0], x1 = CornerPos[i1], x2 = CornerPos[i2];
    real3 nor = lbmmpm_cross(x1-x0,x2-x1);
    real nrm = norm(nor);
    if (nrm > 1.0e-30) nor = nor/nrm;

    int Nx = (int)lbmaux[0].Nx, Ny = (int)lbmaux[0].Ny;
    int nx = (int)(ic % Nx);
    int ny = (int)((ic/Nx) % Ny);
    int nz = (int)(ic/((size_t)Nx*Ny));
    real3 xn = make_real3(nx*lbmaux[0].dx,ny*lbmaux[0].dx,nz*lbmaux[0].dx);
    int k; real r,s,t;
    if (lbmmpm_facet_hit(x0,x1,x2,xn,nor,lbmaux,k,r,s,t)) Gamma[ic] = 2.0;
}

// Immersed-boundary forcing over the MPM boundary facets.  The host builds a
// compact list of (facet, candidate cell) pairs from the facet AABBs; one
// thread handles one pair and scans the lattice directions.  Building the list
// on the CPU keeps the GPU work proportional to the actual overlap instead of
// launching one under-parallelised thread for every facet.
__global__ void cudaImprintFacets(real3 const * __restrict__ CornerPos,
        real3 const * __restrict__ CornerVel,
        size_t const * __restrict__ FacVerts, FacPair const * __restrict__ Pairs,
        size_t numPairs, real * __restrict__ F, real const * __restrict__ Rho,
        real * __restrict__ Gamma, real3 * __restrict__ CornerH,
        FLBM::lbm_aux const * __restrict__ lbmaux)
{
    size_t ip = threadIdx.x + blockIdx.x*blockDim.x;
    if (ip >= numPairs) return;
    size_t ibf = Pairs[ip].facet;
    size_t ic  = Pairs[ip].cell;
    size_t i0 = FacVerts[3*ibf+0];
    size_t i1 = FacVerts[3*ibf+1];
    size_t i2 = FacVerts[3*ibf+2];
    real3 x0 = CornerPos[i0], x1 = CornerPos[i1], x2 = CornerPos[i2];
    real3 v0 = CornerVel[i0], v1 = CornerVel[i1], v2 = CornerVel[i2];
    real3 nor = lbmmpm_cross(x1-x0,x2-x1);
    real nrm = norm(nor);
    if (nrm > 1.0e-30) nor = nor/nrm;

    int Nx = (int)lbmaux[0].Nx, Ny = (int)lbmaux[0].Ny, Nz = (int)lbmaux[0].Nz;
    int nx = (int)(ic % Nx);
    int ny = (int)((ic/Nx) % Ny);
    int nz = (int)(ic/((size_t)Nx*Ny));
    real ld = lbmaux[0].dx;
    real3 xn = make_real3(nx*ld,ny*ld,nz*ld);
    int Nc = (int)lbmaux[0].Ncells;
    int NN = (int)lbmaux[0].Nneigh;
    real3 B = x0 - xn;
    for (int k=1; k<NN; k++)
    {
        Mat3 M(ld*lbmaux[0].C[k].x, ld*lbmaux[0].C[k].y, ld*lbmaux[0].C[k].z,
               x0.x-x1.x,           x0.y-x1.y,           x0.z-x1.z,
               x0.x-x2.x,           x0.y-x2.y,           x0.z-x2.z);
        real3 rst;
        if (!M.SolAlt(B,rst,1.0e-20)) continue;
        real r=rst.x, s=rst.y, t=rst.z;
        if (!(r>0.0 && r<=1.0 && s>0.0 && t>0.0 && s+t<1.0)) continue;
        if (!(dotreal3(xn-x0,nor) > 0.0)) continue;
        Gamma[ic] = 2.0;
        int ko = (int)lbmaux[0].Op[k];
        real3 vw = v0 + s*(v1-v0) + t*(v2-v0);
        int nxf=nx+(int)lbmaux[0].C[ko].x;
        int nyf=ny+(int)lbmaux[0].C[ko].y;
        int nzf=nz+(int)lbmaux[0].C[ko].z;
        int nxs=nx+(int)lbmaux[0].C[k].x;
        int nys=ny+(int)lbmaux[0].C[k].y;
        int nzs=nz+(int)lbmaux[0].C[k].z;
        if (nxf<0 || nyf<0 || nzf<0 || nxf>=Nx || nyf>=Ny || nzf>=Nz) continue;
        if (nxs<0 || nys<0 || nzs<0 || nxs>=Nx || nys>=Ny || nzs>=Nz) continue;
        size_t nf  = nxf + (size_t)nyf*Nx + (size_t)nzf*Nx*Ny;
        size_t ns  = nxs + (size_t)nys*Nx + (size_t)nzs*Nx*Ny;
        if (Gamma[ns] < 0.5) continue;
        real rho = Rho[ic];
        real Fk  = F[(size_t)k*Nc+ic];
        real Fnf = F[(size_t)ko*Nc+nf];
        real Fns = F[(size_t)k*Nc+ns];
        real newF = (r*Fnf + (1.0-r)*Fk + r*Fns
                    + 6.0*rho*lbmaux[0].W[k]*dotreal3(lbmaux[0].C[ko],vw)/lbmaux[0].Cs)/(1.0+r);
        if (newF < 0.0) newF = 0.0;
        F[(size_t)ko*Nc+ic] = newF;

        real3 Fv = (ld*ld*lbmaux[0].Cs)*
            (Fns*(lbmaux[0].Cs*lbmaux[0].C[k]-vw)
            -newF*(lbmaux[0].Cs*lbmaux[0].C[ko]-vw)
            -2.0*rho*lbmaux[0].W[k]*lbmaux[0].Cs*lbmaux[0].C[k]);
        real f = norm(Fv);
        if (f > 1.0e-12)
        {
            real3 F0 = (1.0-s-t)*Fv, F1 = s*Fv, F2 = t*Fv;
            atomicAdd(&CornerH[i0].x,F0.x); atomicAdd(&CornerH[i0].y,F0.y); atomicAdd(&CornerH[i0].z,F0.z);
            atomicAdd(&CornerH[i1].x,F1.x); atomicAdd(&CornerH[i1].y,F1.y); atomicAdd(&CornerH[i1].z,F1.z);
            atomicAdd(&CornerH[i2].x,F2.x); atomicAdd(&CornerH[i2].y,F2.y); atomicAdd(&CornerH[i2].z,F2.z);
        }
    }
}

__global__ void cudaVelDen(real const * __restrict__ F, real * __restrict__ Rho,
        real3 * __restrict__ Vel, bool const * __restrict__ IsSolid,
        FLBM::lbm_aux const * __restrict__ lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x*blockDim.x;
    if (ic >= lbmaux[0].Ncells) return;
    Rho[ic] = 0.0;
    Vel[ic] = make_real3(0.0,0.0,0.0);
    if (IsSolid[ic]) return;
    int NN = (int)lbmaux[0].Nneigh;
    size_t Nc = lbmaux[0].Ncells;
    for (int k=0;k<NN;k++)
    {
        real f = F[(size_t)k*Nc+ic];
        Rho[ic] += f;
        Vel[ic] = Vel[ic] + f*lbmaux[0].C[k];
    }
    Vel[ic] = (lbmaux[0].Cs/Rho[ic])*Vel[ic];
}

template<int NN>
__global__ void cudaCollideMPM(real const * __restrict__ F, real * __restrict__ Ftemp,
        real3 const * __restrict__ BForce, real * __restrict__ Rho,
        real3 * __restrict__ Vel, bool const * __restrict__ IsSolid,
        real tau0, real Sc, FLBM::lbm_aux const * __restrict__ lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x*blockDim.x;
    if (ic >= lbmaux[0].Ncells) return;
    size_t Nc = lbmaux[0].Ncells;
    if (IsSolid[ic])
    {
        Rho[ic] = 0.0;
        Vel[ic] = make_real3(0.0,0.0,0.0);
        for (int k=0;k<NN;k++)
            Ftemp[(size_t)k*Nc+ic] = F[(size_t)lbmaux[0].Op[k]*Nc+ic];
        return;
    }
    // Fused cudaVelDen: reconstruct rho/u from F before collision.
    real rho = 0.0;
    real3 vel0 = make_real3(0.0,0.0,0.0);
    real Cs = lbmaux[0].Cs;
    for (int k=0;k<NN;k++)
    {
        real f = F[(size_t)k*Nc+ic];
        rho  += f;
        vel0 = vel0 + f*lbmaux[0].C[k];
    }
    real invRho = 1.0/rho;
    Vel[ic] = (Cs*invRho)*vel0;
    Rho[ic] = rho;

    real NonEq[NN];
    real Q = 0.0;
    real tau = tau0;
    real3 vel = Vel[ic] + lbmaux[0].dt*tau*invRho*BForce[ic];
    real VdotV = dotreal3(vel,vel);
    real invCs = 1.0/Cs;
    real invCs2 = invCs*invCs;
    real V2 = 1.5*VdotV*invCs2;
    for (int k=0;k<NN;k++)
    {
        real cu = dotreal3(vel,lbmaux[0].C[k])*invCs;
        real Feq = lbmaux[0].W[k]*rho*(1.0+3.0*cu+4.5*cu*cu-V2);
        NonEq[k] = F[(size_t)k*Nc+ic] - Feq;
        Q += NonEq[k]*NonEq[k]*lbmaux[0].EEk[k];
    }
    Q = sqrt(2.0*Q);
    tau = 0.5*(tau+sqrt(tau*tau+6.0*Q*Sc*invRho));
    // Same fixed-point iteration as the CPU limiter, but defer the Ftemp
    // stores until the final alpha is known.
    real alpha = 1.0, alphat = 1.0;
    bool valid = true;
    int ntry = 0;
    while (valid && ntry < 16)
    {
        ntry++;
        valid = false;
        alpha = alphat;
        for (int k=0;k<NN;k++)
        {
            real noneq = NonEq[k]/tau;
            real f = F[(size_t)k*Nc+ic] - alpha*noneq;
            if (f < -1.0e-6)
            {
                real temp = F[(size_t)k*Nc+ic]/noneq;
                if (temp < alphat) alphat = temp;
                valid = true;
            }
        }
    }
    for (int k=0;k<NN;k++)
        Ftemp[(size_t)k*Nc+ic] = F[(size_t)k*Nc+ic] - alpha*(NonEq[k]/tau);
}

__global__ void cudaStreamMPM(real const * __restrict__ F, real * __restrict__ Ftemp,
        FLBM::lbm_aux const * __restrict__ lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x*blockDim.x;
    if (ic >= lbmaux[0].Ncells) return;
    int Nx=(int)lbmaux[0].Nx, Ny=(int)lbmaux[0].Ny, Nz=(int)lbmaux[0].Nz;
    int nx=(int)(ic%Nx), ny=(int)((ic/Nx)%Ny), nz=(int)(ic/((size_t)Nx*Ny));
    size_t Nc=lbmaux[0].Ncells, il=0;
    for (int k=0;k<(int)lbmaux[0].Nneigh;k++)
    {
        int mx = nx + (int)lbmaux[0].C[k].x;
        int my = ny + (int)lbmaux[0].C[k].y;
        int mz = nz + (int)lbmaux[0].C[k].z;
        if (mx<0) mx += Nx; else if (mx>=Nx) mx -= Nx;
        if (my<0) my += Ny; else if (my>=Ny) my -= Ny;
        if (mz<0) mz += Nz; else if (mz>=Nz) mz -= Nz;
        size_t dest=mx+(size_t)my*Nx+(size_t)mz*Nx*Ny;
        Ftemp[FIDX(Nc,lbmaux[0].Nneigh,il,k,dest)] = F[FIDX(Nc,lbmaux[0].Nneigh,il,k,ic)];
    }
}

} // namespace LBMMPM
#endif
