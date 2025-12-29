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
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/


#ifndef MECHSYS_MPM_CUH
#define MECHSYS_MPM_CUH
namespace MPM
{

struct mpm_aux
{
    real   Mmin;         ///< The minimun mass of particle
    real   Dx;           ///< Grid size in physical units
    real   Dt;           ///< Time step
    real   Time;         ///< Time clock
    real   Tf;           ///< Final Time
    real   Gn;           ///< Dissipative constant
    size_t iter   = 0;   ///< Iteration clock;
    size_t Nx;           ///< Integer vector with the dimensions of the MPM node mesh
    size_t Ny;           ///< Integer vector with the dimensions of the MPM node mesh
    size_t Nz;           ///< Integer vector with the dimensions of the MPM node mesh
    size_t Nnodes;       ///< Integer vector with the dimensions of the MPM node mesh
    size_t Npart;        ///< Number of MPM particles
    real3  Xlow;         ///< vectors of extreme positions of the MPM underlaying mesh 
    real3  Xhigh;
};

__global__ void cudaResetNode(NodeCU * Nodes, size_t * SelIndexes, size_t * validn, mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=validn[0]) return;
    size_t in = SelIndexes[ic];
    Nodes[in].Reset();
}

__global__ void cudaParticleToNode(ParticleCU * Par, NodeCU * Nodes, bool * ValidMask, mpm_aux * mpmaux)
{   
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=mpmaux[0].Npart) return;
    real3 Xp  = Par[ic].x;
    real3 Lp  = Par[ic].Dx;
    Mat3 Vsp  = -Par[ic].V*Par[ic].S;
    real3 Fex = Par[ic].b;
    real Dx   = mpmaux[0].Dx;
    real Dt   = mpmaux[0].Dt;

    size_t nlx  = (size_t)fmax((trunc((Xp.x-0.5*Lp.x-mpmaux[0].Xlow.x)/Dx)-1),0.0);
    size_t nly  = (size_t)fmax((trunc((Xp.y-0.5*Lp.y-mpmaux[0].Xlow.y)/Dx)-1),0.0);
    size_t nlz  = (size_t)fmax((trunc((Xp.z-0.5*Lp.z-mpmaux[0].Xlow.z)/Dx)-1),0.0);
    size_t nhx =  (size_t)fmin((ceil ((Xp.x+0.5*Lp.x-mpmaux[0].Xlow.x)/Dx)+1),(real)mpmaux[0].Nx-1.0);
    size_t nhy =  (size_t)fmin((ceil ((Xp.y+0.5*Lp.y-mpmaux[0].Xlow.y)/Dx)+1),(real)mpmaux[0].Ny-1.0);
    size_t nhz =  (size_t)fmin((ceil ((Xp.z+0.5*Lp.z-mpmaux[0].Xlow.z)/Dx)+1),(real)mpmaux[0].Nz-1.0);

    Par[ic].Pnodes = 0;

    for (size_t nx=nlx; nx<=nhx; nx++)
    for (size_t ny=nly; ny<=nhy; ny++)
    for (size_t nz=nlz; nz<=nhz; nz++)
    {
        real Sf=0.0; //Shape function value
        real3 Gs = make_real3(0.0,0.0,0.0);
        real3 Xn = make_real3(nx*Dx,ny*Dx,nz*Dx);
        Xn = Xn + mpmaux[0].Xlow;
        GIMP3DCU(Xp,Xn,Dx,Lp,Sf,Gs);
        if (Sf>0.0)
        {
            size_t nc  = nx + ny*mpmaux[0].Nx + nz*mpmaux[0].Nx*mpmaux[0].Ny;
            Par[ic].Nodes[Par[ic].Pnodes] = nc;
            Par[ic].Sfunc[Par[ic].Pnodes] = Sf;
            Par[ic].GSf  [Par[ic].Pnodes] = Gs;
            Par[ic].Pnodes++;
            real ms = Sf*Par[ic].m;
            real3 VspGs = Vsp*Gs;
            real3 dF = Sf*Fex + VspGs;
            atomicAdd(&Nodes[nc].Mass, ms);
            atomicAdd(&Nodes[nc].Mn.x, ms*Par[ic].v.x+Dt*dF.x);
            atomicAdd(&Nodes[nc].Mn.y, ms*Par[ic].v.y+Dt*dF.y);
            atomicAdd(&Nodes[nc].Mn.z, ms*Par[ic].v.z+Dt*dF.z);
            atomicAdd(&Nodes[nc].Fn.x, dF.x);
            atomicAdd(&Nodes[nc].Fn.y, dF.y);
            atomicAdd(&Nodes[nc].Fn.z, dF.z);
            ValidMask[nc] = true;
        }
    }
}

__global__ void cudaSolveNode(NodeCU * Nodes, size_t * SelIndexes, size_t * validn, mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=validn[0]) return;
    size_t in = SelIndexes[ic];
    Nodes[in].UpdateNode(mpmaux[0].Gn,mpmaux[0].Dt,mpmaux[0].Mmin);
    Nodes[in].FixNode();
}

__global__ void cudaNodeToParticle (ParticleCU * Par, NodeCU * Nodes, mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=mpmaux[0].Npart) return;
    if (ic==0)
    {
        mpmaux[0].iter++;
        mpmaux[0].Time += mpmaux[0].Dt;
    }

    Par[ic].L.SetZero();

    for (size_t nn=0; nn<Par[ic].Pnodes; nn++)
    {
        size_t nc = Par[ic].Nodes[nn];
        double sf = Par[ic].Sfunc[nn];
        real3  gs = Par[ic].GSf[nn];
        real3  an = Nodes[nc].Fn/Nodes[nc].Mass;
        Par[ic].v = Par[ic].v + sf*mpmaux[0].Dt*an;
        Par[ic].x = Par[ic].x + sf*mpmaux[0].Dt*Nodes[nc].Vn;
        Mat3 Gv = Dyad(gs,Nodes[nc].Vn);
        Par[ic].L += Gv;
    }

    Par[ic].CalcVol   (mpmaux[0].Dt);
    Par[ic].CalcStress(mpmaux[0].Dt);
    Par[ic].Reset     (mpmaux[0].Dt);
}



}
#endif
