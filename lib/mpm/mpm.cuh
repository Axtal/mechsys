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
    int    ShapeFunction; ///< Particle-grid interpolation function
    size_t iter   = 0;   ///< Iteration clock;
    size_t Nx;           ///< Integer vector with the dimensions of the MPM node mesh
    size_t Ny;           ///< Integer vector with the dimensions of the MPM node mesh
    size_t Nz;           ///< Integer vector with the dimensions of the MPM node mesh
    size_t Nnodes;       ///< Integer vector with the dimensions of the MPM node mesh
    size_t Npart;        ///< Number of MPM particles
    real3  Xlow;         ///< vectors of extreme positions of the MPM underlaying mesh 
    real3  Xhigh;
    bool   PeriodicX;    ///< Periodic particle-grid mapping along x
    bool   PeriodicY;    ///< Periodic particle-grid mapping along y
    bool   PeriodicZ;    ///< Periodic particle-grid mapping along z
};

//Reset only the nodes touched in the previous step, using the compact list
//built by the scatter.  The hot fields live in separate SoA arrays.
__global__ void cudaResetNode(real * Mass, real3 * Mn, real3 * Fn, int * valid,
                              size_t * SelIndexes, unsigned long long * validn, mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=validn[0]) return;
    size_t in = SelIndexes[ic];
    Mass[in] = 0.0;
    Mn  [in] = make_real3(0.0,0.0,0.0);
    Fn  [in] = make_real3(0.0,0.0,0.0);
    valid[in] = 0;
}

//Scatter, and at the same time build the compact list of touched nodes: the
//first contribution to a node appends its index.  This removes the old
//per-step thrust::fill + copy_if + thrust::distance + cudaMemcpy sequence.
__global__ void cudaParticleToNode(ParticleCU * Par, real * Mass, real3 * Mn, real3 * Fn,
                                  int * valid, size_t * SelIndexes, unsigned long long * validn,
                                  mpm_aux * mpmaux)
{   
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=mpmaux[0].Npart) return;
    real3 Xp  = Par[ic].x;
    real3 Lp  = make_real3(fabs(Par[ic].Dx.x),fabs(Par[ic].Dx.y),fabs(Par[ic].Dx.z));
    Mat3 Vsp  = -Par[ic].V*Par[ic].S;
    real3 Fex = Par[ic].b;
    real Dx   = mpmaux[0].Dx;
    real Dt   = mpmaux[0].Dt;

    int nlx = 0;
    int nly = 0;
    int nlz = 0;
    int nhx = 0;
    int nhy = 0;
    int nhz = 0;
    if (mpmaux[0].PeriodicX)
    {
        ShapeNodeRangeRawCU(mpmaux[0].ShapeFunction,Xp.x,mpmaux[0].Xlow.x,Dx,Lp.x,nlx,nhx);
    }
    else
    {
        size_t nl = 0;
        size_t nh = 0;
        ShapeNodeRangeCU(mpmaux[0].ShapeFunction,Xp.x,mpmaux[0].Xlow.x,Dx,Lp.x,mpmaux[0].Nx,nl,nh);
        nlx = (int)nl; nhx = (int)nh;
    }
    if (mpmaux[0].PeriodicY)
    {
        ShapeNodeRangeRawCU(mpmaux[0].ShapeFunction,Xp.y,mpmaux[0].Xlow.y,Dx,Lp.y,nly,nhy);
    }
    else
    {
        size_t nl = 0;
        size_t nh = 0;
        ShapeNodeRangeCU(mpmaux[0].ShapeFunction,Xp.y,mpmaux[0].Xlow.y,Dx,Lp.y,mpmaux[0].Ny,nl,nh);
        nly = (int)nl; nhy = (int)nh;
    }
    if (mpmaux[0].PeriodicZ)
    {
        ShapeNodeRangeRawCU(mpmaux[0].ShapeFunction,Xp.z,mpmaux[0].Xlow.z,Dx,Lp.z,nlz,nhz);
    }
    else
    {
        size_t nl = 0;
        size_t nh = 0;
        ShapeNodeRangeCU(mpmaux[0].ShapeFunction,Xp.z,mpmaux[0].Xlow.z,Dx,Lp.z,mpmaux[0].Nz,nl,nh);
        nlz = (int)nl; nhz = (int)nh;
    }

    Par[ic].Pnodes = 0;

    for (int nx=nlx; nx<=nhx; nx++)
    for (int ny=nly; ny<=nhy; ny++)
    for (int nz=nlz; nz<=nhz; nz++)
    {
        real Sf=0.0; //Shape function value
        real3 Gs = make_real3(0.0,0.0,0.0);
        real3 Xn = make_real3(nx*Dx,ny*Dx,nz*Dx);
        Xn = Xn + mpmaux[0].Xlow;
        MPMShapeFunction3DCU(mpmaux[0].ShapeFunction,Xp,Xn,Dx,Lp,Sf,Gs);
        if (Sf>0.0)
        {
            size_t ix = mpmaux[0].PeriodicX ? PeriodicNodeIndexCU(nx,mpmaux[0].Nx) : (size_t)nx;
            size_t iy = mpmaux[0].PeriodicY ? PeriodicNodeIndexCU(ny,mpmaux[0].Ny) : (size_t)ny;
            size_t iz = mpmaux[0].PeriodicZ ? PeriodicNodeIndexCU(nz,mpmaux[0].Nz) : (size_t)nz;
            size_t nc = ix + iy*mpmaux[0].Nx + iz*mpmaux[0].Nx*mpmaux[0].Ny;
            Par[ic].Nodes[Par[ic].Pnodes] = nc;
            Par[ic].Sfunc[Par[ic].Pnodes] = Sf;
            Par[ic].GSf  [Par[ic].Pnodes] = Gs;
            Par[ic].Pnodes++;
            real ms = Sf*Par[ic].m;
            real3 VspGs = Vsp*Gs;
            real3 dF = Sf*Fex + VspGs;
            atomicAdd(&Mass[nc], ms);
            atomicAdd(&Mn[nc].x, ms*Par[ic].v.x+Dt*dF.x);
            atomicAdd(&Mn[nc].y, ms*Par[ic].v.y+Dt*dF.y);
            atomicAdd(&Mn[nc].z, ms*Par[ic].v.z+Dt*dF.z);
            atomicAdd(&Fn[nc].x, dF.x);
            atomicAdd(&Fn[nc].y, dF.y);
            atomicAdd(&Fn[nc].z, dF.z);
            if (atomicExch(&valid[nc],1)==0)
            {
                unsigned long long k = atomicAdd(validn,1);
                SelIndexes[k] = nc;
            }
        }
    }
}

//Nodal update over the touched list.  The launch covers all nodes and exits
//early, so the host never needs the list size.
__global__ void cudaSolveNode(NodeCU * Nodes, real * Mass, real3 * Mn, real3 * Fn, real3 * Vn,
                              size_t * SelIndexes, unsigned long long * validn, mpm_aux * mpmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=validn[0]) return;
    size_t in = SelIndexes[ic];
    real  m  = Mass[in];
    real3 mn = Mn  [in];
    real3 fn = Fn  [in];
    real3 fdamp = mpmaux[0].Gn*mn;
    fn = fn - fdamp;
    mn = mn - mpmaux[0].Dt*fdamp;
    real3 vn;
    if (m<mpmaux[0].Mmin) vn = make_real3(0.0,0.0,0.0);
    else                  vn = (1.0/m)*mn;
    NodeCU & nd = Nodes[in];
    if (nd.vxf) { fn.x = 0.0; vn.x = nd.Vnf.x; mn.x = m*vn.x; }
    if (nd.vyf) { fn.y = 0.0; vn.y = nd.Vnf.y; mn.y = m*vn.y; }
    if (nd.vzf) { fn.z = 0.0; vn.z = nd.Vnf.z; mn.z = m*vn.z; }
    Mass[in] = m;
    Mn  [in] = mn;
    Fn  [in] = fn;
    Vn  [in] = vn;
}

//Gather.  Mass/Fn/Vn are read from the coalesced SoA arrays.
__global__ void cudaNodeToParticle (ParticleCU * Par, real * Mass, real3 * Fn, real3 * Vn, mpm_aux * mpmaux)
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
        real sf = Par[ic].Sfunc[nn];
        real3  gs = Par[ic].GSf[nn];
        real3  an = Fn[nc]/Mass[nc];
        Par[ic].v = Par[ic].v + sf*mpmaux[0].Dt*an;
        Par[ic].x = Par[ic].x + sf*mpmaux[0].Dt*Vn[nc];
        Mat3 Gv = Dyad(gs,Vn[nc]);
        Par[ic].L += Gv;
    }

    Par[ic].CalcVol   (mpmaux[0].Dt);
    Par[ic].CalcStress(mpmaux[0].Dt);
    Par[ic].Reset     (mpmaux[0].Dt);
    if (mpmaux[0].PeriodicX)
    {
        real L = mpmaux[0].Xhigh.x - mpmaux[0].Xlow.x;
        while (Par[ic].x.x <  mpmaux[0].Xlow .x) {Par[ic].x.x += L; Par[ic].xb.x += L;}
        while (Par[ic].x.x >= mpmaux[0].Xhigh.x) {Par[ic].x.x -= L; Par[ic].xb.x -= L;}
    }
    if (mpmaux[0].PeriodicY)
    {
        real L = mpmaux[0].Xhigh.y - mpmaux[0].Xlow.y;
        while (Par[ic].x.y <  mpmaux[0].Xlow .y) {Par[ic].x.y += L; Par[ic].xb.y += L;}
        while (Par[ic].x.y >= mpmaux[0].Xhigh.y) {Par[ic].x.y -= L; Par[ic].xb.y -= L;}
    }
    if (mpmaux[0].PeriodicZ)
    {
        real L = mpmaux[0].Xhigh.z - mpmaux[0].Xlow.z;
        while (Par[ic].x.z <  mpmaux[0].Xlow .z) {Par[ic].x.z += L; Par[ic].xb.z += L;}
        while (Par[ic].x.z >= mpmaux[0].Xhigh.z) {Par[ic].x.z -= L; Par[ic].xb.z -= L;}
    }
}



}
#endif
