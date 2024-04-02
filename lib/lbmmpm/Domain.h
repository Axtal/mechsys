/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
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


#ifndef MECHSYS_LBMMPM_DOMAIN_H
#define MECHSYS_LBMMPM_DOMAIN_H

// Mechsys
#include <mechsys/flbm/Domain.h>
#include <mechsys/mpm/Domain.h>

namespace LBMMPM
{

struct MtData
{
    Array<size_t > vacant; ///< An array of vacnat nodes
};

class Domain
{
public:

    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    // Constructors
    Domain ();                    ///< Default constructor

    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    double                    nu, ///< Viscosity of the fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Methods
    //Build LBM mesh
    void LBMmesh (LBMethod      Method, ///< Type of array, for example D2Q9
          double                    nu, ///< Viscosity of the fluid
          iVec3_t               Ndim,   ///< Cell divisions per side
          double                dx,     ///< Space spacing
          double                dt);    ///< Time step

    void Reset();                    ///< Reset LBM grid
    void ImprintLattice();              ///< Imprint the MPM material point into the LBM grid
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);            ///< Solve the Domain dynamics
    
    //Variables
    double           dt;                ///< Time Step
    size_t        Nproc;                ///< Number of processors
    size_t      idx_out;                ///< The discrete time step for output
    String      FileKey;                ///< File Key for output files
    void *     UserData;                ///< User Data
    double         Time;                ///< Simulation time variable
    MtData *        MTD;
    FLBM::Domain LBMDOM;                ///< The LBM domain
    MPM::Domain  MPMDOM;                ///< The MPM domain
};

inline Domain::Domain()
{
    Nproc = 1;
    Time = 0.0;
    MPMDOM = MPM::Domain();
}

inline Domain::Domain(LBMethod TheMethod, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Nproc  = 1;
    dt = Thedt;
    Time = 0.0;
    LBMDOM = FLBM::Domain(TheMethod, Thenu, TheNdim, Thedx, Thedt);
    MPMDOM = MPM ::Domain();

    //LBMDOM.Omeis  = new double *** [TheNdim(0)];
    LBMDOM.Inside = new bool   **  [TheNdim(0)];
    LBMDOM.Gamma  = new double **  [TheNdim(0)];
    for (size_t ix=0; ix< TheNdim(0); ix++)
    {
        //LBMDOM.Omeis [ix] = new double ** [TheNdim(1)];
        LBMDOM.Inside[ix] = new bool   *  [TheNdim(1)];
        LBMDOM.Gamma [ix] = new double *  [TheNdim(1)];
        for (size_t iy=0; iy< TheNdim(1); iy++)
        {
            //LBMDOM.Omeis [ix][iy] = new double * [TheNdim(2)];
            LBMDOM.Inside[ix][iy] = new bool     [TheNdim(2)];
            LBMDOM.Gamma [ix][iy] = new double   [TheNdim(2)];
            for (size_t iz=0; iz< TheNdim(2); iz++)
            {
                //LBMDOM.Omeis [ix][iy][iz] = new double [LBMDOM.Nneigh];
                LBMDOM.Inside[ix][iy][iz] = false;
                LBMDOM.Gamma [ix][iy][iz] = 0.0;
            }
        }
    }
}

inline void Domain::LBMmesh(LBMethod TheMethod, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)
{
    dt = Thedt;
    LBMDOM = FLBM::Domain(TheMethod, Thenu, TheNdim, Thedx, Thedt);

    //LBMDOM.Omeis  = new double *** [TheNdim(0)];
    LBMDOM.Inside = new bool   **  [TheNdim(0)];
    LBMDOM.Gamma  = new double **  [TheNdim(0)];
    for (size_t ix=0; ix< TheNdim(0); ix++)
    {
        //LBMDOM.Omeis [ix] = new double ** [TheNdim(1)];
        LBMDOM.Inside[ix] = new bool   *  [TheNdim(1)];
        LBMDOM.Gamma [ix] = new double *  [TheNdim(1)];
        for (size_t iy=0; iy< TheNdim(1); iy++)
        {
            //LBMDOM.Omeis [ix][iy] = new double * [TheNdim(2)];
            LBMDOM.Inside[ix][iy] = new bool     [TheNdim(2)];
            LBMDOM.Gamma [ix][iy] = new double   [TheNdim(2)];
            for (size_t iz=0; iz< TheNdim(2); iz++)
            {
                //LBMDOM.Omeis [ix][iy][iz] = new double [LBMDOM.Nneigh];
                LBMDOM.Inside[ix][iy][iz] = false;
                LBMDOM.Gamma [ix][iy][iz] = 0.0;
            }
        }
    }
}

void Domain::ImprintLattice()
{
    double Cs = LBMDOM.Cs;
    double ld = LBMDOM.dx;
    double Tau = LBMDOM.Tau[0];

    ///*
        //std::cout << "0" << std::endl;
    //Finding the cells inside the solid body
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < MPMDOM.BodyMesh->Cells.Size(); ip++)
    {
        iVec3_t nlow,nhigh;
        Vec3_t xmin = MPMDOM.BodyMesh->Cells[ip]->Xmin;
        Vec3_t xmax = MPMDOM.BodyMesh->Cells[ip]->Xmax;

        //size_t i0 = MPMDOM.BodyMesh->Cells[ip]->V[0]->ID;
        //size_t i1 = MPMDOM.BodyMesh->Cells[ip]->V[1]->ID;
        //size_t i2 = MPMDOM.BodyMesh->Cells[ip]->V[2]->ID;
        //size_t i3 = MPMDOM.BodyMesh->Cells[ip]->V[3]->ID;
//
        //Vec3_t v0 = MPMDOM.Corners[i0]->v;
        //Vec3_t v1 = MPMDOM.Corners[i1]->v;
        //Vec3_t v2 = MPMDOM.Corners[i2]->v;
        //Vec3_t v3 = MPMDOM.Corners[i3]->v;

        nlow (0) = (size_t)std::max((trunc((xmin(0)-ld)/ld)),0.0);
        nlow (1) = (size_t)std::max((trunc((xmin(1)-ld)/ld)),0.0);
        nlow (2) = (size_t)std::max((trunc((xmin(2)-ld)/ld)),0.0);
        nhigh(0) = (size_t)std::min((ceil ((xmax(0)+ld)/ld)),(double)LBMDOM.Ndim(0)-1.0);
        nhigh(1) = (size_t)std::min((ceil ((xmax(1)+ld)/ld)),(double)LBMDOM.Ndim(1)-1.0);
        nhigh(2) = (size_t)std::min((ceil ((xmax(2)+ld)/ld)),(double)LBMDOM.Ndim(2)-1.0);
        if (xmax(0)<0.0) nhigh(0) = 0;
        if (xmax(1)<0.0) nhigh(1) = 0;
        if (xmax(2)<0.0) nhigh(2) = 0;
        for (size_t nx=nlow(0); nx<=nhigh(0); nx++)
        for (size_t ny=nlow(1); ny<=nhigh(1); ny++)
        for (size_t nz=nlow(2); nz<=nhigh(2); nz++)
        {
            Vec3_t xn(nx*ld,ny*ld,nz*ld);
            //Vec3_t rst;
            if (InsideCell(xn,MPMDOM.BodyMesh->Cells[ip]))
            {
                LBMDOM.Gamma [nx][ny][nz] = 1.0;
                LBMDOM.Inside[nx][ny][nz] = true;
                //Vec3_t vin = v0 + rst(0)*(v1-v0) + rst(1)*(v2-v0) + rst(2)*(v3-v0);
                //for (size_t k = 0; k < LBMDOM.Nneigh; k++)
                //{
                    //LBMDOM.Omeis[nx][ny][nz][k] = LBMDOM.F[0][nx][ny][nz][LBMDOM.Op[k]]-LBMDOM.F[0][nx][ny][nz][k]+6.0*LBMDOM.Rho[0][nx][ny][nz]*LBMDOM.W[k]*dot(LBMDOM.C[k],vin)/Cs;
                    //LBMDOM.F[0][nx][ny][nz][k] = LBMDOM.Feq(k,LBMDOM.Rho[0][nx][ny][nz],Vp);
                //}
            }
        }
    }
    //*/
    ///*
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nx=0; nx < LBMDOM.Ndim(0); nx++)
    for (size_t ny=0; ny < LBMDOM.Ndim(1); ny++)
    for (size_t nz=0; nz < LBMDOM.Ndim(2); nz++)
    {
        if (LBMDOM.Gamma[nx][ny][nz]<0.5&&LBMDOM.Inside[nx][ny][nz])
        {
            for (size_t k=0; k<LBMDOM.Nneigh; k++)
            {
                LBMDOM.F[0][nx][ny][nz][k] = 0.0;
            }
            size_t naem = 0;
            for (size_t k=0; k<LBMDOM.Nneigh; k++)
            {
                size_t nxf = (int) nx +   LBMDOM.C[k](0);                    
                size_t nyf = (int) ny +   LBMDOM.C[k](1);                    
                size_t nzf = (int) nz +   LBMDOM.C[k](2);                    
                size_t nxff= (int) nx + 2*LBMDOM.C[k](0);                    
                size_t nyff= (int) ny + 2*LBMDOM.C[k](1);                    
                size_t nzff= (int) nz + 2*LBMDOM.C[k](2);                    
                if (nxf>=LBMDOM.Ndim(0)||nyf>=LBMDOM.Ndim(1)||nzf>=LBMDOM.Ndim(2))
                {
                    continue;
                }
                if (nxff>=LBMDOM.Ndim(0)||nyff>=LBMDOM.Ndim(1)||nzff>=LBMDOM.Ndim(2))
                {
                    nxff = nxf;
                    nyff = nyf;
                    nzff = nzf;
                }
                if (LBMDOM.Gamma[nxf][nyf][nzf]<0.5&&!LBMDOM.Inside[nxf][nyf][nzf])
                {
                    for (size_t kt=0;kt<LBMDOM.Nneigh;kt++)
                    {
                        LBMDOM.F[0][nx][ny][nz][kt] += 2.0*LBMDOM.F[0][nxf][nyf][nzf][kt]-LBMDOM.F[0][nxff][nyff][nzff][kt];
                    }
                    naem ++;
                }
            }
            LBMDOM.Vel[0][nx][ny][nz] = OrthoSys::O;
            LBMDOM.Rho[0][nx][ny][nz] = 0.0;
            for (size_t k=0; k<LBMDOM.Nneigh; k++)
            {
                //LBMDOM.F[0][nx][ny][nz][k] /= naem;
                //if (LBMDOM.F[0][nx][ny][nz][k]<0.0)
                //{
                    //std::cout << "Refilling algorithm gives nagative value" << std::endl;
                    //std::cout << "F   " << LBMDOM.F[0][nx][ny][nz][k]  << std::endl;
                    //std::cout << "Pos " << iVec3_t(nx,ny,nz) << std::endl;
                    //std::cout << "k   " << k << std::endl;
                    //std::cout << "vaem " << vaem << std::endl;
                    //std::cout << "raem " << raem << std::endl;
                    //throw new Fatal("Refilling algorithm gives nagative value");
                //}
                LBMDOM.F[0][nx][ny][nz][k] = fabs(LBMDOM.F[0][nx][ny][nz][k])/naem;
                LBMDOM.Rho[0][nx][ny][nz] += LBMDOM.F[0][nx][ny][nz][k];
                LBMDOM.Vel[0][nx][ny][nz] += LBMDOM.F[0][nx][ny][nz][k]*LBMDOM.C[k];
            }
            LBMDOM.Vel[0][nx][ny][nz] *= LBMDOM.Cs/LBMDOM.Rho[0][nx][ny][nz];
        }
    }

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nx=0; nx < LBMDOM.Ndim(0); nx++)
    for (size_t ny=0; ny < LBMDOM.Ndim(1); ny++)
    for (size_t nz=0; nz < LBMDOM.Ndim(2); nz++)
    {
        if (LBMDOM.Gamma[nx][ny][nz]<0.5&&LBMDOM.Inside[nx][ny][nz])
        {
            LBMDOM.Inside[nx][ny][nz] = false;
        }
    }
    //*/
    ///*
        //std::cout << "1" << std::endl;
    //Executing the IBB forcing scheme
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ibf =0; ibf < MPMDOM.BFac.Size(); ibf++)
    {
        size_t i0 = MPMDOM.BFac[ibf].Ce->V[MPMDOM.BFac[ibf].Facets[0]]->ID;
        size_t i1 = MPMDOM.BFac[ibf].Ce->V[MPMDOM.BFac[ibf].Facets[1]]->ID;
        size_t i2 = MPMDOM.BFac[ibf].Ce->V[MPMDOM.BFac[ibf].Facets[2]]->ID;
        Vec3_t x0 = MPMDOM.BodyMesh->Verts[i0]->C;
        Vec3_t x1 = MPMDOM.BodyMesh->Verts[i1]->C;
        Vec3_t x2 = MPMDOM.BodyMesh->Verts[i2]->C;
        //Vec3_t v0 = MPMDOM.Particles[i0+MPMDOM.Meshidx]->v;
        //Vec3_t v1 = MPMDOM.Particles[i1+MPMDOM.Meshidx]->v;
        //Vec3_t v2 = MPMDOM.Particles[i2+MPMDOM.Meshidx]->v;
        Vec3_t v0 = MPMDOM.Corners[i0]->v;
        Vec3_t v1 = MPMDOM.Corners[i1]->v;
        Vec3_t v2 = MPMDOM.Corners[i2]->v;
        iVec3_t nlow,nhigh;
        Vec3_t xmin = MPMDOM.BFac[ibf].Xmin;
        Vec3_t xmax = MPMDOM.BFac[ibf].Xmax;
        nlow (0) = (size_t)std::max((trunc((xmin(0)-ld)/ld)),0.0);
        nlow (1) = (size_t)std::max((trunc((xmin(1)-ld)/ld)),0.0);
        nlow (2) = (size_t)std::max((trunc((xmin(2)-ld)/ld)),0.0);
        nhigh(0) = (size_t)std::min((ceil ((xmax(0)+ld)/ld)),(double)LBMDOM.Ndim(0)-1.0);
        nhigh(1) = (size_t)std::min((ceil ((xmax(1)+ld)/ld)),(double)LBMDOM.Ndim(1)-1.0);
        nhigh(2) = (size_t)std::min((ceil ((xmax(2)+ld)/ld)),(double)LBMDOM.Ndim(2)-1.0);
        if (xmax(0)<0.0) nhigh(0) = 0;
        if (xmax(1)<0.0) nhigh(1) = 0;
        if (xmax(2)<0.0) nhigh(2) = 0;
        //std::cout << "2" << " " << nlow << " " << nhigh << " " << xmin << " " << xmax << std::endl;
        //std::cout << ibf << " " << x0 << " " << x1 << " " << x2 << std::endl;
        //std::cout << nlow << " " << nhigh << std::endl;
        for (size_t nx=nlow(0); nx<=nhigh(0); nx++)
        for (size_t ny=nlow(1); ny<=nhigh(1); ny++)
        for (size_t nz=nlow(2); nz<=nhigh(2); nz++)
        {
            Vec3_t xn(nx*ld,ny*ld,nz*ld);
            if (LBMDOM.IsSolid[0][nx][ny][nz]) continue;                
            //Vec3_t xn((nx+0.5)*ld,(ny+0.5)*ld,(nz+0.5)*ld);
            Vec3_t B = x0 - xn;
            for (size_t k = 1; k < LBMDOM.Nneigh; k++)
            {
        //std::cout << "3 " << k << std::endl;
                Mat3_t M;
                M(0,0) = ld*LBMDOM.C[k](0); M(0,1) = x0(0) - x1(0); M(0,2) = x0(0) - x2(0);
                M(1,0) = ld*LBMDOM.C[k](1); M(1,1) = x0(1) - x1(1); M(1,2) = x0(1) - x2(1);
                M(2,0) = ld*LBMDOM.C[k](2); M(2,1) = x0(2) - x1(2); M(2,2) = x0(2) - x2(2);
                Vec3_t rst;
                if (!SolAlt(M,B,rst,1.0e-20)) continue;
        //std::cout << "4 " << k << std::endl;
                //if (nz==30&&ny==25&&nx==21) std::cout << M << " " << B << " " << rst << std::endl;
                double r = rst(0), s = rst(1), t = rst(2);
                double tol = 0.0;
                if ((r>0.0-tol&&r<=1.0+tol)&&(s>0.0-tol)&&(t>0.0-tol)&&(s+t<1.0+tol)&&dot(xn-x0,MPMDOM.BFac[ibf].Nor)>0.0-tol)
                {
                    //if (ibf==0&&nz==30&&nx==21) std::cout << M << " " << B << " " << rst << std::endl;
                    LBMDOM.Gamma[nx][ny][nz] = 2.0;
                    size_t ko = LBMDOM.Op[k];
                    Vec3_t vw = v0 + s*(v1-v0) + t*(v2-v0);
                    size_t nxf= (int) nx + LBMDOM.C[ko](0);                    
                    size_t nyf= (int) ny + LBMDOM.C[ko](1);                    
                    size_t nzf= (int) nz + LBMDOM.C[ko](2);                    
                    size_t nxs= (int) nx + LBMDOM.C[k ](0);                    
                    size_t nys= (int) ny + LBMDOM.C[k ](1);                    
                    size_t nzs= (int) nz + LBMDOM.C[k ](2);                    
                    if (LBMDOM.Gamma[nxs][nys][nzs]<0.5)
                    {
                        //std::cout << "Incorrected solid neighbour found" << std::endl;
                        //std::cout << "Position " << iVec3_t(nx,ny,nz) << std::endl;
                        //std::cout << "Solid "    << iVec3_t(nxs,nys,nzs) << std::endl;
                        //std::cout << "Direction "<< k << " " << LBMDOM.C[k] << std::endl;
                        //std::cout << "Cell No "  << MPMDOM.BFac[ibf].Ce->ID << std::endl;
                        //throw new Fatal("LBMMPM::Imprintlattice incorrect solid neighbour found");
                        continue;
                    }
                    if (nxf>=LBMDOM.Ndim(0)||nyf>=LBMDOM.Ndim(1)||nzf>=LBMDOM.Ndim(2))
                    {
                        //std::cout << "Neighbor cells out of bounds" << std::endl;
                        //throw new Fatal("LBMMPM: Imprintlattice");
                        continue;
                    }

                    double rho = LBMDOM.Rho[0][nx][ny][nz];

                    LBMDOM.F[0][nx][ny][nz][ko] = (r*LBMDOM.F[0][nxf][nyf][nzf][ko]+(1.0-r)*LBMDOM.F[0][nx][ny][nz][k] + 
                    r*LBMDOM.F[0][nxs][nys][nzs][k]+6.0*rho*LBMDOM.W[k]*dot(LBMDOM.C[ko],vw)/Cs)/(1+r);


                    if (LBMDOM.F[0][nx][ny][nz][ko]<0.0)
                    {
                        std::cout << "ibb gave negative value" << std::endl; 
                        std::cout << "rst " << rst << std::endl; 
                        std::cout << "vw " << vw << std::endl; 
                        std::cout << "Fs " <<  LBMDOM.F[0][nxs][nys][nzs][k] << std::endl; 
                        std::cout << "F "  << LBMDOM.F[0][nx][ny][nz][ko] << std::endl; 
                        std::cout << "Ff " <<  LBMDOM.F[0][nxf][nyf][nzf][ko] << std::endl;
                        std::cout << "Time " << Time << std::endl; 
                        throw new Fatal("Domain::Imprintlattice ibb gave negative value");
                    }

                    LBMDOM.F[0][nxs][nys][nzs][k] = LBMDOM.F[0][nx][ny][nz][k];

                    Vec3_t F = ld*ld*Cs*((Cs*LBMDOM.C[k]-vw)*LBMDOM.F[0][nxs][nys][nzs][k]-(Cs*LBMDOM.C[ko]-vw)*LBMDOM.F[0][nx][ny][nz][ko]-2.0*rho*LBMDOM.W[k]*Cs*LBMDOM.C[k]);
                    //Vec3_t F = ld*ld*Cs*((Cs*LBMDOM.C[k]-vw)*LBMDOM.F[0][nxs][nys][nzs][k]-(Cs*LBMDOM.C[ko]-vw)*LBMDOM.F[0][nx][ny][nz][ko]);
                    double f = norm(F);
                    if (isnan(f))
                    {
                        std::cout << "nan value detected" << std::endl;
                        std::cout << "F  = " << F  <<std::endl;
                        throw new Fatal("Domain::Imprintlattice nan value found");
                    }
                    if (f>1.0e-12)
                    {
                        Vec3_t F0 = (1.0-s-t)*F;
                        Vec3_t F1 = s*F;
                        Vec3_t F2 = t*F;
                        
                        //omp_set_lock(&MPMDOM.Particles[i0+MPMDOM.Meshidx]->lck);
                        //MPMDOM.Particles[i0+MPMDOM.Meshidx]->h += F0;
                        //omp_unset_lock(&MPMDOM.Particles[i0+MPMDOM.Meshidx]->lck);
                        //omp_set_lock(&MPMDOM.Particles[i1+MPMDOM.Meshidx]->lck);
                        //MPMDOM.Particles[i1+MPMDOM.Meshidx]->h += F1;
                        //omp_unset_lock(&MPMDOM.Particles[i1+MPMDOM.Meshidx]->lck);
                        //omp_set_lock(&MPMDOM.Particles[i2+MPMDOM.Meshidx]->lck);
                        //MPMDOM.Particles[i2+MPMDOM.Meshidx]->h += F2;
                        //omp_unset_lock(&MPMDOM.Particles[i2+MPMDOM.Meshidx]->lck);

                        omp_set_lock(&MPMDOM.Corners[i0]->lck);
                        MPMDOM.Corners[i0]->h += F0;
                        omp_unset_lock(&MPMDOM.Corners[i0]->lck);
                        omp_set_lock(&MPMDOM.Corners[i1]->lck);
                        MPMDOM.Corners[i1]->h += F1;
                        omp_unset_lock(&MPMDOM.Corners[i1]->lck);
                        omp_set_lock(&MPMDOM.Corners[i2]->lck);
                        MPMDOM.Corners[i2]->h += F2;
                        omp_unset_lock(&MPMDOM.Corners[i2]->lck);
                    }
                }
            }
        }
    }
}

void Domain::Reset()
{
    //std::cout << "0" << std::endl;
    //std::cout << LBMDOM.Ndim << " " << Nproc<< " " << LBMDOM.Nneigh << std::endl;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ix=0; ix < LBMDOM.Ndim(0); ix++)
    for (size_t iy=0; iy < LBMDOM.Ndim(1); iy++)
    for (size_t iz=0; iz < LBMDOM.Ndim(2); iz++)
    {
        //std::cout << iVec3_t(ix,iy,iz) << std::endl;
        //LBMDOM.Gamma[ix][iy][iz] = (double) LBMDOM.IsSolid[0][ix][iy][iz];
        LBMDOM.Gamma[ix][iy][iz] = 0.0;
        //for (size_t k=0; k<LBMDOM.Nneigh; k++)
        //{
            //LBMDOM.Omeis[ix][iy][iz][k] = 0.0;
        //}
    }

    //std::cout << "1" << std::endl;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ic=0; ic<MPMDOM.Corners.Size(); ic++)
    {
        MPMDOM.Corners[ic]->h = OrthoSys::O;
    }

    //std::cout << "2" << std::endl;
    MPMDOM.UpdateMesh();
}

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{
    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Nproc = TheNproc;
    LBMDOM.Nproc = TheNproc;
    MPMDOM.Nproc = TheNproc;
    MPMDOM.Dt    = dt;
    MTD = new MtData[Nproc];
    
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1    , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                    , TERM_RST);
    printf("%s  C parameter                      =  %g%s\n"       ,TERM_CLR2, LBMDOM.Cs                             , TERM_RST);
    for (size_t i=0;i<LBMDOM.Nl;i++)
    {
    printf("%s  Tau of Lattice %zd                 =  %g%s\n"       ,TERM_CLR2, i, LBMDOM.Tau[i]                    , TERM_RST);
    }


    //Initializing MPM particles
    MPMDOM.Mmin = MPMDOM.Particles[0]->m;
    for (size_t ip=0;ip<MPMDOM.Particles.Size();ip++)
    {
        MPMDOM.Particles[ip]->xb = MPMDOM.Particles[ip]->x - MPMDOM.Particles[ip]->v*dt;
        if (MPMDOM.Mmin > MPMDOM.Particles[ip]->m) MPMDOM.Mmin = MPMDOM.Particles[ip]->m;
    }
    MPMDOM.Mmin *= 5.0e-12;

    //std::cout << "0" << std::endl;
    MPMDOM.BoundaryMesh();
    Reset();
    ImprintLattice();
    double tout = Time;

    //std::cout << "1" << std::endl;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time >= tout)
        {
            if (TheFileKey!=NULL)
            {
                String fnmpm;
                fnmpm.Printf    ("%s_mpm_%04d", TheFileKey, idx_out);
                String fnlbm;
                fnlbm.Printf    ("%s_lbm_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    //WriteXDMF(fn.CStr());
                    MPMDOM.WriteXDMF(fnmpm.CStr());
                    //LBMDOM.WriteXDMF(fnlbm.CStr());
                    LBMDOM.WriteXDMF_DEM(fnlbm.CStr());
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }
        //std::cout << "1" << std::endl;
        //Coupling between LBM and MPM
        Reset();
        ImprintLattice();

        //std::cout << "2" << std::endl;
        //The LBM dynamics
        LBMDOM.VelDen();
        LBMDOM.CollideMPM();
        //LBMDOM.CollideMRT();
        LBMDOM.StreamMPM();

        //std::cout << "3" << std::endl;
        //The MPM dynamics
        MPMDOM.OneStepCPI();
        //MPMDOM.OneStepUSF();
        //MPMDOM.ParticleToNode();
        //MPMDOM.NodeToParticle();
        //std::cout << "4" << std::endl;
        //Reset domain
        

        Time += dt;
        LBMDOM.Time += dt;
        MPMDOM.Time += dt;
        //std::cout << Time << std::endl;
    }
}
}
#endif
