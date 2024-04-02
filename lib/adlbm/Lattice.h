/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2014 Sergio Galindo                                    *
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

#ifndef MECHSYS_ADLBM_LATTICE_H
#define MECHSYS_ADLBM_LATTICE_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// MechSys
#include <mechsys/adlbm/Cell.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>

class Lattice
{
public:
    //typedefs
    typedef void (*ptFun_t) (Lattice & Lat, void * UserData);

    //Constructors
    Lattice () {};            //Default
    Lattice (LBMethod Method, double Thenu, double TheDif, iVec3_t TheNdim, double Thedx, double Thedt);

    //Methods
    void Collide    (size_t Np = 1);                                  ///< Collision operation                                              ///< Apply the interaction forces and the collision operator
    void Stream1    (size_t Np = 1);                                  ///< Stream the velocity distributions
    void Stream2    (size_t Np = 1);                                  ///< Stream the velocity distributions
    void CalcProps  (size_t Np = 1);                                  ///< Calculate the fluid properties
    Cell * GetCell(iVec3_t const & v);                                ///< Get pointer to cell at v


    //Data
    iVec3_t                                         Ndim;             ///< Dimensions of the lattice
    size_t                                          idx_out;          ///< The discrete time step
    size_t                                          Ncells;           ///< Number of cells per lattice
    double                                          Time;             ///< The current time
    double                                          dx;               ///< grid space
    double                                          dt;               ///< time step
    double                                          Nu;               ///< Real viscosity
    double                                          Dif;              ///< Diffusion coefficient
    double                                          Tau;              ///< Relaxation time
    double                                          Tauc;             ///< Relaxation time for diffusion 
    Cell                                         ** Cells;            ///< Array of pointer cells
    Array <double>                                  EEk;              ///< Diadic velocity tensor trace
    double                                          Sc;               ///< Smagorinsky constant
};

inline Lattice::Lattice(LBMethod TheMethod, double Thenu, double TheDif, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Ndim = TheNdim;
    dx   = Thedx;
    dt   = Thedt;
    Nu   = Thenu;
    Dif  = TheDif;
    Tau  = 3.0*Nu *dt/(dx*dx) + 0.5;
    Tauc = 3.0*Dif*dt/(dx*dx) + 0.5;


    Cells = new Cell * [Ndim[0]*Ndim[1]*Ndim[2]];
    Ncells = Ndim[0]*Ndim[1]*Ndim[2];
    size_t n = 0;
    for (size_t k=0;k<Ndim[2];k++)
    for (size_t j=0;j<Ndim[1];j++)
    for (size_t i=0;i<Ndim[0];i++)
    {
        Cells[n] = new Cell(n,TheMethod,iVec3_t(i,j,k),Ndim,dx/dt,dt);
        Cells[n]->Dif  = Dif;
        Cells[n]->Tauc = 3.0*Dif*dt/(dx*dx) + 0.5;
        n++;
    } 
    Cells[0]->Cs = dx/dt;

    EEk.Resize(Cells[0]->Nneigh);
    for (size_t k=0;k<Cells[0]->Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(Cells[0]->C[k][n]*Cells[0]->C[k][m]);
        }
    }
    Sc = 0.17;
}

inline void Lattice::Collide(size_t Np)
{
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Ncells;i++)
    {
        Cell * c = Cells[i];
        Array<double> NonEq(c->Nneigh);
        double TauS = Tau;
        double Q = 0.0;
        for (size_t k=0;k<c->Nneigh;k++)
        {
            double FDeqn = c->Feq(k);
            NonEq[k] = c->F[k] - FDeqn;
            Q += NonEq[k]*NonEq[k]*EEk[k];
        }
        Q = sqrt(2.0*Q);
        TauS = 0.5*(TauS + sqrt(TauS*TauS + 6.0*Q*Sc/c->Rho));
        for (size_t k=0;k<c->Nneigh;k++)
        {
            if (!c->IsSolid)
            {
                c->Ftemp[k] = c->F[k] - (NonEq[k])/TauS;
            }
            else
            {
                c->Ftemp[k] = c->F[c->Op[k]];
            }
            if (!c->IsNodif)
            {
                c->Gtemp[k] = c->G[k] - (c->G[k] - c->Geq(k))/c->Tauc;
            }
            else
            {
                c->Gtemp[k] = c->G[c->Op[k]];
            }
        }
        for (size_t k=0;k<c->Nneigh;k++)
        {
            c->F[k] = c->Ftemp[k];
            c->G[k] = c->Gtemp[k];
        }
    }   
}

inline void Lattice::Stream1(size_t Np)
{
    // Assign temporal distributions
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Ncells;i++)
    for (size_t j=1;j<Cell::Nneigh;j++)
    {
        Cells[Cells[i]->Neighs[j]]->Ftemp[j] = Cells[i]->F[j];
        Cells[Cells[i]->Neighs[j]]->Gtemp[j] = Cells[i]->G[j];
    }
}

inline void Lattice::Stream2(size_t Np)
{
    //Swap the distribution values
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Ncells;i++)
    {
        double * Ftemp  = Cells[i]->F;
        Cells[i]->F     = Cells[i]->Ftemp;
        Cells[i]->Ftemp = Ftemp;
        double * Gtemp  = Cells[i]->G;
        Cells[i]->G     = Cells[i]->Gtemp;
        Cells[i]->Gtemp = Gtemp;
    }
}

inline void Lattice::CalcProps(size_t Np)
{
    //Calculate fields and densities
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Ncells;i++)
    {
        Cells[i]->CalcProp();
    }
}



inline Cell * Lattice::GetCell(iVec3_t const & v)
{
    return Cells[v[0] + v[1]*Ndim[0] + v[2]*Ndim[0]*Ndim[1]];
}

#endif
