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

#ifndef MECHSYS_EMLBM_LATTICE_H
#define MECHSYS_EMLBM_LATTICE_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// MechSys
#include <mechsys/emlbm2/Cell.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>

class Lattice
{
public:
    //typedefs
    typedef void (*ptFun_t) (Lattice & Lat, void * UserData);

    //Constructors
    Lattice () {};            //Default
    Lattice (iVec3_t Ndim, double dx, double dt);

    //Methods
    void Stream1    (size_t Np);                                      ///< Stream the velocity distributions
    void Stream2    (size_t Np);                                      ///< Stream the velocity distributions
    void CalcField  (size_t Np);                                      ///< Calculate the Electric and Magnetic Fields
    Cell * GetCell(iVec3_t const & v);                     ///< Get pointer to cell at v


    //Data
    iVec3_t                                   Ndim;             // Dimensions of the lattice
    size_t                                    idx_out;          // The discrete time step
    size_t                                    Ncells;           // Number of cells per lattice
    double                                    Time;             // The current time
    double                                    dx;               // grid space
    double                                    dt;               // time step
    Cell                                   ** Cells;            // Array of pointer cells
};

inline Lattice::Lattice(iVec3_t TheNdim, double Thedx, double Thedt)
{
    Ndim = TheNdim;
    dx   = Thedx;
    dt   = Thedt;

    //Cells.Resize(Ndim[0]*Ndim[1]*Ndim[2]);
    Cells = new Cell * [Ndim[0]*Ndim[1]*Ndim[2]];
    Ncells = Ndim[0]*Ndim[1]*Ndim[2];
    size_t n = 0;
    for (size_t k=0;k<Ndim[2];k++)
    for (size_t j=0;j<Ndim[1];j++)
    for (size_t i=0;i<Ndim[0];i++)
    {
        //Cells[n] =  new Cell(n,TheMethod,iVec3_t(i,j,k),Ndim,dx/dt,Tau);
        Cells[n] = new Cell(n,iVec3_t(i,j,k),Ndim,dx/dt,dt);
        n++;
    } 
}

inline void Lattice::Stream1(size_t Np)
{
    // Assign temporal distributions
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Ncells;i++)
    for (size_t j=0;j<Cells[i]->Nneigh;j++)
    {
        for (size_t k=0;k<2;k++)
        {
            Cells[Cells[i]->Neighs[j]]->FEtemp[k][j] = Cells[i]->FE[k][j];
            Cells[Cells[i]->Neighs[j]]->FBtemp[k][j] = Cells[i]->FB[k][j];
        }
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
        double ** FEtemp  = Cells[i]->FE;
        Cells[i]->FE      = Cells[i]->FEtemp;
        Cells[i]->FEtemp  = FEtemp;
        double ** FBtemp  = Cells[i]->FB;
        Cells[i]->FB      = Cells[i]->FBtemp;
        Cells[i]->FBtemp  = FBtemp;
    }
}

inline void Lattice::CalcField(size_t Np)
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
