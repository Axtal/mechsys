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
#include <mechsys/emlbm/Cell.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>

class Lattice
{
public:
    //typedefs
    typedef void (*ptFun_t) (Lattice & Lat, void * UserData);

    //Constructors
    Lattice () {};            //Default
    Lattice (LBMethod Method, double Tau, iVec3_t Ndim, double dx, double dt);

    //Methods
    void Stream1  (size_t n = 0, size_t Np = 1);                                      ///< Stream the velocity distributions
    void Stream2  (size_t n = 0, size_t Np = 1);                                      ///< Stream the velocity distributions
    void CalcField(size_t n = 0, size_t Np = 1);                                      ///< Calculate the magnetic and electric fields
    Cell * GetCell(iVec3_t const & v);                                                ///< Get pointer to cell at v


    //Data
    iVec3_t                                   Ndim;             // Dimensions of the lattice
    size_t                                    idx_out;          // The discrete time step
    size_t                                    Ncells;           // Number of cells per lattice
    double                                    Time;             // The current time
    double                                    dx;               // grid space
    double                                    dt;               // time step
    double                                    Tau;              // Relaxation time
    Cell                                   ** Cells;            // Array of pointer cells
};

inline Lattice::Lattice(LBMethod TheMethod, double TheTau, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Ndim = TheNdim;
    dx   = Thedx;
    dt   = Thedt;
    Tau  = TheTau;

    //Cells.Resize(Ndim[0]*Ndim[1]*Ndim[2]);
    Cells = new Cell * [Ndim[0]*Ndim[1]*Ndim[2]];
    Ncells = Ndim[0]*Ndim[1]*Ndim[2];
    size_t n = 0;
    for (size_t k=0;k<Ndim[2];k++)
    for (size_t j=0;j<Ndim[1];j++)
    for (size_t i=0;i<Ndim[0];i++)
    {
        //Cells[n] =  new Cell(n,TheMethod,iVec3_t(i,j,k),Ndim,dx/dt,Tau);
        Cells[n] = new Cell(n,TheMethod,iVec3_t(i,j,k),Ndim,dx/dt,dt);
        n++;
    } 
}

inline void Lattice::Stream1(size_t n, size_t Np)
{
	size_t Ni = Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Ncells : Fn = (n+1)*Ni;
    // Assign temporal distributions
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    for (size_t j=1;j<Cells[i]->Nneigh;j++)
    {
        Cells[Cells[i]->Neighs[j]]->Htemp[j] = Cells[i]->H[j];
        for (size_t k=0;k<4;k++)
        {
            Cells[Cells[i]->Neighs[j]]->Ftemp[k][j] = Cells[i]->F[k][j];
            Cells[Cells[i]->Neighs[j]]->Gtemp[k][j] = Cells[i]->G[k][j];
        }
    }
}

inline void Lattice::Stream2(size_t n, size_t Np)
{
	size_t Ni = Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Ncells : Fn = (n+1)*Ni;
    //Swap the distribution values
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        double ** Ftemp  = Cells[i]->F;
        Cells[i]->F      = Cells[i]->Ftemp;
        Cells[i]->Ftemp  = Ftemp;
        double ** Gtemp  = Cells[i]->G;
        Cells[i]->G      = Cells[i]->Gtemp;
        Cells[i]->Gtemp  = Gtemp;
        double *  Htemp  = Cells[i]->H;
        Cells[i]->H      = Cells[i]->Htemp;
        Cells[i]->Htemp  = Htemp;
        Cells[i]->CalcProp();
    }
}

inline void Lattice::CalcField(size_t n, size_t Np)
{
	size_t Ni = Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Ncells : Fn = (n+1)*Ni;
    //Calculate fields
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        Cell * c = Cells[i];
        //if (c->Index[0]==Ndim[0]/2&&c->Index[1]==Ndim[1]/2&&c->Index[2]==Ndim[2]/2) std::cout << c->Index << " " << Cells[Cells[i]->Neighs[4]]->Index << " " << 
            //Cells[Cells[i]->Neighs[5]]->Index << std::endl;
        c->B(0) = 0.5*(Cells[c->Neighs[3]]->A[3]-Cells[c->Neighs[4]]->A[3])/dx-0.5*(Cells[c->Neighs[5]]->A[2]-Cells[c->Neighs[6]]->A[2])/dx;
        c->B(1) = 0.5*(Cells[c->Neighs[5]]->A[1]-Cells[c->Neighs[6]]->A[1])/dx-0.5*(Cells[c->Neighs[1]]->A[3]-Cells[c->Neighs[2]]->A[3])/dx;
        c->B(2) = 0.5*(Cells[c->Neighs[1]]->A[2]-Cells[c->Neighs[2]]->A[2])/dx-0.5*(Cells[c->Neighs[3]]->A[1]-Cells[c->Neighs[4]]->A[1])/dx;
        c->E(0) = -0.5*(Cells[c->Neighs[1]]->A[0]-Cells[c->Neighs[2]]->A[0])/dx-0.5*(3.0*c->A[1]-4.0*c->Ap[1]+c->Al[1])/dt;
        c->E(1) = -0.5*(Cells[c->Neighs[3]]->A[0]-Cells[c->Neighs[4]]->A[0])/dx-0.5*(3.0*c->A[2]-4.0*c->Ap[2]+c->Al[2])/dt;
        c->E(2) = -0.5*(Cells[c->Neighs[5]]->A[0]-Cells[c->Neighs[6]]->A[0])/dx-0.5*(3.0*c->A[3]-4.0*c->Ap[3]+c->Al[3])/dt;
        c->E   /= c->Eps;
    }

}



inline Cell * Lattice::GetCell(iVec3_t const & v)
{
    return Cells[v[0] + v[1]*Ndim[0] + v[2]*Ndim[0]*Ndim[1]];
}

#endif
