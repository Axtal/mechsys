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

#ifndef MECHSYS_EMLBM_CELL_H
#define MECHSYS_EMLBM_CELL_H

// Std lib
#ifdef USE_THREAD
    #include <pthread.h>
#endif

//OpenMP
#ifdef USE_OMP
    #include <omp.h>
#endif

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>



enum LBMethod
{
    D3Q7     ///< 3D 7 velocities
};

class Cell
{
public:
	static const double   WEIGHTSD3Q7  [7]; ///< Weights for the equilibrium distribution functions (D3Q7)
	static const Vec3_t   LVELOCD3Q7   [7]; ///< Local velocities (D3Q7)
	static const size_t   OPPOSITED3Q7 [7]; ///< Opposite directions (D3Q7)
   
    //Constructor
    Cell (size_t ID, LBMethod Method, iVec3_t Indexes, iVec3_t Ndim, double Cs, double Dt); ///< Constructor, it receives the grid type, ht einteger position and the total integer dimension of the domain and the spatial and time steps
    
    // Methods
    double       Feq(size_t mu, size_t k);                         ///< Calculate the equilibrium distribution function F
    double       Geq(size_t mu, size_t k);                         ///< Calculate the equilibrium distribution function G
    double       Heq(size_t k);                                    ///< Calculate the equilibrium distribution function H
    void         CalcProp();                                       ///< Calculate the vectorial properties with the new distributions functions
    void         Initialize();                                     ///< Initialize cell with a given velocity and density

#ifdef USE_THREAD
    pthread_mutex_t lck;
    //std::mutex mtex;       ///< to protect variables in multithreading
#endif

    // Data
    //LBMethod     Method;   ///< Is 2D, 3D and how many velocities it has
    bool         IsSolid;  ///< It is a solid node
    //double       Tau;      ///< Relaxation Time
    double       Gs;       ///< Interaction constant between solid and fluid
    size_t       Nneigh;   ///< Number of neighbors
    size_t       ID;       ///< Tag for the particle
    double       Cs;       ///< Velocity of the grid
    double       Eps;      ///< Dialectric constant
    double       Dt;       ///< Time step
    
    iVec3_t      Index;    ///< Vector of indexes

    double const  *  W;     ///< Pointer to array of weights
    Vec3_t const  *  C;     ///< Pointer to velocity constants
    double        ** F;     ///< Distribution functions for vector potentials
    double        ** Ftemp; ///< Temporary distribution functions
    double        ** G;     ///< Distribution functions for auxiliary vectors S
    double        ** Gtemp; ///< Temporary distribution functions
    double        *  H;     ///< Distribution functions for density
    double        *  Htemp; ///< Temporary distribution functions
    double        *  A;     ///< Array of potentials
    double        *  Ap;    ///< Array of previous potentials
    double        *  Al;    ///< Array of late potentials (2 time steps before)
    double        *  J;     ///< 4-current
    double        *  Sig;   ///< Array of currents
    double        ** S;     ///< Array of auxiliary vectors
    size_t        *  Neighs;///< Array of neighbors indexes
    Vec3_t           E;     ///< Electric Field
    Vec3_t           B;     ///< Magnetic Field

};

inline Cell::Cell(size_t TheID, LBMethod TheMethod, iVec3_t TheIndexes, iVec3_t TheNdim, double TheCs, double TheDt)
{
    IsSolid= false;
    ID     = TheID;
    //Method = TheMethod;
    Index  = TheIndexes;
    Cs     = TheCs;
    Dt     = TheDt;
    Gs     = 1.0;
    Eps    = 1.0;
    //Tau    = TheTau;
    if (TheMethod==D3Q7)
    {
        Nneigh = 7;
        W      = WEIGHTSD3Q7;
        C      = LVELOCD3Q7;
    }
    F      = new double *[4];
    Ftemp  = new double *[4];
    G      = new double *[4];
    Gtemp  = new double *[4];
    A      = new double  [4];
    Ap     = new double  [4];
    Al     = new double  [4];
    Sig    = new double  [4];
    J      = new double  [4];
    S      = new double *[4];
    for (size_t i=0;i<4;i++)
    {
        F    [i]  = new double [Nneigh];
        Ftemp[i]  = new double [Nneigh];
        G    [i]  = new double [Nneigh];
        Gtemp[i]  = new double [Nneigh];
        S    [i]  = new double [3];
    }
    H      = new double [Nneigh];
    Htemp  = new double [Nneigh];
    Neighs = new size_t [Nneigh];

    //Set neighbors
    for (size_t k=0;k<Nneigh;k++)
    {
        //iVec3_t nindex = Index + C[k];
        blitz::TinyVector<int,3>   nindex;
        nindex[0] = static_cast<int>(Index[0])+ C[k][0];
        nindex[1] = static_cast<int>(Index[1])+ C[k][1];
        nindex[2] = static_cast<int>(Index[2])+ C[k][2];
        if (nindex[0]==                          -1) nindex[0] = TheNdim[0]-1;
        if (nindex[0]==static_cast<int>(TheNdim[0])) nindex[0] = 0;
        if (nindex[1]==                          -1) nindex[1] = TheNdim[1]-1;
        if (nindex[1]==static_cast<int>(TheNdim[1])) nindex[1] = 0;
        if (nindex[2]==                          -1) nindex[2] = TheNdim[2]-1;
        if (nindex[2]==static_cast<int>(TheNdim[2])) nindex[2] = 0;

        Neighs[k] =  nindex[0] + nindex[1]*TheNdim[0] + nindex[2]*TheNdim[0]*TheNdim[1];
    }

#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#endif
}


inline void Cell::CalcProp()
{
    for (size_t i=0;i<4;i++)
    {
        Al  [i] = Ap[i];
        Ap  [i] =  A[i];
        A   [i] = 0.0;
        Sig [i] = 0.0;
        S[i][0] = 0.0;
        S[i][1] = 0.0;
        S[i][2] = 0.0;
        for (size_t j=0;j<Nneigh;j++)
        {
            A[i]    += F[i][j];
            S[i][0] += F[i][j]*C[j][0];
            S[i][1] += F[i][j]*C[j][1];
            S[i][2] += F[i][j]*C[j][2];
            Sig [i] += G[i][j];
        }
        //A  [i] = A  [i] + 0.5*Dt*Sig[i];
        Sig[i] = Sig[i] + 0.5*Dt*J[i];
        A  [i] = A  [i] + 0.5*Dt*Sig[i];
    }
    J[0] = 0.0;
    for (size_t j=0;j<Nneigh;j++)
    {
        J[0] += H[j];
    }
    J[0] *= Cs*Cs/Eps;
    //J[0] *= Cs*Cs;
}

inline double Cell::Feq(size_t mu,size_t k)
{
    if (k==0) return 4.0*W[k]*(1.0-3.0*Cs*Cs/Eps)*A[mu];
    else      return 4.0*W[k]*(Cs*Cs/Eps*A[mu] + C[k][0]*S[mu][0] + C[k][1]*S[mu][1] + C[k][2]*S[mu][2]);
}
inline double Cell::Geq(size_t mu,size_t k)
{
    return W[k]*Sig[mu];
}
inline double Cell::Heq(size_t k)
{
    return W[k]*(J[0]/(Cs*Cs/Eps)+4.0*(C[k][0]*J[1] + C[k][1]*J[2] + C[k][2]*J[3]));
    //return W[k]*(J[0]/(Cs*Cs)+4.0*(C[k][0]*J[1] + C[k][1]*J[2] + C[k][2]*J[3]));
}

inline void Cell::Initialize()
{
    for (size_t i=0;i<4;i++)
    {
        for (size_t j=0;j<Nneigh;j++)
        {
            F    [i][j] = 0.0;
            Ftemp[i][j] = 0.0;
            G    [i][j] = -W[j]*J[i];
            Gtemp[i][j] = -W[j]*J[i];
        }
    }
    for (size_t j=0;j<Nneigh;j++)
    {
        H    [j] = 0.0;
        Htemp[j] = 0.0;
    }
}



const double Cell::WEIGHTSD3Q7   [ 7] = { 1./4., 1./8., 1./8., 1./8., 1./8, 1./8., 1./8 };
const Vec3_t Cell::LVELOCD3Q7    [ 7] = { {0,0,0}, {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };

#endif // MECHSYS_EMLBM_CELL_H
