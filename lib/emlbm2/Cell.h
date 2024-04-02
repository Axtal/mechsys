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



class Cell
{
public:
	static const Vec3_t  C   [13]; ///< Local velocities (D3Q13)
	static const Vec3_t  D1  [12]; ///< Auxiliary vectors for electric field
	static const Vec3_t  D2  [12]; ///< Auxiliary vectors for electric field
	static const Vec3_t  H1  [12]; ///< Auxiliary vectors for magnetic field
	static const Vec3_t  H2  [12]; ///< Auxiliary vectors for magnetic field
   
    //Constructor
    Cell (size_t ID, iVec3_t Indexes, iVec3_t Ndim, double Cs, double Dt); ///< Constructor, it receives the grid type, ht einteger position and the total integer dimension of the domain and the spatial and time steps
    
    // Methods
    double       FEeq(size_t mu, size_t k);                         ///< Calculate the equilibrium distribution function FE
    double       FBeq(size_t mu, size_t k);                         ///< Calculate the equilibrium distribution function FB
    void         CalcProp();                                       ///< Calculate the vectorial properties with the new distributions functions
    void         Initialize(double TheRho, Vec3_t & TheJ, Vec3_t & TheE, Vec3_t & TheB);                                     ///< Initialize cell with a given velocity and density

    // Data
    size_t       Nneigh;   ///< Number of neighbors
    size_t       ID;       ///< Tag for the particle
    double       Eps;      ///< Dielectric constant
    double       Mu;       ///< Magnetic permitivity
    double       Sig;      ///< Electric conductivity
    double       Dt;       ///< Time step
    
    iVec3_t      Index;    ///< Vector of indexes

    double        *      F0; ///< Distribution functions for vector potentials
    double        *  F0temp; ///< Temporary distribution functions
    double        **     FE; ///< Distribution functions for vector potentials
    double        ** FEtemp; ///< Temporary distribution functions
    double        **     FB; ///< Distribution functions for vector potentials
    double        ** FBtemp; ///< Temporary distribution functions
    size_t        *  Neighs; ///< Array of neighbors indexes
    double              Rho; ///< Charge Density
    double             Rhof; ///< Charge Density of free charges
    Vec3_t                J; ///< Current density
    Vec3_t               Jf; ///< Current free density
    Vec3_t                E; ///< Electric Field
    Vec3_t                B; ///< Magnetic Field

};

inline Cell::Cell(size_t TheID, iVec3_t TheIndexes, iVec3_t TheNdim, double TheCs, double TheDt)
{
    ID      = TheID;
    Index   = TheIndexes;
    Dt      = TheDt;
    Eps     = 1.0;
    Mu      = 1.0;
    Sig     = 0.0;
    Nneigh  = 12;
    Rhof    = 0.0;
    Jf      = OrthoSys::O;
    F0      = new double  [2];
    F0temp  = new double  [2];
    FE      = new double *[2];
    FEtemp  = new double *[2];
    FB      = new double *[2];
    FBtemp  = new double *[2];
    for (size_t i=0;i<2;i++)
    {
        FE    [i]  = new double [Nneigh];
        FEtemp[i]  = new double [Nneigh];
        FB    [i]  = new double [Nneigh];
        FBtemp[i]  = new double [Nneigh];
    }
    
    Initialize(0.0,OrthoSys::O,OrthoSys::O,OrthoSys::O);


    Neighs  = new size_t [Nneigh];

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

}


inline void Cell::CalcProp()
{
    E   = OrthoSys::O;
    B   = OrthoSys::O;
    Rho = F0[0];
    for (size_t i=0;i<Nneigh;i++)
    {
        E   += FE[0][i]*D1[i] + FE[1][i]*D2[i];
        B   += FB[0][i]*H1[i] + FB[1][i]*H2[i];
        Rho += FE[0][i] + FE[1][i];
    }
    //E /= sqrt(2.0)*Eps;
    E /= Eps;
    J  = Sig*E/(1.0+0.25*Sig/Eps);
    //J  = Sig*E;
}

inline double Cell::FEeq(size_t mu,size_t k)
{
    if (mu==0) return (1.0/16.0)*dot(C[k],J+Jf)+(Eps/4.0)*dot(E-1.0/(4.0*Eps)*(J+Jf),D1[k])+(1.0/(8.0*Mu))*dot(B,H1[k]);
    else       return (1.0/16.0)*dot(C[k],J+Jf)+(Eps/4.0)*dot(E-1.0/(4.0*Eps)*(J+Jf),D2[k])+(1.0/(8.0*Mu))*dot(B,H2[k]);
}

inline double Cell::FBeq(size_t mu,size_t k)
{
    if (mu==0) return (1.0/16.0)*dot(C[k],J+Jf)+(1.0/4.0)*dot(E-1.0/(4.0*Eps)*(J+Jf),D1[k])+(1.0/8.0)*dot(B,H1[k]);
    else       return (1.0/16.0)*dot(C[k],J+Jf)+(1.0/4.0)*dot(E-1.0/(4.0*Eps)*(J+Jf),D2[k])+(1.0/8.0)*dot(B,H2[k]);
}

inline void Cell::Initialize(double TheRho, Vec3_t & TheJ, Vec3_t & TheE, Vec3_t & TheB)
{
    Rho = TheRho;
    J   = TheJ;
    E   = TheE;
    B   = TheB;

    F0[0] = F0[1] = TheRho;
    //F0[0] = F0[1] = Rhof;
    for (size_t i=0;i<Nneigh;i++)
    {
        FE[0][i] = FEeq(0,i);
        FE[1][i] = FEeq(1,i);
        FB[0][i] = FBeq(0,i);
        FB[1][i] = FBeq(1,i);
    }
}



const Vec3_t Cell::C   [13] = { { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0}, {-1.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0, 1.0}, { 0.0,-1.0, 1.0}, { 0.0,-1.0,-1.0}, { 0.0, 1.0,-1.0}, {0.0,0.0,0.0} };
const Vec3_t Cell::D1  [12] = { {-0.5, 0.5, 0.0}, {-0.5,-0.5, 0.0}, { 0.5,-0.5, 0.0}, { 0.5, 0.5, 0.0}, {-0.5, 0.0, 0.5}, {-0.5, 0.0,-0.5}, { 0.5, 0.0,-0.5}, { 0.5, 0.0, 0.5}, { 0.0,-0.5, 0.5}, { 0.0,-0.5,-0.5}, { 0.0, 0.5,-0.5}, { 0.0, 0.5, 0.5} }; 
const Vec3_t Cell::D2  [12] = { { 0.5,-0.5, 0.0}, { 0.5, 0.5, 0.0}, {-0.5, 0.5, 0.0}, {-0.5,-0.5, 0.0}, { 0.5, 0.0,-0.5}, { 0.5, 0.0, 0.5}, {-0.5, 0.0, 0.5}, {-0.5, 0.0,-0.5}, { 0.0, 0.5,-0.5}, { 0.0, 0.5, 0.5}, { 0.0,-0.5, 0.5}, { 0.0,-0.5,-0.5} };
const Vec3_t Cell::H1  [12] = { { 0.0, 0.0, 1.0}, { 0.0, 0.0, 1.0}, { 0.0, 0.0, 1.0}, { 0.0, 0.0, 1.0}, { 0.0,-1.0, 0.0}, { 0.0,-1.0, 0.0}, { 0.0,-1.0, 0.0}, { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0} };
const Vec3_t Cell::H2  [12] = { { 0.0, 0.0,-1.0}, { 0.0, 0.0,-1.0}, { 0.0, 0.0,-1.0}, { 0.0, 0.0,-1.0}, { 0.0, 1.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0} };

#endif // MECHSYS_EMLBM_CELL_H
