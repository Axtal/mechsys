/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_LBM_CELL_H
#define MECHSYS_LBM_CELL_H

// Std lib
#ifdef USE_OMP
    #include <omp.h>
#endif

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>

enum LBMethod
{
    D2Q5,     ///< 2D 5 velocities
    D2Q9,     ///< 2D 9 velocities
    D3Q15,    ///< 3D 15 velocities
    D3Q19,    ///< 3D 19 velocities
    D3Q27     ///< 3D 27 velocities
};

class Cell
{
public:
	static const double   WEIGHTSD2Q5   [ 5]; ///< Weights for the equilibrium distribution functions (D2Q5)
	static const double   WEIGHTSD2Q9   [ 9]; ///< Weights for the equilibrium distribution functions (D2Q9)
	static const double   WEIGHTSD3Q15  [15]; ///< Weights for the equilibrium distribution functions (D3Q15)
	static const double   WEIGHTSD3Q19  [19]; ///< Weights for the equilibrium distribution functions (D3Q19)
	//static const double   WEIGHTSD3Q27  [27]; ///< Weights for the equilibrium distribution functions (D3Q27)
	static const Vec3_t   LVELOCD2Q5    [ 5]; ///< Local velocities (D2Q5) 
	static const Vec3_t   LVELOCD2Q9    [ 9]; ///< Local velocities (D2Q9) 
	static const Vec3_t   LVELOCD3Q15   [15]; ///< Local velocities (D3Q15)
	static const Vec3_t   LVELOCD3Q19   [19]; ///< Local velocities (D3Q19)
	//static const Vec3_t   LVELOCD3Q27   [27]; ///< Local velocities (D3Q27)
	static const size_t   OPPOSITED2Q5  [ 5]; ///< Opposite directions (D2Q5) 
	static const size_t   OPPOSITED2Q9  [ 9]; ///< Opposite directions (D2Q9) 
	static const size_t   OPPOSITED3Q15 [15]; ///< Opposite directions (D3Q15)
	static const size_t   OPPOSITED3Q19 [19]; ///< Opposite directions (D3Q19)
	//static const size_t   OPPOSITED3Q27 [27]; ///< Opposite directions (D3Q27)
   
    //Constructor
    Cell (size_t ID, LBMethod Method, iVec3_t Indexes, iVec3_t Ndim, double Cs, double Tau); ///< Constructor, it receives the grid type, ht einteger position and the total integer dimension of the domain and the spatial and time steps
    
    // Methods
    double       Density();                                        ///< Calculate density
    void         Velocity(Vec3_t & V);                             ///< Calculate velocity
    double       VelDen(Vec3_t & V);                               ///< Calculate both (faster)
    double       Feq(size_t k, Vec3_t const & V, double Rho);      ///< Calculate the equilibrium distribution function
    void         Initialize(double Rho, Vec3_t const & V);         ///< Initialize cell with a given velocity and density
    void         BounceBack();                                     ///< Apply the bounceback rule to this cell

#ifdef USE_OMP
    omp_lock_t      lck;             ///< to protect variables in multithreading
#endif

    // Data
    //LBMethod     Method;   ///< Is 2D, 3D and how many velocities it has
    bool         IsSolid;  ///< It is a solid node
    double       Tau;      ///< Relaxation Time
    double       Gamma;    ///< Solid/Fluid ratio
    double       Gammap;   ///< Solid/Fluid ratio on the previous time step
    double       Gs;       ///< Interaction constant between solid and fluid
    size_t       Nneigh;   ///< Number of neighbors
    size_t       ID;       ///< Tag for the particle
    double       Cs;       ///< Velocity of the grid
    double       Rho;      ///< Density of the Cell
    iVec3_t      Index;    ///< Vector of indexes
    Vec3_t       Vel;      ///< Velocity of the Cell
    Vec3_t       VelBC;    ///< Velocity at boundary
    Vec3_t       BForce;   ///< Applied body force
    Vec3_t       BForcef;  ///< Fixed Applied body force
    double       RhoBC;    ///< Density at boundary
    double       Pf;       ///< Percolation parameter for partial porosity calculation

    size_t const  * Op;    ///< Pointer to opposite velocities
    double const  * W;     ///< Pointer to array of weights
    Vec3_t const  * C;     ///< Pointer to velocity constants
    //Array<double>   F;     ///< Distribution functions
    //Array<double>   Ftemp; ///< Temporary distribution functions
    //Array<double>   Omeis; ///< Array of collision operators
    //Array<size_t>   Neighs;///< Array of neighbors indexes
    double        * F;     ///< Distribution functions
    double        * Ftemp; ///< Temporary distribution functions
    double        * Omeis; ///< Array of collision operators
    size_t        * Neighs;///< Array of neighbors indexes
};

inline Cell::Cell(size_t TheID, LBMethod TheMethod, iVec3_t TheIndexes, iVec3_t TheNdim, double TheCs, double TheTau)
{
    IsSolid= false;
    ID     = TheID;
    //Method = TheMethod;
    Index  = TheIndexes;
    Cs     = TheCs;
    Gamma  = 0.0;
    Gammap = 0.0;
    Pf     = 1.0;
    Gs     = 1.0;
    //Tau    = TheTau;
    if (TheMethod==D2Q5)
    {
        Nneigh = 5;
        W      = WEIGHTSD2Q5;
        C      = LVELOCD2Q5;
        Op     = OPPOSITED2Q5;
    }
    if (TheMethod==D2Q9)
    {
        Nneigh = 9;
        W      = WEIGHTSD2Q9;
        C      = LVELOCD2Q9;
        Op     = OPPOSITED2Q9;
    }
    if (TheMethod==D3Q15)
    {
        Nneigh = 15;
        W      = WEIGHTSD3Q15;
        C      = LVELOCD3Q15;
        Op     = OPPOSITED3Q15;
    }
    if (TheMethod==D3Q19)
    {
        Nneigh = 19;
        W      = WEIGHTSD3Q19;
        C      = LVELOCD3Q19;
        Op     = OPPOSITED3Q19;
    }
    //F.     Resize(Nneigh);
    //Ftemp. Resize(Nneigh);
    //Neighs.Resize(Nneigh);
    //Omeis .Resize(Nneigh);
    F      = new double [Nneigh];
    Ftemp  = new double [Nneigh];
    Neighs = new size_t [Nneigh];
    Omeis  = new double [Nneigh];
    BForcef = 0.0,0.0,0.0;

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
#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

inline double Cell::Density()
{
    if (IsSolid) return 0.0;
    double rho = 0.0;
    for (size_t k=0;k<Nneigh;k++) rho += F[k];
    return rho;
}

inline void Cell::Velocity(Vec3_t & V)
{
    V = 0.0, 0.0, 0.0;
    if (IsSolid) return;
    double rho = Density();
    if (rho<1.0e-12) return;
    for (size_t k=0;k<Nneigh;k++) V += F[k]*C[k];
    V *= Cs/rho;
    //V *= Pf*Cs/rho;
}

inline double Cell::VelDen(Vec3_t & V)
{
    V = OrthoSys::O;
    if (IsSolid) return 0.0;
    double rho = 0.0;
    for (size_t k=0;k<Nneigh;k++)
    {
        if (std::isnan(F[k])) F[k] = 1.0e-12;
        V   += F[k]*C[k];
        rho += F[k];
    }
    //if (rho<1.0e-12)
    //{
       //V = OrthoSys::O;
       //return 0.0;
    //}
    if(std::isnan(rho))
    { 
        std::cout << "NaN found in cell: " << Index << std::endl;
        for (size_t k=0;k<Nneigh;k++)
        {
            std::cout << F[k] << std::endl;
        }
        throw new Fatal("NaN found in one of the cells");
    }
    V *= Cs/rho;
    //V *= Pf*Cs/rho;
    //V += 0.5*BForce/rho;
    return rho;
}

inline double Cell::Feq(size_t k, Vec3_t const & V, double TheRho)
{
    double VdotC = dot(V,C[k]);
    double VdotV = dot(V,V);
    return W[k]*TheRho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

inline void Cell::Initialize(double TheRho, Vec3_t const & V)
{
    for (size_t k=0;k<Nneigh;k++) F[k] = Feq(k,V,TheRho);
    Rho = VelDen(Vel);
}

inline void Cell::BounceBack()
{
    for (size_t j = 0;j<Nneigh;j++) Ftemp[j] = F[j];
    for (size_t j = 0;j<Nneigh;j++) F[j]     = Ftemp[Op[j]];
}


const double Cell::WEIGHTSD2Q5   [ 5] = { 2./6., 1./6., 1./6., 1./6., 1./6.};
const double Cell::WEIGHTSD2Q9   [ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
const double Cell::WEIGHTSD3Q15  [15] = { 2./9., 1./9., 1./9., 1./9., 1./9.,  1./9.,  1./9., 1./72., 1./72. , 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
const double Cell::WEIGHTSD3Q19  [19] = { 1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};
const size_t Cell::OPPOSITED2Q5  [ 5] = { 0, 3, 4, 1, 2 };                                                       ///< Opposite directions (D2Q5) 
const size_t Cell::OPPOSITED2Q9  [ 9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };                                           ///< Opposite directions (D2Q9) 
const size_t Cell::OPPOSITED3Q15 [15] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};                     ///< Opposite directions (D3Q15)
const size_t Cell::OPPOSITED3Q19 [19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};     ///< Opposite directions (D3Q19)
const Vec3_t Cell::LVELOCD2Q5  [ 5] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0} };
const Vec3_t Cell::LVELOCD2Q9  [ 9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
const Vec3_t Cell::LVELOCD3Q15 [15] =
{
	{ 0, 0, 0}, { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, 
	{ 0, 0, 1}, { 0, 0,-1}, { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, 
	{-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} 
};
const Vec3_t Cell::LVELOCD3Q19 [19] =
{
	{ 0, 0, 0}, 
    { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1}, 
    { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 1, 0, 1}, {-1, 0,-1},
    { 1, 0,-1}, {-1, 0, 1}, { 0, 1, 1}, { 0,-1,-1}, { 0, 1,-1}, { 0,-1, 1}
};

/*TODO: Go back to the previous version since this is not nice
 *      Changes made for g++ 4.3.2 */
//Vec3_t Cell::LVELOCD2Q9  [ 9]; ///< Local velocities (D2Q9) 
//Vec3_t Cell::LVELOCD3Q15 [15]; ///< Local velocities (D3Q15)
//int __allocate_constants ()
//{
    //Cell::LVELOCD2Q9  [0] =  0 ,  0 , 0;
    //Cell::LVELOCD2Q9  [1] =  1 ,  0 , 0;
    //Cell::LVELOCD2Q9  [2] =  0 ,  1 , 0;
    //Cell::LVELOCD2Q9  [3] = -1 ,  0 , 0;
    //Cell::LVELOCD2Q9  [4] =  0 , -1 , 0;
    //Cell::LVELOCD2Q9  [5] =  1 ,  1 , 0;
    //Cell::LVELOCD2Q9  [6] = -1 ,  1 , 0;
    //Cell::LVELOCD2Q9  [7] = -1 , -1 , 0;
    //Cell::LVELOCD2Q9  [8] =  1 , -1 , 0;
//
    //Cell::LVELOCD3Q15[ 0]= 0, 0, 0;  Cell::LVELOCD3Q15[ 1]= 1, 0, 0;  Cell::LVELOCD3Q15[ 2]=-1, 0, 0;  Cell::LVELOCD3Q15[ 3]= 0, 1, 0;  Cell::LVELOCD3Q15[ 4]= 0,-1, 0;
    //Cell::LVELOCD3Q15[ 5]= 0, 0, 1;  Cell::LVELOCD3Q15[ 6]= 0, 0,-1;  Cell::LVELOCD3Q15[ 7]= 1, 1, 1;  Cell::LVELOCD3Q15[ 8]=-1,-1,-1;  Cell::LVELOCD3Q15[ 9]= 1, 1,-1; 
    //Cell::LVELOCD3Q15[10]=-1,-1, 1;  Cell::LVELOCD3Q15[11]= 1,-1, 1;  Cell::LVELOCD3Q15[12]=-1, 1,-1;  Cell::LVELOCD3Q15[13]= 1,-1,-1;  Cell::LVELOCD3Q15[14]=-1, 1, 1;
//
    //return 0;
//}
//int __dummy__allocate_constants = __allocate_constants();

#endif // MECHSYS_LBM_CELL_H
