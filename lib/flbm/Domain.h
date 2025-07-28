/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2016 Sergio Galindo                                    *
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


#ifndef MECHSYS_FLBM_DOMAIN_H
#define MECHSYS_FLBM_DOMAIN_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#ifdef USE_CUDA
#include <mechsys/flbm/lbm.cuh>
#endif

// Std lib
#ifdef USE_OMP
#include <omp.h>
#endif

//STD
#include<iostream>

// Mechsys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/numstreams.h>

enum LBMethod
{
    D2Q5,     ///< 2D 5 velocities
    D2Q9,     ///< 2D 9 velocities
    D3Q15,    ///< 3D 15 velocities
    D3Q19,    ///< 3D 19 velocities
    //D3Q27     ///< 3D 27 velocities
};

enum LBMSolver
{
    NavierStokes,
    AdvectionDiffusion,
    PhaseFieldIce,
    PhaseField,
    ShallowWater,
};

namespace FLBM
{

inline size_t Pt2idx(iVec3_t iv, iVec3_t & Dim) // Calculates the index of the cell at coordinates iv for a cubic lattice of dimensions Dim
{
    return iv(0) + iv(1)*Dim(0) + iv(2)*Dim(0)*Dim(1);
}

inline void   idx2Pt(size_t n, iVec3_t & iv, iVec3_t & Dim) // Calculates the coordinates from the index
{
    iv(0) = n%Dim(0);
    iv(1) = (n/Dim(0))%(Dim(1));
    iv(2) = n/(Dim(0)*Dim(1));
}

class Domain
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
    static const double   MD2Q5       [5][5]; ///< MRT transformation matrix (D2Q5)
    static const double   MD2Q9       [9][9]; ///< MRT transformation matrix (D2Q9)
    static const double   MD3Q15    [15][15]; ///< MRT transformation matrix (D3Q15)
    static const double   MD3Q19    [19][19]; ///< MRT transformation matrix (D3Q19)
    //static const size_t   MD3Q27    [27][27]; ///< MRT transformation matrix (D3Q27)
    //static const double   SD2Q5          [5]; ///< MRT relaxation time vector (D2Q5)
    //static const double   SD2Q9          [9]; ///< MRT relaxation time vector (D2Q9)
    //static const double   SD3Q15        [15]; ///< MRT relaxation time vector (D3Q15)
    //static const double   SD3Q19        [19]; ///< MRT relaxation time vector (D3Q19)
    //static       double   SD3Q19        [27]; ///< MRT relaxation time vector (D3Q27)
    
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructors
    Domain () {};                 ///< Default constructor
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    Array<double>         nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt,     ///< Time step
    LBMSolver Solver=NavierStokes);    ///< The solver 

    //Special constructor with only one component, the parameters are the same as above
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    double                nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt,     ///< Time step
    LBMSolver Solver=NavierStokes);    ///< The solver 

    //Methods
    void   SolidPhase(Vec3_t const & XC, double RC);                              ///< Same function as SolidDisk producing a solid sphere
    void   SolidCube (Vec3_t const & XC, double RC);                              ///< Add a solid cube with center XC and side LC
    void   ApplyForcesSC();                                                       ///< Apply the molecular forces for the single component case
    void   ApplyForcesSCSS();                                                     ///< Apply the molecular forces for the single component case with solid surfaces
    void   ApplyForcesMP();                                                       ///< Apply the molecular forces for the multiphase case
    void   ApplyForcesSCMP();                                                     ///< Apply the molecular forces for the both previous cases
    void   CollideMRT();                                                          ///< The collide step of LBM with MRT
    void   CollideSC();                                                           ///< The collide step of LBM for single component simulations
    void   CollideSCDEM();                                                        ///< The collide step of LBM for DEM coupling
    void   CollideMPM();                                                          ///< The collide step of LBM for MPM coupling
    void   CollideMP();                                                           ///< The collide step of LBM for multi phase simulations
    void   CollideEFS_SRT();                                                      ///< The collide step of LBM for enhanced forcing scheme
    void   StreamSC();                                                            ///< The stream step of LBM SC
    void   StreamMPM();                                                           ///< The stream step of LBM for MPM coupling
    void   StreamMP();                                                            ///< The stream step of LBM MP
    void   StreamEFS_SRT();                                                       ///< The strean step of LBM EFS
    void   VelDen();                                                              ///< Calculate density and velocity from the functions
    void   Initialize(size_t k, iVec3_t idx, double Rho, Vec3_t & Vel);           ///< Initialize each cell with a given density and velocity
    void   Initialize(iVec3_t idx, double Pre, double H, double Phi
            , Vec3_t & Vel);                                                      ///< Same as previous, but for the Phase Field Ice model
    void   Initialize(iVec3_t idx, double Pre, double Phi
            , Vec3_t & Vel);                                                      ///< Same as previous, but for the Phase Field model
    void   InitializeSW(iVec3_t idx, double H, Vec3_t & Vel);                     ///< Initialize each cell with a given height and velocity for SW
                                                                                  ///solver
    double Feq(size_t k, double Rho, Vec3_t & Vel);                               ///< The equilibrium function
    double FeqSW(size_t k, double H, Vec3_t & Vel);                               ///< The equilibrium function for SW solver
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);            ///< Solve the Domain dynamics
    
    //Writing Methods
    void WriteXDMF         (char const * FileKey);                                ///< Write the domain data in xdmf file
    void WriteXDMF_DEM     (char const * FileKey);                                ///< Write the domain data in xdmf file for DEM simulations
   
    //Methods for CUDA 
    #ifdef USE_CUDA
    void UpLoadDevice (size_t Nc = 1);                                            ///< Upload the buffers into the coprocessor device
    void DnLoadDevice (size_t Nc = 1);                                            ///< Download the buffers from the coprocessor device
    #endif
    
    #ifdef USE_OMP
    omp_lock_t      lck;                      ///< to protect variables in multithreading
    #endif

    //Data
    LBMSolver     Solver;                      ///< What type of solver is used in LBM
    double *****  F;                           ///< The array containing the individual functions with the order of the lattice, the x,y,z coordinates and the order of the function.
    double *****  Ftemp;                       ///< A similar array to hold provitional data
    double ****   Omeis;                       ///< An array for collision detection
    bool   ****   IsSolid;                     ///< An array of bools with an identifier to see if the cell is a solid cell
    bool   ***    Inside;                      ///< An array of bools with an identifier to check if the cell is inside a solid
    int    ***    DEMInside;                   ///< An array of indexes of DEM particles signalign if the cell is inside a particle or not
    Vec3_t ****   Vel;                         ///< The fluid velocities
    Vec3_t ****   BForce;                      ///< Body Force for each cell
    double ****   Rho;                         ///< The fluid densities
    double ***    Gamma;                       ///< Information on the overlapping volume fraction
    double ***    Gammaf;                      ///< Information on the prescribed overlapping volume fraction
    double *      Tau;                         ///< The characteristic time of the lattice
    //Shallow Water parameters
    double        g;                           ///< Gravity for SW
    //Shan Chen model parameters
    double *      G;                           ///< The attractive constant for multiphase simulations
    double *      Gs;                          ///< The attractive constant for solid phase
    double *      Psi;                         ///< Parameters for the Shan Chen pseudo potential
    double *      Rhoref;                      ///< Parameters for the Shan Chen pseudo potential
    double        Gmix;                        ///< The mixing constant for multicomponent simulations
    //Phase Field for Ice parameters 0 solid 1 liquid 2 gas                                             
    double        rho[3];                      ///< Density of the phases
    double        cap[3];                      ///< Heat capcity for each phase
    double        kap[3];                      ///< heat conductivity for each phase    
    double        thick;                       ///< thickness of the phase field interfase;
    double        sigma;                       ///< surface tension of the interfase
    double        Ts;                          ///< Solidus temperature
    double        Tl;                          ///< Liquidus temperature
    double        L;                           ///< Latent heat
    
    size_t const  * Op;                        ///< An array containing the indexes of the opposite direction for bounce back conditions
    double const  *  W;                        ///< An array with the direction weights
    double *      EEk;                         ///< Diadic product of the velocity vectors
    Vec3_t const  *  C;                        ///< The array of lattice velocities
    Mat_t         M;                           ///< Transformation matrix to momentum space for MRT calculations
    Mat_t         Minv;                        ///< Inverse Transformation matrix to momentum space for MRT calculations
    Vec_t         S;                           ///< Vector of relaxation times for MRT
    size_t        Nneigh;                      ///< Number of Neighbors, depends on the scheme
    double        dt;                          ///< Time Step
    double        dx;                          ///< Grid size
    double        Cs;                          ///< Lattice Velocity
    bool          IsFirstTime;                 ///< Bool variable checking if it is the first time function Setup is called
    iVec3_t       Ndim;                        ///< Lattice Dimensions
    size_t        Ncells;                      ///< Number of cells
    size_t        Nproc;                       ///< Number of processors for openmp
    size_t        Nthread;                     ///< Number of CUDA threads
    size_t        idx_out;                     ///< The discrete time step for output
    String        FileKey;                     ///< File Key for output files
    void *        UserData;                    ///< User Data
    size_t        Step;                        ///< Lenght of averaging cube to save data
    double        Time;                        ///< Simulation time variable
    size_t        Nl;                          ///< Number of lattices (fluids)
    double        Sc;                          ///< Smagorinsky constant
    Array<String> Rhonames;                    ///< Array of strings for xdmf files                                          

    //Array for pair calculation
    size_t       NCellPairs;                  ///< Number of cell pairs
    iVec3_t   *  CellPairs;                   ///< Pairs of cells for molecular force calculation
                                              
    //Speed up ignoring deep solid boundaries data set
    #ifdef IGNORESOLID
    Array<size_t> ValidNodes;
    #endif

    #ifdef USE_CUDA
    thrust::device_vector<real>      bF;                      ///< Buffer with the distribution functions
    thrust::device_vector<real>      bFtemp;                  ///< Buffer with the distribution functions temporal
    thrust::device_vector<bool>      bIsSolid;                ///< Buffer with the solid bool information
    thrust::device_vector<real3>     bBForce;                 ///< Buffer with the body forces
    thrust::device_vector<real3>     bVel;                    ///< Buffer with the cell velocities
    thrust::device_vector<real>      bRho;                    ///< Buffer with the cell densities
    thrust::device_vector<uint3>     bCellPairs;              ///< Pairs of cells for molecular force calculation
    //thrust::device_vector<lbm_aux>   blbmaux;                 ///< Buffer with the structure containing generic lbm information
    real                           * pF;                      ///< Pointer to Buffer with the distribution functions
    real                           * pFtemp;                  ///< Pointer to Buffer with the distribution functions temporal
    bool                           * pIsSolid;                ///< Pointer to Buffer with the solid bool information
    real3                          * pBForce;                 ///< Pointer to Buffer with the body forces
    real3                          * pVel;                    ///< Pointer to Buffer with the cell velocities
    real                           * pRho;                    ///< Pointer to Buffer with the cell densities
    uint3                          * pCellPairs;              ///< 
    lbm_aux                        * plbmaux;                 ///< Pointer to Buffer with the structure containing generic lbm information
    #endif
};

inline Domain::Domain(LBMethod TheMethod, Array<double> nu, iVec3_t TheNdim, double Thedx, double Thedt, LBMSolver TheSolver)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (nu.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one fluid please");
    if (TheNdim(2) >1&&(TheMethod==D2Q9 ||TheMethod==D2Q5 ))  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (TheNdim(2)==1&&(TheMethod==D3Q15||TheMethod==D3Q19))  throw new Fatal("LBM::Domain: Ndim(2) is 1. Either change the method to D2Q9 or increase the z-dimension");
   
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
        M.Resize(Nneigh,Nneigh);
        Minv.Resize(Nneigh,Nneigh);
        for (size_t n=0;n<Nneigh;n++)
        for (size_t m=0;m<Nneigh;m++)
        {
            M(n,m) = MD3Q15[n][m];
        }
        Inv(M,Minv);
    }
    if (TheMethod==D3Q19)
    {
        Nneigh = 19;
        W      = WEIGHTSD3Q19;
        C      = LVELOCD3Q19;
        Op     = OPPOSITED3Q19;
    }
    

    Time        = 0.0;
    dt          = Thedt;
    dx          = Thedx;
    Cs          = dx/dt;
    Step        = 1;
    Sc          = 0.17;
    Nl          = nu.Size();
    Ndim        = TheNdim;
    Ncells      = Ndim(0)*Ndim(1)*Ndim(2);
    NCellPairs  = 0;
    Nthread     = 256;
    IsFirstTime = true;
    Gamma       = NULL;
    Gammaf      = NULL;
    Omeis       = NULL;
    Rhonames.Resize(Nl);
    Solver = TheSolver;

    if (Solver==AdvectionDiffusion)
    {
        Rhonames[0].Printf("Density");
        Rhonames[1].Printf("Concentration");
    }

    if (Solver==PhaseField)
    {
        Rhonames[0].Printf("Pressure");
        Rhonames[1].Printf("Phase");
        rho[0] = 1.0;
        rho[1] = 1.0;
        thick  = 4.0;
        sigma  = 1.0;
    }

    if (Solver==PhaseFieldIce)
    {
        Rhonames[0].Printf("Pressure");
        Rhonames[1].Printf("Phase");
        Rhonames[2].Printf("Enthalpy");

        rho[0] = 1.0;
        rho[1] = 1.0;
        rho[2] = 1.0;
        cap[0] = 1.0;
        cap[1] = 1.0;
        cap[2] = 1.0;
        kap[0] = 1.0;
        kap[1] = 1.0;
        kap[2] = 1.0;
        thick  = 4.0;
        sigma  = 1.0;
        Ts     = 0.8;
        Tl     = 1.0;
        L      = 1.0;
    }


    Tau    = new double [Nl];

    if (Solver==NavierStokes)
    {
        G      = new double [Nl];
        Gs     = new double [Nl];
        Rhoref = new double [Nl];
        Psi    = new double [Nl];
        Gmix= 0.0;
    }

    F       = new double **** [Nl];
    Ftemp   = new double **** [Nl];
    Vel     = new Vec3_t ***  [Nl];
    BForce  = new Vec3_t ***  [Nl];
    Rho     = new double ***  [Nl];
    IsSolid = new bool   ***  [Nl];
    for (size_t i=0;i<Nl;i++)
    {
        Tau     [i]    = 3.0*nu[i]*dt/(dx*dx)+0.5;
        if (Solver==NavierStokes) 
        {
            G       [i]    = 0.0;
            Gs      [i]    = 0.0;
            Rhoref  [i]    = 200.0;
            Psi     [i]    = 4.0;
            Rhonames[i].Printf("Density_%d",i);
        }
        F       [i]    = new double *** [Ndim(0)];
        Ftemp   [i]    = new double *** [Ndim(0)];
        Vel     [i]    = new Vec3_t **  [Ndim(0)];
        BForce  [i]    = new Vec3_t **  [Ndim(0)];
        Rho     [i]    = new double **  [Ndim(0)];
        IsSolid [i]    = new bool   **  [Ndim(0)];

        for (size_t nx=0;nx<Ndim(0);nx++)
        {
            F       [i][nx]    = new double ** [Ndim(1)];
            Ftemp   [i][nx]    = new double ** [Ndim(1)];
            Vel     [i][nx]    = new Vec3_t *  [Ndim(1)];
            BForce  [i][nx]    = new Vec3_t *  [Ndim(1)];
            Rho     [i][nx]    = new double *  [Ndim(1)];
            IsSolid [i][nx]    = new bool   *  [Ndim(1)];
            for (size_t ny=0;ny<Ndim(1);ny++)
            {
                F       [i][nx][ny]    = new double * [Ndim(2)];
                Ftemp   [i][nx][ny]    = new double * [Ndim(2)];
                Vel     [i][nx][ny]    = new Vec3_t   [Ndim(2)];
                BForce  [i][nx][ny]    = new Vec3_t   [Ndim(2)];
                Rho     [i][nx][ny]    = new double   [Ndim(2)];
                IsSolid [i][nx][ny]    = new bool     [Ndim(2)];
                for (size_t nz=0;nz<Ndim(2);nz++)
                {
                    F    [i][nx][ny][nz]    = new double [Nneigh];
                    Ftemp[i][nx][ny][nz]    = new double [Nneigh];
                    IsSolid[i][nx][ny][nz]  = false;
                    for (size_t nn=0;nn<Nneigh;nn++)
                    {
                        F    [i][nx][ny][nz][nn] = 0.0;
                        Ftemp[i][nx][ny][nz][nn] = 0.0;
                    }
                }
            }
        }
    }


    EEk = new double [Nneigh];
    for (size_t k=0;k<Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(C[k][n]*C[k][m]);
        }
    }

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Nl*Ncells,TERM_RST);
    printf("%s  Num of cells x = %zd%s\n",TERM_CLR2,Ndim(0)  ,TERM_RST);
    printf("%s  Num of cells y = %zd%s\n",TERM_CLR2,Ndim(1)  ,TERM_RST);
    printf("%s  Num of cells z = %zd%s\n",TERM_CLR2,Ndim(2)  ,TERM_RST);
#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

inline Domain::Domain(LBMethod TheMethod, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt, LBMSolver TheSolver)
{
    Array<double> nu(1);
    nu[0] = Thenu;
    
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (nu.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one fluid please");
    if (TheNdim(2) >1&&(TheMethod==D2Q9 ||TheMethod==D2Q5 ))  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (TheNdim(2)==1&&(TheMethod==D3Q15||TheMethod==D3Q19))  throw new Fatal("LBM::Domain: Ndim(2) is 1. Either change the method to D2Q9 or increase the z-dimension");
   
    if (Nl==1&&Solver==AdvectionDiffusion) throw new Fatal ("FLBM::Solve: AdvectionDiffusion solver needs two LBM layers");


    Time        = 0.0;
    dt          = Thedt;
    dx          = Thedx;
    Cs          = dx/dt;
    Step        = 1;
    Sc          = 0.17;
    Nl          = 1;
    Ndim        = TheNdim;
    Ncells      = Ndim(0)*Ndim(1)*Ndim(2);
    NCellPairs  = 0;
    Nthread     = 256;
    IsFirstTime = true;
    Gamma       = NULL;
    Gammaf      = NULL;
    Omeis       = NULL;
    Rhonames.Resize(Nl);
    Solver = TheSolver;

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
        M.Resize(Nneigh,Nneigh);
        Minv.Resize(Nneigh,Nneigh);
        for (size_t n=0;n<Nneigh;n++)
        for (size_t m=0;m<Nneigh;m++)
        {
            M(n,m) = MD3Q15[n][m];
        }
        Inv(M,Minv);
        double tau = 3.0*Thenu*dt/(dx*dx)+0.5;
        double s   = 8.0*(2.0-1.0/tau)/(8.0-1.0/tau);
        S.Resize(Nneigh);
        S = 0.0,1.0/tau,1.0/tau,0.0,s,0.0,s,0.0,s,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,s;
    }
    if (TheMethod==D3Q19)
    {
        Nneigh = 19;
        W      = WEIGHTSD3Q19;
        C      = LVELOCD3Q19;
        Op     = OPPOSITED3Q19;
    }
    
    if (Solver==ShallowWater)
    {
        Rhonames[0].Printf("Height");
        rho[0] = 1.0;
        cap[0] = 1.0;
        g      = 1.0;
        if (TheMethod!=D2Q9) throw new Fatal("ShallowWater solver only works with velocity configuration D2Q9");
    }


    Tau    = new double [Nl];
    G      = new double [Nl];
    Gs     = new double [Nl];
    Rhoref = new double [Nl];
    Psi    = new double [Nl];
    Gmix   = 0.0;

    F       = new double **** [Nl];
    Ftemp   = new double **** [Nl];
    Vel     = new Vec3_t ***  [Nl];
    BForce  = new Vec3_t ***  [Nl];
    Rho     = new double ***  [Nl];
    IsSolid = new bool   ***  [Nl];

    for (size_t i=0;i<Nl;i++)
    {
        Tau     [i]    = 3.0*nu[i]*dt/(dx*dx)+0.5;
        F       [i]    = new double *** [Ndim(0)];
        Ftemp   [i]    = new double *** [Ndim(0)];
        Vel     [i]    = new Vec3_t **  [Ndim(0)];
        BForce  [i]    = new Vec3_t **  [Ndim(0)];
        Rho     [i]    = new double **  [Ndim(0)];
        IsSolid [i]    = new bool   **  [Ndim(0)];
        if (Solver==NavierStokes)
        {
            G       [i]    = 0.0;
            Gs      [i]    = 0.0;
            Rhoref  [i]    = 200.0;
            Psi     [i]    = 4.0;
            Rhonames[i].Printf("Density_%d",i);
        }
        for (size_t nx=0;nx<Ndim(0);nx++)
        {
            F       [i][nx]    = new double ** [Ndim(1)];
            Ftemp   [i][nx]    = new double ** [Ndim(1)];
            Vel     [i][nx]    = new Vec3_t *  [Ndim(1)];
            BForce  [i][nx]    = new Vec3_t *  [Ndim(1)];
            Rho     [i][nx]    = new double *  [Ndim(1)];
            IsSolid [i][nx]    = new bool   *  [Ndim(1)];
            for (size_t ny=0;ny<Ndim(1);ny++)
            {
                F       [i][nx][ny]    = new double * [Ndim(2)];
                Ftemp   [i][nx][ny]    = new double * [Ndim(2)];
                Vel     [i][nx][ny]    = new Vec3_t   [Ndim(2)];
                BForce  [i][nx][ny]    = new Vec3_t   [Ndim(2)];
                Rho     [i][nx][ny]    = new double   [Ndim(2)];
                IsSolid [i][nx][ny]    = new bool     [Ndim(2)];
                for (size_t nz=0;nz<Ndim(2);nz++)
                {
                    F    [i][nx][ny][nz]    = new double [Nneigh];
                    Ftemp[i][nx][ny][nz]    = new double [Nneigh];
                    IsSolid[i][nx][ny][nz]  = false;
                    for (size_t nn=0;nn<Nneigh;nn++)
                    {
                        F    [i][nx][ny][nz][nn] = 0.0;
                        Ftemp[i][nx][ny][nz][nn] = 0.0;
                    }
                }
            }
        }
    }


    EEk = new double [Nneigh];
    for (size_t k=0;k<Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(C[k][n]*C[k][m]);
        }
    }

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Nl*Ncells,TERM_RST);
    printf("%s  Num of cells x = %zd%s\n",TERM_CLR2,Ndim(0)  ,TERM_RST);
    printf("%s  Num of cells y = %zd%s\n",TERM_CLR2,Ndim(1)  ,TERM_RST);
    printf("%s  Num of cells z = %zd%s\n",TERM_CLR2,Ndim(2)  ,TERM_RST);
#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
}

inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Ndim[0]/Step;
    size_t  Ny = Ndim[1]/Step;
    size_t  Nz = Ndim[2]/Step;

    for (size_t j=0;j<Nl;j++)
    {
        // Creating data sets
        double * Density   = new double[  Nx*Ny*Nz];
        double * Gammap    = new double[  Nx*Ny*Nz];
        double * Vvec      = new double[3*Nx*Ny*Nz];

        size_t i=0;
        for (size_t m=0;m<Ndim(2);m+=Step)
        for (size_t l=0;l<Ndim(1);l+=Step)
        for (size_t n=0;n<Ndim(0);n+=Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            Vec3_t vel    = OrthoSys::O;

            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho    += Rho    [j][n+ni][l+li][m+mi];
                gamma  += IsSolid[j][n+ni][l+li][m+mi] ? 1.0: 0.0;
                vel    += Vel    [j][n+ni][l+li][m+mi];
            }
            rho  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            vel  /= Step*Step*Step;
            Density [i]  = (double) rho;
            Gammap  [i]  = (double) gamma;
            Vvec[3*i  ]  = (double) vel(0);
            Vvec[3*i+1]  = (double) vel(1);
            Vvec[3*i+2]  = (double) vel(2);
            i++;
        }
        
        //Writing data to h5 file
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname = Rhonames[j];
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Density );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Gammap  );
        }
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vvec    );
        dims[0] = 1;
        int N[1];
        N[0] = Nx;
        dsname.Printf("Nx");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Ny;
        dsname.Printf("Ny");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Nz;
        dsname.Printf("Nz");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);

        delete [] Density ;
        delete [] Gammap  ;
        delete [] Vvec    ;
    }

    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    //std::cout << "2" << std::endl;

    if (Ndim(2)==1)
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Ndim(1) << " " << Ndim(0) << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"" << Rhonames[j] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/" << Rhonames[j] << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    
    else
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"LBM_Mesh\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*dx << " " << Step*dx  << " " << Step*dx  << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"" << Rhonames[j] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/" << Rhonames[j] << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

inline void Domain::WriteXDMF_DEM(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Ndim[0]/Step;
    size_t  Ny = Ndim[1]/Step;
    size_t  Nz = Ndim[2]/Step;

    for (size_t j=0;j<Nl;j++)
    {
        // Creating data sets
        double * Density   = new double[  Nx*Ny*Nz];
        double * Gammap    = new double[  Nx*Ny*Nz];
        double * Vvec      = new double[3*Nx*Ny*Nz];

        size_t i=0;
        for (size_t m=0;m<Ndim(2);m+=Step)
        for (size_t l=0;l<Ndim(1);l+=Step)
        for (size_t n=0;n<Ndim(0);n+=Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            Vec3_t vel    = OrthoSys::O;

            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho    += (1.0-Gamma[n+ni][l+li][m+mi])*Rho    [j][n+ni][l+li][m+mi];
                //rho    += Rho    [j][n+ni][l+li][m+mi];
                gamma  += Gamma[n+ni][l+li][m+mi];
                vel    += (1.0-Gamma[n+ni][l+li][m+mi])*Vel    [j][n+ni][l+li][m+mi];
                //vel    += Vel    [j][n+ni][l+li][m+mi];
            }
            rho  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            vel  /= Step*Step*Step;
            Density [i]  = (double) rho;
            Gammap  [i]  = (double) gamma;
            Vvec[3*i  ]  = (double) vel(0);
            Vvec[3*i+1]  = (double) vel(1);
            Vvec[3*i+2]  = (double) vel(2);
            i++;
        }
        
        //Writing data to h5 file
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Density );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Gammap  );
        }
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vvec    );
        dims[0] = 1;
        int N[1];
        N[0] = Nx;
        dsname.Printf("Nx");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Ny;
        dsname.Printf("Ny");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Nz;
        dsname.Printf("Nz");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);

        delete [] Density ;
        delete [] Gammap  ;
        delete [] Vvec    ;
    }

    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    //std::cout << "2" << std::endl;

    if (Ndim(2)==1)
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Ndim(1) << " " << Ndim(0) << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    
    else
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"LBM_Mesh\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*dx << " " << Step*dx  << " " << Step*dx  << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

inline double Domain::Feq(size_t k, double Rho, Vec3_t & V)
{
    double VdotC = dot(V,C[k]);
    double VdotV = dot(V,V);
    return W[k]*Rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

inline double Domain::FeqSW(size_t k, double h, Vec3_t & V)
{
    double VdotC = dot(V,C[k]);
    double VdotV = dot(V,V);
    if (k==0)
    {
        return h - 5.0/6.0*g*h*h/(Cs*Cs) - 2.0/3.0*h*VdotV;
    }
    else
    {
        return W[k]*h*(1.5*g*h/(Cs*Cs) + 3.0*VdotC + 4.5*VdotC*VdotC - 1.5*VdotV);
    }
}

inline void Domain::Initialize(size_t il, iVec3_t idx, double TheRho, Vec3_t & TheVel)
{
    size_t ix = idx(0);
    size_t iy = idx(1);
    size_t iz = idx(2);

    BForce[il][ix][iy][iz] = OrthoSys::O;

    for (size_t k=0;k<Nneigh;k++)
    {
        F[il][ix][iy][iz][k] = Feq(k,TheRho,TheVel);
    }

    if (!IsSolid[il][ix][iy][iz])
    {
        Vel[il][ix][iy][iz] = TheVel;
        Rho[il][ix][iy][iz] = TheRho;
    }
    else
    {
        Vel[il][ix][iy][iz] = OrthoSys::O;
        Rho[il][ix][iy][iz] = 0.0;
    }
}

inline void Domain::InitializeSW(iVec3_t idx, double TheH, Vec3_t & TheVel)
{
    size_t ix = idx(0);
    size_t iy = idx(1);
    size_t iz = idx(2);

    BForce[0][ix][iy][iz] = OrthoSys::O;

    for (size_t k=0;k<Nneigh;k++)
    {
        F[0][ix][iy][iz][k] = FeqSW(k,TheH,TheVel);
    }

    if (!IsSolid[0][ix][iy][iz])
    {
        Vel[0][ix][iy][iz] = TheVel;
        Rho[0][ix][iy][iz] = TheH;
    }
    else
    {
        Vel[0][ix][iy][iz] = OrthoSys::O;
        Rho[0][ix][iy][iz] = 0.0;
    }
}

inline void Domain::Initialize(iVec3_t idx, double ThePre, double TheH, double ThePhi, Vec3_t & TheVel)
{
    size_t ix = idx(0);
    size_t iy = idx(1);
    size_t iz = idx(2);
    for (size_t il=0;il<3;il++)
    {
        BForce[il][ix][iy][iz] = OrthoSys::O;
    }

    double Hs = Ts*cap[0];
    double Hl = Tl*cap[1];
    double fl = 0.0;
    if      (TheH>=Hs&&TheH<=Hl) fl = (TheH-Hs)/(Hl-Hs);
    else if (TheH>Hl)            fl = 1.0;
    Vel[2][ix][iy][0](0) = fl;
    Vel[2][ix][iy][0](1) = fl;
    Vel[2][ix][iy][0](2) = 0.0;
    double fs   = 1.0-fl;

    double Cp   = ThePhi*(fs*cap[0] + fl*cap[1])+fl*(1.0-ThePhi)*cap[2];
    double Temp = TheH/Cp;
    if      (TheH>=Hs&&TheH<=Hl) Temp = Ts + (TheH-Hs)/(Hl-Hs)*(Tl - Ts);
    else if (TheH>Hl)            Temp = Tl + (TheH-Hl)/Cp;

    double rhof = ThePhi*(fs*rho[0] + fl*rho[1])+fl*(1.0-ThePhi)*rho[2];
    Rho[0][ix][iy][iz] = 0.0;
    double VdotV = dot(TheVel,TheVel);
    for (size_t k=0;k<Nneigh;k++)
    {
        // Pressure
        double VdotC = dot(TheVel,C[k]);
        double sk    = 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs);
        F[0][ix][iy][iz][k] = W[k]*(3.0*ThePre/(Cs*Cs) + rhof*sk);
        if (k==0) F[0][ix][iy][iz][k] += -3.0*ThePre/(Cs*Cs);
        if (k!=0) Rho[0][ix][iy][iz]  += F[0][ix][iy][iz][k];

        //Phase Field
        F[1][ix][iy][iz][k] = W[k]*ThePhi*(1.0 + 3.0*VdotC/Cs);
        Rho[1][ix][iy][iz] += F[1][ix][iy][iz][k];

        //Enthalpy
        F[2][ix][iy][iz][k] = W[k]*Cp*Temp*(1.0 + sk);
        if (k==0) F[2][ix][iy][iz][k] += TheH - Cp*Temp;
        Rho[2][ix][iy][iz] += F[2][ix][iy][iz][k];
    }

    Rho[0][ix][iy][iz] *= Cs*Cs/3.0/(1.0-W[0]);


    Vel[0][ix][iy][iz] = TheVel;
    Vel[1][ix][iy][iz] = ThePhi*TheVel;



}

inline void Domain::Initialize(iVec3_t idx, double ThePre, double ThePhi, Vec3_t & TheVel)
{
    size_t ix = idx(0);
    size_t iy = idx(1);
    size_t iz = idx(2);
    for (size_t il=0;il<2;il++)
    {
        BForce[il][ix][iy][iz] = OrthoSys::O;
    }

    double rhof = ThePhi*rho[0]+(1.0-ThePhi)*rho[1];
    Rho[0][ix][iy][iz] = 0.0;
    double VdotV = dot(TheVel,TheVel);
    for (size_t k=0;k<Nneigh;k++)
    {
        // Pressure
        double VdotC = dot(TheVel,C[k]);
        double sk    = 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs);
        F[0][ix][iy][iz][k] = W[k]*(3.0*ThePre/(Cs*Cs) + rhof*sk);
        if (k==0) F[0][ix][iy][iz][k] += -3.0*ThePre/(Cs*Cs);
        if (k!=0) Rho[0][ix][iy][iz]  += F[0][ix][iy][iz][k];

        //Phase Field
        F[1][ix][iy][iz][k] = W[k]*ThePhi*(1.0 + 3.0*VdotC/Cs);
        Rho[1][ix][iy][iz] += F[1][ix][iy][iz][k];

    }

    Rho[0][ix][iy][iz] *= Cs*Cs/3.0/(1.0-W[0]);


    Vel[0][ix][iy][iz] = TheVel;
    Vel[1][ix][iy][iz] = ThePhi*TheVel;


}

inline void Domain::SolidPhase(Vec3_t const & XC, double RC)
{
    for (size_t m=0; m<Ndim(0); m++)
    for (size_t n=0; n<Ndim(1); n++)
    for (size_t l=0; l<Ndim(2); l++)  
    {
        Vec3_t XX(m*dx,n*dx,l*dx);
        for (size_t il=0;il<Nl;il++)
        {
            if (norm(XX-XC) < RC)    IsSolid [il][m][n][l] = true ;
        }
    }
}

inline void Domain::SolidCube(Vec3_t const & XC, double LC)
{
    for (size_t m=0; m<Ndim(0); m++)
    for (size_t n=0; n<Ndim(1); n++)
    for (size_t l=0; l<Ndim(2); l++)  
    {
        Vec3_t XX(m*dx,n*dx,l*dx);
        for (size_t il=0;il<Nl;il++)
        {
            if ((fabs(XX(0)-XC(0)) < 0.5*LC)&&(fabs(XX(1)-XC(1)) < 0.5*LC)&&(fabs(XX(2)-XC(2)) < 0.5*LC))    IsSolid [il][m][n][l] = true ;
        }
    }
}

inline void Domain::ApplyForcesSC()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<NCellPairs;n++)
    {
        iVec3_t idxc,idxn;
        size_t k = CellPairs[n](2);
        idx2Pt(CellPairs[n](0),idxc,Ndim);
        idx2Pt(CellPairs[n](1),idxn,Ndim);

        size_t ixc = idxc(0);
        size_t iyc = idxc(1);
        size_t izc = idxc(2);
        size_t ixn = idxn(0);
        size_t iyn = idxn(1);
        size_t izn = idxn(2);

        double psic = 0.0;
        double psin = 0.0;

        IsSolid[0][ixc][iyc][izc] ? psic = 0.0 : psic = Psi[0]*exp(-Rhoref[0]/Rho[0][ixc][iyc][izc]);
        IsSolid[0][ixn][iyn][izn] ? psin = 0.0 : psin = Psi[0]*exp(-Rhoref[0]/Rho[0][ixn][iyn][izn]);
//
        Vec3_t bforce = -G[0]*W[k]*C[k]*psic*psin;

        BForce[0][ixc][iyc][izc] += bforce;
        BForce[0][ixn][iyn][izn] -= bforce;
        //
        //IsSolid[0][ixc][iyc][izc] ? psic = 0.0 : psic = 0.5*Psi[0]*Psi[0]*exp(-2.0*Rhoref[0]/Rho[0][ixc][iyc][izc]);
        //IsSolid[0][ixn][iyn][izn] ? psin = 0.0 : psin = 0.5*Psi[0]*Psi[0]*exp(-2.0*Rhoref[0]/Rho[0][ixn][iyn][izn]);
//
        //Vec3_t bforcec = -G[0]*W[k]*C[k]*psin;
        //Vec3_t bforcen = -G[0]*W[k]*C[k]*psic;
//
        //BForce[0][ixc][iyc][izc] += bforcec;
        //BForce[0][ixn][iyn][izn] -= bforcen;
    }
}

inline void Domain::ApplyForcesSCSS()
{
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t n=0;n<NCellPairs;n++)
    {
        iVec3_t idxc,idxn;
        size_t k = CellPairs[n](2);
        idx2Pt(CellPairs[n](0),idxc,Ndim);
        idx2Pt(CellPairs[n](1),idxn,Ndim);

        size_t ixc = idxc(0);
        size_t iyc = idxc(1);
        size_t izc = idxc(2);
        size_t ixn = idxn(0);
        size_t iyn = idxn(1);
        size_t izn = idxn(2);

        double psic  = 0.0;
        double psin  = 0.0;
        bool   solid = false;

        IsSolid[0][ixc][iyc][izc] ? psic = 1.0, solid = true : psic = Rho[0][ixc][iyc][izc];
        IsSolid[0][ixn][iyn][izn] ? psin = 1.0, solid = true : psin = Rho[0][ixn][iyn][izn];

        Vec3_t bforce = OrthoSys::O;
        if (solid) bforce = -Gs[0]*W[k]*C[k]*psic*psin;

        BForce[0][ixc][iyc][izc] += bforce;
        BForce[0][ixn][iyn][izn] -= bforce;
    }
}

inline void Domain::ApplyForcesMP()
{
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t n=0;n<NCellPairs;n++)
    {
        iVec3_t idxc,idxn;
        size_t k = CellPairs[n](2);
        idx2Pt(CellPairs[n](0),idxc,Ndim);
        idx2Pt(CellPairs[n](1),idxn,Ndim);

        size_t ixc = idxc(0);
        size_t iyc = idxc(1);
        size_t izc = idxc(2);
        size_t ixn = idxn(0);
        size_t iyn = idxn(1);
        size_t izn = idxn(2);

        double psic = 0.0;
        double psin = 0.0;

        double Gt   = Gmix;

        IsSolid[0][ixc][iyc][izc] ? psic = 1.0, Gt = Gs[1] : psic = Rho[0][ixc][iyc][izc];
        IsSolid[1][ixn][iyn][izn] ? psin = 1.0, Gt = Gs[0] : psin = Rho[1][ixn][iyn][izn];

        Vec3_t bforce = -Gt*W[k]*C[k]*psic*psin;

        BForce[0][ixc][iyc][izc] += bforce;
        BForce[1][ixn][iyn][izn] -= bforce;

        Gt          = Gmix;

        IsSolid[1][ixc][iyc][izc] ? psic = 1.0, Gt = Gs[0] : psic = Rho[1][ixc][iyc][izc];
        IsSolid[0][ixn][iyn][izn] ? psin = 1.0, Gt = Gs[1] : psin = Rho[0][ixn][iyn][izn];

        bforce      = -Gt*W[k]*C[k]*psic*psin;

        BForce[1][ixc][iyc][izc] += bforce;
        BForce[0][ixn][iyn][izn] -= bforce;
    }
}

inline void Domain::ApplyForcesSCMP()
{
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t n=0;n<NCellPairs;n++)
    {
        iVec3_t idxc,idxn;
        size_t k = CellPairs[n](2);
        idx2Pt(CellPairs[n](0),idxc,Ndim);
        idx2Pt(CellPairs[n](1),idxn,Ndim);

        size_t ixc = idxc(0);
        size_t iyc = idxc(1);
        size_t izc = idxc(2);
        size_t ixn = idxn(0);
        size_t iyn = idxn(1);
        size_t izn = idxn(2);

        double psic = 0.0;
        double psin = 0.0;

        IsSolid[0][ixc][iyc][izc] ? psic = 0.0 : psic = Psi[0]*exp(-Rhoref[0]/Rho[0][ixc][iyc][izc]);
        IsSolid[0][ixn][iyn][izn] ? psin = 0.0 : psin = Psi[0]*exp(-Rhoref[0]/Rho[0][ixn][iyn][izn]);

        Vec3_t bforce = -G[0]*W[k]*C[k]*psic*psin;

        BForce[0][ixc][iyc][izc] += bforce;
        BForce[0][ixn][iyn][izn] -= bforce;

        IsSolid[1][ixc][iyc][izc] ? psic = 0.0 : psic = Psi[1]*exp(-Rhoref[1]/Rho[1][ixc][iyc][izc]);
        IsSolid[1][ixn][iyn][izn] ? psin = 0.0 : psin = Psi[1]*exp(-Rhoref[1]/Rho[1][ixn][iyn][izn]);

        bforce        = -G[1]*W[k]*C[k]*psic*psin;

        BForce[1][ixc][iyc][izc] += bforce;
        BForce[1][ixn][iyn][izn] -= bforce;

        double Gt   = Gmix;

        IsSolid[0][ixc][iyc][izc] ? psic = 1.0, Gt = Gs[1] : psic = Rho[0][ixc][iyc][izc];
        IsSolid[1][ixn][iyn][izn] ? psin = 1.0, Gt = Gs[0] : psin = Rho[1][ixn][iyn][izn];

        bforce      = -Gt*W[k]*C[k]*psic*psin;

        BForce[0][ixc][iyc][izc] += bforce;
        BForce[1][ixn][iyn][izn] -= bforce;

        Gt          = Gmix;

        IsSolid[1][ixc][iyc][izc] ? psic = 1.0, Gt = Gs[0] : psic = Rho[1][ixc][iyc][izc];
        IsSolid[0][ixn][iyn][izn] ? psin = 1.0, Gt = Gs[1] : psin = Rho[0][ixn][iyn][izn];

        bforce      = -Gt*W[k]*C[k]*psic*psin;

        BForce[1][ixc][iyc][izc] += bforce;
        BForce[0][ixn][iyn][izn] -= bforce;
    }
}

inline void Domain::CollideSC()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #ifdef IGNORESOLID
    for (size_t nv=0;nv<ValidNodes.Size();nv++)
    #else
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    #endif
    {
        #ifdef IGNORESOLID
        iVec3_t idx;
        idx2Pt(ValidNodes[nv],idx,Ndim);
        size_t ix = idx(0);
        size_t iy = idx(1);
        size_t iz = idx(2);
        #endif
        if (!IsSolid[0][ix][iy][iz])
        {
            double NonEq[Nneigh];
            double Q = 0.0;
            double tau = Tau[0];
            double rho = Rho[0][ix][iy][iz];
            Vec3_t vel = Vel[0][ix][iy][iz]+dt*tau*BForce[0][ix][iy][iz]/rho;
            double VdotV = dot(vel,vel);
            for (size_t k=0;k<Nneigh;k++)
            {
                double VdotC = dot(vel,C[k]);
                double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                NonEq[k] = F[0][ix][iy][iz][k] - Feq;
                Q +=  NonEq[k]*NonEq[k]*EEk[k];
            }
            Q = sqrt(2.0*Q);
            tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));

            bool valid = true;
            double alpha = 1.0;
            while (valid)
            {
                valid = false;
                for (size_t k=0;k<Nneigh;k++)
                {
                    Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][k] - alpha*NonEq[k]/tau;
                    if (Ftemp[0][ix][iy][iz][k]<-1.0e-12)
                    {
                        //std::cout << Ftemp[0][ix][iy][iz][k] << std::endl;
                        double temp =  tau*F[0][ix][iy][iz][k]/NonEq[k];
                        if (temp<alpha) alpha = temp;
                        valid = true;
                    }
                    if (std::isnan(Ftemp[0][ix][iy][iz][k]))
                    {
                        std::cout << "CollideSC: Nan found, resetting" << std::endl;
                        std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                        throw new Fatal("Domain::CollideSC: Distribution funcitons gave nan value, check parameters");
                    }
                }
            }
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][Op[k]];
            }
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::CollideSCDEM()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[0][ix][iy][iz])
        {
            double NonEq[Nneigh];
            double Q = 0.0;
            double tau = Tau[0];
            double gamma = Gamma[ix][iy][iz];
            double Bn = (gamma*(tau-0.5))/((1.0-gamma)+(tau-0.5));
            //double Bn = gamma;
            double rho = Rho[0][ix][iy][iz];
            Vec3_t vel = Vel[0][ix][iy][iz]+dt*tau*BForce[0][ix][iy][iz]/rho;
            double VdotV = dot(vel,vel);
            for (size_t k=0;k<Nneigh;k++)
            {
                double VdotC = dot(vel,C[k]);
                double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                NonEq[k] = F[0][ix][iy][iz][k] - Feq;
                Q +=  NonEq[k]*NonEq[k]*EEk[k];
            }
            Q = sqrt(2.0*Q);
            tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));

            bool valid = true;
            double alpha = 1.0;
            while (valid)
            {
                valid = false;
                for (size_t k=0;k<Nneigh;k++)
                {
                    double Ome = Omeis[ix][iy][iz][k];
                    double noneq = (1.0 - Bn)*NonEq[k]/tau-Bn*Ome;
                    Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][k] - alpha*(noneq);
                    if (Ftemp[0][ix][iy][iz][k]<-1.0e-12)
                    {
                        //std::cout << Ftemp[0][ix][iy][iz][k] << std::endl;
                        double temp =  F[0][ix][iy][iz][k]/noneq;
                        if (temp<alpha) alpha = temp;
                        valid = true;
                    }
                    if (std::isnan(Ftemp[0][ix][iy][iz][k]))
                    {
                        std::cout << "CollideSCDEM: Nan found, resetting" << std::endl;
                        std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << Ome << " " << Bn << " " << NonEq[k] << " " << tau << " " << F[0][ix][iy][iz][k] << " " << Rho[0][ix][iy][iz] <<  std::endl;
                        throw new Fatal("Domain::CollideMPM: Distribution functions gave nan value, check parameters");
                    }
                }
            }
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][Op[k]];
            }
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::CollideMPM()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[0][ix][iy][iz])
        {
            double NonEq[Nneigh];
            double Q = 0.0;
            double tau = Tau[0];
            //double gamma = Gamma[ix][iy][iz];
            //double Bn = (gamma*(tau-0.5))/((1.0-gamma)+(tau-0.5));
            double rho = Rho[0][ix][iy][iz];
            Vec3_t vel = Vel[0][ix][iy][iz]+dt*tau*BForce[0][ix][iy][iz]/rho;
            double VdotV = dot(vel,vel);
            for (size_t k=0;k<Nneigh;k++)
            {
                double VdotC = dot(vel,C[k]);
                double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                NonEq[k] = F[0][ix][iy][iz][k] - Feq;
                Q +=  NonEq[k]*NonEq[k]*EEk[k];
            }
            Q = sqrt(2.0*Q);
            tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));
            
            bool valid = true;
            double alpha = 1.0;
            double alphat= 1.0;
            //double Bn = gamma;
            //double Bn = 0.0;
            size_t ntry = 0;
            while (valid)
            {
                ntry++;
                valid = false;
                alpha = alphat;
                for (size_t k=0;k<Nneigh;k++)
                {
                    //double Ome = Omeis[ix][iy][iz][k];
                    //if (gamma>1.5)
                    //{ 
                        //Bn = 0.0;
                    //}
                    //double noneq = (1.0 - Bn)*NonEq[k]/tau-Ome;
                    //double noneq = (1.0 - Bn)*NonEq[k]/tau;
                    double noneq = NonEq[k]/tau;
                    Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][k] - alpha*(noneq);
                    if (Ftemp[0][ix][iy][iz][k]<-1.0e-6)
                    {
                        //std::cout << Ftemp[0][ix][iy][iz][k] << std::endl;
                        double temp =  F[0][ix][iy][iz][k]/noneq;
                        if (temp<alphat) alphat = temp;
                        valid = true;
                        if (ntry>3)
                        {
                            std::cout << "Too many tries to fix distribution function" << std::endl;
                            std::cout << "idx " << iVec3_t(ix,iy,iz) << std::endl;
                            std::cout << "k   " << k                 << C[k] << std::endl;
                            std::cout << "F+k " << Ftemp[0][ix][iy][iz][k] << std::endl;
                            std::cout << "Fk "  << F    [0][ix][iy][iz][k] << std::endl;
                            //std::cout << "Ome "  << Ome                     << std::endl;
                            std::cout << "Rho "  << Rho[0][ix][iy][iz]      << std::endl;
                            std::cout << "Vel "  << Vel[0][ix][iy][iz]      << std::endl;
                            std::cout << "Gamma "  << Gamma[ix][iy][iz]      << std::endl;
                            std::cout << "Time "   << Time      << std::endl;
                            throw new Fatal("Domain::CollideMPM: Distribution functions gave negative value");
                        }
                    }
                    if (std::isnan(Ftemp[0][ix][iy][iz][k]))
                    {
                        std::cout << "CollideMPM: Nan found, resetting" << std::endl;
                        //std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << " " << Ome << std::endl;
                        std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                        throw new Fatal("Domain::CollideMPM: Distribution funcitons gave nan value, check parameters");
                    }
                }
            }
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][Op[k]];
            }
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::CollideMRT()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[0][ix][iy][iz])
        {
            //double NonEq[Nneigh];
            //double Q = 0.0;
            //double tau = Tau[0];
            double rho = Rho[0][ix][iy][iz];
            Vec3_t vel = Vel[0][ix][iy][iz];
            //Vec3_t vel = Vel[0][ix][iy][iz]+dt*BForce[0][ix][iy][iz]/rho;
            //double VdotV = dot(vel,vel);
            //for (size_t k=0;k<Nneigh;k++)
            //{
                //double VdotC = dot(vel,C[k]);
                //double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                //NonEq[k] = F[0][ix][iy][iz][k] - Feq;
                //Q +=  NonEq[k]*NonEq[k]*EEk[k];
            //}
            //Q = sqrt(2.0*Q);
            //tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));
            double *f = F[0][ix][iy][iz];
            double *ft= Ftemp[0][ix][iy][iz];
            double fneq[Nneigh];
            int n=Nneigh,m=1;
            double a = 1.0,b = 0.0;
            dgemv_("N",&n,&n,&a,M.data,&n,f,&m,&b,ft,&m);

            ft[ 0] = 0.0; 
            ft[ 1] = S( 1)*(ft[ 1] + rho - rho*dot(vel,vel)/(Cs*Cs));
            ft[ 2] = S( 2)*(ft[ 2] - rho);
            ft[ 3] = 0.0;
            ft[ 4] = S( 4)*(ft[ 4] + 7.0/3.0*rho*vel(0)/Cs); 
            ft[ 5] = 0.0;
            ft[ 6] = S( 6)*(ft[ 6] + 7.0/3.0*rho*vel(1)/Cs); 
            ft[ 7] = 0.0;
            ft[ 8] = S( 8)*(ft[ 8] + 7.0/3.0*rho*vel(2)/Cs); 
            ft[ 9] = S( 9)*(ft[ 9] - rho*(2.0*vel(0)*vel(0)-vel(1)*vel(1)-vel(2)*vel(2))/(Cs*Cs));
            ft[10] = S(10)*(ft[10] - rho*(vel(1)*vel(1)-vel(2)*vel(2))/(Cs*Cs));
            ft[11] = S(11)*(ft[11] - rho*(vel(0)*vel(1))/(Cs*Cs));
            ft[12] = S(12)*(ft[12] - rho*(vel(1)*vel(2))/(Cs*Cs));
            ft[13] = S(13)*(ft[13] - rho*(vel(0)*vel(2))/(Cs*Cs));
            ft[14] = S(14)* ft[14];
            
            dgemv_("N",&n,&n,&a,Minv.data,&n,ft,&m,&b,fneq,&m);


            bool valid = true;
            double alpha = 1.0, alphat = 1.0;
            while (valid)
            {
                alpha = alphat;
                valid = false;
                for (size_t k=0;k<Nneigh;k++)
                {
                    Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][k] - alpha*(fneq[k] - 3.0*W[k]*dot(C[k],BForce[0][ix][iy][iz])*dt*dt/dx);
                    if (Ftemp[0][ix][iy][iz][k]<-1.0e-12)
                    {
                        //std::cout << Ftemp[0][ix][iy][iz][k] << std::endl;
                        valid = true;
                        double temp =  F[0][ix][iy][iz][k]/(fneq[k] - 3.0*W[k]*dot(C[k],BForce[0][ix][iy][iz])*dt*dt/dx);
                        if (temp<alphat) alphat = temp;
                    }
                    //if (std::isnan(Ftemp[0][ix][iy][iz][k]))
                    //{
                        //std::cout << "CollideSC: Nan found, resetting" << std::endl;
                        //std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                        //throw new Fatal("Domain::CollideSC: Distribution funcitons gave nan value, check parameters");
                    //}
                }
            }
        }
        else
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][Op[k]];
            }
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::CollideMP ()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t Vmix = OrthoSys::O;
        double den  = 0.0;
        for (size_t il=0;il<Nl;il++)
        {
            Vmix += Rho[il][ix][iy][iz]*Vel[il][ix][iy][iz]/Tau[il];
            den  += Rho[il][ix][iy][iz]/Tau[il];
        }
        Vmix /= den;

        for (size_t il=0;il<Nl;il++)
        {
            if (!IsSolid[il][ix][iy][iz])
            {
                double rho = Rho[il][ix][iy][iz];
                Vec3_t vel = Vmix + dt*Tau[il]*BForce[il][ix][iy][iz]/rho;
                //Vec3_t vel = Vel[il][ix][iy][iz] + dt*Tau[il]*BForce[il][ix][iy][iz]/rho;
                double VdotV = dot(vel,vel);
                bool valid = true;
                double alpha = 1.0;
                while (valid)
                {
                    valid = false;
                    for (size_t k=0;k<Nneigh;k++)
                    {
                        double VdotC = dot(vel,C[k]);
                        double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                        Ftemp[il][ix][iy][iz][k] = F[il][ix][iy][iz][k] - alpha*(F[il][ix][iy][iz][k] - Feq)/Tau[il];
                        if (Ftemp[il][ix][iy][iz][k]<-1.0e-12)
                        {
                            //std::cout << Ftemp[il][ix][iy][iz][k] << std::endl;
                            double temp =  Tau[il]*F[il][ix][iy][iz][k]/(F[il][ix][iy][iz][k] - Feq);
                            if (temp<alpha) alpha = temp;
                            valid = true;
                        }
                        if (std::isnan(Ftemp[il][ix][iy][iz][k]))
                        {
                            std::cout << "CollideMP: Nan found, resetting" << std::endl;
                            std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                            throw new Fatal("Domain::CollideMP: Distribution funcitons gave nan value, check parameters");
                        }
                    }
                }
            }
            else
            {
                for (size_t k=0;k<Nneigh;k++)
                {
                    Ftemp[il][ix][iy][iz][k] = F[il][ix][iy][iz][Op[k]];
                }
            }
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::CollideEFS_SRT()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[0][ix][iy][iz])
        {
            double tau   = Tau[0];
            double rho   = Rho[0][ix][iy][iz];
            Vec3_t vel   = Vel[0][ix][iy][iz];

            for (size_t k=0;k<Nneigh;k++)
            {
                double Feqtmp = Feq(k,rho,vel);
                double Ftmp   = 3.0/rho * Feqtmp*dot(BForce[0][ix][iy][iz],(C[k]-vel))/(Cs*Cs);
                Ftemp[0][ix][iy][iz][k] = 1.0/tau*Feqtmp + (1.0-0.5/tau)*dt*Ftmp + (1.0-1.0/tau)*F[0][ix][iy][iz][k];   
            }
        }

        else
        //if    (IsSolid[0][ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][Op[k]];                                        // bounce-back BC scheme --- 20170228
            }
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::StreamSC()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #ifdef IGNORESOLID
    for (size_t nv=0;nv<ValidNodes.Size();nv++)
    #else
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    #endif
    {
        #ifdef IGNORESOLID
        iVec3_t idx;
        idx2Pt(ValidNodes[nv],idx,Ndim);
        size_t ix = idx(0);
        size_t iy = idx(1);
        size_t iz = idx(2);
        #endif
        for (size_t k=0;k<Nneigh;k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);
            Ftemp[0][nix][niy][niz][k] = F[0][ix][iy][iz][k];
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #ifdef IGNORESOLID
    for (size_t nv=0;nv<ValidNodes.Size();nv++)
    #else
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    #endif
    {
        #ifdef IGNORESOLID
        iVec3_t idx;
        idx2Pt(ValidNodes[nv],idx,Ndim);
        size_t ix = idx(0);
        size_t iy = idx(1);
        size_t iz = idx(2);
        #endif
        BForce[0][ix][iy][iz] = OrthoSys::O;
        Vel   [0][ix][iy][iz] = OrthoSys::O;
        Rho   [0][ix][iy][iz] = 0.0;
        if (!IsSolid[0][ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Rho[0][ix][iy][iz] +=  F[0][ix][iy][iz][k];
                Vel[0][ix][iy][iz] +=  F[0][ix][iy][iz][k]*C[k];
            }
            Vel[0][ix][iy][iz] *= Cs/Rho[0][ix][iy][iz];
        }
    }
}

inline void Domain::StreamMPM()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        for (size_t k=0;k<Nneigh;k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);
            //if (fabs(Gamma[ix][iy][iz]-1.0)>0.5&&fabs(Gamma[nix][niy][niz]-1.0)>0.5)
            //if (fabs(Gamma[ix][iy][iz]-1.0)>0.5)
            //{
                Ftemp[0][nix][niy][niz][k] = F[0][ix][iy][iz][k];
            //}
            //else
            //{
                 //Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][k];
            //}
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}

inline void Domain::StreamMP()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    for (size_t il=0;il<Nl;il++)
    {
        #ifdef USE_OMP
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        #endif
        for (size_t ix=0;ix<nx;ix++)
        for (size_t iy=0;iy<ny;iy++)
        for (size_t iz=0;iz<nz;iz++)
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
                size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
                size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);
                Ftemp[il][nix][niy][niz][k] = F[il][ix][iy][iz][k];
            }
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        for (size_t il=0;il<Nl;il++)
        {
            Vel   [il][ix][iy][iz] = OrthoSys::O;
            Rho   [il][ix][iy][iz] = 0.0;
            if (!IsSolid[il][ix][iy][iz])
            {
                for (size_t k=0;k<Nneigh;k++)
                {
                    Rho[il][ix][iy][iz] +=  F[il][ix][iy][iz][k];
                    Vel[il][ix][iy][iz] +=  F[il][ix][iy][iz][k]*C[k];
                }
                Vel[il][ix][iy][iz] *= Cs/Rho[il][ix][iy][iz];
            }
            BForce[il][ix][iy][iz] = OrthoSys::O;
        }
    }
}

inline void Domain::StreamEFS_SRT()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        for (size_t k=0;k<Nneigh;k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);                          // periodic BC scheme 
            Ftemp[0][nix][niy][niz][k] = F[0][ix][iy][iz][k];                                            //
        }
    }

    double ***** tmp = F;
    F  = Ftemp;
    Ftemp = tmp;

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vel   [0][ix][iy][iz] = OrthoSys::O;                                                               // in fact the solid nodes are specified zero values. Good.   --- 20170302 
        Rho   [0][ix][iy][iz] = 0.0;
        if (!IsSolid[0][ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Rho[0][ix][iy][iz] +=  F[0][ix][iy][iz][k];                                               // calculating density
                Vel[0][ix][iy][iz] +=  F[0][ix][iy][iz][k]*C[k];
            }
            Vel[0][ix][iy][iz] += 0.5*dt*BForce[0][ix][iy][iz];  // calculating velocity
            Vel[0][ix][iy][iz] *= Cs/Rho[0][ix][iy][iz];
        }
        BForce[0][ix][iy][iz] = OrthoSys::O; 
    }

}

inline void Domain::VelDen()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        //BForce[0][ix][iy][iz] = OrthoSys::O;
        Vel   [0][ix][iy][iz] = OrthoSys::O;
        Rho   [0][ix][iy][iz] = 0.0;
        if (!IsSolid[0][ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Rho[0][ix][iy][iz] +=  F[0][ix][iy][iz][k];
                Vel[0][ix][iy][iz] +=  F[0][ix][iy][iz][k]*C[k];
            }
            Vel[0][ix][iy][iz] *= Cs/Rho[0][ix][iy][iz];
        }
    }
}

#ifdef USE_CUDA
inline void Domain::UpLoadDevice(size_t Nc)
{

    lbm_aux lbmaux;

    lbmaux.Nx = Ndim(0);
    lbmaux.Ny = Ndim(1);
    lbmaux.Nz = Ndim(2);
    lbmaux.Ncells    = Ncells;
    lbmaux.Nneigh    = Nneigh;
    lbmaux.NCPairs   = NCellPairs;
    lbmaux.Nl        = Nl;
    lbmaux.Gmix      = Gmix;
    lbmaux.Cs        = Cs;
    lbmaux.Sc        = Sc;
    lbmaux.dt        = dt;
    lbmaux.dx        = dx;
    lbmaux.iter      = 0;
    lbmaux.Time      = Time;
    lbmaux.g         = g;
    lbmaux.thick     = thick;
    lbmaux.sigma     = sigma;
    lbmaux.Ts        = Ts;
    lbmaux.Tl        = Tl;
    lbmaux.L         = L;

    for (size_t nn=0;nn<Nneigh;nn++)
    {
       lbmaux.C  [nn].x    = C  [nn](0); 
       lbmaux.C  [nn].y    = C  [nn](1); 
       lbmaux.C  [nn].z    = C  [nn](2); 
       lbmaux.EEk[nn]      = EEk[nn]   ; 
       lbmaux.W  [nn]      = W  [nn]   ;
       lbmaux.Op [nn]      = Op [nn]   ;
       if (S.Size()==0) continue;
       lbmaux.S  [nn]      = S  (nn)   ;
       for (size_t mm=0;mm<Nneigh;mm++)
       {
           lbmaux.M [mm + nn*Nneigh] = M   (nn,mm);
           lbmaux.Mi[mm + nn*Nneigh] = Minv(nn,mm);
       }
    }
    
    thrust::host_vector<real>  hF             (Nl*Ncells*Nneigh);
    thrust::host_vector<real>  hFtemp         (Nl*Ncells*Nneigh);
    thrust::host_vector<bool>  hIsSolid       (Nl*Ncells); 
    thrust::host_vector<real3> hBForce        (Nl*Ncells); 
    thrust::host_vector<real3> hVel           (Nl*Ncells); 
    thrust::host_vector<real>  hRho           (Nl*Ncells); 
    thrust::host_vector<uint3> hCellsPairs    (NCellPairs); 



    for (size_t il=0;il<Nl;il++)
    {
        lbmaux.Tau     [il]    = Tau     [il];
        if (Solver==NavierStokes)
        {
            lbmaux.G       [il]    = G       [il];
            lbmaux.Gs      [il]    = Gs      [il];
            lbmaux.Rhoref  [il]    = Rhoref  [il];
            lbmaux.Psi     [il]    = Psi     [il];
        }
        else if (Solver==PhaseField)
        {
            lbmaux.rho     [il]    = rho     [il];
        }
        else if (Solver==PhaseFieldIce)
        {
            lbmaux.rho     [il]    = rho     [il];
            lbmaux.cap     [il]    = cap     [il];
            lbmaux.kap     [il]    = kap     [il];
        }

        for (size_t nz=0;nz<Ndim(2);nz++)
        {
            #pragma omp parallel for schedule(static) num_threads(Nc)
            for (size_t ny=0;ny<Ndim(1);ny++)
            {
                for (size_t nx=0;nx<Ndim(0);nx++)
                {
                    size_t Nm = nx + ny*Ndim(0) + nz*Ndim(1)*Ndim(0) + il*Ncells;
                    hIsSolid[Nm]      = IsSolid  [il][nx][ny][nz];
                    hRho[Nm]          = Rho      [il][nx][ny][nz];
                    hVel[Nm]          = make_real3(Vel[il][nx][ny][nz][0],Vel[il][nx][ny][nz][1],Vel[il][nx][ny][nz][2]);
                    hBForce[Nm]       = make_real3(BForce[il][nx][ny][nz][0],BForce[il][nx][ny][nz][1],BForce[il][nx][ny][nz][2]);
                    //Nm++;
                    for (size_t nn=0;nn<Nneigh;nn++)
                    {
                        size_t Nn = nn + nx*Nneigh + ny*Ndim(0)*Nneigh + nz*Ndim(1)*Ndim(0)*Nneigh + il*Ncells*Nneigh; 
                        hF    [Nn] = F    [il][nx][ny][nz][nn];
                        hFtemp[Nn] = Ftemp[il][nx][ny][nz][nn];
                        //Nn++;
                    }
                }
            }
        }
    }

    for (size_t ic=0;ic<NCellPairs;ic++)
    {
        hCellsPairs[ic] = make_uint3(CellPairs[ic](0),CellPairs[ic](1),CellPairs[ic](2));
    }

    //blbmaux.push_back(lbmaux);

    bF         = hF         ;
    bFtemp     = hFtemp     ;
    bIsSolid   = hIsSolid   ;
    bBForce    = hBForce    ;
    bVel       = hVel       ;
    bRho       = hRho       ;
    bCellPairs = hCellsPairs;


    pF         = thrust::raw_pointer_cast(bF.data());
    pFtemp     = thrust::raw_pointer_cast(bFtemp.data());
    pIsSolid   = thrust::raw_pointer_cast(bIsSolid.data());
    pBForce    = thrust::raw_pointer_cast(bBForce.data());
    pVel       = thrust::raw_pointer_cast(bVel.data());
    pRho       = thrust::raw_pointer_cast(bRho.data());
    pCellPairs = thrust::raw_pointer_cast(bCellPairs.data());
    //plbmaux    = thrust::raw_pointer_cast(blbmaux.data());

    cudaMalloc(&plbmaux, sizeof(lbm_aux));
    cudaMemcpy(plbmaux, &lbmaux, sizeof(lbm_aux), cudaMemcpyHostToDevice);
    

    //cudaCheckUpLoad<<<1,1>>>(plbmaux);

}

inline void Domain::DnLoadDevice(size_t Nc)
{
    thrust::host_vector<real3> hVel;
    thrust::host_vector<real>  hRho;

    hVel = bVel;
    hRho = bRho;

    for (size_t il=0;il<Nl;il++)
    {
        for (size_t nz=0;nz<Ndim(2);nz++)
        {
            #pragma omp parallel for schedule(static) num_threads(Nc)
            for (size_t ny=0;ny<Ndim(1);ny++)
            {
                for (size_t nx=0;nx<Ndim(0);nx++)
                {
                    size_t Nm = nx + ny*Ndim(0) + nz*Ndim(1)*Ndim(0) + il*Ncells;
                    Rho[il][nx][ny][nz]    = hRho[Nm];
                    Vel[il][nx][ny][nz][0] = static_cast<real3>(hVel[Nm]).x;
                    Vel[il][nx][ny][nz][1] = static_cast<real3>(hVel[Nm]).y;
                    Vel[il][nx][ny][nz][2] = static_cast<real3>(hVel[Nm]).z;
                }
            }
        }
    }

    

}
#endif

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{
    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Nproc = TheNproc;

    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1    , TERM_RST);
    printf("%s  Number of cores                  =  %zd%s\n"      ,TERM_CLR2, Nproc                                 , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                    , TERM_RST);
    printf("%s  Grid size                        =  %g%s\n"       ,TERM_CLR2, dx                                    , TERM_RST);
    printf("%s  C parameter                      =  %g%s\n"       ,TERM_CLR2, dx/dt                                 , TERM_RST);
    printf("%s  Number of LBM cells              =  %zd%s\n"      ,TERM_CLR2, Ncells                                , TERM_RST);
    for (size_t i=0;i<Nl;i++)
    {
    printf("%s  Tau of Lattice %zd                 =  %g%s\n"       ,TERM_CLR2, i, Tau[i]                             , TERM_RST);
    }
    printf("%s  Simulated Time                   =  %g%s\n"       ,TERM_CLR2, Tf                                    , TERM_RST);
    #ifdef USE_CUDA
    cudaDeviceProp prop;
    cudaGetDeviceProperties (&prop,0);
    std::cout 
        << TERM_CLR2 
        << "  Using GPU:                       =  " << prop.name << TERM_RST << std::endl;
    #endif


    if (Solver==NavierStokes)
    {
        Array<iVec3_t> CPP(0);

        size_t nx = Ndim(0);
        size_t ny = Ndim(1);
        size_t nz = Ndim(2);

        for (size_t iz=0;iz<nz;iz++)
        for (size_t iy=0;iy<ny;iy++)
        for (size_t ix=0;ix<nx;ix++)
        {
            size_t nc = Pt2idx(iVec3_t(ix,iy,iz),Ndim);
            #ifdef IGNORESOLID
            if (!IsSolid[0][ix][iy][iz]) ValidNodes.Push(nc);
            bool maybevalid = false;
            #endif
            for (size_t k=1;k<Nneigh;k++)
            {
                size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
                size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
                size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);
                size_t nb  = Pt2idx(iVec3_t(nix,niy,niz),Ndim);
                #ifdef IGNORESOLID
                if (!IsSolid[0][nix][niy][niz]) maybevalid = true;
                #endif
                if ((fabs(G[0])>1.0e-12||(fabs(Gs[0])>1.0e-12)||(Nl>1))&&(nb>nc))
                {
                    CPP.Push(iVec3_t(nc,nb,k));
                }
            }
            #ifdef IGNORESOLID
            if (IsSolid[0][ix][iy][iz]&&maybevalid) ValidNodes.Push(nc);
            #endif
        }

        NCellPairs = CPP.Size();
        CellPairs = new iVec3_t [NCellPairs];
        for (size_t n=0;n<NCellPairs;n++)
        {
            CellPairs[n] = CPP[n];
        }
    }

    #ifdef IGNORESOLID
    std::cout 
        << TERM_CLR2 
        << "  Number of valid nodes            =  " << ValidNodes.Size() << TERM_RST << std::endl;
    #endif

    #ifdef USE_CUDA
    UpLoadDevice();
    #endif
    
    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time >= tout)
        {
            //std::cout << "3" << std::endl;
            #ifdef USE_CUDA
            DnLoadDevice();
            #endif
            //std::cout << "4" << std::endl;
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    WriteXDMF(fn.CStr());
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }
        #ifdef USE_CUDA
        if (Solver==NavierStokes)
        {
            if (Nl==1)
            {
                if      (fabs(G[0]) >0.0) 
                {
                    cudaApplyForcesSC<<<NCellPairs/Nthread+1,Nthread>>>(pCellPairs,pIsSolid,pBForce,pRho,plbmaux);
                    //cudaDeviceSynchronize();
                }
                //cudaCollideSC_MRT<<<Nl*Ncells/256+1,256>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
                cudaCollideSC<<<Nl*Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
                //cudaDeviceSynchronize();
            }
            else 
            {
                if ((fabs(G[0])>1.0e-12)||(fabs(G[1])>1.0e-12)||(fabs(Gmix)>1.0e-12))
                {
                    cudaApplyForcesSCMP<<<NCellPairs/Nthread+1,Nthread>>>(pCellPairs,pIsSolid,pBForce,pRho,plbmaux);
                    //cudaDeviceSynchronize();
                }

                cudaCollideMP<<<Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
                //cudaDeviceSynchronize();
            }
        }
        else if (Solver==AdvectionDiffusion)
        {
            cudaCollideAD<<<Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
        }
        else if (Solver==PhaseFieldIce)
        {
            cudaCollidePFI<<<Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
        }
        else if (Solver==PhaseField)
        {
            cudaCollidePF<<<Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
        }
        else if (Solver==ShallowWater)
        {
            cudaCollideSW<<<Nl*Ncells/Nthread+1,Nthread>>>(pF,pFtemp,pBForce,pVel,pRho,plbmaux);
        }

        real * tmp = pF;
        pF = pFtemp;
        pFtemp = tmp;
        cudaStream1<<<Nl*Ncells/Nthread+1,Nthread>>>(pF,pFtemp,pBForce,plbmaux);
        //cudaDeviceSynchronize();

        tmp = pF;
        pF = pFtemp;
        pFtemp = tmp;
        if (Solver==PhaseFieldIce)
        {
            cudaStreamPFI2<<<Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
        }
        else if (Solver==PhaseField)
        {
            cudaStreamPF2<<<Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
        }
        else if (Solver==ShallowWater)
        {
            cudaStreamSW2<<<Nl*Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
        }
        else
        {
            cudaStream2<<<Nl*Ncells/Nthread+1,Nthread>>>(pIsSolid,pF,pFtemp,pBForce,pVel,pRho,plbmaux);
        }
        //cudaDeviceSynchronize();
        #else
        if (Nl==1)
        {
            if      (fabs(G[0]) >1.0e-12) ApplyForcesSC();
            else if (fabs(Gs[0])>1.0e-12) ApplyForcesSCSS();
            CollideSC();
            //CollideMRT();
            //CollideEFS_SRT();
            StreamSC();
            //StreamEFS_SRT();
        }
        else
        {
            if ((fabs(G[0])>1.0e-12)||(fabs(G[1])>1.0e-12)) ApplyForcesSCMP();
            else ApplyForcesMP();
            CollideMP();
            StreamMP();
        }
        #endif

        Time += dt;
        //std::cout << Time << std::endl;
    }
}

const double Domain::WEIGHTSD2Q5   [ 5] = { 2./6., 1./6., 1./6., 1./6., 1./6 };
const double Domain::WEIGHTSD2Q9   [ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const double Domain::WEIGHTSD3Q15  [15] = { 2./9., 1./9., 1./9., 1./9., 1./9.,  1./9.,  1./9., 1./72., 1./72. , 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
const double Domain::WEIGHTSD3Q19  [19] = { 1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};
const size_t Domain::OPPOSITED2Q5  [ 5] = { 0, 3, 4, 1, 2 };                                                       ///< Opposite directions (D2Q5) 
const size_t Domain::OPPOSITED2Q9  [ 9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };                                           ///< Opposite directions (D2Q9) 
const size_t Domain::OPPOSITED3Q15 [15] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};                     ///< Opposite directions (D3Q15)
const size_t Domain::OPPOSITED3Q19 [19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};     ///< Opposite directions (D3Q19)


const Vec3_t Domain::LVELOCD2Q5  [ 5] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0} };
const Vec3_t Domain::LVELOCD2Q9  [ 9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
const Vec3_t Domain::LVELOCD3Q15 [15] =
{
	{ 0, 0, 0}, { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, 
	{ 0, 0, 1}, { 0, 0,-1}, { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, 
	{-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} 
};
const Vec3_t Domain::LVELOCD3Q19 [19] =
{
	{ 0, 0, 0}, 
    { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1}, 
    { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 1, 0, 1}, {-1, 0,-1},
    { 1, 0,-1}, {-1, 0, 1}, { 0, 1, 1}, { 0,-1,-1}, { 0, 1,-1}, { 0,-1, 1}
};
const double Domain::MD3Q15 [15][15]  = { { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          {-2.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          {16.0,-4.0,-4.0,-4.0,-4.0,-4.0,-4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                                          { 0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                          { 0.0,-4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 1.0,-1.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0,-4.0, 4.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0,-4.0, 4.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0},
                                          { 0.0, 2.0, 2.0,-1.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                          { 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0,-1.0},
                                          { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0,-1.0} };
}
#endif
