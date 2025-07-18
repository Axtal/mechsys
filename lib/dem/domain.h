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

/** @file dem/domain.h .*/

#ifndef MECHSYS_DEM_DOMAIN_H
#define MECHSYS_DEM_DOMAIN_H

// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <list>
#include <utility>

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// Voro++
//#include "src/voro++.cc"
#include "voro++.hh"

// VTK
#ifdef USE_VTK
//#include <vtkCellType.h>
//#include <vtkPolyhedron.h>
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkPolyData.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#endif // USE_VTK

// MechSys
#include <mechsys/dem/interacton.h>
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/tree.h>
#ifdef USE_CUDA
#include <mechsys/dem/dem.cuh>
#endif

namespace DEM
{

struct MtData;

class Domain
{
public:
    // typedefs
    typedef void (*ptFun_t) (Domain & Dom, void * UserData);

    // Constructor
    Domain(void * UserData=NULL,size_t ContactLaw=0);

    // Destructor
    ~Domain();

    // Particle generation
    void GenSpheres      (int Tag, double L, size_t N, double rho, char const * Type,
                          size_t Randomseed, double fraction, double RminFraction = 1.0);                                        ///< General spheres
    void GenSpheresBox (int Tag, Vec3_t const & X0, Vec3_t const & X1,                                                           ///< Generate spheres within a rectangular box defined by the vectors X0 and X1
                        double R, double rho, char const * Type, size_t Randomseed, double fraction, double RminFraction);
    void GenRice         (int Tag, double L, size_t N, double R, double rho, size_t Randomseed, double fraction);                ///< General rices
    void GenBox          (int InitialTag, Vec3_t & Xmin, Vec3_t & Xmax, double R, double Cf, bool Cohesion=false);               ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBox          (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, bool Cohesion=false);            ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenOpenBox      (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf);                                 ///< Generate five walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBoundingBox  (int InitialTag, double R, double Cf,bool Cohesion=false);                                              ///< Generate o bounding box enclosing the previous included particles.
    void GenBoundingPlane(int InitialTag, int GroupTag, double R, double Cf,bool Cohesion=false);                                ///< Same as GenBounding but only generates one pair of planes surrounding particles with a given Tag.
    void GenFromMesh     (Mesh::Generic & M, double R, double rho, bool cohesion=false, bool MC=true, double thickness = 0.0);   ///< Generate particles from a FEM mesh generator
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bool Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                        ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni
    void AddVoroPack(int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bVec3_t Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                     ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni, Periodic conditions are chosen for each particle
    // Single particle addition
    void AddSphere   (int Tag, Vec3_t const & X, double R, double rho);                                                          ///< Add sphere
    void AddCube     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddRecBox   (int Tag, Vec3_t const & X, Vec3_t const & L, double R, double rho, double Angle=0, Vec3_t * Axis=NULL);    ///< Add a rectangular box with dimensions given by the vector L
    void AddTetra    (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a tetrahedron at position X with spheroradius R, side of length L and density rho
    void AddDrill    (int Tag, Vec3_t const & X, double R, double Lt, double Ll, double rho);                                    ///< A drill made as a combination of a cube and a pyramid.
    void AddRice     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a rice at position X with spheroradius R, side of length L and density rho
    void AddPlane    (int Tag, Vec3_t const & X, double R, double Lx,double Ly, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddVoroCell (int Tag, voro::voronoicell & VC, double R, double rho, bool Erode, Vec3_t nv = iVec3_t(1.0,1.0,1.0));      ///< Add a single voronoi cell, it should be built before tough
    void AddTorus    (int Tag, Vec3_t const & X, Vec3_t const & N, double Rmax, double R, double rho);                           ///< Add a single torus at position X with a normal N, circunference Rmax and spheroradius R
    void AddCylinder (int Tag, Vec3_t const & X0, double R0, Vec3_t const & X1, double R1, double R, double rho);                ///< Add a cylinder formed by the connection of two circles at positions X0 and X1 and radii R0 and R1
    void AddFromJson (int Tag, char const * Filename, double R, double rho, double scale,bool Erode = false);                    ///< Add a particle generated from Json mesh
    void AddFromOBJ  (int Tag, char const * Filename, double R, double rho, double scale,bool Erode = false);                    ///< Add a particle generated from OBJ  mesh file


    // 
    //void AddParticle (DEM::Particle * Pa);                                                                                       ///< Add a particle as an exact copy of particle Pa

    // Methods
    void SetProps          (Dict & D);                                                                          ///< Set the properties of individual grains by dictionaries
    void Initialize        (double dt=0.0);                                                                     ///< Set the particles to a initial state and asign the possible insteractions
    void Solve             (double tf, double dt, double dtOut, ptFun_t ptSetup=NULL, ptFun_t ptReport=NULL,
                            char const * FileKey=NULL, bool Render=true, size_t Nproc=1,double minEkin=0.0);       ///< Run simulation the simulation up to time tf, with dt and dtOut the time and report steps. The funstion Setup and Report are used to control the workflow form outside, filekey is used to name the report files. VOut has the options 0 no visualization, 1 povray, 2 xmdf and 3 both. minEkin is a minimun of kinetic energy before the simulation stops
#ifdef USE_HDF5    
    void WriteBF           (char const * FileKey);                                                              ///< Save a h5 with branch and force information
    void WriteFrac         (char const * FileKey);                                                              ///< Save a xdmf file for fracture visualization
    void WriteXDMF         (char const * FileKey);                                                              ///< Save a xdmf file for visualization
    void Save              (char const * FileKey);                                                              ///< Save the current domain
    void Load              (char const * FileKey);                                                              ///< Load the domain form a file
#endif

    void UpdateLinkedCells ();                                                                                  ///< Update the linked cells
    void BoundingBox       (Vec3_t & minX, Vec3_t & maxX);                                                      ///< Defines the rectangular box that encloses the free particles.
    void BoundingBoxTag    (Vec3_t & minX, Vec3_t & maxX, int Tag);                                             ///< Defines the rectangular box that encloses the particles of a given Tag.
    void BoundingBoxAll    (Vec3_t & minX, Vec3_t & maxX);                                                      ///< Defines the rectangular box that encloses all the particles particles.
    void Center            (Vec3_t C = Vec3_t(0.0,0.0,0.0));                                                    ///< Centers the domain around C
    void ClearInteractons  ();                                                                                  ///< Reset the interactons
    void ResetInteractons  ();                                                                                  ///< Reset the interactons
    void ResetDisplacements();                                                                                  ///< Reset the displacements
    double MaxDisplacement ();                                                                                  ///< Calculate maximun displacement
    double MaxDim          ();                                                                                  ///< Calculate maximun size of free particles
    void ResetContacts     ();                                                                                  ///< Reset the contact list
    void UpdateContacts    ();                                                                                  ///< Update all the contact lists
    void EnergyOutput      (size_t IdxOut, std::ostream & OutFile);                                             ///< Output of the energy variables
    void GetGSD            (Array<double> & X, Array<double> & Y, Array<double> & D, size_t NDiv=10) const;     ///< Get the Grain Size Distribution
    void Clusters          ();                                                                                  ///< Check the bounded particles in the domain and how many connected clusters are still present
    void DelParticles      (Array<int> const & Tags);                                                           ///< Delete particle
    void DelParticlesIdx   (Array<int> const & idxs);                                                           ///< Delete particle by index

    // Access methods
    Particle       * GetParticle  (int Tag, bool Check=true);       ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    Particle const & GetParticle  (int Tag, bool Check=true) const; ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    void             GetParticles (int Tag, Array<Particle*> & P);  ///< Find all particles with Tag

    // Auxiliar methods
    void   LinearMomentum  (Vec3_t & L);                    ///< Return total momentum of the system
    void   AngularMomentum (Vec3_t & L);                    ///< Return total angular momentum of the system
    double CalcEnergy      (double & Ekin, double & Epot);  ///< Return total energy of the system
    double CriticalDt      ();                              ///< Calculate critical time step from existing particles


#ifdef USE_OMP
    omp_lock_t                                        lck;                         ///< to protect variables in multithreading
#endif
    Array<std::pair<size_t, size_t> >                 ListPosPairs;                ///< List of all possible particles pairs
    iVec3_t                                           LCellDim;                    ///< Dimensions of the linked cell array
    Array<Array <size_t> >                            LinkedCell;                  ///< Linked Cell array for optimization.
    Vec3_t                                            LCxmin;                      ///< Bounding box low   limit for the linked cell array
    Vec3_t                                            LCxmax;                      ///< Bounding box upper limit for the linked cell array
    Vec3_t                                            DxLC;                        ///< box sizes of the linked cells in each dimension
    std::map<std::pair<int,int>, double>              FricCoeff;                   ///< Friction coefficient map defined by user

    // Data
    bool                                              Initialized;                 ///< System (particles and interactons) initialized ?
    bool                                              RotPar;                      ///< Check if particles should be rotated, useful if particle displacements are small
    bool                                              Finished;                    ///< Has the simulation finished
    bool                                              Dilate;                      ///< True if eroded particles should be dilated for visualization
    Array<size_t>                                     FreePar;                     ///< Particles that are free
    Array<size_t>                                     NoFreePar;                   ///< Particles that are not free
    Array<Particle*>                                  Particles;                   ///< All particles in domain
    Array<Interacton*>                                Interactons;                 ///< All interactons
    Array<CInteracton*>                               CInteractons;                ///< Contact interactons
    Array<BInteracton*>                               BInteractons;                ///< Cohesion interactons
    double                                            Time;                        ///< Current time
    double                                            Dt;                          ///< Time step
    double                                            Evis;                        ///< Energy dissipated by the viscosity of the grains
    double                                            Efric;                       ///< Energy dissipated by friction
    double                                            Wext;                        ///< Work done by external forces
    double                                            Vs;                          ///< Volume occupied by the grains
    double                                            Ms;                          ///< Total mass of the particles
    double                                            Alpha;                       ///< Verlet distance
    double                                            Beta;                        ///< Binmultiplier
    double                                            Xmax;                        ///< Maximun distance along the X axis (Periodic Boundary)
    double                                            Xmin;                        ///< Minimun distance along the X axis (Periodic Boundary)
    double                                            Ymax;                        ///< Maximun distance along the Y axis (Periodic Boundary)
    double                                            Ymin;                        ///< Minimun distance along the Y axis (Periodic Boundary)
    double                                            Zmax;                        ///< Maximun distance along the Z axis (Periodic Boundary)
    double                                            Zmin;                        ///< Minimun distance along the Z axis (Periodic Boundary)
    Vec3_t                                            Per;                         ///< Periodic distance vector
    double                                            MaxDmax;                     ///< Maximun value for the radious of the spheres surronding each particle
    void *                                            UserData;                    ///< Some user data
    String                                            FileKey;                     ///< File Key for output files
    size_t                                            Nproc;                       ///< Number of cores for multithreading
    size_t                                            idx_out;                     ///< Index of output
    size_t                                            iter;                        ///< Iteration counter
    size_t                                            ContactLaw;                  ///< Contact law index                                                                                  
    std::unordered_map<size_t,size_t>                 PairtoCInt;                  ///< A map to identify which interacton has a given pair
    Array<Array <int> >                               Listofclusters;              ///< List of particles belonging to bounded clusters (applies only for cohesion simulations)
    MtData *                                          MTD;                         ///< Multithread data

    // Some utilities when the interactions are mainly between spheres
    bool                                              MostlySpheres;               ///< If the simulation is mainly between spheres this should be true
    FrictionMap_t                                     FricSpheres;                 ///< The friction value for spheres only
    FrictionMap_t                                     RollSpheres;                 ///< Map storing the rolling resistance between spheres
    void     CalcForceSphere();                                                    ///< Calculate force between only spheres spheres
    
    //CUDA Implementation
#ifdef USE_CUDA
    size_t Nthread = 256;                                                   ///< Number of GPU threads 
    void UpLoadDevice(size_t Nc=1, bool first = true);                      ///< Upload the domain to cuda device
    void DnLoadDevice(size_t Nc=1, bool force = true);                      ///< Download the key data from device
    void UpdateContactsDevice();                                            ///< Update contacts lists in device
    thrust::device_vector<ParticleCU>      bParticlesCU;                    ///< device vector of particles
    thrust::device_vector<DynParticleCU>   bDynParticlesCU;                 ///< device vector of particles
    thrust::device_vector<real>            bMaxDCU;                         ///< device vector of maximun displacements
    thrust::device_vector<real3>           bVertsCU;                        ///< device vector of particle vertices
    thrust::device_vector<real3>           bVertsoCU;                       ///< device vector of particle vertices original positions
    thrust::device_vector<size_t>          bEdgesCU;                        ///< device vector of edges indexes to vertices
    thrust::device_vector<size_t>          bFacidCU;                        ///< device vector of ordered vertex indexes belonging to a face
    thrust::device_vector<size_t>          bFacesCU;                        ///< device vector of faces indexes to vertices
    thrust::device_vector<InteractonCU>    bInteractons;                    ///< Arrays of possible interaction pairs
    thrust::device_vector<ComInteractonCU> bComInteractons;                 
    thrust::device_vector<DynInteractonCU> bDynInteractonsVV;
    thrust::device_vector<DynInteractonCU> bDynInteractonsEE;
    thrust::device_vector<DynInteractonCU> bDynInteractonsVF;
    thrust::device_vector<DynInteractonCU> bDynInteractonsFV;
    dem_aux                                demaux;                          ///< structure with auxiliary data
    // Pointers to the GPU arrays
    ParticleCU       * pParticlesCU;                                         
    DynParticleCU    * pDynParticlesCU;                                         
    real             * pMaxDCU;                                             
    real3            * pVertsCU;                                             
    real3            * pVertsoCU;                                             
    size_t           * pEdgesCU;                                             
    size_t           * pFacidCU;                                             
    size_t           * pFacesCU;                                             
    InteractonCU     * pInteractons;                  
    ComInteractonCU  * pComInteractons;                  
    DynInteractonCU  * pDynInteractonsVV;
    DynInteractonCU  * pDynInteractonsEE;
    DynInteractonCU  * pDynInteractonsVF;
    DynInteractonCU  * pDynInteractonsFV;
    dem_aux          * pdemaux;
    // Pointers to the main function calculating forces and dynamics
    ForceVV_ptr_t      pForceVV;                                            ///< pointer to the force calculation for vertex vertex
    ForceEE_ptr_t      pForceEE;                                            ///< pointer to the force calculation for edge edge
    ForceVF_ptr_t      pForceVF;                                            ///< pointer to the force calculation for vertex face
    ForceFV_ptr_t      pForceFV;                                            ///< pointer to the force calculation for face vertex
    Translate_ptr_t    pTranslate;                                          ///< pointer to the translation function
    Rotate_ptr_t       pRotate;                                             ///< pointer to the rotation function
    Reset_ptr_t        pReset;                                              ///< pointer to the reset function
    void             * pExtraParams;                                        ///< pointer to a structure of extra parameters for force calculation
#endif
    
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

struct MtData   /// A structure for the multi-thread data
{
    size_t                        ProcRank; ///< Rank of the thread
    size_t                          N_Proc; ///< Total number of threads
    DEM::Domain *                      Dom; ///< Pointer to the lbm domain
    double                             Dmx; ///< Maximun displacement
    Array<std::pair<size_t,size_t> >    LC; ///< A temporal list of new contacts
    Array<size_t>                      LCI; ///< A temporal array of posible Cinteractions
    Array<size_t>                      LCB; ///< A temporal array of posible Binteractions
    Array<std::pair<size_t,size_t> >  CLCI; ///< An array of confirmed pairs
    Array<std::pair<size_t,size_t> >   LPC; ///< A temporal list of new contacts for periodic boundary conditions
    Array<size_t>                     LPCI; ///< A temporal array of possible Cinteractions for periodic boundary conditions
    Array<std::pair<size_t,size_t> > CLPCI; ///< An array of confirmed Cinteractions for periodic boundary conditions
    Array<size_t>                      LBP; ///< A temporal array of possible boundary particles
    Array<std::pair<iVec3_t,size_t> >  LLC; ///< A temporal array of possible linked cells locations
    Array<std::pair<size_t,size_t> >   LPP; ///< A temporal array of possible particle types
};

// Constructor & Destructor

inline Domain::Domain (void * UD, size_t contactlaw)
{
    Initialized = false;
    Dilate = false;
    RotPar = true;
    Time = 0.0;
    iter = 0;
    Alpha = 0.05;
    Beta  = 2.0;
    UserData = UD;
    MostlySpheres = false;
    Xmax = Xmin = Ymax = Ymin = Zmax = Zmin =0.0;
    ContactLaw = contactlaw;
#ifdef USE_OMP
    omp_init_lock(&lck);
#endif
#ifdef USE_CUDA
    if (contactlaw==0) pForceVV   = CalcForceVV;
    if (contactlaw==1) pForceVV   = CalcForceVV_Hertz;
    pForceEE   = CalcForceEE;
    pForceVF   = CalcForceVF;
    pForceFV   = CalcForceFV;
    pTranslate = Translate;
    pRotate    = Rotate;
    pReset     = Reset;
    pExtraParams = NULL;
#endif
}

inline Domain::~Domain ()
{
    //for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    //for (size_t i=0; i<CInteractons.Size(); ++i) if (CInteractons[i]!=NULL) delete CInteractons[i];
    //for (size_t i=0; i<BInteractons.Size(); ++i) if (BInteractons[i]!=NULL) delete BInteractons[i];
}

//All the methods for particle generation

#include<mechsys/dem/dompargen.h>

// Methods

inline void Domain::SetProps (Dict & D)
{
    for (size_t i =0 ; i<Particles.Size(); i++)
    {
        for (size_t j=0; j<D.Keys.Size(); ++j)
        {
            int tag = D.Keys[j];
            if (tag==Particles[i]->Tag)
            {
                SDPair const & p = D(tag);
                if (p.HasKey("Gn"))
                {
                    Particles[i]->Props.Gn = p("Gn");
                }
                if (p.HasKey("Gt"))
                {
                    Particles[i]->Props.Gt = p("Gt");
                }
                if (p.HasKey("Gv"))
                {
                    Particles[i]->Props.Gv = p("Gv");
                }
                if (p.HasKey("Gm"))
                {
                    Particles[i]->Props.Gm = p("Gm");
                }
                if (p.HasKey("Kn"))
                {
                    Particles[i]->Props.Kn = p("Kn");
                }
                if (p.HasKey("Kt"))
                {
                    Particles[i]->Props.Kt = p("Kt");
                }
                if (p.HasKey("Bn"))
                {
                    Particles[i]->Props.Bn = p("Bn");
                }
                if (p.HasKey("Bt"))
                {
                    Particles[i]->Props.Bt = p("Bt");
                }
                if (p.HasKey("Bm"))
                {
                    Particles[i]->Props.Bm = p("Bm");
                }
                if (p.HasKey("Mu"))
                {
                    Particles[i]->Props.Mu = p("Mu");
                }
                if (p.HasKey("Eps"))
                {
                    Particles[i]->Props.eps = p("Eps");
                }
                if (p.HasKey("Beta"))
                {
                    Particles[i]->Props.Beta = p("Beta");
                }
                if (p.HasKey("Eta"))
                {
                    Particles[i]->Props.Eta = p("Eta");
                }
            }
        }
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        BInteractons[i]->UpdateParameters();
    }
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        CInteractons[i]->UpdateParameters();
    }
}

inline void Domain::Initialize (double dt)
{
    if (!Initialized)
    {
        // initialize all particles
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->Initialize(i);
            Particles[i]->InitializeVelocity(dt);
        }
        //Initializing the energies
        Evis = 0.0;
        Efric = 0.0;
        Wext = 0.0;

        // initialize
        //ResetInteractons();
        // info
        Util::Stopwatch stopwatch;
        printf("\n%s--- Initializing particles ------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
        // set flag
        Initialized = true;

        // info
        double Ekin, Epot, Etot;
        Etot = CalcEnergy (Ekin, Epot);
        printf("%s  Kinematic energy   = %g%s\n",TERM_CLR4, Ekin, TERM_RST);
        printf("%s  Potential energy   = %g%s\n",TERM_CLR4, Epot, TERM_RST);
        printf("%s  Total energy       = %g%s\n",TERM_CLR2, Etot, TERM_RST);
    }
    else
    {
        for (size_t i=0; i<Particles.Size(); i++)
        {
            if (Particles[i]->vxf) Particles[i]->xb(0) = Particles[i]->x(0) - Particles[i]->v(0)*dt;
            if (Particles[i]->vyf) Particles[i]->xb(1) = Particles[i]->x(1) - Particles[i]->v(1)*dt;
            if (Particles[i]->vzf) Particles[i]->xb(2) = Particles[i]->x(2) - Particles[i]->v(2)*dt;
        }
    }

}

inline void Domain::Solve (double tf, double dt, double dtOut, ptFun_t ptSetup, ptFun_t ptReport, char const * TheFileKey, bool Render, size_t TheNproc, double minEkin)
{
    // Assigning some domain particles especifically to the output
    FileKey.Printf("%s",TheFileKey);
    idx_out = 0;
    
    //Assigning the value for the time step to the domain variable
    Dt = dt;
    Nproc = TheNproc;

    // initialize particles
    Initialize (Dt);


    // calc the total volume of particles (solids)
    FreePar.Resize(0);
    NoFreePar.Resize(0);
    Vs = 0.0;
    Ms = 0.0;
    MaxDmax        =  0.0;
    double MaxKn   =  0.0;
    double MaxBn   =  0.0;
    double MinDmax = -1.0;
    double MinMass = -1.0;

    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        if (Particles[i]->IsFree())
        {
            Vs += Particles[i]->Props.V;
            Ms += Particles[i]->Props.m;
            if (Particles[i]->Dmax     > MaxDmax) MaxDmax = Particles[i]->Dmax;
            if (Particles[i]->Props.Kn > MaxKn  ) MaxKn   = Particles[i]->Props.Kn;
            if (Particles[i]->Dmax     < MinDmax||(MinDmax<0.0)) MinDmax = Particles[i]->Dmax;
            if (ContactLaw==0)
            {
                if (Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = Particles[i]->Props.m;
            }
            else if (ContactLaw==1)
            {
                if (Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = Particles[i]->Props.m/Particles[i]->Dmax;
            }
            FreePar.Push(i);
        }
        else NoFreePar.Push(i);
        Particles[i]->Index = i;
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        double pbn = std::max(BInteractons[i]->Bn/BInteractons[i]->L0,BInteractons[i]->Bt/BInteractons[i]->L0);
        if (pbn > MaxBn) MaxBn = pbn;
    }


    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    printf("%s  Total mass   of free particles            =  %g%s\n"        ,TERM_CLR4, Ms                                   , TERM_RST);
    printf("%s  Total volume of free particles            =  %g%s\n"        ,TERM_CLR4, Vs                                   , TERM_RST);
    printf("%s  Total number of particles                 =  %zd%s\n"       ,TERM_CLR2, Particles.Size()                     , TERM_RST);
    printf("%s  Time step                                 =  %g%s\n"        ,TERM_CLR2, dt                                   , TERM_RST);
    printf("%s  Simulated Time                            =  %g%s\n"        ,TERM_CLR2, tf                                   , TERM_RST);
    printf("%s  Verlet distance                           =  %g%s\n"        ,TERM_CLR2, Alpha                                , TERM_RST);
    printf("%s  Suggested Time Step                       =  %g%s\n"        ,TERM_CLR5, 0.1*sqrt(MinMass/(MaxKn+MaxBn))      , TERM_RST);
    if (Xmax-Xmin>1.0e-12)
    printf("%s  Periodic Boundary conditions in X between =  %g and %g%s\n" ,TERM_CLR5, Xmin, Xmax                           , TERM_RST);
    if (Ymax-Ymin>1.0e-12)
    printf("%s  Periodic Boundary conditions in Y between =  %g and %g%s\n" ,TERM_CLR5, Ymin, Ymax                           , TERM_RST);
    if (Zmax-Zmin>1.0e-12)
    printf("%s  Periodic Boundary conditions in Y between =  %g and %g%s\n" ,TERM_CLR5, Zmin, Zmax                           , TERM_RST);


    if (Alpha > 2.0*Beta*MaxDmax)
    {
        Alpha = 2.0*Beta*MaxDmax;
        printf("%s  Verlet distance changed to                =  %g%s\n"   ,TERM_CLR2, Alpha                                    , TERM_RST);
    }
    fflush(stdout); 

    // solve
    double t0   = Time;     // initial time
    double tout = t0; // time position for output

    Finished = false;

    Per = Vec3_t(Xmax-Xmin,Ymax-Ymin,Zmax-Zmin);

    // string to output energy data, if user gives the FileKey
    std::ostringstream oss_energy; 
    EnergyOutput (idx_out, oss_energy);

    MTD = new DEM::MtData[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = this;
        MTD[i].Dmx      = 0.0;
    }
    
    UpdateContacts();
#ifdef USE_CUDA
    //if (iter==0) UpLoadDevice(Nproc,true);
    //else         UpLoadDevice(Nproc,false);
    UpLoadDevice(Nproc,true);
    cudaDeviceProp prop;
    cudaGetDeviceProperties (&prop,0);
    std::cout 
        << TERM_CLR2 
        << "  Using GPU:                                =  " << prop.name       << TERM_RST << std::endl
        << "  Using Number of CUDA threads:             =  " << Nthread         << TERM_RST << std::endl
        << "  Number of vertices                        =  " << demaux.nverts   << TERM_RST << std::endl
        << "  Number of edges                           =  " << demaux.nedges/2 << TERM_RST << std::endl
        << "  Number of faces                           =  " << demaux.nfaces/2 << TERM_RST << std::endl;
#endif
    size_t iter_b = iter;
    size_t iter_t = 0;
    size_t numup  = 0;
    // run
    while (Time<tf)
    {

        // output
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time>=tout)
        {
#ifdef USE_CUDA
            DnLoadDevice(Nproc,true);
#endif
            double Ekin,Epot;
            CalcEnergy(Ekin,Epot);
            if (Ekin<minEkin&&Time>0.1*tf)
            {
                printf("\n%s--- Minimun energy reached ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
                break;
            }
            if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            if (TheFileKey!=NULL)
            {
                String fn,fb;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                fb.Printf    ("%s_bf_%04d", TheFileKey, idx_out);
                WriteXDMF    (fn.CStr());
                WriteBF      (fb.CStr());
                //EnergyOutput (idx_out, oss_energy);
            }
            //if (BInteractons.Size()>3) Clusters();
            idx_out++;
            tout += dtOut;
        }
#ifdef USE_CUDA
        //std::cout << "1" << std::endl;
        //Initialize particles
        pReset<<<(demaux.nparts+demaux.ncoint)/Nthread+1,Nthread>>>(pParticlesCU,pDynParticlesCU,pInteractons,pComInteractons,pdemaux,pExtraParams);
        //cudaDeviceSynchronize();
        //Calculate forces
        pForceVV<<<demaux.nvvint/Nthread+1,Nthread>>>(pInteractons, pComInteractons, pDynInteractonsVV, pParticlesCU, pDynParticlesCU, pdemaux, pExtraParams);
        //std::cout << "2" << std::endl;
        //cudaDeviceSynchronize();
        pForceEE<<<demaux.neeint/Nthread+1,Nthread>>>(pEdgesCU,pVertsCU,pInteractons, pComInteractons, pDynInteractonsEE, pParticlesCU, pDynParticlesCU, pdemaux, pExtraParams);
        //std::cout << "3" << std::endl;
        //cudaDeviceSynchronize();
        pForceVF<<<demaux.nvfint/Nthread+1,Nthread>>>(pFacesCU,pFacidCU,pVertsCU,pInteractons, pComInteractons, pDynInteractonsVF, pParticlesCU, pDynParticlesCU, pdemaux, pExtraParams);
        //std::cout << "4" << std::endl;
        //cudaDeviceSynchronize();
        pForceFV<<<demaux.nfvint/Nthread+1,Nthread>>>(pFacesCU,pFacidCU,pVertsCU,pInteractons, pComInteractons, pDynInteractonsFV, pParticlesCU, pDynParticlesCU, pdemaux, pExtraParams);
        //std::cout << "5" << std::endl;
        //cudaDeviceSynchronize();
        //Move Particles
        pTranslate<<<demaux.nparts/Nthread+1,Nthread>>>(pVertsCU, pParticlesCU, pDynParticlesCU, pdemaux, pExtraParams);
        //std::cout << "6" << std::endl;
        //cudaDeviceSynchronize();
        pRotate   <<<demaux.nparts/Nthread+1,Nthread>>>(pVertsCU, pParticlesCU, pDynParticlesCU, pdemaux, pExtraParams);
        //std::cout << "7" << std::endl;
        //cudaDeviceSynchronize();
        MaxD     <<<demaux.nverts/Nthread+1,Nthread>>>(pVertsCU, pVertsoCU, pMaxDCU, pdemaux);
        //std::cout << "8" << std::endl;
        //cudaDeviceSynchronize();
        real maxdis = 0.0;
        thrust::device_vector<real>::iterator it=thrust::max_element(bMaxDCU.begin(),bMaxDCU.end());
        maxdis = *it;
        //std::cout << "9" << std::endl;
        //std::cout << iter << std::endl;
        //Update Pair lists
        demaux.iter++;
        demaux.Time += Dt;
        if (maxdis>Alpha)
        {
            //std::cout << iter << std::endl;
            numup++;
            iter_t+= iter - iter_b;
            iter_b = iter;
            UpdateContactsDevice();
            //cudaDeviceSynchronize();
        }
#else
        //Initialize particles
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<Particles.Size(); i++)
        {
            // set the force and torque to the fixed values
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Tf;

            //Particles[i]->Bdry = false;
        }

        //Calculate forces
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<Interactons.Size(); i++)
        {
            //std::cout << Interactons[i]->I1 << " " << Interactons[i]->I2 << std::endl;
            //std::cout << Interactons[i]->P1->Index << " " << Interactons[i]->P2->Index << std::endl;
            
            //DEM::Interacton * p = Interactons[i];
            //std::cout << p->I1 << " " << p->I2 << std::endl;
		    if (Interactons[i]->CalcForce(Dt,Per,iter,ContactLaw))
            {
                String f_error(FileKey+"_error");
                Save     (f_error.CStr());
                WriteXDMF(f_error.CStr());
                std::cout << "Maximun overlap detected between particles at time " << Time << std::endl;
                std::cout << "Iteration number                                   " << iter << std::endl;
                sleep(1);
                throw new Fatal("Maximun overlap detected between particles");
            }
            omp_set_lock  (&Interactons[i]->P1->lck);
            Interactons[i]->P1->F += Interactons[i]->F1;
            Interactons[i]->P1->T += Interactons[i]->T1;
            omp_unset_lock(&Interactons[i]->P1->lck);
            omp_set_lock  (&Interactons[i]->P2->lck);
            Interactons[i]->P2->F += Interactons[i]->F2;
            Interactons[i]->P2->T += Interactons[i]->T2;
            omp_unset_lock(&Interactons[i]->P2->lck);

        }

        if(MostlySpheres) CalcForceSphere();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Nproc;i++)
        {
            MTD[i].Dmx = 0.0;
        }

        if (RotPar)
        {
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<Particles.Size(); i++)
            {
		        Particles[i]->Translate(Dt);
		        Particles[i]->Rotate(Dt);
                if (Particles[i]->MaxDisplacement()>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = Particles[i]->MaxDisplacement();
            }
        }
        else
        {
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<Particles.Size(); i++)
            {
		        Particles[i]->Translate(Dt);
                if (Particles[i]->MaxDisplacement()>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = Particles[i]->MaxDisplacement();
            }
        }

        double maxdis = 0.0;
        for (size_t i=0;i<Nproc;i++)
        {
            if (maxdis<MTD[i].Dmx) maxdis = MTD[i].Dmx;
        }

        if (maxdis>Alpha)
        {
            UpdateContacts();
        }
#endif
        Time += Dt;
        iter++;
    }

    // last output
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);
#ifdef USE_CUDA
    DnLoadDevice(Nproc,true);
#endif

    // save energy data
    //if (TheFileKey!=NULL)
    //{
        //String fn;
        //fn.Printf("%s_energy.res",TheFileKey);
        //std::ofstream fe(fn.CStr());
        //fe << oss_energy.str();
        //fe.close();
    //}

    // info
    double Ekin, Epot, Etot;
    Etot = CalcEnergy (Ekin, Epot);
    printf("%s  Total number of iterations                                = %zd%s\n",TERM_CLR4,iter, TERM_RST);
    printf("%s  Average number of iterations between contact list updates = %zd%s\n",TERM_CLR4,iter_t/numup, TERM_RST);
    printf("%s  Kinematic energy                                          = %g%s\n",TERM_CLR4, Ekin, TERM_RST);
    printf("%s  Potential energy                                          = %g%s\n",TERM_CLR4, Epot, TERM_RST);
    printf("%s  Total energy                                              = %g%s\n",TERM_CLR2, Etot, TERM_RST);
}

inline void Domain::WriteBF (char const * FileKey)
{

    size_t n_fn = 0;
    size_t n_rl = 0;

    for (size_t i=0;i<CInteractons.Size();i++)
    {
        //if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree())) n_fn++;
        if (norm(CInteractons[i]->Fnet)>1.0e-12)
        {
            n_fn++;
            if (CInteractons[i]->P1->Verts.Size()==1&&CInteractons[i]->P2->Verts.Size()==1) n_rl++;
        }
    }

    if (n_fn==0) return;
    
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    float  *  Fnnet = new float[3*n_fn];
    float  *  Ftnet = new float[3*n_fn];
    float  *  Froll = new float[3*n_fn];
    float  * Branch = new float[3*n_fn];
    float  *   Orig = new float[3*n_fn];
    int    *    ID1 = new   int[  n_fn];
    int    *    ID2 = new   int[  n_fn];

    size_t idx = 0;

    // Saving Collision forces
    for (size_t i=0;i<CInteractons.Size();i++)
    {
        //if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree()))
        if (norm(CInteractons[i]->Fnet)>1.0e-12)
        {
            Fnnet [3*idx  ] = float(CInteractons[i]->Fnet  (0));
            Fnnet [3*idx+1] = float(CInteractons[i]->Fnet  (1));
            Fnnet [3*idx+2] = float(CInteractons[i]->Fnet  (2));
            Ftnet [3*idx  ] = float(CInteractons[i]->Ftnet (0));
            Ftnet [3*idx+1] = float(CInteractons[i]->Ftnet (1));
            Ftnet [3*idx+2] = float(CInteractons[i]->Ftnet (2));
            if (n_rl>0)
            {
            Froll [3*idx  ] = float(CInteractons[i]->Fn    (0));
            Froll [3*idx+1] = float(CInteractons[i]->Fn    (1));
            Froll [3*idx+2] = float(CInteractons[i]->Fn    (2));
            }
            //Branch[3*idx  ] = float(CInteractons[i]->P1->x(0)-CInteractons[i]->P2->x(0));
            //Branch[3*idx+1] = float(CInteractons[i]->P1->x(1)-CInteractons[i]->P2->x(1)); 
            //Branch[3*idx+2] = float(CInteractons[i]->P1->x(2)-CInteractons[i]->P2->x(2)); 
            Branch[3*idx  ] = float(CInteractons[i]->Branch(0));
            Branch[3*idx+1] = float(CInteractons[i]->Branch(1)); 
            Branch[3*idx+2] = float(CInteractons[i]->Branch(2)); 
            //Orig  [3*idx  ] = 0.5*float(CInteractons[i]->P1->x(0)+CInteractons[i]->P2->x(0));
            //Orig  [3*idx+1] = 0.5*float(CInteractons[i]->P1->x(1)+CInteractons[i]->P2->x(1)); 
            //Orig  [3*idx+2] = 0.5*float(CInteractons[i]->P1->x(2)+CInteractons[i]->P2->x(2)); 
            Orig  [3*idx  ] = float(CInteractons[i]->P2->x(0));
            Orig  [3*idx+1] = float(CInteractons[i]->P2->x(1)); 
            Orig  [3*idx+2] = float(CInteractons[i]->P2->x(2)); 
            ID1   [idx]     = int  (CInteractons[i]->P1->Index);
            ID2   [idx]     = int  (CInteractons[i]->P2->Index);
            idx++;
        }
    }

    hsize_t dims[1];
    dims[0] = 3*n_fn;
    String dsname;
    dsname.Printf("Normal");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Fnnet );
    dsname.Printf("Tangential");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ftnet );
    if (n_rl>0)
    {
    dsname.Printf("Rolling");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Froll );
    }
    dsname.Printf("Branch");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Branch);
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Orig);
    dims[0] = n_fn;
    dsname.Printf("ID1");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,ID1   );
    dsname.Printf("ID2");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,ID2   );


    delete [] Fnnet;
    delete [] Ftnet;
    delete [] Froll;
    delete [] Branch;
    delete [] Orig;
    delete [] ID1;
    delete [] ID2;

    //Saving Cohesive forces
    if (BInteractons.Size()>0)
    {
    float  *  Bnnet = new float[3*BInteractons.Size()];
    float  *  Btnet = new float[3*BInteractons.Size()];
    float  *  BOrig = new float[3*BInteractons.Size()];
    int    *  BVal  = new   int[  BInteractons.Size()];
    int    *  BID1  = new   int[  BInteractons.Size()];
    int    *  BID2  = new   int[  BInteractons.Size()];

    idx = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        Bnnet [3*idx  ] = float(BInteractons[i]->Fnet  (0));
        Bnnet [3*idx+1] = float(BInteractons[i]->Fnet  (1));
        Bnnet [3*idx+2] = float(BInteractons[i]->Fnet  (2));
        Btnet [3*idx  ] = float(BInteractons[i]->Ftnet (0));
        Btnet [3*idx+1] = float(BInteractons[i]->Ftnet (1));
        Btnet [3*idx+2] = float(BInteractons[i]->Ftnet (2));
        BOrig [3*idx  ] = float(BInteractons[i]->xnet(0));
        BOrig [3*idx+1] = float(BInteractons[i]->xnet(1)); 
        BOrig [3*idx+2] = float(BInteractons[i]->xnet(2)); 
        BVal  [idx]     = int  (BInteractons[i]->valid);
        BID1  [idx]     = int  (BInteractons[i]->P1->Index);
        BID2  [idx]     = int  (BInteractons[i]->P2->Index);
        idx++;

    }
    hsize_t dims[1];
    dims[0] = 3*BInteractons.Size();
    String dsname;
    dsname.Printf("BNormal");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Bnnet );
    dsname.Printf("BTangential");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Btnet );
    dsname.Printf("BPosition");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,BOrig );
    dims[0] = BInteractons.Size();
    dsname.Printf("BVal");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BVal  );
    dsname.Printf("BID1");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BID1  );
    dsname.Printf("BID2");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BID2  );

    delete [] Bnnet;
    delete [] Btnet;
    delete [] BOrig;
    delete [] BVal ;
    delete [] BID1 ;
    delete [] BID2 ;


    }


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"BranchForce\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << n_fn << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << n_fn << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Normal\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Normal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tangential\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tangential \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    if (n_rl>0)
    {
    oss << "     <Attribute Name=\"Rolling\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Rolling \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    }
    oss << "     <Attribute Name=\"Branch\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Branch \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ID1\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ID1 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ID2\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ID2 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    if (BInteractons.Size()>0)
    {
    oss << "   <Grid Name=\"CohesiveForce\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << BInteractons.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << BInteractons.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/BPosition \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"BNormal\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BNormal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BTangential\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BTangential \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BVal\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BVal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BID1\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BID1 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BID2\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BID2 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::WriteXDMF (char const * FileKey)
{
    size_t N_Faces = 0;
    size_t N_Verts = 0;
    size_t N_Edges = 0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        if (Particles[i]->Edges.Size()==1) N_Edges++;
        for (size_t j=0;j<Particles[i]->Faces.Size();j++)
        {
            N_Faces += Particles[i]->Faces[j]->Edges.Size();
        }
        N_Verts += Particles[i]->Verts.Size() + Particles[i]->Faces.Size();
    }

    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (N_Faces>0)
    {

        //Geometric information
        float  * Verts   = new float [3*N_Verts];
        int    * FaceCon = new int   [3*N_Faces];
        
        //Atributes
        int    * Tags    = new int   [  N_Faces];
        int    * Clus    = new int   [  N_Faces];
        float  * Vel     = new float [  N_Faces];
        float  * Ome     = new float [  N_Faces];
        //float  * Stress  = new float [9*N_Faces];

        size_t n_verts = 0;
        size_t n_faces = 0;
        size_t n_attrs = 0;
        //size_t n_attrv = 0;
        //size_t n_attrt = 0;
        for (size_t i=0;i<Particles.Size();i++)
        {
            Particle * Pa = Particles[i];
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(Pa->Verts.Size());
            Array<Vec3_t> Vres (Pa->Verts.Size());
            for (size_t j=0;j<Pa->Verts.Size();j++)
            {
                Vtemp[j] = *Pa->Verts[j];
                Vres [j] = *Pa->Verts[j];
            }
            double multiplier = 0.0;
            if (Dilate&&Pa->Eroded&&Pa->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,Pa->EdgeCon,Pa->FaceCon,Vres,Pa->Props.R);
                multiplier = 1.0;
            }
            for (size_t j=0;j<Pa->Verts.Size();j++)
            {
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(0);
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(1);
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(2);
                Verts[n_verts++] = float(Vres[j](0));
                Verts[n_verts++] = float(Vres[j](1));
                Verts[n_verts++] = float(Vres[j](2));
            }
            size_t n_reff = n_verts/3;
            for (size_t j=0;j<Pa->FaceCon.Size();j++)
            {
                Vec3_t C,N;
                Pa->Faces[j]->Centroid(C);
                Pa->Faces[j]->Normal(N);
                Verts[n_verts++] = float(C(0) + multiplier*Pa->Props.R*N(0));
                Verts[n_verts++] = float(C(1) + multiplier*Pa->Props.R*N(1));
                Verts[n_verts++] = float(C(2) + multiplier*Pa->Props.R*N(2));
                //Verts[n_verts++] = (float) C(0);
                //Verts[n_verts++] = (float) C(1);
                //Verts[n_verts++] = (float) C(2);
                for (size_t k=0;k<Pa->FaceCon[j].Size();k++)
                {
                    size_t nin = Pa->FaceCon[j][k];
                    size_t nen = Pa->FaceCon[j][(k+1)%Pa->FaceCon[j].Size()];
                    FaceCon[n_faces++] = int(n_reff + j);  
                    FaceCon[n_faces++] = int(n_refv + nin);
                    FaceCon[n_faces++] = int(n_refv + nen);

                    //Writing the attributes
                    Tags  [n_attrs] = int(Pa->Tag);
                    Clus  [n_attrs] = size_t(Pa->Cluster);
                    Vel   [n_attrs] = float(norm(Pa->v));
                    Ome   [n_attrs] = float(norm(Pa->w));
                    n_attrs++;

                    //Vel [n_attrv  ] = (float) Pa->v(0);
                    //Vel [n_attrv+1] = (float) Pa->v(1);
                    //Vel [n_attrv+2] = (float) Pa->v(2);
                    //n_attrv += 3;

                    //Stress[n_attrt  ] = (float) Pa->M(0,0);
                    //Stress[n_attrt+1] = (float) Pa->M(1,0);
                    //Stress[n_attrt+2] = (float) Pa->M(2,0);
                    //Stress[n_attrt+3] = (float) Pa->M(0,1);
                    //Stress[n_attrt+4] = (float) Pa->M(1,1);
                    //Stress[n_attrt+5] = (float) Pa->M(2,1);
                    //Stress[n_attrt+6] = (float) Pa->M(0,2);
                    //Stress[n_attrt+7] = (float) Pa->M(1,2);
                    //Stress[n_attrt+8] = (float) Pa->M(2,2);
                    //n_attrt += 9;
                }
            }
        }

        //Write the data
        hsize_t dims[1];
        String dsname;
        dims[0] = 3*N_Verts;
        dsname.Printf("Verts");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Verts);
        dims[0] = 3*N_Faces;
        dsname.Printf("FaceCon");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,FaceCon);
        dims[0] = N_Faces;
        dsname.Printf("Tag");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tags   );
        dims[0] = N_Faces;
        dsname.Printf("Cluster");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Clus   );
        dims[0] = N_Faces;
        dsname.Printf("Velocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vel);
        dims[0] = N_Faces;
        dsname.Printf("AngVelocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ome);
        
        //Erasing the data
        delete [] Verts;
        delete [] FaceCon;
        delete [] Tags;
        delete [] Clus;
        delete [] Vel;
        delete [] Ome;
        //delete [] Stress;
    }
    // Storing center of mass data
    
    float * Radius = new float[  Particles.Size()];
    float * Posvec = new float[3*Particles.Size()];
    float * Rod    = new float[3*Particles.Size()];
    float * Velvec = new float[3*Particles.Size()];
    float * Omevec = new float[3*Particles.Size()];
    float * Aomvec = new float[3*Particles.Size()];
    float * Amovec = new float[3*Particles.Size()];
    float * Ufovec = new float[3*Particles.Size()];
    float * Hfovec = new float[3*Particles.Size()];
    float * Inertm = new float[6*Particles.Size()];
    float * Ekin   = new float[  Particles.Size()];
    int   * Tag    = new int  [  Particles.Size()];

    for (size_t i=0;i<Particles.Size();i++)
    {
        Mat3_t Inertia,Inertiar,R,Rt,t;
        Inertia(0,0) = Particles[i]->I(0); Inertia(0,1) = 0.0; Inertia(0,2) = 0.0;
        Inertia(1,1) = Particles[i]->I(1); Inertia(1,0) = 0.0; Inertia(1,2) = 0.0;
        Inertia(2,2) = Particles[i]->I(2); Inertia(2,0) = 0.0; Inertia(2,1) = 0.0;

        RotationMatrix(Particles[i]->Q,R);
        Rt = ~R;

        Mult(R,Inertia,t);
        Mult(t,Rt,Inertiar);

        Vec3_t Ao,Ome,L,t1,t2;
        Rotation(Particles[i]->w,Particles[i]->Q,Ome);
        Rotation(Particles[i]->wa,Particles[i]->Q,Ao);
        t1 = Particles[i]->I(0)*Particles[i]->w(0),Particles[i]->I(1)*Particles[i]->w(1),Particles[i]->I(2)*Particles[i]->w(2);
        Rotation (t1,Particles[i]->Q,t2);
        L = Particles[i]->Props.m*cross(Particles[i]->x,Particles[i]->v)+t2;


        Particles[i]->Verts.Size()<=2? Radius[i] = float(Particles[i]->Props.R) : Radius[i] = 0.0;
        Posvec[3*i  ] = float(Particles[i]->x(0));
        Posvec[3*i+1] = float(Particles[i]->x(1));
        Posvec[3*i+2] = float(Particles[i]->x(2));
        if (Particles[i]->Edges.Size()==1)
        {
        Vec3_t V1 = *Particles[i]->Verts[1];
        Vec3_t V0 = *Particles[i]->Verts[0];
        Rod   [3*i  ] = float(V1(0) - V0(0));
        Rod   [3*i+1] = float(V1(1) - V0(1));
        Rod   [3*i+2] = float(V1(2) - V0(2));
        }
        else
        {
        Rod   [3*i  ] = 0.0;
        Rod   [3*i+1] = 0.0;
        Rod   [3*i+2] = 0.0;
        }
        Velvec[3*i  ] = float(Particles[i]->v(0));
        Velvec[3*i+1] = float(Particles[i]->v(1));
        Velvec[3*i+2] = float(Particles[i]->v(2));
        Omevec[3*i  ] = float(Ome(0));
        Omevec[3*i+1] = float(Ome(1)); 
        Omevec[3*i+2] = float(Ome(2)); 
        Aomvec[3*i  ] = float(Ao(0));
        Aomvec[3*i+1] = float(Ao(1)); 
        Aomvec[3*i+2] = float(Ao(2)); 
        Amovec[3*i  ] = float(L(0));
        Amovec[3*i+1] = float(L(1)); 
        Amovec[3*i+2] = float(L(2)); 
        Ufovec[3*i  ] = (float) Particles[i]->F(0);
        Ufovec[3*i+1] = (float) Particles[i]->F(1); 
        Ufovec[3*i+2] = (float) Particles[i]->F(2); 
        Hfovec[3*i  ] = (float) Particles[i]->Flbm(0);
        Hfovec[3*i+1] = (float) Particles[i]->Flbm(1); 
        Hfovec[3*i+2] = (float) Particles[i]->Flbm(2); 
        Inertm[6*i  ] = (float) Inertiar(0,0);
        Inertm[6*i+1] = (float) Inertiar(0,1);
        Inertm[6*i+2] = (float) Inertiar(0,2);
        Inertm[6*i+3] = (float) Inertiar(1,1);
        Inertm[6*i+4] = (float) Inertiar(1,2);
        Inertm[6*i+5] = (float) Inertiar(2,2);
        Ekin  [i]     = float(Particles[i]->Ekin+Particles[i]->Erot);
        Tag   [i]     = int  (Particles[i]->Tag);  
    }

    hsize_t dims[1];
    dims[0] = 6*Particles.Size();
    String dsname;
    dsname.Printf("Inertia");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Inertm);
    dims[0] = 3*Particles.Size();
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("Rod");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Rod);
    dsname.Printf("PVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dsname.Printf("PAngVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Omevec);
    dsname.Printf("PAngacceleration");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Aomvec);
    dsname.Printf("PAngMomentum");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Amovec);
    dsname.Printf("PUForce");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ufovec);
    dsname.Printf("PHForce");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Hfovec);
    dims[0] = Particles.Size();
    dsname.Printf("Radius");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
    dsname.Printf("PEkin");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ekin);
    dsname.Printf("PTag");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tag   );


    delete [] Radius;
    delete [] Posvec;
    delete [] Rod;
    delete [] Velvec;
    delete [] Omevec;
    delete [] Aomvec;
    delete [] Amovec;
    delete [] Ufovec;
    delete [] Hfovec;
    delete [] Inertm;
    delete [] Ekin;
    delete [] Tag;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
    

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    if(N_Faces>0)
    {
    oss << "   <Grid Name=\"DEM_Faces\">\n";
    oss << "     <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << N_Faces << "\">\n";
    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << N_Faces << " 3\">\n";
    oss << "        " << fn.CStr() <<":/FaceCon \n";
    oss << "       </DataItem>\n";
    oss << "     </Topology>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << N_Verts << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Verts \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Cluster\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Cluster \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/AngVelocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << "   <Grid Name=\"DEM_Center\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"PRod\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Rod\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Ekin\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PEkin \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVel\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Angacc\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngacceleration\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngMom\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngMomentum\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Hforce\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PHForce\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Uforce\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PUForce\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Inertia\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Inertia\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";


    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::WriteFrac (char const * FileKey)
{

    // Counting the number of non valid cohesive interactons
    size_t N_Faces = 0;
    size_t N_Verts = 0;
    size_t nvbi = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (!BInteractons[i]->valid) 
        {
            Particle * P1 = BInteractons[i]->P1;
            Particle * P2 = BInteractons[i]->P2;
            Face     * F1 = P1->Faces[BInteractons[i]->IF1];
            Face     * F2 = P2->Faces[BInteractons[i]->IF2];
            nvbi++;
            N_Faces += F1->Edges.Size();
            N_Faces += F2->Edges.Size();
            N_Verts += F1->Edges.Size() + 1;
            N_Verts += F2->Edges.Size() + 1;
        }
    }

    //std::cout << "1 " << nvbi << std::endl;

    if (nvbi==0) return;

    //Geometric information
    float  * Verts   = new float [3*N_Verts];
    int    * FaceCon = new int   [3*N_Faces];

    //Atributes
    int    * Tags    = new int   [  N_Faces];
       
    size_t n_verts = 0;
    size_t n_faces = 0;
    size_t n_attrs = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (BInteractons[i]->valid) continue;
        Particle * P1 = BInteractons[i]->P1;
        Particle * P2 = BInteractons[i]->P2;
        size_t    IF1 = BInteractons[i]->IF1;
        size_t    IF2 = BInteractons[i]->IF2;
        //std::cout << P1 << " " << P2 <<std::endl;

        //For P1
        {
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(P1->Verts.Size());
            Array<Vec3_t> Vres (P1->Verts.Size());
            for (size_t j=0;j<P1->Verts.Size();j++)
            {
                Vtemp[j] = *P1->Verts[j];
                Vres [j] = *P1->Verts[j];
            }
            double multiplier = 0.0;
            if (P1->Eroded&&P1->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,P1->EdgeCon,P1->FaceCon,Vres,P1->Props.R);
                multiplier = 1.0;
            }
            for (size_t j=0;j<P1->FaceCon[IF1].Size();j++)
            {
                size_t k = P1->FaceCon[IF1][j];
                Verts[n_verts++] = float(Vres[k](0));
                Verts[n_verts++] = float(Vres[k](1));
                Verts[n_verts++] = float(Vres[k](2));
            }
            size_t n_reff = n_verts/3;
            Vec3_t C,N;
            P1->Faces[IF1]->Centroid(C);
            P1->Faces[IF1]->Normal(N);
            Verts[n_verts++] = float(C(0) + multiplier*P1->Props.R*N(0));
            Verts[n_verts++] = float(C(1) + multiplier*P1->Props.R*N(1));
            Verts[n_verts++] = float(C(2) + multiplier*P1->Props.R*N(2));
            for (size_t j=0;j<P1->FaceCon[IF1].Size();j++)
            {
                FaceCon[n_faces++] = int(n_reff);  
                FaceCon[n_faces++] = int(n_refv + j);
                FaceCon[n_faces++] = int(n_refv + (j+1)%(P1->FaceCon[IF1].Size()));
                Tags   [n_attrs]   = int(P1->Tag);
                n_attrs++;
            }
        }
        //std::cout << "2" << std::endl;
        //For P2
        {
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(P2->Verts.Size());
            Array<Vec3_t> Vres (P2->Verts.Size());
            for (size_t j=0;j<P2->Verts.Size();j++)
            {
                Vtemp[j] = *P2->Verts[j];
                Vres [j] = *P2->Verts[j];
            }
            //std::cout << "3" << std::endl;
            double multiplier = 0.0;
            if (P2->Eroded&&P2->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,P2->EdgeCon,P2->FaceCon,Vres,P2->Props.R);
                multiplier = 1.0;
            }
            //std::cout << "4" << std::endl;
            for (size_t j=0;j<P2->FaceCon[IF2].Size();j++)
            {
                size_t k = P2->FaceCon[IF2][j];
                Verts[n_verts++] = float(Vres[k](0));
                Verts[n_verts++] = float(Vres[k](1));
                Verts[n_verts++] = float(Vres[k](2));
            }
            //std::cout << "5" << std::endl;
            size_t n_reff = n_verts/3;
            Vec3_t C,N;
            P2->Faces[IF2]->Centroid(C);
            P2->Faces[IF2]->Normal(N);
            Verts[n_verts++] = float(C(0) + multiplier*P2->Props.R*N(0));
            Verts[n_verts++] = float(C(1) + multiplier*P2->Props.R*N(1));
            Verts[n_verts++] = float(C(2) + multiplier*P2->Props.R*N(2));
            //std::cout << "6" << std::endl;
            for (size_t j=0;j<P2->FaceCon[IF2].Size();j++)
            {
                FaceCon[n_faces++] = int(n_reff);  
                FaceCon[n_faces++] = int(n_refv + j);
                FaceCon[n_faces++] = int(n_refv + (j+1)%(P2->FaceCon[IF2].Size()));
                Tags   [n_attrs]   = int(P2->Tag);
                n_attrs++;
            }
            //std::cout << "7" << std::endl;
        }
    }
    //std::cout << n_faces << " " << N_Faces << std::endl;
    //Write the data
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[1];
    String dsname;
    dims[0] = 3*N_Verts;
    dsname.Printf("Verts");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Verts);
    dims[0] = 3*N_Faces;
    dsname.Printf("FaceCon");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,FaceCon);
    dims[0] = N_Faces;
    dsname.Printf("Tag");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tags   );

    //Erasing the data
    delete [] Verts;
    delete [] FaceCon;
    delete [] Tags;

    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"DEM_Faces\">\n";
    oss << "     <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << N_Faces << "\">\n";
    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << N_Faces << " 3\">\n";
    oss << "        " << fn.CStr() <<":/FaceCon \n";
    oss << "       </DataItem>\n";
    oss << "     </Topology>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << N_Verts << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Verts \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";

    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

inline void Domain::Save (char const * FileKey)
{

    // Opening the file for writing
    String fn(FileKey);
    fn.append(".hdf5");
    //if (Util::FileExists(fn))
    //{
        //String command;
        //command.Printf("rm %s",fn.CStr());
        //system(command.CStr());
    //}
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Storing the number of particles in the domain
    int data[1];
    data[0]=Particles.Size();
    hsize_t dims[1];
    dims[0]=1;
    H5LTmake_dataset_int(file_id,"/NP",1,dims,data);

    // Storing information about periodic BCs
    double perdat[1];
    perdat[0]=Xmax;
    H5LTmake_dataset_double(file_id,"/Xmax",1,dims,perdat);
    perdat[0]=Xmin;
    H5LTmake_dataset_double(file_id,"/Xmin",1,dims,perdat);
    perdat[0]=Ymax;
    H5LTmake_dataset_double(file_id,"/Ymax",1,dims,perdat);
    perdat[0]=Ymin;
    H5LTmake_dataset_double(file_id,"/Ymin",1,dims,perdat);
    perdat[0]=Zmax;
    H5LTmake_dataset_double(file_id,"/Zmax",1,dims,perdat);
    perdat[0]=Zmin;
    H5LTmake_dataset_double(file_id,"/Zmin",1,dims,perdat);
    


    for (size_t i=0; i<Particles.Size(); i++)
    {
        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gcreate(file_id, par.CStr(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


        // Storing some scalar variables
        double dat[1];
        dat[0] = Particles[i]->Props.R;
        H5LTmake_dataset_double(group_id,"SR",1,dims,dat);
        dat[0] = Particles[i]->Props.rho;
        H5LTmake_dataset_double(group_id,"Rho",1,dims,dat);
        dat[0] = Particles[i]->Props.m;
        H5LTmake_dataset_double(group_id,"m",1,dims,dat);
        dat[0] = Particles[i]->Props.V;
        H5LTmake_dataset_double(group_id,"V",1,dims,dat);
        dat[0] = Particles[i]->Diam;
        H5LTmake_dataset_double(group_id,"Diam",1,dims,dat);
        dat[0] = Particles[i]->Dmax;
        H5LTmake_dataset_double(group_id,"Dmax",1,dims,dat);
        int datint[1];
        datint[0] = Particles[i]->Index;
        H5LTmake_dataset_int(group_id,"Index",1,dims,datint);


        int tag[1];
        tag[0] = Particles[i]->Tag;
        H5LTmake_dataset_int(group_id,"Tag",1,dims,tag);

        // Storing vectorial variables
        double cd[3];
        hsize_t dd[1];
        dd[0] = 3;

        cd[0]=Particles[i]->x(0);
        cd[1]=Particles[i]->x(1);
        cd[2]=Particles[i]->x(2);
        H5LTmake_dataset_double(group_id,"x",1,dd,cd);

        cd[0]=Particles[i]->xb(0);
        cd[1]=Particles[i]->xb(1);
        cd[2]=Particles[i]->xb(2);
        H5LTmake_dataset_double(group_id,"xb",1,dd,cd);

        cd[0]=Particles[i]->v(0);
        cd[1]=Particles[i]->v(1);
        cd[2]=Particles[i]->v(2);
        H5LTmake_dataset_double(group_id,"v",1,dd,cd);

        cd[0]=Particles[i]->w(0);
        cd[1]=Particles[i]->w(1);
        cd[2]=Particles[i]->w(2);
        H5LTmake_dataset_double(group_id,"w",1,dd,cd);

        cd[0]=Particles[i]->wb(0);
        cd[1]=Particles[i]->wb(1);
        cd[2]=Particles[i]->wb(2);
        H5LTmake_dataset_double(group_id,"wb",1,dd,cd);

        cd[0]=Particles[i]->I(0);
        cd[1]=Particles[i]->I(1);
        cd[2]=Particles[i]->I(2);
        H5LTmake_dataset_double(group_id,"I",1,dd,cd);

        double cq[4];
        dd[0] = 4;
        cq[0]=Particles[i]->Q(0);
        cq[1]=Particles[i]->Q(1);
        cq[2]=Particles[i]->Q(2);
        cq[3]=Particles[i]->Q(3);
        H5LTmake_dataset_double(group_id,"Q",1,dd,cq);




        // Storing the number of vertices of each particle
        data[0] = Particles[i]->Verts.Size();
        H5LTmake_dataset_int(group_id,"n_vertices",1,dims,data);
        hid_t gv_id;
        gv_id = H5Gcreate(group_id,"Verts", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Storing each vertex 
        for (size_t j=0;j<Particles[i]->Verts.Size();j++)
        {
            String parv;
            parv.Printf("Verts_%08d",j);
            double cod[3];
            cod[0]=(*Particles[i]->Verts[j])(0);
            cod[1]=(*Particles[i]->Verts[j])(1);
            cod[2]=(*Particles[i]->Verts[j])(2);
            hsize_t dim[1];
            dim[0]=3;
            H5LTmake_dataset_double(gv_id,parv.CStr(),1,dim,cod);
        }

        // Number of edges of the particle
        data[0] = Particles[i]->Edges.Size();
        H5LTmake_dataset_int(group_id,"n_edges",1,dims,data);
        gv_id = H5Gcreate(group_id,"Edges", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Edges
        for (size_t j=0;j<Particles[i]->Edges.Size();j++)
        {
            String parv;
            parv.Printf("Edges_%08d",j);
            int co[2];
            co[0] = Particles[i]->EdgeCon[j][0];
            co[1] = Particles[i]->EdgeCon[j][1];
            hsize_t dim[1];
            dim[0] =2;
            H5LTmake_dataset_int(gv_id,parv.CStr(),1,dim,co);
        }
        
        // Number of faces of the particle
        data[0] = Particles[i]->Faces.Size();
        H5LTmake_dataset_int(group_id,"n_faces",1,dims,data);
        gv_id = H5Gcreate(group_id,"Faces", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Number of cylinders of the particle
        data[0] = Particles[i]->Cylinders.Size();
        H5LTmake_dataset_int(group_id,"n_cylinders",1,dims,data);
        
        // Faces
        for (size_t j=0;j<Particles[i]->Faces.Size();j++)
        {
            String parv;
            parv.Printf("Faces_%08d",j);
            int co[Particles[i]->FaceCon[j].Size()];
            hsize_t dim[1];
            dim[0]= Particles[i]->FaceCon[j].Size();
            for (size_t k=0;k<Particles[i]->FaceCon[j].Size();k++)
            {
                co[k]=Particles[i]->FaceCon[j][k];
            }
            H5LTmake_dataset_int(gv_id,parv.CStr(),1,dim,co);
        }
        
    }

    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
    //sleep(5);
}

inline void Domain::Load (char const * FileKey)
{

    // Opening the file for reading
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Number of particles in the domain
    int data[1];
    H5LTread_dataset_int(file_id,"/NP",data);
    size_t NP = data[0];

    // Reading information about periodic BCs
    double perdat[1];
    H5LTread_dataset_double(file_id,"/Xmax",perdat);
    Xmax=perdat[0];
    H5LTread_dataset_double(file_id,"/Xmin",perdat);
    Xmin=perdat[0];
    H5LTread_dataset_double(file_id,"/Ymax",perdat);
    Ymax=perdat[0];
    H5LTread_dataset_double(file_id,"/Ymin",perdat);
    Ymin=perdat[0];
    H5LTread_dataset_double(file_id,"/Zmax",perdat);
    Zmax=perdat[0];
    H5LTread_dataset_double(file_id,"/Zmin",perdat);
    Zmin=perdat[0];

    // Loading the particles
    for (size_t i=0; i<NP; i++)
    {

        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gopen(file_id, par.CStr(),H5P_DEFAULT);

        // Finding the particle's position for the domain decomposition
        double X[3];
        H5LTread_dataset_double(group_id,"x",X);


        // Loading the Vertices
        H5LTread_dataset_int(group_id,"n_vertices",data);
        size_t nv = data[0];
        hid_t gv_id;
        gv_id = H5Gopen(group_id,"Verts", H5P_DEFAULT);
        Array<Vec3_t> V;

        for (size_t j=0;j<nv;j++)
        {
            String parv;
            parv.Printf("Verts_%08d",j);
            double cod[3];
            H5LTread_dataset_double(gv_id,parv.CStr(),cod);
            V.Push(Vec3_t(cod[0],cod[1],cod[2]));
        }
        
        // Loading the edges
        H5LTread_dataset_int(group_id,"n_edges",data);
        size_t ne = data[0];
        gv_id = H5Gopen(group_id,"Edges", H5P_DEFAULT);
        Array<Array <int> > E;

        for (size_t j=0;j<ne;j++)
        {
            String parv;
            parv.Printf("Edges_%08d",j);
            int cod[2];
            H5LTread_dataset_int(gv_id,parv.CStr(),cod);
            Array<int> Ep(2);
            Ep[0]=cod[0];
            Ep[1]=cod[1];
            E.Push(Ep);
        }

        // Loading the faces

        // Number of faces of the particle
        H5LTread_dataset_int(group_id,"n_faces",data);
        size_t nf = data[0];
        gv_id = H5Gopen(group_id,"Faces", H5P_DEFAULT);
        Array<Array <int> > F;
        
        // Faces
        for (size_t j=0;j<nf;j++)
        {
            String parv;
            parv.Printf("Faces_%08d",j);
            hsize_t dim[1];
            H5LTget_dataset_info(gv_id,parv.CStr(),dim,NULL,NULL);
            size_t ns = (size_t)dim[0];
            int co[ns];
            Array<int> Fp(ns);

            H5LTread_dataset_int(gv_id,parv.CStr(),co);
            
            for (size_t k=0;k<ns;k++)
            {
                Fp[k] = co[k];
            }

            F.Push(Fp);

        }

        // Number of cylinders
        H5LTread_dataset_int(group_id,"n_cylinders",data);
        size_t nc = data[0];

        Particles.Push (new Particle(-1,V,E,F,OrthoSys::O,OrthoSys::O,0.1,1.0));

        // Loading cylinder data if applicable
        if (nc>0)
        {   
            Vec3_t X0 = 0.5*(*Particles[Particles.Size()-1]->Verts[0] + *Particles[Particles.Size()-1]->Verts[2]);
            Vec3_t X1 = 0.5*(*Particles[Particles.Size()-1]->Verts[3] + *Particles[Particles.Size()-1]->Verts[5]);
            Particles[Particles.Size()-1]->Tori.Push     (new Torus(&X0,Particles[Particles.Size()-1]->Verts[0],Particles[Particles.Size()-1]->Verts[1]));
            Particles[Particles.Size()-1]->Tori.Push     (new Torus(&X1,Particles[Particles.Size()-1]->Verts[3],Particles[Particles.Size()-1]->Verts[4]));
            Particles[Particles.Size()-1]->Cylinders.Push(new Cylinder(Particles[Particles.Size()-1]->Tori[0],Particles[Particles.Size()-1]->Tori[1],Particles[Particles.Size()-1]->Verts[2],Particles[Particles.Size()-1]->Verts[5]));
        }

        // Loading vectorial variables
        Particles[Particles.Size()-1]->x = Vec3_t(X[0],X[1],X[2]);
        double cd[3];
        H5LTread_dataset_double(group_id,"xb",cd);
        Particles[Particles.Size()-1]->xb = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"v",cd);
        Particles[Particles.Size()-1]->v = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"w",cd);
        Particles[Particles.Size()-1]->w = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"wb",cd);
        Particles[Particles.Size()-1]->wb = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"I",cd);
        Particles[Particles.Size()-1]->I = Vec3_t(cd[0],cd[1],cd[2]);

        double cq[4];
        H5LTread_dataset_double(group_id,"Q",cq);
        Particles[Particles.Size()-1]->Q = Quaternion_t(cq[0],cq[1],cq[2],cq[3]);
    
        // Loading the scalar quantities of the particle
        double dat[1];
        H5LTread_dataset_double(group_id,"SR",dat);
        Particles[Particles.Size()-1]->Props.R = dat[0];
        H5LTread_dataset_double(group_id,"Rho",dat);
        Particles[Particles.Size()-1]->Props.rho = dat[0];
        H5LTread_dataset_double(group_id,"m",dat);
        Particles[Particles.Size()-1]->Props.m = dat[0];
        H5LTread_dataset_double(group_id,"V",dat);
        Particles[Particles.Size()-1]->Props.V = dat[0];
        H5LTread_dataset_double(group_id,"Diam",dat);
        Particles[Particles.Size()-1]->Diam = dat[0];
        H5LTread_dataset_double(group_id,"Dmax",dat);
        Particles[Particles.Size()-1]->Dmax = dat[0];
        int datint[1];
        H5LTread_dataset_int(group_id,"Index",datint);
        //Particles[Particles.Size()-1]->Index = datint[0];
        Particles[Particles.Size()-1]->Index = Particles.Size()-1;
        int tag[1];
        H5LTread_dataset_int(group_id,"Tag",tag);
        Particles[Particles.Size()-1]->Tag = tag[0];
        Particles[Particles.Size()-1]->PropsReady = true;

    }


    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);
}

inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    if (Particles.Size()==0) throw new Fatal("DEM::Domain::BoundingBox: There are no particles to build the bounding box");
    //minX = Vec3_t(Particles[0]->MinX(), Particles[0]->MinY(), Particles[0]->MinZ());
    //maxX = Vec3_t(Particles[0]->MaxX(), Particles[0]->MaxY(), Particles[0]->MaxZ());
    minX = Vec3_t( 1.6e308, 1.6e308, 1.6e308);
    maxX = Vec3_t(-1.6e308,-1.6e308,-1.6e308);
    for (size_t i=0; i<Particles.Size(); i++)
    {
        if (minX(0)>Particles[i]->MinX()&&Particles[i]->IsFree()) minX(0) = Particles[i]->MinX();
        if (minX(1)>Particles[i]->MinY()&&Particles[i]->IsFree()) minX(1) = Particles[i]->MinY();
        if (minX(2)>Particles[i]->MinZ()&&Particles[i]->IsFree()) minX(2) = Particles[i]->MinZ();
        if (maxX(0)<Particles[i]->MaxX()&&Particles[i]->IsFree()) maxX(0) = Particles[i]->MaxX();
        if (maxX(1)<Particles[i]->MaxY()&&Particles[i]->IsFree()) maxX(1) = Particles[i]->MaxY();
        if (maxX(2)<Particles[i]->MaxZ()&&Particles[i]->IsFree()) maxX(2) = Particles[i]->MaxZ();
    }
    if (Xmax-Xmin>1.0e-12)
    {
        minX(0) = Xmin;
        maxX(0) = Xmax;
    }
    if (Ymax-Ymin>1.0e-12)
    {
        minX(1) = Ymin;
        maxX(1) = Ymax;
    }
    if (Zmax-Zmin>1.0e-12)
    {
        minX(2) = Zmin;
        maxX(2) = Zmax;
    }
}

inline void Domain::BoundingBoxTag(Vec3_t & minX, Vec3_t & maxX, int Tag)
{
    if (Particles.Size()==0) throw new Fatal("DEM::Domain::BoundingBox: There are no particles to build the bounding box");
    //minX = Vec3_t(Particles[0]->MinX(), Particles[0]->MinY(), Particles[0]->MinZ());
    //maxX = Vec3_t(Particles[0]->MaxX(), Particles[0]->MaxY(), Particles[0]->MaxZ());
    minX = Vec3_t( 1.6e308, 1.6e308, 1.6e308);
    maxX = Vec3_t(-1.6e308,-1.6e308,-1.6e308);
    for (size_t i=0; i<Particles.Size(); i++)
    {
        if (minX(0)>Particles[i]->MinX()&&Particles[i]->Tag==Tag) minX(0) = Particles[i]->MinX();
        if (minX(1)>Particles[i]->MinY()&&Particles[i]->Tag==Tag) minX(1) = Particles[i]->MinY();
        if (minX(2)>Particles[i]->MinZ()&&Particles[i]->Tag==Tag) minX(2) = Particles[i]->MinZ();
        if (maxX(0)<Particles[i]->MaxX()&&Particles[i]->Tag==Tag) maxX(0) = Particles[i]->MaxX();
        if (maxX(1)<Particles[i]->MaxY()&&Particles[i]->Tag==Tag) maxX(1) = Particles[i]->MaxY();
        if (maxX(2)<Particles[i]->MaxZ()&&Particles[i]->Tag==Tag) maxX(2) = Particles[i]->MaxZ();
    }
    if (Xmax-Xmin>1.0e-12)
    {
        minX(0) = Xmin;
        maxX(0) = Xmax;
    }
    if (Ymax-Ymin>1.0e-12)
    {
        minX(1) = Ymin;
        maxX(1) = Ymax;
    }
    if (Zmax-Zmin>1.0e-12)
    {
        minX(2) = Zmin;
        maxX(2) = Zmax;
    }
}

inline void Domain::BoundingBoxAll(Vec3_t & minX, Vec3_t & maxX)
{
    if (Particles.Size()==0) throw new Fatal("DEM::Domain::BoundingBox: There are no particles to build the bounding box");
    minX = Vec3_t(Particles[0]->MinX(), Particles[0]->MinY(), Particles[0]->MinZ());
    maxX = Vec3_t(Particles[0]->MaxX(), Particles[0]->MaxY(), Particles[0]->MaxZ());
    for (size_t i=1; i<Particles.Size(); i++)
    {
        if (minX(0)>Particles[i]->MinX()) minX(0) = Particles[i]->MinX();
        if (minX(1)>Particles[i]->MinY()) minX(1) = Particles[i]->MinY();
        if (minX(2)>Particles[i]->MinZ()) minX(2) = Particles[i]->MinZ();
        if (maxX(0)<Particles[i]->MaxX()) maxX(0) = Particles[i]->MaxX();
        if (maxX(1)<Particles[i]->MaxY()) maxX(1) = Particles[i]->MaxY();
        if (maxX(2)<Particles[i]->MaxZ()) maxX(2) = Particles[i]->MaxZ();
    }
    if (Xmax-Xmin>1.0e-12)
    {
        minX(0) = Xmin;
        maxX(0) = Xmax;
    }
    if (Ymax-Ymin>1.0e-12)
    {
        minX(1) = Ymin;
        maxX(1) = Ymax;
    }
    if (Zmax-Zmin>1.0e-12)
    {
        minX(2) = Zmin;
        maxX(2) = Zmax;
    }
}

inline void Domain::Center(Vec3_t C)
{
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    Vec3_t Transport(-0.5*(maxX+minX));
    Transport += C;
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Translate(Transport);
}

inline void Domain::ClearInteractons()
{
    // delete old interactors
    for (size_t i=0; i<CInteractons.Size(); ++i)
    {
        if (CInteractons[i]!=NULL) delete CInteractons[i];
    }
    CInteractons.Resize(0);
    for (size_t i=0; i<BInteractons.Size(); ++i)
    {
        if (BInteractons[i]!=NULL) delete BInteractons[i];
    }
    BInteractons.Resize(0);
    Interactons.Resize(0);

    PairtoCInt.clear();
}

inline void Domain::ResetInteractons()
{
    // delete old interactors
    for (size_t i=0; i<CInteractons.Size(); ++i)
    {
        if (CInteractons[i]!=NULL) delete CInteractons[i];
    }

    // new interactors
    CInteractons.Resize(0);
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            bool pj_has_vf = !Particles[j]->IsFree();


            // if both particles have any component specified or they are far away, don't create any intereactor
            bool close = (Distance(Particles[i]->x,Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close ) continue;
            size_t hash = HashFunction(i,j);
            //Listofpairs.insert(hash);
            PairtoCInt[hash] = CInteractons.Size();

            // if both particles are spheres (just one vertex)
            if (Particles[i]->Verts.Size()==1 && Particles[j]->Verts.Size()==1)
            {
                CInteractons.Push (new CInteractonSphere(Particles[i],Particles[j]));
            }

            // normal particles
            else
            {
                CInteractons.Push (new CInteracton(Particles[i],Particles[j]));
            }

            std::pair<int,int> p (Particles[i]->Tag,Particles[j]->Tag);
            if (FricCoeff.count(p)==1)
            {
                CInteractons[CInteractons.Size()-1]->Mu = FricCoeff[p];                    
            }
        }
    }
}

inline void Domain::ResetDisplacements()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].Dmx = 0.0;
        MTD[i].LLC.Resize(0);
    }
    bool px = (Xmax-Xmin)>Alpha;
    bool py = (Ymax-Ymin)>Alpha;
    bool pz = (Zmax-Zmin)>Alpha;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Particles.Size();i++)
    {
        Particles[i]->ResetDisplacements();
        if(Particles[i]->IsFree())
        {
            //#ifndef USE_CUDA
            if (px&&Particles[i]->x(0)>=Xmax) Particles[i]->Translate(Vec3_t(Xmin-Xmax,0.0,0.0));
            if (px&&Particles[i]->x(0)< Xmin) Particles[i]->Translate(Vec3_t(Xmax-Xmin,0.0,0.0));
            if (py&&Particles[i]->x(1)>=Ymax) Particles[i]->Translate(Vec3_t(0.0,Ymin-Ymax,0.0));
            if (py&&Particles[i]->x(1)< Ymin) Particles[i]->Translate(Vec3_t(0.0,Ymax-Ymin,0.0));
            if (pz&&Particles[i]->x(2)>=Zmax) Particles[i]->Translate(Vec3_t(0.0,0.0,Zmin-Zmax));
            if (pz&&Particles[i]->x(2)< Zmin) Particles[i]->Translate(Vec3_t(0.0,0.0,Zmax-Zmin));
            //#endif
            iVec3_t idx;
            idx(0) = floor((Particles[i]->x(0)-LCxmin(0))/DxLC(0));
            idx(1) = floor((Particles[i]->x(1)-LCxmin(1))/DxLC(1));
            idx(2) = floor((Particles[i]->x(2)-LCxmin(2))/DxLC(2));
            MTD[omp_get_thread_num()].LLC.Push(std::make_pair(idx,i));
        }
    }
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LLC.Size();j++)
        {
            size_t idx = Pt2idx(MTD[i].LLC[j].first,LCellDim);
            LinkedCell[idx].Push(MTD[i].LLC[j].second);
        }
    }
}

inline void Domain::UpdateLinkedCells()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].LPP.Resize(0);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<FreePar.Size()  ;i++)
    for (size_t j=0;j<NoFreePar.Size();j++)
    {
        size_t i1 = std::min(FreePar[i],NoFreePar[j]);
        size_t i2 = std::max(FreePar[i],NoFreePar[j]);
        MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
    }
    bool px = (Xmax-Xmin)>Alpha;
    bool py = (Ymax-Ymin)>Alpha;
    bool pz = (Zmax-Zmin)>Alpha;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t idx=0;idx<LinkedCell.Size();idx++)
    {
        if (LinkedCell[idx].Size()==0) continue;
        iVec3_t index;
        idx2Pt(idx,index,LCellDim);
        for (size_t n=0  ;n<LinkedCell[idx].Size()-1;n++)
        {
            //std::cout << LinkedCell[idx][n] << " ";
            for (size_t m=n+1;m<LinkedCell[idx].Size()  ;m++)
            {
                size_t i1 = std::min(LinkedCell[idx][n],LinkedCell[idx][m]);
                size_t i2 = std::max(LinkedCell[idx][n],LinkedCell[idx][m]);
                MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
            }
        }
        //std::cout << std::endl;
        int i = index(0);
        int j = index(1);
        int k = index(2);
        std::set<size_t> UniqueLC;
        for (int knb=(pz ? k-1 : std::max(0,k-1));knb<=(pz ? k+1 : std::min(k+1,(int)LCellDim(2)));knb++)
        for (int jnb=(py ? j-1 : std::max(0,j-1));jnb<=(py ? j+1 : std::min(j+1,(int)LCellDim(1)));jnb++)
        for (int inb=(px ? i-1 : std::max(0,i-1));inb<=(px ? i+1 : std::min(i+1,(int)LCellDim(0)));inb++)
        {
            size_t ix = (inb+(int)LCellDim(0))%LCellDim(0);
            size_t iy = (jnb+(int)LCellDim(1))%LCellDim(1);
            size_t iz = (knb+(int)LCellDim(2))%LCellDim(2);
            iVec3_t Ptnb(ix,iy,iz);
            size_t idxnb = Pt2idx(Ptnb,LCellDim);
            if (UniqueLC.count(idxnb)!=0) continue;
            UniqueLC.insert(idxnb);
            if (idxnb>idx)
            {
                for (size_t n=0;n<LinkedCell[idx].Size()  ;n++)
                {
                    for (size_t m=0;m<LinkedCell[idxnb].Size()  ;m++)
                    {
                        size_t i1 = std::min(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        size_t i2 = std::max(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
                    }
                }
            }
        }
    }
    size_t Npp = 0;
    for (size_t i=0;i<Nproc;i++)
    {
        Npp += MTD[i].LPP.Size();
    }
    ListPosPairs.Resize(Npp);
    size_t idx = 0;
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LPP.Size();j++)
        {
            ListPosPairs[idx] = MTD[i].LPP[j];
            idx++;
        }
    }
}

inline double Domain::MaxDisplacement()
{
    double md = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        double mpd = Particles[i]->MaxDisplacement();
        if (mpd > md) md = mpd;
    }
    return md;
}

inline double Domain::MaxDim()
{
    double md = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {  
        if (!Particles[i]->IsFree()) continue;
        double mpd = Particles[i]->Dmax;
        if (mpd > md) md = mpd;
    }
    return md;
}

inline void Domain::ResetContacts()
{   
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].LC.Resize(0);
        MTD[i].LCI.Resize(0);
        MTD[i].CLCI.Resize(0);
        MTD[i].LCB.Resize(0);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<ListPosPairs.Size();n++)
    {
        size_t i = ListPosPairs[n].first;
        size_t j = ListPosPairs[n].second;
        bool pi_has_vf = !Particles[i]->IsFree();
        bool pj_has_vf = !Particles[j]->IsFree();
        bool close = (Distance(Particles[i]->x,Particles[j]->x,Per)<=Particles[i]->Dmax+Particles[j]->Dmax+2.0*Alpha);
        if ((pi_has_vf && pj_has_vf) || !close) continue;
        MTD[omp_get_thread_num()].CLCI.Push(std::make_pair(i,j));
        size_t hash = HashFunction(i,j);
        if (PairtoCInt.count(hash)!=0) continue;
        //std::set<size_t>::iterator it = Listofpairs.find(hash);
        //if (it != Listofpairs.end())
        //{
            //continue;
        //}
        MTD[omp_get_thread_num()].LC.Push(std::make_pair(i,j));
    }
    Array<std::pair<size_t,size_t> > CPairs;
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LC.Size();j++)
        {
            size_t n = MTD[i].LC[j].first;
            size_t m = MTD[i].LC[j].second;
            size_t hash = HashFunction(n,m);
            //Listofpairs.insert(hash);
            PairtoCInt[hash] = CInteractons.Size();
            if (Particles[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
            {
                if (!MostlySpheres) CInteractons.Push (new CInteractonSphere(Particles[n],Particles[m],ContactLaw));
            }
            else
            {
                CInteractons.Push (new CInteracton(Particles[n],Particles[m],ContactLaw));
            }
            std::pair<int,int> p (Particles[n]->Tag,Particles[m]->Tag);
            if (FricCoeff.count(p)==1)
            {
                CInteractons[CInteractons.Size()-1]->Mu = FricCoeff[p];                    
            }
        }
        for (size_t j=0;j<MTD[i].CLCI.Size();j++)
        {
            size_t n = MTD[i].CLCI[j].first;
            size_t m = MTD[i].CLCI[j].second;
            CPairs.Push(std::make_pair(n,m));
        }
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<CPairs.Size();n++)
    {
        size_t i = CPairs[n].first;
        size_t j = CPairs[n].second;
        size_t hash = HashFunction(i,j);
        size_t nc = PairtoCInt[hash];
        if(CInteractons[nc]->UpdateContacts(Alpha,Per,iter)) MTD[omp_get_thread_num()].LCI.Push(nc);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<BInteractons.Size();n++)
    {
        if(BInteractons[n]->UpdateContacts(Alpha,Per)) MTD[omp_get_thread_num()].LCB.Push(n);
    }
     
    Interactons.Resize(0);
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LCI.Size();j++)
        {
            Interactons.Push(CInteractons[MTD[i].LCI[j]]);
        }
        for (size_t j=0;j<MTD[i].LCB.Size();j++)
        {
            Interactons.Push(BInteractons[MTD[i].LCB[j]]);
        }
    }
}

inline void Domain::UpdateContacts()
{
    LinkedCell.Resize(0);

    if (FreePar.Size()!=0)
    {
        BoundingBox(LCxmin,LCxmax);
        LCellDim(0) = (LCxmax(0) - LCxmin(0))/(2.0*Beta*MaxDmax);
        LCellDim(1) = (LCxmax(1) - LCxmin(1))/(2.0*Beta*MaxDmax);
        LCellDim(2) = (LCxmax(2) - LCxmin(2))/(2.0*Beta*MaxDmax);
        if(LCellDim(0)==0) LCellDim(0) = 1;
        if(LCellDim(1)==0) LCellDim(1) = 1;
        if(LCellDim(2)==0) LCellDim(2) = 1;
        DxLC    (0) = (LCxmax(0) - LCxmin(0))/LCellDim(0);
        DxLC    (1) = (LCxmax(1) - LCxmin(1))/LCellDim(1);
        DxLC    (2) = (LCxmax(2) - LCxmin(2))/LCellDim(2);

        LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));
    }

    ResetDisplacements();

    UpdateLinkedCells();

    ResetContacts();
}

inline void Domain::CalcForceSphere()
{
    //std::cout << "Pairs size = " << ListPosPairs.Size() << std::endl;
    //std::set<std::pair<Particle *,Particle *> >::iterator it;
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
    //for (it=Listofpairs.begin();it!=Listofpairs.end();++it)
    for (size_t np=0;np<ListPosPairs.size();np++)
    {
        //std::cout << "1" << std::endl;
        size_t i = ListPosPairs[np].first;
        size_t j = ListPosPairs[np].second;
        //std::cout << i << " " << j << std::endl;
        DEM::Particle * P1 = Particles[i];
        DEM::Particle * P2 = Particles[j];
        //DEM::Particle * P1 = it->first;
        //DEM::Particle * P2 = it->second;
        //size_t i = P1->Index;
        //size_t j = P2->Index;
        bool pi = (P1->Verts.Size()==1);
        bool pj = (P2->Verts.Size()==1);
        //bool close = (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*Alpha);
        if (!(pi&&pj)) continue;
        //std::cout << i << " " << j << std::endl;

        Vec3_t xi = P1->x;
        Vec3_t xf = P2->x;
        double dist = norm(P1->x - P2->x);
        double delta = P1->Props.R + P2->Props.R - dist;
        if (delta>0)
        {
            //std::cout << "2" << std::endl;
            Vec3_t n = (xf-xi)/dist;
            Vec3_t x = xi+n*((P1->Props.R*P1->Props.R-P2->Props.R*P2->Props.R+dist*dist)/(2*dist));
            Vec3_t t1,t2,x1,x2;
            Rotation(P1->w,P1->Q,t1);
            Rotation(P2->w,P2->Q,t2);
            x1 = x - P1->x;
            x2 = x - P2->x;
            Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
            Vec3_t vt = vrel - dot(n,vrel)*n;

            double Kn = ReducedValue(P1->Props.Kn,P2->Props.Kn);
            double Kt = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
            double me = ReducedValue(P1->Props.m ,P2->Props.m );
            double Gn = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn);
            double Gt = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt);
            double beta = 2*ReducedValue(P1->Props.Beta,P2->Props.Beta);
            double eta  = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);
            //double Bn   = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn);
            //double Bt   = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt);
            //double eps  = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
            double Mu;

            if (P1->Props.Mu>1.0e-12&&P2->Props.Mu>1.0e-12)
            {
                Mu          = std::max(P1->Props.Mu,P2->Props.Mu);
            }
            else 
            {
                Mu          = 0.0;
            }

            if (Gn < -0.001)
            {
                if (fabs(Gn)>1.0) throw new Fatal("CInteractonSphere the restitution coefficient is greater than 1");
                Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
                Gt = 0.0;
            }
            Gn *= me;
            Gt *= me;

            Vec3_t Fn = Kn*delta*n;
            //std::pair<int,int> p;
            //p = std::make_pair(i,j);
            size_t p = HashFunction(i,j);
            if (FricSpheres.count(p)==0) 
            {             
#ifdef USE_OMP
                omp_set_lock  (&lck);
#endif
                FricSpheres[p] = OrthoSys::O;
#ifdef USE_OMP
                omp_unset_lock(&lck);
#endif
            }
            FricSpheres[p] += vt*Dt;
            FricSpheres[p] -= dot(FricSpheres[p],n)*n;
            Vec3_t tan = FricSpheres[p];
            if (norm(tan)>0.0) tan/=norm(tan);
            if (norm(FricSpheres[p])>Mu*norm(Fn)/Kt)
            {
                FricSpheres[p] = Mu*norm(Fn)/Kt*tan;
            }
            Vec3_t F = Fn + Kt*FricSpheres[p] + Gn*dot(n,vrel)*n + Gt*vt;
            //Vec3_t F = Fn + P1->Props.Gn*dot(n,vrel)*n + P1->Props.Gt*vt;
            Vec3_t F1   = -F;
            Vec3_t F2   =  F;

            Vec3_t T, Tt;
            Tt = cross (x1,F);
            Quaternion_t q;
            Conjugate (P1->Q,q);
            Rotation  (Tt,q,T);
            Vec3_t T1 = -T;

            Tt = cross (x2,F);
            Conjugate (P2->Q,q);
            Rotation  (Tt,q,T);
            Vec3_t T2 =  T;
            
            //std::cout << "3" << std::endl;


            //rolling resistance
            if (dot(F,n)<0) F-=dot(F,n)*n;
            Vec3_t Normal = Fn/norm(Fn);
            Vec3_t Vr = P1->Props.R*P2->Props.R*cross(Vec3_t(t1 - t2),Normal)/(P1->Props.R+P2->Props.R);
            if (RollSpheres.count(p)==0) 
            {
#ifdef USE_OMP
                omp_set_lock  (&lck);
#endif
                RollSpheres[p] = OrthoSys::O;
#ifdef USE_OMP
                omp_unset_lock(&lck);
#endif
            }
            RollSpheres[p] += Vr*Dt;
            RollSpheres[p] -= dot(RollSpheres[p],Normal)*Normal;
            tan = RollSpheres[p];
            if (norm(tan)>0.0) tan/=norm(tan);
            double Kr = beta*Kt;
            if (norm(RollSpheres[p])>eta*Mu*norm(Fn)/Kr)
            {
                RollSpheres[p] = eta*Mu*norm(Fn)/Kr*tan;
            }
            Vec3_t Ft = -Kr*RollSpheres[p];

            //std::cout << "4" << std::endl;

            Tt = P1->Props.R*cross(Normal,Ft);
            Conjugate (P1->Q,q);
            Rotation  (Tt,q,T);
            T1 += T;
//
            Tt = P2->Props.R*cross(Normal,Ft);
            Conjugate (P2->Q,q);
            Rotation  (Tt,q,T);
            T2 -= T;

#ifdef USE_OMP
            omp_set_lock  (&P1->lck);
#endif
            P1->F += F1;
            P1->T += T1;
#ifdef USE_OMP
            omp_unset_lock(&P1->lck);
            omp_set_lock  (&P2->lck);
#endif
            P2->F += F2;
            P2->T += T2;
#ifdef USE_OMP
            omp_unset_lock(&P2->lck);
#endif
        }
    }
    
}

// Auxiliar methods

inline void Domain::LinearMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        L += Particles[i]->Props.m*Particles[i]->v;
    }
}

inline void Domain::AngularMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Vec3_t t1,t2;
        t1 = Particles[i]->I(0)*Particles[i]->w(0),Particles[i]->I(1)*Particles[i]->w(1),Particles[i]->I(2)*Particles[i]->w(2);
        Rotation (t1,Particles[i]->Q,t2);
        L += Particles[i]->Props.m*cross(Particles[i]->x,Particles[i]->v)+t2;
    }
}

inline double Domain::CalcEnergy (double & Ekin, double & Epot)
{
    // kinematic energy
    Ekin = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Ekin += 0.5*Particles[i]->Props.m*dot(Particles[i]->v,Particles[i]->v)
                + 0.5*(Particles[i]->I(0)*Particles[i]->w(0)*Particles[i]->w(0)
                      +Particles[i]->I(1)*Particles[i]->w(1)*Particles[i]->w(1)
                      +Particles[i]->I(2)*Particles[i]->w(2)*Particles[i]->w(2));
    }

    // potential energy
    Epot = 0.0;
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        Epot += CInteractons[i]->Epot;
    }

    // total energy
    return Ekin + Epot;
}

inline double Domain::CriticalDt ()
{
    if (ContactLaw==0)
    {
        double MaxKn   =  0.0;
        double MaxBn   =  0.0;
        double MinMass = -1.0;
        for (size_t i=0; i<Particles.Size(); i++) 
        { 
            if (Particles[i]->IsFree())
            {
                if (Particles[i]->Props.Kn > MaxKn  ) MaxKn   = Particles[i]->Props.Kn;
                if (Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = Particles[i]->Props.m;
            }
        }
        for (size_t i=0; i<BInteractons.Size(); i++)
        {
            double pbn = std::max(BInteractons[i]->Bn/BInteractons[i]->L0,BInteractons[i]->Bt/BInteractons[i]->L0);
            if (pbn > MaxBn) MaxBn = pbn;
        }

        return 0.1*sqrt(MinMass/(MaxKn+MaxBn));
    }
    else if (ContactLaw==1)
    {
        double MaxKn   =  0.0;
        double MinMass = -1.0;
        for (size_t i=0; i<Particles.Size(); i++) 
        { 
            if (Particles[i]->IsFree())
            {
                if (Particles[i]->Props.Kn > MaxKn  ) MaxKn   = Particles[i]->Props.Kn;
                if (Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = Particles[i]->Props.m/Particles[i]->Dmax;
            }
        }

        return 0.1*sqrt(MinMass/MaxKn);
    }
}

inline void Domain::EnergyOutput (size_t IdxOut, std::ostream & OF)
{
    // header
    if (IdxOut==0)
    {
        OF << Util::_10_6 << "Time" << Util::_8s << "Ekin" << Util::_8s << "Epot" << Util::_8s << "Evis" << Util::_8s << "Efric" << Util::_8s << "Wext" << std::endl;
    }
    double Ekin,Epot;
    CalcEnergy(Ekin,Epot);
    OF << Util::_10_6 << Time << Util::_8s << Ekin << Util::_8s << Epot << Util::_8s << Evis << Util::_8s << Efric << Util::_8s << Wext << std::endl;
}

inline void Domain::GetGSD (Array<double> & X, Array<double> & Y, Array<double> & D, size_t NDiv) const
{
    // calc GSD information
    Array<double> Vg;
    double Vs = 0.0;

    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particle * P = Particles[i];
        double Diam = sqrt((P->MaxX()-P->MinX())*(P->MaxX()-P->MinX())+(P->MaxY()-P->MinY())*(P->MaxY()-P->MinY())+(P->MaxZ()-P->MinZ())*(P->MaxX()-P->MinX()));
        Vs += Particles[i]->Props.V;
        Vg.Push(Particles[i]->Props.V);
        D.Push(Diam);
    }
    double Dmin  = D[D.TheMin()]; // minimum diameter
    double Dmax  = D[D.TheMax()]; // maximum diameter
    double Dspan = (Dmax-Dmin)/NDiv;
    for (size_t i=0; i<=NDiv; i++)
    {
        X.Push (i*Dspan+Dmin);
        double cumsum = 0;
        for (size_t j=0; j<D.Size(); j++)
        {
            if (D[j]<=i*Dspan+Dmin) cumsum++;
        }
        Y.Push (cumsum/Particles.Size());
    }
}

inline void Domain::Clusters ()
{
    Array<int> connections;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (BInteractons[i]->valid)
        {
            connections.Push(BInteractons[i]->P1->Index);
            connections.Push(BInteractons[i]->P2->Index);
        }
    }

    Util::Tree tree(connections);
    tree.GetClusters(Listofclusters);
    for (size_t i=0;i<Listofclusters.Size();i++)
    {
        for (size_t j=0;j<Listofclusters[i].Size();j++)
        {
            Particles[Listofclusters[i][j]]->Cluster = i;
        }
    }
    //std::cout << Listofclusters.Size() << std::endl;
    //std::cout << BInteractons.Size() << std::endl;
}

inline void Domain::DelParticles (Array<int> const & Tags)
{
    Array<int> idxs; // indices to be deleted
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        for (size_t j=0; j<Tags.Size(); ++j)
        {
            if (Particles[i]->Tag==Tags[j]) idxs.Push(i);
        }
    }
    //if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particle to be deleted");
    if (idxs.Size()<1) std::cout<< "Warning: Domain::DelParticles: Could not find any particle to be deleted" << std::endl;
    Particles.DelItems (idxs);
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        Particles[i]->Index = i;
    }
}

inline void Domain::DelParticlesIdx (Array<int> const & idxs)
{
    //if (idxs.Size()<1) throw new Fatal("Domain::DelParticlesIdx: no particles to be deleted");
    if (idxs.Size()<1) std::cout<< "Warning: Domain::DelParticles: Could not find any particle to be deleted" << std::endl;
    Particles.DelItems (idxs);
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        Particles[i]->Index = i;
    }
}

inline Particle * Domain::GetParticle (int Tag, bool Check)
{
    size_t idx   = 0;
    size_t count = 0;
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag)
        {
            if (!Check) return Particles[i];
            idx = i;
            count++;
        }
    }
    if      (count==0) throw new Fatal("Domain::GetParticle: Could not find Particle with Tag==%d",Tag);
    else if (count>1)  throw new Fatal("Domain::GetParticle: There are more than one particle with Tag==%d",Tag);
    return Particles[idx];
}

inline Particle const & Domain::GetParticle (int Tag, bool Check) const
{
    size_t idx   = 0;
    size_t count = 0;
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag)
        {
            if (!Check) return (*Particles[i]);
            idx = i;
            count++;
        }
    }
    if      (count==0) throw new Fatal("Domain::GetParticle: Could not find Particle with Tag==%d",Tag);
    else if (count>1)  throw new Fatal("Domain::GetParticle: There are more than one particle with Tag==%d",Tag);
    return (*Particles[idx]);
}

inline void Domain::GetParticles (int Tag, Array<Particle*> & P)
{
    P.Resize(0);
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag) P.Push(Particles[i]);
    }
    if (P.Size()==0) throw new Fatal("Domain::GetParticles: Could not find any Particle with Tag==%d",Tag);
}

#ifdef USE_CUDA
//////////////////////////////////////////// CUDA IMPLEMENTATION ///////////////////////////////////////////
inline void Domain::UpLoadDevice(size_t Nc, bool first)
{
    if (first)
    {
        demaux.Time = Time;
        demaux.dt   = Dt;
        demaux.Per.x  = Per(0);
        demaux.Per.y  = Per(1);
        demaux.Per.z  = Per(2);
        demaux.Xmin   = Xmin;
        demaux.Ymin   = Ymin;
        demaux.Zmin   = Zmin;
        demaux.Xmax   = Xmax;
        demaux.Ymax   = Ymax;
        demaux.Zmax   = Zmax;
        demaux.px     = (Xmax-Xmin)>Alpha;
        demaux.py     = (Ymax-Ymin)>Alpha;
        demaux.pz     = (Zmax-Zmin)>Alpha;
        demaux.nparts = Particles.Size();

        demaux.nverts = 0;
        demaux.nedges = 0;
        demaux.nfaces = 0;
        demaux.nfacid = 0;

        for (size_t ip=0;ip<Particles.Size();ip++)
        {
            demaux.nverts +=   Particles[ip]->Verts.Size();
            demaux.nedges += 2*Particles[ip]->Edges.Size();
            demaux.nfaces += 2*Particles[ip]->Faces.Size();
            for (size_t idf=0;idf<Particles[ip]->Faces.Size();idf++)
            {
                demaux.nfacid += Particles[ip]->Faces[idf]->Edges.Size();
            }
        }

        thrust::host_vector<ParticleCU>    hParticlesCU   (Particles.Size());
        thrust::host_vector<DynParticleCU> hDynParticlesCU(Particles.Size());
        thrust::host_vector<real3>         hVertsCU       (demaux.nverts   );
        thrust::host_vector<size_t>        hEdgesCU       (demaux.nedges   );
        thrust::host_vector<size_t>        hFacidCU       (demaux.nfacid   );
        thrust::host_vector<size_t>        hFacesCU       (demaux.nfaces   );

        size_t idv = 0;
        size_t ide = 0;
        size_t idf = 0;
        size_t idg = 0;
        for (size_t ip=0;ip<Particles.Size();ip++)
        {
            hParticlesCU[ip] .Nvi = idv;
            Particles   [ip]->Nvi = idv;
            for (size_t iv=0;iv<Particles[ip]->Verts.Size();iv++)
            {
                hVertsCU[idv].x = (*Particles[ip]->Verts[iv])(0);
                hVertsCU[idv].y = (*Particles[ip]->Verts[iv])(1);
                hVertsCU[idv].z = (*Particles[ip]->Verts[iv])(2);
                idv++;
            }
            hParticlesCU[ip].Nvf = idv;

            hParticlesCU[ip] .Nei = ide;
            Particles   [ip]->Nei = ide;
            for (size_t ie=0;ie<Particles[ip]->EdgeCon.Size();ie++)
            {
                hEdgesCU[2*ide  ] = Particles[ip]->EdgeCon[ie][0]+hParticlesCU[ip].Nvi;
                hEdgesCU[2*ide+1] = Particles[ip]->EdgeCon[ie][1]+hParticlesCU[ip].Nvi;
                ide++;
            }
            hParticlesCU[ip].Nef = ide;

            hParticlesCU[ip] .Nfi = idf;
            Particles   [ip]->Nfi = idf;
            for (size_t ic=0;ic<Particles[ip]->FaceCon.Size();ic++)
            {
                hFacesCU[2*idf  ] = idg;
                for (size_t ig=0;ig<Particles[ip]->FaceCon[ic].Size();ig++)
                {
                    hFacidCU[idg] = Particles[ip]->FaceCon[ic][ig]+hParticlesCU[ip].Nvi;
                    idg++;
                }
                hFacesCU[2*idf+1] = idg;
                idf++;                
            }
            hParticlesCU[ip].Nff = idf;
        }

        //std::cout << "2" << std::endl;
        #pragma omp parallel for schedule(static) num_threads(Nc)
        for (size_t ip=0;ip<Particles.Size();ip++)
        {
            UploadParticle(hDynParticlesCU[ip],hParticlesCU[ip],*Particles[ip]);
        }
        
        
        //std::cout << "3" << std::endl;
        bMaxDCU        .resize(demaux.nverts);
        bParticlesCU   .resize(hParticlesCU   .size());
        bDynParticlesCU.resize(hDynParticlesCU.size());
        bVertsCU       .resize(hVertsCU       .size());
        bVertsoCU      .resize(hVertsCU       .size());
        bEdgesCU       .resize(hEdgesCU       .size());
        bFacidCU       .resize(hFacidCU       .size());
        bFacesCU       .resize(hFacesCU       .size());
        thrust::copy(hParticlesCU   .begin(), hParticlesCU   .end(), bParticlesCU   .begin());
        thrust::copy(hDynParticlesCU.begin(), hDynParticlesCU.end(), bDynParticlesCU.begin());
        thrust::copy(hVertsCU       .begin(), hVertsCU       .end(), bVertsCU       .begin());
        thrust::copy(hVertsCU       .begin(), hVertsCU       .end(), bVertsoCU      .begin());
        thrust::copy(hEdgesCU       .begin(), hEdgesCU       .end(), bEdgesCU       .begin());
        thrust::copy(hFacidCU       .begin(), hFacidCU       .end(), bFacidCU       .begin());
        thrust::copy(hFacesCU       .begin(), hFacesCU       .end(), bFacesCU       .begin());
        pMaxDCU          = thrust::raw_pointer_cast(bMaxDCU        .data());
        pParticlesCU     = thrust::raw_pointer_cast(bParticlesCU   .data());
        pDynParticlesCU  = thrust::raw_pointer_cast(bDynParticlesCU.data());
        pVertsCU         = thrust::raw_pointer_cast(bVertsCU       .data());
        pVertsoCU        = thrust::raw_pointer_cast(bVertsoCU      .data());
        pEdgesCU         = thrust::raw_pointer_cast(bEdgesCU       .data());
        pFacidCU         = thrust::raw_pointer_cast(bFacidCU       .data());
        pFacesCU         = thrust::raw_pointer_cast(bFacesCU       .data());

    }
   
    demaux.nvvint = 0;
    demaux.neeint = 0;
    demaux.nvfint = 0;
    demaux.nfvint = 0;
    demaux.ncoint = Interactons.Size();

    size_t idvv = 0;
    size_t idee = 0;
    size_t idvf = 0;
    size_t idfv = 0;

    size_t * Ivv = new size_t[2*demaux.ncoint];
    size_t * Iee = new size_t[2*demaux.ncoint];
    size_t * Ivf = new size_t[2*demaux.ncoint];
    size_t * Ifv = new size_t[2*demaux.ncoint];

    for (size_t ii=0;ii<Interactons.Size();ii++)
    {
        size_t i1 = Interactons[ii]->I1;
        size_t i2 = Interactons[ii]->I2;
        size_t hash = HashFunction(i1,i2);
        DEM::CInteracton * Ci = CInteractons[PairtoCInt[hash]];
        Ivv[2*ii  ] = idvv;
        Iee[2*ii  ] = idee;
        Ivf[2*ii  ] = idvf;
        Ifv[2*ii  ] = idfv;
        if (Particles[i1]->Verts.Size()==1 && Particles[i2]->Verts.Size()==1)
        {
            demaux.nvvint++;
            idvv++;
        }
        else
        {
            demaux.neeint += Ci->Lee.Size();
            demaux.nvfint += Ci->Lvf.Size();
            demaux.nfvint += Ci->Lfv.Size();
            idee += Ci->Lee.Size();
            idvf += Ci->Lvf.Size();
            idfv += Ci->Lfv.Size();
        }
        Ivv[2*ii+1] = idvv;
        Iee[2*ii+1] = idee;
        Ivf[2*ii+1] = idvf;
        Ifv[2*ii+1] = idfv;
    }

    thrust::host_vector<InteractonCU>    hInteractons     (demaux.ncoint);
    thrust::host_vector<ComInteractonCU> hComInteractons  (demaux.ncoint);
    thrust::host_vector<DynInteractonCU> hDynInteractonsVV(demaux.nvvint);
    thrust::host_vector<DynInteractonCU> hDynInteractonsEE(demaux.neeint);
    thrust::host_vector<DynInteractonCU> hDynInteractonsVF(demaux.nvfint);
    thrust::host_vector<DynInteractonCU> hDynInteractonsFV(demaux.nfvint);

    //size_t idvv = 0;
    //size_t idee = 0;
    //size_t idvf = 0;
    //size_t idfv = 0;

    #pragma omp parallel for schedule(static) num_threads(Nc)
    for (size_t ii=0;ii<Interactons.Size();ii++)
    {
        size_t i1 = Interactons[ii]->I1;
        size_t i2 = Interactons[ii]->I2;
        size_t hash = HashFunction(i1,i2);
        real   r1   = Particles[i1]->Props.R;
        real   r2   = Particles[i2]->Props.R;
        DEM::CInteracton * Ci = CInteractons[PairtoCInt[hash]];
        InteractonCU     Icu;
        ComInteractonCU CIcu;
        Icu.BothFree= Ci->BothFree;
        Icu.Kn      = Ci->Kn;
        Icu.Kt      = Ci->Kt;
        Icu.Gn      = Ci->Gn;
        Icu.Gt      = Ci->Gt;
        Icu.Mu      = Ci->Mu;
        Icu.Fnf     = make_real3(0.0,0.0,0.0);
        Icu.Ftf     = make_real3(0.0,0.0,0.0);
        CIcu.I1     = i1;
        CIcu.I2     = i2;
        if (Particles[i1]->Verts.Size()==1 && Particles[i2]->Verts.Size()==1)
        {

                DEM::CInteractonSphere * Cis = static_cast<DEM::CInteractonSphere *>(Ci);
                
                Icu.Beta    = Cis->beta;
                Icu.Eta     = Cis->eta;

                DynInteractonCU DIcu;
                DIcu.Dmax1  = r1;
                DIcu.Dmax2  = r2;
                DIcu.Idx    = ii;
                Vec3_t Ft   = OrthoSys::O;
                if (norm(Cis->Fdvv)>0.0) Ft = Cis->Fdvv*Cis->Kt;
                DIcu.Ft     = make_real3(Ft(0),Ft(1),Ft(2));
                Ft          = OrthoSys::O;
                if (norm(Cis->Fdr) >0.0) Ft = Cis->Fdr*Cis->Kt*Cis->beta;
                DIcu.Fr     = make_real3(Ft(0),Ft(1),Ft(2));

                hDynInteractonsVV[Ivv[2*ii]] = DIcu;
                //hDynInteractonsVV[idvv] = DIcu;
                //idvv++;
        }
        else
        {
            for (size_t iee=0;iee<Ci->Lee.Size();iee++)
            {
                size_t if1  = Ci->Lee[iee].first ;
                size_t if2  = Ci->Lee[iee].second;
                DynInteractonCU DIcu;
                DIcu.Dmax1   = Particles[i1]->Edges[if1]->Dmax + r1;
                DIcu.Dmax2   = Particles[i2]->Edges[if2]->Dmax + r2;
                DIcu.Idx     = ii;
                DIcu.IF1     = if1 + Particles[i1]->Nei;
                DIcu.IF2     = if2 + Particles[i2]->Nei;

                size_t p    = HashFunction(if1,if2);
                Vec3_t Ft   = OrthoSys::O;
                if (Ci->Fdee.count(p)!=0) Ft = Ci->Fdee[p]*Ci->Kt;
                DIcu.Ft     = make_real3(Ft(0),Ft(1),Ft(2));

                hDynInteractonsEE[Iee[2*ii]+iee] = DIcu;
                //hDynInteractonsEE[idee] = DIcu;
                //idee ++;
            }
            for (size_t ivf=0;ivf<Ci->Lvf.Size();ivf++)
            {
                size_t if1  = Ci->Lvf[ivf].first ;
                size_t if2  = Ci->Lvf[ivf].second;
                DynInteractonCU DIcu;
                DIcu.Dmax1   = r1;
                DIcu.Dmax2   = Particles[i2]->Faces[if2]->Dmax + r2;
                DIcu.Idx     = ii;
                DIcu.IF1     = if1 + Particles[i1]->Nvi;
                DIcu.IF2     = if2 + Particles[i2]->Nfi;

                size_t p    = HashFunction(if1,if2);
                Vec3_t Ft   = OrthoSys::O;
                if (Ci->Fdvf.count(p)!=0) Ft = Ci->Fdvf[p]*Ci->Kt;
                DIcu.Ft     = make_real3(Ft(0),Ft(1),Ft(2));

                hDynInteractonsVF[Ivf[2*ii]+ivf] = DIcu;
                //hDynInteractonsVF[idvf] = DIcu;
                //idvf ++;
            }
            for (size_t ifv=0;ifv<Ci->Lfv.Size();ifv++)
            {
                size_t if1  = Ci->Lfv[ifv].first ;
                size_t if2  = Ci->Lfv[ifv].second;
                DynInteractonCU DIcu;
                DIcu.Dmax1   = Particles[i1]->Faces[if1]->Dmax + r1;
                DIcu.Dmax2   = r2;
                DIcu.Idx     = ii;
                DIcu.IF1     = if1 + Particles[i1]->Nfi;
                DIcu.IF2     = if2 + Particles[i2]->Nvi;

                size_t p    = HashFunction(if1,if2);
                Vec3_t Ft   = OrthoSys::O;
                if (Ci->Fdfv.count(p)!=0) Ft = Ci->Fdfv[p]*Ci->Kt;
                DIcu.Ft     = make_real3(Ft(0),Ft(1),Ft(2));

                hDynInteractonsFV[Ifv[2*ii]+ifv] = DIcu;
                //hDynInteractonsFV[idfv] = DIcu;
                //idfv ++;
            }
        }
        hInteractons     [ii  ] = Icu;
        hComInteractons  [ii  ] = CIcu;
    } 

    delete [] Ivv;
    delete [] Iee;
    delete [] Ivf;
    delete [] Ifv;
    
    bInteractons     .resize(hInteractons     .size());
    bComInteractons  .resize(hComInteractons  .size());
    bDynInteractonsVV.resize(hDynInteractonsVV.size());
    bDynInteractonsEE.resize(hDynInteractonsEE.size());
    bDynInteractonsVF.resize(hDynInteractonsVF.size());
    bDynInteractonsFV.resize(hDynInteractonsFV.size());

    thrust::copy(hInteractons     .begin(),hInteractons     .end(),bInteractons     .begin());
    thrust::copy(hComInteractons  .begin(),hComInteractons  .end(),bComInteractons  .begin());
    thrust::copy(hDynInteractonsVV.begin(),hDynInteractonsVV.end(),bDynInteractonsVV.begin());
    thrust::copy(hDynInteractonsEE.begin(),hDynInteractonsEE.end(),bDynInteractonsEE.begin());
    thrust::copy(hDynInteractonsVF.begin(),hDynInteractonsVF.end(),bDynInteractonsVF.begin());
    thrust::copy(hDynInteractonsFV.begin(),hDynInteractonsFV.end(),bDynInteractonsFV.begin());
   
    pInteractons      = thrust::raw_pointer_cast(bInteractons     .data());
    pComInteractons   = thrust::raw_pointer_cast(bComInteractons  .data());
    pDynInteractonsVV = thrust::raw_pointer_cast(bDynInteractonsVV.data());
    pDynInteractonsEE = thrust::raw_pointer_cast(bDynInteractonsEE.data());
    pDynInteractonsVF = thrust::raw_pointer_cast(bDynInteractonsVF.data());
    pDynInteractonsFV = thrust::raw_pointer_cast(bDynInteractonsFV.data());

    //if (!first) cudaFree  (pdemaux);
    if (iter!=0) cudaFree  (pdemaux);
    cudaMalloc(&pdemaux, sizeof(dem_aux));
    cudaMemcpy(pdemaux, &demaux, sizeof(dem_aux), cudaMemcpyHostToDevice);
    //std::cout << "4" << std::endl;
    //
}

inline void Domain::DnLoadDevice(size_t Nc, bool force)
{
    thrust::host_vector<DynParticleCU> hDynParticlesCU = bDynParticlesCU;
    thrust::host_vector<real3>         hVertsCU        = bVertsCU;
   
    #pragma omp parallel for schedule(static) num_threads(Nc)
    for (size_t ip=0;ip<Particles.Size();ip++)
    {
        DnloadParticle(hDynParticlesCU[ip],*Particles[ip]);
        for (size_t iv=0;iv<Particles[ip]->Verts.Size();iv++)
        {
            size_t idv = Particles[ip]->Nvi + iv;
            (*Particles[ip]->Verts[iv])(0) = hVertsCU[idv].x;
            (*Particles[ip]->Verts[iv])(1) = hVertsCU[idv].y;
            (*Particles[ip]->Verts[iv])(2) = hVertsCU[idv].z;
        }
        for (size_t ie=0;ie<Particles[ip]->Edges.Size();ie++)
        {
            Particles[ip]->Edges[ie]->UpdatedL();
        }
        for (size_t ic=0;ic<Particles[ip]->Faces.Size();ic++)
        {
            Particles[ip]->Faces[ic]->UpdatedL();
        }
    }

    if (force)
    {
        thrust::host_vector<ComInteractonCU> hComInteractons   = bComInteractons;
        thrust::host_vector<DynInteractonCU> hDynInteractonsVV = bDynInteractonsVV;
        thrust::host_vector<DynInteractonCU> hDynInteractonsEE = bDynInteractonsEE;
        thrust::host_vector<DynInteractonCU> hDynInteractonsVF = bDynInteractonsVF;
        thrust::host_vector<DynInteractonCU> hDynInteractonsFV = bDynInteractonsFV;

        #pragma omp parallel for schedule(static) num_threads(Nc)
        for (size_t ivv=0;ivv<hDynInteractonsVV.size();ivv++)
        {
            size_t id = hDynInteractonsVV[ivv].Idx;
            size_t i1 = hComInteractons  [id ].I1;
            size_t i2 = hComInteractons  [id ].I2;
            size_t hash = HashFunction(i1,i2);
            DEM::CInteracton * Ci = CInteractons[PairtoCInt[hash]];
            DEM::CInteractonSphere * Cis = static_cast<DEM::CInteractonSphere *>(Ci);
            Cis->Fdvv(0) = hDynInteractonsVV[ivv].Ft.x/Cis->Kt;
            Cis->Fdvv(1) = hDynInteractonsVV[ivv].Ft.y/Cis->Kt;
            Cis->Fdvv(2) = hDynInteractonsVV[ivv].Ft.z/Cis->Kt;
            Cis->Fdr (0) = hDynInteractonsVV[ivv].Fr.x/Cis->Kt;
            Cis->Fdr (1) = hDynInteractonsVV[ivv].Fr.y/Cis->Kt;
            Cis->Fdr (2) = hDynInteractonsVV[ivv].Fr.z/Cis->Kt;
            if (fabs(Cis->beta)>0.0)
            {
                Cis->Fdr/=Cis->beta;
            }
        }

        #pragma omp parallel for schedule(static) num_threads(Nc)
        for (size_t iee=0;iee<hDynInteractonsEE.size();iee++)
        {
            size_t id = hDynInteractonsEE[iee].Idx;
            size_t i1 = hComInteractons  [id ].I1;
            size_t i2 = hComInteractons  [id ].I2;
            size_t f1 = hDynInteractonsEE[iee].IF1-Particles[i1]->Nei;
            size_t f2 = hDynInteractonsEE[iee].IF2-Particles[i2]->Nei;
            size_t hash = HashFunction(i1,i2);
            DEM::CInteracton * Ci = CInteractons[PairtoCInt[hash]];
            hash = HashFunction(f1,f2);
            Ci->Fdee[hash](0) = hDynInteractonsEE[iee].Ft.x/Ci->Kt;
            Ci->Fdee[hash](1) = hDynInteractonsEE[iee].Ft.y/Ci->Kt;
            Ci->Fdee[hash](2) = hDynInteractonsEE[iee].Ft.z/Ci->Kt;
        }

        #pragma omp parallel for schedule(static) num_threads(Nc)
        for (size_t ivf=0;ivf<hDynInteractonsVF.size();ivf++)
        {
            size_t id = hDynInteractonsVF[ivf].Idx;
            size_t i1 = hComInteractons  [id ].I1;
            size_t i2 = hComInteractons  [id ].I2;
            size_t f1 = hDynInteractonsVF[ivf].IF1-Particles[i1]->Nvi;
            size_t f2 = hDynInteractonsVF[ivf].IF2-Particles[i2]->Nfi;
            size_t hash = HashFunction(i1,i2);
            DEM::CInteracton * Ci = CInteractons[PairtoCInt[hash]];
            hash = HashFunction(f1,f2);
            Ci->Fdvf[hash](0) = hDynInteractonsVF[ivf].Ft.x/Ci->Kt;
            Ci->Fdvf[hash](1) = hDynInteractonsVF[ivf].Ft.y/Ci->Kt;
            Ci->Fdvf[hash](2) = hDynInteractonsVF[ivf].Ft.z/Ci->Kt;
        }

        #pragma omp parallel for schedule(static) num_threads(Nc)
        for (size_t ifv=0;ifv<hDynInteractonsFV.size();ifv++)
        {
            size_t id = hDynInteractonsFV[ifv].Idx;
            size_t i1 = hComInteractons  [id ].I1;
            size_t i2 = hComInteractons  [id ].I2;
            size_t f1 = hDynInteractonsFV[ifv].IF1-Particles[i1]->Nfi;
            size_t f2 = hDynInteractonsFV[ifv].IF2-Particles[i2]->Nvi;
            size_t hash = HashFunction(i1,i2);
            DEM::CInteracton * Ci = CInteractons[PairtoCInt[hash]];
            hash = HashFunction(f1,f2);
            Ci->Fdfv[hash](0) = hDynInteractonsFV[ifv].Ft.x/Ci->Kt;
            Ci->Fdfv[hash](1) = hDynInteractonsFV[ifv].Ft.y/Ci->Kt;
            Ci->Fdfv[hash](2) = hDynInteractonsFV[ifv].Ft.z/Ci->Kt;
        }

        #pragma omp parallel for schedule(static) num_threads(Nc)
        for (size_t ii=0;ii<hComInteractons.size();ii++)
        {
            DEM::CInteracton * Ci = static_cast<DEM::CInteracton *>(Interactons[ii]);
            Ci->Fnet(0)  = hComInteractons[ii].Fnnet.x;
            Ci->Fnet(1)  = hComInteractons[ii].Fnnet.y;
            Ci->Fnet(2)  = hComInteractons[ii].Fnnet.z;
            Ci->Ftnet(0) = hComInteractons[ii].Ftnet.x;
            Ci->Ftnet(1) = hComInteractons[ii].Ftnet.y;
            Ci->Ftnet(2) = hComInteractons[ii].Ftnet.z;
            if (Ci->BothFree) // I had to do this due to the two definitions for BranchVec function
            {
                BranchVec(Ci->P2->x,Ci->P1->x,Ci->Branch,Per);
            }
            else
            {
                Ci->Branch = Ci->P1->x-Ci->P2->x;
            }
            //if (norm(Ci->Branch)<0.5*(Ci->P1->Props.R+Ci->P2->Props.R))
            //{
                //std::cout << Ci->P1->Index << " " << Ci->P2->Index << " " << Ci->P1->x << " " << Ci->P2->x << " " <<std::endl;
            //}
        }
    }
}

inline void Domain::UpdateContactsDevice()
{
    ResetMaxD<<<demaux.nparts/Nthread+1,Nthread>>>(pVertsCU, pVertsoCU, pMaxDCU, pParticlesCU, pDynParticlesCU, pdemaux);
    DnLoadDevice(Nproc,true);
    UpdateContacts();
    UpLoadDevice(Nproc,false);
}
#endif
}; // namespace DEM

#endif // MECHSYS_DEM_DOMAIN_H
