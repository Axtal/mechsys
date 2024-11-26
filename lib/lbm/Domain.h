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


#ifndef MECHSYS_LBM_DOMAIN_H
#define MECHSYS_LBM_DOMAIN_H

// STD
#include <map>
#include <vector>
#include <utility>
#include <set>
#include <ctime>
#include <ratio>
#include <chrono>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/lbm/Lattice.h>
#include <mechsys/lbm/Interacton.h>
//#include <mechsys/mesh/mesh.h>
//#include <mechsys/util/util.h>
//#include <mechsys/util/maps.h>
#include <mechsys/mesh/unstructured.h>
//#include <mechsys/util/tree.h>

using std::set;
using std::map;
using std::pair;
using std::make_pair;

namespace LBM
{

struct ParticleCellPair
{
    size_t ICell;         ///< Index of the cell
    size_t IPar;          ///< Index of the particle
    Array<size_t> IGeo;   ///< Array of index of the geometric feature
};

struct MtData;

class Domain
{
public:
    // Transformation matrices for MRT LBM
    static const double   MD2Q5       [5][5]; ///< MRT transformation matrix (D2Q5)
    static const double   MD2Q9       [9][9]; ///< MRT transformation matrix (D2Q9)
    static const double   MD3Q15    [15][15]; ///< MRT transformation matrix (D3Q15)
    static const double   MD3Q19    [19][19]; ///< MRT transformation matrix (D3Q19)
    
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructors
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    Array<double>         nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Special constructor with only one component, the parameters are the same as above
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    double                nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step
    
    //Utility methods
    void BoundingBox       (Vec3_t & minX, Vec3_t & maxX);                                                      ///< Defines the rectangular box that encloses the particles.
    void Center            (Vec3_t C = Vec3_t(0.0,0.0,0.0));                                                    ///< Centers the domain around C
    void SetProps          (Dict & D);                                                                          ///< Set the properties of individual grains by dictionaries
    void Clusters          ();                                                                                  ///< Check the bounded particles in the domain and how many connected clusters are still present
    
    //Methods for adding particles
    void AddSphere       (int Tag, Vec3_t const & X, double R, double rho);                                                                              ///< Add sphere
    void GenSpheresBox   (int Tag, Vec3_t const & X0, Vec3_t const & X1, double R, double rho, size_t Randomseed, double fraction, double RminFraction); ///< Create an array of spheres
    void AddCube         (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);                                ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddRecBox       (int Tag, Vec3_t const & X, Vec3_t const & L, double R, double rho, double Angle=0, Vec3_t * Axis=NULL);                        ///< Add a rectangular box with dimensions given by the vector L
    void AddOcta         (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);                                ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddTetra        (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);                                ///< Add a octahedron at position X with spheroradius R, side of length L and density rho
    void AddPlane        (int Tag, Vec3_t const & X, double R, double Lx,double Ly , double rho, double Angle=0, Vec3_t * Axis=NULL);                    ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void GenBox          (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, double rho, bool Cohesion=false);                        ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBox          (int InitialTag, Vec3_t const & X0, Vec3_t const & X1, double R, double Cf, double rho);                                        ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBoundingBox  (int InitialTag, double R, double Cf, double rho, bool Cohesion=false);                                                         ///< Generate o bounding box enclosing the previous included particles.
    void GenFromMesh     (Mesh::Generic & M, double R, double rho, bool cohesion=false, bool MC=true, double thickness = 0.0);                           ///< Generate particles from a FEM mesh generator
    void AddVoroCell     (int Tag, voro::voronoicell & VC, double R, double rho, bool Erode, Vec3_t nv = iVec3_t(1.0,1.0,1.0));                          ///< Add a single voronoi cell, it should be built before tough
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bool Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                                                ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bVec3_t Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                                             ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni, Periodic conditions are chosen for each particle
    void AddFromJson (int Tag, char const * Filename, double R, double rho, double scale,bool Erode = false);                                            ///< Add a particle generated from Json mesh
    void AddDisk(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt);                      ///< Add a disk element in 2D
    
    // Access methods
    DEM::Particle       * GetParticle  (int Tag, bool Check=true);       ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    DEM::Particle const & GetParticle  (int Tag, bool Check=true) const; ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    void                 DelParticles  (Array<int> const & Tags);        ///< Delete particle
    //Methods
#ifdef USE_HDF5
    void WriteBF           (char const * FileKey);  ///< Save a h5 with branch and force information
    void WriteFrac         (char const * FileKey);  ///< Save a xdmf file for fracture visualization
    void WriteXDMF         (char const * FileKey);  ///< Write the domain data in xdmf file
    void WriteXDMF_D       (char const * FileKey);  ///< Write the domain data in xdmf file with double precision
    void Load              (char const * FileKey);  ///< Load particle data from Mechsys DEM
#endif
    void UpdateLinkedCells ();                                                                                  ///< Update the linked cells

    void Initialize       (double dt=0.0);                                                                                              ///< Set the particles to a initial state and asign the possible insteractions
    void ApplyForce       (size_t n = 0, size_t Np = 1, bool MC=false);                                                                 ///< Apply the interaction forces
    void CollideMRT       (size_t n = 0, size_t Np = 1);                                                                                ///< Apply the collision operator with DEM particles in the case of single component fluid
    void CollideSC        (size_t n = 0, size_t Np = 1);                                                                                ///< Apply the collision operator with DEM particles in the case of single component fluid
    void CollideMC        (size_t n = 0, size_t Np = 1);                                                                                ///< Apply the collision operator with DEM particles in the case of multiple component fluid
    void CollideNoPar     (size_t n = 0, size_t Np = 1);                                                                                ///< Apply the collision operator for the case of no DEM particles
    void ImprintLatticeSC (size_t n = 0, size_t Np = 1);                                                                                ///< Imprint the DEM particles into the lattices when there is a single fluid component
    void ImprintLatticeMC (size_t n = 0, size_t Np = 1);                                                                                ///< Imprint the DEM particles into the lattices when there are multiple fluids
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                                ///< Solve the Domain dynamics
    void ResetContacts();                                                                                                             ///< Reset contacts for verlet method DEM
    void ResetDisplacements();                                                                                                        ///< Reset the displacements for the verlet method DEM
    double  MaxDisplacement();                                                                                                        ///< Give the maximun displacement of DEM particles

#ifdef USE_OMP
    Array<pair<size_t, size_t> >                ListPosPairs;         ///< List of all possible particles pairs
    Array<pair<size_t, size_t> >                List2DPairs;          ///< List of all possible particles pairs in 2D
    iVec3_t                                     LCellDim;             ///< Dimensions of the linked cell array
    Array<Array <size_t> >                      LinkedCell;           ///< Linked Cell array for optimization.
    Vec3_t                                      LCxmin;               ///< Bounding box low   limit for the linked cell array
    Vec3_t                                      LCxmax;               ///< Bounding box upper limit for the linked cell array
#endif
    //Data
    bool                                         Initialized;         ///< System (particles and interactons) initialized ?
    bool                                              RotPar;         ///< Check if particles should be rotated, useful if particle displacements are small
    bool                                              PrtVec;         ///< Print Vector data into the xdmf-h5 files
    bool                                              PrtDou;         ///< Use double precision in h5 output files
    bool                                              PrtPer;         ///< Print percolation parameter when applicable
    bool                                            Finished;         ///< Has the simulation finished
    bool                                              Dilate;         ///< True if eroded particles should be dilated for visualization
    Array<size_t>                                    FreePar;         ///< Particles that are free
    Array<size_t>                                  NoFreePar;         ///< Particles that are not free
    String                                           FileKey;         ///< File Key for output files
    Array <Lattice>                                      Lat;         ///< Fluid Lattices
    Array <DEM::Particle *>                        Particles;         ///< Array of Particles
    Array <LBM::Disk     *>                            Disks;         ///< Array of Disks for 2D calculation
    Array <DEM::Interacton *>                    Interactons;         ///< Array of interactons
    Array <LBM::DiskPair   *>                      DiskPairs;         ///< Array of interactons for 2D calculation
    Array <DEM::CInteracton *>                  CInteractons;         ///< Array of valid  collision interactons
    Array <DEM::BInteracton *>                  BInteractons;         ///< Cohesion interactons
    Array <iVec3_t>                                CellPairs;         ///< Pairs of cells
    Array <ParticleCellPair>                    ParCellPairs;         ///< Pairs of cells and particles
    set<pair<DEM::Particle *, DEM::Particle *> > Listofpairs;         ///< List of pair of particles associated per interacton for memory optimization
    set<pair<LBM::Disk *, LBM::Disk *> >     ListofDiskPairs;         ///< List of pair of disks associated per interacton for memory optimization
    double                                              Time;         ///< Time of the simulation
    double                                                dt;         ///< Timestep
    double                                             dtdem;         ///< Timestep for DEM
    double                                             Alpha;         ///< Verlet distance
    double                                              Beta;         ///< Binmultiplier
    double                                              Gmix;         ///< Interaction constant for the mixture
    double                                            Voltot;         ///< toala particle volume
    double                                           MaxDmax;         ///< Maximun value for the radious of the spheres surronding each particle
    double                                                Ms;         ///< Total mass of the particles
    double                                                Vs;         ///< Volume occupied by the grains
    double                                                Sc;         ///< Smagorinsky constant
    double                                             Fconv;         ///< Force conversion factor
    Array <double>                                       EEk;         ///< Diadic velocity tensor trace
    void *                                          UserData;         ///< User Data
    Mat_t                                                  M;         ///< Transformation matrix to momentum space for MRT calculations
    Mat_t                                               Minv;         ///< Inverse Transformation matrix to momentum space for MRT calculations
    Vec_t                                                  S;         ///< Vector of relaxation times for MRT
    size_t                                             Nproc;         ///< Number of cores for multithreading
    size_t                                           idx_out;         ///< The discrete time step
    size_t                                              Step;         ///< The space step to reduce the size of the h5 file for visualization
    Array<Array <int> >                       Listofclusters;         ///< List of particles belonging to bounded clusters (applies only for cohesion simulations)
    MtData *                                             MTD;         ///< Multithread data
};

struct MtData
{
    size_t                        ProcRank; ///< Rank of the thread
    size_t                          N_Proc; ///< Total number of threads
    LBM::Domain *                      Dom; ///< Pointer to the lbm domain
    double                             Dmx; ///< Maximun displacement
    double                              dt; ///< Time step
    Array<pair<size_t,size_t> >         LC; ///< A temporal list of new contacts
    Array<size_t>                      LCI; ///< A temporal array of posible Cinteractions
    Array<size_t>                      LCB; ///< A temporal array of posible Binteractions
    Array<ParticleCellPair>            LPC; ///< A temporal array of possible particle cell contacts
    Array<std::pair<iVec3_t,size_t> >  LLC; ///< A temporal array of possible linked cells locations
    Array<std::pair<size_t,size_t> >   LPP; ///< A temporal array of possible partcle types
};

inline Domain::Domain(LBMethod Method, Array<double> nu, iVec3_t Ndim, double dx, double Thedt)
{
    Initialized = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (nu.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one fluid please");
    if (Ndim(2) >1&&Method==D2Q9)  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (Ndim(2)==1&&Method==D3Q15) throw new Fatal("LBM::Domain: Ndim(2) is 1. Either change the method to D2Q9 or increase the z-dimension");
    for (size_t i=0;i<nu.Size();i++)
    {
        Lat.Push(Lattice(Method,nu[i],Ndim,dx,Thedt));
    }
    Time   = 0.0;
    dt     = Thedt;
    dtdem  = 0.0;
    Alpha  = 10.0;
    Beta   = 2.0;
    Step   = 1;
    if (nu.Size()>1) Sc = 0.0;
    else             Sc = 0.17;
    PrtVec = true;
    PrtDou = false;
    RotPar = true;
    Fconv  = 1.0;


    EEk.Resize(Lat[0].Cells[0]->Nneigh);
    for (size_t k=0;k<Lat[0].Cells[0]->Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(Lat[0].Cells[0]->C[k][n]*Lat[0].Cells[0]->C[k][m]);
        }
    }

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Size()*Lat[0].Ncells,TERM_RST);
}

inline Domain::Domain(LBMethod Method, double nu, iVec3_t Ndim, double dx, double Thedt)
{
    Initialized = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    Lat.Push(Lattice(Method,nu,Ndim,dx,Thedt));
    if (Ndim(2) >1&&Method==D2Q9)  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (Ndim(2)==1&&Method==D3Q15) throw new Fatal("LBM::Domain: Ndim(2) is greater than 1. Either change the method to D2Q9 or increse the z-dimension");
    Time   = 0.0;
    dt     = Thedt;
    dtdem  = 0.0;
    Alpha  = 10.0;
    Beta   = 2.0;
    Step   = 1;
    Sc     = 0.17;
    PrtVec = true;
    PrtDou = false;
    RotPar = true;
    Fconv  = 1.0;

    EEk.Resize(Lat[0].Cells[0]->Nneigh);
    for (size_t k=0;k<Lat[0].Cells[0]->Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(Lat[0].Cells[0]->C[k][n]*Lat[0].Cells[0]->C[k][m]);
        }
    }

    if (Method==D3Q15)
    {
        size_t Nneigh = 15;
        M.Resize(Nneigh,Nneigh);
        Minv.Resize(Nneigh,Nneigh);
        for (size_t n=0;n<Nneigh;n++)
        for (size_t m=0;m<Nneigh;m++)
        {
            M(n,m) = MD3Q15[n][m];
        }
        Inv(M,Minv);
        double tau = 3.0*nu*dt/(dx*dx)+0.5;
        double s   = 8.0*(2.0-1.0/tau)/(8.0-1.0/tau);
        S.Resize(Nneigh);
        S = 0.0,1.0/tau,1.0/tau,0.0,s,0.0,s,0.0,s,1.0/tau,1.0/tau,1.0/tau,1.0/tau,1.0/tau,s;
    }

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Size()*Lat[0].Ncells,TERM_RST);
}

#ifdef USE_HDF5

inline void Domain::WriteBF (char const * FileKey)
{

    size_t n_fn = 0;

    for (size_t i=0;i<CInteractons.Size();i++)
    {
        //if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree())) n_fn++;
        if (norm(CInteractons[i]->Fnet)>1.0e-12) n_fn++;
    }

    if (n_fn==0) return;
    
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    float  *  Fnnet = new float[3*n_fn];
    float  *  Ftnet = new float[3*n_fn];
    float  * Branch = new float[3*n_fn];
    float  *   Orig = new float[3*n_fn];
    float  *    Fn  = new float[  n_fn];
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
            Branch[3*idx  ] = float(CInteractons[i]->P1->x(0)-CInteractons[i]->P2->x(0));
            Branch[3*idx+1] = float(CInteractons[i]->P1->x(1)-CInteractons[i]->P2->x(1)); 
            Branch[3*idx+2] = float(CInteractons[i]->P1->x(2)-CInteractons[i]->P2->x(2)); 
            //Orig  [3*idx  ] = 0.5*float(CInteractons[i]->P1->x(0)+CInteractons[i]->P2->x(0));
            //Orig  [3*idx+1] = 0.5*float(CInteractons[i]->P1->x(1)+CInteractons[i]->P2->x(1)); 
            //Orig  [3*idx+2] = 0.5*float(CInteractons[i]->P1->x(2)+CInteractons[i]->P2->x(2)); 
            Orig  [3*idx  ] = float(CInteractons[i]->P2->x(0));
            Orig  [3*idx+1] = float(CInteractons[i]->P2->x(1)); 
            Orig  [3*idx+2] = float(CInteractons[i]->P2->x(2)); 
            Fn    [idx]     = float(norm(CInteractons[i]->Fnet));
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
    dsname.Printf("Branch");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Branch);
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Orig);
    dims[0] = n_fn;
    dsname.Printf("Fn");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Fn);
    dsname.Printf("ID1");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,ID1   );
    dsname.Printf("ID2");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,ID2   );


    delete [] Fnnet;
    delete [] Ftnet;
    delete [] Branch;
    delete [] Orig;
    delete [] Fn;
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
    oss << "     <Attribute Name=\"Branch\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Branch \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Fn\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Fn \n";
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
            DEM::Particle * P1 = BInteractons[i]->P1;
            DEM::Particle * P2 = BInteractons[i]->P2;
            DEM::Face     * F1 = P1->Faces[BInteractons[i]->IF1];
            DEM::Face     * F2 = P2->Faces[BInteractons[i]->IF2];
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
        DEM::Particle * P1 = BInteractons[i]->P1;
        DEM::Particle * P2 = BInteractons[i]->P2;
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

inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Lat[0].Ndim[0]/Step;
    size_t  Ny = Lat[0].Ndim[1]/Step;
    size_t  Nz = Lat[0].Ndim[2]/Step;
    for (size_t j=0;j<Lat.Size();j++)
    {
        // Creating data sets
        float * Density   = new float[  Nx*Ny*Nz];
        float * Gamma     = new float[  Nx*Ny*Nz];
        float * Per       = new float[  Nx*Ny*Nz];
        float * Vvec      = new float[3*Nx*Ny*Nz];

        size_t i=0;
        for (size_t m=0;m<Lat[0].Ndim(2);m+=Step)
        for (size_t l=0;l<Lat[0].Ndim(1);l+=Step)
        for (size_t n=0;n<Lat[0].Ndim(0);n+=Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            double per    = 0.0;
            Vec3_t vel    = OrthoSys::O;

            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Rho;
                gamma+= std::max(Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Gamma,(double) Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->IsSolid);
                per  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Pf;
                vel  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel;
                vel  += dt*0.5*Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->BForce;
            }
            rho  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            per  /= Step*Step*Step;
            vel  /= Step*Step*Step;
            //Gamma   [i]  = (float) Lat[j].Cells[i]->IsSolid? 1.0: gamma;
            Gamma   [i]  = gamma;
            Per     [i]  = per;
            //Density [i]  = (float) rho;
            //Vvec[3*i  ]  = (float) vel(0);
            //Vvec[3*i+1]  = (float) vel(1);
            //Vvec[3*i+2]  = (float) vel(2);
            Density [i]  = (float) rho   *(1.0-Gamma[i]);
            Vvec[3*i  ]  = (float) vel(0)*(1.0-Gamma[i]);
            Vvec[3*i+1]  = (float) vel(1)*(1.0-Gamma[i]);
            Vvec[3*i+2]  = (float) vel(2)*(1.0-Gamma[i]);
            i++;
        }

        //Write the data
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Gamma   );
            if (PrtPer)
            {
                dsname.Printf("Per");
                H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Per );
            }
        }
        if (PrtVec)
        {
            dims[0] = 3*Nx*Ny*Nz;
            dsname.Printf("Velocity_%d",j);
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vvec    );
        }
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
        delete [] Gamma   ;
        delete [] Per     ;
        delete [] Vvec    ;
    }


    size_t N_Faces = 0;
    size_t N_Verts = 0;
    //Writing particle data
    if (Particles.Size()>0)
    {
        for (size_t i=0; i<Particles.Size(); i++) 
        { 
            for (size_t j=0;j<Particles[i]->Faces.Size();j++)
            {
                N_Faces += Particles[i]->Faces[j]->Edges.Size();
            }
            N_Verts += Particles[i]->Verts.Size() + Particles[i]->Faces.Size();
        }
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
                DEM::Particle * Pa = Particles[i];
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
                    Verts[n_verts++] = (float) Vres[j](0);
                    Verts[n_verts++] = (float) Vres[j](1);
                    Verts[n_verts++] = (float) Vres[j](2);
                }
                size_t n_reff = n_verts/3;
                for (size_t j=0;j<Pa->FaceCon.Size();j++)
                {
                    Vec3_t C,N;
                    Pa->Faces[j]->Centroid(C);
                    Pa->Faces[j]->Normal(N);
                    Verts[n_verts++] = (float) C(0) + multiplier*Pa->Props.R*N(0);
                    Verts[n_verts++] = (float) C(1) + multiplier*Pa->Props.R*N(1);
                    Verts[n_verts++] = (float) C(2) + multiplier*Pa->Props.R*N(2);
                    //Verts[n_verts++] = (float) C(0);
                    //Verts[n_verts++] = (float) C(1);
                    //Verts[n_verts++] = (float) C(2);
                    for (size_t k=0;k<Pa->FaceCon[j].Size();k++)
                    {
                        size_t nin = Pa->FaceCon[j][k];
                        size_t nen = Pa->FaceCon[j][(k+1)%Pa->FaceCon[j].Size()];
                        FaceCon[n_faces++] = (int) n_reff + j;  
                        FaceCon[n_faces++] = (int) n_refv + nin;
                        FaceCon[n_faces++] = (int) n_refv + nen;

                        //Writing the attributes
                        Tags  [n_attrs] = (int)   Pa->Tag;
                        Clus  [n_attrs] = size_t(Pa->Cluster);
                        Vel   [n_attrs] = (float) norm(Pa->v);
                        Ome   [n_attrs] = (float) norm(Pa->w);
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
            //dims[0] = 9*N_Faces;
            //dsname.Printf("Stress");
            //H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Stress);
            
            //Erasing the data
            delete [] Verts;
            delete [] FaceCon;
            delete [] Tags;
            delete [] Clus;
            delete [] Vel;
            delete [] Ome;
            //delete [] Stress;
        }
        //Creating data sets
        float * Radius = new float[  Particles.Size()];
        float * Posvec = new float[3*Particles.Size()];
        float * Velvec = new float[3*Particles.Size()];
        float * Omevec = new float[3*Particles.Size()];
        float * Aomvec = new float[3*Particles.Size()];
        float * Amovec = new float[3*Particles.Size()];
        float * Ufovec = new float[3*Particles.Size()];
        float * Forvec = new float[3*Particles.Size()];
        float * Inertm = new float[6*Particles.Size()];
        float * Ekin   = new float[  Particles.Size()];
        int   * Tags   = new int  [  Particles.Size()];

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


            Particles[i]->Verts.Size()==1 ? Radius[i] = float(Particles[i]->Dmax) : Radius[i] = 0.0;
            Posvec[3*i  ] = (float) Particles[i]->x(0);
            Posvec[3*i+1] = (float) Particles[i]->x(1);
            Posvec[3*i+2] = (float) Particles[i]->x(2);
            Velvec[3*i  ] = (float) Particles[i]->v(0);
            Velvec[3*i+1] = (float) Particles[i]->v(1);
            Velvec[3*i+2] = (float) Particles[i]->v(2);
            Omevec[3*i  ] = (float) Ome(0);
            Omevec[3*i+1] = (float) Ome(1);
            Omevec[3*i+2] = (float) Ome(2);
            Aomvec[3*i  ] = float(Ao(0));
            Aomvec[3*i+1] = float(Ao(1)); 
            Aomvec[3*i+2] = float(Ao(2)); 
            Amovec[3*i  ] = float(L(0));
            Amovec[3*i+1] = float(L(1)); 
            Amovec[3*i+2] = float(L(2)); 
            Ufovec[3*i  ] = (float) Particles[i]->F(0);
            Ufovec[3*i+1] = (float) Particles[i]->F(1); 
            Ufovec[3*i+2] = (float) Particles[i]->F(2); 
            Forvec[3*i  ] = (float) Particles[i]->Flbm(0);
            Forvec[3*i+1] = (float) Particles[i]->Flbm(1);
            Forvec[3*i+2] = (float) Particles[i]->Flbm(2);
            Inertm[6*i  ] = (float) Inertiar(0,0);
            Inertm[6*i+1] = (float) Inertiar(0,1);
            Inertm[6*i+2] = (float) Inertiar(0,2);
            Inertm[6*i+3] = (float) Inertiar(1,1);
            Inertm[6*i+4] = (float) Inertiar(1,2);
            Inertm[6*i+5] = (float) Inertiar(2,2);
            Ekin  [i]     = float(Particles[i]->Ekin+Particles[i]->Erot);
            Tags  [i]     = (int)   Particles[i]->Tag;
        }

        hsize_t dims[1];
        String dsname;
        dims[0] = 6*Particles.Size();
        dsname.Printf("Inertia");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Inertm);
        dims[0] = 3*Particles.Size();
        dsname.Printf("Position");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
        dsname.Printf("PVelocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
        dsname.Printf("PAngVel");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Omevec);
        dsname.Printf("PAngacceleration");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Aomvec);
        dsname.Printf("PAngMomentum");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Amovec);
        dsname.Printf("PUForce");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ufovec);
        dsname.Printf("PForce");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Forvec);
        dims[0] = Particles.Size();
        dsname.Printf("Radius");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
        dsname.Printf("PEkin");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ekin);
        dsname.Printf("PTag");
        H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tags  );

        delete [] Radius;
        delete [] Posvec;
        delete [] Velvec;
        delete [] Omevec;
        delete [] Aomvec;
        delete [] Amovec;
        delete [] Ufovec;
        delete [] Forvec;
        delete [] Inertm;
        delete [] Ekin;
        delete [] Tags;
    }
    
    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    //std::cout << "2" << std::endl;

    if (Lat[0].Ndim(2)==1)
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Lat[0].Ndim(1) << " " << Lat[0].Ndim(0) << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Lat.Size();j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        if (PrtVec)
        {
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        }
        if (PrtPer)
        {
        oss << "     <Attribute Name=\"Per\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Per\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
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
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*Lat[0].dx << " " << Step*Lat[0].dx  << " " << Step*Lat[0].dx  << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Lat.Size();j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        if (PrtVec)
        {
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        }
        if (PrtPer)
        {
        oss << "     <Attribute Name=\"Per\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Per\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        if(Particles.Size()>0)
        {
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
        //oss << "     </Attribute>\n";
        //oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Cell\">\n";
        //oss << "       <DataItem Dimensions=\"" << N_Faces << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        //oss << "        " << fn.CStr() <<":/AngVelocity \n";
        //oss << "       </DataItem>\n";
        //oss << "     </Attribute>\n";
        //oss << "     <Attribute Name=\"Stress\" AttributeType=\"Tensor\" Center=\"Cell\">\n";
        //oss << "       <DataItem Dimensions=\"" << N_Faces << " 3 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        //oss << "        " << fn.CStr() <<":/Stress \n";
        //oss << "       </DataItem>\n";
        //oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        }
        oss << "   <Grid Name=\"DEM_Center\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
        oss << "     <Geometry GeometryType=\"XYZ\">\n";
        oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
        oss << "        " << fn.CStr() <<":/Position \n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
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
        oss << "        " << fn.CStr() <<":/PAngVel\n";
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
        oss << "     <Attribute Name=\"Uforce\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/PUForce\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Force\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/PForce\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Inertia\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Inertia\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        }
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::WriteXDMF_D(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Lat[0].Ndim[0]/Step;
    size_t  Ny = Lat[0].Ndim[1]/Step;
    size_t  Nz = Lat[0].Ndim[2]/Step;
    for (size_t j=0;j<Lat.Size();j++)
    {
        // Creating data sets
        double * Density   = new double[  Nx*Ny*Nz];
        double * Gamma     = new double[  Nx*Ny*Nz];
        double * Per       = new double[  Nx*Ny*Nz];
        double * Vvec      = new double[3*Nx*Ny*Nz];

        size_t i=0;
        for (size_t m=0;m<Lat[0].Ndim(2);m+=Step)
        for (size_t l=0;l<Lat[0].Ndim(1);l+=Step)
        for (size_t n=0;n<Lat[0].Ndim(0);n+=Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            double per    = 0.0;
            Vec3_t vel    = OrthoSys::O;

            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Rho;
                gamma+= std::max(Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Gamma,(double) Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->IsSolid);
                per  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Pf;
                vel  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel;
                vel  += dt*0.5*Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->BForce;
            }
            rho  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            per  /= Step*Step*Step;
            vel  /= Step*Step*Step;
            //Gamma   [i]  = (double) Lat[j].Cells[i]->IsSolid? 1.0: gamma;
            Gamma   [i]  = gamma;
            Per     [i]  = per;
            //Density [i]  = (double) rho;
            //Vvec[3*i  ]  = (double) vel(0);
            //Vvec[3*i+1]  = (double) vel(1);
            //Vvec[3*i+2]  = (double) vel(2);
            Density [i]  = (double) rho   *(1.0-Gamma[i]);
            Vvec[3*i  ]  = (double) vel(0)*(1.0-Gamma[i]);
            Vvec[3*i+1]  = (double) vel(1)*(1.0-Gamma[i]);
            Vvec[3*i+2]  = (double) vel(2)*(1.0-Gamma[i]);
            i++;
        }

        //Write the data
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Density );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Gamma   );
            if (PrtPer)
            {
                dsname.Printf("Per");
                H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Per );
            }
        }
        if (PrtVec)
        {
            dims[0] = 3*Nx*Ny*Nz;
            dsname.Printf("Velocity_%d",j);
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vvec    );
        }
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
        delete [] Gamma   ;
        delete [] Per     ;
        delete [] Vvec    ;
    }


    size_t N_Faces = 0;
    size_t N_Verts = 0;
    //Writing particle data
    if (Particles.Size()>0)
    {
        for (size_t i=0; i<Particles.Size(); i++) 
        { 
            for (size_t j=0;j<Particles[i]->Faces.Size();j++)
            {
                N_Faces += Particles[i]->Faces[j]->Edges.Size();
            }
            N_Verts += Particles[i]->Verts.Size() + Particles[i]->Faces.Size();
        }
        if (N_Faces>0)
        {

            //Geometric information
            double  * Verts   = new double [3*N_Verts];
            int    * FaceCon = new int   [3*N_Faces];
            
            //Atributes
            int    * Tags    = new int   [  N_Faces];
            int    * Clus    = new int   [  N_Faces];
            double  * Vel     = new double [  N_Faces];
            double  * Ome     = new double [  N_Faces];
            //double  * Stress  = new double [9*N_Faces];

            size_t n_verts = 0;
            size_t n_faces = 0;
            size_t n_attrs = 0;
            //size_t n_attrv = 0;
            //size_t n_attrt = 0;
            for (size_t i=0;i<Particles.Size();i++)
            {
                DEM::Particle * Pa = Particles[i];
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
                    //Verts[n_verts++] = (double) (*Pa->Verts[j])(0);
                    //Verts[n_verts++] = (double) (*Pa->Verts[j])(1);
                    //Verts[n_verts++] = (double) (*Pa->Verts[j])(2);
                    Verts[n_verts++] = (double) Vres[j](0);
                    Verts[n_verts++] = (double) Vres[j](1);
                    Verts[n_verts++] = (double) Vres[j](2);
                }
                size_t n_reff = n_verts/3;
                for (size_t j=0;j<Pa->FaceCon.Size();j++)
                {
                    Vec3_t C,N;
                    Pa->Faces[j]->Centroid(C);
                    Pa->Faces[j]->Normal(N);
                    Verts[n_verts++] = (double) C(0) + multiplier*Pa->Props.R*N(0);
                    Verts[n_verts++] = (double) C(1) + multiplier*Pa->Props.R*N(1);
                    Verts[n_verts++] = (double) C(2) + multiplier*Pa->Props.R*N(2);
                    //Verts[n_verts++] = (double) C(0);
                    //Verts[n_verts++] = (double) C(1);
                    //Verts[n_verts++] = (double) C(2);
                    for (size_t k=0;k<Pa->FaceCon[j].Size();k++)
                    {
                        size_t nin = Pa->FaceCon[j][k];
                        size_t nen = Pa->FaceCon[j][(k+1)%Pa->FaceCon[j].Size()];
                        FaceCon[n_faces++] = (int) n_reff + j;  
                        FaceCon[n_faces++] = (int) n_refv + nin;
                        FaceCon[n_faces++] = (int) n_refv + nen;

                        //Writing the attributes
                        Tags  [n_attrs] = (int)   Pa->Tag;
                        Clus  [n_attrs] = size_t(Pa->Cluster);
                        Vel   [n_attrs] = (double) norm(Pa->v);
                        Ome   [n_attrs] = (double) norm(Pa->w);
                        n_attrs++;

                        //Vel [n_attrv  ] = (double) Pa->v(0);
                        //Vel [n_attrv+1] = (double) Pa->v(1);
                        //Vel [n_attrv+2] = (double) Pa->v(2);
                        //n_attrv += 3;

                        //Stress[n_attrt  ] = (double) Pa->M(0,0);
                        //Stress[n_attrt+1] = (double) Pa->M(1,0);
                        //Stress[n_attrt+2] = (double) Pa->M(2,0);
                        //Stress[n_attrt+3] = (double) Pa->M(0,1);
                        //Stress[n_attrt+4] = (double) Pa->M(1,1);
                        //Stress[n_attrt+5] = (double) Pa->M(2,1);
                        //Stress[n_attrt+6] = (double) Pa->M(0,2);
                        //Stress[n_attrt+7] = (double) Pa->M(1,2);
                        //Stress[n_attrt+8] = (double) Pa->M(2,2);
                        //n_attrt += 9;
                    }
                }
            }

            //Write the data
            hsize_t dims[1];
            String dsname;
            dims[0] = 3*N_Verts;
            dsname.Printf("Verts");
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Verts);
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
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vel);
            dims[0] = N_Faces;
            dsname.Printf("AngVelocity");
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Ome);
            //dims[0] = 9*N_Faces;
            //dsname.Printf("Stress");
            //H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Stress);
            
            //Erasing the data
            delete [] Verts;
            delete [] FaceCon;
            delete [] Tags;
            delete [] Clus;
            delete [] Vel;
            delete [] Ome;
            //delete [] Stress;
        }
        //Creating data sets
        double * Radius = new double[  Particles.Size()];
        double * Posvec = new double[3*Particles.Size()];
        double * Velvec = new double[3*Particles.Size()];
        double * Omevec = new double[3*Particles.Size()];
        double * Forvec = new double[3*Particles.Size()];
        double * Torvec = new double[3*Particles.Size()];
        int   * Tags   = new int  [  Particles.Size()];
        for (size_t i=0;i<Particles.Size();i++)
        {
            Vec3_t Ome,Tor;
            Rotation(Particles[i]->w,Particles[i]->Q,Ome);
            Rotation(Particles[i]->T,Particles[i]->Q,Tor);
            Particles[i]->Verts.Size()==1 ? Radius[i] = double(Particles[i]->Dmax) : Radius[i] = 0.0;
            Posvec[3*i  ] = (double) Particles[i]->x(0);
            Posvec[3*i+1] = (double) Particles[i]->x(1);
            Posvec[3*i+2] = (double) Particles[i]->x(2);
            Velvec[3*i  ] = (double) Particles[i]->v(0);
            Velvec[3*i+1] = (double) Particles[i]->v(1);
            Velvec[3*i+2] = (double) Particles[i]->v(2);
            Omevec[3*i  ] = (double) Ome(0);
            Omevec[3*i+1] = (double) Ome(1);
            Omevec[3*i+2] = (double) Ome(2);
            Forvec[3*i  ] = (double) Particles[i]->F(0);
            Forvec[3*i+1] = (double) Particles[i]->F(1);
            Forvec[3*i+2] = (double) Particles[i]->F(2);
            Torvec[3*i  ] = (double) Tor(0);
            Torvec[3*i+1] = (double) Tor(1);
            Torvec[3*i+2] = (double) Tor(2);
            Tags  [i]     = (int)   Particles[i]->Tag;
        }

        hsize_t dims[1];
        dims[0] = 3*Particles.Size();
        String dsname;
        dsname.Printf("Position");
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Posvec);
        dsname.Printf("PVelocity");
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Velvec);
        dsname.Printf("PAngVel");
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Omevec);
        dsname.Printf("PForce");
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Forvec);
        dsname.Printf("PTorque");
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Torvec);
        dims[0] = Particles.Size();
        dsname.Printf("Radius");
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Radius);
        dsname.Printf("PTag");
        H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tags  );


        delete [] Radius;
        delete [] Posvec;
        delete [] Velvec;
        delete [] Omevec;
        delete [] Forvec;
        delete [] Torvec;
        delete [] Tags  ;
    }
    
    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    //std::cout << "2" << std::endl;

    if (Lat[0].Ndim(2)==1)
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Lat[0].Ndim(1) << " " << Lat[0].Ndim(0) << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Lat.Size();j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        if (PrtVec)
        {
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        }
        if (PrtPer)
        {
        oss << "     <Attribute Name=\"Per\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Per\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
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
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*Lat[0].dx << " " << Step*Lat[0].dx  << " " << Step*Lat[0].dx  << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Lat.Size();j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        if (PrtVec)
        {
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        }
        if (PrtPer)
        {
        oss << "     <Attribute Name=\"Per\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Per\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        if(Particles.Size()>0)
        {
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
        //oss << "     </Attribute>\n";
        //oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Cell\">\n";
        //oss << "       <DataItem Dimensions=\"" << N_Faces << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        //oss << "        " << fn.CStr() <<":/AngVelocity \n";
        //oss << "       </DataItem>\n";
        //oss << "     </Attribute>\n";
        //oss << "     <Attribute Name=\"Stress\" AttributeType=\"Tensor\" Center=\"Cell\">\n";
        //oss << "       <DataItem Dimensions=\"" << N_Faces << " 3 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        //oss << "        " << fn.CStr() <<":/Stress \n";
        //oss << "       </DataItem>\n";
        //oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        }
        oss << "   <Grid Name=\"DEM_Center\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
        oss << "     <Geometry GeometryType=\"XYZ\">\n";
        oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
        oss << "        " << fn.CStr() <<":/Position \n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Radius \n";
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
        oss << "        " << fn.CStr() <<":/PAngVel\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Force\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/PForce\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Torque\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/PTorque\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        }
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
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

        Particles.Push (new DEM::Particle(-1,V,E,F,OrthoSys::O,OrthoSys::O,0.1,1.0));

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

#endif

inline void Domain::ApplyForce(size_t n, size_t Np, bool MC)
{
    size_t Ni = CellPairs.Size()/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = CellPairs.Size() : Fn = (n+1)*Ni;

    if (MC)
    {
#ifdef USE_OMP
        In = 0;
        Fn = CellPairs.Size();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
        for (size_t i=In;i<Fn;i++)
        {
            size_t ind1 = CellPairs[i](0);
            size_t ind2 = CellPairs[i](1);
            size_t vec  = CellPairs[i](2);
            Cell * c = Lat[0].Cells[ind1];
            Cell *nb = Lat[1].Cells[ind2];
            double nb_psi,psi,G=Gmix/(dt*dt);
            bool check = false;
            if (c->IsSolid||fabs(c->Gamma-1.0)<1.0e-12)
            //if (c->IsSolid||c->Gamma>1.0e-12)
            {
                psi    = 1.0;
                G      = Lat[1].Gs*c->Gs;
                check  = true;
            }
            else psi   = c ->Rho;
            if (nb->IsSolid||fabs(nb->Gamma-1.0)<1.0e-12)
            //if (nb->IsSolid||nb->Gamma>1.0e-12)
            {
                nb_psi = 1.0;
                G      = Lat[0].Gs*nb->Gs;
                // this is to ignore forces where both are solid nodes
                if (check) G = 0.0;
            }
            else nb_psi = nb->Rho;
            Vec3_t  BF = -G*psi*nb_psi*c->W[vec]*c->C[vec];

#ifdef USE_OMP
            omp_set_lock        (&c->lck);
#endif
            c ->BForce += BF;
#ifdef USE_OMP
            omp_unset_lock      (&c ->lck);
            omp_set_lock        (&nb->lck);
#endif
            nb->BForce -= BF;
#ifdef USE_OMP
            omp_unset_lock      (&nb->lck);
#endif

            c  = Lat[1].Cells[ind1];
            nb = Lat[0].Cells[ind2];
            G  = Gmix/(dt*dt);
            check = false;
            if (c->IsSolid||fabs(c->Gamma-1.0)<1.0e-12)
            //if (c->IsSolid||c->Gamma>1.0e-12)
            {
                psi    = 1.0;
                G      = Lat[0].Gs;
                check  = true;
            }
            else psi   = c ->Rho;
            if (nb->IsSolid||fabs(nb->Gamma-1.0)<1.0e-12)
            //if (nb->IsSolid||nb->Gamma>1.0e-12)
            {
                nb_psi = 1.0;
                G      = Lat[1].Gs;
                if (check) G = 0.0;
            }
            else nb_psi = nb->Rho;
            BF = -G*psi*nb_psi*c->W[vec]*c->C[vec];

#ifdef USE_OMP
            omp_set_lock        (&c->lck);
#endif
            c ->BForce += BF;
#ifdef USE_OMP
            omp_unset_lock      (&c ->lck);
            omp_set_lock        (&nb->lck);
#endif
            nb->BForce -= BF;
#ifdef USE_OMP
            omp_unset_lock      (&nb->lck);
#endif
        }
    }
    else
    {
#ifdef USE_OMP
        In = 0;
        Fn = CellPairs.Size();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
        for (size_t i=In;i<Fn;i++)
        {
            size_t ind1 = CellPairs[i](0);
            size_t ind2 = CellPairs[i](1);
            size_t vec  = CellPairs[i](2);
            for (size_t j=0;j<Lat.Size();j++)
            {
                bool solid = false;
                Cell * c = Lat[j].Cells[ind1];
                double psi;
                if (c->IsSolid||fabs(c->Gamma-1.0)<1.0e-12)
                //if (c->IsSolid||c->Gamma>0.0)
                {
                    psi   = 1.0;
                    solid = true;
                }
                else fabs(Lat[j].G)>1.0e-12 ? psi = Lat[j].Psi(c->Rho) : psi = c->Rho;
                for (size_t k=0;k<Lat.Size();k++)
                {
                    Cell * nb = Lat[k].Cells[ind2];
                    double nb_psi;
                    if (nb->IsSolid||fabs(nb->Gamma-1.0)<1.0e-12)
                    //if (nb->IsSolid||nb->Gamma>0.0)
                    {
                        nb_psi = 1.0;
                        solid  = true;
                    }
                    else fabs(Lat[j].G)>1.0e-12 ? nb_psi = Lat[k].Psi(nb->Rho) : nb_psi = nb->Rho;
                    double G;
                    solid ? G = Lat[j].Gs*2.0*ReducedValue(nb->Gs,c->Gs) : G = Lat[j].G; 
                    Vec3_t BF(OrthoSys::O);
                    if (j==k)
                    {
                        BF += -G*psi*nb_psi*c->W[vec]*c->C[vec];
                    }
                    else if(!solid)
                    {
                        BF += -Gmix*c->Rho*nb->Rho*c->W[vec]*c->C[vec];
                    }
#ifdef USE_THREAD
                    pthread_mutex_lock  (&c ->lck);
#elif USE_OMP
                    omp_set_lock        (&c ->lck);
#endif
                    c ->BForce += BF;
#ifdef USE_THREAD
                    pthread_mutex_unlock(&c ->lck);
                    pthread_mutex_lock  (&nb->lck);
#elif USE_OMP
                    omp_unset_lock      (&c ->lck);
                    omp_set_lock        (&nb->lck);
#endif
                    nb->BForce -= BF;
#ifdef USE_THREAD
                    pthread_mutex_unlock(&nb->lck);
#elif USE_OMP
                    omp_unset_lock      (&nb->lck);
#endif
                }
            }
        }
    }
}

void Domain::CollideMRT (size_t n, size_t Np)
{
    //std::cout << "SCP" << std::endl;
	size_t Ni = Lat[0].Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Lat[0].Ncells : Fn = (n+1)*Ni;
#ifdef USE_OMP
    In = 0;
    Fn = Lat[0].Ncells;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        Cell * c = Lat[0].Cells[i];
        double rho = c->Rho;
        Vec3_t vel = c->Vel;
        if (!c->IsSolid)
        {
            double *f = c->F;
            double *ft= c->Ftemp;
            double fneq[c->Nneigh];

            int n=c->Nneigh,m=1;
            double a = 1.0,b = 0.0,Cs = c->Cs;
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
            
            double Bn;
            //rho<10e-12 ? Bn =0.0 : Bn = (c->Gamma*(Lat[j].Tau-0.5))/((1.0-c->Gamma)+(Lat[j].Tau-0.5));
            //rho<10e-12 ? Bn =0.0 : Bn = c->Gamma;
            Bn = (c->Gamma*(Lat[0].Tau-0.5))/((1.0-c->Gamma)+(Lat[0].Tau-0.5));
            //Bn = c->Gamma;
            //Bn = floor(c->Gamma);

            bool valid  = true;
            double alphal = 1.0;
            double alphat = 1.0;
            size_t num = 0;
            while (valid)
            {
                num++;
                valid = false;
                alphal  = alphat;
                for (size_t k=0;k<c->Nneigh;k++)
                {
                    //double FDeqn = c->Feq(k,DV,rho);
                    //c->Ftemp[k] = c->F[k] - alphal*((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]);
                    c->Ftemp[k] = c->F[k] - alphal*((1.0 - Bn)*fneq[k] - Bn*c->Omeis[k] - dt*3.0*c->W[k]*dot(c->BForce,c->C[k])/Cs);
                    if (c->Ftemp[k]<-1.0e-12&&num<2)
                    {
                        //double temp = fabs(c->F[k]/((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]));
                        double temp = fabs(c->F[k]/((1.0 - Bn)*fneq[k] - Bn*c->Omeis[k] - dt*3.0*c->W[k]*dot(c->BForce,c->C[k])/Cs));
                        if (temp<alphat) alphat = temp;
                        valid = true;
                    }
                }
            }
            for (size_t k=0;k<c->Nneigh;k++)
            {
                if (std::isnan(c->Ftemp[k]))
                {
                    //c->Gamma = 2.0;
                    #ifdef USE_HDF5
                    //WriteXDMF("error");
                    #endif
                    //std::cout << "CollideSC: Nan found, resetting" << std::endl;
                    std::cout << c->Density() << " " << c->BForce << " " << num << " " << alphat << " " << c->Index << " " << c->IsSolid << " " << c->Gamma << " " << k << " " << std::endl;
                    c->Ftemp[k] = c->F[k];
                    //throw new Fatal("Domain::Collide: Body force gives nan value, check parameters");
                }
                c->F[k] = fabs(c->Ftemp[k]);
                //c->F[k] = fabs(c->Ftemp[k])*c->Rho/newrho;
                //std::cout << newrho << std::endl;
            }
        }
        else
        {
            for (size_t j = 1;j<c->Nneigh;j++)
            {
                c->Ftemp[j] = c->F[j];
            }
            for (size_t j = 1;j<c->Nneigh;j++)
            {
                c->F[j]     = c->Ftemp[c->Op[j]];
            }
        }
    }
}

void Domain::CollideSC (size_t n, size_t Np)
{
    //std::cout << "SCP" << std::endl;
	size_t Ni = Lat[0].Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Lat[0].Ncells : Fn = (n+1)*Ni;
#ifdef USE_OMP
    In = 0;
    Fn = Lat[0].Ncells;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        Cell * c = Lat[0].Cells[i];
        double rho = c->Rho;
        //if (rho<1.0e-12) continue;
        //if (c->IsSolid||rho<1.0e-12) continue;
        //if (fabs(c->Gamma-1.0)<1.0e-12&&(fabs(Lat[j].G)>1.0e-12||Gmix>1.0e-12)) continue;
        //if (fabs(c->Gamma-1.0)<1.0e-12) continue;
        if (!c->IsSolid)
        {
            //Calculate Smagorinsky LES model
            Array<double> NonEq(c->Nneigh);
            //Array<double> FEq  (c->Nneigh);
            double Tau = Lat[0].Tau;
            Vec3_t DV  = c->Vel + c->BForce*Tau*dt/rho;
            //Vec3_t DV  = Vmix + c->BForce*Tau/rho;
            double Q = 0.0;
            for (size_t k=0;k<c->Nneigh;k++)
            {
                double FDeqn = c->Feq(k,DV,rho);
                NonEq[k] = c->F[k] - FDeqn;
                //FEq  [k] = FDeqn;
                Q += NonEq[k]*NonEq[k]*EEk[k];
                //if(c->Gamma>1.0e-12) c->Omeis[k] /= c->Gamma;
            }
            Q = sqrt(2.0*Q);
            Tau = 0.5*(Tau + sqrt(Tau*Tau + 6.0*Q*Sc/rho));
            c->Tau = Tau;
            //std::cout << Tau << std::endl;

            double Bn;
            //rho<10e-12 ? Bn =0.0 : Bn = (c->Gamma*(Lat[j].Tau-0.5))/((1.0-c->Gamma)+(Lat[j].Tau-0.5));
            //rho<10e-12 ? Bn =0.0 : Bn = c->Gamma;
            Bn = (c->Gamma*(Lat[0].Tau-0.5))/((1.0-c->Gamma)+(Lat[0].Tau-0.5));
            //Bn = (c->Gamma*(Tau-0.5))/((1.0-c->Gamma)+(Tau-0.5));
            //Bn = c->Gamma;
            //Bn = floor(c->Gamma);
            bool valid  = true;
            double alphal = 1.0;
            double alphat = 1.0;
            size_t num = 0;
            while (valid)
            {
                num++;
                valid = false;
                alphal  = alphat;
                for (size_t k=0;k<c->Nneigh;k++)
                {
                    //double FDeqn = c->Feq(k,DV,rho);
                    //c->Ftemp[k] = c->F[k] - alphal*((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]);
                    c->Ftemp[k] = c->F[k] - alphal*((1.0 - Bn)*(NonEq[k]/Tau - (1.0 - c->Pf)*(c->F[c->Op[k]] - c->F[k] + NonEq[k]/Tau)) - Bn*c->Omeis[k]);
                    if (c->Ftemp[k]<-1.0e-12&&num<2)
                    {
                        //double temp = fabs(c->F[k]/((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]));
                        double temp = fabs(c->F[k]/((1.0 - Bn)*(NonEq[k]/Tau - (1.0 - c->Pf)*(c->F[c->Op[k]] - c->F[k] + NonEq[k]/Tau)) - Bn*c->Omeis[k]));
                        if (temp<alphat) alphat = temp;
                        valid = true;
                    }
                }
            }
            for (size_t k=0;k<c->Nneigh;k++)
            {
                if (std::isnan(c->Ftemp[k]))
                {
                    //c->Gamma = 2.0;
                    #ifdef USE_HDF5
                    //WriteXDMF("error");
                    #endif
                    std::cout << "CollideSC: Nan found, resetting" << std::endl;
                    std::cout << c->Density() << " " << c->BForce << " " << num << " " << alphat << " " << c->Index << " " << c->IsSolid << " " << c->Gamma << " " << k << " " << c->Omeis[k] << std::endl;
                    //c->Ftemp[k] = c->F[k];
                    throw new Fatal("Domain::Collide: Body force gives nan value, check parameters");
                }
                c->F[k] = fabs(c->Ftemp[k]);
                //c->F[k] = fabs(c->Ftemp[k])*c->Rho/newrho;
                //std::cout << newrho << std::endl;
            }
        }
        else
        {
            for (size_t j = 1;j<c->Nneigh;j++)
            {
                c->Ftemp[j] = c->F[j];
            }
            for (size_t j = 1;j<c->Nneigh;j++)
            {
                c->F[j]     = c->Ftemp[c->Op[j]];
            }
        }
    }
}

void Domain::CollideMC (size_t n, size_t Np)
{
    //std::cout << "MCP" << std::endl;
	size_t Ni = Lat[0].Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Lat[0].Ncells : Fn = (n+1)*Ni;
#ifdef USE_OMP
    In = 0;
    Fn = Lat[0].Ncells;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        Vec3_t num(0.0,0.0,0.0);
        double den = 0.0;
        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double tau = Lat[j].Tau;
            num += c->Vel*c->Rho/tau;
            den += c->Rho/tau;
        }
        Vec3_t Vmix = num/den;


        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double rho = c->Rho;
            //if (rho<1.0e-12) continue;
            //if (c->IsSolid||rho<1.0e-12) continue;
            //if (fabs(c->Gamma-1.0)<1.0e-12&&(fabs(Lat[j].G)>1.0e-12||Gmix>1.0e-12)) continue;
            //if (fabs(c->Gamma-1.0)<1.0e-12) continue;
            if (!c->IsSolid)
            {
                //Array<double> FEq  (c->Nneigh);
                double Tau = Lat[j].Tau;
                Vec3_t DV  = Vmix + c->BForce*Tau*dt/rho;
                //Vec3_t DV  = Vmix + c->BForce*dt/rho;

                double Bn;
                //rho<10e-12 ? Bn =0.0 : Bn = (c->Gamma*(Lat[j].Tau-0.5))/((1.0-c->Gamma)+(Lat[j].Tau-0.5));
                //rho<10e-12 ? Bn =0.0 : Bn = c->Gamma;
                //Bn = (c->Gamma*(Lat[j].Tau-0.5))/((1.0-c->Gamma)+(Lat[j].Tau-0.5));
                Bn = floor(c->Gamma);
                bool valid  = true;
                double alphal = 1.0;
                double alphat = 1.0;
                size_t num = 0;
                while (valid)
                {
                    num++;
                    valid = false;
                    alphal  = alphat;
                    for (size_t k=0;k<c->Nneigh;k++)
                    {
                        if (num==1&&c->Omeis[k]>1.0e-12) c->Omeis[k] /= c->Gamma;
                        double FDeqn = c->Feq(k,DV,rho);
                        c->Ftemp[k] = c->F[k] - alphal*((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]);
                        //c->Ftemp[k] = c->F[k] - alphal*((1 - Bn)*(NonEq[k])/Tau - Bn*c->Omeis[k]);
                        //newrho += c->Ftemp[k];
                        if (c->Ftemp[k]<-1.0e-12&&num<2)
                        {
                            double temp = fabs(c->F[k]/((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]));
                            //double temp = fabs(c->F[k]/((1 - Bn)*(NonEq[k])/Tau - Bn*c->Omeis[k]));
                            if (temp<alphat) alphat = temp;
                            valid = true;
                        }
                    }
                }
                for (size_t k=0;k<c->Nneigh;k++)
                {
                    if (std::isnan(c->Ftemp[k]))
                    {
                        //c->Gamma = 2.0;
                        #ifdef USE_HDF5
                        //WriteXDMF("error");
                        #endif
                        //std::cout << "CollideMC: Nan found, resetting" << std::endl;
                        //std::cout << c->Density() << " " << c->BForce << " " << num << " " << alphat << " " << c->Index << " " << c->IsSolid << " " << j << " " << k << " " << std::endl;
                        c->Ftemp[k] = c->F[k];
                        //throw new Fatal("Domain::Collide: Body force gives nan value, check parameters");
                    }
                    c->F[k] = fabs(c->Ftemp[k]);
                    //c->F[k] = fabs(c->Ftemp[k])*c->Rho/newrho;
                    //std::cout << newrho << std::endl;
                }
            }
            else
            {
                for (size_t j = 1;j<c->Nneigh;j++)
                {
                    c->Ftemp[j] = c->F[j];
                }
                for (size_t j = 1;j<c->Nneigh;j++)
                {
                    c->F[j]     = c->Ftemp[c->Op[j]];
                }
            }
        }
    }   
}

void Domain::CollideNoPar (size_t n, size_t Np)
{
    //std::cout << "CNP" << std::endl;
	size_t Ni = Lat[0].Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Lat[0].Ncells : Fn = (n+1)*Ni;
#ifdef USE_OMP
    In = 0;
    Fn = Lat[0].Ncells;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        Vec3_t num(0.0,0.0,0.0);
        double den = 0.0;
        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double tau = Lat[j].Tau;
            num += c->Vel*c->Rho/tau;
            den += c->Rho/tau;
        }
        Vec3_t Vmix = num/den;


        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double rho = c->Rho;
            //if (rho<1.0e-12) continue;
            //if (c->IsSolid||rho<1.0e-12) continue;
            //if (fabs(c->Gamma-1.0)<1.0e-12&&(fabs(Lat[j].G)>1.0e-12||Gmix>1.0e-12)) continue;
            //if (fabs(c->Gamma-1.0)<1.0e-12) continue;
            if (!c->IsSolid)
            {
                //Calculate Smagorinsky LES model
                Array<double> NonEq(c->Nneigh);
                //Array<double> FEq  (c->Nneigh);
                double Tau = Lat[j].Tau;
                Vec3_t DV  = Vmix + c->BForce*Tau*dt/rho;
                //Vec3_t DV  = Vmix + c->BForce*Tau/rho;
                double Q = 0.0;
                for (size_t k=0;k<c->Nneigh;k++)
                {
                    double FDeqn = c->Feq(k,DV,rho);
                    NonEq[k] = c->F[k] - FDeqn;
                    //FEq  [k] = FDeqn;
                    Q += NonEq[k]*NonEq[k]*EEk[k];
                }
                Q = sqrt(2.0*Q);
                Tau = 0.5*(Tau + sqrt(Tau*Tau + 6.0*Q*Sc/rho));
                c->Tau = Tau;
                //std::cout << Tau << std::endl;

                bool valid  = true;
                double alphal = 1.0;
                double alphat = 1.0;
                size_t num = 0;
                while (valid)
                {
                    num++;
                    valid = false;
                    alphal  = alphat;
                    for (size_t k=0;k<c->Nneigh;k++)
                    {
                        c->Ftemp[k] = c->F[k] - alphal*((NonEq[k])/Tau - (1.0 - c->Pf)*(c->F[c->Op[k]] - c->F[k] + NonEq[k]/Tau));
                        //newrho += c->Ftemp[k];
                        if (c->Ftemp[k]<-1.0e-12&&num<2)
                        {
                            double temp = fabs(c->F[k]/(NonEq[k]/Tau - (1.0 - c->Pf)*(c->F[c->Op[k]] - c->F[k] + NonEq[k]/Tau)));
                            if (temp<alphat) alphat = temp;
                            valid = true;
                        }
                    }
                }
                for (size_t k=0;k<c->Nneigh;k++)
                {
                    if (std::isnan(c->Ftemp[k]))
                    {
                        //c->Gamma = 2.0;
                        #ifdef USE_HDF5
                        //WriteXDMF("error");
                        #endif
                        std::cout << "CollideNoPar: Nan found, resetting" << std::endl;
                        std::cout << c->Density() << " " << c->BForce << " " << num << " " << alphat << " " << c->Index << " " << c->IsSolid << " " << j << " " << k << " " << std::endl;
                        c->Ftemp[k] = c->F[k];
                        throw new Fatal("Domain::CollideNoPar: Body force gives nan value, check parameters");
                    }
                    c->F[k] = fabs(c->Ftemp[k]);
                    //c->F[k] = fabs(c->Ftemp[k])*c->Rho/newrho;
                    //std::cout << newrho << std::endl;
                }
            }
            else
            {
                for (size_t j = 1;j<c->Nneigh;j++)
                {
                    c->Ftemp[j] = c->F[j];
                }
                for (size_t j = 1;j<c->Nneigh;j++)
                {
                    c->F[j]     = c->Ftemp[c->Op[j]];
                }
            }
        }
    }   
}

void Domain::ImprintLatticeSC (size_t n,size_t Np)
{
    
	size_t Ni = ParCellPairs.Size()/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = ParCellPairs.Size() : Fn = (n+1)*Ni;
    //std::cout << "Im proccess = " << n << std::endl;
    // 2D imprint
    if (Lat[0].Ndim(2)==1)
    {
    #ifdef USE_OMP
        In = 0;
        Fn = ParCellPairs.Size();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
        for (size_t i = In;i<Fn;i++)
        {
            LBM::Disk      * Pa   = Disks[ParCellPairs[i].IPar];
            Cell           * cell = Lat[0].Cells[ParCellPairs[i].ICell];
            double x              = Lat[0].dx*(cell->Index(0));
            double y              = Lat[0].dx*(cell->Index(1));
            double z              = Lat[0].dx*(cell->Index(2));
            Vec3_t  C(x,y,z);
            double len = DEM::DiskSquare(Pa->X,C,Pa->R,Lat[0].dx);
            if (fabs(len)<1.0e-12) continue;
            double Tau = Lat[0].Tau;
            cell = Lat[0].Cells[ParCellPairs[i].ICell];
            double gamma  = len/(4.0*Lat[0].dx);
            cell->Gamma   = std::min(gamma+cell->Gamma,1.0);
            Vec3_t B      = C - Pa->X;
            Vec3_t tmp;
            Rotation(Pa->W,Pa->Q,tmp);
            Vec3_t VelP   = Pa->V + cross(tmp,B);
            double rho = cell->Rho;
            double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
            //double Bn  = gamma;
            //double Bn  = floor(gamma);
            size_t ncells = cell->Nneigh;
            Vec3_t Flbm = OrthoSys::O;
            for (size_t k=0;k<ncells;k++)
            {
                double Fvpp    = cell->Feq(cell->Op[k],VelP,rho);
                double Fvp     = cell->Feq(k          ,VelP,rho);
                double Omega   = cell->F[cell->Op[k]] - Fvpp - (cell->F[k] - Fvp);

                //cell->Omeis[k] += gamma*Omega;
                cell->Omeis[k] = Omega;
                Flbm += -Fconv*Bn*Omega*cell->Cs*cell->Cs*Lat[0].dx*Lat[0].dx*cell->C[k];
            }
            Vec3_t T,Tt;
            Tt =           cross(B,Flbm);
            Quaternion_t q;
            Conjugate    (Pa->Q,q);
            Rotation     (Tt,q,T);
            //std::cout << "1" << std::endl;
    #ifdef USE_OMP
            omp_set_lock      (&Pa->lck);
    #endif
            Pa->F          += Flbm;
            Pa->Flbm       += Flbm;
            Pa->T          += T;
    #ifdef USE_OMP
            omp_unset_lock    (&Pa->lck);
    #endif

        }
    }

    //3D imprint
    else
    {
    #ifdef USE_OMP
        In = 0;
        Fn = ParCellPairs.Size();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
        for (size_t i = In;i<Fn;i++)
        {
            DEM::Particle  * Pa   = Particles[ParCellPairs[i].IPar];
            Cell           * cell = Lat[0].Cells[ParCellPairs[i].ICell];
            double x              = Lat[0].dx*(cell->Index(0));
            double y              = Lat[0].dx*(cell->Index(1));
            double z              = Lat[0].dx*(cell->Index(2));
            Vec3_t  C(x,y,z);
            Vec3_t  Xtemp,Xs,Xstemp;
            double len,minl = Pa->Dmax;
    
            //std::cout << "1" << std::endl;
            //if (Pa->IsInsideFaceOnly(C)) len = 12.0*Lat[0].dx;
            //else if (ParCellPairs[i].IGeo.Size()==0) continue;
            if (norm(C-Pa->x)>Pa->Dmax) continue;
            len = 12.0*Lat[0].dx;
            Vec3_t Nor = OrthoSys::O;
            //else
            if (ParCellPairs[i].IGeo.Size()>0) 
            {
                if (Pa->Faces.Size()>0)
                {
                    DEM::Distance(C,*Pa->Faces[ParCellPairs[i].IGeo[0]],Xtemp,Xs);
                    minl = norm(Xtemp-Xs);
                    Nor = Pa->Faces[ParCellPairs[i].IGeo[0]]->Nor;
                    for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                    {
                        DEM::Distance(C,*Pa->Faces[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                        if (norm(Xtemp-Xstemp) < minl)
                        {
                            minl = norm(Xtemp-Xstemp);
                            Xs   = Xstemp;
                            Nor = Pa->Faces[ParCellPairs[i].IGeo[j]]->Nor;
                        }
                    }
                }
                else if (Pa->Edges.Size()>0)
                {
                    DEM::Distance(C,*Pa->Edges[ParCellPairs[i].IGeo[0]],Xtemp,Xs);
                    minl = norm(Xtemp-Xs);
                    for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                    {
                        DEM::Distance(C,*Pa->Edges[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                        if (norm(Xtemp-Xstemp) < minl)
                        {
                            minl = norm(Xtemp-Xstemp);
                            Xs   = Xstemp;
                        }
                    }
                }
                else if (Pa->Verts.Size()>0)
                {
                    DEM::Distance(C,*Pa->Verts[ParCellPairs[i].IGeo[0]],Xtemp,Xs);
                    minl = norm(Xtemp-Xs);
                    for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                    {
                        DEM::Distance(C,*Pa->Verts[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                        if (norm(Xtemp-Xstemp) < minl)
                        {
                            minl = norm(Xtemp-Xstemp);
                            Xs   = Xstemp;
                        }
                    }
                }
                double dotpro = dot(C-Xs,Nor);
                if (dotpro>0.0||fabs(dotpro)<0.95*minl||Pa->Faces.Size()<4) 
                {
                    len = DEM::SphereCube(Xs,C,Pa->Props.R,Lat[0].dx);
                }
            }
            //std::cout << "2" << std::endl;
            if (fabs(len)<1.0e-12) continue;
            double Tau = Lat[0].Tau;
            cell = Lat[0].Cells[ParCellPairs[i].ICell];
            double gamma  = len/(12.0*Lat[0].dx);
            cell->Gamma = gamma;
            //cell->Gamma   = std::max(gamma,cell->Gamma);
            //cell->Gamma   = std::min(gamma+cell->Gamma,1.0);
            //if (fabs(cell->Gamma-1.0)<1.0e-12)
            //if (fabs(cell->Gamma-1.0)<1.0e-12&&(fabs(Lat[0].G)>1.0e-12||Gmix>1.0e-12)) 
            Vec3_t B      = C - Pa->x;
            Vec3_t tmp;
            Rotation(Pa->w,Pa->Q,tmp);
            Vec3_t VelP   = Pa->v + cross(tmp,B);
            double rho = cell->Rho;
            double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
            //double Bn  = gamma;
            //double Bn  = floor(gamma);
            size_t ncells = cell->Nneigh;
            Vec3_t Flbm = OrthoSys::O;
            for (size_t k=0;k<ncells;k++)
            {
                double Fvpp     = cell->Feq(cell->Op[k],VelP,rho);
                double Fvp      = cell->Feq(k          ,VelP,rho);
                double Omega    = cell->F[cell->Op[k]] - Fvpp - (cell->F[k] - Fvp);
                //cell->Omeis[k] += Omega;
                //cell->Omeis[k] += gamma*Omega;
                cell->Omeis[k] = Omega;
                Flbm += -Fconv*Bn*Omega*cell->Cs*cell->Cs*Lat[0].dx*Lat[0].dx*cell->C[k];
            }
            Vec3_t T,Tt;
            Tt =           cross(B,Flbm);
            Quaternion_t q;
            Conjugate    (Pa->Q,q);
            Rotation     (Tt,q,T);
            //std::cout << "1" << std::endl;
    #ifdef USE_OMP
            omp_set_lock      (&Pa->lck);
    #endif
            Pa->F          += Flbm;
            Pa->Flbm       += Flbm;
            Pa->T          += T;
    #ifdef USE_OMP
            omp_unset_lock    (&Pa->lck);
    #endif
            //std::cout << "3" << std::endl;
        }
    }
}

void Domain::ImprintLatticeMC (size_t n,size_t Np)
{
    
	size_t Ni = ParCellPairs.Size()/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = ParCellPairs.Size() : Fn = (n+1)*Ni;
    //std::cout << "Im proccess = " << n << std::endl;
    // 2D imprint
    if (Lat[0].Ndim(2)==1)
    {
    #ifdef USE_OMP
        In = 0;
        Fn = ParCellPairs.Size();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
        for (size_t i = In;i<Fn;i++)
        {
            LBM::Disk      * Pa   = Disks[ParCellPairs[i].IPar];
            Cell           * cell = Lat[0].Cells[ParCellPairs[i].ICell];
            double x              = Lat[0].dx*(cell->Index(0));
            double y              = Lat[0].dx*(cell->Index(1));
            double z              = Lat[0].dx*(cell->Index(2));
            Vec3_t  C(x,y,z);
            double len = DEM::DiskSquare(Pa->X,C,Pa->R,Lat[0].dx);
            if (fabs(len)<1.0e-12) continue;
            for (size_t j=0;j<Lat.Size();j++)
            {
                cell = Lat[j].Cells[ParCellPairs[i].ICell];
                double gamma  = len/(4.0*Lat[0].dx);
                cell->Gamma   = std::min(gamma+cell->Gamma,1.0);
                //if (fabs(cell->Gamma-1.0)<1.0e-12) 
                Vec3_t B      = C - Pa->X;
                Vec3_t tmp;
                Rotation(Pa->W,Pa->Q,tmp);
                Vec3_t VelP   = Pa->V + cross(tmp,B);
                double rho = cell->Rho;
                //double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
                //double Bn  = gamma;
                double Bn  = floor(gamma);
                size_t ncells = cell->Nneigh;
                Vec3_t Flbm = OrthoSys::O;
                for (size_t k=0;k<ncells;k++)
                {
                    double Fvpp    = cell->Feq(cell->Op[k],VelP,rho);
                    double Fvp     = cell->Feq(k          ,VelP,rho);
                    double Omega   = cell->F[cell->Op[k]] - Fvpp - (cell->F[k] - Fvp);

                    cell->Omeis[k] += gamma*Omega;
                    Cell *nb        = Lat[j].Cells[cell->Neighs[k]];

                    double Fcontact;
                    (cell->Gammap>0.0)&&(cell->Gammap<1.0) ? Fcontact = Lat[j].Gs*cell->W[k]*cell->Rho*floor(nb->Gammap) : Fcontact = 0.0;
                    Flbm += -(Fconv*Bn*Omega*cell->Cs*cell->Cs*Lat[0].dx*Lat[0].dx - Fcontact)*cell->C[k];

                }
                Vec3_t T,Tt;
                Tt =           cross(B,Flbm);
                Quaternion_t q;
                Conjugate    (Pa->Q,q);
                Rotation     (Tt,q,T);
                //std::cout << "1" << std::endl;
    #ifdef USE_OMP
                omp_set_lock      (&Pa->lck);
    #endif
                Pa->F          += Flbm;
                Pa->Flbm       += Flbm;
                Pa->T          += T;
    #ifdef USE_OMP
                omp_unset_lock    (&Pa->lck);
    #endif
            }

        }
    }

    //3D imprint
    else
    {
    #ifdef USE_OMP
        In = 0;
        Fn = ParCellPairs.Size();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
        for (size_t i = In;i<Fn;i++)
        {
            DEM::Particle  * Pa   = Particles[ParCellPairs[i].IPar];
            Cell           * cell = Lat[0].Cells[ParCellPairs[i].ICell];
            double x              = Lat[0].dx*(cell->Index(0));
            double y              = Lat[0].dx*(cell->Index(1));
            double z              = Lat[0].dx*(cell->Index(2));
            Vec3_t  C(x,y,z);
            Vec3_t  Xtemp,Xs,Xstemp;
            double len,minl = Pa->Dmax;
    
            //std::cout << "1" << std::endl;
            //if (Pa->IsInsideFaceOnly(C)) len = 12.0*Lat[0].dx;
            //else if (ParCellPairs[i].IGeo.Size()==0) continue;
            if (norm(C-Pa->x)>Pa->Dmax) continue;
            len = 12.0*Lat[0].dx;
            Vec3_t Nor = OrthoSys::O;
            //else
            if (ParCellPairs[i].IGeo.Size()>0) 
            {
                if (Pa->Faces.Size()>0)
                {
                    DEM::Distance(C,*Pa->Faces[ParCellPairs[i].IGeo[0]],Xtemp,Xs);
                    minl = norm(Xtemp-Xs);
                    Nor = Pa->Faces[ParCellPairs[i].IGeo[0]]->Nor;
                    for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                    {
                        DEM::Distance(C,*Pa->Faces[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                        if (norm(Xtemp-Xstemp) < minl)
                        {
                            minl = norm(Xtemp-Xstemp);
                            Xs   = Xstemp;
                            Nor = Pa->Faces[ParCellPairs[i].IGeo[j]]->Nor;
                        }
                    }
                }
                else if (Pa->Edges.Size()>0)
                {
                    DEM::Distance(C,*Pa->Edges[ParCellPairs[i].IGeo[0]],Xtemp,Xs);
                    minl = norm(Xtemp-Xs);
                    for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                    {
                        DEM::Distance(C,*Pa->Edges[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                        if (norm(Xtemp-Xstemp) < minl)
                        {
                            minl = norm(Xtemp-Xstemp);
                            Xs   = Xstemp;
                        }
                    }
                }
                else if (Pa->Verts.Size()>0)
                {
                    DEM::Distance(C,*Pa->Verts[ParCellPairs[i].IGeo[0]],Xtemp,Xs);
                    minl = norm(Xtemp-Xs);
                    for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                    {
                        DEM::Distance(C,*Pa->Verts[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                        if (norm(Xtemp-Xstemp) < minl)
                        {
                            minl = norm(Xtemp-Xstemp);
                            Xs   = Xstemp;
                        }
                    }
                }
                double dotpro = dot(C-Xs,Nor);
                if (dotpro>0.0||fabs(dotpro)<0.95*minl||Pa->Faces.Size()<4) 
                {
                    len = DEM::SphereCube(Xs,C,Pa->Props.R,Lat[0].dx);
                }
            }
            //std::cout << "2" << std::endl;
            if (fabs(len)<1.0e-12) continue;
    
            for (size_t j=0;j<Lat.Size();j++)
            {
                cell = Lat[j].Cells[ParCellPairs[i].ICell];
                //cell->Gamma   = std::max(len/(12.0*Lat[0].dx),cell->Gamma);
                double gamma  = len/(12.0*Lat[0].dx);
                cell->Gamma   = std::min(gamma+cell->Gamma,1.0);
                //if (fabs(cell->Gamma-1.0)<1.0e-12)
                //if (fabs(cell->Gamma-1.0)<1.0e-12&&(fabs(Lat[0].G)>1.0e-12||Gmix>1.0e-12)) 
                Vec3_t B      = C - Pa->x;
                Vec3_t tmp;
                Rotation(Pa->w,Pa->Q,tmp);
                Vec3_t VelP   = Pa->v + cross(tmp,B);
                double rho = cell->Rho;
                double Bn  = (gamma*(Lat[j].Tau-0.5))/((1.0-gamma)+(Lat[j].Tau-0.5));
                //double Bn  = gamma;
                //double Bn  = floor(gamma);
                size_t ncells = cell->Nneigh;
                Vec3_t Flbm = OrthoSys::O;
                for (size_t k=0;k<ncells;k++)
                {
                    double Fvpp     = cell->Feq(cell->Op[k],VelP,rho);
                    double Fvp      = cell->Feq(k          ,VelP,rho);
                    double Omega    = cell->F[cell->Op[k]] - Fvpp - (cell->F[k] - Fvp);
                    //cell->Omeis[k] += Omega;
                    cell->Omeis[k] += gamma*Omega;
                    Cell *nb        = Lat[j].Cells[cell->Neighs[k]];

                    double Fcontact;
                    (cell->Gammap>0.0)&&(cell->Gammap<1.0) ? Fcontact = Lat[j].Gs*cell->W[k]*cell->Rho*floor(nb->Gammap) : Fcontact = 0.0;
                    Flbm += -(Fconv*Bn*Omega*cell->Cs*cell->Cs*Lat[0].dx*Lat[0].dx - Fcontact)*cell->C[k];
                }
                Vec3_t T,Tt;
                Tt =           cross(B,Flbm);
                Quaternion_t q;
                Conjugate    (Pa->Q,q);
                Rotation     (Tt,q,T);
                //std::cout << "1" << std::endl;
    #ifdef USE_OMP
                omp_set_lock      (&Pa->lck);
    #endif
                Pa->F          += Flbm;
                Pa->Flbm       += Flbm;
                Pa->T          += T;
    #ifdef USE_OMP
                omp_unset_lock    (&Pa->lck);
    #endif
                //std::cout << "3" << std::endl;
                
            }
        }
    }
}

inline void Domain::ResetDisplacements()
{
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].Dmx = 0.0;
        MTD[i].LLC.Resize(0);
    }
    //std::cout << "1" << std::endl;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Particles.Size();i++)
    {
        Particles[i]->ResetDisplacements();
        if(Particles[i]->IsFree())
        {
            iVec3_t idx = (Particles[i]->x - LCxmin)/(2.0*Beta*MaxDmax);
            MTD[omp_get_thread_num()].LLC.Push(std::make_pair(idx,i));
        }
    }
    //std::cout << "2" << std::endl;
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LLC.Size();j++)
        {
            //std::cout << i << " " << j << std::endl;
            size_t idx = DEM::Pt2idx(MTD[i].LLC[j].first,LCellDim);
            LinkedCell[idx].Push(MTD[i].LLC[j].second);
        }
    }

    //Only for 2D Disks
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Disks.Size();i++)
    {
        Disks[i]->X0 = Disks[i]->X;
    }
    //std::cout << "3" << std::endl;
#else
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->ResetDisplacements();
    }
#endif
}

inline double Domain::MaxDisplacement()
{
    double md = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        double mdp = Particles[i]->MaxDisplacement();
        if (mdp>md) md = mdp;
    }
    return md;
}

#ifdef USE_OMP
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
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t idx=0;idx<LinkedCell.Size();idx++)
    {
        if (LinkedCell[idx].Size()==0) continue;
        iVec3_t index;
        DEM::idx2Pt(idx,index,LCellDim);
        for (size_t n=0  ;n<LinkedCell[idx].Size()-1;n++)
        for (size_t m=n+1;m<LinkedCell[idx].Size()  ;m++)
        {
            size_t i1 = LinkedCell[idx][n];
            size_t i2 = LinkedCell[idx][m];
            if (i1==i2) continue;
            MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
        }
        size_t i = index(0);
        size_t j = index(1);
        size_t k = index(2);
        for (size_t knb=std::max(0,int(k)-1);knb<=std::min(LCellDim(2)-1,k+1);knb++)
        for (size_t jnb=std::max(0,int(j)-1);jnb<=std::min(LCellDim(1)-1,j+1);jnb++)
        for (size_t inb=std::max(0,int(i)-1);inb<=std::min(LCellDim(0)-1,i+1);inb++)
        {
            iVec3_t Ptnb(inb,jnb,knb);
            size_t idxnb = DEM::Pt2idx(Ptnb,LCellDim);
            if (idxnb>idx)
            {
                for (size_t n=0;n<LinkedCell[idx].Size()  ;n++)
                {
                    for (size_t m=0;m<LinkedCell[idxnb].Size()  ;m++)
                    {
                        size_t i1 = std::min(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        size_t i2 = std::max(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        if (i1==i2) continue;
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
#endif

inline void Domain::ResetContacts()
{
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].LC.Resize(0);
        MTD[i].LCI.Resize(0);
        MTD[i].LCB.Resize(0);
        MTD[i].LPC.Resize(0);
    }

    //2D case
    if (Lat[0].Ndim(2)==1)
    {
        //2D case
        #pragma omp parallel for schedule(static) num_threads(Nproc)
	    for (size_t i=0;i<Disks.Size();i++)
        {
            LBM::Disk * Pa = Disks[i];
            for (size_t n=std::max(0.0,double(Pa->X(0)-Pa->R-2.0*Alpha-Lat[0].dx)/Lat[0].dx);n<=std::min(double(Lat[0].Ndim(0)-1),double(Pa->X(0)+Pa->R+2.0*Alpha+Lat[0].dx)/Lat[0].dx);n++)
            for (size_t m=std::max(0.0,double(Pa->X(1)-Pa->R-2.0*Alpha-Lat[0].dx)/Lat[0].dx);m<=std::min(double(Lat[0].Ndim(1)-1),double(Pa->X(1)+Pa->R+2.0*Alpha+Lat[0].dx)/Lat[0].dx);m++)
            {
                Cell  * cell = Lat[0].GetCell(iVec3_t(n,m,0));
                double x     = Lat[0].dx*(cell->Index(0));
                double y     = Lat[0].dx*(cell->Index(1));
                double z     = Lat[0].dx*(cell->Index(2));
                Vec3_t  C(x,y,z);
                if ((norm(C-Pa->X)>2.0*Alpha+2.0*Lat[0].dx+Pa->R)) continue;
                ParticleCellPair NewPCP;
                NewPCP.IPar = i;
                NewPCP.ICell= cell->ID;

                MTD[omp_get_thread_num()].LPC.Push(NewPCP);
            }
        }
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t n=0;n<List2DPairs.Size();n++)
        {
            size_t i = List2DPairs[n].first;
            size_t j = List2DPairs[n].second;
            bool pi_has_vf = !Disks[i]->IsFree();
            bool pj_has_vf = !Disks[j]->IsFree();
            bool close = (DEM::Distance(Disks[i]->X,Disks[j]->X)<=Disks[i]->R+Disks[j]->R+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close) continue;
            set<pair<Disk *, Disk *> >::iterator it = ListofDiskPairs.find(make_pair(Disks[i],Disks[j]));
            if (it != ListofDiskPairs.end())
            {
                continue;
            }
            MTD[omp_get_thread_num()].LC.Push(std::make_pair(i,j));
        }
        for (size_t i=0;i<Nproc;i++)
        {
            //std::cout << MTD[i].LC.Size() << std::endl;
            for (size_t j=0;j<MTD[i].LC.Size();j++)
            {
            //std::cout << MTD[i].LC.Size() << std::endl;
                size_t n = MTD[i].LC[j].first;
                size_t m = MTD[i].LC[j].second;
                ListofDiskPairs.insert(std::make_pair(Disks[n],Disks[m]));
                DiskPairs.Push(new DiskPair(Disks[n],Disks[m]));
            }
        }
    }

    else
    {
        //std::cout << "1" << std::endl;
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t n=0;n<ListPosPairs.Size();n++)
        {
            size_t i = ListPosPairs[n].first;
            size_t j = ListPosPairs[n].second;
            bool pi_has_vf = !Particles[i]->IsFree();
            bool pj_has_vf = !Particles[j]->IsFree();
            bool close = (DEM::Distance(Particles[i]->x,Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close) continue;
            std::set<std::pair<DEM::Particle *, DEM::Particle *> >::iterator it = Listofpairs.find(std::make_pair(Particles[i],Particles[j]));
            if (it != Listofpairs.end())
            {
                continue;
            }
            MTD[omp_get_thread_num()].LC.Push(std::make_pair(i,j));
        }
        //std::cout << "2" << std::endl;
        for (size_t i=0;i<Nproc;i++)
        {
            //std::cout << MTD[i].LC.Size() << std::endl;
            for (size_t j=0;j<MTD[i].LC.Size();j++)
            {
            //std::cout << MTD[i].LC.Size() << std::endl;
                size_t n = MTD[i].LC[j].first;
                size_t m = MTD[i].LC[j].second;
                Listofpairs.insert(std::make_pair(Particles[n],Particles[m]));
                if (Particles[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
                {
                    CInteractons.Push (new DEM::CInteractonSphere(Particles[n],Particles[m]));
                }
                else
                {
                    CInteractons.Push (new DEM::CInteracton(Particles[n],Particles[m]));
                }
            }
        }
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t n=0;n<CInteractons.Size();n++)
        {
            if(CInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LCI.Push(n);
        }
        //std::cout << "3" << std::endl;
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t n=0;n<BInteractons.Size();n++)
        {
            if(BInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LCB.Push(n);
        }
        //std::cout << "4" << std::endl;
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

        //std::cout << "5" << std::endl;
        #pragma omp parallel for schedule(static) num_threads(Nproc)
	    for (size_t i=0;i<Particles.Size();i++)
        {
            //Cell  * cell = dat.Dom->Lat[0].Cells[i];
            DEM::Particle * Pa = Particles[i];
            if (Pa->Bdry) continue;
            //std::cout << std::max(0.0,double(Pa->x(0)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " 
                      //<< std::max(0.0,double(Pa->x(1)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " 
                      //<< std::max(0.0,double(Pa->x(2)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " << std::endl;
            //std::cout << std::min(double(dat.Dom->Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " "
                      //<< std::min(double(dat.Dom->Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " "
                      //<< std::min(double(dat.Dom->Lat[0].Ndim(2)-1),double(Pa->x(2)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " << std::endl;
            for (size_t n=std::max(0.0,double(Pa->x(0)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);n<=std::min(double(Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);n++)
            for (size_t m=std::max(0.0,double(Pa->x(1)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);m<=std::min(double(Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);m++)
            for (size_t l=std::max(0.0,double(Pa->x(2)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);l<=std::min(double(Lat[0].Ndim(2)-1),double(Pa->x(2)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);l++)
            //for (size_t n=0;n<dat.Dom->Particles.Size();n++)
            {
                //DEM::Particle * Pa = dat.Dom->Particles[n];
                Cell  * cell = Lat[0].GetCell(iVec3_t(n,m,l));
                double x     = Lat[0].dx*(cell->Index(0));
                double y     = Lat[0].dx*(cell->Index(1));
                double z     = Lat[0].dx*(cell->Index(2));
                Vec3_t  C(x,y,z);
                if ((norm(C-Pa->x)>2.0*Alpha+2.0*Lat[0].dx+Pa->Dmax)||cell->IsSolid) continue;
                ParticleCellPair NewPCP;
                NewPCP.IPar = i;
                //NewPCP.IPar = n;
                NewPCP.ICell= cell->ID;
                
                bool valid = true;
                Vec3_t Nor = OrthoSys::O;
                if (Pa->Faces.Size()>0)
                {
                    double minl = 2.0*Pa->Dmax;
                    Vec3_t Xs = OrthoSys::O;
                    for (size_t j=0;j<Pa->Faces.Size();j++)
                    {
                        Vec3_t Xstemp,Xtemp;
                        DEM::Distance(C,*Pa->Faces[j],Xtemp,Xstemp);
                        double dist = norm(Xtemp - Xstemp);
                        if (dist<minl)
                        {
                            Nor  = Pa->Faces[j]->Nor;
                            minl = dist;
                            Xs   = Xstemp;
                        }
                        if (dist<2.0*Alpha+Pa->Props.R)
                        {
                            if (Pa->Faces[j]->Area()<2.0*M_PI*Pa->Props.R*Pa->Props.R)
                            {
                                continue;
                            }
                            NewPCP.IGeo.Push(j);
                        }
                    }
                    double dotpro = dot(C-Xs,Nor);
                    if ((dotpro>0.0||fabs(dotpro)<0.95*minl)&&Pa->Faces.Size()>3&&minl>2.0*Alpha+Pa->Props.R) valid = false;
                    else if (minl>2.0*Alpha+Pa->Props.R&&Pa->Faces.Size()<4)               valid = false;
                    else if (Pa->Faces.Size()>3&&NewPCP.IGeo.Size()==0&&!Pa->IsInsideFaceOnly(C)) valid = false;
                }
                else if (Pa->Edges.Size()>0)
                {
                    for (size_t j=0;j<Pa->Edges.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Edges[j])<2.0*Alpha+Pa->Props.R) 
                        {
                            NewPCP.IGeo.Push(j);
                        }
                        else valid = false;
                    }
                }
                else if (Pa->Verts.Size()>0)
                {
                    for (size_t j=0;j<Pa->Verts.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Verts[j])<2.0*Alpha+Pa->Props.R)
                        {
                            NewPCP.IGeo.Push(j);
                        }
                        else valid = false;
                    }
                }
                if (valid) MTD[omp_get_thread_num()].LPC.Push(NewPCP);
                //bool valid = false;
                //if (Pa->Faces.Size()>0)
                //{
                    //for (size_t j=0;j<Pa->Faces.Size();j++)
                    //{
                        //if (DEM::Distance(C,*Pa->Faces[j])<2.0*Alpha+Pa->Props.R)
                        //{
                            //if (Pa->Faces[j]->Area()<2.0*M_PI*Pa->Props.R*Pa->Props.R)
                            //{
                                //continue;
                            //}
                            //NewPCP.IGeo.Push(j);
                            //valid = true;
                        //}
                    //}
                //}
                //else if (Pa->Edges.Size()>0)
                //{
                    //for (size_t j=0;j<Pa->Edges.Size();j++)
                    //{
                        //if (DEM::Distance(C,*Pa->Edges[j])<2.0*Alpha+Pa->Props.R) 
                        //{
                            //NewPCP.IGeo.Push(j);
                            //valid = true;
                        //}
                    //}
                //}
                //else if (Pa->Verts.Size()>0)
                //{
                    //for (size_t j=0;j<Pa->Verts.Size();j++)
                    //{
                        //if (DEM::Distance(C,*Pa->Verts[j])<2.0*Alpha+Pa->Props.R)
                        //{
                            //NewPCP.IGeo.Push(j);
                            //valid = true;
                        //}
                    //}
                //}
                //if (Pa->IsInsideFaceOnly(C)) valid = true;
                //if (valid) MTD[omp_get_thread_num()].LPC.Push(NewPCP);
            }
        }
    }
    //std::cout << "6" << std::endl;

    ParCellPairs.Resize(0);
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LPC.Size();j++)
        {
            ParCellPairs.Push(MTD[i].LPC[j]);
        }
    }
#else
	if (Particles.Size()==0) return;
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            bool pj_has_vf = !Particles[j]->IsFree();

            bool close = (norm(Particles[i]->x-Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close) continue;
            
            // checking if the interacton exist for that pair of particles
            set<std::pair<DEM::Particle *, DEM::Particle *> >::iterator it = Listofpairs.find(std::make_pair(Particles[i],Particles[j]));
            if (it != Listofpairs.end())
            {
                continue;
            }
            Listofpairs.insert(std::make_pair(Particles[i],Particles[j]));
            if (Particles[i]->Verts.Size()==1 && Particles[j]->Verts.Size()==1)
            {
                CInteractons.Push (new DEM::CInteractonSphere(Particles[i],Particles[j]));
            }
            else
            {
                CInteractons.Push (new DEM::CInteracton(Particles[i],Particles[j]));
            }
        }
    }

    Interactons.Resize(0);
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        if(CInteractons[i]->UpdateContacts(Alpha)) Interactons.Push(CInteractons[i]);
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        if(BInteractons[i]->UpdateContacts(Alpha)) Interactons.Push(BInteractons[i]);
    }

    ParCellPairs.Resize(0);

    for (size_t i=0; i<Particles.Size(); i++)
    {
        DEM::Particle * Pa = Particles[i];
        for (size_t n=std::max(0.0,double(Pa->x(0)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);n<=std::min(double(Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);n++)
        for (size_t m=std::max(0.0,double(Pa->x(1)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);m<=std::min(double(Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);m++)
        for (size_t l=std::max(0.0,double(Pa->x(2)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);l<=std::min(double(Lat[0].Ndim(2)-1),double(Pa->x(2)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);l++)
        {
            Cell  * cell = Lat[0].GetCell(iVec3_t(n,m,l));
            double x     = Lat[0].dx*(cell->Index(0));
            double y     = Lat[0].dx*(cell->Index(1));
            double z     = Lat[0].dx*(cell->Index(2));
            Vec3_t  C(x,y,z);
            ParticleCellPair NewPCP;
            NewPCP.IPar = i;
            NewPCP.ICell= cell->ID;
            bool valid = false;
            if (Pa->Faces.Size()>0)
            {
                for (size_t j=0;j<Pa->Faces.Size();j++)
                {
                    if (DEM::Distance(C,*Pa->Faces[j])<2.0*Alpha+Pa->Props.R)                        
                    {
                        if (Pa->Faces[j]->Area()<2.0*M_PI*Pa->Props.R*Pa->Props.R)
                        {
                            //std::cout << "it did" << std::endl;
                            continue;
                        }
                        NewPCP.IGeo.Push(j);
                        valid = true;
                    }
                }
            }
            else if (Pa->Edges.Size()>0)
            {
                for (size_t j=0;j<Pa->Edges.Size();j++)
                {
                    if (DEM::Distance(C,*Pa->Edges[j])<2.0*Alpha+Pa->Props.R)
                    {
                        NewPCP.IGeo.Push(j);
                        valid = true;
                    }
                }
            }
            else if (Pa->Verts.Size()>0)
            {
                for (size_t j=0;j<Pa->Verts.Size();j++)
                {
                    if (DEM::Distance(C,*Pa->Verts[j])<2.0*Alpha+Pa->Props.R)
                    {
                        NewPCP.IGeo.Push(j);
                        valid = true;
                    }
                }
            }
            if (Pa->IsInsideFaceOnly(C)) valid = true;
            if (valid) ParCellPairs.Push(NewPCP);
        }
    }
#endif
}

//Utility methods
inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    if (Particles.Size()==0&&Disks.Size()==0)
    {
        minX = OrthoSys::O;
        maxX = OrthoSys::O;
    }
    else if (Particles.Size()>0)
    {
        minX = Vec3_t(Particles[0]->MinX(), Particles[0]->MinY(), Particles[0]->MinZ());
        maxX = Vec3_t(Particles[0]->MaxX(), Particles[0]->MaxY(), Particles[0]->MaxZ());
        for (size_t i=1; i<Particles.Size(); i++)
        {
            if (minX(0)>Particles[i]->MinX()&&Particles[i]->IsFree()) minX(0) = Particles[i]->MinX();
            if (minX(1)>Particles[i]->MinY()&&Particles[i]->IsFree()) minX(1) = Particles[i]->MinY();
            if (minX(2)>Particles[i]->MinZ()&&Particles[i]->IsFree()) minX(2) = Particles[i]->MinZ();
            if (maxX(0)<Particles[i]->MaxX()&&Particles[i]->IsFree()) maxX(0) = Particles[i]->MaxX();
            if (maxX(1)<Particles[i]->MaxY()&&Particles[i]->IsFree()) maxX(1) = Particles[i]->MaxY();
            if (maxX(2)<Particles[i]->MaxZ()&&Particles[i]->IsFree()) maxX(2) = Particles[i]->MaxZ();
        }
    }
    else
    {
        minX = Vec3_t(Disks[0]->X(0)-Disks[0]->R, Disks[0]->X(1)-Disks[0]->R, Disks[0]->X(2)-Disks[0]->R);
        maxX = Vec3_t(Disks[0]->X(0)+Disks[0]->R, Disks[0]->X(1)+Disks[0]->R, Disks[0]->X(2)+Disks[0]->R);
        for (size_t i=1; i<Disks.Size(); i++)
        {
            if ((minX(0)>Disks[i]->X(0)-Disks[i]->R)&&Disks[i]->IsFree()) minX(0) = Disks[i]->X(0)-Disks[i]->R;
            if ((minX(1)>Disks[i]->X(1)-Disks[i]->R)&&Disks[i]->IsFree()) minX(1) = Disks[i]->X(1)-Disks[i]->R;
            if ((minX(2)>Disks[i]->X(2)-Disks[i]->R)&&Disks[i]->IsFree()) minX(2) = Disks[i]->X(2)-Disks[i]->R;
            if ((maxX(0)<Disks[i]->X(0)+Disks[i]->R)&&Disks[i]->IsFree()) maxX(0) = Disks[i]->X(0)+Disks[i]->R;
            if ((maxX(1)<Disks[i]->X(1)+Disks[i]->R)&&Disks[i]->IsFree()) maxX(1) = Disks[i]->X(1)+Disks[i]->R;
            if ((maxX(2)<Disks[i]->X(2)+Disks[i]->R)&&Disks[i]->IsFree()) maxX(2) = Disks[i]->X(2)+Disks[i]->R;
        }                                                                                                
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
}

//DEM particle methods
#include <mechsys/lbm/Dompargen.h>


//Dynamic methods
inline void Domain::Initialize (double dt)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing particles ------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    // set flag
    if (!Initialized)
    {
        Initialized = true;
        // initialize all particles
        Voltot = 0.0;
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->Initialize(i);
            Particles[i]->InitializeVelocity(dt);
            Voltot+=Particles[i]->Props.V;
        }
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

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{

    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Finished = false;
    Nproc = TheNproc;

    if (dtdem<1.0e-12||dtdem>dt) dtdem = dt;

    // initialize particles
    Initialize (dtdem);
    // calc the total volume of particles (solids)
    FreePar.Resize(0);
    NoFreePar.Resize(0);
    Vs = 0.0;
    Ms = 0.0;
    MaxDmax        =  0.0;
    double MaxKn   = -1.0;
    double MaxBn   = -1.0;
    double MinDmax = -1.0;
    double MinMass = -1.0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        if (Particles[i]->IsFree())
        {
            Vs += Particles[i]->Props.V;
            Ms += Particles[i]->Props.m;
            if (Particles[i]->Dmax     > MaxDmax) MaxDmax = Particles[i]->Dmax;
            if (Particles[i]->Props.Kn > MaxKn  ||(MinDmax<0.0)) MaxKn   = Particles[i]->Props.Kn;
            if (Particles[i]->Dmax     < MinDmax||(MinDmax<0.0)) MinDmax = Particles[i]->Dmax;
            if (Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = Particles[i]->Props.m;
            FreePar.Push(i);
        }
        else NoFreePar.Push(i);
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        double pbn = BInteractons[i]->Bn/BInteractons[i]->L0;
        if (pbn > MaxBn||(MinDmax<0.0)) MaxBn = pbn;
    }


    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    printf("%s  Porosity                         =  %g%s\n"       ,TERM_CLR4, 1.0 - Lat[0].SolidFraction()         , TERM_RST);
    printf("%s  Lattice velocity C = dx/dt       =  %g%s\n"       ,TERM_CLR4, Lat[0].dx/Lat[0].dt                  , TERM_RST);
    printf("%s  Total mass   of free particles   =  %g%s\n"       ,TERM_CLR4, Ms                                   , TERM_RST);
    printf("%s  Total volume of free particles   =  %g%s\n"       ,TERM_CLR4, Vs                                   , TERM_RST);
    printf("%s  Total number of particles        =  %zd%s\n"      ,TERM_CLR2, Particles.Size()+Disks.Size()        , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                   , TERM_RST);
    printf("%s  Time step for DEM                =  %g%s\n"       ,TERM_CLR2, dtdem                                , TERM_RST);
    printf("%s  Verlet distance                  =  %g%s\n"       ,TERM_CLR2, Alpha                                , TERM_RST);
    for (size_t i=0;i<Lat.Size();i++)
    {
    printf("%s  Tau of Lattice %zd                 =  %g%s\n"       ,TERM_CLR2, i, Lat[i].Tau                        , TERM_RST);
    }

    if (FreePar.Size()>0)
    {
    printf("%s  Suggested Time Step for DEM      =  %g%s\n"       ,TERM_CLR5, 0.1*sqrt(MinMass/(MaxKn+MaxBn))      , TERM_RST);
    printf("%s  Suggested Verlet distance        =  %g or %g%s\n" ,TERM_CLR5, 0.5*MinDmax, 0.25*(MinDmax + MaxDmax), TERM_RST);

    if (Alpha > MinDmax&&MinDmax>0.0)
    {
        Alpha = MinDmax;
        printf("%s  Verlet distance changed to       =  %g%s\n"   ,TERM_CLR2, Alpha                                    , TERM_RST);
    }

    if (Alpha < 0.0) throw new Fatal("Verlet distance cannot be negative");
    }

     


    //std::cout << "1" << std::endl;
    // Creates pair of cells to speed up body force calculation
    for (size_t i=0;i<Lat[0].Ncells;i++)
    {
        Cell * c = Lat[0].Cells[i];
        if (fabs(1.0-c->Pf) > 1.0e-8) PrtPer = true;
        for (size_t j=1;j<c->Nneigh;j++)
        {
            Cell * nb = Lat[0].Cells[c->Neighs[j]];
            if (nb->ID>c->ID) 
            {
                //if (!c->IsSolid||!nb->IsSolid) CellPairs.Push(iVec3_t(i,nb->ID,j));
                CellPairs.Push(iVec3_t(i,nb->ID,j));
            }
        }
    }
    //std::cout << "2" << std::endl;
    
    MTD = new LBM::MtData[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = this;
        MTD[i].Dmx      = 0.0;
        MTD[i].dt       = Lat[0].dt;
    }
    //std::cout << "3" << std::endl;
#ifdef USE_OMP
    if (Disks.Size()>0)
    {
        List2DPairs.Resize(0);
        for (size_t i=0;i<Disks.Size()-1;i++)
        for (size_t j=i+1;j<Disks.Size();j++)
        {
            List2DPairs.Push(make_pair(i,j));
        }
    }

    LinkedCell.Resize(0);
    BoundingBox(LCxmin,LCxmax);
    LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax) + iVec3_t(1,1,1);
    LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));
    //std::cout << "1" << std::endl;
    ResetDisplacements();

    //std::cout << "2" << std::endl;
    UpdateLinkedCells();

    //std::cout << "3" << std::endl;
    ResetContacts();

    //std::cout << "4" << std::endl;
    ImprintLatticeSC(0,Nproc);    
#else

    //No longer used, OpenMp should be used by default

#endif
    double tout = Time;
    double tlbm = Time;

    //std::cout << "4" << std::endl;
    while (Time < Tf)
    {
        //std::cout << Interactons.Size() << " " << CInteractons.Size() << " " << BInteractons.Size() << " " << ParCellPairs.Size() << " " << Particles.Size() << std::endl;
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time >= tout)
        {
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    if (PrtDou) WriteXDMF_D(fn.CStr());
                    else        WriteXDMF  (fn.CStr());
                    #else
                    //WriteVTK (fn.CStr());
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            if (BInteractons.Size()>0) Clusters();
            tout += dtOut;
            idx_out++;
        }


#ifdef USE_OMP 
        //std::chrono::high_resolution_clock::time_point ti1 = std::chrono::high_resolution_clock::now();
        //Initialize all the particles and cells
        //std::cout << "0" <<std::endl;
        for (size_t i=0;i<Lat.Size();i++)
        {
            Lat[i].SetZeroGamma(0,Nproc);
        }
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for(size_t i=0;i<Particles.Size();i++)
        {
            Particles[i]->F    = Particles[i]->Ff;
            Particles[i]->Flbm = OrthoSys::O;
            Particles[i]->T    = Particles[i]->Tf;
        }
        //2D case
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for(size_t i=0;i<Disks.Size();i++)
        {
            Disks[i]->F    = Disks[i]->Ff;
            Disks[i]->Flbm = OrthoSys::O;
            Disks[i]->T    = Disks[i]->Tf;
        }

        //std::cout << "1" <<std::endl;
        
        //Imprint the particles into the lattice
        if (Particles.Size()>0||Disks.Size()>0)
        {
            if (Lat.Size()>1||fabs(Lat[0].Gs)>0.0)
            {
                ImprintLatticeMC(0,Nproc);
            }
            else
            {
                ImprintLatticeSC(0,Nproc);
            }
        }

        //std::chrono::high_resolution_clock::time_point ti2 = std::chrono::high_resolution_clock::now();
//
        //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( ti2 - ti1 ).count();
//
        //std::cout << "Duration of DEM/LBM part " << duration << std::endl;
        //
        //std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //std::cout << "2" <<std::endl;
        //Calculate interparticle forces
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<Interactons.Size(); i++)
        {
		    if (Interactons[i]->CalcForce(dtdem))
            {
                WriteXDMF("error");
                std::cout << "Maximun overlap detected between particles at time " << Time << std::endl;
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

        //std::cout << "3" <<std::endl;
        //2D case
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<DiskPairs.Size(); i++)
        {
		    DiskPairs[i]->CalcForce(dtdem);
            omp_set_lock  (&DiskPairs[i]->P1->lck);
            DiskPairs[i]->P1->F += DiskPairs[i]->F1;
            DiskPairs[i]->P1->T += DiskPairs[i]->T1;
            omp_unset_lock(&DiskPairs[i]->P1->lck);
            omp_set_lock  (&DiskPairs[i]->P2->lck);
            DiskPairs[i]->P2->F += DiskPairs[i]->F2;
            DiskPairs[i]->P2->T += DiskPairs[i]->T2;
            omp_unset_lock(&DiskPairs[i]->P2->lck);
        }
        //if (DiskPairs.Size()>0)
        //if (norm(DiskPairs[0]->F1)>1.0e10) 
        //{
            //std::cout << DiskPairs.Size() << std::endl;
            //std::cout << DiskPairs[0]->F1(0) << std::endl;
            //std::cout << Disks[0]->W << " " << Disks[1]->W << std::endl;
            //std::cout << Disks[0]->V << " " << Disks[1]->V << std::endl;
            //std::cout << Disks[0]->X << " " << Disks[1]->X << std::endl;
            //throw new Fatal("Bad");
        //}
        //if (isnan(norm(Disks[1]->F))) std::cout << Disks[0]->W << " " << Disks[1]->W << std::endl;
        //Checking if particles have moved beyond the verlet distance
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
		        Particles[i]->Translate(dtdem);
		        Particles[i]->Rotate(dtdem);
                if (Particles[i]->MaxDisplacement()>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = Particles[i]->MaxDisplacement();
            }
        }
        else
        {
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<Particles.Size(); i++)
            {
		        Particles[i]->Translate(dtdem);
                if (Particles[i]->MaxDisplacement()>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = Particles[i]->MaxDisplacement();
            }

        }

        //2D case

        if (RotPar)
        {
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<Disks.Size(); i++)
            {
		        Disks[i]->Translate(dtdem);
		        Disks[i]->Rotate(dtdem);
                if (norm(Disks[i]->X0-Disks[i]->X)>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = norm(Disks[i]->X0-Disks[i]->X);
            }
        }
        else
        {
            #pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i=0; i<Disks.Size(); i++)
            {
		        Disks[i]->Translate(dtdem);
                if (norm(Disks[i]->X0-Disks[i]->X)>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = norm(Disks[i]->X0-Disks[i]->X);
            }
        }

        double maxdis = 0.0;
        for (size_t i=0;i<Nproc;i++)
        {
            if (maxdis<MTD[i].Dmx) maxdis = MTD[i].Dmx;
        }
         
        //std::cout << "4 " << maxdis << std::endl;
        if (maxdis>Alpha)
        {
            LinkedCell.Resize(0);
            BoundingBox(LCxmin,LCxmax);
            LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax) + iVec3_t(1,1,1);
            LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));

            //std::cout << "4a" <<std::endl;
            ResetDisplacements();

            //std::cout << "4b" <<std::endl;
            UpdateLinkedCells();

            //std::cout << "4c" <<std::endl;
            ResetContacts();

        }
        //std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        //auto duration12 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

        //std::cout << "Duration of DEM part " << duration12 << std::endl;
        
        //std::chrono::high_resolution_clock::time_point t3 = std::chrono::high_resolution_clock::now();
        if (Time>=tlbm)
        {
            //Apply molecular forces
            if (Lat.Size()>1||(fabs(Lat[0].G)+fabs(Lat[0].Gs)>1.0e-12))
            {
                bool MC = false;
                if (Lat.Size()==2)
                {
                    if (fabs(Lat[0].G)<1.0e-9&&fabs(Lat[1].G)<1.0e-9) MC = true;
                }
                ApplyForce(0,Nproc,MC);
            }

            //Apply collision operator
            if (Particles.Size()>0||Disks.Size()>0)
            {
                if (Lat.Size()>1)
                {
                    CollideMC(0,Nproc);
                }
                else
                {
                    CollideSC(0,Nproc);
                    //CollideMRT(0,Nproc);
                }
            }
            else
            {
                CollideNoPar(0,Nproc);
                //CollideMRT(0,Nproc);
            }

            //Stream the distribution functions
            for (size_t i=0;i<Lat.Size();i++)
            {
                Lat[i].Stream(0,Nproc);
            }
            tlbm += dt;
        }
        
        //std::chrono::high_resolution_clock::time_point t4 = std::chrono::high_resolution_clock::now();
        //auto duration43 = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
//
        //std::cout << "Duration of LBM part " << duration43 << std::endl;
        //std::cout << "5" <<std::endl;
#else
        //No longer used, OpenMp should be used by default

#endif

        Time += dtdem;
        //std::cout << Time << " " << tlbm << std::endl;
    }
    // last output
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

//Matrices definitions

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
