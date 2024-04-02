/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2017 Maziar Gholami                                    *
 * Copyright (C) 2020 Mario Trujillo                                    *
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

#ifndef SPH_DOMAIN_H
#define SPH_DOMAIN_H

#include <stdio.h>    // for NULL
#include <algorithm>  // for min,max

#include <hdf5.h>
#include <hdf5_hl.h>

#include <omp.h>

#include <mechsys/sph/Particle.h>
#include <mechsys/sph/Functions.h>
#include <mechsys/sph/Boundary_Condition.h>

//Include DEM functions
#include <mechsys/dem/interacton.h>


//**********************CPU TIME COUNTER********************
#include <chrono>
using namespace std::chrono; 
//**********************CPU TIME COUNTER********************

//***************************COLORS*************************
#ifndef _COLORS_
#define _COLORS_
/* FOREGROUND */
#define RST   "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST
#endif  
/* _COLORS_ */
#include <iostream>
//***************************COLORS*************************


//C++ Enum used for easiness of coding in the input files
enum Kernels_Type { Qubic_Spline=0, Quintic=1, Quintic_Spline=2 };
enum Viscosity_Eq_Type { Morris=0, Shao=1, Incompressible_Full=2, Takeda=3 };
enum Gradient_Type { Squared_density=0, Multiplied_density=1 };    

namespace SPH {

class Domain
{
public:
	typedef void (*PtVel) (Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry);
	typedef void (*PtOut) (Particle * Particles, double & Prop1, double & Prop2,  double & Prop3);
	typedef void (*PtDom) (Domain & dom);
    // Constructor
    Domain();

    // Destructor
    ~Domain();

    // Domain Part
    
    // Particle generation
    void AddSingleParticle	(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed);		//Add one particle
    void AddBoxLength				(int tag, Vec3_t const &V, double Lx, double Ly, double Lz,double r, double Density,
																	double h,int type, int rotation, bool random, bool Fixed);									//Add a cube of particles with a defined dimensions
    void AddBoxNo						(int tag, Vec3_t const &V, size_t nx, size_t ny, size_t nz,double r, double Density,
																	double h,int type, int rotation, bool random, bool Fixed);									//Add a cube of particles with a defined numbers

    //DEM Particle Generation
    void AddSphere      (int Tag,Vec3_t const & X, double R, double rho);  ///< Add a DEM sphere or Disk
    void AddSegment     (int Tag, Vec3_t const & X0, Vec3_t const & X1, double R, double rho);            ///< Add a rice with limits at X0 and X1 with spheroradius R, side of length L and density rho
    void AddPlane       (int Tag, Vec3_t const & X, double R, double Lx,double Ly, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho
        
    void DelParticles				(int const & Tags);					//Delete particles by tag
    void CheckParticleLeave	();													//Check if any particles leave the domain, they will be deleted

    void YZPlaneCellsNeighbourSearch(int q1);						//Create pairs of particles in cells of XZ plan
    void MainNeighbourSearch				();									//Create pairs of particles in the whole domain
    void StartAcceleration					(Vec3_t const & a = Vec3_t(0.0,0.0,0.0));	//Add a fixed acceleration such as the Gravity
    void PrimaryComputeAcceleration	();									//Compute the solid boundary properties
    void LastComputeAcceleration		();								//Compute the acceleration due to the other particles
    void CalculateForceDEM		();									    //Compute the acceleration due to DEM Particles
    void CalcForce11		(Particle * P1, Particle * P2);	//Calculates the contact force between fluid-fluid particles
    void CalcForce2233	    (Particle * P1, Particle * P2);	//Calculates the contact force between soil-soil/solid-solid particles
    void CalcForce12		(Particle * P1, Particle * P2);	//Calculates the contact force between fluid-solid particles
    void CalcForce13		(Particle * P1, Particle * P2);	//Calculates the contact force between fluid-soil particles
    void KGC                (Particle * P1, Particle * P2); //Kernel Gradient Correction
    void HomoAcc            ();                             //Homogenize the accelaration for fixed groups
    void Move						(double dt);										//Move particles
    void MoveDEM    				(double dt);										//Move DEM particles

    void Solve					(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);		///< The solving function

    void CellInitiate		();															//Find the size of the domain as a cube, make cells and HOCs
    void ListGenerate		();															//Generate linked-list
    void CellReset			();															//Reset HOCs and particles' LL to initial value of -1

    void WriteXDMF			(char const * FileKey);					//Save a XDMF file for the visualization


    void InFlowBCLeave	();
    void InFlowBCFresh	();
    void WholeVelocity	();

    void Kernel_Set				(Kernels_Type const & KT);
    void Viscosity_Eq_Set		(Viscosity_Eq_Type const & VQ);
    void Gradient_Approach_Set	(Gradient_Type const & GT);

    // Data
    Array <Particle*>		    Particles; 	    ///< Array of particles
    Array <DEM::Particle*>	    DEMParticles; 	///< Array of DEM particles
    Array <DEM::CInteracton*>   DEMInteractons; ///< DEM particles interactions
    Array<Array <Vec3_t> >      FrictionMap;    ///< Map containing friction information between SPH and DEM particles
    double                      DEMstiff;  ///< Stiffness constant for DEM interaction
    double                      DEMdiss;   ///< Dissipation constant A for DEM interaction
    double                      DEMalpha;  ///< Dissipation constant B for DEM interaction
    double                      DEMDamp;   ///< Dissipation constant C=f(A,B) for DEM interaction
    double                      DEML;      ///< Characteristic length for DEM interaction

    double					    R;			///< Particle Radius in addrandombox

    double					    sqrt_h_a;	//Coefficient for determining Time Step based on acceleration (can be defined by user)

    int 					    Dimension;  ///< Dimension of the problem

    double					    MuMax;		///< Max Dynamic viscosity for calculating the timestep
    double					    CsMax;		///< Max speed of sound for calculating the timestep

    Vec3_t					    Gravity;	///< Gravity acceleration

    bool                        KGcorrection;

    Vec3_t                 	    TRPR;		///< Top right-hand point at rear of the domain as a cube
    Vec3_t                 	    BLPF;       ///< Bottom left-hand point at front of the domain as a cube
    Vec3_t                 	    CellSize;   ///< Calculated cell size according to (cell size >= 2h)
    int		               	    CellNo[3];  ///< No. of cells for linked list
    double 					    hmax;		///< Max of h for the cell size  determination
    Vec3_t                 	    DomSize;	///< Each component of the vector is the domain size in that direction if periodic boundary condition is defined in that direction as well
    double					    rhomax;

    int						    *** HOC;	///< Array of "Head of Chain" for each cell

    size_t					    SWIType;	///< Selecting variable to choose Soil-Water Interaction type
    bool					    FSI;		///< Selecting variable to choose Fluid-Structure Interaction

    double 					    XSPH;		///< Velocity correction factor

    double                      DeltaSPH;	     ///< Factor of Delta SPH that is between 0 and 1 depending on the DSPH term
    double  					DeltaStress;   ///< Factor of Delta Stresses that is between 0 and 10
    bool 						Shifting;      

    double 					    InitialDist;	///< Initial distance of particles for Inflow BC

    double					    AvgVelocity;	///< Average velocity of the last two column for x periodic constant velocity

    size_t					    Nproc;		///< No of threads which are going to use in parallel calculation
    omp_lock_t 				    dom_lock;	///< Open MP lock to lock Interactions array
    Boundary				    BC;
    PtOut					    UserOutput;
    PtVel 					    InCon;
    PtVel 					    OutCon;
    PtVel 					    AllCon;
    Vec3_t					    DomMax;
    Vec3_t					    DomMin;
    PtDom					    GeneralBefore;	///< Pointer to a function: to modify particles properties before CalcForce function
    PtDom					    GeneralAfter;	///< Pointer to a function: to modify particles properties after CalcForce function
    size_t					    Scheme;		    ///< Integration scheme: 0 = Modified Verlet, 1 = Leapfrog

    Array<Array<std::pair<size_t,size_t> > >	SMPairs;            ///< Same material pairs
    Array<Array<std::pair<size_t,size_t> > >	NSMPairs;           ///< Not same material pairs
    Array<Array<std::pair<size_t,size_t> > >	FSMPairs;           ///< Fixed same materials pairs
    Array< size_t > 				            FixedParticles;     ///< Array with the indexes of Fixed particles
    Array< size_t >				                FreeFSIParticles;
    Array< Array<size_t> >                      HomoGroups;         ///< Group lists for particle acceleration homogenization

    Array<std::pair<size_t,size_t> >		Initial;
    Mat3_t I;
    String					OutputName[3];


    //private:
    void Periodic_X_Correction	(Vec3_t & x, double const & h, Particle * P1, Particle * P2);		//Corrects xij for the periodic boundary condition
    void AdaptiveTimeStep				();		                                                    //Uses the minimum time step to smoothly vary the time step

    void PrintInput			(char const * FileKey);		//Print out some initial parameters as a file
    void InitialChecks	();		//Checks some parameter before proceeding to the solution
    void TimestepCheck	();		//Checks the user time step with CFL approach

    size_t					VisEq;					//Choose viscosity Eq based on different SPH discretisation
    size_t					KernelType;				//Choose a kernel
    size_t					GradientType;			//Choose a Gradient approach 1/Rho i^2 + 1/Rho j^2 or 1/(Rho i * Rho j)
    size_t					SoilCMType;				//Choose the constitutive model type for soil, elastoplastic or hypoplastic
    double 					Cellfac;				//Define the compact support of a kernel

    double					Time;    				//Current time of simulation at each solving step
    double					deltat;					//Time Step
    double					deltatmin;			//Minimum Time Step
    double					deltatint;			//Initial Time Step

};

void General(Domain & dom)
{
}

void OutPut(Particle * Particles, double & Prop1, double & Prop2,  double & Prop3)
{
	Prop1 = 0.0;
	Prop2 = 0.0;
	Prop3 = 0.0;
}

void InFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.inv;
	Den = bdry.inDensity;
}

void OutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.outv;
	Den = bdry.outDensity;
}

void AllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.allv;
	Den = bdry.allDensity;
}

// Constructor
inline Domain::Domain ()
{
    std::cout <<           "  " << std::endl;
    std::cout << BOLD(FBLU("|--------------------------------------------------------------------|")) << std::endl;
    std::cout << BOLD(FBLU("| ****************************************************************** |")) << std::endl;
    std::cout << BOLD(FBLU("| * -------------------------------------------------------------- * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |                                                            | * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |             SMOOTH PARTICLE HYDRODYNAMICS (SPH)            | * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |  A library from Mechsys: Multi-Physics Simulation Library  | * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |                                                           .| * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |                           .uu****.               .uu****u.*| * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |                          *  '**uuuu*         .uuuuuuuuuuuu*| * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |                      O         :uuuu*     .*uuuu.*.uuuuuuuu| * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |                               .*uuuuu*..*uuuuuu*  *uuuuuuuu| * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |***..                       ..**uuuuuuu**uuuuuu*    .uuuuuuu| * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu*..u*uuuuuuu| * |")) << std::endl;
    std::cout << BOLD(FBLU("| * |uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu| * |")) << std::endl;
    std::cout << BOLD(FBLU("| *  ------------------------------------------------------------  * |")) << std::endl;
    std::cout << BOLD(FBLU("| ****************************************************************** |")) << std::endl;
    std::cout << BOLD(FBLU("|--------------------------------------------------------------------|")) << std::endl;
    std::cout <<           "  " << std::endl;
    const time_t loctime = time(0);  // Get the local time 
    std::cout << "++++++++++++++ Started at the local time:  " << asctime(localtime(&loctime)) << std::endl; // Print the local time


    OutputName[0] = "Property1";
    OutputName[1] = "Property2";
    OutputName[2] = "Property3";
    Time    = 0.0;

    Dimension = 2;
    DomSize	= 0.0,0.0,0.0;

    Gravity	= 0.0,0.0,0.0;

    Cellfac = 2.0;

    KernelType	= 0;
    SWIType	= 0;
    FSI		= false;
    VisEq	= 0;
    Scheme	= 0;
    GradientType = 0;
    KGcorrection = false;

    DEMstiff = 0.1;
    DEMdiss  = 0.65;
    DEMalpha = 50.0;
    DEMDamp  = -(log(DEMdiss))/(DEMalpha*sqrt( log(DEMdiss)*log(DEMdiss) + M_PI*M_PI ) );
    DEML     = 0.0;

    XSPH	= 0.0;
    DeltaSPH = 0.0;
    DeltaStress = 0.0;
    Shifting = false;
    InitialDist = 0.0;

    AvgVelocity = 0.0;
    hmax	= 0.0;

    omp_init_lock (&dom_lock);
    Nproc	= 1;

    deltat	= 0.0;
    deltatint	= 0.0;
    deltatmin	= 0.0;
    sqrt_h_a = 0.0025;

    TRPR = 0.0;
    BLPF = 0.0;

    InCon = & InFlowCon;
    OutCon = & OutFlowCon;
    AllCon = & AllFlowCon;
    GeneralBefore = & General;
    GeneralAfter = & General;
    UserOutput = & OutPut;

    DomMax = -100000000000.0;
    DomMin = 100000000000.0;
    I = OrthoSys::I;
}

inline Domain::~Domain ()
{
	size_t Max = Particles.Size();
	for (size_t i=1; i<=Max; i++)  Particles.DelItem(Max-i);
}

//All the methods for DEM particle generation

#include<mechsys/sph/Dompargen.h>

inline void Domain::Periodic_X_Correction(Vec3_t & x, double const & h, Particle * P1, Particle * P2)
{
	if (DomSize(0)>0.0) {if (x(0)>2*Cellfac*h || x(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? x(0) -= DomSize(0) : x(0) += DomSize(0);}}
	if (DomSize(1)>0.0) {if (x(1)>2*Cellfac*h || x(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? x(1) -= DomSize(1) : x(1) += DomSize(1);}}
	if (DomSize(2)>0.0) {if (x(2)>2*Cellfac*h || x(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? x(2) -= DomSize(2) : x(2) += DomSize(2);}}
}

inline void Domain::Kernel_Set(Kernels_Type const & KT)
{
	KernelType = KT;
	if (KernelType==2) Cellfac = 3.0; else Cellfac = 2.0;
}

inline void Domain::Viscosity_Eq_Set(Viscosity_Eq_Type const & VQ)
{
	VisEq = VQ;
}

inline void Domain::Gradient_Approach_Set(Gradient_Type const & GT)
{
	GradientType = GT;
}

inline void Domain::AdaptiveTimeStep()
{
	if (deltatint>deltatmin)
	{
		if (deltat<deltatmin)
			deltat		= 2.0*deltat*deltatmin/(deltat+deltatmin);
		else
			deltat		= deltatmin;
	}
	else
	{
		if (deltatint!=deltat)
			deltat		= 2.0*deltat*deltatint/(deltat+deltatint);
		else
			deltat		= deltatint;
	}

	if (deltat<(deltatint/1.0e5))
		throw new Fatal("Too large time step, please choose a smaller time step initially to make the simulation more stable");
}

inline void Domain::AddSingleParticle(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed)
{
    if (DEML==0.0 || 2.0*h/Cellfac < DEML) DEML = 2.0*h/Cellfac;

   	Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),Mass,Density,h,Fixed));
}

inline void Domain::AddBoxNo(int tag, Vec3_t const & V, size_t nx, size_t ny, size_t nz, double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
{
    if (DEML==0.0 || 2.0*r < DEML) DEML = 2.0*r;

    if (!(type==0 || type==1))
    {
	   	std::cout << "Packing Type is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Hexagonal Close Packing" << std::endl;
		std::cout << "1 => Cubic Packing" << std::endl;
	    abort();
    }

    if (!(rotation==0 || rotation==90))
    {
	   	std::cout << "Packing Rotation Angle is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => " << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << std::endl;
		std::cout << "90 =>" << std::endl;
		std::cout << "  0   0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0   0  " << std::endl;
		abort();
    }

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by AddBoxNo with defined numbers of particles--------------" << std::endl;

    size_t PrePS = Particles.Size();

    double x,y;

    double qin = 0.03;
    srand(100);

    if (Dimension==3)
    {
   		double z;

    	if (type==0)
    	{
    		//Hexagonal close packing
 		    for (size_t k=0; k<nz; k++)
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				if ((k%2!=0) && (j%2!=0)) x = V(0) + (2*i+(j%2)+(k%2)-1)*r; else x = V(0) + (2*i+(j%2)+(k%2)+1)*r;
				y = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
				z = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
			}
    	}
    	else
    	{
    		//Cubic packing
		    for (size_t k=0; k<nz; k++)
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				x = V(0) + (2.0*i+1)*r;
				y = V(1) + (2.0*j+1)*r;
				z = V(2) + (2.0*k+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
			}
    	}

        //Calculate particles' mass in 3D
        Vec3_t temp, Max=V;
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS);

		#pragma omp parallel for num_threads(Nproc)
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			Particles[i]->Mass = Mass;
		}
    }

    if (Dimension==2)
    {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		if (rotation==0)
    		{
    		    for (size_t j=0; j<ny; j++)
    		    for (size_t i=0; i<nx; i++)
    		    {
					x = V(0) + (2*i+(j%2)+1)*r;
					y = V(1) + (sqrt(3.0)*j+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
				}
			}
    		else
    		{
    		    for (size_t i=0; i<nx; i++)
    		    for (size_t j=0; j<ny; j++)
    		    {
					x = V(0) + (sqrt(3.0)*i+1)*r;
					y = V(1) + (2*j+(i%2)+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
				}
    		}
    	}
    	else
    	{
    		//Cubic packing
		    for (size_t j=0; j<ny; j++)
		    for (size_t i=0; i<nx; i++)
		    {
				x = V(0) + (2*i+1)*r;
				y = V(1) + (2*j+1)*r;
				if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
					else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed));
			}

    	}
    }

	R = r;
}

inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
{
    if (DEML==0.0 || 2.0*r < DEML) DEML = 2.0*r;

    if (!(type==0 || type==1))
    {
	   	std::cout << "Packing Type is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Hexagonal Close Packing" << std::endl;
		std::cout << "1 => Cubic Packing" << std::endl;
	    abort();
    }

    if (!(rotation==0 || rotation==90))
    {
	   	std::cout << "Packing Rotation Angle is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => " << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << std::endl;
		std::cout << "90 => Cubic Close Packing" << std::endl;
		std::cout << "  0   0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0   0  " << std::endl;
		abort();
    }

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by AddBoxLength with defined length of particles-----------" << std::endl;

    size_t PrePS = Particles.Size();

    double x,y,xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);

    if (Dimension==3)
    {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		double z,zp;
			size_t k=0;
			zp = V(2);

			while (zp <= (V(2)+Lz-r))
			{
				j = 0;
				yp = V(1);
				while (yp <= (V(1)+Ly-r))
				{
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						if ((k%2!=0) && (j%2!=0)) x = V(0) + (2*i+(j%2)+(k%2)-1)*r; else x = V(0) + (2*i+(j%2)+(k%2)+1)*r;
						y = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
						z = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
						i++;
						if ((k%2!=0) && (j%2!=0)) xp = V(0) + (2*i+(j%2)+(k%2)-1)*r; else xp = V(0) + (2*i+(j%2)+(k%2)+1)*r;
					}
					j++;
					yp = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
				}
				k++;
				zp = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
			}
    	}
    	else
    	{
    		//Cubic packing
    		double z,zp;
			size_t k=0;
			zp = V(2);

			while (zp <= (V(2)+Lz-r))
			{
				j = 0;
				yp = V(1);
				while (yp <= (V(1)+Ly-r))
				{
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						x = V(0) + (2.0*i+1)*r;
						y = V(1) + (2.0*j+1)*r;
						z = V(2) + (2.0*k+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
						i++;
						xp = V(0) + (2*i+1)*r;
					}
					j++;
					yp = V(1) + (2.0*j+1)*r;
				}
				k++;
				zp = V(2) + (2.0*k+1)*r;
			}
    	}

        //Calculate particles' mass in 3D
        Vec3_t temp, Max=V;
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS);

		#pragma omp parallel for num_threads(Nproc)
		for (size_t i=PrePS; i<Particles.Size(); i++)
		{
			Particles[i]->Mass = Mass;
		}
    }

    if (Dimension==2)
    {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		if (rotation==0)
    		{
				j = 0;
				yp = V(1);

				while (yp <= (V(1)+Ly-r))
				{
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						x = V(0) + (2*i+(j%2)+1)*r;
						y = V(1) + (sqrt(3.0)*j+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						i++;
						xp = V(0) + (2*i+(j%2)+1)*r;
					}
					j++;
					yp = V(1) + (sqrt(3.0)*j+1)*r;
				}
			}
    		else
    		{
				i = 0;
				xp = V(0);

				while (xp <= (V(0)+Lx-r))
				{
					j = 0;
					yp = V(1);
					while (yp <= (V(1)+Ly-r))
					{
						x = V(0) + (sqrt(3.0)*i+1)*r;
						y = V(1) + (2*j+(i%2)+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						j++;
						yp = V(1) + (2*j+(i%2)+1)*r;
					}
					i++;
					xp = V(0) + (sqrt(3.0)*i+1)*r;
				}
    		}
    	}
    	else
    	{
    		//Cubic packing
    		j = 0;
			yp = V(1);

			while (yp <= (V(1)+Ly-r))
			{
				i = 0;
				xp = V(0);
				while (xp <= (V(0)+Lx-r))
				{
					x = V(0) + (2*i+1)*r;
					y = V(1) + (2*j+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed));
					i++;
					xp = V(0) + (2*i+1)*r;
				}
				j++;
				yp = V(1) + (2*j+1)*r;
			}

    	}
    }

	R = r;
}

inline void Domain::DelParticles (int const & Tags)
{
    Array<int> idxs; // indices to be deleted

	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->ID==Tags)
		{
			omp_set_lock(&dom_lock);
        	idxs.Push(i);
			omp_unset_lock(&dom_lock);
		}
    }
    if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particles to delete");
    Particles.DelItems (idxs);

    std::cout << "\n" << "Particle(s) with Tag No. " << Tags << " has been deleted" << std::endl;
}

inline void Domain::CheckParticleLeave ()
{
	Array <int> DelParticles;

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
    {
		if ((Particles[i]->x(0) > TRPR(0)) || (Particles[i]->x(1) > TRPR(1)) || (Particles[i]->x(2) > TRPR(2)) ||
			(Particles[i]->x(0) < BLPF(0)) || (Particles[i]->x(1) < BLPF(1)) || (Particles[i]->x(2) < BLPF(2)))
		{
			omp_set_lock(&dom_lock);
			DelParticles.Push(i);
			omp_unset_lock(&dom_lock);
		}
    }

	if (DelParticles.Size()>0)
	{
		std::cout<< DelParticles.Size()<< " particle(s) left the Domain"<<std::endl;
		for (size_t i=0; i<DelParticles.Size(); i++)
		{
			std::cout<<""<<std::endl;
			std::cout<<"Particle Number   = "<<DelParticles[i]<<std::endl;
			std::cout<<"Particle Material = "<<Particles[DelParticles[i]]->Material<<std::endl;
			std::cout<<"x = "<<Particles[DelParticles[i]]->x<<std::endl;
			std::cout<<"v = "<<Particles[DelParticles[i]]->v<<std::endl;
			std::cout<<"a = "<<Particles[DelParticles[i]]->a<<std::endl;
		}
		Particles.DelItems(DelParticles);
	}
}

inline void Domain::CellInitiate ()
{
	if (!(norm(TRPR)>0.0) && !(norm(BLPF)>0.0))
	{
		// Calculate Domain Size
		BLPF = Particles[0]->x;
		TRPR = Particles[0]->x;
		hmax = Particles[0]->h;
		rhomax = Particles[0]->Density;

		for (size_t i=0; i<Particles.Size(); i++)
		{
			if (Particles[i]->x(0) > TRPR(0)) TRPR(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > TRPR(1)) TRPR(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > TRPR(2)) TRPR(2) = Particles[i]->x(2);

			if (Particles[i]->x(0) < BLPF(0)) BLPF(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) < BLPF(1)) BLPF(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) < BLPF(2)) BLPF(2) = Particles[i]->x(2);

			if (Particles[i]->h > hmax) hmax=Particles[i]->h;
			if (Particles[i]->Density > rhomax) rhomax=Particles[i]->Density;
			if (Particles[i]->Mu > MuMax) MuMax=Particles[i]->Mu;
			if (Particles[i]->Cs > CsMax) CsMax=Particles[i]->Cs;
		}
	}

	// Override the calculated domain size
	if (DomMax(0)>TRPR(0)) TRPR(0) = DomMax(0);
	if (DomMax(1)>TRPR(1)) TRPR(1) = DomMax(1);
	if (DomMax(2)>TRPR(2)) TRPR(2) = DomMax(2);
	if (DomMin(0)<BLPF(0)) BLPF(0) = DomMin(0);
	if (DomMin(1)<BLPF(1)) BLPF(1) = DomMin(1);
	if (DomMin(2)<BLPF(2)) BLPF(2) = DomMin(2);


	//Because of Hexagonal close packing in x direction domain is modified
	if (!BC.Periodic[0]) {TRPR(0) += hmax/2;	BLPF(0) -= hmax/2;}else{TRPR(0) += R; BLPF(0) -= R;}
	if (!BC.Periodic[1]) {TRPR(1) += hmax/2;	BLPF(1) -= hmax/2;}else{TRPR(1) += R; BLPF(1) -= R;}
	if (!BC.Periodic[2]) {TRPR(2) += hmax/2;	BLPF(2) -= hmax/2;}else{TRPR(2) += R; BLPF(2) -= R;}

    // Calculate Cells Properties
	switch (Dimension)
	{case 2:
		if (double (ceil(((TRPR(0)-BLPF(0))/(Cellfac*hmax)))-((TRPR(0)-BLPF(0))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[0] = int(ceil((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		else
			CellNo[0] = int(floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(1)-BLPF(1))/(Cellfac*hmax)))-((TRPR(1)-BLPF(1))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[1] = int(ceil((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
		else
			CellNo[1] = int(floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));

		CellNo[2] = 1;

		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],0.0);
		break;

	case 3:
		if (double (ceil(((TRPR(0)-BLPF(0))/(Cellfac*hmax)))-((TRPR(0)-BLPF(0))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[0] = int(ceil((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		else
			CellNo[0] = int(floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(1)-BLPF(1))/(Cellfac*hmax)))-((TRPR(1)-BLPF(1))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[1] = int(ceil((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
		else
			CellNo[1] = int(floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(2)-BLPF(2))/(Cellfac*hmax)))-((TRPR(2)-BLPF(2))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[2] = int(ceil((TRPR(2)-BLPF(2))/(Cellfac*hmax)));
		else
			CellNo[2] = int(floor((TRPR(2)-BLPF(2))/(Cellfac*hmax)));

		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],(TRPR(2)-BLPF(2))/CellNo[2]);
		break;

	default:
    	std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
		break;
	}

	// Periodic BC modifications
	if (BC.Periodic[0]) CellNo[0] += 2;
    if (BC.Periodic[1]) CellNo[1] += 2;
    if (BC.Periodic[2]) CellNo[2] += 2;

    if (BC.Periodic[0]) DomSize[0] = (TRPR(0)-BLPF(0));
    if (BC.Periodic[1]) DomSize[1] = (TRPR(1)-BLPF(1));
    if (BC.Periodic[2]) DomSize[2] = (TRPR(2)-BLPF(2));

    // Initiate Head of Chain array for Linked-List
    HOC = new int**[(int) CellNo[0]];
    for(int i =0; i<CellNo[0]; i++){
       HOC[i] = new int*[CellNo[1]];
       for(int j =0; j<CellNo[1]; j++){
           HOC[i][j] = new int[CellNo[2]];
           for(int k = 0; k<CellNo[2];k++){
              HOC[i][j][k] = -1;
           }
       }
    }
    // Initiate Pairs array for neibour searching
    for(size_t i=0 ; i<Nproc ; i++)
    {
	    SMPairs.Push(Initial);
	    NSMPairs.Push(Initial);
	    FSMPairs.Push(Initial);
    }
}

inline void Domain::ListGenerate ()
{
	int i, j, k, temp=0;
	switch (Dimension)
	{case 2:
		for (size_t a=0; a<Particles.Size(); a++)
		{
			i= (int) (floor((Particles[a]->x(0) - BLPF(0)) / CellSize(0)));
			j= (int) (floor((Particles[a]->x(1) - BLPF(1)) / CellSize(1)));

			if (i<0)
            {
                    if ((BLPF(0) - Particles[a]->x(0)) <= hmax) i=0;
                            else std::cout<<"Leaving i<0"<<std::endl;
            }
            if (j<0)
            {
                    if ((BLPF(1) - Particles[a]->x(1)) <= hmax) j=0;
                            else std::cout<<"Leaving j<0"<<std::endl;
            }
			if (i>=CellNo[0])
			{
					if ((Particles[a]->x(0) - TRPR(0)) <= hmax) i=CellNo[0]-1;
							else std::cout<<"Leaving i>=CellNo"<<std::endl;
			}
            if (j>=CellNo[1])
            {
                    if ((Particles[a]->x(1) - TRPR(1)) <= hmax) j=CellNo[1]-1;
                            else std::cout<<"Leaving j>=CellNo"<<std::endl;
            }

			temp = HOC[i][j][0];
			HOC[i][j][0] = a;
			Particles[a]->LL = temp;
			Particles[a]->CC[0] = i;
			Particles[a]->CC[1] = j;
			Particles[a]->CC[2] = 0;
			if (!Particles[a]->IsFree) FixedParticles.Push(a);
		}
		break;

	case 3:
		for (size_t a=0; a<Particles.Size(); a++)
		{
			i= (int) (floor((Particles[a]->x(0) - BLPF(0)) / CellSize(0)));
			j= (int) (floor((Particles[a]->x(1) - BLPF(1)) / CellSize(1)));
			k= (int) (floor((Particles[a]->x(2) - BLPF(2)) / CellSize(2)));

            if (i<0)
            {
                    if ((BLPF(0) - Particles[a]->x(0))<=hmax) i=0;
                            else std::cout<<"Leaving"<<std::endl;
            }
            if (j<0)
            {
                    if ((BLPF(1) - Particles[a]->x(1))<=hmax) j=0;
                            else std::cout<<"Leaving"<<std::endl;
            }
            if (k<0)
            {
                    if ((BLPF(2) - Particles[a]->x(2))<=hmax) k=0;
                            else std::cout<<"Leaving"<<std::endl;
            }
			if (i>=CellNo[0])
			{
					if ((Particles[a]->x(0) - TRPR(0))<=hmax) i=CellNo[0]-1;
							else std::cout<<"Leaving"<<std::endl;
			}
            if (j>=CellNo[1])
            {
                    if ((Particles[a]->x(1) - TRPR(1))<=hmax) j=CellNo[1]-1;
                            else std::cout<<"Leaving"<<std::endl;
            }
            if (k>=CellNo[2])
            {
                    if ((Particles[a]->x(2) - TRPR(2))<=hmax) k=CellNo[2]-1;
                            else std::cout<<"Leaving"<<std::endl;
            }

            temp = HOC[i][j][k];
			HOC[i][j][k] = a;
			Particles[a]->LL = temp;
			Particles[a]->CC[0] = i;
			Particles[a]->CC[1] = j;
			Particles[a]->CC[2] = k;
			if (!Particles[a]->IsFree) FixedParticles.Push(a);
		}
		break;

	default:
    	std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
		break;
	}

	if (BC.Periodic[0])
	{
	   for(int j =0; j<CellNo[1]; j++)
	   for(int k =0; k<CellNo[2]; k++)
	   {
		  HOC[CellNo[0]-1][j][k] =  HOC[1][j][k];
		  HOC[CellNo[0]-2][j][k] =  HOC[0][j][k];
	   }
	}
	if (BC.Periodic[1])
	{
	   for(int i =0; i<CellNo[0]; i++)
	   for(int k =0; k<CellNo[2]; k++)
	   {
		  HOC[i][CellNo[1]-1][k] =  HOC[i][1][k];
		  HOC[i][CellNo[1]-2][k] =  HOC[i][0][k];
	   }
	}
	if (BC.Periodic[2])
	{
	   for(int i =0; i<CellNo[0]; i++)
	   for(int j =0; j<CellNo[1]; j++)
	   {
		  HOC[i][j][CellNo[2]-1] =  HOC[i][j][1];
		  HOC[i][j][CellNo[2]-2] =  HOC[i][j][0];
	   }
	}
}

inline void Domain::CellReset ()
{

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for(int i =0; i<CellNo[0]; i++)
    {
		for(int j =0; j<CellNo[1]; j++)
		for(int k =0; k<CellNo[2];k++)
		{
			HOC[i][j][k] = -1;
		}
    }

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t a=0; a<Particles.Size(); a++)
    {
    	Particles[a]->LL = -1;
    }

    FixedParticles.Clear();
}

inline void Domain::MainNeighbourSearch()
{
    int q1;

    if (BC.Periodic[0])
    {
	#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
	for (q1=1;q1<(CellNo[0]-1); q1++)	YZPlaneCellsNeighbourSearch(q1);
    }
    else
    {
	#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
    	for (q1=0;q1<CellNo[0]; q1++)	YZPlaneCellsNeighbourSearch(q1);
    }
}

inline void Domain::YZPlaneCellsNeighbourSearch(int q1)
{
	int q3,q2;
	size_t T = omp_get_thread_num();

	for (BC.Periodic[2] ? q3=1 : q3=0;BC.Periodic[2] ? (q3<(CellNo[2]-1)) : (q3<CellNo[2]); q3++)
	for (BC.Periodic[1] ? q2=1 : q2=0;BC.Periodic[1] ? (q2<(CellNo[1]-1)) : (q2<CellNo[1]); q2++)
	{
		if (HOC[q1][q2][q3]==-1) continue;
		else
		{
			int temp1, temp2;
			temp1 = HOC[q1][q2][q3];

			while (temp1 != -1)
			{
				// The current cell  => self cell interactions
				temp2 = Particles[temp1]->LL;
				while (temp2 != -1)
				{
					if (Particles[temp1]->IsFree || Particles[temp2]->IsFree)
					{
						if (Particles[temp1]->Material == Particles[temp2]->Material)
						{
							if (Particles[temp1]->IsFree&&Particles[temp2]->IsFree)
								SMPairs[T].Push(std::make_pair(temp1, temp2));
							else
								FSMPairs[T].Push(std::make_pair(temp1, temp2));

						}
						else if (Particles[temp1]->IsFree && Particles[temp2]->IsFree)
                        {
							NSMPairs[T].Push(std::make_pair(temp1, temp2));
                        }
					}
					temp2 = Particles[temp2]->LL;
				}

				// (q1 + 1, q2 , q3)
				if (q1+1< CellNo[0])
				{
					temp2 = HOC[q1+1][q2][q3];
					while (temp2 != -1)
					{
						if (Particles[temp1]->IsFree || Particles[temp2]->IsFree)
						{
							if (Particles[temp1]->Material == Particles[temp2]->Material)
							{
								if (Particles[temp1]->IsFree&&Particles[temp2]->IsFree)
									SMPairs[T].Push(std::make_pair(temp1, temp2));
								else
									FSMPairs[T].Push(std::make_pair(temp1, temp2));

							}
						    else if (Particles[temp1]->IsFree && Particles[temp2]->IsFree)
                            {
								NSMPairs[T].Push(std::make_pair(temp1, temp2));
                            }
						}
						temp2 = Particles[temp2]->LL;
					}
				}

				// (q1 + a, q2 + 1, q3) & a[-1,1]
				if (q2+1< CellNo[1])
				{
					for (int i = q1-1; i <= q1+1; i++)
					{
						if (i<CellNo[0] && i>=0)
						{
							temp2 = HOC[i][q2+1][q3];
							while (temp2 != -1)
							{
								if (Particles[temp1]->IsFree || Particles[temp2]->IsFree)
								{
									if (Particles[temp1]->Material == Particles[temp2]->Material)
									{
										if (Particles[temp1]->IsFree&&Particles[temp2]->IsFree)
											SMPairs[T].Push(std::make_pair(temp1, temp2));
										else
											FSMPairs[T].Push(std::make_pair(temp1, temp2));

									}
						            else if (Particles[temp1]->IsFree && Particles[temp2]->IsFree)
                                    {
							        	NSMPairs[T].Push(std::make_pair(temp1, temp2));
                                    }
								}
								temp2 = Particles[temp2]->LL;
							}
						}
					}
				}

				// (q1 + a, q2 + b, q3 + 1) & a,b[-1,1] => all 9 cells above the current cell
				if (q3+1< CellNo[2])
				{
					for (int j=q2-1; j<=q2+1; j++)
					for (int i=q1-1; i<=q1+1; i++)
					{
						if (i<CellNo[0] && i>=0 && j<CellNo[1] && j>=0)
						{
							temp2 = HOC[i][j][q3+1];
							while (temp2 != -1)
							{
								if (Particles[temp1]->IsFree || Particles[temp2]->IsFree)
								{
									if (Particles[temp1]->Material == Particles[temp2]->Material)
									{
										if (Particles[temp1]->IsFree&&Particles[temp2]->IsFree)
											SMPairs[T].Push(std::make_pair(temp1, temp2));
										else
											FSMPairs[T].Push(std::make_pair(temp1, temp2));

									}
						            else if (Particles[temp1]->IsFree && Particles[temp2]->IsFree)
                                    {
							        	NSMPairs[T].Push(std::make_pair(temp1, temp2));
                                    }
								}
								temp2 = Particles[temp2]->LL;
							}
						}
					}
				}
				temp1 = Particles[temp1]->LL;
			}
		}
	}
}

inline void Domain::StartAcceleration (Vec3_t const & a)
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
	{
	    	if (Particles[i]->IsFree)
    		{
			// Tensile Instability for all soil and solid particles
			if (Particles[i]->Material > 1 && Particles[i]->TI > 0.0)
        		{
				// XY plane must be used, It is very slow in 3D
				if (Dimension == 2)
				{
					double teta, Sigmaxx, Sigmayy, C, S;

					if ((Particles[i]->Sigma(0,0)-Particles[i]->Sigma(1,1))!=0.0)
						teta = 0.5*atan(2.0*Particles[i]->Sigma(0,1)/(Particles[i]->Sigma(0,0)-Particles[i]->Sigma(1,1)));
					else
						teta = M_PI/4.0;

					C = cos(teta);
					S = sin(teta);
					Sigmaxx = C*C*Particles[i]->Sigma(0,0) + 2.0*C*S*Particles[i]->Sigma(0,1) + S*S*Particles[i]->Sigma(1,1);
					Sigmayy = S*S*Particles[i]->Sigma(0,0) - 2.0*C*S*Particles[i]->Sigma(0,1) + C*C*Particles[i]->Sigma(1,1);
					if (Sigmaxx>0) Sigmaxx = -Particles[i]->TI * Sigmaxx/(Particles[i]->Density*Particles[i]->Density); else Sigmaxx = 0.0;
					if (Sigmayy>0) Sigmayy = -Particles[i]->TI * Sigmayy/(Particles[i]->Density*Particles[i]->Density); else Sigmayy = 0.0;
					Particles[i]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					Particles[i]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					Particles[i]->TIR(0,1) = Particles[i]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				}
				else
				{
					Mat3_t Vec,Val,VecT,temp;

					MatRotation(Particles[i]->Sigma,Vec,VecT,Val);
					if (Val(0,0)>0) Val(0,0) = -Particles[i]->TI * Val(0,0)/(Particles[i]->Density*Particles[i]->Density); else Val(0,0) = 0.0;
					if (Val(1,1)>0) Val(1,1) = -Particles[i]->TI * Val(1,1)/(Particles[i]->Density*Particles[i]->Density); else Val(1,1) = 0.0;
					if (Val(2,2)>0) Val(2,2) = -Particles[i]->TI * Val(2,2)/(Particles[i]->Density*Particles[i]->Density); else Val(2,2) = 0.0;
					Mult(Vec,Val,temp);
					Mult(temp,VecT,Particles[i]->TIR);
				}
			}
	    	}
	    	else
	    	{
	       		// Reset the pressure and the induced velocity for solid boundaries
	    		Particles[i]->NSv = 0.0;
	    		Particles[i]->Pressure = 0.0;
	        	set_to_zero(Particles[i]->Sigma);
	        	set_to_zero(Particles[i]->ShearStress);
	    	}

        Particles[i]->Lcorr.Resize(Dimension,Dimension);
        Particles[i]->Lcorr.SetValues(0.0);

		//Reset to zero for all particles
		Particles[i]->a		= a;
		Particles[i]->SatCheck	= false;
		Particles[i]->dDensity	= 0.0;
		Particles[i]->VXSPH	= 0.0;
		Particles[i]->ZWab	= 0.0;
		Particles[i]->SumDen	= 0.0;
		Particles[i]->SumKernel	= 0.0;
		Particles[i]->FSISumKernel	= 0.0;
		Particles[i]->FSINSv		= 0.0;
		Particles[i]->FSIPressure	= 0.0;
		Particles[i]->qmin		= 1.0e16;
//	        set_to_zero(Particles[i]->FSISigma);
		if (Dimension == 2) Particles[i]->v(2) = 0.0;
		set_to_zero(Particles[i]->StrainRate);
		set_to_zero(Particles[i]->RotationRate);
        set_to_zero(Particles[i]->FilStress);
        set_to_zero(Particles[i]->Shift);
		Particles[i]->S		= 0.0;
	}

	#pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t j=0;j<DEMParticles.Size();j++)
    {
        DEMParticles[j]->F = DEMParticles[j]->Ff;
        DEMParticles[j]->T = DEMParticles[j]->Tf;
    }
}

inline void Domain::PrimaryComputeAcceleration ()
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t k=0; k<Nproc;k++)
	{
		size_t P1,P2;
		Vec3_t xij;
		double h,K;
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		for (size_t a=0; a<FSMPairs[k].Size();a++)
		{
			P1	= FSMPairs[k][a].first;
			P2	= FSMPairs[k][a].second;
			xij	= Particles[P1]->x-Particles[P2]->x;
			h	= (Particles[P1]->h+Particles[P2]->h)/2.0;

			Periodic_X_Correction(xij, h, Particles[P1], Particles[P2]);


			K	= Kernel(Dimension, KernelType, norm(xij)/h, h);

			if (!Particles[P1]->IsFree)
			{
				omp_set_lock(&Particles[P1]->my_lock);
				Particles[P1]->SumKernel+= K;
				if (Particles[P1]->Material < 3)	Particles[P1]->Pressure	+= Particles[P2]->Pressure * K + dot(Gravity,xij)*Particles[P2]->Density*K;
				if (Particles[P1]->Material > 1)	Particles[P1]->Sigma 	 = Particles[P1]->Sigma + K * Particles[P2]->Sigma;
				if (Particles[P1]->NoSlip)		Particles[P1]->NSv 	+= Particles[P2]->v * K;
				omp_unset_lock(&Particles[P1]->my_lock);
			}
			else
			{
				omp_set_lock(&Particles[P2]->my_lock);
				Particles[P2]->SumKernel+= K;
				if (Particles[P2]->Material < 3)	Particles[P2]->Pressure	+= Particles[P1]->Pressure * K + dot(Gravity,xij)*Particles[P1]->Density*K;
				if (Particles[P2]->Material > 1)	Particles[P2]->Sigma	 = Particles[P2]->Sigma + K * Particles[P1]->Sigma;
				if (Particles[P2]->NoSlip)		Particles[P2]->NSv 	+= Particles[P1]->v * K;
				omp_unset_lock(&Particles[P2]->my_lock);
			}
		}
		if (SWIType < 2)
		{
			for (size_t a=0; a<NSMPairs[k].Size();a++)
			{
				P1 = NSMPairs[k][a].first;
				P2 = NSMPairs[k][a].second;
				if (Particles[P1]->Material == 3 && Particles[P1]->Material*Particles[P2]->Material == 3)
				{
					if (!Particles[P1]->SatCheck)
						if (Particles[P2]->CC[1] >= Particles[P1]->CC[1])
							if (Particles[P2]->x(1) >= Particles[P1]->x(1))
							{
								omp_set_lock(&Particles[P1]->my_lock);
									Particles[P1]->SatCheck = true;
								omp_unset_lock(&Particles[P1]->my_lock);
							}
				}
				if (Particles[P2]->Material == 3  && Particles[P1]->Material*Particles[P2]->Material == 3)
				{
					if (!Particles[P2]->SatCheck)
						if (Particles[P1]->CC[1] >= Particles[P2]->CC[1])
							if (Particles[P1]->x(1) >= Particles[P2]->x(1))
							{
								omp_set_lock(&Particles[P2]->my_lock);
									Particles[P2]->SatCheck = true;
								omp_unset_lock(&Particles[P2]->my_lock);
							}
				}
			}
		}
		if (FSI)
		{
			for (size_t a=0; a<NSMPairs[k].Size();a++)
			{
				P1 = NSMPairs[k][a].first;
				P2 = NSMPairs[k][a].second;
				if (Particles[P1]->Material*Particles[P2]->Material == 2 && (Particles[P1]->IsFree && Particles[P2]->IsFree))
				{

					xij	= Particles[P1]->x-Particles[P2]->x;
					h	= (Particles[P1]->h+Particles[P2]->h)/2.0;

					Periodic_X_Correction(xij, h, Particles[P1], Particles[P2]);

					K	= Kernel(Dimension, KernelType, norm(xij)/h, h);

					if (Particles[P1]->Material == 1)
					{
						omp_set_lock(&dom_lock);
							FreeFSIParticles.Push(P2);
						omp_unset_lock(&dom_lock);

						omp_set_lock(&Particles[P2]->my_lock);
							Particles[P2]->FSISumKernel	+= K;
							Particles[P2]->NSv 		+= Particles[P1]->v * K;
							Particles[P2]->FSIPressure	+= Particles[P1]->Pressure * K + dot(Gravity,xij)*Particles[P1]->Density*K;
						omp_unset_lock(&Particles[P2]->my_lock);

//						Particles[P1]->FSISumKernel	+= K;
//						Particles[P1]->NSv 		+= Particles[P2]->v * K;
//						Particles[P1]->FSIPressure	+= Particles[P2]->Pressure * K + dot(Gravity,xij)*Particles[P2]->Density*K;
//						Particles[P1]->FSISigma	 	 = Particles[P1]->FSISigma + K * Particles[P2]->Sigma;
					}
					else
					{
						omp_set_lock(&dom_lock);
							FreeFSIParticles.Push(P1);
						omp_unset_lock(&dom_lock);

//						Particles[P2]->FSISumKernel	+= K;
//						Particles[P2]->FSINSv 		+= Particles[P1]->v * K;
//						Particles[P2]->FSIPressure	+= Particles[P1]->Pressure * K + dot(Gravity,xij)*Particles[P1]->Density*K;
//						Particles[P2]->FSISigma	 	 = Particles[P2]->FSISigma + K * Particles[P1]->Sigma;

						omp_set_lock(&Particles[P1]->my_lock);
							Particles[P1]->FSISumKernel	+= K;
							Particles[P1]->FSINSv 		+= Particles[P2]->v * K;
							Particles[P1]->FSIPressure	+= Particles[P2]->Pressure * K + dot(Gravity,xij)*Particles[P2]->Density*K;
						omp_unset_lock(&Particles[P1]->my_lock);
					}

				}
			}
		}

	}

	if (FSI)
	{
		// Calculateing the finala value of the smoothed pressure, velocity and stress for fixed particles
		#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (size_t i=0; i<FreeFSIParticles.Size(); i++)
		{
			size_t a = FreeFSIParticles[i];
			Particles[a]->FSIPressure	= Particles[a]->FSIPressure/Particles[a]->FSISumKernel;
//			Particles[a]->FSISigma		= 1.0/Particles[a]->FSISumKernel*Particles[a]->FSISigma;
			Particles[a]->FSINSv		= Particles[a]->FSINSv/Particles[a]->FSISumKernel;

		}
		FreeFSIParticles.Clear();
	}


	// Calculateing the finala value of the smoothed pressure, velocity and stress for fixed particles
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<FixedParticles.Size(); i++)
		if (Particles[FixedParticles[i]]->SumKernel!= 0.0)
		{
			size_t a = FixedParticles[i];
			if (Particles[a]->Material < 3)	Particles[a]->Pressure	= Particles[a]->Pressure/Particles[a]->SumKernel;
			if (Particles[a]->Material > 1) Particles[a]->Sigma	= 1.0/Particles[a]->SumKernel*Particles[a]->Sigma;
			if (Particles[a]->NoSlip)	Particles[a]->NSv	= Particles[a]->NSv/Particles[a]->SumKernel;

			// Tensile Instability for fixed soil and solid particles
			if (Particles[a]->Material > 1 && Particles[a]->TI > 0.0)
			{
				// XY plane must be used, It is very slow in 3D
				if (Dimension == 2)
				{
					double teta, Sigmaxx, Sigmayy, C, S;
					if ((Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1))!=0.0)
						teta = 0.5*atan(2.0*Particles[a]->Sigma(0,1)/(Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1)));
					else
						teta = M_PI/4.0;

					C = cos(teta);
					S = sin(teta);
					Sigmaxx = C*C*Particles[a]->Sigma(0,0) + 2.0*C*S*Particles[a]->Sigma(0,1) + S*S*Particles[a]->Sigma(1,1);
					Sigmayy = S*S*Particles[a]->Sigma(0,0) - 2.0*C*S*Particles[a]->Sigma(0,1) + C*C*Particles[a]->Sigma(1,1);
					if (Sigmaxx>0) Sigmaxx = -Particles[a]->TI * Sigmaxx/(Particles[a]->Density*Particles[a]->Density); else Sigmaxx = 0.0;
					if (Sigmayy>0) Sigmayy = -Particles[a]->TI * Sigmayy/(Particles[a]->Density*Particles[a]->Density); else Sigmayy = 0.0;
					Particles[a]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					Particles[a]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					Particles[a]->TIR(0,1) = Particles[a]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				}
				else
				{
					Mat3_t Vec,Val,VecT,temp;
					MatRotation(Particles[a]->Sigma,Vec,VecT,Val);
					if (Val(0,0)>0) Val(0,0) = -Particles[a]->TI * Val(0,0)/(Particles[a]->Density*Particles[a]->Density); else Val(0,0) = 0.0;
					if (Val(1,1)>0) Val(1,1) = -Particles[a]->TI * Val(1,1)/(Particles[a]->Density*Particles[a]->Density); else Val(1,1) = 0.0;
					if (Val(2,2)>0) Val(2,2) = -Particles[a]->TI * Val(2,2)/(Particles[a]->Density*Particles[a]->Density); else Val(2,2) = 0.0;
					Mult(Vec,Val,temp);
					Mult(temp,VecT,Particles[a]->TIR);
				}
			}
		}


	if (SWIType < 2)
	{
		#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (size_t i=0; i<Particles.Size(); i++)
		{
			if (Particles[i]->Material == 3)
			{
				if (Particles[i]->SatCheck && !Particles[i]->IsSat)
				{
					Particles[i]->Mass		= Particles[i]->V*(Particles[i]->RefDensity - Particles[i]->RhoF);
					Particles[i]->Density		= Particles[i]->Density - Particles[i]->RhoF;
					Particles[i]->Densityb		= Particles[i]->Densityb - Particles[i]->RhoF;
					Particles[i]->RefDensity	= Particles[i]->RefDensity - Particles[i]->RhoF;
					Particles[i]->IsSat		= true;
				}
				if (!Particles[i]->SatCheck && Particles[i]->IsSat)
				{
					Particles[i]->Mass		= Particles[i]->V*(Particles[i]->RefDensity + Particles[i]->RhoF);
					Particles[i]->Density		= Particles[i]->Density + Particles[i]->RhoF;
					Particles[i]->Densityb		= Particles[i]->Densityb + Particles[i]->RhoF;
					Particles[i]->RefDensity	= Particles[i]->RefDensity + Particles[i]->RhoF;
					Particles[i]->IsSat		= false;
				}
			}
		}
	}

        //std::cout << "3 " << Particles[14966]->Sigma << std::endl;
}

inline void Domain::KGC(Particle * P1, Particle * P2)
{
    // This function obtain the kernel gradient correction matrix L for each SPH particle
    double h    = (P1->h+P2->h)/2;
    Vec3_t xij  = P1->x - P2->x;
    double rij  = norm(xij);
    double GK   = GradKernel(Dimension, KernelType, rij/h, h);
    double di=0.0,dj=0.0,mi=0.0,mj=0.0;

    if (P1->Material*P2->Material == 9){
            if (!P1->IsFree) {di = P2->Density; mi = P1->FPMassC * P2->Mass;} else {di = P1->Density;   mi = P1->Mass;}
            if (!P2->IsFree) {dj = P1->Density; mj = P2->FPMassC * P1->Mass;} else {dj = P2->Density;   mj = P2->Mass;}
    }
    else{
        if (!P1->IsFree){
            di = DensitySolid(P2->PresEq, P2->Cs, P2->P0,P1->Pressure, P2->RefDensity);
            mi = P1->FPMassC * P2->Mass;
        }
        else{
            di = P1->Density;
            mi = P1->Mass;
        }
        if (!P2->IsFree){
            dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->Pressure, P1->RefDensity);
            mj = P2->FPMassC * P1->Mass;
        }
        else{
            dj = P2->Density;
            mj = P2->Mass;
        }
    }

    switch(Dimension){
        case 2: {
            Mat_t tempCorr;
            tempCorr.Resize(Dimension,Dimension);
            tempCorr.SetValues(0.0);

            tempCorr(0,0) = xij(0)*xij(0);
            tempCorr(0,1) = xij(0)*xij(1);
            tempCorr(1,0) = xij(1)*xij(0);
            tempCorr(1,1) = xij(1)*xij(1);

            omp_set_lock(&P1->my_lock);
            P1->Lcorr += (mj/dj)*GK*tempCorr;
            omp_unset_lock(&P1->my_lock);

            omp_set_lock(&P2->my_lock);
            P2->Lcorr += (mi/di)*GK*tempCorr;
            omp_unset_lock(&P2->my_lock);
            break;
        }
        case 3: {
            Mat_t tempCorr;
            tempCorr.Resize(Dimension,Dimension);
            tempCorr.SetValues(0.0);

            tempCorr(0,0) = xij(0)*xij(0);
            tempCorr(0,1) = xij(0)*xij(1);
            tempCorr(0,2) = xij(0)*xij(2);
            tempCorr(1,0) = xij(1)*xij(0);
            tempCorr(1,1) = xij(1)*xij(1);
            tempCorr(1,2) = xij(1)*xij(2);
            tempCorr(2,0) = xij(2)*xij(0);
            tempCorr(2,1) = xij(2)*xij(1);
            tempCorr(2,2) = xij(2)*xij(2);

            omp_set_lock(&P1->my_lock);
            P1->Lcorr += (mj/dj)*GK*tempCorr;
            omp_unset_lock(&P1->my_lock);

            omp_set_lock(&P2->my_lock);
            P2->Lcorr += (mi/di)*GK*tempCorr;
            omp_unset_lock(&P2->my_lock);
            break;

        }
    }

}

inline void Domain::LastComputeAcceleration ()
{
    if (KGcorrection){
        #pragma omp parallel for schedule (static) num_threads(Nproc)
        for (size_t k=0; k<Nproc;k++){
            for (size_t i=0; i<SMPairs[k].Size();i++){
                if (Particles[SMPairs[k][i].first]->Material == 3)
                    KGC(Particles[SMPairs[k][i].first],Particles[SMPairs[k][i].second]);
            }
        }
    }

	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t k=0; k<Nproc;k++)
	{
		for (size_t i=0; i<SMPairs[k].Size();i++)
			if (Particles[SMPairs[k][i].first]->Material == 1)
				CalcForce11(Particles[SMPairs[k][i].first],Particles[SMPairs[k][i].second]);
			else
            {
                CalcForce2233(Particles[SMPairs[k][i].first],Particles[SMPairs[k][i].second]);
            }

		for (size_t i=0; i<FSMPairs[k].Size();i++)
			if (Particles[FSMPairs[k][i].first]->Material == 1)
				CalcForce11(Particles[FSMPairs[k][i].first],Particles[FSMPairs[k][i].second]);
			else
				CalcForce2233(Particles[FSMPairs[k][i].first],Particles[FSMPairs[k][i].second]);
	}

	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t k=0; k<NSMPairs.Size();k++)
	{
		for (size_t i=0; i<NSMPairs[k].Size();i++)
		{
			if (Particles[NSMPairs[k][i].first]->Material*Particles[NSMPairs[k][i].second]->Material == 3)
				CalcForce13(Particles[NSMPairs[k][i].first],Particles[NSMPairs[k][i].second]);
			else if (Particles[NSMPairs[k][i].first]->Material*Particles[NSMPairs[k][i].second]->Material == 2)
				CalcForce12(Particles[NSMPairs[k][i].first],Particles[NSMPairs[k][i].second]);
			else
			{
				std::cout<<"Out of Interaction types"<<std::endl;
				abort();
			}
		}
	}

	for (size_t i=0 ; i<Nproc ; i++)
	{
		SMPairs[i].Clear();
		FSMPairs[i].Clear();
		NSMPairs[i].Clear();
	}

		//Min time step check based on the acceleration
		double test	= 0.0;
		deltatmin	= deltatint;
		#pragma omp parallel for schedule (static) private(test) num_threads(Nproc)
		for (size_t i=0; i<Particles.Size(); i++)
		{
			if (Particles[i]->IsFree)
			{
				test = sqrt(Particles[i]->h/norm(Particles[i]->a));
				if (deltatmin > (sqrt_h_a*test))
				{
					omp_set_lock(&dom_lock);
					deltatmin = sqrt_h_a*test;
					omp_unset_lock(&dom_lock);
				}
			}
		}
}

inline void Domain::CalculateForceDEM ()
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
    {
        for (size_t j=0;j<DEMParticles.Size();j++)
        {
            if (DEM::Distance(Particles[i]->x,DEMParticles[j]->x)>Cellfac*Particles[i]->h + DEMParticles[j]->Dmax) continue;
            Vec3_t xi, xf;
            if (DEMParticles[j]->Verts.Size()==1)
            {
                double contactR = DEMParticles[j]->Props.R;
                xf = Particles[i]->x;
                xi = DEMParticles[j]->x + contactR*(xf-DEMParticles[j]->x)/norm(xf-DEMParticles[j]->x);
            }
            else
            {
                DEM::Distance(Particles[i]->x,*DEMParticles[j]->Faces[0],xf,xi);
            }
            double dist  = norm(xf-xi);
            double delta = DEMParticles[j]->Props.eps - dist;
            //std::cout << dist << " " << delta << std::endl;
            if (Particles[i]->Material == 3 || Particles[i]->Material == 2) //Interaccion between solid particles and DEM particles
            {
                if (delta > 0.0)
                {
                    Vec3_t n = (xf-xi)/dist;
                    Vec3_t tw,arm;
                    arm = xi - DEMParticles[j]->x;
                    Rotation(DEMParticles[j]->w,DEMParticles[j]->Q,tw);
                    Vec3_t vrel = Particles[i]->v - DEMParticles[j]->v - cross(tw,arm);
                    Vec3_t vn = dot(vrel,n)*n;
                    Vec3_t vt = vrel-vn;
                    Vec3_t Fn = DEMParticles[j]->Props.Kn*delta*n; // linear normal contact force SPH-DEM
                    //Vec3_t Fn = DEMParticles[j]->Props.Kn*(DEMParticles[j]->Props.eps/dist)*delta*n - DEMDamp*(Particles[i]->Mass*vn/deltat); // non-linear (Sheikh et al, 2020)
                    FrictionMap[j][i] += deltat*vt;
                    Vec3_t t = FrictionMap[j][i];
                    if (norm(t)>0.0) t/=norm(t);
                    //double Mu = tan(Particles[i]->phi);
                    double Mu = DEMParticles[j]->Props.Mu;
                    //if (norm(FrictionMap[j][i])>DEMParticles[j]->Props.Mu*norm(Fn)/DEMParticles[j]->Props.Kn)
                    if (norm(FrictionMap[j][i]) > Mu*norm(Fn)/DEMParticles[j]->Props.Kn)
                    {
                        FrictionMap[j][i] = Mu*norm(Fn)/DEMParticles[j]->Props.Kn*t;
                    }
                    //Particles[i]->a = 1.47/2.0*(Particles[i]->a - Gravity);
                    //Particles[i]->a += Fn/Particles[i]->Mass - DEMParticles[j]->Props.Kn*FrictionMap[j][i]/Particles[i]->Mass - DEMParticles[j]->Props.Gn*Particles[i]->Mass*vt + Gravity;
                    //Particles[i]->a += Fn/Particles[i]->Mass - DEMParticles[j]->Props.Kn*FrictionMap[j][i]/Particles[i]->Mass - DEMParticles[j]->Props.Gn*Particles[i]->Mass*vt;
                    //Vec3_t F = Fn - DEMParticles[j]->Props.Kn*FrictionMap[j][i] - DEMParticles[j]->Props.Gn*dot(n,vrel)*n;
                    Vec3_t F = Fn - DEMParticles[j]->Props.Kn*FrictionMap[j][i];
                    if (Dimension == 2) F(2) = 0.0;
                    Vec3_t T, Tt;
                    Vec3_t x = xi - DEMParticles[j]->x;
                    Tt = cross(x,F);
                    Quaternion_t q;
                    Conjugate(DEMParticles[j]->Q,q);
                    Rotation(Tt,q,T);
                    Particles[i]->a += F/Particles[i]->Mass;
                    omp_set_lock  (&DEMParticles[j]->lck);
                    DEMParticles[j]->F -= F;
                    DEMParticles[j]->T -= T;
                    omp_unset_lock  (&DEMParticles[j]->lck);
                }
            }
            else //Interaction for fluid particles with DEM particles
            {
                double h    = Particles[i]->h;
                if (dist < Cellfac*h)
                {
                    double mass     = Particles[i]->Mass;
                    //double rho      = mass/(dist*area);
                    //double rhoref   = mass/(l*area);
                    //double Pressure = EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0, rho, rhoref);
                    Vec3_t n = (xf-xi)/dist;
                    Vec3_t xij = xf-xi;
                    Vec3_t tw,arm;
                    arm = xi - DEMParticles[j]->x;
                    Rotation(DEMParticles[j]->w,DEMParticles[j]->Q,tw);
                    Vec3_t vrel = Particles[i]->v - DEMParticles[j]->v - cross(tw,arm);
                    Vec3_t vn = dot(vrel,n)*n;
                    Vec3_t Fn = OrthoSys::O;
                    if (delta > 0.0)  Fn = DEMParticles[j]->Props.Kn*delta*n; // linear normal contact force SPH-DEM
                    //if (delta > 0.0) Fn = DEMParticles[j]->Props.Kn*(DEMParticles[j]->Props.eps/dist)*delta*n - DEMDamp*(Particles[i]->Mass*vn/deltat); // non-linear (Sheikh et al, 2020)
		            double GK	= 2.0/(3.0*pow(h,Dimension))*GradKernel(Dimension, KernelType, dist/h, h)/Kernel(Dimension, KernelType, 0.0, h);
		            //double GK	= GradKernel(Dimension, KernelType, dist/h, h);
		            Mat3_t StrainRate;
		            set_to_zero(StrainRate);
		            double Mu = Particles[i]->Mu;
		            if ((Particles[i]->T0>1.0e-12) || (Particles[i]->LES)){
						StrainRate =	2.0*vrel(0)*xij(0) , vrel(0)*xij(1)+vrel(1)*xij(0) , vrel(0)*xij(2)+vrel(2)*xij(0) ,
										vrel(0)*xij(1)+vrel(1)*xij(0) , 2.0*vrel(1)*xij(1) , vrel(1)*xij(2)+vrel(2)*xij(1) ,
						 				vrel(0)*xij(2)+vrel(2)*xij(0) , vrel(1)*xij(2)+vrel(2)*xij(1) , 2.0*vrel(2)*xij(2) ;
						StrainRate = 	-GK*StrainRate;
						Particles[i]->StrainRate		= Particles[i]->StrainRate + mass/Particles[i]->Density*StrainRate;
					}

		            Vec3_t VI = 0.0;
		            Viscous_Force(VisEq, VI, Mu, Particles[i]->Density, Particles[i]->Density, GK, vrel, Dimension, KernelType, dist, h, xij, vrel);
                    Vec3_t F = Fn + mass*mass*VI;
                    if (Dimension == 2) F(2) = 0.0;
                    Vec3_t T, Tt;
                    Vec3_t x = xi - DEMParticles[j]->x;
                    Tt = cross(x,F);
                    Quaternion_t q;
                    Conjugate(DEMParticles[j]->Q,q);
                    Rotation(Tt,q,T);

                    Particles[i]->a += F/mass;
                    //Particles[i]->a += mass*(-2.0*Particles[i]->Pressure/(Particles[i]->Density*Particles[i]->Density)*GK*xij+VI);
                    //Particles[i]->dDensity += mass*dot(vrel,GK*xij);
                    omp_set_lock  (&DEMParticles[j]->lck);
                    DEMParticles[j]->F -= F;
                    DEMParticles[j]->T -= T;
                    omp_unset_lock  (&DEMParticles[j]->lck);
                }
            }
        }
    }

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0; i<DEMInteractons.Size(); i++)
    {
	    if (DEMInteractons[i]->CalcForce(deltat))
        {
            std::cout << "Maximun overlap detected between particles at time " << Time << std::endl;
            sleep(1);
            throw new Fatal("Maximun overlap detected between particles");
        }
        omp_set_lock  (&DEMInteractons[i]->P1->lck);
        DEMInteractons[i]->P1->F += DEMInteractons[i]->F1;
        DEMInteractons[i]->P1->T += DEMInteractons[i]->T1;
        omp_unset_lock(&DEMInteractons[i]->P1->lck);
        omp_set_lock  (&DEMInteractons[i]->P2->lck);
        DEMInteractons[i]->P2->F += DEMInteractons[i]->F2;
        DEMInteractons[i]->P2->T += DEMInteractons[i]->T2;
        omp_unset_lock(&DEMInteractons[i]->P2->lck);
    }
}

inline void Domain::MoveDEM (double dt)
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0; i<DEMParticles.Size(); i++)
    {
        if (Dimension == 2)
        {
            DEMParticles[i]->F(2) = 0.0;
            DEMParticles[i]->T(0) = 0.0;
            DEMParticles[i]->T(1) = 0.0;
        }
        //std::cout << "1" << std::endl;
        DEMParticles[i]->Translate(dt);
        //std::cout << "2" << std::endl;
        DEMParticles[i]->Rotate(dt);
    }
}

inline void Domain::HomoAcc()
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t ng=0;ng<HomoGroups.Size();ng++)
    {
        Vec3_t acc  = OrthoSys::O;
        double mass = 0.0;
        for (size_t np=0;np<HomoGroups[ng].Size();np++)
        {
            acc  += Particles[HomoGroups[ng][np]]->a*Particles[HomoGroups[ng][np]]->Mass;
            mass += Particles[HomoGroups[ng][np]]->Mass;
        }
        acc /= mass;

        for (size_t np=0;np<HomoGroups[ng].Size();np++)
        {
            Particles[HomoGroups[ng][np]]->a = acc;
        }
    }
}

inline void Domain::Move (double dt)
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
		if (Particles[i]->IsFree)
		{
			if (Particles[i]->InOut>0)
			{
				Particles[i]->a = 0.0;
				if (Particles[i]->InOut == 1)
				{
					Particles[i]->dDensity = 0.0;
					Particles[i]->ZWab = 0.0;
				}
				else
				{
					if (BC.outDensity>0.0)
					{
						Particles[i]->dDensity = 0.0;
						Particles[i]->ZWab = 0.0;
					}
				}
			}
		    Particles[i]->Move(dt,DomSize,TRPR,BLPF,Scheme,I);
		}
}

inline void Domain::InFlowBCLeave()
{
	size_t a,b;
	Array <int> DelPart,TempPart;
	Array<std::pair<Vec3_t,size_t> > AddPart;

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
		if ((Particles[i]->x(0) > TRPR(0)) || (Particles[i]->x(1) > TRPR(1)) || (Particles[i]->x(2) > TRPR(2)) ||
				(Particles[i]->x(0) < BLPF(0)) || (Particles[i]->x(1) < BLPF(1)) || (Particles[i]->x(2) < BLPF(2)))
		{
			Particles[i]->InOut	= 0;
			omp_set_lock(&dom_lock);
			DelPart.Push(i);
			omp_unset_lock(&dom_lock);
		}

	if (BC.InOutFlow==1 || BC.InOutFlow==3)
	{
		for (size_t i=0 ; i<BC.InPart.Size() ; i++)
			if(Particles[BC.InPart[i]]->x(0) > BC.InFlowLoc1)
			{
				Vec3_t temp1	 = Particles[BC.InPart[i]]->x;
				temp1(0)	-= (BC.InFlowLoc3-BC.InFlowLoc2+InitialDist);
				Particles[BC.InPart[i]]->InOut		= 0;
				Particles[BC.InPart[i]]->FirstStep	= false;
				Particles[BC.InPart[i]]->ShepardCounter	= 0;
				AddPart.Push(std::make_pair(temp1,BC.InPart[i]));
				TempPart.Push(i);
			}
		BC.InPart.DelItems(TempPart);
		TempPart.Clear();
	}

	if (AddPart.Size() >= DelPart.Size())
	{
		for (size_t i=0 ; i<DelPart.Size() ; i++)
		{
			a = DelPart[i];
			b = AddPart[i].second;
			Particles[a]->x 		= AddPart[i].first;
			Particles[a]->Material		= 1;
			Particles[a]->InOut		= 1;
			Particles[a]->FirstStep		= false;

			Particles[a]->P0		= Particles[b]->P0;
			Particles[a]->PresEq		= Particles[b]->PresEq;
			Particles[a]->Cs		= Particles[b]->Cs;

			Particles[a]->Alpha		= Particles[b]->Alpha;
			Particles[a]->Beta		= Particles[b]->Beta;
			Particles[a]->Mu		= Particles[b]->Mu;
			Particles[a]->MuRef		= Particles[b]->MuRef;
			Particles[a]->T0		= 0.0; // Inflow is not capable of injecting non-Newtonian fluid

			Particles[a]->Mass 		= Particles[b]->Mass;
			Particles[a]->h			= Particles[b]->h;

			Particles[a]->ID 		= Particles[b]->ID;

			Particles[a]->TI		= 0.0; // Inflow is not capable of condidering the tensile instability

			Particles[a]->RefDensity	= Particles[b]->RefDensity; // The density for inflow must always be defined

			Particles[a]->ct		= Particles[b]->ct;

			Particles[a]->Shepard		= Particles[b]->Shepard;
			Particles[a]->ShepardStep	= Particles[b]->ShepardStep;
			Particles[a]->ShepardCounter	= Particles[b]->ShepardCounter;

			Particles[a]->LES		= Particles[b]->LES;
			Particles[a]->CSmag		= Particles[b]->CSmag;
			BC.InPart.Push(a);
		}

		if (AddPart.Size() != DelPart.Size())
			for (size_t i=DelPart.Size() ; i<AddPart.Size() ; i++)
			{
				b = AddPart[i].second;

				Particles.Push(new Particle(Particles[b]->ID,AddPart[i].first,Particles[b]->v,Particles[b]->Mass,Particles[b]->RefDensity,Particles[b]->h,false));

				a = Particles.Size()-1;
				Particles[a]->Material		= 1;
				Particles[a]->InOut		= 1;
				Particles[a]->FirstStep		= false;

				Particles[a]->P0		= Particles[b]->P0;
				Particles[a]->PresEq		= Particles[b]->PresEq;
				Particles[a]->Cs		= Particles[b]->Cs;

				Particles[a]->Alpha		= Particles[b]->Alpha;
				Particles[a]->Beta		= Particles[b]->Beta;
				Particles[a]->Mu		= Particles[b]->Mu;
				Particles[a]->MuRef		= Particles[b]->MuRef;
				Particles[a]->T0		= 0.0; // Inflow is not capable of injecting non-Newtonian fluid

				Particles[a]->TI		= 0.0; // Inflow is not capable of condidering the tensile instability

				Particles[a]->ct		= Particles[b]->ct;

				Particles[a]->Shepard		= Particles[b]->Shepard;
				Particles[a]->ShepardStep	= Particles[b]->ShepardStep;
				Particles[a]->ShepardCounter	= Particles[b]->ShepardCounter;

				Particles[a]->LES		= Particles[b]->LES;
				Particles[a]->CSmag		= Particles[b]->CSmag;
				BC.InPart.Push(a);
			}

		DelPart.Clear();
		AddPart.Clear();
	}
	else
	{
		for (size_t i=0 ; i<AddPart.Size() ; i++)
		{
			a = DelPart[i];
			b = AddPart[i].second;
			Particles[a]->x 		= AddPart[i].first;
			Particles[a]->Material		= 1;
			Particles[a]->InOut		= 1;
			Particles[a]->FirstStep		= false;

			Particles[a]->P0		= Particles[b]->P0;
			Particles[a]->PresEq		= Particles[b]->PresEq;
			Particles[a]->Cs		= Particles[b]->Cs;

			Particles[a]->Alpha		= Particles[b]->Alpha;
			Particles[a]->Beta		= Particles[b]->Beta;
			Particles[a]->Mu		= Particles[b]->Mu;
			Particles[a]->MuRef		= Particles[b]->MuRef;
			Particles[a]->T0		= 0.0; // Inflow is not capable of injecting non-Newtonian fluid

			Particles[a]->Mass 		= Particles[b]->Mass;
			Particles[a]->h			= Particles[b]->h;

			Particles[a]->ID 		= Particles[b]->ID;

			Particles[a]->TI		= 0.0; // Inflow is not capable of condidering the tensile instability

			Particles[a]->RefDensity	= Particles[b]->RefDensity; // The density for inflow must always be defined

			Particles[a]->ct		= Particles[b]->ct;

			Particles[a]->Shepard		= Particles[b]->Shepard;
			Particles[a]->ShepardStep	= Particles[b]->ShepardStep;
			Particles[a]->ShepardCounter	= Particles[b]->ShepardCounter;

			Particles[a]->LES		= Particles[b]->LES;
			Particles[a]->CSmag		= Particles[b]->CSmag;
			BC.InPart.Push(a);
		}
		for (size_t i=AddPart.Size() ; i<DelPart.Size() ; i++)
		{
			TempPart.Push(DelPart[i]);
		}
		Particles.DelItems(TempPart);
		BC.inoutcounter = 1;
		DelPart.Clear();
		AddPart.Clear();
	}
}

inline void Domain::InFlowBCFresh()
{
	int temp, temp1;
	int q1,q2,q3;
	if (BC.inoutcounter == 0)
	{
		if (BC.InOutFlow==1 || BC.InOutFlow==3)
		{
			if (!(norm(BC.inv)>0.0) || !(BC.inDensity>0.0))
			{
				std::cout<< "BC.inv or BC.inDensity are not defined, please define them and run the code again."<<std::endl;
				abort();
			}

			BC.InPart.Clear();
			BC.InFlowLoc1  = BLPF(0) + BC.cellfac*hmax;
			temp1 = (int) (floor((BC.InFlowLoc1 - BLPF(0)) / CellSize(0)));

			for (q2=0; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
			for (q3=0; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
			for (q1=0; q1<(temp1 + 1)                                      ; q1++)
			{
				if (HOC[q1][q2][q3]!=-1)
				{
					temp = HOC[q1][q2][q3];
					while (temp != -1)
					{
						if (Particles[temp]->IsFree && (Particles[temp]->x(0) <= BC.InFlowLoc1) )
						{
							BC.InPart.Push(temp);
							Particles[temp]->InOut = 1;
						}
						temp = Particles[temp]->LL;
					}
				}
			}
			BC.InFlowLoc2  = Particles[BC.InPart[0]]->x(0);
			BC.InFlowLoc3  = Particles[BC.InPart[0]]->x(0);
			#pragma omp parallel for schedule(static) num_threads(Nproc)
			for (size_t i=0 ; i<BC.InPart.Size() ; i++)
			{
				if (Particles[BC.InPart[i]]->x(0) < BC.InFlowLoc2) BC.InFlowLoc2  = Particles[BC.InPart[i]]->x(0);
				if (Particles[BC.InPart[i]]->x(0) > BC.InFlowLoc3) BC.InFlowLoc3  = Particles[BC.InPart[i]]->x(0);
			}
		}

		if (BC.InOutFlow==2 || BC.InOutFlow==3)
			BC.OutFlowLoc = TRPR(0) - BC.cellfac*hmax;

		BC.inoutcounter = 2;
	}

	if (BC.inoutcounter == 1)
	{
		BC.InPart.Clear();
		temp1 = (int) (floor((BC.InFlowLoc1 - BLPF(0)) / CellSize(0)));

		for (q2=0; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
		for (q3=0; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
		for (q1=0; q1<(temp1 + 1)                                      ; q1++)
		{
			if (HOC[q1][q2][q3]!=-1)
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree && (Particles[temp]->x(0) <= BC.InFlowLoc1) && Particles[temp]->InOut==1)
						BC.InPart.Push(temp);
					temp = Particles[temp]->LL;
				}
			}
		}
		BC.inoutcounter = 2;
	}


	if (BC.InOutFlow==2 || BC.InOutFlow==3)
	{
		BC.OutPart.Clear();
		temp1 = (int) (floor((BC.OutFlowLoc - BLPF(0)) / CellSize(0)));

		for (q2=0     ; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
		for (q3=0     ; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
		for (q1=temp1 ; q1<CellNo[0]                                        ; q1++)
		{
			if (HOC[q1][q2][q3]!=-1)
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree && (Particles[temp]->x(0) >= BC.OutFlowLoc) )
					{
						BC.OutPart.Push(temp);
						Particles[temp]->InOut = 2;
					}
					temp = Particles[temp]->LL;
				}
			}
		}
	}

	Vec3_t vel = 0.0;
	double den = 0.0;
	if (BC.InPart.Size()>0)
		#pragma omp parallel for schedule(static) private(vel,den) num_threads(Nproc)
		for (size_t i=0 ; i<BC.InPart.Size() ; i++)
		{
			size_t a = BC.InPart[i];
			InCon(Particles[a]->x,vel,den,BC);
			Particles[a]->v  = vel;
			Particles[a]->vb = vel;
			Particles[a]->va = vel;
			Particles[a]->Density  = den;
			Particles[a]->Densityb = den;
			Particles[a]->Densitya = den;
    		Particles[a]->Pressure = EOS(Particles[a]->PresEq, Particles[a]->Cs, Particles[a]->P0,Particles[a]->Density, Particles[a]->RefDensity);
		}

	double temp11;
	if (BC.MassConservation)
		temp11 = BC.InPart.Size()*1.0/(BC.OutPart.Size()*1.0);
	else
		temp11 = 1.0;

	if (BC.OutPart.Size()>0)
		#pragma omp parallel for schedule(static) private(vel,den) num_threads(Nproc)
		for (size_t i=0 ; i<BC.OutPart.Size() ; i++)
		{
			size_t a = BC.OutPart[i];
			OutCon(Particles[a]->x,vel,den,BC);
			if (norm(BC.outv)>0.0 || temp11 != 1.0)
			{
				Particles[a]->v  = temp11*vel;
				Particles[a]->vb = temp11*vel;
				Particles[a]->va = temp11*vel;
			}
			if (BC.outDensity>0.0)
			{
				Particles[a]->Density  = den;
				Particles[a]->Densityb = den;
 				Particles[a]->Densitya = den;
    				Particles[a]->Pressure = EOS(Particles[a]->PresEq, Particles[a]->Cs, Particles[a]->P0,Particles[a]->Density, Particles[a]->RefDensity);
			}
		}
}

inline void Domain::WholeVelocity()
{
    //Apply a constant velocity to all particles in the initial time step
    if (norm(BC.allv)>0.0 || BC.allDensity>0.0)
    {
    	Vec3_t vel = 0.0;
    	double den = 0.0;

	#pragma omp parallel for schedule (static) private(vel,den) num_threads(Nproc)
    	for (size_t i=0 ; i<Particles.Size() ; i++)
    	{
	    	AllCon(Particles[i]->x,vel,den,BC);
    		if (Particles[i]->IsFree && norm(BC.allv)>0.0 && Particles[i]->Material == 1)
    		{
		    	Particles[i]->v		= vel;
 	    	}
    		if (Particles[i]->IsFree && BC.allDensity>0.0 && Particles[i]->Material == 1)
    		{
		    	Particles[i]->Density	= den;
		    	Particles[i]->Pressure	= EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0,Particles[i]->Density, Particles[i]->RefDensity);
    		}
    	}
    }
}

inline void Domain::InitialChecks()
{
	//initializing identity matrix
	if (Dimension == 2) I(2,2) = 0;

	if (Dimension<=1 || Dimension>3)
	{
		std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
	}

	if (BC.InOutFlow>0 && BC.Periodic[0])
		throw new Fatal("Periodic BC in the X direction cannot be used with In/Out-Flow BC simultaneously");

    //Setting DEM eps distance to zero initially
    FrictionMap.Resize(DEMParticles.Size());
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<DEMParticles.Size(); i++)
    {
        DEMParticles[i]->Props.Kn  = 0.0;
        DEMParticles[i]->InitializeVelocity(deltat);
        FrictionMap[i].Resize(Particles.Size());
    }

	for (size_t i=0; i<DEMParticles.Size(); i++)
	//#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (size_t i=0; i<Particles.Size(); i++)
	{

		        
        //Initialising parameters for Non-Newtonian fluids
		if ((Particles[i]->Material == 1) && (Particles[i]->T0>1.0e-12)){
			Particles[i]->Mumax = Particles[i]->MuRef + Particles[i]->T0*Particles[i]->m;
        	Particles[i]->gamma0 = Particles[i]->T0/(Particles[i]->Mumax-Particles[i]->MuRef);
            Particles[i]->Mu0=1.0e3*Particles[i]->MuRef;  
            Particles[i]->Kcross = Particles[i]->Mu0/Particles[i]->T0;}
        	//std::cout<<"numax "<<Particles[i]->Mumax<<std::endl;


		//Initialising pressure of solid and fluid particles
		if (Particles[i]->Material < 3)
			Particles[i]->Pressure = EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0,Particles[i]->Density, Particles[i]->RefDensity);

		// Initialising the permeability for soil particles
		if (Particles[i]->Material == 3)
		{
			switch(Particles[i]->SeepageType)
			{
				case 0:
					break;

				case 1:
					Particles[i]->k = Particles[i]->n0*Particles[i]->n0*Particles[i]->n0*Particles[i]->d*Particles[i]->d/(180.0*(1.0-Particles[i]->n0)*(1.0-Particles[i]->n0));
					break;

				case 2:
					Particles[i]->k = Particles[i]->n0*Particles[i]->n0*Particles[i]->n0*Particles[i]->d*Particles[i]->d/(150.0*(1.0-Particles[i]->n0)*(1.0-Particles[i]->n0));
					Particles[i]->k2= 1.75*(1.0-Particles[i]->n0)/(Particles[i]->n0*Particles[i]->n0*Particles[i]->n0*Particles[i]->d);
					break;

				case 3:
					Particles[i]->k = Particles[i]->n0*Particles[i]->n0*Particles[i]->n0*Particles[i]->d*Particles[i]->d/(150.0*(1.0-Particles[i]->n0)*(1.0-Particles[i]->n0));
					Particles[i]->k2= 0.4/(Particles[i]->n0*Particles[i]->n0*Particles[i]->d);
					break;

				default:
					std::cout << "Seepage Type No is out of range. Please correct it and run again" << std::endl;
					std::cout << "0 => Darcy's Law" << std::endl;
					std::cout << "1 => Darcy's Law & KozenyCarman Eq" << std::endl;
					std::cout << "2 => The Forchheimer Eq & Ergun Coeffs" << std::endl;
					std::cout << "3 => The Forchheimer Eq & Den Adel Coeffs" << std::endl;
					abort();
					break;
			}

			Particles[i]->n = Particles[i]->n0;

		}
        //This part updates eps so the minimun distance at the beginning is equal to eps
        for (size_t j = 0;j<DEMParticles.Size();j++){
            double Kn = DEMstiff*Particles[i]->Mass/(deltat*deltat);
            if (Kn<DEMParticles[j]->Props.Kn || DEMParticles[j]->Props.Kn == 0.0) DEMParticles[j]->Props.Kn = Kn;
        }
	}

    //Creating interacton array for DEM particles
	for (size_t i=0  ; i<DEMParticles.Size()-1; i++)
	for (size_t j=i+1; j<DEMParticles.Size()  ; j++)
    {
        if (!DEMParticles[i]->IsFree()&&!DEMParticles[j]->IsFree()) continue;
        if (DEMParticles[i]->Verts.Size()==1 && DEMParticles[j]->Verts.Size()==1)
        {
            DEMInteractons.Push (new DEM::CInteractonSphere(DEMParticles[i],DEMParticles[j]));
        }
        else
        {
            DEMInteractons.Push (new DEM::CInteracton(DEMParticles[i],DEMParticles[j]));
        }
    }

    for (size_t i=0; i<DEMInteractons.Size(); i++)
    {
        DEMInteractons[i]->UpdateContacts(1.0e12);
    }
        
}

inline void Domain::TimestepCheck ()
{
	// Check the time step
	double t1,t2;
	t1 = 0.25*hmax/(CsMax);
	if (MuMax>0.0) t2 = 0.125*hmax*hmax*rhomax/MuMax; else t2 =1000000.0;

	std::cout << "Max allowable time step using CFL = "<< std::min(t1,t2) << " S" << std::endl;
	std::cout << "User Time Step = "<< deltatint  << " S" << std::endl;

	if (deltatint > std::min(t1,t2))
	throw new Fatal("Please decrease the time step to the allowable range");
}

inline void Domain::Solve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx)
{
	size_t idx_out = 1;
	double tout = Time;

	//Initializing adaptive time step variables
	deltat = deltatint = deltatmin	= dt;

	InitialChecks();
	CellInitiate();
	ListGenerate();
	PrintInput(TheFileKey);
	TimestepCheck();
	WholeVelocity();


	//Initial model output
	if (TheFileKey!=NULL)
	{
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\n--------------Initial Condition has been generated--------------\n" << std::endl;
	}


	std::cout <<           " " << std::endl;
	std::cout << BOLD(FGRN("----------------------------------------------------------------------")) << std::endl;
	std::cout << BOLD(FGRN("|                               Solving                              |")) << std::endl; // Print the local time
	std::cout << BOLD(FGRN("----------------------------------------------------------------------")) << std::endl;
	std::cout <<           " " << std::endl;

	auto start = high_resolution_clock::now(); // Get starting timepoint 

	while (Time<tf && idx_out<=maxidx)
	{
		StartAcceleration(Gravity);
		if (BC.InOutFlow>0) InFlowBCFresh();
		MainNeighbourSearch();
		GeneralBefore(*this); 
		PrimaryComputeAcceleration();
		LastComputeAcceleration();
        CalculateForceDEM();
		GeneralAfter(*this);

		// output
		if (Time>=tout)
		{
			if (TheFileKey!=NULL)
			{
				String fn;
				fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				WriteXDMF    (fn.CStr());
				std::cout << "\nOutput No. " << idx_out << ", at t = " << Time << "s has been generated" << std::endl;
				std::cout << "Current time-step, dt = " <<deltat<<std::endl;
			}
			idx_out++;
			tout += dtOut;
		}

		AdaptiveTimeStep();
		HomoAcc();
		Move(deltat);
        MoveDEM(deltat);
		Time += deltat;
		if (BC.InOutFlow>0) InFlowBCLeave(); else CheckParticleLeave ();
		CellReset();
		ListGenerate();

	}
	
	// This part obtain the total computational cost of the simulation
	auto stop = high_resolution_clock::now();  // Get ending timepoint 
	auto duration = duration_cast<seconds>(stop - start); // Estimating the duration of the computational time

	std::cout << " " << std::endl;
  	std::cout << "++++++++++++++ Total computational time: " << duration.count() << " seconds ++++++++++++++" << std::endl;
	std::cout << " " << std::endl;

	std::cout <<           " " << std::endl;
	std::cout << BOLD(FGRN("----------------------------------------------------------------------")) << std::endl;
	std::cout <<           "|                The simulation finished successfully                |"  << std::endl;
	std::cout << BOLD(FGRN("----------------------------------------------------------------------")) << std::endl;
	std::cout <<           "\n" << std::endl;

}

inline void Domain::PrintInput(char const * FileKey)
{
	//type definition to shorten coding
	std::ostringstream oss;

    //Writing Inputs in a Log file
	String fn(FileKey);
    oss << "___________________________________________________________";
    oss << "\nComputational setup:\n";
    oss << "___________________________________________________________\n";
    
	oss << "\nDimension = "<< Dimension << "D\n";

    oss << "\nTime integration scheme = ";
    switch (Scheme)
    {
        case 0:
            oss << "Modified Verlet\n";
            break;
        case 1:
            oss << "Leapfrog\n";
            break;
    }
    
	oss << "\nKernel Type = ";
	switch (KernelType)
	{
		case 0:
		    oss << "Qubic Spline\n";
	    	break;
		case 1:
		    oss << "Quintic\n";
		    break;
		case 2:
		    oss << "Quintic Spline\n";
		    break;
	}

	oss << "\nViscosity Equation = ";
	switch (VisEq)
	{
		case 0:
			oss << "0 => Morris et al 1997\n";
			break;
		case 1:
			oss << "1 => Shao et al 2003\n";
			break;
		case 2:
			oss << "2 => Real viscosity for incompressible fluids\n";
			break;
		case 3:
			oss << "3 => Takeda et al 1994 (Real viscosity for compressible fluids)\n";
			break;
	}

    oss << "\nSoil-water interaction = ";
    switch (SWIType)
    {
        case 0:
            oss << "0 => Seepage + buoyant\n";
            break;
        case 1:
            oss << "1 => Seepage + surface erosion + buoyant\n";
            break;
        case 2:
            oss << "2 => Seepage + pore water pressure\n";
            break;
        case 3:
            oss << "3 => No fluid-soil interaction\n";
            break;
    }
    
	oss << "\nComputational domain size:\n";
	oss << "Bottom Left-Corner Front = " << BLPF <<" m\n";
	oss << "Top Right-Corner Rear    = " << TRPR <<" m\n";

	oss << "\nMax of the smoothing lengths, h = " << hmax << " m\n";
	oss << "Cell factor in Linked List (based on kernels) = " << Cellfac << "\n";

	oss << "\nCell Size in XYZ Directions = " << CellSize <<" m\n";
	oss << "No of Cells in XYZ Directions = ( " << CellNo[0] << " , " << CellNo[1] << " , " << CellNo[2] <<" )\n" ;

	oss << "\nInitial number of particles = " << Particles.Size() << "\n";

	oss << "\nInitial time-step = "<< deltatint << " s \n";
	oss << "\nInitial dx = "<< InitialDist << " m \n";

	oss << "\nExternal acceleration (Gravity enabled)= "<<Gravity<< " m/s2\n";

	oss << "\nNumber of threads = "<<Nproc<<"\n";

	oss << "\nPeriodic Boundary Condition X dir= " << (BC.Periodic[0] ? "True" : "False") << "\n";
	oss << "Periodic Boundary Condition Y dir= " << (BC.Periodic[1] ? "True" : "False") << "\n";
	oss << "Periodic Boundary Condition Z dir= " << (BC.Periodic[2] ? "True" : "False") << "\n";

	oss << "\nGroups with an homogeneous acceleration: "<<HomoGroups.Size()<< "\n";

    oss << "\n___________________________________________________________";
    oss << "\nParametric setup:";
    oss << "\n___________________________________________________________\n";

    oss << "\nTo see the physical parameters check the input file '  .inp'\n";
    
    //oss << "\nViscosity " << Particles[0]->Mu << "\n";

	fn = FileKey;
	fn.append("_log.dat");
	std::ofstream of(fn.CStr(), std::ios::out);
	of << oss.str();
	of.close();
}

inline void Domain::WriteXDMF (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);



    
    //Storing DEM data
    size_t N_Faces = 0;
    size_t N_Verts = 0;
    if (DEMParticles.Size()>0)
    {
        for (size_t i=0; i<DEMParticles.Size(); i++) 
        { 
            for (size_t j=0;j<DEMParticles[i]->Faces.Size();j++)
            {
                N_Faces += DEMParticles[i]->Faces[j]->Edges.Size();
            }
            N_Verts += DEMParticles[i]->Verts.Size() + DEMParticles[i]->Faces.Size();
        }

        if (N_Faces>0)
        {

            //Geometric information
            float  * Verts   = new float [3*N_Verts];
            int    * FaceCon = new int   [3*N_Faces];
            
            //Atributes
            int    * Tags    = new int   [  N_Faces];

            size_t n_verts = 0;
            size_t n_faces = 0;
            size_t n_attrs = 0;
            for (size_t i=0;i<DEMParticles.Size();i++)
            {
                DEM::Particle * Pa = DEMParticles[i];
                size_t n_refv = n_verts/3;
                Array<Vec3_t> Vtemp(Pa->Verts.Size());
                Array<Vec3_t> Vres (Pa->Verts.Size());
                for (size_t j=0;j<Pa->Verts.Size();j++)
                {
                    Vtemp[j] = *Pa->Verts[j];
                    Vres [j] = *Pa->Verts[j];
                }
                double multiplier = 0.0;
                //if (Dilate&&Pa->Eroded&&Pa->Faces.Size()>=4)
                //{
                    //DEM::Dilation(Vtemp,Pa->EdgeCon,Pa->FaceCon,Vres,Pa->Props.R);
                    //multiplier = 1.0;
                //}
                for (size_t j=0;j<Pa->Verts.Size();j++)
                {
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
                    for (size_t k=0;k<Pa->FaceCon[j].Size();k++)
                    {
                        size_t nin = Pa->FaceCon[j][k];
                        size_t nen = Pa->FaceCon[j][(k+1)%Pa->FaceCon[j].Size()];
                        FaceCon[n_faces++] = int(n_reff + j);  
                        FaceCon[n_faces++] = int(n_refv + nin);
                        FaceCon[n_faces++] = int(n_refv + nen);

                        Tags  [n_attrs] = int(Pa->Tag);
                        n_attrs++;
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
            dsname.Printf("DFTag");
            H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tags   );
            
            //Erasing the data
            delete [] Verts;
            delete [] FaceCon;
            delete [] Tags;
        }

        float * DEMRadius = new float[  DEMParticles.Size()];
        float * DEMPosvec = new float[3*DEMParticles.Size()];
        float * DEMVelvec = new float[3*DEMParticles.Size()];
        float * DEMOmevec = new float[3*DEMParticles.Size()];
        float * DEMForvec = new float[3*DEMParticles.Size()];
        float * DEMTorvec = new float[3*DEMParticles.Size()];
        int   * DEMTag    = new int  [  DEMParticles.Size()];

        for (size_t i=0;i<DEMParticles.Size();i++)
        {
            Vec3_t Ome,Tor;
            Rotation(DEMParticles[i]->w,DEMParticles[i]->Q,Ome);
            Rotation(DEMParticles[i]->T,DEMParticles[i]->Q,Tor);
            DEMParticles[i]->Verts.Size()==1 ? DEMRadius[i] = float(DEMParticles[i]->Dmax) : DEMRadius[i] = 0.0;
            DEMPosvec[3*i  ] = float(DEMParticles[i]->x(0));
            DEMPosvec[3*i+1] = float(DEMParticles[i]->x(1));
            DEMPosvec[3*i+2] = float(DEMParticles[i]->x(2));
            DEMVelvec[3*i  ] = float(DEMParticles[i]->v(0));
            DEMVelvec[3*i+1] = float(DEMParticles[i]->v(1));
            DEMVelvec[3*i+2] = float(DEMParticles[i]->v(2));
            DEMOmevec[3*i  ] = float(Ome(0));
            DEMOmevec[3*i+1] = float(Ome(1)); 
            DEMOmevec[3*i+2] = float(Ome(2)); 
            DEMForvec[3*i  ] = float(DEMParticles[i]->F(0));
            DEMForvec[3*i+1] = float(DEMParticles[i]->F(1));
            DEMForvec[3*i+2] = float(DEMParticles[i]->F(2));
            DEMTorvec[3*i  ] = float(Tor(0));
            DEMTorvec[3*i+1] = float(Tor(1)); 
            DEMTorvec[3*i+2] = float(Tor(2)); 
            DEMTag   [i]     = int  (DEMParticles[i]->Tag);  
        }

        hsize_t dims[1];
        dims[0] = 3*DEMParticles.Size();
        String dsname;
        dsname.Printf("DEMPosition");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,DEMPosvec);
        dsname.Printf("DEMVelocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,DEMVelvec);
        dsname.Printf("DEMAngVelocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,DEMOmevec);
        dsname.Printf("DEMForce");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,DEMForvec);
        dsname.Printf("DEMTorque");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,DEMTorvec);
        dims[0] = DEMParticles.Size();
        dsname.Printf("DEMRadius");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,DEMRadius);
        dsname.Printf("DEMTag");
        H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,DEMTag   );

        delete [] DEMRadius;
        delete [] DEMPosvec;
        delete [] DEMVelvec;
        delete [] DEMOmevec;
        delete [] DEMForvec;
        delete [] DEMTorvec;
        delete [] DEMTag;
    }


    //Storing SPH data
    float * Posvec	= new float[3*Particles.Size()];
    float * Velvec	= new float[3*Particles.Size()];
    float * ACCvec	= new float[3*Particles.Size()];
    float * Pressure= new float[  Particles.Size()];
    float * Prop1	= new float[  Particles.Size()];
    float * Prop2	= new float[  Particles.Size()];
    float * Prop3	= new float[  Particles.Size()];
    float * Density	= new float[  Particles.Size()];
    float * Mass	= new float[  Particles.Size()];
    float * sh		= new float[  Particles.Size()];
    int   * Tag		= new int  [  Particles.Size()];
    float * Sigma	= new float[6*Particles.Size()];
    float * Strain	= new float[6*Particles.Size()];

	double P1,P2,P3;

    #pragma omp parallel for schedule (static) private(P1,P2,P3) num_threads(Nproc)
    for (size_t i=0;i<Particles.Size();i++)
    {
        Posvec  [3*i  ] = float(Particles[i]->x(0));
        Posvec  [3*i+1] = float(Particles[i]->x(1));
        Posvec  [3*i+2] = float(Particles[i]->x(2));
        Velvec  [3*i  ] = float(Particles[i]->v(0));
        Velvec  [3*i+1] = float(Particles[i]->v(1));
        Velvec  [3*i+2] = float(Particles[i]->v(2));
        ACCvec  [3*i  ] = float(Particles[i]->a(0));
        ACCvec  [3*i+1] = float(Particles[i]->a(1));
        ACCvec  [3*i+2] = float(Particles[i]->a(2));
       	Pressure[i    ] = float(Particles[i]->Pressure);
        Density [i    ] = float(Particles[i]->Density);
        Mass	[i    ] = float(Particles[i]->Mass);
        sh	    [i    ] = float(Particles[i]->h);
        Tag     [i    ] = int  (Particles[i]->ID);
        Sigma   [6*i  ] = float(Particles[i]->Sigma(0,0));
        Sigma   [6*i+1] = float(Particles[i]->Sigma(0,1));
        Sigma   [6*i+2] = float(Particles[i]->Sigma(0,2));
        Sigma   [6*i+3] = float(Particles[i]->Sigma(1,1));
        Sigma   [6*i+4] = float(Particles[i]->Sigma(1,2));
        Sigma   [6*i+5] = float(Particles[i]->Sigma(2,2));
        Strain  [6*i  ] = float(Particles[i]->Strain(0,0));
        Strain  [6*i+1] = float(Particles[i]->Strain(0,1));
        Strain  [6*i+2] = float(Particles[i]->Strain(0,2));
        Strain  [6*i+3] = float(Particles[i]->Strain(1,1));
        Strain  [6*i+4] = float(Particles[i]->Strain(1,2));
        Strain  [6*i+5] = float(Particles[i]->Strain(2,2));

	    UserOutput(Particles[i],P1,P2,P3);
        Prop1	[i    ] = float(P1);
        Prop2	[i    ] = float(P2);
        Prop3	[i    ] = float(P3);
   }
    
    int data[1];
    String dsname;
    hsize_t dims[1];
    dims[0]=1;
    data[0]=Particles.Size();
    dsname.Printf("/NP");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,data);
    dims[0] = 3*Particles.Size();
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("Velocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dsname.Printf("Acceleration");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,ACCvec);
    dims[0] = Particles.Size();
    dsname.Printf("Tag");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tag);
    dsname.Printf("Pressure");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Pressure);
    dsname.Printf("Density");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density);
    dsname.Printf(OutputName[0]);
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop1);
    dsname.Printf(OutputName[1]);
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop2);
    dsname.Printf(OutputName[2]);
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Prop3);
    dsname.Printf("Mass");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Mass);
    dsname.Printf("h");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,sh);
    dims[0] = 6*Particles.Size();
    dsname.Printf("Sigma");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma);
    dsname.Printf("Strain");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Strain);



    delete [] Posvec;
    delete [] Velvec;
    delete [] ACCvec;
    delete [] Pressure;
    delete [] Prop1;
    delete [] Prop2;
    delete [] Prop3;
    delete [] Density;
    delete [] Mass;
    delete [] sh;
    delete [] Tag;
    delete [] Sigma;
    delete [] Strain;

   //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    if (DEMParticles.Size()>0)
    {
    if (N_Faces>0)
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
    oss << "        " << fn.CStr() <<":/DFTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << "   <Grid Name=\"DEM_Center\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << DEMParticles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << DEMParticles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/DEMPosition \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DEMParticles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/DEMRadius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DEMParticles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/DEMTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DEMParticles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/DEMVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVel\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DEMParticles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/DEMAngVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Force\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DEMParticles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/DEMForce\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Torque\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << DEMParticles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/DEMTorque\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << "   <Grid Name=\"SPHCenter\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"10\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Position\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Acceleration\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Acceleration \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Density \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Pressure \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"" << OutputName[0] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/" << OutputName[0] << " \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"" << OutputName[1] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/" << OutputName[1] << " \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"" << OutputName[2] << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"10\"  Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/" << OutputName[2] << " \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Sigma\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Sigma \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Strain\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Strain \n";
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

//////////////////////////// FORCE CALCULATION METHODS ///////////////////////////

inline void Domain::CalcForce11(Particle * P1, Particle * P2)
{
    // This function computes the contact force between fluid-fluid paricles 
	double h		= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	Periodic_X_Correction(xij, h, P1, P2);
	double rij	= norm(xij);

	if ((rij/h)<=Cellfac)
	{
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;
		Vec3_t vij	= P1->v - P2->v;


		if (!P1->IsFree)
		{
			di = DensitySolid(P2->PresEq, P2->Cs, P2->P0,P1->Pressure, P2->RefDensity);
			mi = P1->FPMassC * P2->Mass;
		}
		else
		{
			di = P1->Density;
			mi = P1->Mass;
		}

		if (!P2->IsFree)
		{
			dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->Pressure, P1->RefDensity);
			mj = P2->FPMassC * P1->Mass;
		}
		else
		{
			dj = P2->Density;
			mj = P2->Mass;
		}

		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h);

		// Artificial Viscosity
		double PIij = 0.0;
		if (Alpha!=0.0 || Beta!=0.0){
			double Ci,Cj;
			if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);										///<(2.75) Li, Liu Book
			if (dot(vij,xij)<0) PIij = (-Alpha*0.5*(Ci+Cj)*MUij+Beta*MUij*MUij)/(0.5*(di+dj));		///<(2.74) Li, Liu Book
		}

		// Tensile Instability
		double TIij = 0.0;
		if (P1->TI > 0.0 || P2->TI > 0.0){
			double Ri,Rj;
			Ri = 0.0;
			Rj = 0.0;
			if ((P1->Pressure+P2->Pressure) < 0.0)
			{
				if (P1->Pressure < 0.0) Ri = -P1->Pressure/(di*di);
				if (P2->Pressure < 0.0) Rj = -P2->Pressure/(dj*dj);
			}
			TIij = (P1->TI*Ri + P2->TI*Rj)*pow((K/Kernel(Dimension, KernelType, (P1->TIInitDist + P2->TIInitDist)/(2.0*h), h)),(P1->TIn+P2->TIn)/2.0);
		}
        // Delta-SPH filter for the density and pressure - see Antuono et al (2012).
		double delfil = 0.0;
		if (DeltaSPH>1.0e-12){
			double Ci,Cj;
			if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
               delfil = DeltaSPH*2.0*h*(0.5*(Ci+Cj))*(dj-di)*(dot( xij, GK*xij )/(rij*rij + 0.01*h*h)); // delta-SPH from Molteni and Colagrossi (2009).
               //delfil = DeltaSPH*(0.5*(Ci+Cj))*(dj-di)*(dot( xij, GK*xij )/(rij + 0.1*h)); // delta-SPH from Ferrari et al. (2009).
		}

		// Real Viscosity
		Mat3_t StrainRate;
		set_to_zero(StrainRate);
		Vec3_t VI = 0.0;
		Vec3_t vab=0.0;
		if ((P1->NoSlip || P2->NoSlip) || (P1->IsFree&&P2->IsFree))
		{
			double Mu = 0.0;
			if (P1->IsFree&&P2->IsFree)
			{
				vab	= vij;
				Mu	= 2.0*P1->Mu*P2->Mu/(P1->Mu+P2->Mu);
			}
			else
			{
				// No-Slip velocity correction
				if (P1->IsFree)	vab = P1->v - (2.0*P2->v-P2->NSv); else vab = (2.0*P1->v-P1->NSv) - P2->v;
				if (!P1->IsFree) Mu = P2->Mu;
				if (!P2->IsFree) Mu = P1->Mu;
			}
			if (P1->T0>0.0 || P2->T0>0.0 || (P1->LES&&P2->LES))
			{
				StrainRate =	2.0*vab(0)*xij(0) , vab(0)*xij(1)+vab(1)*xij(0) , vab(0)*xij(2)+vab(2)*xij(0) ,
								vab(0)*xij(1)+vab(1)*xij(0) , 2.0*vab(1)*xij(1) , vab(1)*xij(2)+vab(2)*xij(1) ,
						 		vab(0)*xij(2)+vab(2)*xij(0) , vab(1)*xij(2)+vab(2)*xij(1) , 2.0*vab(2)*xij(2) ;
				StrainRate = 	-GK * StrainRate;
			}

			Viscous_Force(VisEq, VI, Mu, di, dj, GK, vab, Dimension, KernelType, rij, h, xij, vij);
		}

		// XSPH Monaghan
		if (XSPH != 0.0 && (P1->IsFree&&P2->IsFree))
		{
			omp_set_lock(&P1->my_lock);
			P1->VXSPH		+= XSPH*mj/(0.5*(di+dj))*K*-vij;
			omp_unset_lock(&P1->my_lock);

			omp_set_lock(&P2->my_lock);
			P2->VXSPH		+= XSPH*mi/(0.5*(di+dj))*K*vij;
			omp_unset_lock(&P2->my_lock);
		}

		// Calculating the forces for the particle 1 & 2
		Vec3_t temp	= 0.0;
		double temp1	= 0.0;

		if (GradientType == 0)
			temp		= -1.0*( P1->Pressure/(di*di) + P2->Pressure/(dj*dj) + PIij + TIij ) * GK*xij + VI;
		else
			temp		= -1.0*( (P1->Pressure + P2->Pressure)/(di*dj)       + PIij + TIij ) * GK*xij + VI;

		if (Dimension == 2) temp(2) = 0.0;
		temp1		= dot( vij , GK*xij );

		omp_set_lock(&P1->my_lock);
		P1->a					+= mj * temp;
		P1->dDensity	+= mj * (di/dj) * temp1 + mj*(-delfil/dj);

		if (P1->IsFree)
		{
			if (P1->T0>0.0 || P1->LES)	P1->StrainRate		= P1->StrainRate + mj/dj*StrainRate;
			if (SWIType == 1)						P1->S							= P1->S + mj/dj*vab(0)*xij(1)*-GK;
			P1->ZWab	+= mj/dj* K;
		}
		else
			P1->ZWab	= 1.0;

		if (P1->Shepard)
		if (P1->ShepardCounter == P1->ShepardStep)
		P1->SumDen += mj*K;
		omp_unset_lock(&P1->my_lock);


		omp_set_lock(&P2->my_lock);
		P2->a					-= mi * temp;
		P2->dDensity	+= mi * (dj/di) * temp1 + mi*(delfil/di);

		if (P2->IsFree)
		{
			if (P2->T0>0.0 || P2->LES)	P2->StrainRate		= P2->StrainRate + mi/di*StrainRate;
			if (SWIType ==1)						P2->S		 					= P2->S + mi/di*vab(0)*xij(1)*-GK;
			P2->ZWab	+= mi/di* K;
		}
		else
			P2->ZWab	= 1.0;

		if (P2->Shepard)
		if (P2->ShepardCounter == P2->ShepardStep)
		P2->SumDen += mi*K;
		omp_unset_lock(&P2->my_lock);
    }
}

inline void Domain::CalcForce2233(Particle * P1, Particle * P2)
{
    // This function computes the contact force between soil-soil particles or solid-solid paricles 
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

	if ((rij/h)<=Cellfac)
	{
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;

		if (P1->Material*P2->Material == 9)
		{
			if (!P1->IsFree) {di = P2->Density;	mi = P1->FPMassC * P2->Mass;} else {di = P1->Density;	mi = P1->Mass;}
			if (!P2->IsFree) {dj = P1->Density;	mj = P2->FPMassC * P1->Mass;} else {dj = P2->Density;	mj = P2->Mass;}
		}
		else
		{
			if (!P1->IsFree)
			{
				di = DensitySolid(P2->PresEq, P2->Cs, P2->P0,P1->Pressure, P2->RefDensity);
				mi = P1->FPMassC * P2->Mass;
			}
			else
			{
				di = P1->Density;
				mi = P1->Mass;
			}
			if (!P2->IsFree)
			{
				dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->Pressure, P1->RefDensity);
				mj = P2->FPMassC * P1->Mass;
			}
			else
			{
				dj = P2->Density;
				mj = P2->Mass;
			}
		}

		Vec3_t vij	= P1->v - P2->v;
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h);


		Mat3_t Sigmaj,Sigmai;
		set_to_zero(Sigmaj);
		set_to_zero(Sigmai);
		Sigmai = P1->Sigma;
		Sigmaj = P2->Sigma;

		Mat3_t Strainj,Straini;
		set_to_zero(Strainj);
		set_to_zero(Straini);
		Straini = P1->Strain;
		Strainj = P2->Strain;

//		if (P1->IsFree) Sigmai = P1->Sigma; else  Sigmai = P2->Sigma;
//		if (P2->IsFree) Sigmaj = P2->Sigma; else  Sigmaj = P1->Sigma;

		// Artificial Viscosity
		Mat3_t PIij;
		set_to_zero(PIij);
		if (Alpha!=0.0 || Beta!=0.0)
		{
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);					///<(2.75) Li, Liu Book
			double Cij;
			if (P1->Material*P2->Material == 9)
				Cij = 0.5*(P1->Cs+P2->Cs);
			else
			{
				double Ci,Cj;
				if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
				if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
				Cij = 0.5*(Ci+Cj);
			}
			if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;		///<(2.74) Li, Liu Book
			//PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;		///<(2.74) Li, Liu Book
		}

		// Tensile Instability
		Mat3_t TIij;
		set_to_zero(TIij);
		if (P1->TI > 0.0 || P2->TI > 0.0) TIij = pow((K/Kernel(Dimension, KernelType, (P1->TIInitDist + P2->TIInitDist)/(2.0*h), h)),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);


        // Delta-SPH filter for the density and pressure - see Antuono et al (2012).
		double delfil = 0.0;
		if (DeltaSPH>1.0e-12)
		{
			double Ci,Cj;
			if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
               delfil = DeltaSPH*2.0*h*(0.5*(Ci+Cj))*(dj-di)*(dot( xij, GK*xij )/(rij*rij + 0.01*h*h)); // delta-SPH from Molteni and Colagrossi (2009).
               //delfil = DeltaSPH*(0.5*(Ci+Cj))*(dj-di)*(dot( xij, GK*xij )/(rij + 0.1*h)); // delta-SPH from Ferrari et al. (2009).
               //delfil = DeltaSPH(P1->DSPH, DeltaSPH, 0.5*(Ci+Cj), di, dj, xij, rij, h, GK); 
		}

        // Similar Delta-SPH filter for stress tensor solely (NO stress gradient)
		Mat3_t delstress;
		set_to_zero(delstress);
		if (DeltaStress>1.0e-12)
		{
			double Ci,Cj;
			if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
               delstress = DeltaStress*h*(0.5*(Ci+Cj))*(dot( xij, GK*xij ))/((0.5*(di+dj))*(rij*rij + 0.01*h*h))*(Sigmaj - Sigmai); // 1.0 delta-SPH from Molteni and Colagrossi (2009).
               //delstress = DeltaStress*2.*h*(0.5*(Ci+Cj))*(dot( xij, GK*xij ))/((1*(di*dj))*(rij*rij + 0.01*h*h))*(Strainj - Straini); // 1.0 delta-SPH from Molteni and Colagrossi (2009).
               //delstress = DeltaSPH*(0.5*(Ci+Cj))*(dot( xij, GK*xij ))/((0.5*(di+dj))*(rij + 0.1*h))*(Sigmaj - Sigmai); // delta-SPH from Ferrari et al. (2009).
               //delfil = DeltasSPH(P1->DSPH, P1->DeltaSPH, 0.5*(Ci+Cj), di, dj, xij, rij, h, GK); 
		}

		// Shitf correction term
		Vec3_t TransC;
		set_to_zero(TransC);
		if (Shifting){
			double Cij = 0.5*(P1->Cs+P2->Cs);
			//TransC = -pow((2.*h),2)*( 1. + 0.2*pow( (K/Kernel(Dimension, KernelType, InitialDist/h, h)), 4) )*(1./(di+dj))*(GK*xij);
			TransC = -0.5*2.*h*Cij*( 1. + 0.2*pow( (K/Kernel(Dimension, KernelType, InitialDist/h, h)), 4.) )*(GK*xij);	
		} 


		// NoSlip BC velocity correction
		Vec3_t vab = 0.0;
		if (P1->IsFree&&P2->IsFree)
		{
			vab = vij;
		}
		else
		{
			if (P1->NoSlip || P2->NoSlip)
			{
				// No-Slip velocity correction
				if (P1->IsFree)	vab = P1->v - (2.0*P2->v-P2->NSv); else vab = (2.0*P1->v-P1->NSv) - P2->v;
			}
			// Please check
			if (!(P1->NoSlip || P2->NoSlip))
			{
				if (P1->IsFree) vab = P1->v - P2->vb; else vab = P1->vb - P2->v;
//				if (P1->IsFree) vab(0) = P1->v(0) + P2->vb(0); else vab(0) = -P1->vb(0) - P2->v(0);
			}
		}

		Mat3_t StrainRate,RotationRate;
		set_to_zero(StrainRate);
		set_to_zero(RotationRate);

		// Calculation strain rate tensor
		StrainRate(0,0) = 2.0*vab(0)*xij(0);
		StrainRate(0,1) = vab(0)*xij(1)+vab(1)*xij(0);
		StrainRate(0,2) = vab(0)*xij(2)+vab(2)*xij(0);
		StrainRate(1,0) = StrainRate(0,1);
		StrainRate(1,1) = 2.0*vab(1)*xij(1);
		StrainRate(1,2) = vab(1)*xij(2)+vab(2)*xij(1);
		StrainRate(2,0) = StrainRate(0,2);
		StrainRate(2,1) = StrainRate(1,2);
		StrainRate(2,2) = 2.0*vab(2)*xij(2);
		StrainRate	    = -0.5 * GK * StrainRate;

		// Calculation rotation rate tensor
		RotationRate(0,1) = vab(0)*xij(1)-vab(1)*xij(0);
		RotationRate(0,2) = vab(0)*xij(2)-vab(2)*xij(0);
		RotationRate(1,2) = vab(1)*xij(2)-vab(2)*xij(1);
		RotationRate(1,0) = -RotationRate(0,1);
		RotationRate(2,0) = -RotationRate(0,2);
		RotationRate(2,1) = -RotationRate(1,2);
		RotationRate	  = -0.5 * GK * RotationRate;

		// XSPH Monaghan
		if (XSPH != 0.0  && (P1->IsFree&&P2->IsFree)){
			omp_set_lock(&P1->my_lock);
			P1->VXSPH += XSPH*mj/(0.5*(di+dj))*K*-vij;
			omp_unset_lock(&P1->my_lock);

			omp_set_lock(&P2->my_lock);
			P2->VXSPH += XSPH*mi/(0.5*(di+dj))*K*vij;
			omp_unset_lock(&P2->my_lock);
		}



		// Calculating the stress gradient bewtween particles 1 and 2
		Vec3_t tempP1 = 0.0;
		Vec3_t tempP2 = 0.0;
		if(!KGcorrection){
			switch(GradientType){
				case 0:{
					Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , tempP1);
					tempP2=tempP1;
					break;}
				case 1:{
					Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , tempP1);
					tempP2=tempP1;
					break;}
			}
		}
		else{ // This is employed to correct the gradient kernel
			Vec3_t tcorrGKP1 = 0.0;
        	Vec3_t tcorrGKP2 = 0.0;
            Vec_t   RSVcorr,SOLcorrP1,SOLcorrP2;
            RSVcorr.Resize(Dimension);
            SOLcorrP1.Resize(Dimension);
            SOLcorrP2.Resize(Dimension);
            RSVcorr.SetValues(0.0);
            SOLcorrP1.SetValues(0.0);
            SOLcorrP2.SetValues(0.0);
            switch(Dimension){
              case 2:{
 				RSVcorr(0) = GK*xij(0);
                RSVcorr(1) = GK*xij(1);
                //if (P1->ID==-100 && P2->ID==-100) std::cout << RSVcorr << " =1= \n" << P1->Lcorr << " " << std::endl;
                Sol(-1.0*P1->Lcorr,RSVcorr,SOLcorrP1);
                tcorrGKP1(0) = SOLcorrP1(0);
                tcorrGKP1(1) = SOLcorrP1(1);
                tcorrGKP1(2) = 0.0;
                //if (P1->ID==-100 && P2->ID==-100) std::cout << RSVcorr << " =2= \n" << P2->Lcorr << " " << std::endl;
                Sol(-1.0*P2->Lcorr,RSVcorr,SOLcorrP2);
                tcorrGKP2(0) = SOLcorrP2(0);
                tcorrGKP2(1) = SOLcorrP2(1);
                tcorrGKP2(2) = 0.0;
                break;}
              case 3:{
                RSVcorr(0) = GK*xij(0);
                RSVcorr(1) = GK*xij(1);
                RSVcorr(2) = GK*xij(2);
                Sol(-1.0*P1->Lcorr,RSVcorr,SOLcorrP1);
                tcorrGKP1(0) = SOLcorrP1(0);
                tcorrGKP1(1) = SOLcorrP1(1);
                tcorrGKP1(2) = SOLcorrP1(2);
                Sol(-1.0*P2->Lcorr,RSVcorr,SOLcorrP2);
                tcorrGKP2(0) = SOLcorrP2(0);
                tcorrGKP2(1) = SOLcorrP2(1);
                tcorrGKP2(2) = SOLcorrP2(2);
                break;}
            }

			// Calculating the forces for the particle 1 & 2
			switch(GradientType){
				case 0:{
					Mult( tcorrGKP1 , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , tempP1);
					Mult( tcorrGKP2 , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , tempP2);
					break;}
				case 1:{
					Mult( tcorrGKP1 , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , tempP1);
					Mult( tcorrGKP2 , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , tempP2);
					break;}
			}
            /*if (P1->ID== -100){
            std::cout << "VectorO: " << GK*xij << std::endl;
            std::cout << "Vector1: " << tcorrGKP1 << std::endl;
            std::cout << "Vector2: " << tcorrGKP2 << std::endl;}*/
            
        }


		if (Dimension == 2) tempP1(2) = tempP2(2) = 0.0;

		double temp1 = 0.0;
		temp1 = dot( vij , GK*xij );

		// Locking the particle 1 for updating the properties
		omp_set_lock(&P1->my_lock);
		P1->a			+= mj * tempP1;
		P1->dDensity	+= mj * (di/dj) * temp1 + mj*(-delfil/dj);

		if (P1->IsFree)
		{
			P1->Shift 		+= TransC*mj/dj;
			P1->FilStress 	= P1->FilStress - mj/dj*delstress;
			P1->ZWab		+= mj/dj* K;
			P1->StrainRate	 = P1->StrainRate + mj/dj*StrainRate;
			P1->RotationRate = P1->RotationRate + mj/dj*RotationRate;
			if (SWIType ==1) P1->S = P1->S + mj/dj*vab(0)*xij(1)*-GK;
		}
		else
			P1->ZWab	= 1.0;

		if (P1->Shepard)
			if (P1->ShepardCounter == P1->ShepardStep)
				P1->SumDen += mj*K;
		omp_unset_lock(&P1->my_lock);

		// Locking the particle 2 for updating the properties
		omp_set_lock(&P2->my_lock);
		P2->a			-= mi * tempP2;
		P2->dDensity	+= mi * (dj/di) * temp1 + mi*(delfil/di);
		if (P2->IsFree)
		{
			P2->Shift 		+= -TransC*mi/di;
			P2->FilStress 	= P2->FilStress + mi/di*delstress;
			P2->ZWab		+= mi/di* K;
			P2->StrainRate	 = P2->StrainRate + mi/di*StrainRate;
			P2->RotationRate = P2->RotationRate + mi/di*RotationRate;
			if (SWIType ==1) P2->S = P2->S + mi/di*vab(0)*xij(1)*-GK;
		}
		else
			P2->ZWab	= 1.0;

		if (P2->Shepard)
			if (P2->ShepardCounter == P2->ShepardStep)
				P2->SumDen += mi*K;

		omp_unset_lock(&P2->my_lock);
	}
}

inline void Domain::CalcForce12(Particle * P1, Particle * P2)
{
    // This function computes the contact force between fluid-solid paricles 
	double h	= (P1->h+P2->h)/2;
//	double h	= std::max(P1->h,P2->h);
	Vec3_t xij	= P1->x - P2->x;

	Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

	if ((rij/h)<=Cellfac)
	{
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		Vec3_t vij	= P1->v - P2->v;
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);

		if (P1->Material == 1)
		{
			di = P1->Density;
			mi = P1->Mass;
			dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->FSIPressure, P1->RefDensity);
			mj = P1->Mass;
		}
		else
		{
			di = DensitySolid(P2->PresEq, P2->Cs, P2->P0,P1->FSIPressure, P2->RefDensity);
			mi = P2->Mass;
			dj = P2->Density;
			mj = P2->Mass;
		}

		// Artificial Viscosity
		double PIij	= 0.0;
		double Alpha	= 0.0;
		double Beta	= 0.0;
		if (Alpha!=0.0 || Beta!=0.0)
		{
			double Ci,Cj;
			if (P1->Material == 1)
			{
				Alpha	= P1->Alpha;
				Beta	= P1->Beta;
				Ci	= SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
				Cj	= SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity);
			}
			else
			{
				Alpha	= P2->Alpha;
				Beta	= P2->Beta;
				Ci 	= SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity);
				Cj 	= SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
			}
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);						///<(2.75) Li, Liu Book
			if (dot(vij,xij)<0) PIij = (-Alpha*0.5*(Ci+Cj)*MUij+Beta*MUij*MUij)/(0.5*(di+dj));		///<(2.74) Li, Liu Book
		}

		// Real Viscosity
		Mat3_t StrainRate;
		set_to_zero(StrainRate);
		Vec3_t VI = 0.0;
		Vec3_t vab= 0.0;
		double Mu = 0.0;
		if (P1->Material == 1)
		{
			vab = P1->v - (2.0*P2->v-P2->FSINSv);
			Mu = P1->Mu;
		}
		else
		{
			vab = (2.0*P1->v-P1->FSINSv) - P2->v;
			Mu = P2->Mu;
		}
		if ((P1->T0>0.0 || P2->T0>0.0) || (P1->LES || P2->LES))
		{
			StrainRate =	2.0*vab(0)*xij(0)            , vab(0)*xij(1)+vab(1)*xij(0) , vab(0)*xij(2)+vab(2)*xij(0) ,
					vab(0)*xij(1)+vab(1)*xij(0)  , 2.0*vab(1)*xij(1)           , vab(1)*xij(2)+vab(2)*xij(1) ,
					 vab(0)*xij(2)+vab(2)*xij(0) , vab(1)*xij(2)+vab(2)*xij(1) , 2.0*vab(2)*xij(2)           ;
			StrainRate = -GK * StrainRate;
		}

		Viscous_Force(VisEq, VI, Mu, di, dj, GK, vab, Dimension, KernelType, rij, h, xij, vij);


		// Calculating the forces for the particle 1 & 2
		Vec3_t temp	= 0.0;
		double temp1	= 0.0;
		if (P1->Material == 1)
		{
			if (GradientType == 0)
				temp		= -1.0*( P1->Pressure/(di*di) + P2->FSIPressure/(dj*dj) + PIij ) * GK*xij + VI;
			else
				temp		= -1.0*( (P1->Pressure + P2->FSIPressure)/(di*dj)       + PIij ) * GK*xij + VI;

			if (Dimension == 2) temp(2) = 0.0;
			temp1		= dot( vij , GK*xij );

			omp_set_lock(&P1->my_lock);
				P1->a		+= mj * temp;
				P1->dDensity	+= mj * (di/dj) * temp1;
				if (P1->T0>0.0 || P1->LES)	P1->StrainRate	 = P1->StrainRate + mj/dj*StrainRate;
			omp_unset_lock(&P1->my_lock);

			omp_set_lock(&P2->my_lock);
				P2->a		-= mi * temp;
			omp_unset_lock(&P2->my_lock);
		}
		else
		{
			if (GradientType == 0)
				temp		= -1.0*( P1->FSIPressure/(di*di) + P2->Pressure/(dj*dj) + PIij ) * GK*xij + VI;
			else
				temp		= -1.0*( (P1->FSIPressure + P2->Pressure)/(di*dj)       + PIij ) * GK*xij + VI;
			if (Dimension == 2) temp(2) = 0.0;
			temp1		= dot( vij , GK*xij );

			omp_set_lock(&P1->my_lock);
				P1->a		+= mj * temp;
			omp_unset_lock(&P1->my_lock);


			omp_set_lock(&P2->my_lock);
				P2->a		-= mi * temp;
				P2->dDensity	+= mi * (dj/di) * temp1;
				if (P2->T0>0.0 || P2->LES)	P2->StrainRate	 = P2->StrainRate + mi/di*StrainRate;
			omp_unset_lock(&P2->my_lock);
		}

    }
}

inline void Domain::CalcForce13(Particle * P1, Particle * P2)
{
    // This function computes the contact force between fluid-soil paricles 
	double h	= std::min(P1->h,P2->h);
	Vec3_t xij	= P1->x - P2->x;

	Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

	if ((rij/h)<=Cellfac)
	{
		double K = Kernel(Dimension, KernelType, rij/h, h)/(P1->Density*P2->Density);
		double SF1=0.0,SF2=0.0;
		Vec3_t SFt=0.0,v=0.0;
		switch(SWIType)
		{
			case 0:
				if (P1->Material == 3 )
				{
					v = P2->v-P1->v;
					Seepage(P1->SeepageType, P1->k, P1->k2, P2->MuRef, P2->RefDensity, SF1, SF2);
					SFt = (SF1*v + SF2*norm(v)*v) *K;
					if (Dimension == 2) SFt(2) = 0.0; 

					omp_set_lock(&P1->my_lock);
						P1->a += P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a -= P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				else
				{
					v = P1->v-P2->v;
					Seepage(P2->SeepageType, P2->k, P2->k2, P1->MuRef, P1->RefDensity, SF1, SF2);
					SFt = (SF1*v + SF2*norm(v)*v) *K;
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a -= P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a += P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				break;
			case 1:
				if (P1->Material == 3 )
				{
					v = P2->v-P1->v;
					if (P1->ZWab<0.25)
					{
						double Cd = 24.0*(P2->MuRef/P2->RefDensity)/(P1->d*norm(v)+0.01*h*h) + 2.0;
						SFt = (3.0/(4.0*P1->d)*P2->RefDensity*(1.0-P1->n0)*Cd*norm(v)*v) *K;
						SFt(1) += (P2->RefDensity*(1.0-P1->n0)*norm(v)*fabs(P2->S-P1->S)) *K;
					}
					else
					{
						Seepage(P1->SeepageType, P1->k, P1->k2, P2->MuRef, P2->RefDensity, SF1, SF2);
						SFt = (SF1*v + SF2*norm(v)*v) *K;
					}
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a += P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a -= P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				else
				{
					v = P1->v-P2->v;
					if (P2->ZWab<0.25)
					{
						double Cd = 24.0*(P1->MuRef/P1->RefDensity)/(P2->d*norm(v)+0.01*h*h) + 2.0;
						SFt = (3.0/(4.0*P2->d)*P1->RefDensity*(1.0-P2->n0)*Cd*norm(v)*v) *K;
						SFt(1) += (P1->RefDensity*(1.0-P2->n0)*norm(v)*fabs(P1->S-P2->S)) *K;;
					}
					else
					{
						Seepage(P2->SeepageType, P2->k, P2->k2, P1->MuRef, P1->RefDensity, SF1, SF2);
						SFt = (SF1*v + SF2*norm(v)*v) *K;
					}
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
						P1->a -= P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
						P2->a += P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				break;
			case 2:
				if (P1->Material == 3 )
				{
					double GK	= GradKernel(Dimension, KernelType, rij/h, h)/(P1->Density*P2->Density);
					v = P2->v-P1->v;
					Seepage(P1->SeepageType, P1->k, P1->k2, P2->MuRef, P2->RefDensity, SF1, SF2);
					SFt = (SF1*v + SF2*norm(v)*v) *K;
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
					P1->a += P2->Mass*SFt - P2->Mass*P2->Pressure*GK*xij;
                    P1->qmin = std::min(P1->qmin, rij/P1->SatR );
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
				    P2->a -= P1->Mass*SFt;
					omp_unset_lock(&P2->my_lock);
				}
				else
				{
					double GK	= GradKernel(Dimension, KernelType, rij/h, h)/(P1->Density*P2->Density);
					v = P1->v-P2->v;
					Seepage(P2->SeepageType, P2->k, P2->k2, P1->MuRef, P1->RefDensity, SF1, SF2);
					SFt = (SF1*v + SF2*norm(v)*v) *K;
					if (Dimension == 2) SFt(2) = 0.0;

					omp_set_lock(&P1->my_lock);
					P1->a -= P2->Mass*SFt;
					omp_unset_lock(&P1->my_lock);

					omp_set_lock(&P2->my_lock);
	        		P2->a += P1->Mass*SFt + P1->Mass*P1->Pressure*GK*xij;
                    P2->qmin = std::min(P2->qmin, rij/P2->SatR );
					omp_unset_lock(&P2->my_lock);
				}
				break;
			case 3:
				break;
			default:
				std::cout << "Soil-Water Interaction Type No is out of range. Please correct it and run again" << std::endl;
				std::cout << "0 => The seepage force + The bouyant unit weight of soil" << std::endl;
				std::cout << "1 => The seepage force + The surface erosion(Lift+Drag)) + The bouyant unit weight of soil" << std::endl;
				std::cout << "2 => The seepage force + The pore water pressure from water particles " << std::endl;
				std::cout << "3 => Zero interaction force" << std::endl;
				abort();
				break;
		}
	}
}

}; // namespace SPH

#endif // SPH_DOMAIN_H
