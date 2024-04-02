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

#ifndef SPH_PARTICLE_H
#define SPH_PARTICLE_H

#include <mechsys/sph/Functions.h>

namespace SPH {

	class Particle
	{
	public:
		// Shepard density correction
		bool   	Shepard;	///< Shepard Filter for the density
		size_t	ShepardCounter;	///< Count number of contributing particles
		size_t	ShepardStep;	///< Cycle number for shepard correction
		double	ZWab;		///< Summation of mb/db*Wab for neighbour particles of the particle a (for Shepard filter)
		double	SumDen;		///< Summation of mb*Wab for neighbour particles of the particle a (for Shepard filter)
        double  SumMin;		///< Minimum value of summation Zwab for Shepard filter
        Mat_t   Lcorr;		///< Kernel Gradient Correction matrix
        Mat3_t  FilStress;	///< Delta-SPH type filter for the soil stress tensor
        Vec3_t  Shift;      ///< Vector to modify distribution of the SPH particles

		bool   	IsFree;		///< Check the particle if it is free to move or not
		size_t	InOut;		///< Check the particle if it is in-flow or out-flow or not
		bool   	IsSat;		///< Check the particle if it is Saturated or not
		bool   	SatCheck;	///< Check the particle Saturation at each time step
		bool   	NoSlip;		///< No-Slip BC

		int    	ID;		    ///< an Integer value to identify the particle set
		int     Material;	///< an Integer value to identify the particle material type: 1 = Fluid, 2 = Solid, 3 = Soil
        int     Constitutive;   ///< an Integer value to identify the particle constitutive model type: 11-15= Fluid, 2 = Solid, 31-32 = Soil

		Vec3_t	x;		///< Position of the particle n
		Vec3_t	vb;		///< Velocity of the particle n-1 (Modified Verlet)
		Vec3_t	va;		///< Velocity of the particle n+1/2 (Leapfrog)
		Vec3_t	v;		///< Velocity of the particle n+1
		Vec3_t	NSv;		///< Velocity of the fixed particle for no-slip BC
		Vec3_t	FSINSv;		///< Velocity of the fixed particle for no-slip BC in FSI
		Vec3_t	VXSPH;		///< Mean Velocity of neighbor particles for updating the particle position (XSPH)
		Vec3_t	a;		///< Acceleration of the particle n

		size_t	PresEq;		///< Selecting variable to choose an equation of state
		double	Cs;		///< Speed of sound
		double	P0;		///< background pressure for equation of state
		double 	Pressure;	///< Pressure of the particle n+1
		double 	FSIPressure;	///< Pressure of the particle n+1 in FSI

		double	Density;	///< Density of the particle n+1
		double 	Densitya;	///< Density of the particle n+1/2 (Leapfrog)
		double 	Densityb;	///< Density of the particle n-1 (Modified Verlet)
		double 	dDensity;	///< Rate of density change in time based on state equations n
		double 	RefDensity;	///< Reference Density of Particle
		double 	FPMassC;	///< Mass coefficient for fixed particles to avoid leaving particles
		double 	Mass;		///< Mass of the particle

		Mat3_t	StrainRate;	///< Global shear Strain rate tensor n
		Mat3_t	RotationRate;	///< Global rotation tensor n
		double	ShearRate;	///< Global shear rate for fluids
		double	SBar;		///< shear component for LES

		Mat3_t	ShearStress;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1
		Mat3_t	ShearStressa;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1/2 (Leapfrog)
		Mat3_t	ShearStressb;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n-1 (Modified Verlet)

		Mat3_t	Sigma;		///< Cauchy stress tensor (Total Stress) n+1
		//    Mat3_t	FSISigma;	///< Cauchy stress tensor (Total Stress) n+1 in FSI
		Mat3_t	Sigmaa;		///< Cauchy stress tensor (Total Stress) n+1/2 (Leapfrog)
		Mat3_t	Sigmab;		///< Cauchy stress tensor (Total Stress) n-1 (Modified Verlet)

		Mat3_t	Strain;		///< Total Strain n+1
		Mat3_t	Straina;	///< Total Strain n+1/2 (Leapfrog)
		Mat3_t	Strainb;	///< Total Strain n-1 (Modified Verlet)

		Mat3_t	TIR;		///< Tensile Instability stress tensor R
		double	TI;		///< Tensile instability factor
		double	TIn;		///< Tensile instability power
		double 	TIInitDist;	///< Initial distance of particles for calculation of tensile instability

		double 	Alpha;		///< Dynamic viscosity coefficient of the fluid particle
		double 	Beta;		///< Dynamic viscosity coefficient of the fluid particle
		double 	Mu;			///< Dynamic viscosity coefficient of the fluid particle
		double 	MuRef;		///< Reference Dynamic viscosity coefficient
        double  Mumax;      ///< Highest viscosity for Bilinear Bingham model
		double 	T0;			///< Yield stress for Bingham fluids
		double 	m;			///< Regularisation parameter for Bingham fluids
        double  t0mu;       ///< Time parameter for Herschel-Bulkley-Kee Model
        double  gamma0;     ///< Constant for the Bilinear Bingham model
        double  Mu0;        ///< Constant parameter of viscosity for Cross model
        double  Kcross;     ///< Constant consistency index for Cross model
        
		size_t	VisM;		///< Non-Newtonian viscosity method
		bool	LES;		///< Large eddy simulation using sub-particle scale
		double	CSmag;		///< Coefficient of Smagorinsky-Lilly model

		double 	G;			///< Shear modulus
		double 	K;			///< Bulk modulus
		double	Sigmay;		///< Tensile yield stress
		size_t	Fail;		///< Failure criteria

		double	c;		///< Cohesion
		double	phi;	///< Friction angel
		double	psi;	///< Dilation angel
		double	n;		///< Porosity
		double	n0;		///< Initial Prosity
		double	k;		///< Permeability
		double	k2;		///< Second Permeability for the Forchheimer Eq
		double	Saturation;	///< Saturation of the soil
		double	d;		///< effective particle size
		double	V;		///< Volume of a particle
		double	RhoF;		    ///< Density of water or any other fluids
		bool	VarPorosity;	///< If yes, it will calculate porosity and permeability based on new calculated porosity
		size_t	SeepageType;	///< Selecting variable to choose a Seepage method
		double	S;		       ///< Velocity derivative for surface erosion

        double  c1hp;   ///< Hypoplastic dimensionless parameter c1
        double  c2hp;   ///< Hypoplastic dimensionless parameter c2
        double  c3hp;   ///< Hypoplastic dimensionless parameter c3
        double  c4hp;   ///< Hypoplastic dimensionless parameter c4

        double  SatR;    ///< Radius to calculate saturation
		double 	h;		///< Smoothing length of the particle
		double  qmin;		///< Min unity distance "q=rij/h" among multiple interacting particles
		int    	LL;		///< Linked-List variable to show the next particle in the list of a cell
		int    	CC[3];		///< Current cell No for the particle (linked-list)
		int	    ct;		///< Correction step for the Modified Verlet Algorithm
		double	SumKernel;	///< Summation of the kernel value for neighbour particles
		double	FSISumKernel;	///< Summation of the kernel value for neighbour particles in FSI
		bool	FirstStep;	///< to initialize the integration scheme

		omp_lock_t my_lock;		///< Open MP lock


		// Constructor
		Particle				(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0, bool Fixed=false);

		// Methods
		void Move				(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin,size_t Scheme, Mat3_t I);	///< Update the important quantities of a particle
		void Move_MVerlet		(Mat3_t I, double dt);										///< Update the important quantities of a particle
		void Move_Leapfrog      (Mat3_t I, double dt);										///< Update the important quantities of a particle
		void translate			(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin);
		void Mat1Newton         (double dt);
        void Mat1BP             (double dt);
        void Mat1BB             (double dt);
        void Mat1CSL            (double dt);
        void Mat1HBK            (double dt);
		void Mat2MVerlet		(double dt);
		void Mat3epMVerlet		(Mat3_t I, double dt);
        void Mat3hpMVerlet      (Mat3_t I, double dt);
		void Mat2Leapfrog		(double dt);
		void Mat3Leapfrog		(Mat3_t I, double dt);
        void Mat3hpMLeapfrog    (Mat3_t I, double dt);
		void ScalebackMat3      (size_t Dimension, size_t Scheme);
	};

    inline Particle::Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0,bool Fixed)
    {
    	ct = 0;
    	a = 0.0;
        x = x0;
        n = 0.0;
        n0 = 0.0;
        k = 0.0;
        k2 = 0.0;
    	
        Cs		= 0.0;
        P0		= 0.0;
        PresEq	= 0;

        Alpha	= 0.0;
        Beta	= 0.0;
    	
        va = 0.0;
        vb = 0.0;
        NSv = 0.0;
        FSINSv = 0.0;
        v = v0;
        VXSPH = 0.0;
        TI		= 0.0;
        TIn		= 4.0;
        TIInitDist  = 0.0;
    
        Densitya = 0.0;
        Densityb = 0.0;
        Density = Density0;
        RefDensity = Density0;
    
        Mass = Mass0;
        FPMassC = 1.0;
        IsFree = !Fixed;
        h = h0;
        SatR = h0;
        qmin = 1.0e16;
        Saturation = 0.0;
        Pressure=0.0;
        FSIPressure=0.0;
        ID = Tag;
        CC[0]= CC[1] = CC[2] = 0;
        LL=0;
        ZWab = 0.0;
        SumDen = 0.0;
        dDensity=0.0;
        ShearRate = 0.0;
        MuRef = Mu = 0.0;
    	VisM = 0;
        T0 = 0.0;
        m = 300.0;
        Mumax = 0.0;  // Value redefined in Domain -InitialChecks-
        gamma0 = 0.0; // Value redefined in Domain -InitialChecks-
        Mu0 = 0.0;    // Value redefined in Domain -InitialChecks-
        Kcross = 0.0; // Value redefined in Domain -InitialChecks-
        t0mu = 1.0e-4;
        SumKernel = 0.0;
        FSISumKernel = 0.0;
        SumMin = 0.6;
        G = 0.0;
        K = 0.0;
        Material = 0;
        Constitutive = 0;
        Fail = 0;
        c = 0.0;
        phi = 0.0;
        psi = 0.0;
        d =0.0;
        Sigmay = 0.0;
        NoSlip = false;
        Shepard = false;
        InOut = 0;
        FirstStep = true;
        V = Mass/RefDensity;
        RhoF = 0.0;
        IsSat = false;
        SatCheck = false;
        ShepardStep = 40;
        ShepardCounter = 0;
        S = 0.0;
        VarPorosity = false;
        SeepageType = 0;
        S = 0;
    	LES = false;
    	SBar = 0.0;
    	CSmag = 0.17;

        c1hp = 0.0;
        c2hp = 0.0;
        c3hp = 0.0;
        c4hp = 0.0;
    
        set_to_zero(Shift);
        set_to_zero(FilStress);

        set_to_zero(Strainb);
        set_to_zero(Strain);
        set_to_zero(Sigmab);
        set_to_zero(Sigma);
    //    set_to_zero(FSISigma);
        set_to_zero(Sigmaa);
        set_to_zero(ShearStress);
        set_to_zero(ShearStressb);
        set_to_zero(TIR);
        set_to_zero(StrainRate);
        set_to_zero(RotationRate);
        omp_init_lock(&my_lock);
    
    }
    
    inline void Particle::Move(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, size_t Scheme, Mat3_t I)
    { 
    	if (Scheme == 0)
    		Move_MVerlet(I, dt);
    	else
    		Move_Leapfrog(I, dt);
    
    
    	//Periodic BC particle position update
    	if (Domainsize(0)>0.0)
    	{
    		(x(0)>(domainmax(0))) ? x(0) -= Domainsize(0) : x(0);
    		(x(0)<(domainmin(0))) ? x(0) += Domainsize(0) : x(0);
    	}
    	if (Domainsize(1)>0.0)
    	{
    		(x(1)>(domainmax(1))) ? x(1) -= Domainsize(1) : x(1);
    		(x(1)<(domainmin(1))) ? x(1) += Domainsize(1) : x(1);
    	}
    	if (Domainsize(2)>0.0)
    	{
    		(x(2)>(domainmax(2))) ? x(2) -= Domainsize(2) : x(2);
    		(x(2)<(domainmin(2))) ? x(2) += Domainsize(2) : x(2);
    	}
    
    }
    
    
    inline void Particle::Move_MVerlet (Mat3_t I, double dt)
    {
    	if (FirstStep)
    	{
    		ct = 30;
    		FirstStep = false;
    	}
    
    	x += dt*(v+VXSPH) + 0.5*dt*dt*a;
    
    	if (ct == 30)
    	{
    		if (Shepard && ShepardCounter == ShepardStep)
    		{
    			if (ZWab>SumMin)
    			{
    				Densityb	= SumDen/ZWab;
    //				Densityb	= Density;
    				Density		= SumDen/ZWab;
    			}
    			else
    			{
    				Densityb	= Density;
    				Density		+=dt*dDensity;
    			}
    		}
    		else
    		{
    			Densityb		= Density;
    			Density			+=dt*dDensity;
    		}
    
    		vb	= v;
    		v	+=dt*a;
    	}
    	else
    	{
    		if (Shepard && ShepardCounter == ShepardStep)
    		{
    			if (ZWab>SumMin)
    			{
    				Densityb	= SumDen/ZWab;
    //				Densityb	= Density;
    				Density		= SumDen/ZWab;
    			}
    			else
    			{
    				double dens	= Density;
    				Density		= Densityb + 2.0*dt*dDensity;
    				Densityb	= dens;
    			}
    		}
    		else
    		{
    			double dens	= Density;
    			Density		= Densityb + 2.0*dt*dDensity;
    			Densityb	= dens;
    		}
    
    		Vec3_t temp;
    		temp	= v;
    		v		= vb + 2*dt*a;
    		vb		= temp;
    	}
    
    	switch (Constitutive){
        case 11:
            Mat1Newton(dt);
    	    break;
        case 12:
            Mat1BP(dt);
            break;
        case 13:
            Mat1BB(dt);
            break;
        case 14:
            Mat1CSL(dt);
            break;
        case 15:
            Mat1HBK(dt);
            break;
        case 2:
            Mat2MVerlet(dt);
            break;
        case 31:
            Mat3epMVerlet(I,dt);
            break;
        case 32:
            Mat3hpMVerlet(I,dt);
            break;
        default:
            std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
            std::cout << "11 => Fluid:Newtonian" << std::endl;
            std::cout << "12 => Fluid:Bingham-Papanastasiou" << std::endl;
            std::cout << "13 => Fluid:Bingham_Bilinear" << std::endl;
            std::cout << "14 => Fluid:Cross_Shao-Lo" << std::endl;
            std::cout << "15 => Fluid:Herschel-Bulkley-Papanastasiou_Kee" << std::endl;
    	    std::cout << "2 => Solid" << std::endl;
            std::cout << "31 => Soil:Elastoplastic" << std::endl;
            std::cout << "32 => Soil:Hypoplastic" << std::endl;
    	    abort();
    	    break;
        }
    	if (ct == 30) ct = 0; else ct++;
    	if (ShepardCounter == ShepardStep) ShepardCounter = 0; else ShepardCounter++;
    }


    inline void Particle::Move_Leapfrog(Mat3_t I, double dt)
    {
        if (FirstStep){
            Densitya = Density - 0.5*dt*dDensity;
            va = v - 0.5*dt*a;
        }
        Densityb = Densitya;
        Densitya += dt*dDensity;
        Density = 0.5*(Densitya+Densityb);
        vb = va;
        va += dt*a;
        v = (va + vb)/2.0;
        x += dt*va;
    
        switch (Constitutive){
        case 11:
            Mat1Newton(dt);
            break;
        case 12:
            Mat1BP(dt);
            break;
        case 13:
            Mat1BB(dt);
            break;
        case 14:
            Mat1CSL(dt);
            break;
        case 15:
            Mat1HBK(dt);
            break;
        case 2:
            Mat2Leapfrog(dt);
            break;
        case 31:
            Mat3Leapfrog(I,dt);
            break;
        case 32:
            Mat3hpMLeapfrog(I,dt);// LOOOOOOOOKKKKK correct
            break;
        default:
            std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
            std::cout << "11 => Fluid:Newtonian" << std::endl;
            std::cout << "12 => Fluid:Bingham-Papanastasiou" << std::endl;
            std::cout << "13 => Fluid:Bingham_Bilinear" << std::endl;
            std::cout << "14 => Fluid:Cross_Shao-Lo" << std::endl;
            std::cout << "15 => Fluid:Herschel-Bulkley-Papanastasiou_Kee" << std::endl;
            std::cout << "2  => Solid" << std::endl;
            std::cout << "31 => Soil:Elastoplastic" << std::endl;
            std::cout << "32 => Soil:Hypoplastic" << std::endl;
            abort();
            break;
        }
        if (FirstStep) FirstStep = false;
    }

    inline void Particle::Mat1Newton(double dt){
        Pressure    = EOS(PresEq, Cs, P0,Density, RefDensity);
        // LES model
        if (LES){
            double temp =  ((StrainRate(0,0)*StrainRate(0,0) + StrainRate(1,1)*StrainRate(1,1) + StrainRate(2,2)*StrainRate(2,2)) +
                           2.0*(StrainRate(0,1)*StrainRate(1,0) + StrainRate(0,2)*StrainRate(2,0) + StrainRate(1,2)*StrainRate(2,1)));

            SBar        = sqrt(2.0*temp);
            Mu  = MuRef + RefDensity*pow((CSmag*h),2.0)*SBar;
        }  
    }

    inline void Particle::Mat1BP(double dt){
        Pressure    = EOS(PresEq, Cs, P0,Density, RefDensity);
        double temp =  ((StrainRate(0,0)*StrainRate(0,0) + StrainRate(1,1)*StrainRate(1,1) + StrainRate(2,2)*StrainRate(2,2)) +
                       2.0*(StrainRate(0,1)*StrainRate(1,0) + StrainRate(0,2)*StrainRate(2,0) + StrainRate(1,2)*StrainRate(2,1)));

        ShearRate = sqrt(0.5*temp);

        // LES model
        if (LES){
            SBar = sqrt(2.0*temp);
            Mu   = MuRef + RefDensity*pow((CSmag*h),2.0)*SBar;
        }
    
        // Non-Newtonian viscosities computation
        // Regularised Bingham-Papanastasiou model
        if ((m*ShearRate) > 1.0e-12){
            Mu = MuRef + T0*(-expm1(-m*ShearRate))/ShearRate;}
        else{
            Mu = Mumax;}
    }

    inline void Particle::Mat1BB(double dt){
        Pressure    = EOS(PresEq, Cs, P0,Density, RefDensity);
        double temp =  ((StrainRate(0,0)*StrainRate(0,0) + StrainRate(1,1)*StrainRate(1,1) + StrainRate(2,2)*StrainRate(2,2)) +
                       2.0*(StrainRate(0,1)*StrainRate(1,0) + StrainRate(0,2)*StrainRate(2,0) + StrainRate(1,2)*StrainRate(2,1)));

        ShearRate = sqrt(0.5*temp);

        // LES model
        if (LES){
            SBar = sqrt(2.0*temp);
            Mu   = MuRef + RefDensity*pow((CSmag*h),2.0)*SBar;
        }

        // Non-Newtonian viscosities computation
        // Bingham Bilinear model
        if (ShearRate > 1.0e-12){
            Mu = MuRef + T0/ShearRate - ((T0*gamma0)/(ShearRate+gamma0))/ShearRate;}
        else{
            Mu = Mumax;}
    }

    inline void Particle::Mat1CSL(double dt){
        Pressure    = EOS(PresEq, Cs, P0,Density, RefDensity);
        double temp =  ((StrainRate(0,0)*StrainRate(0,0) + StrainRate(1,1)*StrainRate(1,1) + StrainRate(2,2)*StrainRate(2,2)) +
                       2.0*(StrainRate(0,1)*StrainRate(1,0) + StrainRate(0,2)*StrainRate(2,0) + StrainRate(1,2)*StrainRate(2,1)));

        ShearRate = sqrt(0.5*temp);

        // LES model
        if (LES){
            SBar = sqrt(2.0*temp);
            Mu   = MuRef + RefDensity*pow((CSmag*h),2.0)*SBar;
        }

        // Non-Newtonian viscosities computation
        // Cross model with Shao-Lo parameters
        Mu = (Mu0 + MuRef*Kcross*ShearRate)/(1.+Kcross*ShearRate);
    }

   inline void Particle::Mat1HBK(double dt){
        Pressure    = EOS(PresEq, Cs, P0,Density, RefDensity);
        double temp =  ((StrainRate(0,0)*StrainRate(0,0) + StrainRate(1,1)*StrainRate(1,1) + StrainRate(2,2)*StrainRate(2,2)) +
                       2.0*(StrainRate(0,1)*StrainRate(1,0) + StrainRate(0,2)*StrainRate(2,0) + StrainRate(1,2)*StrainRate(2,1)));

        ShearRate = sqrt(0.5*temp);

        // LES model
        if (LES){
            SBar = sqrt(2.0*temp);
            Mu   = MuRef + RefDensity*pow((CSmag*h),2.0)*SBar;
        }
    
        // Non-Newtonian viscosities computation
        // Regularised Herschel-Bulkley-Kee model (power-law)
        if ((t0mu*ShearRate) > 1.0e-12){
            Mu = MuRef*(1.+expm1(-t0mu*ShearRate)) + T0*(-expm1(-m*ShearRate))/ShearRate;}
        else{
            Mu = Mumax;}
    }

    inline void Particle::Mat2MVerlet(double dt){
    	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);
    
    	// Jaumann rate terms
    	Mat3_t RotationRateT, Stress,SRT,RS;
    	Trans(RotationRate,RotationRateT);
    	Mult(ShearStress,RotationRateT,SRT);
    	Mult(RotationRate,ShearStress,RS);
    
    	// Elastic prediction step (ShearStress_e n+1)
    	Stress			= ShearStress;
    	if (ct == 30)
    		ShearStress	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;
    	else
    		ShearStress	= 2.0*dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressb;
    	ShearStressb	= Stress;
    
    	if (Fail == 1)
    	{
    		double J2	= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                        + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    		//Scale back
    		ShearStress	= std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStress;
    	}
    
    	Sigma = -Pressure * OrthoSys::I + ShearStress;
    
    	Stress	= Strain;
    	if (ct == 30)
    		Strain	= dt*StrainRate + Strain;
    	else
    		Strain	= 2.0*dt*StrainRate + Strainb;
    	Strainb	= Stress;
    
    
    	if (Fail > 1)
    	{
    		std::cout<<"Undefined failure criteria for solids"<<std::endl;
    		abort();
    	}
    }
    
    inline void Particle::Mat2Leapfrog(double dt)
    {
    	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);
    
    	// Jaumann rate terms
    	Mat3_t RotationRateT,SRT,RS;
    	Trans(RotationRate,RotationRateT);
    	Mult(ShearStress,RotationRateT,SRT);
    	Mult(RotationRate,ShearStress,RS);
    
    	// Elastic prediction step (ShearStress_e n+1)
    	if (FirstStep)
    		ShearStressa	= -dt/2.0*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;
    
    	ShearStressb	= ShearStressa;
    	ShearStressa	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressa;
    
    	if (Fail == 1)
    	{
    		double J2	= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                        + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    		//Scale back
    		ShearStressa= std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStressa;
    	}
    	ShearStress	= 1.0/2.0*(ShearStressa+ShearStressb);
    
    	Sigma = -Pressure * OrthoSys::I + ShearStress;
    
    	if (FirstStep)
    		Straina	= -dt/2.0*StrainRate + Strain;
    	Strainb	= Straina;
    	Straina	= dt*StrainRate + Straina;
    	Strain	= 1.0/2.0*(Straina+Strainb);
    
    
    	if (Fail > 1)
    	{
    		std::cout<<"Undefined failure criteria for solids"<<std::endl;
    		abort();
    	}
    }
    

    inline void Particle::Mat3hpMVerlet(Mat3_t I, double dt)
    {// Hypoplatic contitutive model as shown in Peng et al 2015. A SPH approach...
    	Mat3_t RotationRateT, Stress, Straintmp, SRT,RS;
    	double I1,Ee,I1strain,trOe;

    	// Objective Zaremba-Jaumann rates terms
    	Trans(RotationRate,RotationRateT);
    	Mult(Sigma,RotationRateT,SRT);
    	Mult(RotationRate,Sigma,RS);

    	// Cohesion addition if greater than zero
    	Stress = Sigma;
    	if(c > 1.e-12){
           Stress(0,0) -= c;
           Stress(1,1) -= c;
           Stress(2,2) -= c;
    	}

    	// To avoid division by zero when the trace of the stress is used
    	I1 = Stress(0,0) + Stress(1,1) + Stress(2,2);
        if(fabs(I1) < 4.0){
            I1 = -4.0;
            Stress(0,0) = -1.0;
            Stress(1,1) = -1.0;
            Stress(2,2) = -2.0;
            Stress(0,1) = 0.0;
            Stress(0,2) = 0.0;
            Stress(1,2) = 0.0;
            Stress(1,0) = 0.0;
            Stress(2,0) = 0.0;
            Stress(2,1) = 0.0;
        }

	// Deviatoric stress tensor S
	//ShearStress = Stress - 1.0/3.0* I1 *OrthoSys::I;
        ShearStress = Stress;
	    ShearStress(0,0) = ShearStress(0,0) - I1/3.;
        ShearStress(1,1) = ShearStress(1,1) - I1/3.;
        ShearStress(2,2) = ShearStress(2,2) - I1/3.;

    	I1strain = StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2); // Trace of the strain tensor
        Ee	 = sqrt(indpAA(StrainRate));         // The Frobenius or Euclidean norm of the strain tensor, sqrt(D^2)
    	trOe = indpAB(Stress,StrainRate);        // Inner double product of stress and strain rate, tr(TD)

   
    	// Prediction of the stress tensor    
    	Stress	= Sigma;
    	if (ct == 30)
    		Sigma	= dt*( SRT + RS + (c1hp*I1)*StrainRate + (c2hp*I1strain + c3hp*trOe/I1)*Stress + (c4hp*Ee)*(Stress + ShearStress) + FilStress) + Sigma;
    	else
    		Sigma	= 2.0*dt*( SRT + RS + (c1hp*I1)*StrainRate + (c2hp*I1strain + c3hp*trOe/I1)*Stress + (c4hp*Ee)*(Stress + ShearStress) + FilStress) + Sigmab;
    	Sigmab	= Stress;


	double press = (Sigma(0,0) + Sigma(1,1) + Sigma(2,2))/3.0;
	if(press > c){
	    Sigma(0,0) -= (press - c);
	    Sigma(1,1) -= (press - c);
	    Sigma(2,2) -= (press - c);
	}

    	// Prediction of the strain tensor
    	Straintmp	= Strain;
    	if (ct == 30)
    		Strain	= dt*StrainRate + Strain;
    	else
    		Strain	= 2.0*dt*StrainRate + Strainb;
    	Strainb	= Straintmp;

    	// Notes:
    	// FilStress: filter to smooth spurious numerical oscillations in the stress rate tensor
    }

    inline void Particle::Mat3hpMLeapfrog(Mat3_t I, double dt)
    {// Hypoplatic contitutive model as shown in Peng et al 2015. A SPH approach...
        Mat3_t RotationRateT, Stress, Straintmp, SRT,RS;
        double I1,Ee,I1strain,trOe;

        // Objective Zaremba-Jaumann rates terms
        Trans(RotationRate,RotationRateT);
        Mult(Sigma,RotationRateT,SRT);
        Mult(RotationRate,Sigma,RS);

        // Cohesion addition if greater than zero
        Stress = Sigma;
        if(c > 1.e-12){
           Stress(0,0) -= c;
           Stress(1,1) -= c;
           Stress(2,2) -= c;
        }

        // To avoid division by zero when the trace of the stress is used
        I1 = Stress(0,0) + Stress(1,1) + Stress(2,2);
        if(fabs(I1) < 4.0){
            I1 = -4.0;
            Stress(0,0) = -1.0;
            Stress(1,1) = -1.0;
            Stress(2,2) = -2.0;
            Stress(0,1) = 0.0;
            Stress(0,2) = 0.0;
            Stress(1,2) = 0.0;
            Stress(1,0) = 0.0;
            Stress(2,0) = 0.0;
            Stress(2,1) = 0.0;
        }

        // Deviatoric stress tensor S
        //ShearStress = Stress - 1.0/3.0* I1 *OrthoSys::I;
        ShearStress = Stress;
        ShearStress(0,0) = ShearStress(0,0) - I1/3.;
        ShearStress(1,1) = ShearStress(1,1) - I1/3.;
        ShearStress(2,2) = ShearStress(2,2) - I1/3.;

        I1strain = StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2); // Trace of the strain tensor
        Ee   = sqrt(indpAA(StrainRate));         // The Frobenius or Euclidean norm of the strain tensor, sqrt(D^2)
        trOe = indpAB(Stress,StrainRate);        // Inner double product of stress and strain rate, tr(TD)

   
        // Prediction of the stress tensor    
        // Elastic prediction step (Sigma_e n+1)
        if (FirstStep){ Sigmaa = -dt/2.0*( SRT + RS + (c1hp*I1)*StrainRate + (c2hp*I1strain + c3hp*trOe/I1)*Stress + (c4hp*Ee)*(Stress + ShearStress) + FilStress) + Sigma;}

        Sigmab  = Sigmaa;
        Sigmaa  = dt*( SRT + RS + (c1hp*I1)*StrainRate + (c2hp*I1strain + c3hp*trOe/I1)*Stress + (c4hp*Ee)*(Stress + ShearStress) + FilStress) + Sigmaa;
        Sigma = 0.5*(Sigmaa+Sigmab);

        double press = (Sigma(0,0) + Sigma(1,1) + Sigma(2,2))/3.0;
        if(press > c){
            Sigma(0,0) -= (press - c);
            Sigma(1,1) -= (press - c);
            Sigma(2,2) -= (press - c);
        }

        // Prediction of the strain tensor
        if (FirstStep){ Straina = - dt/2.0*StrainRate + Strain;}
        
        Strainb = Straina;
        Straina = Straina + dt*StrainRate;
        Strain  = 0.5*(Straina+Strainb);
        // Notes:
        // FilStress: filter to smooth spurious numerical oscillations in the stress rate tensor
    }



    inline void Particle::Mat3epMVerlet(Mat3_t I, double dt)
    {
    	Mat3_t RotationRateT, Stress, SRT,RS;
    	double I1,J2,alpha,kf,I1strain;
    
    	// Jaumann rate terms
    	Trans(RotationRate,RotationRateT);
    	Mult(Sigma,RotationRateT,SRT);
    	Mult(RotationRate,Sigma,RS);

	// Volumetric strain
    	I1strain = StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2);
    
    	// Elastic prediction step (Sigma_e n+1)
    	Stress	= Sigma;
    	if (ct == 30)
    		Sigma	= dt*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I) + SRT + RS + FilStress) + Sigma;
    	else
    		Sigma	= 2.0*dt*( I1strain*K*OrthoSys::I + 2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I) + SRT + RS + FilStress) + Sigmab;
    	Sigmab	= Stress;
    
    	if (Fail>1)
    	{
    		if (I(2,2)==0.0)
    		{
    			// Drucker-Prager failure criterion for plane strain
    			alpha	= tan(phi) / sqrt(9.0+12.0*tan(phi)*tan(phi));
    			kf		= 3.0 * c  / sqrt(9.0+12.0*tan(phi)*tan(phi));
    		}
    		else
    		{
    			// Drucker-Prager failure criterion for 3D
    			alpha	= (2.0*  sin(phi)) / (sqrt(3.0)*(3.0-sin(phi)));
    			kf		= (6.0*c*cos(phi)) / (sqrt(3.0)*(3.0-sin(phi)));
    		}

    
    
    		// Bring back stress to the apex of the failure criteria
    		I1		= Sigma(0,0) + Sigma(1,1) + Sigma(2,2);
    		if ((kf-alpha*I1)<0.0)
    		{
    			double Ratio;
    			if (alpha == 0.0) Ratio =0.0; else Ratio = kf/alpha;
    			Sigma(0,0) -= 1.0/3.0*(I1-Ratio);
    			Sigma(1,1) -= 1.0/3.0*(I1-Ratio);
    			Sigma(2,2) -= 1.0/3.0*(I1-Ratio);
    			I1 			= Ratio;
    		}
    
    		// Shear stress based on the elastic assumption (S_e n+1)
    		ShearStress = Sigma - 1.0/3.0* I1 *OrthoSys::I;
    		J2 			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                        + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    
    
    		// Check the elastic prediction step by the failure criteria
    		if ((sqrt(J2)+alpha*I1-kf)>0.0)
    		{
    			// Shear stress based on the existing stress (S n)
    			ShearStress = Stress - 1.0/3.0*(Stress(0,0)+Stress(1,1)+Stress(2,2))*OrthoSys::I;
    			J2 			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                            + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    
    			if (sqrt(J2)>0.0)
    			{
    				Mat3_t temp, Plastic;
    				double sum,dLanda;
    
    				// calculating the plastic term based on the existing shear stress and strain rate
    				temp	= abab(ShearStress,StrainRate);
    				sum		= temp(0,0)+temp(0,1)+temp(0,2)+temp(1,0)+temp(1,1)+temp(1,2)+temp(2,0)+temp(2,1)+temp(2,2);
    				switch (Fail)
    				{
    				case 2:
    					dLanda	= 1.0/(9.0*alpha*alpha*K+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
    					Plastic	= 3.0*alpha*K*I + G/sqrt(J2)*ShearStress;
    					break;
    				case 3:
    					dLanda	= 1.0/(9.0*alpha*K*3.0*sin(psi)+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
    					Plastic	= 3.0*3.0*sin(psi)*K*I + G/sqrt(J2)*ShearStress;
    					break;
    				default:
    					std::cout << "Failure Type No is out of range. Please correct it and run again" << std::endl;
    					std::cout << "2 => Associated flow rule" << std::endl;
    					std::cout << "3 => non-associated flow rule" << std::endl;
    					abort();
    					break;
    				}
    				// Apply the plastic term
                    if (ct == 30)
                        Sigma = Sigma - dt*(dLanda*Plastic);
                    else
                        Sigma = Sigma - 2.0*dt*(dLanda*Plastic);
    			}
    
    			//Scale back
    			I1			= Sigma(0,0) + Sigma(1,1) + Sigma(2,2);
    			if ((kf-alpha*I1)<0.0)
    			{
    				double Ratio;
    				if (alpha == 0.0) Ratio =0.0; else Ratio = kf/alpha;
    				Sigma(0,0) -= 1.0/3.0*(I1-Ratio);
    				Sigma(1,1) -= 1.0/3.0*(I1-Ratio);
    				Sigma(2,2) -= 1.0/3.0*(I1-Ratio);
    				I1 			= Ratio;
    			}
    			ShearStress	= Sigma - 1.0/3.0* I1 *OrthoSys::I;
    			J2			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                            + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    
    			if ((sqrt(J2)+alpha*I1-kf)>0.0 && sqrt(J2)>0.0) Sigma = I1/3.0*OrthoSys::I + (kf-alpha*I1)/sqrt(J2) * ShearStress;
    		}
    	}
    
    	Stress	= Strain;
    	if (ct == 30)
    		Strain	= dt*StrainRate + Strain;
    	else
    		Strain	= 2.0*dt*StrainRate + Strainb;
    	Strainb	= Stress;
    
    	if (VarPorosity)
    	{
    		if (IsFree)
    		{
    			double ev = (Strain(0,0)+Strain(1,1)+Strain(2,2));
    			n = (n0+ev)/(1.0+ev);
    			switch(SeepageType)
    			{
    				case 0:
    					break;
    				case 1:
    					k = n*n*n*d*d/(180.0*(1.0-n)*(1.0-n));
    					break;
    				case 2:
    					k = n*n*n*d*d/(150.0*(1.0-n)*(1.0-n));
    					k2= 1.75*(1.0-n)/(n*n*n*d);
    					break;
    				case 3:
    					k = n*n*n*d*d/(150.0*(1.0-n)*(1.0-n));
    					k2= 0.4/(n*n*d);
    					break;
    				default:
    					std::cout << "Seepage Type No is out of range. Please correct it and run again" << std::endl;
    					std::cout << "0 => Darcy's Law" << std::endl;
    					std::cout << "1 => Darcy's Law & Kozeny–Carman Eq" << std::endl;
    					std::cout << "2 => The Forchheimer Eq & Ergun Coeffs" << std::endl;
    					std::cout << "3 => The Forchheimer Eq & Den Adel Coeffs" << std::endl;
    					abort();
    					break;
    			}
    		}
    		else
    			n = n0;
    	}
    	else
    		n = n0;
    
    
    }
    
    inline void Particle::ScalebackMat3(size_t Dimension,size_t Scheme)
    {
    	double I1,J2,alpha,kf;
    
    	if (Dimension==0.0)
    	{
    		// Drucker-Prager failure criterion for plane strain
    		alpha	= tan(phi) / sqrt(9.0+12.0*tan(phi)*tan(phi));
    		kf		= 3.0 * c  / sqrt(9.0+12.0*tan(phi)*tan(phi));
    	}
    	else
    	{
    		// Drucker-Prager failure criterion for 3D
    		alpha	= (2.0*  sin(phi)) / (sqrt(3.0)*(3.0-sin(phi)));
    		kf		= (6.0*c*cos(phi)) / (sqrt(3.0)*(3.0-sin(phi)));
    	}
    
    	// Bring back stressb to the apex of the failure criteria
    	I1	= Sigmab(0,0) + Sigmab(1,1) + Sigmab(2,2);
    	if ((kf-alpha*I1)<0.0)
    	{
    		double Ratio;
    		if (alpha == 0.0) Ratio =0.0; else Ratio = kf/alpha;
    		Sigmab(0,0) -= 1.0/3.0*(I1-Ratio);
    		Sigmab(1,1) -= 1.0/3.0*(I1-Ratio);
    		Sigmab(2,2) -= 1.0/3.0*(I1-Ratio);
    		I1 	     = Ratio;
    	}
    
    	ShearStress	= Sigmab - 1.0/3.0* I1 *OrthoSys::I;
    	J2			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                    + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    
    	if ((sqrt(J2)+alpha*I1-kf)>0.0 && sqrt(J2)>0.0) Sigmab = I1/3.0*OrthoSys::I + (kf-alpha*I1)/sqrt(J2) * ShearStress;
    
    	// Bring back stress to the apex of the failure criteria
    	if (Scheme == 0)
    	{
    		I1	= Sigma(0,0) + Sigma(1,1) + Sigma(2,2);
    		if ((kf-alpha*I1)<0.0)
    		{
    			double Ratio;
    			if (alpha == 0.0) Ratio =0.0; else Ratio = kf/alpha;
    			Sigma(0,0) -= 1.0/3.0*(I1-Ratio);
    			Sigma(1,1) -= 1.0/3.0*(I1-Ratio);
    			Sigma(2,2) -= 1.0/3.0*(I1-Ratio);
    			I1 	    = Ratio;
    		}
    
    		ShearStress	= Sigma - 1.0/3.0* I1 *OrthoSys::I;
    		J2			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                        + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    
    		if ((sqrt(J2)+alpha*I1-kf)>0.0 && sqrt(J2)>0.0) Sigma = I1/3.0*OrthoSys::I + (kf-alpha*I1)/sqrt(J2) * ShearStress;
    	}
    	else
    	{
    		I1	= Sigmaa(0,0) + Sigmaa(1,1) + Sigmaa(2,2);
    		if ((kf-alpha*I1)<0.0)
    		{
    			double Ratio;
    			if (alpha == 0.0) Ratio =0.0; else Ratio = kf/alpha;
    			Sigmaa(0,0) -= 1.0/3.0*(I1-Ratio);
    			Sigmaa(1,1) -= 1.0/3.0*(I1-Ratio);
    			Sigmaa(2,2) -= 1.0/3.0*(I1-Ratio);
    			I1 			= Ratio;
    		}
    		ShearStress	= Sigmaa - 1.0/3.0* I1 *OrthoSys::I;
    		J2			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                        + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );

    		if ((sqrt(J2)+alpha*I1-kf)>0.0 && sqrt(J2)>0.0) Sigmaa = I1/3.0*OrthoSys::I + (kf-alpha*I1)/sqrt(J2) * ShearStress;
    	}
    }
    
    
    inline void Particle::Mat3Leapfrog(Mat3_t I, double dt)
    {
    	Mat3_t RotationRateT, Stress, SRT,RS;
    	double I1,J2,alpha,kf,I1strain;
    
    	// Jaumann rate terms
    	Trans(RotationRate,RotationRateT);
    	Mult(Sigma,RotationRateT,SRT);
    	Mult(RotationRate,Sigma,RS);
    
    	// Volumetric strain
    	I1strain = StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2);
    
    	// Elastic prediction step (Sigma_e n+1)
    	if (FirstStep)
    		Sigmaa	= -dt/2.0*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigma;
    
    	Sigmab	= Sigmaa;
    	Sigmaa	= dt*(I1strain*K*OrthoSys::I + 2.0*G*(StrainRate-1.0/3.0*I1strain*OrthoSys::I) + SRT + RS) + Sigmaa;
    
    	if (Fail>1)
    	{
    		if (I(2,2)==0.0)
    		{
    			// Drucker-Prager failure criterion for plane strain
    			alpha	= tan(phi) / sqrt(9.0+12.0*tan(phi)*tan(phi));
    			kf		= 3.0 * c  / sqrt(9.0+12.0*tan(phi)*tan(phi));
    		}
    		else
    		{
    			// Drucker-Prager failure criterion for 3D
    			alpha	= (2.0*  sin(phi)) / (sqrt(3.0)*(3.0-sin(phi)));
    			kf		= (6.0*c*cos(phi)) / (sqrt(3.0)*(3.0-sin(phi)));
    		}
    
    		// Bring back stress to the apex of the failure criteria
    		I1		= Sigmaa(0,0) + Sigmaa(1,1) + Sigmaa(2,2);
    		if ((kf-alpha*I1)<0.0)
    		{
    			double Ratio;
    			if (alpha == 0.0) Ratio =0.0; else Ratio = kf/alpha;
    			Sigmaa(0,0) -= 1.0/3.0*(I1-Ratio);
    			Sigmaa(1,1) -= 1.0/3.0*(I1-Ratio);
    			Sigmaa(2,2) -= 1.0/3.0*(I1-Ratio);
    			I1 			= Ratio;
    		}
    
    		// Shear stress based on the elastic assumption (S_e n+1)
    		ShearStress = Sigmaa - 1.0/3.0* I1 *OrthoSys::I;
    		J2 			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                        + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    
    
    		// Check the elastic prediction step by the failure criteria
    		if ((sqrt(J2)+alpha*I1-kf)>0.0)
    		{
    			// Shear stress based on the existing stress (S n)
    			ShearStress = Sigma - 1.0/3.0*(Sigma(0,0)+Sigma(1,1)+Sigma(2,2))*OrthoSys::I;
    			J2 			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                            + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );
    
    			if (sqrt(J2)>0.0)
    			{
    				Mat3_t temp, Plastic;
    				double sum,dLanda;
    
    				// calculating the plastic term based on the existing shear stress and strain rate
    				temp	= abab(ShearStress,StrainRate);
    				sum		= temp(0,0)+temp(0,1)+temp(0,2)+temp(1,0)+temp(1,1)+temp(1,2)+temp(2,0)+temp(2,1)+temp(2,2);
    				switch (Fail)
    				{
    				case 2:
    					dLanda	= 1.0/(9.0*alpha*alpha*K+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
    					Plastic	= 3.0*alpha*K*I + G/sqrt(J2)*ShearStress;
    					break;
    				case 3:
    					dLanda	= 1.0/(9.0*alpha*K*3.0*sin(psi)+G)*( (3.0*alpha*K*I1strain) + (G/sqrt(J2))*sum );
    					Plastic	= 3.0*3.0*sin(psi)*K*I + G/sqrt(J2)*ShearStress;
    					break;
    				default:
    					std::cout << "Failure Type No is out of range. Please correct it and run again" << std::endl;
    					std::cout << "2 => Associated flow rule" << std::endl;
    					std::cout << "3 => non-associated flow rule" << std::endl;
    					abort();
    					break;
    				}
    				Sigmaa = Sigmaa - dt*(dLanda*Plastic);
    			}
    
    			I1	= Sigmaa(0,0) + Sigmaa(1,1) + Sigmaa(2,2);
    			if ((kf-alpha*I1)<0.0)
    			{
    				double Ratio;
    				if (alpha == 0.0) Ratio =0.0; else Ratio = kf/alpha;
    				Sigmaa(0,0) -= 1.0/3.0*(I1-Ratio);
    				Sigmaa(1,1) -= 1.0/3.0*(I1-Ratio);
    				Sigmaa(2,2) -= 1.0/3.0*(I1-Ratio);
    				I1 			= Ratio;
    			}
    			ShearStress	= Sigmaa - 1.0/3.0* I1 *OrthoSys::I;
    			J2			= 0.5*( (ShearStress(0,0)*ShearStress(0,0) + ShearStress(1,1)*ShearStress(1,1) + ShearStress(2,2)*ShearStress(2,2)) 
                            + 2.0*(ShearStress(0,1)*ShearStress(1,0) + ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,2)*ShearStress(2,1)) );

    			if ((sqrt(J2)+alpha*I1-kf)>0.0 && sqrt(J2)>0.0) Sigmaa = I1/3.0*OrthoSys::I + (kf-alpha*I1)/sqrt(J2) * ShearStress;
    		}
    	}
    	Sigma = 1.0/2.0*(Sigmaa+Sigmab);
    
    	if (FirstStep)
    		Straina	= -dt/2.0*StrainRate + Strain;
    	Strainb	= Straina;
    	Straina	= dt*StrainRate + Straina;
    	Strain	= 1.0/2.0*(Straina+Strainb);
    
    	if (VarPorosity)
    	{
    		if (IsFree)
    		{
    			double ev = (Strain(0,0)+Strain(1,1)+Strain(2,2));
    			n = (n0+ev)/(1.0+ev);
    			switch(SeepageType)
    			{
    				case 0:
    					break;
    				case 1:
    					k = n*n*n*d*d/(180.0*(1.0-n)*(1.0-n));
    					break;
    				case 2:
    					k = n*n*n*d*d/(150.0*(1.0-n)*(1.0-n));
    					k2= 1.75*(1.0-n)/(n*n*n*d);
    					break;
    				case 3:
    					k = n*n*n*d*d/(150.0*(1.0-n)*(1.0-n));
    					k2= 0.4/(n*n*d);
    					break;
    				default:
    					std::cout << "Seepage Type No is out of range. Please correct it and run again" << std::endl;
    					std::cout << "0 => Darcy's Law" << std::endl;
    					std::cout << "1 => Darcy's Law & Kozeny–Carman Eq" << std::endl;
    					std::cout << "2 => The Forchheimer Eq & Ergun Coeffs" << std::endl;
    					std::cout << "3 => The Forchheimer Eq & Den Adel Coeffs" << std::endl;
    					abort();
    					break;
    			}
    		}
    		else
    			n = n0;
    	}
    	else
    		n = n0;
    
    
    }
    

    inline void Particle::translate(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin)
    {
    	x = x + dt*v + 0.5*dt*dt*a;
    
    	// Evolve velocity
    	Vec3_t temp;
    	temp = v;
    	v = vb + 2*dt*a;
    	vb = temp;
    
    	//Periodic BC particle position update
    	if (Domainsize(0)>0.0)
    	{
    		(x(0)>(domainmax(0))) ? x(0) -= Domainsize(0) : x(0);
    		(x(0)<(domainmin(0))) ? x(0) += Domainsize(0) : x(0);
    	}
    	if (Domainsize(1)>0.0)
    	{
    		(x(1)>(domainmax(1))) ? x(1) -= Domainsize(1) : x(1);
    		(x(1)<(domainmin(1))) ? x(1) += Domainsize(1) : x(1);
    	}
    	if (Domainsize(2)>0.0)
    	{
    		(x(2)>(domainmax(2))) ? x(2) -= Domainsize(2) : x(2);
    		(x(2)<(domainmin(2))) ? x(2) += Domainsize(2) : x(2);
    	}
    }
    
}; // namespace SPH

#endif //SPH_PARTICLE_H
