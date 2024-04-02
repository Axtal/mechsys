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

///////////////////////////// test 08 SPH-DEM test case in 2D drag coefficient

#include <mechsys/sph/Domain.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

//Simulation time
double Tf = 5.0e3;
//Average Velocity
double U = 0.0;
void After(SPH::Domain & dom)
{
    dom.BC.inv = Vec3_t(std::min(10.0*dom.Time/Tf*U,U),0.0,0.0);
}

int main(int argc, char **argv) try
{
    size_t Nproc = omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);

    SPH::Domain		dom;
    dom.Nproc	   	    = Nproc;
	dom.Dimension		= 2;
	dom.BC.Periodic[1]	= true;
	dom.BC.Periodic[0]	= true;
	dom.VisEq		    = 0;
	dom.KernelType		= 0;
	dom.Scheme		    = 1;
    dom.Gravity		    = 1.0e-5,0.0,0.0;
    //dom.GradientType    = 1;
    //dom.GeneralAfter    = &After;


	double Mu		= 1.0e-1;
	double rho		= 1000.0;
    double H	    = 0.2;
    double L	    = H;
    double chi      = 0.2;
	double Rc		= H*sqrt(chi/M_PI);
    double n	    = 200.0;
	double dx		= H/n;
	double h		= dx*1.1;
	//double Re		= 1.0;
	//U		        = Re*Mu/(rho*2.0*Rc);
	double Cs		= 0.07;
	//double P0		= Cs*Cs*rho*1.0;
	double P0		= 0.0;
	double timestep	= (0.05*h/Cs);

	//dom.BC.InOutFlow	= 3;
	//dom.BC.allv		    = U,0.0,0.0;
	//dom.BC.inDensity	= rho;
	//dom.BC.inv		    = U,0.0,0.0;
//	dom.BC.outv		= U,0.0,0.0;
//	dom.BC.outDensity	= rho;
	dom.InitialDist 	= dx;
//	dom.BC.MassConservation = true;

	//std::cout<<"Re = "<<Re<<std::endl;
	std::cout<<"V  = "<<U<<std::endl;
	std::cout<<"Cs = "<<Cs<<std::endl;
	std::cout<<"P0 = "<<P0<<std::endl;
	std::cout<<"Time Step = "<<timestep<<std::endl;
	std::cout<<"Resolution = "<<(2.0*Rc/dx)<<std::endl;

	dom.AddBoxLength(1 ,OrthoSys::O, L , H  , 0.0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->Shepard	= true;
		//dom.Particles[a]->Alpha	    = 0.1;
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->P0		= P0;
		dom.Particles[a]->Mu		= Mu;
		dom.Particles[a]->MuRef		= Mu;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Material	= 1;
		double xb=dom.Particles[a]->x(0)-0.5*L;
		double yb=dom.Particles[a]->x(1)-0.5*H;
		if ((xb*xb+yb*yb)<((Rc+dx)*(Rc+dx)))
		{
			dom.Particles[a]->ID=4;
		}
	}
	dom.DelParticles(4);


    dom.AddSphere(-1,Vec3_t(0.5*L,0.5*H,0.0),Rc,rho);
    dom.DEMParticles[0]->Props.eps = dx;
    dom.DEMParticles[0]->FixVeloc();

	dom.Solve(/*tf*/Tf,/*dt*/timestep,/*dtOut*/Tf/100.0,"test08",10000000);
	return 0;
}
MECHSYS_CATCH
