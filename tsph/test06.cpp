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

///////////////////////////// test 06 SPH-DEM test case in 2D for fluids Poiselle flow

#include <mechsys/sph/Domain.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

//Simulation time
double Tf = 200.0;

int main(int argc, char **argv) try
{

        size_t Nproc = omp_get_max_threads();

        if (argc>1) Nproc = atoi(argv[1]);

        SPH::Domain		dom;
        dom.Dimension	= 2;
        dom.Nproc		= Nproc;
        dom.VisEq = 0;
    	dom.KernelType	= 0;
    	dom.Gravity		= 0.0001,0.0,0.0;
    	dom.Scheme		= 1;

    	double H	= 0.1;
    	double n	= 40.0;
    	double L	= 0.1;
    	double rho	= 1000.0;
        double Mu   = 1.0e-1;
        double Cs   = 0.07;
        //double P0   = Cs*Cs*rho;
        double P0   = 0.0;
       

	  	double dx	= H / n;
    	double h	= dx*1.1;

    	dom.InitialDist	= dx;
        double timestep		= (0.1*h/(Cs));
        dom.BC.Periodic[0] = true;

        cout<<"Time Step = "<<timestep<<endl;
        cout<<"Cs = " << Cs <<endl;
        cout<<"h  = " << h  <<endl;
        cout<<"dx = " << dx <<endl;

     	dom.AddBoxLength(1 ,Vec3_t (0.0,  0.0, 0.0), L, H,  0.0, dx/2.0, rho, h, 1, 0, false, false);

     	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		    dom.Particles[a]->Material	= 1;
                dom.Particles[a]->PresEq    = 0;
                dom.Particles[a]->Mu        = Mu;
                dom.Particles[a]->MuRef     = Mu;
    		    dom.Particles[a]->Cs		= Cs;
    		    dom.Particles[a]->P0		= P0;
    		    //dom.Particles[a]->Shepard   = true;
    	}

        //std::cout << dom.Particles[0]->x << std::endl;
        //std::cout << dom.Particles[dom.Particles.Size()-1]->x << std::endl;

        dom.AddSegment(-1,Vec3_t(-L-dx,-0.5*dx,0.0),Vec3_t(2.0*L+dx,-0.5*dx,0.0),dx,1000.0);
        dom.AddSegment(-2,Vec3_t(-L-dx,H+0.5*dx,0.0),Vec3_t(2.0*L+dx,H+0.5*dx,0.0),dx,1000.0);
        //dom.AddSegment(-3,Vec3_t(-L-dx,-1.5*dx,0.0),Vec3_t(2.0*L+dx,-1.5*dx,0.0),dx,1000.0);
        //dom.AddSegment(-4,Vec3_t(-L-dx,H+1.5*dx,0.0),Vec3_t(2.0*L+dx,H+1.5*dx,0.0),dx,1000.0);
        //dom.AddSegment(-5,Vec3_t(-L-dx,-2.5*dx,0.0),Vec3_t(2.0*L+dx,-2.5*dx,0.0),dx,1000.0);
        //dom.AddSegment(-6,Vec3_t(-L-dx,H+2.5*dx,0.0),Vec3_t(2.0*L+dx,H+2.5*dx,0.0),dx,1000.0);
        
        for (size_t j = 0;j<dom.DEMParticles.Size();j++)
        {
            dom.DEMParticles[j]->Props.eps = 0.0;
    	    for (size_t i=0; i<dom.Particles.Size(); i++)
    	    {
                double dist = DEM::Distance(dom.Particles[i]->x,*dom.DEMParticles[j]->Faces[0]);
                if (dist<dom.DEMParticles[j]->Props.eps || dom.DEMParticles[j]->Props.eps == 0.0) dom.DEMParticles[j]->Props.eps = dist;
                
            }
            dom.DEMParticles[j]->FixVeloc();
        }

        //dom.XSPH = 0.1;
        dom.GradientType = 1;
    	//dom.Solve(/*tf*/10.0*timestep,/*dt*/timestep,/*dtOut*/Tf/1000,"test06",10000000);
    	dom.Solve(/*tf*/Tf,/*dt*/timestep,/*dtOut*/Tf/100,"test06",10000000);
        return 0;
}
MECHSYS_CATCH
