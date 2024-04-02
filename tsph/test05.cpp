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

///////////////////////////// test 05 SPH-DEM test case in 3D

#include <mechsys/sph/Domain.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

void UserOutput (SPH::Particle * Pa, double & Prop1, double & Prop2, double & Prop3)
{
    Prop1 = -(Pa->Sigma(0,0)+Pa->Sigma(1,1)+Pa->Sigma(2,2))/3.0;
    Prop2 = 0.0;
    Prop3 = 0.0;
}

int main(int argc, char **argv) try
{

        size_t Nproc = omp_get_max_threads();

        if (argc>1) Nproc = atoi(argv[1]);

        SPH::Domain		dom;
        dom.UserOutput = &UserOutput;
        dom.OutputName[0] = "PressureSoil";

        dom.Dimension	= 3;
        dom.Nproc		= Nproc;
    	dom.KernelType	= 0;
    	dom.Gravity		= 0.0, -9.81, 0.0;
    	//dom.Gravity		= 0.0, 0.0,-9.81;
    	dom.Scheme		= 0;

    	double H	= 0.1;
    	double W	= 0.1;
    	double L	= 0.1;
    	double n	= 40.0;
    	double rho	= 2650.0;
    	double K	= 10.0e6;
    	double Phi	= 20.0;
    	double c	= 5.0e3;
    	double Psi	= 0.0;
    	double Nu	= 0.3;
	
        double Mu 	= 1.002e-3;
	    double m	= 30.0;
	    double T0 	= 4.0e-4;

    	double dx	= H / n;
    	double h	= dx*1.2;
    	double E	= (3.0*(1.0-2.0*Nu))*K;
    	double G	= E/(2.0*(1.0+Nu));
        double Cs	= sqrt(E/rho);

    	dom.InitialDist	= dx;
        double timestep		= (0.2*h/(Cs));
        

        cout<<"Time Step = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"G = "<<G<<endl;
        cout<<"E = "<<E<<endl;
        cout<<"K = "<<K<<endl;
        cout<<"h = "<<h<<endl;

     	dom.AddBoxLength(1 ,Vec3_t (0.0,  0.0, 0.0), L, W, H, dx/2.0, rho, h, 1, 0, false, false);

     	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
			    dom.Particles[a]->Alpha		= 0.5;
			    dom.Particles[a]->Beta		= 0.5;
    		    dom.Particles[a]->Cs		= Cs;
    		    dom.Particles[a]->G			= G;
    		    dom.Particles[a]->K			= K;
    		    dom.Particles[a]->Material	= 3;
    		    dom.Particles[a]->Fail		= 3;
			    dom.Particles[a]->TI		= 0.0;
			    dom.Particles[a]->TIn		= 0.0;
    		    dom.Particles[a]->c			= c;
    		    dom.Particles[a]->phi		= Phi/180.0*M_PI;
    		    dom.Particles[a]->psi		= Psi/180.0*M_PI;
    	}

        dom.DomMax = Vec3_t(L+dx,W+dx,H+dx);
        dom.DomMin = Vec3_t(-dx,-dx,-dx);

        //dom.AddPlane(-1,Vec3_t(0.5*L,0.5*W,-dx),dx,2.0*L+2.0*dx,2.0*W+2.0*dx,1000.0,0.0,&OrthoSys::e0);
        dom.AddPlane(-1,Vec3_t(0.5*L,-dx,0.5*H),dx,2.0*L+2.0*dx,2.0*W+2.0*dx,1000.0,M_PI/2.0,&OrthoSys::e0);
        
        for (size_t j = 0;j<dom.DEMParticles.Size();j++)
        {
            dom.DEMParticles[j]->Props.eps = 0.0;
    	    for (size_t i=0; i<dom.Particles.Size(); i++)
    	    {
                double dist = DEM::Distance(dom.Particles[i]->x,*dom.DEMParticles[j]->Faces[0]);
                if (dist<dom.DEMParticles[j]->Props.eps || dom.DEMParticles[j]->Props.eps == 0.0) dom.DEMParticles[j]->Props.eps = dist;
                
            }
        }

    	dom.Solve(/*tf*/1000.0,/*dt*/timestep,/*dtOut*/0.005,"test05",200);
        return 0;
}
MECHSYS_CATCH
