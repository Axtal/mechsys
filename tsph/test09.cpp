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

///////////////////////////// test 09 SPH-DEM embedded pipe in 2d

#include <mechsys/sph/Domain.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

//Simulation time
double Tf = 2.0;

//Damping coefficient (changing from DampS0 to DampS1) in time from T0 to T1
double DampS  = 0.0;
double DampS0 = 0.0;
double DampS1 = 0.0;
double T0     = 0.2*Tf;
double T1     = 0.3*Tf;



void UserOutput (SPH::Particle * Pa, double & Prop1, double & Prop2, double & Prop3)
{
    Prop1 = -(Pa->Sigma(0,0)+Pa->Sigma(1,1)+Pa->Sigma(2,2))/3.0;
    Prop2 = 0.0;
    Prop3 = 0.0;
}

void After(SPH::Domain & dom)
{
    //if (dom.Time > T0)
    //{
        //DampS = std::max(DampS0*(T1-dom.Time)/(T1-T0)+DampS1*(dom.Time-T0)/(T1-T0),DampS1);
    //}

	//#pragma omp parallel for schedule (static) num_threads(dom.Nproc)
	//for (size_t i=0; i<dom.Particles.Size(); i++)
    //{
        //dom.Particles[i]->a -= DampS*dom.Particles[i]->v;
    //}
}

int main(int argc, char **argv) try
{

        size_t Nproc = omp_get_max_threads();

        if (argc>1) Nproc = atoi(argv[1]);

        SPH::Domain		dom;
        dom.UserOutput = &UserOutput;
        dom.OutputName[0] = "PressureSoil";

        dom.Dimension	= 2;
        dom.Nproc		= Nproc;
    	dom.KernelType	= 0;
    	dom.Gravity		= 0.0, -9.81, 0.0;
    	dom.Scheme		= 1;

    	double H	= 0.2;
    	double n	= 100.0;
    	double L	= 0.2;
    	double rho	= 2650.0;
    	double K	= 1.0e6;
    	double Phi	= 30.0;
    	double c	= 0.0e3;
    	double Psi	= 0.0;
    	double Nu	= 0.3;
        double Rc   = 0.1*H;
	

    	double dx	= H / n;
    	double h	= dx*1.2;
    	double E	= (3.0*(1.0-2.0*Nu))*K;
    	double G	= E/(2.0*(1.0+Nu));
        double Cs	= sqrt(E/rho);

    	dom.InitialDist	= dx;
        double timestep		= (0.1*h/(Cs));
        DampS	= 0.2*sqrt(E/(rho*h*h));
        DampS0  = 0.2*sqrt(E/(rho*h*h));
        DampS1  = 0.2*sqrt(E/(rho*h*h));

        cout<<"Time Step = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"G = "<<G<<endl;
        cout<<"E = "<<E<<endl;
        cout<<"K = "<<K<<endl;
        cout<<"h = "<<h<<endl;

     	dom.AddBoxLength(1 ,Vec3_t (0.0,  0.0, 0.0), L, H,  0.0, dx/2.0, rho, h, 1, 0, false, false);
        dom.AddSegment(-1,Vec3_t(-L-dx,-dx,0.0),Vec3_t(2.0*L+dx,-dx,0.0),dx,1000.0);
        dom.AddSegment(-2,Vec3_t(-dx,-dx,0.0),Vec3_t(-dx,H+dx,0.0),dx,1000.0);
        dom.AddSegment(-3,Vec3_t(L+dx,-dx,0.0),Vec3_t(L+dx,H+dx,0.0),dx,1000.0);
        dom.AddSegment(-4,Vec3_t(-dx,H+3.0*Rc,0.0),Vec3_t(0.4*L,H+2.5*Rc,0.0),dx,1000.0);
        dom.AddSphere (1, Vec3_t(0.5*L,1.0*H+Rc+2.0*dx,0.0),Rc,rho);
        dom.AddSphere (2, Vec3_t(0.2*L,1.0*H+4.0*Rc+2.0*dx,0.0),Rc,rho);
        dom.DEMParticles[0]->FixVeloc();
        dom.DEMParticles[1]->FixVeloc();
        dom.DEMParticles[2]->FixVeloc();
        dom.DEMParticles[3]->FixVeloc();

     	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
			    dom.Particles[a]->Alpha		= 5.0;
			    dom.Particles[a]->Beta		= 0.0;
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
    		    dom.Particles[a]->Shepard   = true;
		        //double xb=dom.Particles[a]->x(0)-0.5*L;
		        //double yb=dom.Particles[a]->x(1)-0.5*H;
		        //if ((xb*xb+yb*yb)<((Rc+dx)*(Rc+dx)))
		        //{
		        	//dom.Particles[a]->ID=4;
		        //}
    	}
	    //dom.DelParticles(4);


        dom.DomMax = Vec3_t(2.0*L+dx,H+dx,0.0);
        dom.DomMin = Vec3_t(-L-dx,-2.0*dx,0.0);

        for (size_t j = 0;j<dom.DEMParticles.Size();j++)
        {
            //dom.DEMParticles[j]->Props.Mu  = tan(Phi/180.0*M_PI);
            //dom.DEMParticles[j]->Props.Mu  = 0.1;
            dom.DEMParticles[j]->Ff = dom.DEMParticles[j]->Props.m*dom.Gravity;
            if (dom.DEMParticles[j]->Tag >= 1) continue;
            dom.DEMParticles[j]->Props.eps = 0.0;
    	    for (size_t i=0; i<dom.Particles.Size(); i++)
    	    {
                double dist = DEM::Distance(dom.Particles[i]->x,*dom.DEMParticles[j]->Faces[0]);
                if (dist<dom.DEMParticles[j]->Props.eps || dom.DEMParticles[j]->Props.eps == 0.0) dom.DEMParticles[j]->Props.eps = dist;
            }
        }

        dom.DEMParticles[0]->Props.Mu  = 0.0;
        dom.DEMParticles[1]->Props.Mu  = 0.0;
        dom.DEMParticles[2]->Props.Mu  = 0.0;
        dom.DEMParticles[3]->Props.Mu  = tan(Phi/180.0*M_PI);
        dom.DEMParticles[4]->Props.Mu  = tan(Phi/180.0*M_PI);
        dom.DEMParticles[4]->Props.eps = dx;
        dom.DEMParticles[4]->Props.Gn  = -0.2;
        dom.DEMParticles[5]->Props.Mu  = tan(Phi/180.0*M_PI);
        dom.DEMParticles[5]->Props.eps = dx;
        dom.DEMParticles[5]->Props.Gn  = -0.2;

        //dom.XSPH = 0.1;
        dom.GradientType = 1;
        dom.GeneralAfter = &After;
    	dom.Solve(/*tf*/Tf,/*dt*/timestep,/*dtOut*/Tf/200,"test09",10000000);
        return 0;
}
MECHSYS_CATCH
