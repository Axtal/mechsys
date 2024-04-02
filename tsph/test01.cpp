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

///////////////////////////// test 01 elastoplastic soil collapse

#include <mechsys/sph/Domain.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

// psi*sqrt(E/(rho*h*h))
double check = 0.0;
double Damp;
double t = 0.02;
void UserDamping(SPH::Domain & domi)
{
	if (domi.Time<t)
    #pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
		if (domi.Particles[i]->IsFree) domi.Particles[i]->a -= Damp * domi.Particles[i]->v;
	if (domi.Time>t && domi.Time<(t+1.0*domi.deltat) && check == 0.0)
	{
//		domi.TI	= 0.5;
//		domi.TIn= 2.55;
		std::cout<<"yes"<<std::endl;
		domi.DelParticles(5);
		check = 1.0;
	}
}

int main(int argc, char **argv) try
{

        size_t Nproc = omp_get_max_threads();
        String filekey  (argv[1]);
        String filename (filekey+".inp");
        if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
        std::ifstream infile(filename.CStr());

        if (argc>2) Nproc = atoi(argv[2]);

        SPH::Domain		dom;

        dom.Dimension	= 2;
        dom.Nproc		= Nproc;
    	dom.KernelType	= 0;
    	dom.Gravity		= 0.0, -9.81, 0.0;
    	dom.Scheme		= 0;

        double dx,h,rho,K,G,timestep,E,Nu,H,n,L,Phi,Psi,c,Cs,Mu,m,T0;

    	H	= 0.1;
    	n	= 40.0;
    	L	= 0.2;
    	rho	= 2650.0;
    	K	= 0.7e4;
    	Phi	= 19.8;
    	c	= 0.0;
    	Psi	= 0.0;
    	Nu	= 0.3;
	
        Mu 	= 1.002e-3;
	    m	= 30.0;
	    T0 	= 4.0e-4;
//        dom.Cs			= sqrt(E/rho);

        
        size_t Case = 0;
        infile >> Case;       infile.ignore(200,'\n');
        infile >> Mu;         infile.ignore(200,'\n');
        infile >> T0;         infile.ignore(200,'\n');
        infile >> Phi;        infile.ignore(200,'\n');
        infile >> c;          infile.ignore(200,'\n');
        infile >> K;          infile.ignore(200,'\n');
        infile >> Nu;         infile.ignore(200,'\n');
        infile >> n;          infile.ignore(200,'\n');
        infile >> H;          infile.ignore(200,'\n');
        infile >> L;          infile.ignore(200,'\n');

    	dx	= H / n;
    	h	= dx*1.2;
    	E	= (3.0*(1.0-2.0*Nu))*K;
    	G	= E/(2.0*(1.0+Nu));
        Cs	= sqrt(E/rho);
    	dom.InitialDist	= dx;
        timestep		= (0.2*h/(Cs));
        

//    	Damp 			= 0.02*sqrt(E/(rho*h*h));
//        dom.GeneralAfter= & UserDamping;

        cout<<"Time Step = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"G = "<<G<<endl;
        cout<<"E = "<<E<<endl;
        cout<<"K = "<<K<<endl;
        cout<<"h = "<<h<<endl;

     	dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx ,  0.0    , 0.0 ), 2.5*L + 6.0*dx + dx/10.0 , 1.1*H  + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
     	dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx , -3.0*dx , 0.0 ), 2.5*L + 3.0*dx + dx/10.0 , 3.0*dx + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );

		double K0 = 1-sin(Phi/180.0*M_PI);

     	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
            if (Case == 1)
            {
			    dom.Particles[a]->Alpha		= 0.1;
			    dom.Particles[a]->Beta		= 0.1;
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

            if (Case == 0)
            {
		        dom.Particles[a]->Cs		= Cs;
		        dom.Particles[a]->Mu 		= Mu;
		        dom.Particles[a]->MuRef 	= Mu;
		        dom.Particles[a]->T0 		= T0;
   		        dom.Particles[a]->m 		= m;
   		        dom.Particles[a]->VisM 		= 0;
		        //dom.Particles[a]->P0		= Cs*Cs*rho*0.01;
		        dom.Particles[a]->P0		= 0.0;
		        dom.Particles[a]->PresEq	= 0;
		        dom.Particles[a]->Material	= 1;
            }
    		if (dom.Particles[a]->x(0)>L&&dom.Particles[a]->x(0)<2.45*L&&dom.Particles[a]->x(1)>0.0)
    		//if (dom.Particles[a]->x(0)>L&&dom.Particles[a]->x(1)>0.0)
    		{
    			dom.Particles[a]->NoSlip	= true;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 5;
    		}
    		if (dom.Particles[a]->x(1)<0.0)
    		{
    			dom.Particles[a]->NoSlip	= true;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 3;
    		}
    		if (dom.Particles[a]->x(0)<0.0 && dom.Particles[a]->x(1)>0.0)
    		{
    			dom.Particles[a]->NoSlip	= true;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 3;
    		}
    		if (dom.Particles[a]->x(0)>2.45*L)
    		{
    			dom.Particles[a]->NoSlip	= true;
    			dom.Particles[a]->IsFree	= false;
    			dom.Particles[a]->ID		= 3;
    		}
    		//if (dom.Particles[a]->x(1)>1.2*H)
    		//{
    			//dom.Particles[a]->NoSlip	= true;
    			//dom.Particles[a]->IsFree	= false;
    			//dom.Particles[a]->ID		= 3;
    		//}
    		if (dom.Particles[a]->x(1)>H && dom.Particles[a]->ID==1)
    		{
    			dom.Particles[a]->ID		= 6;
    		}
         	for (size_t i=0; i<10; i++)
    			if (dom.Particles[a]->ID == 1)
    			{
    				if (dom.Particles[a]->x(0)>(7.0*dx + i*8.0*dx + dx/2.1) && dom.Particles[a]->x(0)<=(8.0*dx + i*8.0*dx + dx/2.1))
    					dom.Particles[a]->ID = 2;
    			}

         	for (size_t i=0; i<5; i++)
    			if (dom.Particles[a]->ID == 1)
    			{
    				if (dom.Particles[a]->x(1)>(7.0*dx + i*8.0*dx + dx/2.1) && dom.Particles[a]->x(1)<=(8.0*dx + i*8.0*dx + dx/2.1))
    					dom.Particles[a]->ID = 2;
    			}

//			if (dom.Particles[a]->ID == 1 || dom.Particles[a]->ID == 2)
//			{
//				dom.Particles[a]->Sigma(1,1)	= rho * -9.81 *(H-dom.Particles[a]->x(1));
//				dom.Particles[a]->Sigma(0,0)	= K0 * rho * -9.81 *(H-dom.Particles[a]->x(1));
//				dom.Particles[a]->Sigma(2,2)	= K0 * rho * -9.81 *(H-dom.Particles[a]->x(1));
//			}
    	}
     	dom.DelParticles(5);
     	dom.DelParticles(6);

    	dom.Solve(/*tf*/1000.0,/*dt*/timestep,/*dtOut*/0.005,"test01",999);
        return 0;
}
MECHSYS_CATCH
