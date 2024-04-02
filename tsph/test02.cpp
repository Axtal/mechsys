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

///////////////////////////// test 02 3D dam brake

#include <mechsys/sph/Domain.h>

int main(int argc, char **argv) try
{
    size_t Nproc    = 1;
    if (argc>1) Nproc = atoi(argv[1]);

    SPH::Domain		dom;

    dom.Dimension   = 3;
    dom.Nproc       = Nproc;
    dom.Scheme		= 0;
	dom.Viscosity_Eq_Set(Morris);
	dom.Kernel_Set(Qubic_Spline);
	dom.Gradient_Approach_Set(Squared_density);


    double xb,yb,zb,dx,Cs,h,Rho,H,L,TL,TH,t,res,Mu,g,TW,W;

    H   = 0.5;
    L   = H;
    TH  = 1.25*H;
    TL  = 1.5*H;
    W   = H/2.0;
    TW  = H;
    res = 40.0;

    g   = 9.81;
    Rho = 998.21;
    Mu	= 1.002e-3;
    dx  = H/res;
    h   = dx*1.2;
    Cs	= 20.0 * sqrt(g*H);
    t   = (0.2*h/Cs);

    dom.InitialDist 	= dx;
    dom.Gravity	      = 0.0, -g ,0.0 ;

    dom.AddBoxLength(1 ,Vec3_t ( -3.0*dx , -3.0*dx , -3.0*dx ), 7.0*dx + TL + dx/10.0 , 3.0*dx + TH + dx/10.0,  6.0*dx + TW + dx/10.0 , dx/2.0 ,Rho, h, 1 , 0 , false, false );

    for (size_t a=0; a<dom.Particles.Size(); a++)
    {
        xb=dom.Particles[a]->x(0);
        yb=dom.Particles[a]->x(1);
        zb=dom.Particles[a]->x(2);

        dom.Particles[a]->Cs      = Cs;
        dom.Particles[a]->PresEq	= 1;
     	dom.Particles[a]->Mu			= Mu;
     	dom.Particles[a]->MuRef		= Mu;
     	dom.Particles[a]->Material= 1;
        dom.Particles[a]->Shepard = false;


        if (xb<0.0)
        {
        		dom.Particles[a]->ID			= 2;
        		dom.Particles[a]->IsFree	= false;
        		dom.Particles[a]->NoSlip	= true;
        }

        if (yb<0.0)
        {
        		dom.Particles[a]->ID			= 2;
        		dom.Particles[a]->IsFree	= false;
        		dom.Particles[a]->NoSlip	= true;
        }

        if (zb<0.0)
        {
        		dom.Particles[a]->ID			= 2;
        		dom.Particles[a]->IsFree	= false;
        		dom.Particles[a]->NoSlip	= true;
        }

        if (xb>TL)
        {
        		dom.Particles[a]->ID			= 2;
        		dom.Particles[a]->IsFree	= false;
        		dom.Particles[a]->NoSlip	= true;
        }

        if (zb>TW)
        {
        	dom.Particles[a]->ID			= 2;
        	dom.Particles[a]->IsFree	= false;
        	dom.Particles[a]->NoSlip	= true;
        }

        if (yb>H && dom.Particles[a]->ID==1)
          dom.Particles[a]->ID=11;

        if (xb>L && dom.Particles[a]->ID==1)
          dom.Particles[a]->ID=11;

        if (zb>W && dom.Particles[a]->ID==1)
          dom.Particles[a]->ID=11;

        if (dom.Particles[a]->ID==1)
        {
          dom.Particles[a]->Density	  = Rho*pow((1+7.0*g*(H-yb)/(Cs*Cs)),(1.0/7.0));
          dom.Particles[a]->Densityb	= Rho*pow((1+7.0*g*(H-yb)/(Cs*Cs)),(1.0/7.0));
        }

    }

    dom.DelParticles(11);

    dom.Solve(/*tf*/50.0,/*dt*/t,/*dtOut*/0.05,"test02",999);
    return 0;
}
MECHSYS_CATCH
