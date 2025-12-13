/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang
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
// Elastic beam vibration.


// MechSys
#include <mechsys/mpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    Array<MPM::Particle *> EndBeamPar; //Particles at the end of the beam
    double                       Fb; //Maximun applied body force
    double                       Nc; //Number of cycles for the vibration
    double                       Tf; //Final simulation time
    size_t                     ndiv;
    size_t                       m;
    double                       Lx;
    double                       E;
    double                       Lz;
    double                       Lt;
    double                       B1;
    double                       w1;
    double                       Bc;
    double                       a;
    std::ofstream                oss_ss1; // The output file for position comparison
    std::ofstream                oss_ss2; // The output file for Energy
};

void Report (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    double x    = dat.Lt;
    double W1   =(cosh(dat.B1*x)-cos(dat.B1*x))+(cos(dat.B1*dat.Lt)+cosh(dat.B1*dat.Lt))/(sin(dat.B1*dat.Lt)+sinh(dat.B1*dat.Lt))*(sin(dat.B1*x)-sinh(dat.B1*x));
    double A1   = 0.017/(dat.w1*W1*5.089877);
    


    if (dom.idx_out==0)
    {
        String fs1;
        fs1.Printf("tmpm02_Comparison.res");
        dat.oss_ss1.open(fs1.CStr());
        dat.oss_ss1 << Util::_10_6 << "Time" << Util::_8s << "x" << Util::_8s << "v" << Util::_8s << "v(mpm)" << "\n";
        
        String fs2;
        fs2.Printf("tmpm02_Energy.res");
        dat.oss_ss2.open(fs2.CStr());
        dat.oss_ss2 << Util::_10_6 << "t" << Util::_8s << "k" << Util::_8s << "p" << Util::_8s << "to" << "\n";
    }

    else 
    {
        double v = 0.0;
        double x = 0.0;
        for (size_t ip=0;ip<dat.EndBeamPar.Size();ip++)
        {
            x += dat.EndBeamPar[ip]->x(0);
            v += dat.EndBeamPar[ip]->x(2)-0.1;
        }
        x /= dat.EndBeamPar.Size();
        v /= dat.EndBeamPar.Size();
        double vt = A1*W1*sin(5*dat.w1*dom.Time);
        dat.oss_ss1 << Util::_10_6 << dom.Time << Util::_8s << x-dat.Bc+dat.a*dat.Bc << Util::_8s << vt << Util::_8s << v << "\n";
        
        double Ek = 0; //Kinetics 
        double Ep = 0; //Potential
        for (size_t ip=0; ip < dom.Particles.Size(); ip++)
        {
            Ek += 0.5*dom.Particles[ip]->m*dot(dom.Particles[ip]->v,dom.Particles[ip]->v);
            for (int i=0;i<3;i++)
            {
                for (int j=0;j<3;j++)
                {
                    Ep += 0.5*dom.Particles[ip]->S(i,j)*dom.Particles[ip]->E(i,j)*dom.Particles[ip]->V;
                }
            }

        }
        double Et = Ek+Ep; //Total Energy
        dat.oss_ss2 << Util::_10_6 << dom.Time << Util::_8s << Ek << Util::_8s << Ep << Util::_8s << Et << "\n";
    }
}

int main(int argc, char **argv) try
{
    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
   
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);
    size_t ndiv = 20; //number of divisions per x lenght 
    //size_t ndiv = 30; //number of divisions per x lenght 
    double K	= 10.0e4; //Bulk modulus
    double Nu	= 0.3; //Poisson ratio
    double	E	= (3.0*(1.0-2.0*Nu))*K; //Young modulus
    double	G	= E/(2.0*(1.0+Nu)); //Shear modulus
    double rho  = 3000.0; //density
    double Cs	= sqrt(E/rho); //Speed of sound
    double h    = 1.0/ndiv; // 
    double dt   = 0.1*h/Cs;
    double Lx   = 1.0; //length of the beam in x
    double Ly   = 0.2;
    double Lz   = 0.2;
    double Dx   = 3.0*Lx/ndiv;
    double Bc   = Dx;
    double a    = 0.348;
    double Lt   = Lx-Bc+a*Bc;

    double x;
    double m    = rho*Ly*Lz*1.0;
    double I    = pow(Ly*Lz,3.0)/12.0;
    double B1   = 0.59686*M_PI/Lt;
    double w1   = B1*B1*sqrt(E*I/m);

    dat.ndiv    = ndiv;
    dat.Lx      = Lx; 
    dat.m       = m;
    dat.E       = E;
    dat.Lz      = Lz;
    dat.Lt      = Lt;
    dat.B1      = B1;
    dat.w1      = w1;
    dat.Bc      = Bc;
    dat.a       = a;

    dom.AddRectangularBeam(-1, Vec3_t(0.0,0.0,0.0), Vec3_t(Lx,Ly,Lz), rho, ndiv); //Adding the beam
    //dom.AddRectangularBeamMesh(-1, Vec3_t(0.0,0.0,0.0), Vec3_t(Lx,Ly,Lz), rho, ndiv); //Adding the beam
    dom.ResizeDomain(Vec3_t(-0.5*Lx,-0.5*Lx,-0.5*Lx),Vec3_t(1.5*Lx,1.5*Lx,1.5*Lx), Dx); // Mesh size


    //Setting properties for the material points
    for (size_t ip=0; ip < dom.Particles.Size(); ip++)
    {
        dom.Particles[ip]->K = K;
        dom.Particles[ip]->G = G;
        if (dom.Particles[ip]->x(0)<Bc)
        {
            dom.Particles[ip]->FixVeloc();
            dom.Particles[ip]->Tag = -2;
        }
        if (dom.Particles[ip]->x(0)>Lx*(ndiv-0.5)/ndiv)
        {
            dom.Particles[ip]->Tag = -3;
            dat.EndBeamPar.Push(dom.Particles[ip]);
        }
        dom.Particles[ip]->v(2) = 0.01*((cosh(B1*(dom.Particles[ip]->x(0)-Bc+a*Bc))-cos(B1*(dom.Particles[ip]->x(0)-Bc+a*Bc)))+(cos(B1*Lt)+cosh(B1*Lt))/(sin(B1*Lt)+sinh(B1*Lt))*(sin(B1*(dom.Particles[ip]->x(0)-Bc+a*Bc))-sinh(B1*(dom.Particles[ip]->x(0)-Bc+a*Bc))));
    }

    //Setting properties for nodes
    for (size_t in=0; in < dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        if (xn(0)<Bc)
        {
            dom.Nodes[in].FixVeloc();
        }
    }

    dom.Gn = 0.0;
    double Tf = 80.0; // Final Time
    dat.Tf = Tf;     
    dat.Nc = 5.0; // number of cycles for vibration

    dom.Solve(Tf,dt,Tf/300,NULL,Report,"tmpm02",true,Nproc);

}
MECHSYS_CATCH

