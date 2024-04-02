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

// Cylinder drag coefficient

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    Array<Cell *> Left;
    Array<Cell *> Right;
    Array<double> Vel;
    double        vmax;
    double        rho;
    double        R;
    double        nu;
    double        Tf;
};

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_force.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "Tx" << Util::_8s << "Ty" << Util::_8s << "Tz" << Util::_8s << "Vx" << Util::_8s << "Fx" << Util::_8s << "F" << Util::_8s << "Rho" << Util::_8s << "Re" << Util::_8s << "CD" << Util::_8s << "CDsphere \n";
    }
    if (!dom.Finished) 
    {
        double M    = 0.0;
        double Vx   = 0.0;
        size_t nc   = 0;
        Vec3_t Flux = OrthoSys::O;
        for (size_t i=0;i<dom.Lat[0].Ncells;i++)
        {
            Cell * c = dom.Lat[0].Cells[i];
            if (c->IsSolid||c->Gamma>1.0e-8) continue;
            Vx   += (1.0 - c->Gamma)*c->Vel(0);
            Flux += c->Rho*c->Vel;
            M    += c->Rho;
            nc++;
        }
        Vx  /=dom.Lat[0].Ncells;
        Flux/=M;
        M   /=nc;
        double CD  = 2.0*dom.Disks[0]->F(0)/(Flux(0)*Flux(0)/M*2.0*dat.R);
        double Re  = 2.0*Flux(0)*dat.R/(M*dat.nu);
        double CDt = 10.7*pow(Re,-0.7239)+0.9402;
        if (Flux(0)<1.0e-12) 
        {
            Flux = OrthoSys::O;
            CD = 0.0;
            Re = 0.0;
        }

        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dom.Disks[0]->F(0) << Util::_8s << dom.Disks[0]->F(1) << Util::_8s << dom.Disks[0]->F(2) << Util::_8s << dom.Disks[0]->T(0) << Util::_8s << dom.Disks[0]->T(1) << Util::_8s << dom.Disks[0]->T(2) << Util::_8s << Vx << Util::_8s << Flux(0) << Util::_8s << norm(Flux) << Util::_8s << M << Util::_8s << Re << Util::_8s << CD << Util::_8s << CDt << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //for (size_t i=0;i<dat.Left.Size();i++)
    //{
        //dat.Left [i]->Initialize(dat.Left [i]->RhoBC,dat.Left [i]->VelBC);
        //dat.Right[i]->Initialize(dat.Right[i]->RhoBC,dat.Right[i]->VelBC);
    //}

    double vel = std::min(10.0*dat.vmax*dom.Time/dat.Tf,dat.vmax);
	// Cells with prescribed velocity
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
	for (size_t i=0; i<dat.Left.Size(); ++i)
	{
		Cell * c = dat.Left[i];
		if (c->IsSolid) continue;
		double rho = (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[3]+c->F[6]+c->F[7]))/(1.0-vel);
		//c->F[1] = c->F[3] + (2.0/3.0)*rho*dat.Vel[i];
		//c->F[5] = c->F[7] + (1.0/6.0)*rho*dat.Vel[i] - 0.5*(c->F[2]-c->F[4]);
		//c->F[8] = c->F[6] + (1.0/6.0)*rho*dat.Vel[i] + 0.5*(c->F[2]-c->F[4]);
		c->F[1] = c->F[3] + (2.0/3.0)*rho*vel;
		c->F[5] = c->F[7] + (1.0/6.0)*rho*vel - 0.5*(c->F[2]-c->F[4]);
		c->F[8] = c->F[6] + (1.0/6.0)*rho*vel + 0.5*(c->F[2]-c->F[4]);
        c->Rho = c->VelDen(c->Vel);
	}

	// Cells with prescribed density
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
	for (size_t i=0; i<dat.Right.Size(); ++i)
	{
		Cell * c = dat.Right[i];
		if (c->IsSolid) continue;
		//double vx = -1.0 + (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[1]+c->F[5]+c->F[8]))/c->RhoBC;
		//c->F[3] = c->F[1] - (2.0/3.0)*c->RhoBC*vx; 
		//c->F[7] = c->F[5] - (1.0/6.0)*c->RhoBC*vx + 0.5*(c->F[2]-c->F[4]);
		//c->F[6] = c->F[8] - (1.0/6.0)*c->RhoBC*vx - 0.5*(c->F[2]-c->F[4]);
		double rho = (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[1]+c->F[5]+c->F[8]))/(1.0+vel);
		c->F[3] = c->F[1] - (2.0/3.0)*rho*vel;
		c->F[7] = c->F[5] - (1.0/6.0)*rho*vel + 0.5*(c->F[2]-c->F[4]);
		c->F[6] = c->F[8] - (1.0/6.0)*rho*vel - 0.5*(c->F[2]-c->F[4]);
        c->Rho = c->VelDen(c->Vel);
	}
}

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc>=2) nproc=atoi(argv[1]);
    double u_max  = 0.1;                // Poiseuille's maximum velocity
    double Re     = 5000.0;                  // Reynold's number
    if (argc>=3) Re=atof(argv[2]);
    size_t nx = 2000;
    size_t ny = 2000;
    int radius = 21;           // radius of inner circle (obstacle)
    double nu     = u_max*(2*radius)/Re; // viscocity
    //double nu     = 1.0/6.0; // viscocity
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    Dom.Fconv = 1.0;
    UserData dat;
    Dom.UserData = &dat;
    dat.R = radius;
    dat.nu = nu;
    dat.Tf = 80000.0;
    Dom.Sc = 0.05;

    dat.vmax = u_max;
    //Assigning the left and right cells
    for (size_t i=0;i<ny;i++)
    {
        dat.Left .Push(Dom.Lat[0].GetCell(iVec3_t(0   ,i,0)));
        dat.Right.Push(Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0)));
        
        // set parabolic profile
        //double L  = ny - 2;                       // channel width in cell units
        //double yp = i - 1.5;                      // ordinate of cell
        //double vx = dat.vmax*4/(L*L)*(L*yp - yp*yp); // horizontal velocity
        double vx = dat.vmax; // horizontal velocity
        dat.Vel.Push(vx);
        dat.Left [i]->RhoBC = 1.0;
        //dat.Right[i]->VelBC = 0.08,0.0,0.0;
        dat.Right[i]->RhoBC = 1.0;

    }
    dat.rho  = 1.0;

	// set inner obstacle
	int obsX   = nx/4;   // x position
	int obsY   = ny/2+3; // y position
    //for (size_t i=0;i<nx;i++)
    //for (size_t j=0;j<ny;j++)
    //{
        //if ((i-obsX)*(i-obsX)+(j-obsY)*(j-obsY)<radius*radius)
        //{
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->IsSolid = true;
        //}
    //}

    Dom.AddDisk(-1,Vec3_t(obsX,obsY,0),Vec3_t(0.0,0.0,0.0),Vec3_t(0.0,0.0,0.0),   3.0,radius,1.0);
    Dom.Disks[0]->FixVeloc();
    //Dom.AddDisk(0,Vec3_t(nx/2.0,obsY,0),Vec3_t(0.0,0.0,0.0),Vec3_t(0.0,0.0,0.0),30.0,2*radius,1.0);



    //Assigning solid boundaries at top and bottom
    //for (size_t i=0;i<nx;i++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    //}

    double rho0 = 1.0;
    Vec3_t v0(0.08,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        //Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
        Dom.Lat[0].Cells[i]->Initialize(rho0, OrthoSys::O);
    }

    //Solving
    Dom.Time = 0.0;
    Dom.Solve(dat.Tf,dat.Tf/200.0,Setup,Report,"tlbm02",true,nproc);
 
}
MECHSYS_CATCH


