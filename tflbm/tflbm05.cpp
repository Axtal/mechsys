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

// Single Fracture flow


// MechSys
#include <mechsys/flbm/Domain.h>


struct UserData
{
    double        rhomax;
    double        rhomin;
    double            dt;
    double            dx;
    std::ofstream oss_ss;      
    #ifdef USE_OCL
    cl::Buffer        bBCRho;
    cl::Program       UserProgram;
    #endif
};

void Setup (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    
    // Cells with prescribed velocity
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
	for (size_t i=0; i<dom.Ndim(1); ++i)
	for (size_t j=0; j<dom.Ndim(2); ++j)
	{
        if (dom.IsSolid[0][0][i][j]) continue;
        double * f = dom.F[0][0][i][j];
        
        f[1] = 1.0/3.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-f[2]-2*f[3]-2*f[4]-2*f[5]-2*f[6]-4*f[8]+2*dat.rhomax);
        f[7] = 1.0/24.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-4*f[2] +f[3]-5*f[4]  +f[5]-5*f[6]+20*f[8]+2*dat.rhomax);
        f[9] = 1.0/24.0*(-2*f[0]+20*f[10]-4*f[12]-4*f[14]-4*f[2]+f[3]-5*f[4]-5*f[5]+f[6]-4*f[8]+2*dat.rhomax);
        f[11]= 1.0/24.0*(-2*f[0]-4*f[10]+20*f[12]-4*f[14]-4*f[2]-5*f[3]+f[4]  +f[5]-5*f[6]-4*f[8]+2*dat.rhomax);
        f[13]= 1.0/24.0*(-2*f[0]-4*f[10]-4 *f[12]+20*f[14]-4*f[2]-5*f[3]+  f[4]-5*f[5]+f[6]-4*f[8]+2*dat.rhomax);

        dom.Vel[0][0][i][j] = OrthoSys::O;
        dom.Rho[0][0][i][j] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][0][i][j] +=  dom.F[0][0][i][j][k];
            dom.Vel[0][0][i][j] +=  dom.F[0][0][i][j][k]*dom.C[k];
        }
        dom.Vel[0][0][i][j] *= dom.Cs/dom.Rho[0][0][i][j];
	}

	// Cells with prescribed density
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
	for (size_t i=0; i<dom.Ndim(1); ++i)
	for (size_t j=0; j<dom.Ndim(2); ++j)
	{
        if (dom.IsSolid[0][dom.Ndim(0)-1][i][j]) continue;
        double * f = dom.F[0][dom.Ndim(0)-1][i][j];

        f[2] = 1/3.0* (-2*f[0]-f[1]-2*(2*f[11]+2*f[13]+f[3]+f[4]+f[5]+f[6]+2*f[7]+2*f[9]-dat.rhomin));
        f[8] = 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] - 5*f[5] + f[6] +20*f[7] - 4*f[9] + 2*dat.rhomin);
        f[10]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] + f[5] - 5*f[6] - 4*f[7] + 20*f[9] + 2*dat.rhomin) ;
        f[12]= 1/24.0*(-2*f[0] - 4*f[1] + 20*f[11] - 4*f[13] + f[3] - 5*f[4] - 5*f[5] + f[6] -  4*f[7] - 4*f[9] + 2*dat.rhomin);
        f[14]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] + 20*f[13] + f[3] - 5*f[4] + f[5] - 5*f[6] -  4*f[7] - 4*f[9] + 2*dat.rhomin);
        
        dom.Vel[0][dom.Ndim(0)-1][i][j] = OrthoSys::O;
        dom.Rho[0][dom.Ndim(0)-1][i][j] = 0.0;
        for (size_t k=0;k<dom.Nneigh;k++)
        {
            dom.Rho[0][dom.Ndim(0)-1][i][j] +=  dom.F[0][dom.Ndim(0)-1][i][j][k];
            dom.Vel[0][dom.Ndim(0)-1][i][j] +=  dom.F[0][dom.Ndim(0)-1][i][j][k]*dom.C[k];
        }
        dom.Vel[0][dom.Ndim(0)-1][i][j] *= dom.Cs/dom.Rho[0][dom.Ndim(0)-1][i][j];
	}
}

void Report (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.IsFirstTime)
    {
        dom.IsFirstTime = false;
        String fs;
        fs.Printf("%s_per.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "V_ave" << Util::_8s << "V_th \n";
    }
    double vave = 0.0;
    for (size_t k=1;k<dom.Ndim(2)-1;k++)
    {
        vave += dom.Vel[0][dom.Ndim(0)/2][dom.Ndim(1)/2][k](0);
    }
    vave /= dom.Ndim(2)-2;
    
    double vth = dat.dx/dat.dt*(dom.Ndim(2)-2)*(dom.Ndim(2)-2)/(12.0*(dom.Tau[0]-0.5))*(dat.rhomax-dat.rhomin)/(dom.Ndim(0));

    dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << vave << Util::_8s << vth << std::endl;
    //std::cout << Util::_10_6 << dom.Time << Util::_8s << vave << Util::_8s << vth << std::endl;
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1;
    double nu    = 1.0/6.0;
    double dx    = 1.0;
    double dt    = 1.0;
    if (argc>=2) Nproc=atoi(argv[1]);
    if (argc>=3) nu   =atof(argv[2]);
    if (argc>=4) dx   =atof(argv[3]);
    if (argc>=5) dt   =atof(argv[4]);
    size_t nx = 100;
    size_t ny = 6;
    size_t nz = 10;
    //size_t nx = 1;
    //size_t ny = 1;
    //size_t nz = 2;
    FLBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    
    UserData dat;
    Dom.UserData = &dat;
    Dom.Sc = 0.0;

    dat.rhomin  = 1.0;
    dat.rhomax  = 1.03;
    dat.dx      = dx;
    dat.dt      = dt;

    //Assigning solid boundaries at top and bottom
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        Dom.IsSolid[0][i][j][0]    = true;
        Dom.IsSolid[0][i][j][nz-1] = true;
    }


    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        Dom.Initialize(0,idx,1.0,v);
    }  

    Dom.Solve(1.0e4*dt,100.0*dt,Setup,Report,"single",true,Nproc);
    //Dom.Solve(1.0,80.0,NULL,NULL,"single",true,Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH
