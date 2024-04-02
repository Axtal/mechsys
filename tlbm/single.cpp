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

// Single fracture flow


// MechSys
#include <mechsys/lbm/Domain.h>


struct UserData
{
    double        rhomax;
    double        rhomin;
    double            pf;
    double            dt;
    double            dx;
    std::ofstream oss_ss;      
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
	for (size_t i=0; i<dom.Lat[0].Ndim(1); ++i)
	for (size_t j=0; j<dom.Lat[0].Ndim(2); ++j)
	{
        Cell * c = dom.Lat[0].GetCell(iVec3_t(0,i,j));
        if (c->IsSolid) continue;
        double * f = c->F;
        
        f[1] = 1.0/3.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-f[2]-2*f[3]-2*f[4]-2*f[5]-2*f[6]-4*f[8]+2*dat.rhomax);
        f[7] = 1.0/24.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-4*f[2] +f[3]-5*f[4]  +f[5]-5*f[6]+20*f[8]+2*dat.rhomax);
        f[9] = 1.0/24.0*(-2*f[0]+20*f[10]-4*f[12]-4*f[14]-4*f[2]+f[3]-5*f[4]-5*f[5]+f[6]-4*f[8]+2*dat.rhomax);
        f[11]= 1.0/24.0*(-2*f[0]-4*f[10]+20*f[12]-4*f[14]-4*f[2]-5*f[3]+f[4]  +f[5]-5*f[6]-4*f[8]+2*dat.rhomax);
        f[13]= 1.0/24.0*(-2*f[0]-4*f[10]-4 *f[12]+20*f[14]-4*f[2]-5*f[3]+  f[4]-5*f[5]+f[6]-4*f[8]+2*dat.rhomax);
        c->Rho = c->VelDen(c->Vel);            
	}

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    #endif
	for (size_t i=0; i<dom.Lat[0].Ndim(1); ++i)
	for (size_t j=0; j<dom.Lat[0].Ndim(2); ++j)
	{
        Cell * c = dom.Lat[0].GetCell(iVec3_t(dom.Lat[0].Ndim(0)-1,i,j));
        if (c->IsSolid) continue;
        double * f = c->F;
        f[2] = 1/3.0* (-2*f[0]-f[1]-2*(2*f[11]+2*f[13]+f[3]+f[4]+f[5]+f[6]+2*f[7]+2*f[9]-dat.rhomin));
        f[8] = 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] - 5*f[5] + f[6] +20*f[7] - 4*f[9] + 2*dat.rhomin);
        f[10]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] + f[5] - 5*f[6] - 4*f[7] + 20*f[9] + 2*dat.rhomin) ;
        f[12]= 1/24.0*(-2*f[0] - 4*f[1] + 20*f[11] - 4*f[13] + f[3] - 5*f[4] - 5*f[5] + f[6] -  4*f[7] - 4*f[9] + 2*dat.rhomin);
        f[14]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] + 20*f[13] + f[3] - 5*f[4] + f[5] - 5*f[6] -  4*f[7] - 4*f[9] + 2*dat.rhomin);
        c->Rho = c->VelDen(c->Vel);            
	}
}

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_per.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "V_ave" << Util::_8s << "V_th" << Util::_8s << "K \n";
    }
    double vave = 0.0;
    for (size_t k=1;k<dom.Lat[0].Ndim(2)-1;k++)
    {
        vave += dom.Lat[0].GetCell(iVec3_t(dom.Lat[0].Ndim(0)/2,dom.Lat[0].Ndim(1)/2,k))->Vel(0);
    }
    vave /= dom.Lat[0].Ndim(2)-2;

    double kave = dat.pf*dat.dx*dat.dt*(dom.Lat[0].Tau-0.5)/(dat.rhomax-dat.rhomin)*vave*(dom.Lat[0].Ndim(0));
    //double kave = dat.pf/(dat.rhomax-dat.rhomin)*vave*(dom.Lat[0].Ndim(0));
    
    double vth = (dom.Lat[0].Ndim(2)-2)*(dom.Lat[0].Ndim(2)-2)/(12.0*(dom.Lat[0].Tau-0.5))*(dat.rhomax-dat.rhomin)/(dom.Lat[0].Ndim(0));

    dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << vave << Util::_8s << vth << Util::_8s << kave << std::endl;
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    double Pf    = 1.0;
    double nu    = 0.16;
    double dx    = 1.0;
    double dt    = 1.0;
    if (argc>=2) Nproc=atoi(argv[1]);
    if (argc>=3) Pf   =atof(argv[2]);
    if (argc>=4) nu   =atof(argv[3]);
    if (argc>=5) dx   =atof(argv[4]);
    if (argc>=6) dt   =atof(argv[5]);
    size_t nx = 100;
    size_t ny = 6;
    size_t nz = 6;
    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    
    UserData dat;
    Dom.UserData = &dat;

    dat.rhomin  = 1.0;
    dat.rhomax  = 1.03;
    //dat.rhomax  = 1.1;
    dat.pf      = Pf;
    dat.dt      = dt;
    dat.dx      = dx;

    //Assigning solid boundaries at top and bottom
    //for (size_t i=0;i<nx;i++)
    //for (size_t j=0;j<ny;j++)
    //{
      //  Dom.Lat[0].GetCell(iVec3_t(i,j,0   ))->IsSolid = true;
      //  Dom.Lat[0].GetCell(iVec3_t(i,j,nz-1))->IsSolid = true;
    //}


    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        iVec3_t idx(ix,iy,iz);
        Dom.Lat[0].GetCell(idx)->Initialize(1.0,OrthoSys::O);
        Dom.Lat[0].GetCell(idx)->Pf = dat.pf;
    }  
    
    Dom.Solve(40000.0*dt,400.0*dt,Setup,Report,"single",true,Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH
