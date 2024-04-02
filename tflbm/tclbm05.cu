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
    thrust::device_vector<real> rhobc; 
    real                     * prhobc; 
    double            dt;
    double            dx;
    std::ofstream oss_ss;      
};

__global__ void Left_BC(real * PBC, bool * IsSolid, real * F, real3 * Vel, real * Rho, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ny*lbmaux[0].Nz) return;
    
    size_t ib  = ic*lbmaux[0].Nx;
    if (!IsSolid[ib])
    {
        size_t iv = ib*lbmaux[0].Nneigh;
        real * f = F + iv;
        f[1] = 1.0/3.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-f[2]-2*f[3]-2*f[4]-2*f[5]-2*f[6]-4*f[8]+2*PBC[0]);
        f[7] = 1.0/24.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-4*f[2] +f[3]-5*f[4]  +f[5]-5*f[6]+20*f[8]+2*PBC[0]);
        f[9] = 1.0/24.0*(-2*f[0]+20*f[10]-4*f[12]-4*f[14]-4*f[2]+f[3]-5*f[4]-5*f[5]+f[6]-4*f[8]+2*PBC[0]);
        f[11]= 1.0/24.0*(-2*f[0]-4*f[10]+20*f[12]-4*f[14]-4*f[2]-5*f[3]+f[4]  +f[5]-5*f[6]-4*f[8]+2*PBC[0]);
        f[13]= 1.0/24.0*(-2*f[0]-4*f[10]-4 *f[12]+20*f[14]-4*f[2]-5*f[3]+  f[4]-5*f[5]+f[6]-4*f[8]+2*PBC[0]);

        Rho   [ib] = 0.0;
        Vel   [ib] = make_real3(0.0,0.0,0.0);
        for(size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ib] += F[iv + k];
            Vel[ib] = Vel[ib] + F[iv + k]*lbmaux[0].C[k];
        }
        Vel[ib] = (lbmaux[0].Cs/Rho[ib])*Vel[ib];
    }
}

__global__ void Right_BC(real * PBC, bool * IsSolid, real * F, real3 * Vel, real * Rho, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Ny*lbmaux[0].Nz) return;
    
    size_t ib  = (ic+1)*lbmaux[0].Nx-1;
    if (!IsSolid[ib])
    {
        size_t iv = ib*lbmaux[0].Nneigh;
        real * f = F + iv;

        f[2] = 1/3.0* (-2*f[0]-f[1]-2*(2*f[11]+2*f[13]+f[3]+f[4]+f[5]+f[6]+2*f[7]+2*f[9]-PBC[1]));
        f[8] = 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] - 5*f[5] + f[6] +20*f[7] - 4*f[9] + 2*PBC[1]);
        f[10]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] + f[5] - 5*f[6] - 4*f[7] + 20*f[9] + 2*PBC[1]) ;
        f[12]= 1/24.0*(-2*f[0] - 4*f[1] + 20*f[11] - 4*f[13] + f[3] - 5*f[4] - 5*f[5] + f[6] -  4*f[7] - 4*f[9] + 2*PBC[1]);
        f[14]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] + 20*f[13] + f[3] - 5*f[4] + f[5] - 5*f[6] -  4*f[7] - 4*f[9] + 2*PBC[1]);

        Rho   [ib] = 0.0;
        Vel   [ib] = make_real3(0.0,0.0,0.0);
        for(size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ib] += F[iv + k];
            Vel[ib] = Vel[ib] + F[iv + k]*lbmaux[0].C[k];
        }
        Vel[ib] = (lbmaux[0].Cs/Rho[ib])*Vel[ib];
    }
}

void Setup (FLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    Left_BC<<<dom.Ndim(1)*dom.Ndim(2)/256+1,256>>>(dat.prhobc,dom.pIsSolid,dom.pF,dom.pVel,dom.pRho,dom.plbmaux);
    cudaDeviceSynchronize();
    Right_BC<<<dom.Ndim(1)*dom.Ndim(2)/256+1,256>>>(dat.prhobc,dom.pIsSolid,dom.pF,dom.pVel,dom.pRho,dom.plbmaux);
    cudaDeviceSynchronize();
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
    
    double vth = dat.dx/dat.dt*(dom.Ndim(2)-2)*(dom.Ndim(2)-2)/(12.0*(dom.Tau[0]-0.5))*(dat.rhobc[0]-dat.rhobc[1])/(dom.Ndim(0));

    dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << vave << Util::_8s << vth << std::endl;
    //std::cout << Util::_10_6 << dom.Time << Util::_8s << vave << Util::_8s << vth << std::endl;
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1;
    double nu    = 1.0/60.0;
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

    dat.rhobc.resize(2);
    dat.rhobc[1]  = 1.0;
    dat.rhobc[0]  = 1.03;
    dat.prhobc = thrust::raw_pointer_cast(dat.rhobc.data());
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

    Dom.Solve(1.0e4*dt,1.0e2*dt,Setup,Report,"single",true,Nproc);
    //Dom.Solve(1.0,80.0,NULL,NULL,"single",true,Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH
