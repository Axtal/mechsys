/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2021 Sergio Galindo                                    *
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
// Testing real bouyancy


// MechSys
#include <mechsys/lbmdem/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    thrust::device_vector<real3> g; 
    real3 *    pg;
    thrust::device_vector<size_t> idxp;
    size_t *   pidxp;
    thrust::device_vector<size_t> imax;
    size_t *   pimax;
    bool       change = false;
    double     gr;
    double     rhof;
    double     Tf;
};

__global__ void ApplyGravity(real3 const * pg, real3 * BForce, real * Rho, FLBM::lbm_aux * lbmaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=lbmaux[0].Nl*lbmaux[0].Ncells) return;
    BForce[ic] = Rho[ic]*pg[0];
}

__global__ void FreePar(size_t const * idx, size_t const * imax, DEM::ParticleCU * Par)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=imax[0]) return;
    Par[idx[ic]].vxf = false;
    Par[idx[ic]].vyf = false;
    Par[idx[ic]].vzf = false;
    Par[idx[ic]].wxf = false;
    Par[idx[ic]].wyf = false;
    Par[idx[ic]].wzf = false;
}

void Setup (LBMDEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    ApplyGravity<<<dom.LBMDOM.Ncells*dom.LBMDOM.Nl/dom.Nthread+1,dom.Nthread>>>(dat.pg,dom.LBMDOM.pBForce,dom.LBMDOM.pRho,dom.LBMDOM.plbmaux);

    //if (dom.Time>0.5*dat.Tf&&!dat.change)
    //{
        //FreePar<<<dat.idxp.size()/dom.Nthread+1,dom.Nthread>>>(dat.pidxp,dat.pimax,dom.DEMDOM.pParticlesCU);
        //dat.change = true;
    //}
}

void Report (LBMDEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_force.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "Fz" << Util::_8s << "Ffz" << Util::_8s << "Ft" << Util::_8s << "Z" << std::endl;
    }
    if (!dom.Finished)
    {
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dom.DEMDOM.Particles[0]->F(2) << Util::_8s << dom.DEMDOM.Particles[0]->Flbm(2) << Util::_8s <<
            dom.DEMDOM.Particles[0]->Props.V*dat.rhof*dat.gr << Util::_8s << dom.DEMDOM.Particles[0]->x(2)  << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

int main(int argc, char **argv) try
{
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>=2) Nproc = atoi(argv[1]);

        
    double nu  = 1.5;
    size_t nx  = 100;
    size_t ny  = 100;
    size_t nz  = 100;
    double dx  = 0.4;
    double dt  = 1.6e-2;
    double R   = 7.2;
    double Cs2 = dx*dx/(3.0*dt*dt);
    LBMDEM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;
    dat.gr       = 1.0e-5*Cs2/(nz*dx);
    dat.rhof     = 2.0;
    dat.g.resize(1);
    dat.g[0]     = make_real3(0.0,0.0,-dat.gr);
    dat.pg       = thrust::raw_pointer_cast(dat.g.data());
    
    dom.DEMDOM.AddSphere(-1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.5*dx*nz),R,1.0);
    //double e = pow(M_PI/6.0,1.0/3.0)*2*R;
    //dom.DEMDOM.AddCube(-1,Vec3_t(0.95*dx*nx,0.95*dx*ny,0.95*dx*nz),0.05*e,e,1.0);
    //dom.DEMDOM.GetParticle(-1)->Ff = dom.DEMDOM.GetParticle(-1)->Props.m*Vec3_t(0.0,0.0,-dat.gr);
    //dom.DEMDOM.AddCube(-2,Vec3_t(0.1*dx*nx,0.1*dx*ny,0.1*dx*nz),0.05*e,e,1.0);
    dom.DEMDOM.GetParticle(-1)->FixVeloc();

    double rho0 = dat.rhof*dat.gr*(nz*dx)/(Cs2*(exp(nz*dx*dat.gr/Cs2)-1));
    //std::cout << dat.gr << " " << dx*nz << " " << Cs2 << std::endl;
    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        double rho = rho0*exp((nz*dx - dx*(iz+1))*dat.gr/Cs2);
        if (((ix*dx-0.5*dx*nx)*(ix*dx-0.5*dx*nx)+(iy*dx-0.5*dx*ny)*(iy*dx-0.5*dx*ny)+(iz*dx-0.5*dx*nz)*(iz*dx-0.5*dx*nz))<R*R)
        {
            dom.LBMDOM.Initialize(0,idx,dat.rhof,OrthoSys::O);
        }
        else
        {
            dom.LBMDOM.Initialize(0,idx,rho,v);
        }
        //dom.LBMDOM.Initialize(0,idx,dat.rhof,OrthoSys::O);
        //if (ix==0||iy==0||iz==0) dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
        if (iz==0) dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
    }   

    dom.Alpha = 2.0*dx;
    //dom.PeriodicX= true;
    //dom.PeriodicY= true;
    //dom.PeriodicZ= true;


    thrust::host_vector<size_t> hidxp;
    hidxp.push_back(0);
    dat.idxp = hidxp;

    dat.imax.push_back(dat.idxp.size());

    dat.pidxp = thrust::raw_pointer_cast(dat.idxp.data());
    dat.pimax = thrust::raw_pointer_cast(dat.imax.data());
    
    double Tf = 1.0e3;
    dat.Tf = Tf;
    dom.Solve(Tf,Tf/200,Setup,Report,"tlbmdem_cu_03",true,Nproc);
}
MECHSYS_CATCH

