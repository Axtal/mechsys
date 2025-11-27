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
// Periodic boundary conditions


// MechSys
#include <mechsys/lbmdem/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    Vec3_t                acc;
    double                 nu;
    double                  R;
};

void Setup (LBMDEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //#pragma omp parallel for schedule(static) num_threads(dom.Nproc)
	//for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
	//for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
	//for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
    //{
        //dom.LBMDOM.BForce[0][ix][iy][iz] = dom.LBMDOM.Rho[0][ix][iy][iz]*dat.acc;
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
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "Fx" << Util::_8s << "Vx" << Util::_8s << "Fa \n";
    }
    if (!dom.Finished)
    {
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << norm(dom.DEMDOM.Particles[0]->Flbm) << Util::_8s << norm(dom.DEMDOM.Particles[0]->v) << Util::_8s << 6.0*M_PI*dat.nu*dat.R*norm(dom.DEMDOM.Particles[0]->v) << std::endl;
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
    size_t nthread = 256;
    if (argc>=3) nthread = atoi(argv[2]);

        
    double nu = 1.6;
    size_t nx = 200;
    size_t ny = 200;
    size_t nz = 200;
    double dx = 0.4;
    double dt = 1.6e-2;
    double R = 1.8;
    LBMDEM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;
    dat.R  = R;
    dat.nu = nu;
    //dat.acc = Vec3_t(1.0e-2,1.0e-2,1.0e-2);
    dat.acc = Vec3_t(1.0e-2,0.0,0.0);
    dom.DEMDOM.AddSphere(-1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.5*dx*nz),R,1.0);
    //dom.DEMDOM.AddSphere(-1,Vec3_t(0.95*dx*nx,0.5*dx*ny,0.5*dx*nz),R,1.0);
    //dom.DEMDOM.AddSphere(-1,Vec3_t(0.95*dx*nx,0.95*dx*ny,0.95*dx*nz),R,1.0);
    //double e = pow(M_PI/6.0,1.0/3.0)*2*R;
    //dom.DEMDOM.AddCube(-1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.5*dx*nz),0.05*e,e,1.0);
    //dom.DEMDOM.AddRecBox(-1,Vec3_t(0.5*dx*nx,0.5*dx*ny,0.5*dx*nz),Vec3_t(e,2.0*e,3.0*e),0.05*e,1.0,0.0,&OrthoSys::e0);
    dom.DEMDOM.GetParticle(-1)->Ff = dom.DEMDOM.GetParticle(-1)->Props.m*dat.acc;
    //dom.DEMDOM.AddRecBox(-2,Vec3_t(0.8*dx*nx,0.5*dx*ny,0.5*dx*nz),Vec3_t(e,2.0*e,3.0*e),0.05*e,1.0);
    //dom.DEMDOM.GetParticle(-2)->FixVeloc();
    //dom.DEMDOM.GetParticle(-2)->FixFree = true;
    //dom.DEMDOM.AddCube(-2,Vec3_t(0.1*dx*nx,0.1*dx*ny,0.1*dx*nz),0.05*e,e,1.0);


    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,1.0/*rho*/,v);
        if ((iz==0)||(iz==nz-1)||(iy==0)||(iy==ny-1)) dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
    }   

    dom.Alpha = 2.0*dx;
    dom.PeriodicX= true;
    dom.PeriodicY= true;
    dom.PeriodicZ= true;
    
    double Tf = 2.0e3;
    dom.Solve(Tf,Tf/200,Setup,Report,"tlbmdem_cu_02",true,Nproc);
}
MECHSYS_CATCH

