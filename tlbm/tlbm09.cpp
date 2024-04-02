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
// Compression test

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    Array<Cell *>         Top;
    double                 Tf;
    double              rhobc;
    double                  g;
    double                 dx;
    double                 Cs;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0;i<dat.Top.Size();i++)
    {
        Cell * c = dat.Top[i];
        if (c->IsSolid) continue;
        /*D3Q15*/
        //c->F[6]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]- 2*c->F[2]- 2*c->F[3]- 2*c->F[4]- c->F[5]- 4*c->F[7]+ 2*dat.rhobc);
        //c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]+ c->F[2]- 5*c->F[3]+ c->F[4]- 4*c->F[5]+ 20*c->F[7]+ 2*dat.rhobc);
        //c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[14]- 5*c->F[2]+ c->F[3]- 5*c->F[4]- 4*c->F[5]- 4*c->F[7]+ 2*dat.rhobc);
        //c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[14]+ c->F[2]+ c->F[3]- 5*c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*dat.rhobc);
        //c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[14]- 5*c->F[2]- 5*c->F[3]+ c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*dat.rhobc);
        /*D3Q19*/
        c->F[6 ] = 1.0/3.0*(-c->F[0]-c->F[1]-c->F[2]-c->F[3]-c->F[4]+c->F[5]-c->F[7]-c->F[8]-c->F[9]-c->F[10]-2.0*c->F[11]-2.0*c->F[14]-2.0*c->F[15]-2.0*c->F[18]+dat.rhobc);
        c->F[12] = 1.0/6.0*(-c->F[0]-c->F[1]-c->F[2]-c->F[3]-c->F[4]-2.0*c->F[5]-c->F[7]-c->F[8]-c->F[9]-c->F[10]+4.0*c->F[11]-2.0*c->F[14]-2.0*c->F[15]-2.0*c->F[18]+dat.rhobc);
        c->F[13] = 1.0/6.0*(-c->F[0]-c->F[1]-c->F[2]-c->F[3]-c->F[4]-2.0*c->F[5]-c->F[7]-c->F[8]-c->F[9]-c->F[10]-2.0*c->F[11]+4.0*c->F[14]-2.0*c->F[15]-2.0*c->F[18]+dat.rhobc);
        c->F[16] = 1.0/6.0*(-c->F[0]-c->F[1]-c->F[2]-c->F[3]-c->F[4]-2.0*c->F[5]-c->F[7]-c->F[8]-c->F[9]-c->F[10]-2.0*c->F[11]-2.0*c->F[14]+4.0*c->F[15]-2.0*c->F[18]+dat.rhobc);
        c->F[17] = 1.0/6.0*(-c->F[0]-c->F[1]-c->F[2]-c->F[3]-c->F[4]-2.0*c->F[5]-c->F[7]-c->F[8]-c->F[9]-c->F[10]-2.0*c->F[11]-2.0*c->F[14]-2.0*c->F[15]+4.0*c->F[18]+dat.rhobc);
        c->Rho = c->VelDen(c->Vel);
    }

}

void Report(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_pressure.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "P_up" << Util::_8s << "P_down" << Util::_8s << "P_th" << Util::_8s << "Z" << Util::_8s << "VZ" << std::endl;
    }
    if (!dom.Finished) 
    {
        size_t Ndown = dom.Particles[0]->x(2)/dat.dx - dom.Lat[0].Ndim(0)/2 - 2;
        double P_down = 0.0;
        size_t ncells = 0;
        for (size_t i=0;i<dom.Lat[0].Ndim(0);i++)
        for (size_t j=0;j<dom.Lat[0].Ndim(1);j++)
        for (size_t k=0;k<Ndown;k++)
        {
            P_down += dat.Cs*dat.Cs*(dom.Lat[0].GetCell(iVec3_t(i,j,k))->Rho-1.0);
            ncells++;
        }
        P_down /= ncells;

        size_t Nup  = dom.Particles[0]->x(2)/dat.dx + dom.Lat[0].Ndim(0)/2 + 2;
        double P_up = 0.0;
        ncells      = 0;
        for (size_t i=0;i<dom.Lat[0].Ndim(0);i++)
        for (size_t j=0;j<dom.Lat[0].Ndim(1);j++)
        for (size_t k=Nup+1;k<dom.Lat[0].Ndim(2);k++)
        {
            P_up += dat.Cs*dat.Cs*(dom.Lat[0].GetCell(iVec3_t(i,j,k))->Rho-1.0);
            ncells++;
        }
        P_up /= ncells;

        double P_th = dom.Particles[0]->Props.m*dat.g/(dom.Lat[0].Ndim(0)*dat.dx*dom.Lat[0].Ndim(1)*dat.dx);
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << P_up << Util::_8s << P_down << Util::_8s << P_th << Util::_8s << dom.Particles[0]->x(2) << Util::_8s << dom.Particles[0]->v(2) << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    if (argc>=2) Nproc=atoi(argv[1]);

    bool   Render = true;
    size_t nx = 50;
    size_t ny = 50;
    size_t nz = 250;
    double nu = 0.16;
    double dx = 1.0;
    double dt = 1.0;
    double R  = 10.0;
    double Tf = 100000.0;
    double g  = 0.333333333e-6;
    
    if (argc>=3) g =atof(argv[2]);
    if (argc>=4) nu=atof(argv[3]);
    

    LBM::Domain Dom(D3Q19, nu, iVec3_t(nx,ny,nz), dx, dt);
    Dom.Step     = 1;
    //Dom.Sc       = 0.0;
    Dom.Alpha    = dx;
    UserData dat;
    Dom.UserData = &dat;
    dat.Tf       = Tf;
    dat.rhobc    = 1.0;
    dat.g        = g;
    dat.dx       = dx;
    dat.Cs       = dx/(dt*sqrt(3.0));

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(1.0,OrthoSys::O);
        if (k==nz-1) dat.Top.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,k)));
        if (k==0)    Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
    }

    Dom.AddCube(-1,Vec3_t(0.5*nx*dx,0.5*ny*dx,0.5*nz*dx),dx,nx*dx,2000.0,0.0,&OrthoSys::e2);
    Dom.Particles[0]->FixVeloc();
    Dom.Particles[0]->vzf = false;
    Dom.Particles[0]->Ff  = Dom.Particles[0]->Props.m*Vec3_t(0.0,0.0,-g);
    Dom.Particles[0]->Props.Gv = 0.0001;
    
    //Solving
    Dom.Solve(Tf,0.01*Tf,Setup,Report,"tlbm09",Render,Nproc);
}
MECHSYS_CATCH

