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
// Stokes law

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    Vec3_t                acc;
    double                 nu;
    double                  R;
    double                 Tf;
    Array<Cell *>        xmin;
    Array<Cell *>        xmax;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c   = dom.Lat[0].Cells[i];
        c->BForcef = c->Rho*dat.acc;
    }

/*    
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0;i<dat.xmin.Size();i++)
    {
        Cell   * c = dat.xmin[i];
        //double * f = c->F;
        double  Cs = c->Cs;
        double vel;
        if (3.0*dom.Time<dat.Tf) vel = dat.acc(0)*3.0*dom.Time/dat.Tf;
        else                 vel = dat.acc(0);
        double rho = (c->F[0] + c->F[3] + c->F[4] + c->F[5] + c->F[6] + 2.0*(c->F[2] + c->F[8] + c->F[10] + c->F[12] + c->F[14]))/(1.0 - dat.acc(0)/Cs);
        c->F[ 1] = c->F[ 2] + 2.0/3.0 *rho*vel/Cs;
        c->F[ 7] = c->F[ 8] + 1.0/12.0*rho*vel/Cs - 1.0/4.0*( c->F[3] - c->F[4] + c->F[5] - c->F[6]);
        c->F[ 9] = c->F[10] + 1.0/12.0*rho*vel/Cs - 1.0/4.0*( c->F[3] - c->F[4] - c->F[5] + c->F[6]);
        c->F[11] = c->F[12] + 1.0/12.0*rho*vel/Cs - 1.0/4.0*(-c->F[3] + c->F[4] + c->F[5] - c->F[6]);
        c->F[13] = c->F[14] + 1.0/12.0*rho*vel/Cs - 1.0/4.0*(-c->F[3] + c->F[4] - c->F[5] + c->F[6]);
        //f[ 1] = -(2.0*(f[0]+f[3]+f[4]+f[5]+f[6]+2.0*f[8]+2.0*f[10]+2.0*f[12]+2.0*f[14])*dat.acc(0)/Cs+f[2]*(dat.acc(0)/Cs+3.0))/(3.0*(dat.acc(0)/Cs-1.0));
        //f[ 7] =  (2.0*f[ 8]*(5.0*dat.acc(0)/Cs-6.0)-(f[0]+2.0*f[2]+f[3]+f[4]+f[5]+f[6]+2.0*f[10]+2.0*f[12]+2.0*f[14])*dat.acc(0)/Cs)/(12.0*(dat.acc(0)/Cs-1.0));
        //f[ 9] =  (2.0*f[10]*(5.0*dat.acc(0)/Cs-6.0)-(f[0]+2.0*f[2]+f[3]+f[4]+f[5]+f[6]+2.0*f[ 8]+2.0*f[12]+2.0*f[14])*dat.acc(0)/Cs)/(12.0*(dat.acc(0)/Cs-1.0));
        //f[11] =  (2.0*f[12]*(5.0*dat.acc(0)/Cs-6.0)-(f[0]+2.0*f[2]+f[3]+f[4]+f[5]+f[6]+2.0*f[ 8]+2.0*f[10]+2.0*f[14])*dat.acc(0)/Cs)/(12.0*(dat.acc(0)/Cs-1.0));
        //f[13] =  (2.0*f[14]*(5.0*dat.acc(0)/Cs-6.0)-(f[0]+2.0*f[2]+f[3]+f[4]+f[5]+f[6]+2.0*f[ 8]+2.0*f[10]+2.0*f[12])*dat.acc(0)/Cs)/(12.0*(dat.acc(0)/Cs-1.0));
        c->Rho = c->VelDen(c->Vel);
    }

//#ifdef USE_OMP
    //#pragma omp parallel for schedule(static) num_threads(dom.Nproc)
//#endif
    //for (size_t i=0;i<dat.xmax.Size();i++)
    //{
        //Cell   * c = dat.xmax[i];
        //double * f = c->F;
        //double  Cs = c->Cs;
        //double vel;
        //if (3.0*dom.Time<dat.Tf) vel = dat.acc(0)*3.0*dom.Time/dat.Tf;
        //else                 vel = dat.acc(0);
        //double rho = (c->F[0] + c->F[3] + c->F[4] + c->F[5] + c->F[6] + 2.0*(c->F[1] + c->F[7] + c->F[9] + c->F[11] + c->F[13]))/(1.0 + dat.acc(0)/Cs);
        //c->F[ 2] = c->F[ 1] - 2.0/3.0 *rho*vel/Cs;
        //c->F[ 8] = c->F[ 7] - 1.0/12.0*rho*vel/Cs + 1.0/4.0*( c->F[3] - c->F[4] + c->F[5] - c->F[6]);
        //c->F[10] = c->F[ 9] - 1.0/12.0*rho*vel/Cs + 1.0/4.0*( c->F[3] - c->F[4] - c->F[5] + c->F[6]);
        //c->F[12] = c->F[11] - 1.0/12.0*rho*vel/Cs + 1.0/4.0*(-c->F[3] + c->F[4] + c->F[5] - c->F[6]);
        //c->F[14] = c->F[13] - 1.0/12.0*rho*vel/Cs + 1.0/4.0*(-c->F[3] + c->F[4] - c->F[5] + c->F[6]);
        //f[ 2] =  -(f[1]*(dat.acc(0)/Cs-3.0)+2.0*(f[0]+f[3]+f[4]+f[5]+f[6]+2.0*f[7]+2.0*f[9]+2.0*f[11]+2.0*f[13])*dat.acc(0)/Cs)/(3.0*(dat.acc(0)/Cs+1.0));
        //f[ 8] =   (2.0*f[ 7]*(5.0*dat.acc(0)/Cs+6.0)-(f[0]+2.0*f[1]+f[3]+f[4]+f[5]+f[6]+2.0*f[9]+2.0*f[11]+2.0*f[13])*dat.acc(0)/Cs)/(12.0*(dat.acc(0)/Cs+1.0));
        //f[10] =   (2.0*f[ 9]*(5.0*dat.acc(0)/Cs+6.0)-(f[0]+2.0*f[1]+f[3]+f[4]+f[5]+f[6]+2.0*f[7]+2.0*f[11]+2.0*f[13])*dat.acc(0)/Cs)/(12.0*(dat.acc(0)/Cs+1.0));
        //f[12] =   (2.0*f[11]*(5.0*dat.acc(0)/Cs+6.0)-(f[0]+2.0*f[1]+f[3]+f[4]+f[5]+f[6]+2.0*f[7]+2.0*f[9 ]+2.0*f[13])*dat.acc(0)/Cs)/(12.0*(dat.acc(0)/Cs+1.0));
        //f[14] =   (2.0*f[13]*(5.0*dat.acc(0)/Cs+6.0)-(f[0]+2.0*f[1]+f[3]+f[4]+f[5]+f[6]+2.0*f[7]+2.0*f[9 ]+2.0*f[11])*dat.acc(0)/Cs)/(12.0*(dat.acc(0)/Cs+1.0));
        //c->Rho = c->VelDen(c->Vel);
    //}

    //for (size_t i=0;i<dat.xmin.Size();i++)
    //{
        //Cell * c = dat.xmin[i];
        //if(c->IsSolid) continue;
        //c->F[1] = 1.0/3.0 *(-2*c->F[0] - 4*c->F[10] - 4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->F[7] = 1.0/24.0*(-2*c->F[0] - 4*c->F[10] - 4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
        //c->F[9] = 1.0/24.0*(-2*c->F[0] + 20*c->F[10] - 4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->F[11]= 1.0/24.0*(-2*c->F[0] - 4*c->F[10] + 20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->F[13]= 1.0/24.0*(-2*c->F[0] - 4*c->F[10] - 4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->Rho = c->VelDen(c->Vel);
    //}
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
#endif
    for (size_t i=0;i<dat.xmax.Size();i++)
    {
        Cell * c = dat.xmax[i];
        if(c->IsSolid) continue;
        c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-1.0));
        c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*1.0);
        c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*1.0) ;
        c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*1.0);
        c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*1.0);
        c->Rho = c->VelDen(c->Vel);
    }
*/
}

void Report(LBM::Domain & dom, void * UD)
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
        double CD  = 2.0*dom.Particles[0]->F(0)/(Flux(0)*Flux(0)/M*M_PI*dat.R*dat.R);
        double Re  = 2.0*Flux(0)*dat.R/(M*dat.nu);
        double CDt = 24.0/Re + 6.0/(1.0+sqrt(Re)) + 0.4;
        if (Flux(0)<1.0e-12) 
        {
            Flux = OrthoSys::O;
            CD = 0.0;
            Re = 0.0;
        }

        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dom.Particles[0]->F(0) << Util::_8s << dom.Particles[0]->F(1) << Util::_8s << dom.Particles[0]->F(2) << Util::_8s << dom.Particles[0]->T(0) << Util::_8s << dom.Particles[0]->T(1) << Util::_8s << dom.Particles[0]->T(2) << Util::_8s << Vx << Util::_8s << Flux(0) << Util::_8s << norm(Flux) << Util::_8s << M << Util::_8s << Re << Util::_8s << CD << Util::_8s << CDt << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

int main(int argc, char **argv) try
{
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    String ptype;
    bool   Render = true;
    size_t nx = 100;
    size_t ny = 50;
    size_t nz = 50;
    double nu = 0.01;
    double dx = 1.0;
    double dt = 1.0;
    double Dp = 0.1;
    double R  = 10.0;
    double w  = 0.001;
    double ang= 0.0;
    double Tf = 40000.0;
    {
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> Render;       infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> nu;           infile.ignore(200,'\n');
        infile >> dx;           infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> Dp;           infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> w;            infile.ignore(200,'\n');
        infile >> ang;          infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
    }
    
    

    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    Dom.Step     = 1;
    Dom.Sc       = 0.0;
    Dom.Alpha    = dx;
    UserData dat;
    Dom.UserData = &dat;
    dat.acc      = Vec3_t(Dp,0.0,0.0);
    dat.R        = R;
    dat.nu       = nu;
    dat.Tf       = Tf;

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(1.0,OrthoSys::O);
        if (i==0   ) dat.xmin.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,k)));
        if (i==nx-1) dat.xmax.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,k)));
    }

    if       (ptype=="sphere")  Dom.AddSphere(-1,Vec3_t(0.5*nx*dx,0.5*ny*dx,0.5*nz*dx),R,3.0);
    else if  (ptype=="tetra" )
    { 
        double e = pow(sqrt(2)*M_PI,1.0/3.0)*2*R;
        Dom.AddTetra(-1,Vec3_t(0.5*nx*dx,0.5*ny*dx,0.5*nz*dx),0.05*e,e,3.0,M_PI/4.0,&OrthoSys::e2);
        Quaternion_t q;
        NormalizeRotation(35.26*M_PI/180.0,OrthoSys::e1,q);
        Dom.Particles[0]->Rotate(q,Dom.Particles[0]->x);
        NormalizeRotation(ang*M_PI/180.0,OrthoSys::e2,q);
        Dom.Particles[0]->Rotate(q,Dom.Particles[0]->x);
    }
    else if  (ptype=="cube"  )
    {
        double e = pow(M_PI/6.0,1.0/3.0)*2*R;
        Dom.AddCube(-1,Vec3_t(0.5*nx*dx,0.5*ny*dx,0.5*nz*dx),0.05*e,e,3.0,ang*M_PI/180.0,&OrthoSys::e2);
    }
    else
    {
        Dom.AddFromJson(-1, ptype.CStr(), 0.05*R,3.0,R);
        Vec3_t t(0.5*nx*dx,0.5*ny*dx,0.5*nz*dx);
        Dom.Particles[0]->Position(t);
    }

    Dom.Particles[0]->FixVeloc();
    Dom.Particles[0]->w = Vec3_t(0.0,0.0,w);
    Dom.Alpha = 2.0*dx;

    //Solving
    Dom.Solve(Tf,0.01*Tf,Setup,Report,filekey.CStr(),Render,Nproc);
}
MECHSYS_CATCH

