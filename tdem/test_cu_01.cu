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

/////////////////////// Test 01 sliding bodies

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

int main(int argc, char **argv) try
{
    DEM::Domain dom;
    double dt = 1.0e-5;     // time step
    double g = 9.8;         // gravity acceleration
    double velocity =  10.0; // velocity of the sliding block.
    double inipos = -5.0;   // Initial position in the x axis
    double mu = 0.1;        // friction coefficient
    //dom.Xmax = 25.0;
    //dom.Xmin =-25.0;
    //dom.Alpha = 100.0;
    //dom.Beta  = 4.0;
    //dom.AddCube(/*Tag*/-1,/*initial position*/Vec3_t(inipos,0.0,0.7), /*spheroradius*/0.1,/*side length*/ 1.0,  /*rho*/1.0, /*orientation*/M_PI/2.0,/*axis of orientation*/&OrthoSys::e0);
    dom.AddCube(/*Tag*/-1,/*initial position*/Vec3_t( inipos,0.0,0.7), /*spheroradius*/0.1,/*side length*/ 1.0,  /*rho*/1.0, /*orientation*/M_PI/2.0,/*axis of orientation*/&OrthoSys::e0);
    dom.AddPlane(/*Tag*/-2,OrthoSys::O,0.1,100.0,100.0,1.0);
    //dom.AddSphere(-3,Vec3_t(0.0,0.0,0.6),0.5,1.0);
    //dom.AddSphere(-4,Vec3_t(1.0,0.0,0.6),0.5,1.0);
    //dom.AddCube(/*Tag*/-5,/*initial position*/Vec3_t(-24.4,0.5,0.5), /*spheroradius*/0.1,/*side length*/ 1.0,  /*rho*/1.0, /*orientation*/M_PI/2.0,/*axis of orientation*/&OrthoSys::e0);
    DEM::Particle *p = dom.GetParticle(-2);
    p->FixVeloc();
    //DEM::Particle * p = dom.GetParticle(-1);
    p = dom.GetParticle(-1);
    //p->FixVeloc();
    p->v  = velocity , 0.0, 0.0;
    dom.Initialize(dt);
    //p->Ff = 0.0,0.0,-p->Props.m*g;
    for (size_t ip=0;ip<dom.Particles.Size();ip++)
    {
        dom.Particles[ip]->Ff = Vec3_t(0.0,0.0,-dom.Particles[ip]->Props.m*9.8);
    }


    Dict B;
    B.Set(-1,"Gn Gt Kn Kt Mu",-0.5,0.0,1.0e7,5.0e6,mu);
    B.Set(-2,"Gn Gt Kn Kt Mu",-0.5,0.0,1.0e7,5.0e6,mu);
    B.Set(-3,"Gn Gt Kn Kt Mu",-0.5,0.0,1.0e7,5.0e6,mu);
    B.Set(-4,"Gn Gt Kn Kt Mu",-0.5,0.0,1.0e7,5.0e6,mu);
    B.Set(-5,"Gn Gt Kn Kt Mu",-0.5,0.0,1.0e7,5.0e6,mu);
    dom.SetProps(B);
    //std::pair<int,int> pr (-1/*tag1*/,-2/*tag2*/);
    //dom.FricCoeff[pr] = mu;
    //dom.AddSphere(-3,Vec3_t(0.0,0.0,0.6),0.5,1.0);
    //dom.AddSphere(-4,Vec3_t(1.0,0.0,0.6),0.5,1.0);
    //dom.AddCube(/*Tag*/-5,/*initial position*/Vec3_t(inipos-1.2,0.0,0.7), /*spheroradius*/0.1,/*side length*/ 1.0,  /*rho*/1.0, /*orientation*/M_PI/4.0,/*axis of orientation*/&OrthoSys::e0);
    //dom.AddPlane(/*Tag*/-6,Vec3_t(0.0,0.0,0.2),0.1,100.0,100.0,1.0,3.0*M_PI/4.0,&OrthoSys::e2);

    dom.Solve(/*final time*/10.0,/*time step*/dt,/*Output step*/0.1,NULL,NULL,/*file key*/"test_cu_01",/*Render video?*/2);

    double finpos = p->x(0);
    std::cout << "  Observed L = " << finpos - inipos<< " Theoretical L = " << velocity*velocity/(2*mu*g) << "\n";
    //dom.UpLoadDevice(1);
/*
    thrust::host_vector<DEM::ParticleCU> hParticlesCU = dom.bParticlesCU;
    thrust::host_vector<real3>           hVertsCU     = dom.bVertsCU;
    thrust::host_vector<size_t>          hEdgesCU     = dom.bEdgesCU;
    thrust::host_vector<size_t>          hFacidCU     = dom.bFacidCU;
    thrust::host_vector<size_t>          hFacesCU     = dom.bFacesCU;
    DEM::ParticleCU                    * pParticlesCU = thrust::raw_pointer_cast(hParticlesCU.data());
    real3                              * pVertsCU     = thrust::raw_pointer_cast(hVertsCU    .data());
    size_t                             * pEdgesCU     = thrust::raw_pointer_cast(hEdgesCU    .data());
    size_t                             * pFacidCU     = thrust::raw_pointer_cast(hFacidCU    .data());
    size_t                             * pFacesCU     = thrust::raw_pointer_cast(hFacesCU    .data());

    for (size_t ip=0;ip<hParticlesCU.size();ip++)
    {
        std::cout << "\nParticle " << ip << "\n" << std::endl;
        for (size_t iv=hParticlesCU[ip].Nvi;iv<hParticlesCU[ip].Nvf;iv++)
        {
            std::cout << "Vertex " << iv << " " << hVertsCU[iv].x << " " << hVertsCU[iv].y << " " << hVertsCU[iv].z << std::endl; 
        }
        for (size_t ie=hParticlesCU[ip].Nei;ie<hParticlesCU[ip].Nef;ie++)
        {
            std::cout << "Edge   "   << ie << " " << hEdgesCU[2*ie  ] << " " << hEdgesCU[2*ie+1] << std::endl; 
        }
        for (size_t ic=hParticlesCU[ip].Nfi;ic<hParticlesCU[ip].Nff;ic++)
        {
            std::cout << "Face   " << ic;
            for (size_t ig=hFacesCU[2*ic];ig<hFacesCU[2*ic+1];ig++)
            {
                std::cout << " " << hFacidCU[ig];
            }
            std::cout << std::endl;
        }
    }

    {    
    thrust::host_vector<DEM::InteractonCU>    hInteractons      = dom.bInteractons     ;
    thrust::host_vector<DEM::DynInteractonCU> hDynInteractonsVV = dom.bDynInteractonsVV;
    thrust::host_vector<DEM::DynInteractonCU> hDynInteractonsEE = dom.bDynInteractonsEE;
    thrust::host_vector<DEM::DynInteractonCU> hDynInteractonsVF = dom.bDynInteractonsVF;
    thrust::host_vector<DEM::DynInteractonCU> hDynInteractonsFV = dom.bDynInteractonsFV;
    //
    //for (size_t ivv=0;ivv<hInteractonsVV.size();ivv++)
    //{
        //std::cout << "\nInteracton VV " << ivv << "\n" << std::endl;
        //std::cout << "Particles " << hInteractonsVV[ivv].I1 << " " << hInteractonsVV[ivv].I2 << std::endl;
    //}
//
    for (size_t iee=0;iee<hDynInteractonsEE.size();iee++)
    {
        std::cout << "\nInteracton EE " << iee << "\n" << std::endl;
        //std::cout << "Particles " << hInteractonsEE[iee].I1 << " " << hInteractonsEE[iee].I2 << std::endl;
        std::cout << "Features  " << hDynInteractonsEE[iee].IF1  << " " << hDynInteractonsEE[iee].IF2  << std::endl;
        std::cout << "Force     " << hDynInteractonsEE[iee].Fn.x << " " << hDynInteractonsEE[iee].Ft.x << std::endl;
    }
//
    //for (size_t ivf=0;ivf<hDynInteractonsVF.size();ivf++)
    //{
        //std::cout << "\nInteracton VF " << ivf << "\n" << std::endl;
        //std::cout << "Particles " << hInteractonsVF[ivf].I1 << " " << hInteractonsVF[ivf].I2 << std::endl;
        //std::cout << "Features  " << hInteractonsVF[ivf].IF1<< " " << hInteractonsVF[ivf].IF2<< std::endl;
        //std::cout << "Force     " << norm(hDynInteractonsVF[ivf].F)  << " " << norm(hDynInteractonsVF[ivf].Ft) << std::endl;
    //}
//
    //for (size_t ifv=0;ifv<hDynInteractonsFV.size();ifv++)
    //{
        //std::cout << "\nInteracton FV " << ifv << "\n" << std::endl;
        //std::cout << "Particles " << hInteractonsFV[ifv].I1 << " " << hInteractonsFV[ifv].I2 << std::endl;
        //std::cout << "Features  " << hInteractonsFV[ifv].IF1<< " " << hInteractonsFV[ifv].IF2<< std::endl;
        //std::cout << "Force     " << norm(hDynInteractonsFV[ifv].F)  << " " << norm(hDynInteractonsFV[ifv].Ft) << std::endl;
    //}
    }
    */
    //dom.DnLoadDevice(1);
    
    //dom.WriteXDMF("test_cu_01");

    //size_t e0 = 1;
    //size_t e1 = 15;
    //real3 xi,xf,s;
    //DEM::DistanceEE(pEdgesCU,pVertsCU,e0,e1,xi,xf,s,dom.demaux.Per);
//
    //std::cout << "\nDistance" <<std::endl;
    //std::cout << "Edges " << e0   << " " << e1                  << std::endl;
    //std::cout << "S     " << s .x << " " << s .y << " " << s .z << std::endl;
    //std::cout << "Xi    " << xi.x << " " << xi.y << " " << xi.z << std::endl;
    //std::cout << "Xf    " << xf.x << " " << xf.y << " " << xf.z << std::endl;

    //real3 v  = make_real3( -0.0, -55.0,-0.5);
    //real3 xf,s;
    //DEM::DistanceFV(pFacesCU,pFacidCU,pVertsCU,0,v,xf,s);
    //std::cout << "\nDistance" <<std::endl;
    //std::cout << "S     " << s .x << " " << s .y << " " << s .z << std::endl;
    //std::cout << "Xf    " << xf.x << " " << xf.y << " " << xf.z << std::endl;
   
    //real3 v = make_real3(1.0,0.0,0.0), axis=make_real3(0.0,1.0,1.0), p;
    //real4 q;
    //NormalizeRotation(M_PI/3.0,axis,q);
    //Rotation(v,q,p);
    //std::cout << "P    " << p.x << " " << p.y << " " << p.z << std::endl;
    return 0;
}
MECHSYS_CATCH

