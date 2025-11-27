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

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

struct UserData
{
    bool               StrainCtrl;   ///< Is a failuretest ?
    bool               RenderVideo;  ///< RenderVideo ?
    size_t             InitialIndex; ///< The initial index marking the bounding box
    double             Thf;          ///< Angle in the p=cte plane
    double             Alp;          ///< Angle in the p q plane
    double             dt;           ///< Time step
    double             tspan;        ///< Time span for the different stages
    Vec3_t             Sig;          ///< Current stress state
    Vec3_t             Sig0;         ///< Initial stress state
    Vec3_t             DSig;         ///< Total stress increment to be applied by Solve => after
    bVec3_t            pSig;         ///< Prescribed stress ?
    Vec3_t             L0;           ///< Initial length of the packing
    std::ofstream      oss_ss;       ///< file for stress strain data

    //Constructor
    UserData() {Sig = 0.0,0.0,0.0;}     
};

void SetTxTest (Vec3_t const & Sigf, bVec3_t const & pEps, Vec3_t const & dEpsdt, double theta, double alpha, bool TheStrainCtrl, UserData & UD, DEM::Domain const & D)
{
    // info
    std::cout << "[1;33m\n--- Setting up Triaxial Test -------------------------------------[0m\n";
    double start = std::clock();

    // Store setting up data
    UD.StrainCtrl = TheStrainCtrl;
    UD.Thf = theta;
    UD.Alp = alpha;
    if (TheStrainCtrl) UD.Sig0 = UD.Sig;

    // initialize particles
    for (size_t i=0; i<D.Particles.Size(); i++) D.Particles[i]->Initialize(i);

    // total stress increment
    UD.DSig = Sigf - UD.Sig;

    // assume all strains prescribed by default
    UD.pSig = false, false, false;

    // Eps(0) prescribed ?
    Vec3_t veloc, force;
    if (pEps(0))
    {
        double height = (D.Particles[UD.InitialIndex]->x(0)-D.Particles[UD.InitialIndex+1]->x(0));
        veloc = 0.5*dEpsdt(0)*height, 0.0, 0.0;
        D.Particles[UD.InitialIndex  ]->Ff = 0.0,0.0,0.0;
        D.Particles[UD.InitialIndex  ]->FixVeloc();
        D.Particles[UD.InitialIndex  ]->v  =  veloc;
        D.Particles[UD.InitialIndex+1]->Ff = 0.0,0.0,0.0;
        D.Particles[UD.InitialIndex+1]->FixVeloc();
        D.Particles[UD.InitialIndex+1]->v  = -veloc;
    }
    else // UD.Sig(0) prescribed
    {
        double area = (D.Particles[UD.InitialIndex+2]->x(1)-D.Particles[UD.InitialIndex+3]->x(1))*(D.Particles[UD.InitialIndex+4]->x(2)-D.Particles[UD.InitialIndex+5]->x(2));
        force = UD.Sig(0)*area, 0.0, 0.0;
        D.Particles[UD.InitialIndex  ]->Ff =  force;
        D.Particles[UD.InitialIndex  ]->FixVeloc();
        D.Particles[UD.InitialIndex  ]->vxf = false;
        D.Particles[UD.InitialIndex+1]->Ff = -force;
        D.Particles[UD.InitialIndex+1]->FixVeloc();
        D.Particles[UD.InitialIndex+1]->vxf = false;
        UD.pSig(0) = true;
    }

    // Eps(1) prescribed ?
    if (pEps(1))
    {
        double height = (D.Particles[UD.InitialIndex+2]->x(1)-D.Particles[UD.InitialIndex+3]->x(1));
        veloc = 0.0, 0.5*dEpsdt(1)*height, 0.0;
        D.Particles[UD.InitialIndex+2]->Ff = 0.0,0.0,0.0;
        D.Particles[UD.InitialIndex+2]->FixVeloc();
        D.Particles[UD.InitialIndex+2]->v  =  veloc;
        D.Particles[UD.InitialIndex+3]->Ff = 0.0,0.0,0.0;
        D.Particles[UD.InitialIndex+3]->FixVeloc();
        D.Particles[UD.InitialIndex+3]->v  = -veloc;
    }
    else // UD.Sig(1) presscribed
    {
        double area = (D.Particles[UD.InitialIndex]->x(0)-D.Particles[UD.InitialIndex+1]->x(0))*(D.Particles[UD.InitialIndex+4]->x(2)-D.Particles[UD.InitialIndex+5]->x(2));
        force = 0.0, UD.Sig(1)*area, 0.0;
        D.Particles[UD.InitialIndex+2]->Ff =  force;
        D.Particles[UD.InitialIndex+2]->FixVeloc();
        D.Particles[UD.InitialIndex+2]->vyf = false;
        D.Particles[UD.InitialIndex+3]->Ff = -force;
        D.Particles[UD.InitialIndex+3]->FixVeloc();
        D.Particles[UD.InitialIndex+3]->vyf = false;
        UD.pSig(1) = true;
    }

    // Eps(2) prescribed ?
    if (pEps(2))
    {
        double height = (D.Particles[UD.InitialIndex+4]->x(2)-D.Particles[UD.InitialIndex+5]->x(2));
        veloc = 0.0, 0.0, 0.5*dEpsdt(2)*height;
        D.Particles[UD.InitialIndex+4]->Ff = 0.0,0.0,0.0;
        D.Particles[UD.InitialIndex+4]->FixVeloc();
        D.Particles[UD.InitialIndex+4]->v  =  veloc;
        D.Particles[UD.InitialIndex+5]->Ff = 0.0,0.0,0.0;
        D.Particles[UD.InitialIndex+5]->FixVeloc();
        D.Particles[UD.InitialIndex+5]->v  = -veloc;
    }
    else // UD.Sig(2) presscribed
    {
        double area = (D.Particles[UD.InitialIndex]->x(0)-D.Particles[UD.InitialIndex+1]->x(0))*(D.Particles[UD.InitialIndex+2]->x(1)-D.Particles[UD.InitialIndex+3]->x(1));
        force = 0.0, 0.0, UD.Sig(2)*area;
        D.Particles[UD.InitialIndex+4]->Ff =  force;
        D.Particles[UD.InitialIndex+4]->FixVeloc();
        D.Particles[UD.InitialIndex+4]->vzf = false;
        D.Particles[UD.InitialIndex+5]->Ff = -force;
        D.Particles[UD.InitialIndex+5]->FixVeloc();
        D.Particles[UD.InitialIndex+5]->vzf = false;
        UD.pSig(2) = true;
    }

    // info
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    
}

void ResetEps(DEM::Domain const & dom, UserData & UD)
{
    UD.L0(0) = dom.Particles[UD.InitialIndex  ]->x(0)-dom.Particles[UD.InitialIndex+1]->x(0);
    UD.L0(1) = dom.Particles[UD.InitialIndex+2]->x(1)-dom.Particles[UD.InitialIndex+3]->x(1);
    UD.L0(2) = dom.Particles[UD.InitialIndex+4]->x(2)-dom.Particles[UD.InitialIndex+5]->x(2);
}

void Setup (DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.StrainCtrl)
    {
        if (!dat.pSig(0))
        {
            double area = (dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1))*(dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2));
            double sig = -0.5*(dom.Particles[dat.InitialIndex  ]->F(0)-dom.Particles[dat.InitialIndex+1]->F(0))/area;
            double dsig = sig - dat.Sig0(0);
            double r = dsig/((2.0/3.0)*sin(dat.Alp)*sin(dat.Thf-2.0*Util::PI/3.0)-cos(dat.Alp));
            dat.Sig(0) = dat.Sig0(0)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf-2.0*Util::PI/3.0);
            dat.Sig(1) = dat.Sig0(1)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf);
            dat.Sig(2) = dat.Sig0(2)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf+2.0*Util::PI/3.0);
        }
        if (!dat.pSig(1))
        {
            double area = (dom.Particles[dat.InitialIndex]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0))*(dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2));
            double sig = -0.5*(dom.Particles[dat.InitialIndex+2]->F(1)-dom.Particles[dat.InitialIndex+3]->F(1))/area;
            double dsig = sig - dat.Sig0(1);
            double r = dsig/((2.0/3.0)*sin(dat.Alp)*sin(dat.Thf)-cos(dat.Alp));
            dat.Sig(0) = dat.Sig0(0)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf-2.0*Util::PI/3.0);
            dat.Sig(1) = dat.Sig0(1)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf);
            dat.Sig(2) = dat.Sig0(2)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf+2.0*Util::PI/3.0);
        }
        if (!dat.pSig(2))
        {
            double area = (dom.Particles[dat.InitialIndex]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0))*(dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1));
            double sig = -0.5*(dom.Particles[dat.InitialIndex+4]->F(2)-dom.Particles[dat.InitialIndex+5]->F(2))/area;
            double dsig = sig - dat.Sig0(2);
            double r = dsig/((2.0/3.0)*sin(dat.Alp)*sin(dat.Thf+2.0*Util::PI/3.0)-cos(dat.Alp));
            dat.Sig(0) = dat.Sig0(0)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf-2.0*Util::PI/3.0);
            dat.Sig(1) = dat.Sig0(1)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf);
            dat.Sig(2) = dat.Sig0(2)-r*cos(dat.Alp) + (2.0/3.0)*r*sin(dat.Alp)*sin(dat.Thf+2.0*Util::PI/3.0);
        }
    }
    Vec3_t force;
    bool   update_sig = false;
    if (dat.pSig(0))
    {
        double area = (dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1))*(dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2));
        force = dat.Sig(0)*area, 0.0, 0.0;
        dom.Particles[dat.InitialIndex  ]->Ff =  force;
        dom.Particles[dat.InitialIndex+1]->Ff = -force;
        if (!dat.StrainCtrl) update_sig = true;
    }
    else if (!dat.StrainCtrl)
    {
        double area = (dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1))*(dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2));
        dat.Sig(0) = -0.5*(dom.Particles[dat.InitialIndex  ]->F(0)-dom.Particles[dat.InitialIndex+1]->F(0))/area;
    }
    if (dat.pSig(1))
    {
        double area = (dom.Particles[dat.InitialIndex]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0))*(dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2));
        force = 0.0, dat.Sig(1)*area, 0.0;
        dom.Particles[dat.InitialIndex+2]->Ff =  force;
        dom.Particles[dat.InitialIndex+3]->Ff = -force;
        if (!dat.StrainCtrl) update_sig = true;
    }
    else if (!dat.StrainCtrl)
    {
        double area = (dom.Particles[dat.InitialIndex]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0))*(dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2));
        dat.Sig(1) = -0.5*(dom.Particles[dat.InitialIndex+2]->F(1)-dom.Particles[dat.InitialIndex+3]->F(1))/area;
    }
    if (dat.pSig(2))
    {
        double area = (dom.Particles[dat.InitialIndex]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0))*(dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1));
        force = 0.0, 0.0, dat.Sig(2)*area;
        dom.Particles[dat.InitialIndex+4]->Ff =  force;
        dom.Particles[dat.InitialIndex+5]->Ff = -force;
        if (!dat.StrainCtrl) update_sig = true;
    }
    else if (!dat.StrainCtrl)
    {
        double area = (dom.Particles[dat.InitialIndex]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0))*(dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1));
        dat.Sig(2) = -0.5*(dom.Particles[dat.InitialIndex+4]->F(2)-dom.Particles[dat.InitialIndex+5]->F(2))/area;
    }
    if (update_sig) dat.Sig += dat.dt*dat.DSig/(dat.tspan);
}

void Report (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
/*    if(dat.RenderVideo)
    {
    //Output force vectors for each interacton
        String ff;
        ff.Printf    ("%s_%08d_chains.res",dom.FileKey.CStr(), dom.idx_out);
        std::ofstream FF(ff.CStr());
        FF <<  Util::_10_6 << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" <<  Util::_8s << "Ftx" << Util::_8s << "Fty" << Util::_8s << "Ftz" 
           <<    Util::_8s << "Nc" << Util::_8s << "Nsc"<< Util::_8s << "P1" <<  Util::_8s << "P2" <<"\n";

        for (size_t i=0; i<dom.CInteractons.Size(); i++)
        {
            if ((norm(dom.CInteractons[i]->Fnet)>0.0)&&(dom.CInteractons[i]->P1->IsFree()&&dom.CInteractons[i]->P2->IsFree()))
            {
                FF << Util::_8s << dom.CInteractons[i]->Fnet(0)   << Util::_8s << dom.CInteractons[i]->Fnet(1)   << Util::_8s <<  dom.CInteractons[i]->Fnet(2) 
                   << Util::_8s << dom.CInteractons[i]->Ftnet(0)  << Util::_8s << dom.CInteractons[i]->Ftnet(1)  << Util::_8s <<  dom.CInteractons[i]->Ftnet(2) 
                   << Util::_8s << dom.CInteractons[i]->Nc        << Util::_8s << dom.CInteractons[i]->Nsc
                   << Util::_8s << dom.CInteractons[i]->P1->Index << Util::_8s << dom.CInteractons[i]->P2->Index << "\n";
            }
        }
        FF.close();
        //Output cohesion interactons information
        if (dom.BInteractons.Size()>0)
        {
            String fc;
            fc.Printf    ("%s_%08d_bonds.res",dom.FileKey.CStr(), dom.idx_out);
            std::ofstream FC(fc.CStr());
            FC <<  Util::_10_6 << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" <<  Util::_8s << "Area" << Util::_8s << "valid" << Util::_8s << "P1" <<  Util::_8s << "P2" <<"\n";

            for (size_t i=0; i<dom.BInteractons.Size(); i++)
            {
                FC << Util::_8s << dom.BInteractons[i]->Fnet(0)   << Util::_8s << dom.BInteractons[i]->Fnet(1)   << Util::_8s <<  dom.BInteractons[i]->Fnet(2) 
                   << Util::_8s << dom.BInteractons[i]->Area      << Util::_8s << dom.BInteractons[i]->valid
                   << Util::_8s << dom.BInteractons[i]->P1->Index << Util::_8s << dom.BInteractons[i]->P2->Index << "\n";
            }
            FC.close();

            String fl;
            fl.Printf    ("%s_%08d_clusters.res",dom.FileKey.CStr(), dom.idx_out);
            std::ofstream FL(fl.CStr());
            FL <<  Util::_10_6 << "Cluster" << Util::_8s << "Mass" << "\n";
            FL <<  Util::_8s   << 0         << Util::_8s << dom.Ms << "\n";
            for (size_t i=0;i<dom.Listofclusters.Size();i++)
            {
                double m = 0.0; //mass of the cluster
                for (size_t j=0;j<dom.Listofclusters[i].Size();j++)
                {
                    m+=dom.Particles[dom.Listofclusters[i][j]]->Props.m;
                }
                FL << Util::_8s << i+1 << Util::_8s << m << "\n"; 
            }
            FL.close();
        }
        //Ouput particle status at each time step
        String fv;
        fv.Printf    ("%s_%08d_particles.res",dom.FileKey.CStr(), dom.idx_out);
        std::ofstream FV(fv.CStr());
        FV <<  Util::_10_6 << "PID" << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" <<"\n";

        for (size_t i=0; i<dom.Particles.Size(); i++)
        {
            if (dom.Particles[i]->IsFree())
            {
                FV << Util::_8s << dom.Particles[i]->Index << Util::_8s << dom.Particles[i]->v(0)   << Util::_8s << dom.Particles[i]->v(1) << Util::_8s << dom.Particles[i]->v(2)
                   << Util::_8s << dom.Particles[i]->x(0)  << Util::_8s << dom.Particles[i]->x(1)   << Util::_8s << dom.Particles[i]->x(2) << "\n";
            }
        }
        FV.close();
        // output triaxial test data
    }
*/
  
    // header
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        // Output of the current time, the stress state sx,sy,sz the strain state ex,ey,ez the void ratio, the coordination number the number of
        // contacts and sliding contacts Nc,Nsc and the number of bonds and broken bonds Nb Nbb
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sx" << Util::_8s << "sy" << Util::_8s << "sz";
        dat.oss_ss <<                          Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez";
        dat.oss_ss << Util::_8s   << "e"                         << Util::_8s << "Nc" << Util::_8s << "Nsc";         
        dat.oss_ss <<                                               Util::_8s << "Nb" << Util::_8s << "Nbb" << "\n";
    }
    if (dat.RenderVideo)
    {
        String ff;
        ff.Printf    ("%s_bf_%04d",dom.FileKey.CStr(), dom.idx_out);
        //dom.WriteVTKContacts (ff.CStr());
        dom.WriteBF(ff.CStr());
    }
    if (!dom.Finished) 
    {
        // stress
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dat.Sig(0) << Util::_8s << dat.Sig(1) << Util::_8s << dat.Sig(2);

        // strain
        dat.oss_ss << Util::_8s << (dom.Particles[dat.InitialIndex  ]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0)-dat.L0(0))/dat.L0(0);
        dat.oss_ss << Util::_8s << (dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1)-dat.L0(1))/dat.L0(1);
        dat.oss_ss << Util::_8s << (dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2)-dat.L0(2))/dat.L0(2);

        // void ratio
        double volumecontainer = (dom.Particles[dat.InitialIndex  ]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0)-dom.Particles[dat.InitialIndex  ]->Props.R+dom.Particles[dat.InitialIndex+1]->Props.R)*
                                 (dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1)-dom.Particles[dat.InitialIndex+2]->Props.R+dom.Particles[dat.InitialIndex+3]->Props.R)*
                                 (dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2)-dom.Particles[dat.InitialIndex+4]->Props.R+dom.Particles[dat.InitialIndex+5]->Props.R);

        dat.oss_ss << Util::_8s << (volumecontainer-dom.Vs)/dom.Vs;

        // Number of contacts Nc, number of sliding contacts Nsc and Coordination number Cn
        size_t Nc = 0;
        size_t Nsc = 0;
        for (auto it=dom.PairtoCInt.begin();it!=dom.PairtoCInt.end();++it)
        {
            if(it->second->I2<dat.InitialIndex)
            {
                Nc  += it->second->Nc;
                Nsc += it->second->Nsc;
            }
        }

        size_t Nb  = dom.BInteractons.Size(); //Number of bonds
        size_t Nbb = 0;                   //Number of broken bonds
        for (size_t i=0; i<dom.BInteractons.Size(); i++)
        {
            if(!dom.BInteractons[i]->valid)
            {
                Nbb ++;
            }
        }

        dat.oss_ss << Util::_8s << Nc << Util::_8s << Nsc << Util::_8s << Nb << Util::_8s << Nbb;

        dat.oss_ss << std::endl;
    }
    else
    {
        dat.oss_ss.close();
        String fn;
        fn.Printf("%s_forces.res",dom.FileKey.CStr());
        std::ofstream OF(fn.CStr());
        OF <<  Util::_10_6 << "Fn" << Util::_8s << "Ft" << Util::_8s << "NContacts" << Util::_8s << "Issliding" << "\n";

        String f;
        f.Printf("%s_stress.res",dom.FileKey.CStr());
        std::ofstream SF(f.CStr());
        Mat3_t S,B;
        for (size_t m=0;m<3;m++)
        {
            for (size_t n=0;n<3;n++)
            {
                S(m,n)=0.0;
                B(m,n)=0.0;
            }
        }
        double volumecontainer = (dom.Particles[dat.InitialIndex  ]->x(0)-dom.Particles[dat.InitialIndex+1]->x(0)-dom.Particles[dat.InitialIndex  ]->Props.R-dom.Particles[dat.InitialIndex+1]->Props.R)*
                                 (dom.Particles[dat.InitialIndex+2]->x(1)-dom.Particles[dat.InitialIndex+3]->x(1)-dom.Particles[dat.InitialIndex+2]->Props.R-dom.Particles[dat.InitialIndex+3]->Props.R)*
                                 (dom.Particles[dat.InitialIndex+4]->x(2)-dom.Particles[dat.InitialIndex+5]->x(2)-dom.Particles[dat.InitialIndex+4]->Props.R-dom.Particles[dat.InitialIndex+5]->Props.R);
        for (auto it=dom.PairtoCInt.begin();it!=dom.PairtoCInt.end();++it)
        {
            DEM::CInteracton * CI = it->second;
            //if (CI->Nc>0)
            if (CI->Nc>0&&CI->P1->IsFree()&&CI->P2->IsFree())
            {
                Vec3_t branch    = CI->P2->x - CI->P1->x;
                for (size_t m=0;m<3;m++)
                {
                    for (size_t n=0;n<3;n++)
                    {
                        S(m,n) += (CI->Fnet(m) + CI->Ftnet(m))*branch(n)/volumecontainer;
                    }
                }

            }
        }
        size_t Ncontacts = 0;
        for (auto it=dom.PairtoCInt.begin();it!=dom.PairtoCInt.end();++it)
        {
            if (norm(it->second->Fnet)>0.0&&it->second->P1->IsFree()&&it->second->P2->IsFree())
            {
                //OF << Util::_10_6 << norm(dom.CInteractons[i]->Fnet) << Util::_8s << norm(dom.CInteractons[i]->Ftnet) << Util::_8s <<  dom.CInteractons[i]->Nc << Util::_8s <<  dom.CInteractons[i]->Nsc << "\n";
                OF << Util::_10_6 << it->second->Fnet(0) << Util::_8s << it->second->Fnet(1) << Util::_8s << it->second->Fnet(2) << Util::_8s <<  "\n";
                //for (size_t m=0;m<3;m++)
                //{
                    //for (size_t n=0;n<3;n++)
                    //{
                        //B(m,n)+=dom.CInteractons[i]->B(m,n);
                    //}
                //}
                Ncontacts+=it->second->Nc;
            }
        }
        OF.close();



        for (size_t m=0;m<3;m++)
        {
            for (size_t n=0;n<3;n++)
            {
                SF << Util::_8s << S(m,n);
            }
            SF << std::endl;
        }
        for (size_t m=0;m<3;m++)
        {
            for (size_t n=0;n<3;n++)
            {
                SF << Util::_8s << B(m,n)/Ncontacts;
            }
            SF << std::endl;
        }
        SF.close();
    }
}

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    //number of threads
    size_t Nproc         = 1;
    bool   mostlyspheres = false;
    double binsize       = 2.0;
    if (argc>=3) Nproc         = atoi(argv[2]);
    if (argc>=4) mostlyspheres = atoi(argv[3]);
    if (argc>=5) binsize       = atof(argv[4]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    
    double verlet;      // Verlet distance for optimization
    String ptype;       // Particle type 
    size_t RenderVideo; // Decide is video should be render
    bool   Cohesion;    // Decide if coheison is going to be simulated
    double fraction;    // Fraction of particles to be generated
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Beta;        // Rolling stiffness coefficient (only for spheres)
    double Eta;         // Plastic moment coefficient (only for spheres, 0 if rolling resistance is not used)
    double Bn;          // Cohesion normal stiffness
    double Bt;          // Cohesion tangential stiffness
    double Bm;          // Cohesion torque stiffness
    double Eps;         // Threshold for breking bonds
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt;          // Time step
    double dtOut;       // Time step for output
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double rho;         // rho
    double p0;          // Pressure for the isotropic compression
    double T0;          // Time span for the compression
    bool   isfailure;   // Flag for a failure stress path
    bool   pssrx;       // Prescribed strain rate in X ?
    bool   pssry;       // Prescribed strain rate in Y ?
    bool   pssrz;       // Prescribed strain rate in Z ?
    double srx;         // Final Strain x
    double sry;         // Final Strain y
    double srz;         // Final Strain z
    double pf;          // Final pressure p
    double qf;          // Final deviatoric stress q
    double thf;         // Angle of the stress path alpha
    double alpf;        // Angle of the p q plane
    double Tf;          // Final time for the test
    {
        infile >> verlet;       infile.ignore(200,'\n');
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> Cohesion;     infile.ignore(200,'\n');
        infile >> fraction;     infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu;           infile.ignore(200,'\n');
        infile >> Beta;         infile.ignore(200,'\n');
        infile >> Eta;          infile.ignore(200,'\n');
        infile >> Bn;           infile.ignore(200,'\n');
        infile >> Bt;           infile.ignore(200,'\n');
        infile >> Bm;           infile.ignore(200,'\n');
        infile >> Eps;          infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> seed;         infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> dtOut;        infile.ignore(200,'\n');
        infile >> Lx;           infile.ignore(200,'\n');
        infile >> Ly;           infile.ignore(200,'\n');
        infile >> Lz;           infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
        infile >> p0;           infile.ignore(200,'\n');
        infile >> T0;           infile.ignore(200,'\n');
        infile >> isfailure;    infile.ignore(200,'\n');
        infile >> pssrx;        infile.ignore(200,'\n');
        infile >> pssry;        infile.ignore(200,'\n');
        infile >> pssrz;        infile.ignore(200,'\n');
        infile >> srx;          infile.ignore(200,'\n');
        infile >> sry;          infile.ignore(200,'\n');
        infile >> srz;          infile.ignore(200,'\n');
        infile >> pf;           infile.ignore(200,'\n');
        infile >> qf;           infile.ignore(200,'\n');
        infile >> thf;          infile.ignore(200,'\n');
        infile >> alpf;         infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
    }

    // domain and User data
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha=verlet;
    dom.Beta = binsize;
    dom.MostlySpheres = mostlyspheres;
    dat.dt = dt;
    dat.RenderVideo = (bool) RenderVideo;


    bool load = false;
    // particle
    if      (ptype=="sphere")    dom.GenSpheres  (-1, Lx, nx, rho, "HCP", seed, fraction, Eps);
    else if (ptype=="sphereboxnormal") 
    {
        Vec3_t Xmin(-0.5*Lx,-0.5*Ly,-0.5*Lz);
        Vec3_t Xmax = -Xmin;
        dom.GenSpheresBox (-1, Xmin, Xmax, R, rho, "Normal", seed, fraction, Eps);
    }
    else if (ptype=="sphereboxhcp") 
    {
        Vec3_t Xmin(-0.5*Lx,-0.5*Ly,-0.5*Lz);
        Vec3_t Xmax = -Xmin;
        dom.GenSpheresBox (-1, Xmin, Xmax, R, rho, "HCP",    seed, fraction, Eps);
    }
    else if (ptype=="voronoi")   dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, Cohesion, !Cohesion, seed, fraction);
    else if (ptype=="tetra")
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        dom.GenFromMesh (mesh,/*R*/R,/*rho*/rho,Cohesion,false);
    }
    else if (ptype=="rice") dom.GenRice(-1,Lx,nx,R,rho,seed,fraction);
    else
    {
        dom.Load(ptype.CStr());
        Array<int> DeleteTags(6);
        DeleteTags = -2,-3,-4,-5,-6,-7;
        dom.DelParticles(DeleteTags);
        load = true;
    }
        
     //throw new Fatal("Packing for particle type not implemented yet");
    dat.InitialIndex = dom.Particles.Size();
    dom.GenBoundingBox (/*InitialTag*/-2, 0.02*R, /*Cf*/1.3,Cohesion);
    
    //Fixing some degrees of freedom for the particles, just for the biaxial test, not for the official version of the test
    //Array<DEM::Particle *> Grains;
    //dom.GetParticles(-1,Grains);
    //for (size_t i=0;i<Grains.Size();i++)
    //{
        //Grains[i]->FixVeloc();
        //Grains[i]->vxf = false;
        //Grains[i]->vyf = false;
        //Grains[i]->vzf = false;
    //}



    // properties of particles prior the triaxial test
    Dict B;
    B.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    //B.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,0.0,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    B.Set(-2,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,0.0,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-3,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,0.0,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-4,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,0.0,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-5,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,0.0,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-6,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,0.0,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-7,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,0.0,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    dom.SetProps(B);

    // stage 1: isotropic compresssion  //////////////////////////////////////////////////////////////////////
    String fkey_a(filekey+"_a");
    String fkey_b(filekey+"_b");
    Vec3_t  sigf;                      // final stress state
    bVec3_t peps(false, false, false); // prescribed strain rates ?
    Vec3_t  depsdt(0.0,0.0,0.0);       // strain rate

    sigf =  Vec3_t(-p0,-p0,-p0);
    if (load) dat.Sig = sigf;
    ResetEps  (dom,dat);
    SetTxTest (sigf, peps, depsdt,0,0,false,dat,dom);
    dat.tspan = T0/2.0 - dom.Time;
    dom.Solve  (/*tf*/T0/2.0, /*dt*/dt, /*dtOut*/dtOut, &Setup, &Report, fkey_a.CStr(),RenderVideo,Nproc);
    dom.Save(fkey_a.CStr());
    SetTxTest (sigf, peps, depsdt,0,0,false,dat,dom);
    dat.tspan = T0 - dom.Time;
    //Dict D;
    //D.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    //dom.SetProps(D);
    dom.Solve (/*tf*/T0, /*dt*/dt, /*dtOut*/dtOut, &Setup, &Report, fkey_b.CStr(),RenderVideo,Nproc);
    dom.Save(fkey_b.CStr());



    // stage 2: The proper triaxial test /////////////////////////////////////////////////////////////////////////
    String fkey_c(filekey+"_c");
    Vec3_t lf;
    pqth2L (pf, qf, thf, lf, "cam");
    sigf   = lf(0), lf(1), lf(2);
    peps   = bVec3_t(pssrx, pssry, pssrz);
    depsdt = Vec3_t(srx/(Tf-dom.Time), sry/(Tf-dom.Time), srz/(Tf-dom.Time));
    
    // run
    ResetEps  (dom,dat);
    SetTxTest (sigf, peps, depsdt, thf*M_PI/180, alpf*M_PI/180, isfailure, dat, dom);
    dat.tspan = Tf - dom.Time;
    dom.Solve     (/*tf*/Tf, /*dt*/dt, /*dtOut*/dtOut, &Setup, &Report, fkey_c.CStr(),RenderVideo,Nproc);
    dom.Save(fkey_c.CStr());

    return 0;
}
MECHSYS_CATCH
