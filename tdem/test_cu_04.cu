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
 * along with this program. If not, see <ttp://www.gnu.org/licenses/>  *
 ************************************************************************/

/////////////////////// Test 04 shearbox to test periiodic boundary conditions

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

struct UdCu //User data for cuda
{
    real   Lz;
    real   x0;
    real   sawL;
    real   Tf;
    real   T0;
    real   str;
    real   mtop;
    real   ztop;
    real   zbtop;
    real   Ftop;
    real   Dz;
    real   Vz;
    size_t i1;
    size_t i2;
    size_t nptop;
    size_t npbot;
};

struct UserData
{
    String test;
    double L0;
    UdCu * pUD;                      ///< Pointer to user data strcuture for cuda
    bool * pcrossed;                 ///< a bool varible to see if the top lid moved beyond the threshold
    thrust::device_vector<size_t> TopPar; 
    thrust::device_vector<size_t> BotPar; 
    Array<size_t> TPar;
    Array<size_t> BPar;
    size_t * pTopPar;
    size_t * pBotPar;
    size_t nptop;
    size_t npbot;
    std::ofstream      oss_ss;       ///< file for stress strain data
};

__global__ void AddTopForce(size_t const * TopPar, DEM::DynParticleCU const * DPar, UdCu * UC)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=UC[0].nptop) return;
    size_t ip = TopPar[ic];
    atomicAdd(&UC[0].Ftop,DPar[ip].F.z);
}

__global__ void MoveTop(UdCu * UC,DEM::dem_aux const * demaux)
{
    real za   = 2.0*UC[0].ztop - UC[0].zbtop + UC[0].Ftop*demaux[0].dt*demaux[0].dt/UC[0].mtop;
    UC[0].Dz  = za - UC[0].ztop;
    UC[0].Vz  = 0.5*(za - UC[0].zbtop)/demaux[0].dt;
    UC[0].zbtop = UC[0].ztop;
    UC[0].ztop  = za;
    UC[0].Ftop  = 0.0;
}

__global__ void MoveTopPar(size_t * TopPar,real3 * Verts, real3 * Vertso, DEM::ParticleCU * Par, DEM::DynParticleCU * DPar, UdCu * UC, DEM::dem_aux const * demaux)
{
    size_t ic = threadIdx.x + blockIdx.x * blockDim.x;
    if (ic>=UC[0].nptop) return;
    size_t ip = TopPar[ic];
    real  vel = UC[0].str*(DPar[UC[0].i2].x.z-DPar[UC[0].i1].x.z-UC[0].Lz)*fmin(10.0*(demaux[0].Time-UC[0].T0)/(UC[0].Tf-UC[0].T0),1.0);
    DPar[ip].v.x  = vel;
    DPar[ip].xb.x = DPar[ip].x.x-vel*demaux[0].dt;
    DPar[ip].v.z  = UC[0].Vz;
    real3 trans = make_real3(0.0,0.0,UC[0].Dz);
    Translate(Verts,Vertso,Par[ip],DPar[ip],trans);
}

__global__ void SetupSaw(real3 * Verts, real3 * Vertso, DEM::ParticleCU * Par, DEM::DynParticleCU * DPar, bool * pcrossed, UdCu * UC, DEM::dem_aux const * demaux)
{
    //printf("UC[0] %g %g %g \n",UC[0].x0,UC[0].str,UC[0].Lz);
    //printf("DPar %g %g %g \n",DPar[UC[0].i2].x.x,DPar[UC[0].i2].x.y,DPar[UC[0].i2].x.z);
    real  vel = UC[0].str*(DPar[UC[0].i2].x.z-DPar[UC[0].i1].x.z-UC[0].Lz)*fmin(10.0*(demaux[0].Time-UC[0].T0)/(UC[0].Tf-UC[0].T0),1.0);
    DPar[UC[0].i2].v.x  = vel;
    DPar[UC[0].i2].xb.x = DPar[UC[0].i2].x.x-vel*demaux[0].dt;
    if (DPar[UC[0].i2].x.x > UC[0].x0 + UC[0].sawL)
    {
        pcrossed[0] = true;
        real3 trans = make_real3(-UC[0].sawL,0.0,0.0);
        Translate(Verts,Vertso,Par[UC[0].i2],DPar[UC[0].i2],trans);
    }
    else pcrossed[0] = false;
}

void AddSawPlate(DEM::Domain & dom, int Tag,  Vec3_t & X, double Lx, double Ly, size_t Ntooth, double depth, double rho, double R)
{
     Array<Vec3_t>        V(4*Ntooth+2);
     Array<Array <int> >  E(2*Ntooth+1);
     Array<Array <int> >  F(2*Ntooth);
     double step = Lx/Ntooth;
     for (size_t i=0;i<Ntooth+1;i++)
     {
         V[i             ] = Vec3_t(i*step,-0.5*Ly,0.0);
         V[i + Ntooth + 1] = Vec3_t(i*step, 0.5*Ly,0.0);
         E[i].Push(i);
         E[i].Push(i+Ntooth+1);
     }
     //std::cout << "1" << std::endl;
     for (size_t i=0;i<Ntooth;i++)
     {
         V[i + 2*Ntooth + 2] = Vec3_t(i*step+0.5*step,-0.5*Ly,depth);
         V[i + 3*Ntooth + 2] = Vec3_t(i*step+0.5*step, 0.5*Ly,depth);
         E[i + Ntooth + 1].Push(i + 2*Ntooth + 2);
         E[i + Ntooth + 1].Push(i + 3*Ntooth + 2);
         F[2*i  ].Push(i               );
         F[2*i  ].Push(i +   Ntooth + 1);
         F[2*i  ].Push(i + 3*Ntooth + 2);
         F[2*i  ].Push(i + 2*Ntooth + 2);
         F[2*i+1].Push(i + 1           );
         F[2*i+1].Push(i +   Ntooth + 2);
         F[2*i+1].Push(i + 3*Ntooth + 2);
         F[2*i+1].Push(i + 2*Ntooth + 2);
     }
     //std::cout << "2" << std::endl;
     dom.Particles.Push(new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
     //std::cout << "3" << std::endl;
     dom.Particles[dom.Particles.Size()-1]->Q          = 1.0,0.0,0.0,0.0;
     dom.Particles[dom.Particles.Size()-1]->Props.V    = Lx*Ly*R;
     dom.Particles[dom.Particles.Size()-1]->Props.m    = rho*Lx*Ly*R;
     dom.Particles[dom.Particles.Size()-1]->I          = 1.0,1.0,1.0;
     dom.Particles[dom.Particles.Size()-1]->I         *= dom.Particles[dom.Particles.Size()-1]->Props.m;
     dom.Particles[dom.Particles.Size()-1]->x          = Vec3_t(0.5*Lx,0.0,0.5*depth);
     dom.Particles[dom.Particles.Size()-1]->Ekin       = 0.0;
     dom.Particles[dom.Particles.Size()-1]->Erot       = 0.0;
     dom.Particles[dom.Particles.Size()-1]->Dmax       = 0.5*sqrt(Lx*Lx+Ly*Ly+depth*depth)+R;
     dom.Particles[dom.Particles.Size()-1]->PropsReady = true;
     dom.Particles[dom.Particles.Size()-1]->Index      = dom.Particles.Size()-1;

     dom.Particles[dom.Particles.Size()-1]->Position(X);
}

void SetupSaw (DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    SetupSaw<<<1,1>>>(dom.pVertsCU,dom.pVertsoCU,dom.pParticlesCU,dom.pDynParticlesCU,dat.pcrossed,dat.pUD,dom.pdemaux);
    bool * crossed = new bool[1];
    cudaMemcpy(crossed,dat.pcrossed, sizeof(bool), cudaMemcpyDeviceToHost);
    if (crossed[0])
    {
        dom.UpdateContactsDevice();
    }
}

void ReportSaw (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sx" << Util::_8s << "ey \n";
    }
    if (!dom.Finished) 
    {
        //double Areaz = (dom.GetParticle(-2)->x(1)-dom.GetParticle(-3)->x(1))*(dom.Xmax - dom.Xmin);
        double Areaz = (dom.Ymax-dom.Ymin)*fabs(dom.Xmax - dom.Xmin);
        //double Sx    = 0.5*(dom.GetParticle(-4)->F(0) - dom.GetParticle(-5)->F(0))/Areaz;
        double Sx    = -dom.GetParticle(-5)->F(0)/Areaz;
        double Ey    = (dom.GetParticle(-5)->x(2) - dom.GetParticle(-4)->x(2) - dat.L0)/dat.L0;
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << Sx << Util::_8s << Ey << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

void SetupNor (DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    AddTopForce<<<dat.nptop/dom.Nthread+1,dom.Nthread>>>(dat.pTopPar,dom.pDynParticlesCU,dat.pUD);

    MoveTop<<<1,1>>>(dat.pUD,dom.pdemaux);

    MoveTopPar<<<dat.nptop/dom.Nthread+1,dom.Nthread>>>(dat.pTopPar,dom.pVertsCU,dom.pVertsoCU,dom.pParticlesCU,dom.pDynParticlesCU,dat.pUD,dom.pdemaux);
}

void ReportNor (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sx" << Util::_8s << "ey \n";
    }
    if (!dom.Finished) 
    {
        double Areaz = (dom.Ymax-dom.Ymin)*fabs(dom.Xmax - dom.Xmin);
        double Fx = 0.0;
        for (size_t ip=0;ip<dat.TPar.Size();ip++)
        {
            size_t id = dat.TPar[ip];
            Fx += dom.Particles[id]->F(0);
        }
        double Sx    = -Fx/Areaz;
        double Ey    = (dom.GetParticle(-5)->x(2) - dom.GetParticle(-4)->x(2) - dat.L0)/dat.L0;
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << Sx << Util::_8s << Ey << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    //number of threads
    size_t Nproc = 1; 
    double dthoot= 4.0;     
    double qin   = 0.0;
    if (argc>=3) Nproc=atoi(argv[2]);
    if (argc>=4) dthoot=atof(argv[3]);
    if (argc>=5) qin   =atof(argv[4]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    
    //-----------
    // select contact model for `CInteracton`
    DEM::CInteracton_PointerToContactFunction = &DEM::CInteracton_UseLinearContact;
    DEM::CInteractonSphere_PointerToContactFunction = &DEM::CInteractonSphere_UseLinearContact;
    DEM::CInteracton_CalcForce_PointerToContactFunction = &DEM::CInteractonSphere_CalcForce_UseLinearContact;
    cudaMemcpyFromSymbol(&DEM::deviceFunction_inside_CalcForceVV, DEM::devFuncPtrLinearContact, sizeof(DEM::Device_CalcForceVV_ContactPointer));

    // DEM::CInteracton_PointerToContactFunction = &DEM::CInteracton_UseHertzContact;
    // DEM::CInteractonSphere_PointerToContactFunction = &DEM::CInteractonSphere_UseHertzContact;
    // DEM::CInteracton_CalcForce_PointerToContactFunction = &DEM::CInteractonSphere_CalcForce_UseHertzContact;
    // cudaMemcpyFromSymbol(&DEM::deviceFunction_inside_CalcForceVV, DEM::devFuncPtrHertzContact, sizeof(DEM::Device_CalcForceVV_ContactPointer));

    double verlet;      // Verlet distance for optimization
    String ptype;       // Particle type 
    String test;        // Particle type 
    size_t  RenderVideo;// Decide is video should be render
    bool   Cohesion;    // Decide if coheison is going to be simulated
    double fraction;    // Fraction of particles to be generated
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu0;         // Microscopic friction coefficient at the initial stage
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
    double str;         // Strain rate for shearing
    double T0;          // Time span for the compression
    double Tf;          // Final time for the test
    {
        infile >> verlet;       infile.ignore(200,'\n');
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> test;         infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> Cohesion;     infile.ignore(200,'\n');
        infile >> fraction;     infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu0;          infile.ignore(200,'\n');
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
        infile >> str;          infile.ignore(200,'\n');
        infile >> T0;           infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
    }

    // domain and User data
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha=verlet;
    dom.Dilate = true;

    bool load = false;
    // particle
    if      (ptype=="sphere")  dom.GenSpheres  (-1, Lx, nx, rho, "HCP", seed, fraction);
    else if (ptype=="sphereboxhcp") 
    {
        Vec3_t Xmin(-0.5*Lx,-0.5*Ly,-0.5*Lz);
        Vec3_t Xmax = -Xmin;
        dom.GenSpheresBox (-1, Xmin, Xmax, R, rho, "HCP",    seed, fraction, Eps);
    }
    else if (ptype=="voronoi")
    {
        if (ny==1) dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, Cohesion, bVec3_t(true,false,true), seed, fraction, Vec3_t(qin,1.0,qin));
        else       dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, Cohesion, bVec3_t(true,true ,true), seed, fraction, Vec3_t(qin,qin,qin));
    }
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
        load = true;
    }
    
    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    double a = 0.1;
    double maxd = dom.MaxDim();
    if (test=="normal")
    {
        Vec3_t X0 = Vec3_t(0.0,0.5*(Xmin(1)+Xmax(1)),Xmin(2)-a*maxd);
        Vec3_t X1 = Vec3_t(0.0,0.5*(Xmin(1)+Xmax(1)),Xmax(2)+a*maxd);
        dom.AddPlane(-4,X0,a*maxd,Lx*2.0,Xmax(1)-Xmin(1)+4.0*maxd,3.0,a*maxd);
        dom.AddPlane(-5,X1,a*maxd,Lx*2.0,Xmax(1)-Xmin(1)+4.0*maxd,3.0,a*maxd);
    }
    else if (test=="sawtooth")
    {
        size_t nthoot = 8*size_t(Lx/(4.0*dthoot*maxd));
        Vec3_t X0 = Vec3_t(0.0,0.5*(Xmin(1)+Xmax(1)),Xmin(2)-a*maxd-0.5*dthoot*maxd);
        Vec3_t X1 = Vec3_t(0.0,0.5*(Xmin(1)+Xmax(1)),Xmax(2)+a*maxd+0.5*dthoot*maxd);
        AddSawPlate(dom,-4,X0,Lx*2.0,Xmax(1)-Xmin(1)+8.0*maxd,nthoot,dthoot*maxd,3.0,a*maxd);
        AddSawPlate(dom,-5,X1,Lx*2.0,Xmax(1)-Xmin(1)+8.0*maxd,nthoot,dthoot*maxd,3.0,a*maxd);
        Quaternion_t q;
        NormalizeRotation (M_PI,OrthoSys::e1,q);
        dom.GetParticle(-5)->Rotate(q,dom.GetParticle(-5)->x);
    }


    dom.Xmax =  0.5*Lx;
    dom.Xmin = -0.5*Lx;
    dom.Ymax =  0.5*Ly;
    dom.Ymin = -0.5*Ly;

    Dict B1;
    B1.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt ,Mu0     ,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    //B1.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt ,0.0     ,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    //B1.Set(-2,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    //B1.Set(-3,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B1.Set(-4,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B1.Set(-5,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    dom.SetProps(B1);

    //dom.GetParticle(-2)->FixVeloc();
    //dom.GetParticle(-3)->FixVeloc();
    dom.GetParticle(-4)->FixVeloc();
    dom.GetParticle(-5)->FixVeloc();
    dom.GetParticle(-4)->vzf = false;
    dom.GetParticle(-5)->vzf = false;

    //double Areaz = (dom.GetParticle(-2)->x(1)-dom.GetParticle(-3)->x(1))*Lx;
    double Areaz = Lx*Ly;
    dom.GetParticle(-4)->Ff = Vec3_t(0.0,0.0, p0*Areaz);
    dom.GetParticle(-5)->Ff = Vec3_t(0.0,0.0,-p0*Areaz);

    String fkey_a  (filekey+"_a");
    String fkey_b  (filekey+"_b");
    String fkeybf_a(filekey+"bf_a");
    String fkeybf_b(filekey+"bf_b");

    dom.Solve (/*tf*/T0, /*dt*/dt, /*dtOut*/dtOut, NULL, NULL, fkey_a.CStr(),RenderVideo,Nproc);
    dom.Save(fkey_a.CStr());
    dom.WriteXDMF(fkey_a.CStr());
    dom.WriteBF(fkeybf_a.CStr());

    Dict B2;
    B2.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt ,Mu     ,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    //B2.Set(-2,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    //B2.Set(-3,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B2.Set(-4,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0     ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B2.Set(-5,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0     ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    dom.SetProps(B2);

    //dom.GetParticle(-4)->vzf = true;
    //dom.GetParticle(-5)->vzf = true;

    //double vel = 0.5*str*(dom.GetParticle(-4)->x(2)-dom.GetParticle(-5)->x(2));
    //dom.GetParticle(-4)->v = Vec3_t( vel,0.0,0.0);
    //dom.GetParticle(-5)->v = Vec3_t(-vel,0.0,0.0);
    if (test=="sawtooth")
    {
        dat.L0   = dom.GetParticle(-5)->x(2) - dom.GetParticle(-4)->x(2);
        UdCu UC;
        UC.Lz    = dom.GetParticle(-5)->x(2)-dom.GetParticle(-4)->x(2)-(dom.GetParticle(-5)->MinZ()-dom.GetParticle(-4)->MaxZ());
        UC.x0    = dom.GetParticle(-5)->x(0);
        UC.str   = str;
        UC.T0    = T0;
        UC.Tf    = Tf;
        UC.sawL  = Lx/4.0;
        UC.i1    = dom.GetParticle(-4)->Index;
        UC.i2    = dom.GetParticle(-5)->Index;
        cudaMalloc(&dat.pUD, sizeof(UdCu));
        cudaMemcpy(dat.pUD, &UC, sizeof(UdCu), cudaMemcpyHostToDevice);
        cudaMalloc(&dat.pcrossed, sizeof(bool));

        dom.Solve (/*tf*/T0+Tf, /*dt*/dt, /*dtOut*/dtOut, &SetupSaw, &ReportSaw, fkey_b.CStr(),RenderVideo,Nproc);
    }
    else if (test=="normal")
    {
        dat.L0   = dom.GetParticle(-5)->x(2) - dom.GetParticle(-4)->x(2);
        UdCu UC;
        UC.Lz    = dom.GetParticle(-5)->x(2)-dom.GetParticle(-4)->x(2)-(dom.GetParticle(-5)->MinZ()-dom.GetParticle(-4)->MaxZ());
        UC.i1    = dom.GetParticle(-4)->Index;
        UC.i2    = dom.GetParticle(-5)->Index;
        UC.str   = str;
        UC.T0    = T0;
        UC.Tf    = Tf;
        UC.mtop  = 0.0;
        UC.ztop  = 0.0;
        UC.zbtop = 0.0;
        UC.Ftop  = 0.0;
        
        double maxd = dom.MaxDim();
        double Ztop = dom.GetParticle(-5)->x(2);
        double Zbot = dom.GetParticle(-4)->x(2);
        thrust::host_vector<size_t> hTopPar;
        thrust::host_vector<size_t> hBotPar;

        for (size_t i=0;i<dom.Particles.Size();i++)
        {
            if (!dom.Particles[i]->IsFree()) continue;
            if (Ztop-dom.Particles[i]->x(2)<dthoot*maxd)
            {
                dom.Particles[i]->FixFree = true;
                dom.Particles[i]->FixVeloc();
                dom.Particles[i]->InitializeVelocity(dt);
                dom.Particles[i]->Tag     = -3;
                hTopPar.push_back(i);
                dat.TPar.Push(i);
                UC.mtop += dom.Particles[i]->Props.m;
            }
            if (dom.Particles[i]->x(2)-Zbot<dthoot*maxd)
            {
                dom.Particles[i]->FixFree = true;
                dom.Particles[i]->FixVeloc();
                dom.Particles[i]->InitializeVelocity(dt);
                dom.Particles[i]->Tag     = -2;
                hBotPar.push_back(i);
                dat.BPar.Push(i);
            }
        }
        dat.TopPar  = hTopPar;
        dat.BotPar  = hBotPar;
        dat.nptop   = hTopPar.size();
        dat.npbot   = hBotPar.size();
        UC .nptop   = hTopPar.size();
        UC .npbot   = hBotPar.size();
        dat.pTopPar = thrust::raw_pointer_cast(dat.TopPar.data());
        dat.pBotPar = thrust::raw_pointer_cast(dat.BotPar.data());

        cudaMalloc(&dat.pUD, sizeof(UdCu));
        cudaMemcpy(dat.pUD, &UC, sizeof(UdCu), cudaMemcpyHostToDevice);

        dom.Solve (/*tf*/T0+Tf, /*dt*/dt, /*dtOut*/dtOut, &SetupNor, &ReportNor, fkey_b.CStr(),RenderVideo,Nproc);
    }
    dom.Save(fkey_b.CStr());

    return 0;
}
MECHSYS_CATCH
