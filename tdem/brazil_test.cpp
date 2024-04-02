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
#include <mechsys/linalg/matvec.h>
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;
using std::ofstream;
using DEM::Domain;

struct UserData
{
    DEM::Particle         * p1;  // Upper plane
    DEM::Particle         * p2;  // Lower plane
    Vec3_t                  X1;  //Initial position of first plate
    Vec3_t                  X2;  //Initial position of second plate
    Vec3_t               force;  // Force on planes
    double                  L0;  // Initial vertical separation of the plates
    double                   S;  // Vertical separation of the planes
    double                strf;  // Final strain
    double              strflu;  // strain fluctuation
    double              angle0;  // Angle of the plates
    double                 ome;  // frequency
    double                  Tf;  // final time
    std::ofstream       oss_ss;  // File to store the forces
};

void Setup (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    dat.force = 0.5*(dat.p2->F-dat.p1->F);
    dat.S     = (dat.L0 - dat.p2->x(1) + dat.p1->x(1))/dat.L0;
    Vec3_t    X2 = dat.X2 - 0.5*dat.L0*(dat.strf*dom.Time/dat.Tf +  dat.strflu*sin(dat.ome*dom.Time))*Vec3_t(sin(dat.angle0),cos(dat.angle0),0.0);
    Vec3_t    X1 = dat.X1 + 0.5*dat.L0*(dat.strf*dom.Time/dat.Tf +  dat.strflu*sin(dat.ome*dom.Time))*Vec3_t(sin(dat.angle0),cos(dat.angle0),0.0);
    dat.p2->v        = (X2 - dat.p2->x)/dom.Dt;
    dat.p1->v        = (X1 - dat.p1->x)/dom.Dt;
    dat.p2->Position(X2);
    dat.p1->Position(X1);
}

void Report (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "fx" << Util::_8s << "fy" << Util::_8s << "fz" << Util::_8s << "S" << Util::_8s << "Nc \n";
    }
    if (!dom.Finished) 
    {
        dat.force = 0.5*(dat.p2->F-dat.p1->F);
        dat.S     = (dat.L0 - dat.p2->x(1) + dat.p1->x(1))/dat.L0;
        String ff;
        ff.Printf    ("%s_bf_%04d",dom.FileKey.CStr(), dom.idx_out);
        dom.WriteBF  (ff.CStr());
        ff.Printf    ("%s_frac_%04d",dom.FileKey.CStr(), dom.idx_out);
        dom.WriteFrac(ff.CStr());
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << fabs(dat.force(0)) << Util::_8s << fabs(dat.force(1)) << Util::_8s << fabs(dat.force(2)) << Util::_8s << dat.S << Util::_8s << dom.Listofclusters.Size() << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

void CreateContact(DEM::Domain & dom, int Tag,  double angle, double angle0,size_t Ndiv, double R, double SR, double thick)
{
    Array<Vec3_t>       V(2*Ndiv+2);
    Array<Array <int> > E(3*Ndiv+1);
    Array<Array <int> > F(Ndiv);

    for (size_t i=0; i<=Ndiv;i++) 
    {
        V[i       ] = Vec3_t((R+1.1*SR)*sin(i*angle*M_PI/(180.0*Ndiv)-0.5*angle*M_PI/180.0+angle0*M_PI/180.0),(R+1.1*SR)*cos(i*angle*M_PI/(180.0*Ndiv)-0.5*angle*M_PI/180.0+angle0*M_PI/180.0),-0.5*thick);
        V[i+Ndiv+1] = Vec3_t((R+1.1*SR)*sin(i*angle*M_PI/(180.0*Ndiv)-0.5*angle*M_PI/180.0+angle0*M_PI/180.0),(R+1.1*SR)*cos(i*angle*M_PI/(180.0*Ndiv)-0.5*angle*M_PI/180.0+angle0*M_PI/180.0), 0.5*thick);
        E[i+2*Ndiv].Resize(2);
        E[i+2*Ndiv] = i,i+Ndiv+1;
    }
    for (size_t i=0; i<Ndiv;i++) 
    {
        E[i     ].Resize(2);
        E[i     ] = i  ,i+1;
        E[i+Ndiv].Resize(2);
        E[i+Ndiv] = i+Ndiv+1,i+Ndiv+2;
        F[i].Resize(4);
        F[i]      = i,i+Ndiv+1,i+Ndiv+2,i+1;
    }
    Vec3_t centr(OrthoSys::O);
    for (size_t i=0; i<V.Size();i++) 
    {
        centr += V[i];
    }
    centr/=V.Size();

    dom.Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,SR,1.0));
    dom.Particles[dom.Particles.Size()-1]->Q           = 1.0,0.0,0.0,0.0;
    dom.Particles[dom.Particles.Size()-1]->Props.V     = thick*R*angle*M_PI/180.0;
    dom.Particles[dom.Particles.Size()-1]->Props.m     = 10.0;
    dom.Particles[dom.Particles.Size()-1]->I           = 1.0,1.0,1.0;
    dom.Particles[dom.Particles.Size()-1]->x           = centr;
    dom.Particles[dom.Particles.Size()-1]->Ekin        = 0.0;
    dom.Particles[dom.Particles.Size()-1]->Erot        = 0.0;
    dom.Particles[dom.Particles.Size()-1]->Dmax        = (R+SR)*angle*M_PI/180.0;
    dom.Particles[dom.Particles.Size()-1]->PropsReady  = true;
    dom.Particles[dom.Particles.Size()-1]->Index       = dom.Particles.Size()-1;
}

int main(int argc, char **argv) try
{
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc>2) Nproc = atoi(argv[2]);

    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    UserData dat; 
    Domain d(&dat);
    Mesh::Unstructured mesh(2);
    size_t Render;
    size_t n_divisions = 30;
    double thickness   = 5.0;
    double radius      = 20.0;
    double Kn          = 1.0e6;
    double Kt          = 3.3e5;
    double Bn          = 1.0e6;
    double Bt          = 3.3e5;
    double Bm          = 3.3e5;
    double Gn          = 16.0;
    double Gt          = 8.0;
    double eps         = 0.01;
    double Amax        = 10.0;
    double Alpha       = 0.1;
    double width       = 0.5;
    double SR          = 0.1;
    double angle0      = 0.1;
    double angle       = 0.1;
    double dt          = 0.0001;
    double dtOut       = 0.1;
    double Tf          = 10.0;
    double strf        = 0.01;
    double strflu      = 0.001;
    double ome         = 1000.0;


    infile >> Render;             infile.ignore(200,'\n');
    infile >> n_divisions;        infile.ignore(200,'\n');
    infile >> thickness;          infile.ignore(200,'\n');
    infile >> radius;             infile.ignore(200,'\n');
    infile >> Kn;                 infile.ignore(200,'\n');
    infile >> Kt;                 infile.ignore(200,'\n');
    infile >> Bn;                 infile.ignore(200,'\n');
    infile >> Bt;                 infile.ignore(200,'\n');
    infile >> Bm;                 infile.ignore(200,'\n');
    infile >> Gn;                 infile.ignore(200,'\n');
    infile >> Gt;                 infile.ignore(200,'\n');
    infile >> eps;                infile.ignore(200,'\n');
    infile >> Amax;               infile.ignore(200,'\n');
    infile >> Alpha;              infile.ignore(200,'\n');
    infile >> width;              infile.ignore(200,'\n');
    infile >> SR;                 infile.ignore(200,'\n');
    infile >> angle0;             infile.ignore(200,'\n');
    infile >> angle;              infile.ignore(200,'\n');
    infile >> dt;                 infile.ignore(200,'\n');
    infile >> dtOut;              infile.ignore(200,'\n');
    infile >> Tf;                 infile.ignore(200,'\n');
    infile >> strf;               infile.ignore(200,'\n');
    infile >> strflu;             infile.ignore(200,'\n');
    infile >> ome;                infile.ignore(200,'\n');


    mesh.Set    (n_divisions, n_divisions, 1, 0);  // division, edges, regions and holes
    mesh.SetReg (0,  -1,  Amax,  0.0, 0.0);  // id, tag, max{volume}, x, y, z <<<<<<< regions
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetPnt(i  , 0, radius*cos(2*i*M_PI/n_divisions),radius*sin(2*i*M_PI/n_divisions));
    }
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetSeg(i, 0, i, (i+1)%n_divisions);
    }

    mesh.Generate();
    d.GenFromMesh(mesh,Alpha,3.0,true,false,thickness);
    d.Alpha = Alpha;
    d.Dilate = true;
    d.Center();
    Vec3_t Xmin,Xmax;
    d.BoundingBox(Xmin,Xmax);
    Vec3_t velocity;
    if (angle<1.0e-3)
    {
        d.AddPlane(-2, Vec3_t(0.0,Xmin(1)-0.5*sqrt(Amax/10),0.0), 0.5*sqrt(Amax/10), width*radius, 1.2*thickness, 1.0, M_PI/2.0, &OrthoSys::e0);
        d.AddPlane(-3, Vec3_t(0.0,Xmax(1)+0.5*sqrt(Amax/10),0.0), 0.5*sqrt(Amax/10), width*radius, 1.2*thickness, 1.0, M_PI/2.0, &OrthoSys::e0);
        velocity = Vec3_t(0.0,strf*radius/Tf,0.0);
    }
    else
    {
        CreateContact(d,-2,angle,angle0      ,n_divisions*angle/360.0,1.0*radius,SR,thickness);
        CreateContact(d,-3,angle,angle0+180.0,n_divisions*angle/360.0,1.0*radius,SR,thickness);
        //CreateContact(d,-2,angle,angle0      ,20.0*angle/60.0,radius,SR,thickness);
        //CreateContact(d,-3,angle,angle0+180.0,20.0*angle/60.0,radius,SR,thickness);
        velocity = -strf*radius/Tf*Vec3_t(sin(angle0),cos(angle0),0.0);
    }
    DEM::Particle * p1 = d.GetParticle(-2);
    DEM::Particle * p2 = d.GetParticle(-3);
    p1->FixVeloc();
    p1->v =  velocity;
    p2->FixVeloc();
    p2->v = -velocity;

    dat.p1     = p1;
    dat.p2     = p2;
    dat.X1     = p1->x;
    dat.X2     = p2->x;
    dat.L0     = p2->x(1) - p1->x(1);
    dat.Tf     = Tf;
    dat.ome    = ome*2*M_PI/Tf;
    dat.strf   = strf;
    dat.strflu = strflu;
    dat.angle0 = angle0;
    dat.force = OrthoSys::O;
    dat.S     = 0.0;
    
    // properties of particles prior the brazilian test
    Dict B;
    B.Set(-1,"Bn Bt Bm Gn Gt Eps Kn Kt"   ,Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt);
    B.Set(-2,"Bn Bt Bm Gn Gt Eps Kn Kt Mu",Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt,0.0);
    B.Set(-3,"Bn Bt Bm Gn Gt Eps Kn Kt Mu",Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt,0.0);
    d.SetProps(B);

    //d.Solve(Tf, dt, dtOut, &Setup, &Report, filekey.CStr(),Render,Nproc);
    d.Solve(Tf, dt, dtOut, NULL, &Report, filekey.CStr(),Render,Nproc);


    return 0;
}
MECHSYS_CATCH
