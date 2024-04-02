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

// Std lib
#include <math.h>

// GSL
#include <gsl/gsl_linalg.h>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;

struct UserData
{
    DEM::Particle *    p;            // the array of particles at which the force is to be applied
    Array<Vec3_t  >    vm0;          // value of the vectors close to the middle section
    Array<Vec3_t *>    vm;           // pointers to the vectors close to the middle section
    String             test;         // Type of test vibraiton or tension
    double             A;            // Area of the plate for stress calculation
    double             Am;           // vibration amplitude
    double             ome;          // vibration frequency
    double             sy;           // Stress State
    double             Tf;           // Final time
    Vec3_t             L0;           // Initial dimensions
    std::ofstream      oss_ss;       // file for stress strain data
};

void Setup (DEM::Domain & Dom, void * UD)
{
    // force at particle tagged -3
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.test=="vibration")    dat.p->Ff=0.0,0.0,dat.Am*sin(dat.ome*Dom.Time);
    if (dat.test=="bending")
    {
        if (Dom.Time < 0.5*dat.Tf) dat.p->Ff=0.0,0.0,dat.Am*2*Dom.Time/dat.Tf;
        else                       dat.p->Ff=0.0,0.0,dat.Am;
    }
    dat.sy = dat.p->F(1)/dat.A;
}

void Report (DEM::Domain & Dom, void * UD)
{ 
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.test=="tensile")
    {
        if (Dom.idx_out==0)
        {
            String fs;
            fs.Printf("%s_walls.res",Dom.FileKey.CStr());
            dat.oss_ss.open(fs.CStr());
            dat.sy = dat.p->F(1)/dat.A;
            // Output of the current time, the stress state sx, and the strains ex,ey and ez
            dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sy" << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez" << std::endl;
        }
        if (!Dom.Finished)
        {
            // Measure the strains and stresses
            Vec3_t Xmin, Xmax;
            Dom.BoundingBox(Xmin, Xmax);
            double ex = (Xmax(0)-Xmin(0)-dat.L0(0))/dat.L0(0);
            double ey = (Xmax(1)-Xmin(1)-dat.L0(1))/dat.L0(1);
            double ez = (Xmax(2)-Xmin(2)-dat.L0(2))/dat.L0(2);
            dat.oss_ss << Util::_10_6 << Dom.Time << Util::_8s << dat.sy << Util::_8s << ex << Util::_8s << ey << Util::_8s << ez << std::endl;
        }
        else dat.oss_ss.close();
    }
    if (dat.test=="bending")
    {
        if (Dom.idx_out==0)
        {
            double tol = dat.L0(2)/20.0;
            for (size_t i=0;i<Dom.Particles.Size();i++)
            {
                for (size_t j=0;j<Dom.Particles[i]->Verts.Size();j++)
                {
                    if (fabs((*Dom.Particles[i]->Verts[j])(2))<tol)
                    {
                        dat.vm.Push (Dom.Particles[i]->Verts[j]);
                        dat.vm0.Push(Vec3_t(*Dom.Particles[i]->Verts[j]));
                    }
                }
            }
        }
        String fs;
        fs.Printf("%s_%08d.res",Dom.FileKey.CStr(),Dom.idx_out);
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s <<"ux" << Util::_8s << "uy" << Util::_8s << "uz" << std::endl;
        for (size_t i=0;i<dat.vm.Size();i++)
        {
            dat.oss_ss << Util::_8s <<            dat.vm0[i](0) << Util::_8s <<            dat.vm0[i](1) << Util::_8s <<            dat.vm0[i](2);
            dat.oss_ss << Util::_8s << (*dat.vm[i])(0)-dat.vm0[i](0) << Util::_8s << (*dat.vm[i])(1)-dat.vm0[i](1) << Util::_8s << (*dat.vm[i])(2)-dat.vm0[i](2) << std::endl;
        }
        dat.oss_ss.close();

    }
}

int main(int argc, char **argv) try
{

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    double verlet;      // Verlet distance for optimization
    String ptype;       // Particle type 
    String test;       // Particle type 
    size_t RenderVideo; // Decide is video should be render
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Bn;          // Cohesion normal stiffness
    double Bt;          // Cohesion tangential stiffness
    double Bm;          // Cohesion torque stiffness
    double Eps;         // Threshold for breking bonds
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt;          // Time step
    double dtOut;       // Time step for output
    double Tf;          // Final time for the test
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double rho;         // rho
    double Am;          // vibration force amplitude
    double ome;         // Frequency of vibration
    double ex;          // Final strain for the tensile test (positive extension, negative compression)
    {
        infile >> verlet;       infile.ignore(200,'\n');
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> test;         infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu;           infile.ignore(200,'\n');
        infile >> Bn;           infile.ignore(200,'\n');
        infile >> Bt;           infile.ignore(200,'\n');
        infile >> Bm;           infile.ignore(200,'\n');
        infile >> Eps;          infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> seed;         infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> dtOut;        infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
        infile >> Lx;           infile.ignore(200,'\n');
        infile >> Ly;           infile.ignore(200,'\n');
        infile >> Lz;           infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
        infile >> Am;           infile.ignore(200,'\n');
        infile >> ome;          infile.ignore(200,'\n');
        infile >> ex;           infile.ignore(200,'\n');
    }



    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha = verlet;
    dom.Dilate= true;

    if (ptype=="voronoi") dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, true, false, seed, 1.0);
    else if (ptype=="cube")
    {
        Mesh::Structured mesh(3);
        mesh.GenBox (false, nx, ny, nz, Lx, Ly, Lz);
        dom.GenFromMesh (mesh, R, rho, true, false);
    }
    else if (ptype=="tetra")
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        dom.GenFromMesh (mesh, R, rho, true, false);
    }
    else throw new Fatal("Packing for particle type not implemented yet");

    // Initialize the UserData structure
    dat.test = test;
    dat.A    = Lx*Ly;
    dat.Am   = Am;
    dat.ome  = ome;
    dat.Tf   = Tf;
    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    dat.L0   = Xmax - Xmin;

    dom.GenBoundingPlane(-2,-1,R,1.0,true);

    //identify the moving lid
    dat.p = dom.GetParticle (-3);
    if (test=="tensile")
    {
        dat.p->FixVeloc();
        dat.p->v = 0.0, ex*dat.L0(1)/Tf, 0.0;
    }

    //set the element properties
    Dict B;
    B.Set(-1,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(-2,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(-3,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    dom.SetProps(B);

    // fix -2 particles at the left extreme of the beam
    DEM::Particle * p;
    p = dom.GetParticle (-2);
    p->FixVeloc();

    dom.Solve (Tf,dt,dtOut, &Setup, &Report, filekey.CStr(), RenderVideo, Nproc);
}
MECHSYS_CATCH
