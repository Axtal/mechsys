/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

// Std Lib
#include <cstdlib> // for srand, rand and RAND_MAX
#include <ctime>   // for time

// VTK
#include <vtkPolygon.h>
#include <vtkLightKit.h>

// MechSys
#include <mechsys/dem/visualise.h>
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/vtk/win.h>
#include <mechsys/vtk/axes.h>

using Util::PI;
using namespace DEM;

struct Data
{
     Data () : Vis(NULL) {}
    ~Data () { if (Vis!=NULL) delete Vis;  if (Lgt!=NULL) Lgt->Delete(); }
    void Init (DEM::Domain const & Dom)
    {
        Vis = new Visualise(Dom, /*parts*/Array<int>(-10,true), /*walls*/Array<int>(-1,-2,-3,-4,-5,-6));
        //VTK::Axes axe(1.0, /*hline*/false, /*reverse*/false, /*full*/true);
        Vis->AddTo  (Win);
        //axe .AddTo  (Win);

        //Vis->PartFaces.SetMaterial (0.2, 0.1, 1.0, 200.0);
        //Vis->PartFaces.SetColor    ("cobalt", 1.0);
        
        Win .Render ();
        Win .Render ();
        
        Lgt = vtkLightKit::New();
        Lgt->AddLightsToRenderer (Win.GetRen());


        //Lgt->SetKeyLightAngle (/*elevation*/60., /*azimuth*/-45.0);
        //Lgt->SetKeyToHeadRatio (0.01);
        Lgt->SetKeyToFillRatio (1.2);

        IdxOut = 0;
    }
    Visualise   * Vis;
    VTK::Win      Win;
    int           IdxOut;
    vtkLightKit * Lgt;
};

void Report (DEM::Domain & dom, void * UserData)
{
    Data & dat = (*static_cast<Data*>(UserData));
    String buf;
    buf.Printf ("genpack_%08d.png", dat.IdxOut);
    dat.Vis->Update();
    dat.Win.WritePNG (buf.CStr());
    dat.IdxOut++;
}

int main(int argc, char **argv) try
{
    // input
    int    num   = 4;      // n per side (with 5 spheres, there is overlapping)
    bool   voro  = false;  // Voronoi particles
    bool   vtk   = true;   // show VTK window
    bool   ascii = true;   // write ASCII file instead of HDF5 (spheres only)
    bool   sim   = true;   // run simulation
    double tf    = 20.0;   // final time
    double dt    = 1.0e-4; // time step
    double dtOut = 1.0;    // output time step
    bool   HCP   = true;   // HCP packing
    if (argc>1) num   = atoi(argv[1]);
    if (argc>2) voro  = atoi(argv[2]);
    if (argc>3) vtk   = atoi(argv[3]);
    if (argc>4) ascii = atoi(argv[4]);
    if (argc>5) sim   = atoi(argv[5]);
    if (argc>6) tf    = atof(argv[6]);
    if (argc>7) dt    = atof(argv[7]);
    if (argc>8) dtOut = atof(argv[8]);
    if (argc>9) HCP   = atoi(argv[9]);
    if (ascii && voro) throw new Fatal("ASCII file works only with spheres");

    // user data and domain
    Data dat;
    DEM::Domain dom((vtk ? &dat : NULL));

    // generate particles
    double R        = 0.05;
    double L        = 22.0;//1.0;
    size_t seed     = 123;
    double fraction = 1.0;
    //if (voro) dom.AddVoroPack (-10, R, L,L,L, num,num,num, /*rho*/1.0, /*cohesion*/false,/*periodic*/true, seed, fraction);
    //if (voro) dom.AddVoroPack (-10, R, L,L,L, num,num,4*num, /*rho*/1.0, /*cohesion*/false,/*periodic*/true, seed, fraction, Vec3_t(0.8,0.8,0.8));
    //if (voro) dom.AddVoroPack (-10, R, L,L,L, num,num,10*num, /*rho*/1.0, /*cohesion*/false,/*periodic*/true, seed, fraction, Vec3_t(1.,1.,1.));
    if (voro) dom.AddVoroPack (-10, R, L,L,L, num,num,4.0*num, /*rho*/1.0, /*cohesion*/false,/*periodic*/true, seed, fraction);
    else
    {
        if (HCP) dom.GenSpheres  (-10, L, num, /*rho*/1.0, "HCP",    /*seed*/1000, /*fraction*/1.0);
        else     dom.GenSpheres  (-10, L, num, /*rho*/1.0, "Normal", /*seed*/1000, /*fraction*/1.0);
    }

    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    std::cout << Xmin << Xmax << std::endl;
    // run simulation
    if (sim)
    {
        // generate walls
        double c = 1.0;
        dom.Center   ();
        dom.AddPlane (-1,Vec3_t(  -c, 0.0, 0.0), /*R*/0.2, /*lx*/2.5*c, /*ly*/2.5*c, /*rho*/1.0, PI/2.0, &OrthoSys::e1); // -x face
        dom.AddPlane (-2,Vec3_t(   c, 0.0, 0.0), /*R*/0.2, /*lx*/2.5*c, /*ly*/2.5*c, /*rho*/1.0, PI/2.0, &OrthoSys::e1); // +x face
        dom.AddPlane (-3,Vec3_t( 0.0,  -c, 0.0), /*R*/0.2, /*lx*/2.5*c, /*ly*/2.5*c, /*rho*/1.0, PI/2.0, &OrthoSys::e0); // -y face
        dom.AddPlane (-4,Vec3_t( 0.0,   c, 0.0), /*R*/0.2, /*lx*/2.5*c, /*ly*/2.5*c, /*rho*/1.0, PI/2.0, &OrthoSys::e0); // +y face
        dom.AddPlane (-5,Vec3_t( 0.0, 0.0,  -c), /*R*/0.2, /*lx*/2.5*c, /*ly*/2.5*c, /*rho*/1.0);                        // -z face
        dom.AddPlane (-6,Vec3_t( 0.0, 0.0,   c), /*R*/0.2, /*lx*/2.5*c, /*ly*/2.5*c, /*rho*/1.0);                        // +z face

        // set properties of particles
        double Kn  = 100.0;
        double Kt  = 1.0;
        double Gn  = 1.0;
        double Gt  = 1.0;
        double mu  = 0.0;
        double bet = 0.0;
        double eta = 0.0;
        Dict prps;
        prps.Set ( -1, "Kn Kt Gn Gt Mu Beta Eta", Kn, Kt, Gn, Gt, mu, bet, eta);
        prps.Set ( -2, "Kn Kt Gn Gt Mu Beta Eta", Kn, Kt, Gn, Gt, mu, bet, eta);
        prps.Set ( -3, "Kn Kt Gn Gt Mu Beta Eta", Kn, Kt, Gn, Gt, mu, bet, eta);
        prps.Set ( -4, "Kn Kt Gn Gt Mu Beta Eta", Kn, Kt, Gn, Gt, mu, bet, eta);
        prps.Set ( -5, "Kn Kt Gn Gt Mu Beta Eta", Kn, Kt, Gn, Gt, mu, bet, eta);
        prps.Set ( -6, "Kn Kt Gn Gt Mu Beta Eta", Kn, Kt, Gn, Gt, mu, bet, eta);
        prps.Set (-10, "Kn Kt Gn Gt Mu Beta Eta", Kn, Kt, Gn, Gt, mu, bet, eta);
        dom.SetProps (prps);

        // initialise
        dom.Initialize(dt);

        // initialise window (must be after dom.Initialize())
        if (vtk)
        {
            dat.Init (dom);
            dat.Win.Show ();
        }

        // set random velocity to particles
        srand (time(NULL));
        for (size_t i=0; i<dom.Particles.Size(); ++i)
        {
            if (dom.Particles[i]->Tag==-10)
            {
                dom.Particles[i]->Ff = dom.Particles[i]->Props.m*Vec3_t(0.0,0.0,-0.01);
            }
        }

        // fix walls
        dom.GetParticle(-1)->FixVeloc();
        dom.GetParticle(-2)->FixVeloc();
        dom.GetParticle(-3)->FixVeloc();
        dom.GetParticle(-4)->FixVeloc();
        dom.GetParticle(-5)->FixVeloc();
        dom.GetParticle(-6)->FixVeloc();

        // solve (fall)
        dom.ResetInteractons (); // this is necessary if FixVeloc is applied after dom.Initialize
        dom.Solve (tf, dt, dtOut, /*setup*/NULL, (vtk ? &Report : NULL));

        // visualise
        if (vtk) dat.Win.Show ();

        // move walls
        double v = 0.001;
        dom.GetParticle(-1)->FixVeloc(  v, 0.0, 0.0);
        dom.GetParticle(-2)->FixVeloc( -v, 0.0, 0.0);
        dom.GetParticle(-3)->FixVeloc(0.0,   v, 0.0);
        dom.GetParticle(-4)->FixVeloc(0.0,  -v, 0.0);
        dom.GetParticle(-5)->FixVeloc(0.0, 0.0,   v);
        dom.GetParticle(-6)->FixVeloc(0.0, 0.0,  -v);

        // solve (isotropic compression)
        dom.ResetInteractons (); // this is necessary if FixVeloc is applied after dom.Initialize
        try
        {
            dom.Solve (10*tf, dt, dtOut, /*setup*/NULL, (vtk ? &Report : NULL));
        }
        catch (Fatal *)
        {
            std::cout << "\n>>> Maximum compression achieved <<<\n" << std::endl;
        }

        // visualise
        if (vtk) dat.Win.Show ();

        // delete walls
        dom.DelParticles (Array<int>(-1,-2,-3,-4,-5,-6));
    }

    // only particles
    else
    {
        dom.Initialize(dt);

        // TODO: delete this
        /*
        for (size_t i=0; i<dom.Particles.Size(); ++i)
        {
            double vx = 0.05*static_cast<double>(1.0*rand()/RAND_MAX);
            double vy = 0.05*static_cast<double>(1.0*rand()/RAND_MAX);
            double vz = 0.05*static_cast<double>(1.0*rand()/RAND_MAX);
            dom.Particles[i]->v = vx, vy, vz;
        }
        */

        if (vtk)
        {
            dat.Init         (dom);
            dat.Win.Show     ();
            dat.Win.WritePNG ("dem_genpack");
        }
    }

    // save file
    cout << "Number of particles = " << dom.Particles.Size() << endl;
    if (ascii)
    {
        // header
        Array<String> keys("id", "xc", "yc", "zc", "ra", "vx", "vy", "vz", "ct");
        std::ostringstream oss;
        oss << Util::_6 << keys[0];
        for (size_t i=1; i<keys.Size()-1; ++i) { oss << Util::_8s << keys[i]; }
        oss << Util::_6 << keys.Last() << "\n";

        // values
        for (size_t i=0; i<dom.Particles.Size(); ++i)
        {
            oss << Util::_6  << i;
            oss << Util::_8s << dom.Particles[i]->x(0);
            oss << Util::_8s << dom.Particles[i]->x(1);
            oss << Util::_8s << dom.Particles[i]->x(2);
            oss << Util::_8s << dom.Particles[i]->Props.R;
            oss << Util::_8s << dom.Particles[i]->v(0);
            oss << Util::_8s << dom.Particles[i]->v(1);
            oss << Util::_8s << dom.Particles[i]->v(2);
            oss << Util::_6  << -1;
            oss << "\n";
        }

        // open file and save data
        String buf("xpack-spheres.dat");
        std::ofstream of(buf.CStr(), std::ios::out);
        of << oss.str();
        of.close();
        cout << "file [1;34m<"<<buf<<">[0m written" << endl;
    }
    else
    {
        // save HDF5
        dom.Save("xpack-voro");
        cout << "file [1;34m<xpack-voro.hdf5>[0m written" << endl;
    }

    // end
    return 0;
}
MECHSYS_CATCH
