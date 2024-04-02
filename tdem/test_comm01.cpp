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
#include <iostream>
#include <math.h>

// GSL
#include <gsl/gsl_linalg.h>

// MechSys
#include <mechsys/util/fatal.h>
#include <mechsys/dem/particle.h>
#include <mechsys/dem/domain.h>

using std::cout;
using std::endl;

#define TAG_PARTICLE 10000

int main(int argc, char **argv) try
{
#ifdef USE_MPI
    MECHSYS_CATCH_PARALLEL = true;

    // init
    MPI::Init (argc, argv);
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors

    // create particle data types
    BuildParticleDataType (MPI_Particle_Type);
    BuildParticlePropsDataType (MPI_Part_Props_Type);

    // check
    cout << "hi, I'm processor # " << my_id << endl;

    size_t verts_size = 0;

    if (my_id==0)
    {
        // gen particle
        DEM::Domain dom;
        dom.AddCube (-1, /*x*/Vec3_t(0.5,0.5,0.5), /*R*/0.1, /*L*/1.0, /*rho*/1.0);
        Particle * p = dom.GetParticle (-1);
        p->Index = 10;
        verts_size = p->Verts.Size();

        // send data
        for (int k=1; k<nprocs; ++k)
        {
            p->SendParticle(k,TAG_PARTICLE);
        }

        // draw the cube
        dom.WriteBPY("proc");
        dom.WritePOV("proc");
    }
    else if (my_id==1)
    {
        DEM::Domain dom;
        // dummy particle
        Particle  *p = new Particle();
        p->ReceiveParticle(TAG_PARTICLE);
        Vec3_t trans(-1.5,0.0,0.0);
        p->Translate(trans);
        // include the particle in the domain and draw it
        dom.Particles.Push(p);
        dom.WriteBPY("proc");
        dom.WritePOV("proc");
    }
    else
    {
        DEM::Domain dom;
        // dummy particle
        Particle  *p = new Particle();
        p->ReceiveParticle(TAG_PARTICLE);
        Vec3_t trans(-0.5,0.0,0.0);
        p->Translate(trans);


        // include the particle in the domain and draw it
        dom.Particles.Push(p);
        dom.WriteBPY("proc");
        dom.WritePOV("proc");
        //cout << "processor # " << my_id << " has got the following particle: \n" << *p;
    }

    cout << "processor # " << my_id << " has finished" << endl;

    // end
    MPI::Finalize();
#else
    throw new Fatal("This program needs MPI");
#endif
    return 0;
}
MECHSYS_CATCH
