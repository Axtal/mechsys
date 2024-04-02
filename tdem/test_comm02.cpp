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

int main(int argc, char **argv) try
{
#ifdef USE_MPI
    MECHSYS_CATCH_PARALLEL = true;

    // init
    MPI::Init (argc, argv);
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    //int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors

    // create particle data type
    BuildParticleDataType (MPI_Particle_Type);
    BuildParticlePropsDataType (MPI_Part_Props_Type);
    
    DEM::Domain dom;
    dom.Load("test_dynamics");
    String fn;
    fn.Printf("proc_%d",my_id);
    dom.Save(fn.CStr());

    // Running the simulation
    double dt = 1.0e-5;
    dom.CamPos = 0.0,30.0,0.0;
    dom.Solve(/*tf*/30.0, dt, /*dtOut*/0.5, NULL, NULL, "test_comm02");

    // end
    MPI::Finalize();
#else
    throw new Fatal("This program needs MPI");
#endif
    return 0;
}
MECHSYS_CATCH
