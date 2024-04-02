/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang                                         *
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


#ifndef MECHSYS_MPM_CORNER_H
#define MECHSYS_MPM_CORNER_H

//STD
#include<iostream>

// Mechsys
#include <mechsys/linalg/matvec.h>

namespace MPM
{
class Corner //Class for corner points
{
public:
    //Constructor
    Corner  (Vec3_t              * x0        ///< Initial position
            );    

    //Data
    Vec3_t * x;             ///< Position of particle
    Vec3_t xb;              ///< Previous position of particle
    Vec3_t v;               ///< Velocity of the particle
    Vec3_t vf;              ///< Fixed velocity for dirichlet BCs
    Vec3_t h;               ///< Hydraulic Force
    bool   vxf, vyf, vzf;   ///< Fixed components of velocity
    size_t Nodes[8];        ///< Array of node indexes relevant to this particle
    double Sfunc[8];        ///< Relevant values of the shape function

    //Data for solid properties
    #ifdef USE_OMP
    omp_lock_t      lck;    ///< to protect variables in multithreading
    #endif
    
    //Methods
    void   FixVeloc           (double vx=0.0, double vy=0.0, double vz=0.0);  ///< Fix all velocities
    bool   IsFree             () {return !vxf&&!vyf&&!vzf;};                  ///< Ask if the particle has any constrain in its movement
    void   Reset              (double dt);                                    ///< Reset relevant quantites to zero 
    
};

inline Corner::Corner (Vec3_t * xt0)
{
    x   = xt0;
    xb  = *x;
    set_to_zero(v);
    set_to_zero(vf);
    set_to_zero(h);
    vxf = false;
    vyf = false;
    vzf = false;

    #ifdef USE_OMP
    omp_init_lock(&lck);    ///< to protect variables in multithreading
    #endif
}

inline void Corner::FixVeloc (double vx, double vy, double vz)
{
    v   = vx, vy, vz;
    vf  = vx, vy, vz;
    vxf = true; vyf = true; vzf = true; 
}

inline void Corner::Reset(double dt)
{
    //set_to_zero(h);
    v = (*x-xb)/dt;
    if (vxf)
    {
        v(0) = vf(0);
        (*x)(0) = xb(0) + v(0)*dt;
        //xb(0) = (*x)(0);
    }
    if (vyf)
    {
        v(1) = vf(1);
        (*x)(1) = xb(1) + v(1)*dt;
        //xb(1) = (*x)(1);
    }
    if (vzf)
    {
        v(2) = vf(2);
        (*x)(2) = xb(2) + v(2)*dt;
        //xb(2) = (*x)(2);
    }
    xb = *x;
}
}
#endif //MECHSYS_MPM_CORNER_H
