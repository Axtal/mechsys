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


#ifndef MECHSYS_MPM_NODE_H
#define MECHSYS_MPM_NODE_H

//STD
#include<iostream>

// Mechsys
#include <mechsys/mpm/Particle.h>

namespace MPM
{

class Node
{
public:
    //Constructors
    Node();
    Node(iVec3_t Pt,        ///< coordinates of the node
         size_t idx         ///< Index of the node
        );    

    //Data
    bool       valid;                       ///< flag to see if the node is valid for calculation
    iVec3_t    X;                           ///< index position of the node
    size_t     Idx;                         ///< Index of the node
    double     Mass;                        ///< Mass of each node
    Vec3_t     Vn;                          ///< Velocity at each node
    Vec3_t     Vnf;                         ///< Fixed Velocity at each node for Dirichlet BCs
    Vec3_t     Mn;                          ///< Momentum at each node
    Vec3_t     Fn;                          ///< Force at each node
    Mat3_t     Sn;                          ///< Stress at each node
    bool       vxf, vyf, vzf;               ///< Fixed components of velocity
    bool       fixCor;

    #ifdef USE_OMP
    omp_lock_t      lck;    ///< to protect variables in multithreading
    #endif

    //Methods
    void FixVeloc           (double vx=0.0, double vy=0.0, double vz=0.0);  ///< Fix all velocities
    void Reset();                                                           ///< Reset all relevant quantities to zero
    void UpdateNode(double Gn, double dt, double mmin);                     ///< Update node quantities with dissipation constant, time step and the minimun mass
    void FixNode();                                                         ///< Fix Dirichlet nodes
};

inline Node::Node()
{
    X    = OrthoSys::O;
    Idx  = 0;
    Mass = 0.0;
    Vn = OrthoSys::O;
    Vnf= OrthoSys::O;
    Mn = OrthoSys::O;
    Fn = OrthoSys::O;
    vxf = false;
    vyf = false;
    vzf = false;
    fixCor = false;

    #ifdef USE_OMP
    omp_init_lock(&lck);    ///< to protect variables in multithreading
    #endif
}

inline Node::Node(iVec3_t Pt0, size_t idx0)
{
    X = Pt0;
    Idx = idx0;
    Mass = 0.0;
    Vn = OrthoSys::O;
    Mn = OrthoSys::O;
    Fn = OrthoSys::O;
    #ifdef USE_OMP
    omp_init_lock(&lck);    ///< to protect variables in multithreading
    #endif
}

inline void Node::Reset()
{
    fixCor = false;
    valid = false;
    Mass  = 0.0;
    set_to_zero(Vn);
    set_to_zero(Mn);
    set_to_zero(Fn);
    set_to_zero(Sn);
}

inline void Node::FixVeloc (double vx, double vy, double vz)
{
    Vn  = vx, vy, vz;
    Vnf = vx, vy, vz;
    vxf = true; vyf = true; vzf = true; 
}

inline void Node::UpdateNode (double gn, double dt, double mmin)
{
    Vec3_t fdamp = gn*Mn;
    Fn -= fdamp;
    Mn -= fdamp*dt;
    //Vn = Mn/Mass;
    if (Mass<mmin)
    {
        Vn = OrthoSys::O;
    }
    else
    {
        //Vec3_t fdamp = OrthoSys::O;
        //if (norm(Mn)>1.0e-12) fdamp = gn*norm(Fn)*Mn/norm(Mn);
        Vn = Mn/Mass;
        //Sn = (1.0/Mass)*Sn;
    }
}

inline void Node::FixNode ()
{
    if (vxf)
    {
        Fn(0) = 0.0;
        Vn(0) =     Vnf(0);
        Mn(0) = Mass*Vn(0);
    }
    if (vyf)
    {
        Fn(1) = 0.0;
        Vn(1) =     Vnf(1);
        Mn(1) = Mass*Vn(1);
    }
    if (vzf)
    {
        Fn(2) = 0.0;
        Vn(2) =     Vnf(2);
        Mn(2) = Mass*Vn(2);
    }
}
}
#endif //MECHSYS_MPM_NODE_H
