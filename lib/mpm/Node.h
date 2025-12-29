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
    vxf = false;
    vyf = false;
    vzf = false;
    fixCor = false;
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

#ifdef USE_CUDA
struct NodeCU
{
    bool       valid;                       ///< mark is a node is valid
    int3       X;                           ///< index position of the node
    size_t     Idx;                         ///< Index of the node
    real       Mass;                        ///< Mass of each node
    real3      Vn;                          ///< Velocity at each node
    real3      Vnf;                         ///< Fixed Velocity at each node for Dirichlet BCs
    real3      Mn;                          ///< Momentum at each node
    real3      Fn;                          ///< Force at each node
    Mat3       Sn;                          ///< Stress at the node
    bool       vxf, vyf, vzf;               ///< Fixed components of velocity
    bool       fixCor;


    __host__ __device__ void Reset()
    {
        fixCor = false;
        valid  = false; 
        Mass   = 0.0;
        Vn     = make_real3(0.0,0.0,0.0);
        Mn     = Vn;
        Fn     = Vn;
        Sn.SetZero();
    }

    __host__ __device__ void UpdateNode(real gn, real dt, real mmin)
    {
        real3 fdamp = gn*Mn;
        Fn = Fn - fdamp;
        Mn = Mn - dt*fdamp;
        if (Mass<mmin)
        {
            Vn = make_real3(0.0,0.0,0.0);
        }
        else
        {
            Vn = (1.0/Mass)*Mn;
        }
    }

    __host__ __device__ void FixNode()
    {
        if (vxf)
        {
            Fn.x = 0.0;
            Vn.x =     Vnf.x;
            Mn.x = Mass*Vn.x;
        }
        if (vyf)
        {
            Fn.y = 0.0;
            Vn.y =     Vnf.y;
            Mn.y = Mass*Vn.y;
        }
        if (vzf)
        {
            Fn.z = 0.0;
            Vn.z =     Vnf.z;
            Mn.z = Mass*Vn.z;
        }
    }
};

__host__ void UploadNode(      NodeCU & NCu, const Node & Node)
{
    NCu.valid  = Node.valid  ;
    NCu.X.x    = Node.X(0)   ;
    NCu.X.y    = Node.X(1)   ;
    NCu.X.z    = Node.X(2)   ;
    NCu.Idx    = Node.Idx    ;
    NCu.Mass   = Node.Mass   ;
    NCu.Vn.x   = Node.Vn(0)  ;
    NCu.Vn.y   = Node.Vn(1)  ;
    NCu.Vn.z   = Node.Vn(2)  ;
    NCu.Vnf.x  = Node.Vnf(0) ;
    NCu.Vnf.y  = Node.Vnf(1) ;
    NCu.Vnf.z  = Node.Vnf(2) ;
    NCu.Mn.x   = Node.Mn(0)  ;
    NCu.Mn.y   = Node.Mn(1)  ;
    NCu.Mn.z   = Node.Mn(2)  ;
    NCu.Fn.x   = Node.Fn(0)  ;
    NCu.Fn.y   = Node.Fn(1)  ;
    NCu.Fn.z   = Node.Fn(2)  ;
    NCu.Sn     = Node.Sn     ;
    NCu.vxf    = Node.vxf    ;
    NCu.vyf    = Node.vyf    ;
    NCu.vzf    = Node.vzf    ;
    NCu.fixCor = Node.fixCor ; 
}

__host__ void DnloadNode(const NodeCU & NCu,       Node & Node)
{
    Node.valid   = NCu.valid ;
    Node.X(0)    = NCu.X.x   ;
    Node.X(1)    = NCu.X.y   ;
    Node.X(2)    = NCu.X.z   ;
    Node.Mass    = NCu.Mass  ;
    Node.Vn(0)   = NCu.Vn.x  ;
    Node.Vn(1)   = NCu.Vn.y  ;
    Node.Vn(2)   = NCu.Vn.z  ;
    Node.Vnf(0)  = NCu.Vnf.x ;
    Node.Vnf(1)  = NCu.Vnf.y ;
    Node.Vnf(2)  = NCu.Vnf.z ;
    Node.Mn(0)   = NCu.Mn.x  ;
    Node.Mn(1)   = NCu.Mn.y  ;
    Node.Mn(2)   = NCu.Mn.z  ;
    Node.Fn(0)   = NCu.Fn.x  ;
    Node.Fn(1)   = NCu.Fn.y  ;
    Node.Fn(2)   = NCu.Fn.z  ;
    copyMat3(Node.Sn,NCu.Sn) ;
}

#endif
}
#endif //MECHSYS_MPM_NODE_H
