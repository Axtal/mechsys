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


#ifndef MECHSYS_LBM_DEM_H
#define MECHSYS_LBM_DEM_H

// STD
#include <algorithm>

// Mechsys
#include <mechsys/lbm/Lattice.h>
#include <mechsys/linalg/quaternion.h>

namespace LBM
{

class Disk
{
public:
    Disk(int Tag, Vec3_t const & X, Vec3_t const & V, Vec3_t const & W, double rho, double R, double dt);

    //Methods
    void Translate   (double dt);
    void Rotate      (double dt);
    void FixVeloc    () {vf=true,true,true; wf=true,true,true;};
    bool IsFree () {return !vf(0)&&!vf(1)&&!vf(2)&&!wf(0)&&!wf(1)&&!wf(2);}; ///< Ask if the particle has any constrain in its movement

#ifdef USE_THREAD
    pthread_mutex_t lck;   ///< Lock to protect variables from race conditions.
#elif USE_OMP
    omp_lock_t      lck;             ///< to protect variables in multithreading
#endif
    // Data
    int  Tag;              ///< Id of the particle
    Vec3_t X0;             ///< initial position of the particle
    Vec3_t X;              ///< position of the particle
    Vec3_t Xb;             ///< position of the particle before
    Vec3_t V;              ///< velocity of the CM
    Vec3_t W;              ///< angular velocity
    Vec3_t Wb;             ///< angular velocity
    Vec3_t F;              ///< Force
    Vec3_t Flbm;           ///< Force exceted by LBM fluid (hydraulic force)
    Vec3_t Ff;             ///< fixed Force
    Vec3_t T;              ///< Torque
    Vec3_t Tf;             ///< fixed Torque
    Quaternion_t Q;        ///< The quaternion representing the rotation
    double R;              ///< Disk radius
    double M;              ///< mass of the disk
    double I;              ///< inertia moment of the particle
    double Gn;             ///< dissipation constant for collision
    double Gt;             ///< dissipation constant for collision
    double Kn;             ///< Stiffness constant
    double Kt;             ///< Tangential stiffness constant
    double Mu;             ///< Friction coefficient
    double Beta;           ///< Rolling stiffness coeffcient
    double Eta;            ///< Plastic moment coefficient
    bVec3_t vf;            ///< prescribed velocities
    bVec3_t wf;            ///< prescribed angular velocities
    
};

inline Disk::Disk(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt)
{
    Tag = TheTag;
    X   = TheX;
    V   = TheV;
    W   = TheW;
    R   = TheR;
    M   = M_PI*R*R*Therho;
    I   = 0.5*M*R*R;
    Xb  = X - dt*V;
    Wb  = W;
    F   = 0.0,0.0,0.0;
    Flbm= 0.0,0.0,0.0;
    T   = 0.0,0.0,0.0;
    Ff  = 0.0,0.0,0.0;
    Tf  = 0.0,0.0,0.0;
    vf  = false,false,false;
    wf  = false,false,false;
    Gn  = 8.0;
    Gt  = 0.0;
    Kn  = 1.0e3;
    Kt  = 5.0e2;
    Mu  = 0.4;
    Eta = 1.0;  
    Beta = 0.12; 
    Q    = 1.0,0.0,0.0,0.0;

#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#elif USE_OMP
    omp_init_lock(&lck);
#endif
}

inline void Disk::Translate(double dt)
{
    //std::cout << F(0) << " " << M << " " << V(0) << std::endl;
    Vec3_t Ft = F;
    if (vf(0)) Ft(0) = 0.0;
    if (vf(1)) Ft(1) = 0.0;
    if (vf(2)) Ft(2) = 0.0;
    //if (isnan(norm(F))) 
    //{
        //std::cout << Tag << std::endl;
    //}

    Vec3_t Xa = 2*X - Xb + Ft*(dt*dt/M);
    Vec3_t tp = Xa - X;
    V         = 0.5*(Xa - Xb)/dt;
    Xb        = X;
    X         = Xa;
}

inline void Disk::Rotate (double dt)
{
    double q0,q1,q2,q3,wx,wy,wz;
    q0 = 0.5*Q(0);
    q1 = 0.5*Q(1);
    q2 = 0.5*Q(2);
    q3 = 0.5*Q(3);

    Vec3_t Tt = T;
    if (wf(0)) T(0) = 0.0;
    if (wf(1)) T(1) = 0.0;
    if (wf(2)) T(2) = 0.0;

    Vec3_t Td = Tt/I;
    W = Wb+0.5*dt*Td;
    wx = W(0);
    wy = W(1);
    wz = W(2);
    Quaternion_t dq(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz),qm;

    Wb  = Wb+Td*dt;
    qm  = Q+dq*(0.5*dt);
    q0  = 0.5*qm(0);
    q1  = 0.5*qm(1);
    q2  = 0.5*qm(2);
    q3  = 0.5*qm(3);
    wx  = Wb(0);
    wy  = Wb(1);
    wz  = Wb(2);
    dq  = Quaternion_t(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz);
    Quaternion_t Qd = (qm+dq*0.5*dt),temp;
    Q  = Qd/norm(Qd);
}

}
#endif
