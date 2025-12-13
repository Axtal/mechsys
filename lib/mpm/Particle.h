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


#ifndef MECHSYS_MPM_PARTICLE_H
#define MECHSYS_MPM_PARTICLE_H

//STD
#include<iostream>

// Mechsys
#include <mechsys/linalg/matvec.h>

namespace MPM
{
class Particle
{
public:
    //Constructor
    Particle(int                         Tag,       ///< Tag of the particle
             Vec3_t              const & x0,        ///< Initial position
             Vec3_t              const & v0,        ///< Initial velocity
             double                      mass,      ///< Initial mass of the material
             double                      vol0       ///< Initial volume of the material
             );    

    //Data
    int Tag;                ///< Tag of particles to identify different groups
    double V0;              ///< Initial volume assigned to the particle
    double V;               ///< Volume assigned to the particle
    double Vt;              ///< Volume of the assigned tetrahedron
    double m;               ///< Mass of the particle
    Vec3_t x0;              ///< Initial position of particle
    Vec3_t x;               ///< Position of particle
    Vec3_t xb;              ///< Previous position of particle
    Vec3_t v;               ///< Velocity of the particle
    Vec3_t vf;              ///< Fixed velocity for dirichlet BCs
    Vec3_t b;               ///< Body Force
    Vec3_t h;               ///< Hydraulic Force
    Mat3_t S;               ///< Stress of the particle
    Mat3_t E;               ///< Strain of the particle
    Mat3_t L;               ///< Velocity gradient tensor
    Mat3_t F;               ///< Deformation tensor
    Vec3_t Dx0;             ///< Size of the particle for GIMP initial
    Vec3_t Dx;              ///< Size of the particle for GIMP
    bool   vxf, vyf, vzf;   ///< Fixed components of velocity
    size_t Corners[4];      ///< Indexes of corners associated with this particle
    Vec3_t Gradvec[4];      ///< Vector of geometry for gradients
    Vec3_t * Vcorner[4];    ///< Pointer to corners
    Array<size_t> Nodes;    ///< Array of node indexes relevant to this particle
    Array<double> Sfunc;    ///< Relevant values of the shape function
    Array<Vec3_t> GSf;      ///< Relevant values of the shape function gradient

    //Data for solid properties
	double G;				///< Shear modulus
	double K;				///< bulk modulus


    #ifdef USE_OMP
    omp_lock_t      lck;    ///< to protect variables in multithreading
    #endif

    //Methods
    bool   IsFree             () {return !vxf&&!vyf&&!vzf;};                                  ///< Ask if the particle has any constrain in its movement
    void   FixVeloc           (double vx=0.0, double vy=0.0, double vz=0.0);  ///< Fix all velocities
    void   CalcVol            (double dt);                                    ///< Caculate volume and other geoemtric properties
    void   CalcVolCPI         (double dt);                                    ///< Caculate volume and other geoemtric properties for CPI
    void   CalcStress         (double dt);                                    ///< Calculate Stress tensor and other related qunatities
    void   Reset              (double dt);                                    ///< Reset relevant quantites to zero 
    
};

inline Particle::Particle (int TheTag, Vec3_t const & xt0, Vec3_t const & v0, double mass0, double Vol0)
{
    Tag = TheTag;
    x0  = xt0;
    x   = x0;
    v   = v0;
    vf  = v0;
    m   = mass0;
    V0  = Vol0;
    V   = V0;
    Vt  = V0;
    Dx0 = Vec3_t(pow(V0,1.0/3.0),pow(V0,1.0/3.0),pow(V0,1.0/3.0));
    Dx  = Vec3_t(pow(V0,1.0/3.0),pow(V0,1.0/3.0),pow(V0,1.0/3.0));
    set_to_zero(b);
    set_to_zero(h);
    set_to_zero(S);
    set_to_zero(E);
    set_to_zero(L);
    F = OrthoSys::I;
    vxf = false;
    vyf = false;
    vzf = false;
    K   = 1.0e4;
    G   = 0.3e4;

    #ifdef USE_OMP
    omp_init_lock(&lck);    ///< to protect variables in multithreading
    #endif
}

inline void Particle::FixVeloc (double vx, double vy, double vz)
{
    v   = vx, vy, vz;
    vf  = vx, vy, vz;
    vxf = true; vyf = true; vzf = true; 
}

inline void Particle::CalcVol(double dt)
{
    // Calculate deformation tensor
    F =  (OrthoSys::I + dt*L)*F;
	
    // Calculating the deformation of a single material point 
    //Dx(0) =  Dx0(0)*sqrt(F(0,0)*F(0,0) + F(1,0)*F(1,0) + F(2,0)*F(2,0));
	//Dx(1) =  Dx0(1)*sqrt(F(0,1)*F(0,1) + F(1,1)*F(1,1) + F(2,1)*F(2,1));
	//Dx(2) =  Dx0(2)*sqrt(F(0,2)*F(0,2) + F(1,2)*F(1,2) + F(2,2)*F(2,2));
    Dx(0) =  Dx0(0)*F(0,0);
	Dx(1) =  Dx0(1)*F(1,1);
	Dx(2) =  Dx0(2)*F(2,2);
    V = V0*Det(F);
}

inline void Particle::CalcVolCPI(double dt)
{
    // Calculate deformation tensor
    F =  (OrthoSys::I + dt*L)*F;
    V = V0*Det(F);
	
    // Calculating the deformation of a single material point 
    Vec3_t a14 = *Vcorner[0] - *Vcorner[3];
    Vec3_t a24 = *Vcorner[1] - *Vcorner[3];
    Vec3_t a34 = *Vcorner[2] - *Vcorner[3];
    Vec3_t a31 = *Vcorner[2] - *Vcorner[0];
    Vec3_t a21 = *Vcorner[1] - *Vcorner[0];

    Vt = dot(a34,cross(a24,a14))/6.0;
    V = Vt;

    Gradvec[0] = cross(a34,a24);
    Gradvec[1] = cross(a14,a34);
    Gradvec[2] = cross(a24,a14);
    Gradvec[3] = cross(a21,a31);

    //if (-dot(a21,Gradvec[0]<0.0)) Gradvec[0] *= -1.0;
    //if ( dot(a24,Gradvec[1]<0.0)) Gradvec[1] *= -1.0;
    //if ( dot(a31,Gradvec[2]<0.0)) Gradvec[2] *= -1.0;
    //if (-dot(a14,Gradvec[3]<0.0)) Gradvec[3] *= -1.0;
    

    //Calculate the current position
    x = 0.25*(*Vcorner[0]+*Vcorner[1]+*Vcorner[2]+*Vcorner[3]);
}

inline void Particle::CalcStress(double dt)
{
    // Calculating the constitutive relation (linear elasticity)
    Mat3_t De = 0.5*dt*(L + ~L);
    Mat3_t Dw = 0.5*dt*(L - ~L);
    E = E + De;
    S = S + (Dw*S - S*Dw) + K*Trace(De)*OrthoSys::I + 2.0*G*(De-(1.0/3.0)*Trace(De)*OrthoSys::I);
}

inline void Particle::Reset(double dt)
{
    Nodes.Resize(0);
    Sfunc.Resize(0);
    GSf  .Resize(0);
    if (vxf)
    {
        v(0) = vf(0);
        x(0) = xb(0) + v(0)*dt;
        xb(0) = x(0);
    }
    if (vyf)
    {
        v(1) = vf(1);
        x(1) = xb(1) + v(1)*dt;
        xb(1) = x(1);
    }
    if (vzf)
    {
        v(2) = vf(2);
        x(2) = xb(2) + v(2)*dt;
        xb(2) = x(2);
    }
}

#ifdef USE_CUDA
struct ParticleCU
{

}
#endif

}
#endif //MECHSYS_MPM_PARTICLE_H
