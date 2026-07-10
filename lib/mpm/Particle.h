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
#include <mechsys/mpm/eos.h>
#include <mechsys/mpm/fluid_filter.h>

namespace MPM
{

    inline double DoubleDot(Mat3_t const & A, Mat3_t const & B)
{
    return A(0,0)*B(0,0) + A(0,1)*B(0,1) + A(0,2)*B(0,2)
         + A(1,0)*B(1,0) + A(1,1)*B(1,1) + A(1,2)*B(1,2)
         + A(2,0)*B(2,0) + A(2,1)*B(2,1) + A(2,2)*B(2,2);
}

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

    int    Constitutive;    ///< 0=Elastic, 1=MohrCoulomb, 2=DruckerPrager, 3=Hypoplastic, 4=Newtonian fluid
    //Data for solid properties
	double G;				///< Shear modulus
	double K;				///< bulk modulus
    // Mohr-Coulomb parameters
    double  Cohesion;  ///< Cohesion [Pa]
    double  Friction;  ///< Friction angle [radians]
    double  Dilation;  ///< Dilation angle [radians]

    int    Fail;           ///< 2=Associated flow, 3=Non-associated flow
    int    ct;             ///< Legacy SPH integration flag; not used by MPM DP update
    Mat3_t Sigmab;         ///< Legacy backup stress kept for compatibility
    Mat3_t Strainb;        ///< Legacy backup strain kept for compatibility
    Mat3_t FilStress;      ///< Legacy SPH stress filter; not used by the MPM constitutive update

    double c1hp;            ///< Hypoplastic dimensionless parameter c1
    double c2hp;            ///< Hypoplastic dimensionless parameter c2
    double c3hp;            ///< Hypoplastic dimensionless parameter c3
    double c4hp;            ///< Hypoplastic dimensionless parameter c4
    double HypoPressureMin; ///< Compressive pressure floor [Pa] for hypoplastic zero-stress starts

    double Mu;               ///< Dynamic viscosity for Newtonian fluid [Pa s]
    double Cs;               ///< Artificial sound speed for weakly-compressible fluid [m/s]
    int    EOSType;          ///< Fluid EOS selector: EOS_TAIT or EOS_LINEAR
    double EOSexp;           ///< Tait/Cole EOS exponent [-]
    double AVLin;            /// Artificial Viscosity Linear coefficient
    double AVQuad;           /// Artificial Viscosity Quadratic coefficient
    double FluidDx;          ///< Characteristic grid spacing for the pressure filter [m]
    double Pressure;         ///< Fluid/solid pressure diagnostic [Pa]

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

    // Constitutive models
    void   ElasticModel       (double dt);
    void   MohrCoulombModel   (double dt);
    void   DruckerPragerModel (double dt);
    void   HypoplasticModel   (double dt);
    void   NewtonianFluidModel(double dt);
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

    Constitutive = 0;  // Default to Elastic
    // Initialize linear elasticity parameters (default values)
    K   = 1.0e4;    // Bulk modulus (Pa)
    G   = 0.3e4;    // Shear modulus (Pa)

    // Initialize Mohr-Coulomb parameters (default values)
    Cohesion = 0.0;              // Cohesion (Pa)
    Friction = 28.0*M_PI/180.0;  // Default friction angle [radians]
    Dilation = 0.0;              // Default dilation angle [radians]

    // Initialize Drucker-Prager specific variables
    Fail = 2;          // Default: associated flow
    ct = 0;            // Legacy SPH integration flag; not used by MPM updates
    set_to_zero(Sigmab);
    set_to_zero(Strainb);
    set_to_zero(FilStress);
    // Default dry granular-flow hypoplastic constants from Guo et al. (2016),
    // Table 3, for 1.0 mm dry spherical polystyrene beads.
    c1hp = -50.0;
    c2hp = -832.1;
    c3hp = -832.1;
    c4hp = 2369.9;
    // Legacy member kept for compatibility; the SPH-faithful hypoplastic
    // update below uses the original zero-trace guard directly.
    HypoPressureMin = 0.0;
    Mu              = 1.0e-3; // Water at room temperature [Pa.s]
    Cs              = 20.0;   // Artificial WC fluid sound speed; set per case.
    EOSType         = EOS_TAIT;
    EOSexp          = 7.0;    // Standard water exponent in Tait/Cole EOS.
    AVLin           = 0.0;
    AVQuad          = 0.0;
    FluidDx         = 0.0;
    Pressure        = 0.0;

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
    if (Constitutive == 4)
    {
        // Ren et al. (2025) weakly-compressible MPM water update:
        // V^{n+1} = (1 + dt*tr(D)) V^n, with D the symmetric velocity gradient.
        Mat3_t D = 0.5*(L + ~L);
        double Jinc = 1.0 + dt*Trace(D);
        if (Jinc < 1.0e-12) Jinc = 1.0e-12; // Prevent non-physical negative volume after a failed step.
        V *= Jinc;

        // F is not a fluid state variable here; keep it isotropic and consistent
        // with the updated volume for diagnostics and particle-size bookkeeping.
        double J = (V0 > 1.0e-12) ? V/V0 : 1.0;
        double stretch = (J > 1.0e-12) ? std::pow(J,1.0/3.0) : 1.0e-4;
        F = stretch*OrthoSys::I;
        Dx = Vec3_t(stretch*Dx0(0),stretch*Dx0(1),stretch*Dx0(2));
        return;
    }

    // Calculate deformation tensor for solid/soil constitutive models.
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
    switch (Constitutive)
    {
        case 0:
            ElasticModel(dt);
            break;
        case 1:
            MohrCoulombModel(dt);
            break;
        case 2:
            DruckerPragerModel(dt);
            break;
        case 3:
            HypoplasticModel(dt);
            break;
        case 4:
            NewtonianFluidModel(dt);
            break;
        default:
            ElasticModel(dt);
            break;
    }
}

inline void Particle::ElasticModel(double dt)
{
    // Calculating the constitutive relation (linear pure-elasticity, no failure criterion implemented)
    Mat3_t I = OrthoSys::I;
    Mat3_t De = 0.5*dt*(L + ~L); // Rate of deformation (symmetric part of L)
    Mat3_t Dw = 0.5*dt*(L - ~L); // Spin tensor (skew-symmetric part of L)
    E = E + De;                  // Update elastic strain (small strain assumption)
    S = S + (Dw*S - S*Dw) + K*Trace(De)*I + 2.0*G*(De-(1.0/3.0)*Trace(De)*I);
}

inline void Particle::MohrCoulombModel(double dt)
{
    // References:
    // - Drucker, D. C. and Prager, W. (1952), "Soil mechanics and plastic analysis for limit design".
    // - Abbo, A. J. and Sloan, S. W. (1995), "A smooth hyperbolic approximation to the Mohr-Coulomb yield criterion".
    // - Simo, J. C. and Hughes, T. J. R. (1998), "Computational Inelasticity".
    // Smooth Mohr-Coulomb p-q approximation, not the exact hexagonal
    // principal-stress Mohr-Coulomb surface.
    // Elastic predictor. This code uses the same stress sign convention as
    // ElasticModel: compression is negative stress, so p = -tr(S)/3 below.
    Mat3_t I = OrthoSys::I;
    Mat3_t De = 0.5*dt*(L + ~L); // Rate of deformation (symmetric part of L)
    Mat3_t Dw = 0.5*dt*(L - ~L); // Spin tensor (skew-symmetric part of L)
    E = E + De;                  // Update elastic strain (small strain assumption)

    // Elastic predictor using the Jaumann rate of stress:
    Mat3_t S_trial = S + (Dw*S - S*Dw) + K*Trace(De)*I + 2.0*G*(De-(1.0/3.0)*Trace(De)*I);

    const double tol = 1.0e-12;
 //   const double phi = Friction;
    const double psi = (Fail == 3) ? Dilation : Friction; // If Fail == 3, use non‑associated flow; otherwise associated (psi = phi).
    const double sin_phi = sin(Friction);
    const double cos_phi = cos(Friction);
    const double sin_psi = sin(psi);
    const double mc_denom = 3.0 - sin_phi; // Denominator used in the definition of M and k.

    if (mc_denom <= tol)
    { // If sin_psi >= 3 (physically impossible), we skip plasticity and keep the elastic predictor.
        S = S_trial;
        return;
    }

    const double M = 6.0*sin_phi/mc_denom; // Slope of the yield surface in (p, q) space.
    const double k = 6.0*Cohesion*cos_phi/mc_denom; // Intercept (cohesion term).
    const double psi_denom = 3.0 - sin_psi;

    // Slope used for the flow rule (plastic potential). For associated flow, Mpsi = M;
    // for non‑associated, it uses the dilation angle, psi.
    const double Mpsi = (fabs(psi_denom) > tol) ? 6.0*sin_psi/psi_denom : M;

    // Convert trial stress to p-q invariants. q measures shear intensity;
    // q_yield is the frictional shear strength available at pressure p.
    double I1 = Trace(S_trial); // First invariant (trace)
    double p = -I1/3.0; // Mean pressure (positive in compression, hence the minus sign).
    double p_apex = (M > tol) ? -k/M : -1.0e30; // The pressure at the apex of the cone thus the apex is at p=0 when c=0
    // and when p<papex​ would be tensile and is mapped directly to the apex.

    // Cohesionless or weakly cohesive sand should not carry tensile mean stress.
    // If the trial state lies beyond the apex, project directly to the apex.
    if (p <= p_apex)
    {
        S = -p_apex*I;
        return;
    }

    Mat3_t dev_S = S_trial - (I1/3.0)*I; // Deviatoric part of the trial stress.
    // Second invariant of the deviatoric stress
    double J2 = 0.5*(dev_S(0,0)*dev_S(0,0) + dev_S(1,1)*dev_S(1,1) + dev_S(2,2)*dev_S(2,2))
              + dev_S(0,1)*dev_S(1,0) + dev_S(0,2)*dev_S(2,0) + dev_S(1,2)*dev_S(2,1);
    if (J2 < 0.0 && J2 > -tol) J2 = 0.0;
    double q = sqrt(3.0*J2); // Shear stress intensity
    double q_yield = M*p + k; // Shear strength at the current pressure according to the Mohr‑Coulomb cone.
    if (q <= q_yield + tol)
    { // If the trial stress lies inside or on the yield surface, no plastic correction is needed; the trial stress is accepted.
        S = S_trial;
        return;
    }

    // Return mapping on the smooth p-q cone. Fail==3 uses the dilation angle
    // for plastic flow; otherwise associated flow uses psi = phi.
    double denom = 3.0*G + K*M*Mpsi; // The denominator comes from the consistency condition:
    if (denom <= tol)
    {
        S = S_trial;
        return;
    }

    double dgamma = (q - q_yield)/denom; // Plastic multiplier
    double q_new = q - 3.0*G*dgamma;
    double p_new = p + K*Mpsi*dgamma;

    if (p_new <= p_apex || q_new <= 0.0)
    {// If the corrected state goes beyond the apex or below zero shear, we clamp it to the apex (pure hydrostatic tension).
        p_new = p_apex;
        q_new = 0.0;
    }

    // Reconstruct the stress tensor from the updated invariants:
    // The deviatoric direction is preserved from the trial state.
    double scale = (q > tol) ? q_new/q : 0.0;
    S = -p_new*I + scale*dev_S;
}
// ========================================================================


// ========================================================================
// DRUCKER-PRAGER CONSTITUTIVE MODEL (CPU Version)
// ========================================================================
inline void Particle::DruckerPragerModel(double dt)
{
    // References:
    // - Bui, H. H., Fukagawa, R., Sako, K. and Ohno, S. (2008),
    //   "Lagrangian meshfree particles method (SPH) for large deformation
    //   and failure flows of geomaterial using elastic-plastic soil
    //   constitutive model".
    // - Drucker, D. C. and Prager, W. (1952), "Soil mechanics and plastic
    //   analysis for limit design".
    // Bui et al. use f = sqrt(J2) + alpha*I1 - kf. With this code's
    // compression-negative stress convention, p = -I1/3 and q = sqrt(3*J2),
    // so the same cone is q - M*p - k = 0. The update below keeps Bui's DP
    // yield and non-associated potential, but performs a stable p-q return
    // mapping for the MPM step instead of the old SPH Modified-Verlet branch.
    const double tol = 1.0e-12;
    Mat3_t I = OrthoSys::I;
    Mat3_t S_old = S;

    // Elastic predictor with the Jaumann objective stress rate. De and Dw are
    // increments because dt is already included here.
    Mat3_t De = 0.5*dt*(L + ~L);
    Mat3_t Dw = 0.5*dt*(L - ~L);
    E = E + De;
    Mat3_t S_trial = S + (Dw*S - S*Dw)
                   + K*Trace(De)*I
                   + 2.0*G*(De - (Trace(De)/3.0)*I);

    // Three-dimensional circumscribed Drucker-Prager cone used in Bui et al.
    // The older SPH code also had a plane-strain formula through its I tensor;
    // this MPM particle stress is stored as a full 3D tensor, so we do not
    // infer a transient 2D mode from zero z-stress components.
    const double sin_phi = sin(Friction);
    const double cos_phi = cos(Friction);
    const double denom = sqrt(3.0)*(3.0 - sin_phi);
    if (denom <= tol)
    {
        S = S_trial;
        Sigmab = S_old;
        return;
    }
    const double alpha = 2.0*sin_phi/denom;
    const double kf = 6.0*Cohesion*cos_phi/denom;
    const double M = 3.0*sqrt(3.0)*alpha;
    const double k = sqrt(3.0)*kf;

    // Convert the DP cone to p-q invariants. The apex projection is the
    // tension-cracking treatment used by Bui et al. for non-cohesive soils.
    double I1 = Trace(S_trial);
    double p = -I1/3.0;
    double p_apex = (M > tol) ? -k/M : -1.0e30;
    if (p <= p_apex)
    {
        S = -p_apex*I;
        Sigmab = S_old;
        return;
    }

    Mat3_t dev_S = S_trial - (I1/3.0)*I;
    double J2 = 0.5*DoubleDot(dev_S, dev_S);
    if (J2 < 0.0 && J2 > -tol) J2 = 0.0;
    double q = sqrt(3.0*J2);
    double q_yield = M*p + k;

    if (q <= q_yield + tol)
    {
        S = S_trial;
        Sigmab = S_old;
        return;
    }

    // Associated flow uses the yield slope M. For Fail==3, Bui's plastic
    // potential g = sqrt(J2) + 3*I1*sin(psi) gives Mpsi = 9*sqrt(3)*sin(psi)
    // in the q - Mpsi*p form. Bui's granular collapse examples use psi = 0.
    const double Mpsi = (Fail == 3) ? 9.0*sqrt(3.0)*sin(Dilation) : M;
    const double plastic_denom = 3.0*G + K*M*Mpsi;
    if (plastic_denom <= tol)
    {
        S = S_trial;
        Sigmab = S_old;
        return;
    }

    const double dgamma = (q - q_yield)/plastic_denom;
    double q_new = q - 3.0*G*dgamma;
    double p_new = p + K*Mpsi*dgamma;

    if (p_new <= p_apex || q_new <= 0.0)
    {
        p_new = p_apex;
        q_new = 0.0;
    }

    const double scale = (q > tol) ? q_new/q : 0.0;
    S = -p_new*I + scale*dev_S;
    Sigmab = S_old;
}
// ========================================================================


// ========================================================================
// HYPOPLASTIC CONSTITUTIVE MODEL (CPU Version)
// ========================================================================
inline void Particle::HypoplasticModel(double dt)
{
    // References:
    // - Guo, X., Peng, C., Wu, W. and Wang, Y. (2016),
    //   "A hypoplastic constitutive model for debris materials",
    //   Acta Geotechnica 11:1217-1229. doi:10.1007/s11440-016-0494-0.
    // - Wu, W., Bauer, E. and Kolymbas, D. (1996), "Hypoplastic
    //   constitutive model with critical state for granular materials".
    // This mirrors the legacy SPH hypoplastic invariant terms. SPH-only
    // Modified-Verlet branching and FilStress smoothing are not used here.
    const double tol = 1.0e-12;
    Mat3_t I = OrthoSys::I;

    // D and W are the MPM equivalents of SPH StrainRate and RotationRate.
    // The SPH routine also includes these Zaremba-Jaumann objective terms.
    Mat3_t D = 0.5*(L + ~L);
    // MPM stores L as grad(N) otimes v, i.e. the transpose of the usual
    // velocity gradient. Use the opposite skew part so the Jaumann term has
    // the same spin sign as LOQUAT's hypoplastic update.
    Mat3_t W = 0.5*(~L - L);

    // The SPH model subtracts cohesion from the normal stresses before using
    // the stress in the hypoplastic law. Keep that convention isolated here.
    Mat3_t Stress = S;
    if (Cohesion > tol)
    {
        Stress(0,0) -= Cohesion;
        Stress(1,1) -= Cohesion;
        Stress(2,2) -= Cohesion;
    }

    // Match the SPH zero-trace guard exactly: only the effective stress used
    // for the invariants is reset, while the actual stress S is still advanced.
    double I1 = Trace(Stress);
    if (fabs(I1) < 4.0)
    {
        I1 = -4.0;
        Stress(0,0) = -1.0;
        Stress(1,1) = -1.0;
        Stress(2,2) = -2.0;
        Stress(0,1) = 0.0;
        Stress(0,2) = 0.0;
        Stress(1,2) = 0.0;
        Stress(1,0) = 0.0;
        Stress(2,0) = 0.0;
        Stress(2,1) = 0.0;
    }

    Mat3_t ShearStress = Stress - (I1/3.0)*I;
    double Dnorm2 = DoubleDot(D, D);
    if (Dnorm2 < 0.0 && Dnorm2 > -tol) Dnorm2 = 0.0;
    const double normD = sqrt(Dnorm2); // Frobenius norm of the strain-rate tensor.
    const double trD = Trace(D); // Trace of the strain-rate tensor.
    const double trSD = DoubleDot(Stress, D); // Inner double product of stress and strain rate, tr(TD)
    // The legacy SPH Mat3hpMVerlet routine uses the guarded effective stress
    // only for invariants, then resets Stress = Sigma for the rate terms below.
    Mat3_t StressForRate = S;
    // D and W are rates, so S_rate is a stress rate; dt is applied once below.
    Mat3_t S_rate = (W*S - S*W)
                  + (c1hp*I1)*D
                  + (c2hp*trD + c3hp*trSD/I1)*StressForRate
                  + (c4hp*normD)*(StressForRate + ShearStress);
    
    Mat3_t S_trial = S + dt*S_rate;
    S_trial = 0.5*(S_trial + ~S_trial);

    // Same tensile mean-stress cap as the SPH routine: for c=0 the mean stress
    // is not allowed to become tensile.
    double mean_stress = Trace(S_trial)/3.0;
    if (mean_stress > Cohesion)
    {
        S_trial = S_trial - (mean_stress - Cohesion)*I;
    }

    // The hypoplastic law does not use accumulated strain as a state variable,
    // but MPM writes E as "Strain" in the output, matching the SPH bookkeeping.
    E = E + dt*D;
    S = S_trial;
    Pressure = std::max(0.0,-Trace(S)/3.0);
}
// ========================================================================

// ========================================================================
// NEWTONIAN FLUID MODEL (CPU Version)
// ========================================================================
inline void Particle::NewtonianFluidModel(double dt)
{
    // References:
    // - Chen, Z.-P., Zhang, X., Sze, K.Y., Kan, L., and Qiu, X.-M. (2018).
    //   "v-p material point method for weakly compressible problems".
    //   Computers and Fluids, 176, pp. 170-181. doi:10.1016/j.compfluid.2018.09.005.
    // - Ren, S., Liu, Z., Trujillo-Vela, M.G., Galindo-Torres, S.A., Tian, X.,
    //   and Zhang, P. (2025). "Simulation of solitary wave generation and
    //   wave-structure interactions using an MPM-SDEM coupling scheme".
    //   Ocean Engineering, 340, 122195. doi:10.1016/j.oceaneng.2025.122195.
    // Compression is negative stress in this MPM code, so S = tau - pI.
    const double tol = 1.0e-12;
    Mat3_t I = OrthoSys::I;
    Mat3_t D = 0.5*(L + ~L);
    const double trD = D(0,0) + D(1,1) + D(2,2);   // compute once
    //Mat3_t E_trial = E + dt*D;

    const double rho0 = (V0 > tol) ? m/V0 : 0.0;
    const double rho  = (V  > tol) ? m/V : rho0;
    double p = EOS(EOSType,rho,rho0,Cs,EOSexp);
    p += VNRLPressure(rho,trD,FluidDx,Cs,AVLin,AVQuad);
    // Keep the EOS gauge pressure signed; expansion needs a restoring pressure.
    Pressure = p;

    Mat3_t Ddev = D - (trD/3.0)*I;
    S = 2.0*Mu*Ddev - Pressure*I;

    // The fluid law is rate-based; E is kept only as an output diagnostic.
    //E = E_trial;
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
    size_t Nnodes;          ///< Number of nodes interacting with this particle
    int    Tag;             ///< Tag of particles to identify different groups
    real   V0;              ///< Initial volume assigned to the particle
    real   V;               ///< Volume assigned to the particle
    real   Vt;              ///< Volume of the assigned tetrahedron
    real   m;               ///< Mass of the particle
    real3  x0;              ///< Initial position of particle
    real3  x;               ///< Position of particle
    real3  xb;              ///< Previous position of particle
    real3  v;               ///< Velocity of the particle
    real3  vf;              ///< Fixed velocity for dirichlet BCs
    real3  b;               ///< Body Force
    real3  h;               ///< Hydraulic Force
    Mat3   S;               ///< Stress of the particle
    Mat3   E;               ///< Strain of the particle
    Mat3   L;               ///< Velocity gradient tensor
    Mat3   F;               ///< Deformation tensor
    real3  Dx0;             ///< Size of the particle for GIMP initial
    real3  Dx;              ///< Size of the particle for GIMP
    bool   vxf, vyf, vzf;   ///< Fixed components of velocity
    size_t Pnodes = 0;      ///< number of nodes interacting with the particle
    size_t Nodes[125];      ///< Array of node indexes relevant to this particle
    real   Sfunc[125];      ///< Relevant values of the shape function
    real3  GSf  [125];      ///< Relevant values of the shape function gradient

    int    Constitutive;    ///< 0=Elastic, 1=MohrCoulomb, 2=DruckerPrager, 3=Hypoplastic, 4=Newtonian fluid
    //Data for solid properties
	real G;				///< Shear modulus
	real K;				///< bulk modulus
    // Mohr-Coulomb parameters
    real  Cohesion;  ///< Cohesion [Pa]
    real  Friction;  ///< Friction angle [radians]
    real  Dilation;  ///< Dilation angle [radians]

    // Pre-computed values for efficiency
    real  sin_phi;          ///< sin(Friction angle)
    real  cos_phi;          ///< cos(Friction angle)
    real  sin_psi;          ///< sin(Dilation angle)
    real  cos_psi;          ///< cos(Dilation angle)

    // Drucker-Prager specific members
    int    Fail;           ///< 2=Associated flow, 3=Non-associated flow
    int    ct;             ///< Legacy SPH integration flag; not used by MPM DP update
    Mat3   Sigmab;         ///< Legacy backup stress kept for compatibility
    Mat3   Strainb;        ///< Legacy backup strain kept for compatibility
    Mat3   FilStress;      ///< Legacy SPH stress filter; not used by the MPM constitutive update

    real   c1hp;            ///< Hypoplastic dimensionless parameter c1
    real   c2hp;            ///< Hypoplastic dimensionless parameter c2
    real   c3hp;            ///< Hypoplastic dimensionless parameter c3
    real   c4hp;            ///< Hypoplastic dimensionless parameter c4
    real   HypoPressureMin; ///< Compressive pressure floor [Pa] for hypoplastic zero-stress starts

    real   Mu;              ///< Dynamic viscosity for Newtonian fluid [Pa s]
    real   Cs;              ///< Artificial sound speed for weakly-compressible fluid [m/s]
    int    EOSType;         ///< Fluid EOS selector: EOS_TAIT or EOS_LINEAR
    real   EOSexp;          ///< Tait/Cole EOS exponent [-]
    real   AVLin;           ///< Linear vNRL artificial-viscosity coefficient cL.
    real   AVQuad;          ///< Quadratic vNRL artificial-viscosity coefficient cQ.
    real   FluidDx;         ///< Characteristic grid spacing for the pressure filter [m]
    real   Pressure;        ///< Fluid/solid pressure diagnostic [Pa]

    __host__ __device__ void Reset (real dt)
    {
        Pnodes = 0;
        if (vxf)
        {
            v.x = vf.x;
            x.x = xb.x + v.x*dt;
            xb.x = x.x;
        }
        if (vyf)
        {
            v.y = vf.y;
            x.y = xb.y + v.y*dt;
            xb.y = x.y;
        }
        if (vzf)
        {
            v.z = vf.z;
            x.z = xb.z + v.z*dt;
            xb.z = x.z;
        }
    }

    __host__ __device__ void CalcVol (real dt)
    {
        Mat3 I;
        if (Constitutive == 4)
        {
            Mat3 D = 0.5*(L + L.transpose());
            real Jinc = 1.0 + dt*D.trace();
            if (Jinc < 1.0e-12) Jinc = 1.0e-12;
            V *= Jinc;

            real J = (V0 > 1.0e-12) ? V/V0 : 1.0;
            real stretch = (J > 1.0e-12) ? pow(J,1.0/3.0) : 1.0e-4;
            F = stretch*I;
            Dx.x = stretch*Dx0.x;
            Dx.y = stretch*Dx0.y;
            Dx.z = stretch*Dx0.z;
            return;
        }
        F =  (I + dt*L)*F;
        Dx.x =  Dx0.x*F(0,0);
	    Dx.y =  Dx0.y*F(1,1);
	    Dx.z =  Dx0.z*F(2,2);
        V = V0*F.determinant();
    }

    __host__ __device__ void CalcStress (real dt)
    {
        switch (Constitutive)
        {
            case 0:
                ElasticModel(dt);
                break;
            case 1:
                MohrCoulombModel(dt);
                break;
            case 2:
                DruckerPragerModel(dt);
                break;
            case 3:
                HypoplasticModel(dt);
                break;
            case 4:
                NewtonianFluidModel(dt);
                break;
            default:
                ElasticModel(dt);
                break;
        }
    }

    __host__ __device__ real DoubleDot(const Mat3 & A, const Mat3 & B) const
    {
        return A(0,0)*B(0,0) + A(0,1)*B(0,1) + A(0,2)*B(0,2)
             + A(1,0)*B(1,0) + A(1,1)*B(1,1) + A(1,2)*B(1,2)
             + A(2,0)*B(2,0) + A(2,1)*B(2,1) + A(2,2)*B(2,2);
    }

    __host__ __device__ void ElasticModel (real dt)
    {
        // Calculating the constitutive relation (linear elasticity)
        Mat3 I;
        Mat3 De = 0.5*dt*(L + L.transpose());
        Mat3 Dw = 0.5*dt*(L - L.transpose());
        E = E + De;

        S = S + (Dw*S - S*Dw) + K*De.trace()*I + 2.0*G*(De-(1.0/3.0)*De.trace()*I);
    }

    __host__ __device__ void MohrCoulombModel(real dt)
    {
        // References:
        // - Drucker, D. C. and Prager, W. (1952), "Soil mechanics and plastic analysis for limit design".
        // - Abbo, A. J. and Sloan, S. W. (1995), "A smooth hyperbolic approximation to the Mohr-Coulomb yield criterion".
        // - Simo, J. C. and Hughes, T. J. R. (1998), "Computational Inelasticity".
        // Smooth Mohr-Coulomb p-q approximation, not the exact hexagonal
        // principal-stress Mohr-Coulomb surface.
        // Elastic predictor. Compression is represented by negative stress,
        // so the pressure invariant used for strength is p = -tr(S)/3.
        Mat3 I;
        Mat3 De = 0.5 * dt * (L + L.transpose());
        Mat3 Dw = 0.5 * dt * (L - L.transpose());
        E = E + De;

        Mat3 S_trial = S + (Dw * S - S * Dw) + K * De.trace() * I + 2.0 * G * (De - (1.0/3.0) * De.trace() * I);

        const real tol = 1.0e-12;
        real sin_flow = (Fail == 3) ? sin_psi : sin_phi;
        real mc_denom = 3.0 - sin_phi;
        if (mc_denom <= tol)
        {
            S = S_trial;
            return;
        }

        real M = 6.0*sin_phi/mc_denom;
        real k = 6.0*Cohesion*cos_phi/mc_denom;
        real psi_denom = 3.0 - sin_flow;
        real Mpsi = (fabs(psi_denom) > tol) ? 6.0*sin_flow/psi_denom : M;

        // Convert trial stress to p-q invariants and evaluate frictional strength.
        real I1 = S_trial(0,0) + S_trial(1,1) + S_trial(2,2);
        real p = -I1/3.0;
        real p_apex = (M > tol) ? -k/M : -1.0e30;

        // Enforce the tensile apex so cohesionless sand cannot carry tension.
        if (p <= p_apex)
        {
            S = -p_apex*I;
            return;
        }

        Mat3 dev_S = S_trial - (I1/3.0)*I;
        real J2 = 0.5 * (dev_S(0,0)*dev_S(0,0) + dev_S(1,1)*dev_S(1,1) + dev_S(2,2)*dev_S(2,2))
                + dev_S(0,1)*dev_S(1,0) + dev_S(0,2)*dev_S(2,0) + dev_S(1,2)*dev_S(2,1);
        if (J2 < 0.0 && J2 > -tol) J2 = 0.0;

        real q = sqrt(3.0*J2);
        real q_yield = M*p + k;
        if (q <= q_yield + tol)
        {
            S = S_trial;
            return;
        }

        // Return mapping on the smooth p-q cone. Fail==3 uses dilation;
        // otherwise the flow is associated and uses psi = phi.
        real plastic_denom = 3.0*G + K*M*Mpsi;
        if (plastic_denom <= tol)
        {
            S = S_trial;
            return;
        }

        real dgamma = (q - q_yield)/plastic_denom;
        real q_new = q - 3.0*G*dgamma;
        real p_new = p + K*Mpsi*dgamma;

        if (p_new <= p_apex || q_new <= 0.0)
        {
            p_new = p_apex;
            q_new = 0.0;
        }

        real scale = (q > tol) ? q_new/q : 0.0;
        S = -p_new*I + scale*dev_S;
    }

    __host__ __device__ void DruckerPragerModel(real dt)
    {
        // References:
        // - Bui, H. H., Fukagawa, R., Sako, K. and Ohno, S. (2008),
        //   "Lagrangian meshfree particles method (SPH) for large deformation
        //   and failure flows of geomaterial using elastic-plastic soil
        //   constitutive model".
        // - Drucker, D. C. and Prager, W. (1952), "Soil mechanics and plastic
        //   analysis for limit design".
        // Bui et al. use f = sqrt(J2) + alpha*I1 - kf. With compression-negative
        // stress, p = -I1/3 and q = sqrt(3*J2), giving q - M*p - k = 0.
        // This CUDA path mirrors the CPU p-q return mapping for the MPM step.
        const real tol = 1.0e-12;
        Mat3 I;
        Mat3 S_old = S;

        // Elastic predictor with the Jaumann objective stress rate. De and Dw
        // already include dt, so the trial stress terms are not multiplied again.
        Mat3 De = 0.5*dt*(L + L.transpose());
        Mat3 Dw = 0.5*dt*(L - L.transpose());
        E = E + De;
        Mat3 S_trial = S + (Dw*S - S*Dw)
                     + K*De.trace()*I
                     + 2.0*G*(De - (De.trace()/3.0)*I);

        // Three-dimensional circumscribed Drucker-Prager cone from Bui et al.
        real denom = sqrt(3.0)*(3.0 - sin_phi);
        if (denom <= tol)
        {
            S = S_trial;
            Sigmab = S_old;
            return;
        }
        real alpha = 2.0*sin_phi/denom;
        real kf = 6.0*Cohesion*cos_phi/denom;
        real M = 3.0*sqrt(3.0)*alpha;
        real k = sqrt(3.0)*kf;

        // Convert the DP cone to p-q invariants and apply the apex projection.
        real I1 = S_trial.trace();
        real p = -I1/3.0;
        real p_apex = (M > tol) ? -k/M : -1.0e30;
        if (p <= p_apex)
        {
            S = -p_apex*I;
            Sigmab = S_old;
            return;
        }

        Mat3 dev_S = S_trial - (I1/3.0)*I;
        real J2 = 0.5*DoubleDot(dev_S, dev_S);
        if (J2 < 0.0 && J2 > -tol) J2 = 0.0;
        real q = sqrt(3.0*J2);
        real q_yield = M*p + k;

        if (q <= q_yield + tol)
        {
            S = S_trial;
            Sigmab = S_old;
            return;
        }

        // Bui's non-associated potential gives Mpsi = 9*sqrt(3)*sin(psi).
        real Mpsi = (Fail == 3) ? 9.0*sqrt(3.0)*sin_psi : M;
        real plastic_denom = 3.0*G + K*M*Mpsi;
        if (plastic_denom <= tol)
        {
            S = S_trial;
            Sigmab = S_old;
            return;
        }

        real dgamma = (q - q_yield)/plastic_denom;
        real q_new = q - 3.0*G*dgamma;
        real p_new = p + K*Mpsi*dgamma;

        if (p_new <= p_apex || q_new <= 0.0)
        {
            p_new = p_apex;
            q_new = 0.0;
        }

        real scale = (q > tol) ? q_new/q : 0.0;
        S = -p_new*I + scale*dev_S;
        Sigmab = S_old;
    }

    __host__ __device__ void HypoplasticModel(real dt)
    {
        // References:
        // - Guo, X., Peng, C., Wu, W. and Wang, Y. (2016),
        //   "A hypoplastic constitutive model for debris materials",
        //   Acta Geotechnica 11:1217-1229. doi:10.1007/s11440-016-0494-0.
        // - Wu, W., Bauer, E. and Kolymbas, D. (1996), "Hypoplastic
        //   constitutive model with critical state for granular materials".
        // CUDA mirrors the CPU MPM update and the legacy SPH hypoplastic form;
        // SPH Modified-Verlet and FilStress terms are not used.
        const real tol = 1.0e-12;
        Mat3 I;

        // D and W are the MPM equivalents of SPH StrainRate and RotationRate.
        Mat3 D = 0.5*(L + L.transpose());
        // L is stored as grad(N) otimes v, so use the opposite skew part to
        // match the spin convention in LOQUAT's hypoplastic stress update.
        Mat3 W = 0.5*(L.transpose() - L);

        Mat3 Stress = S;
        if (Cohesion > tol)
        {
            Stress(0,0) -= Cohesion;
            Stress(1,1) -= Cohesion;
            Stress(2,2) -= Cohesion;
        }

        real I1 = Stress.trace();
        if (fabs(I1) < 4.0)
        {
            I1 = -4.0;
            Stress(0,0) = -1.0;
            Stress(1,1) = -1.0;
            Stress(2,2) = -2.0;
            Stress(0,1) = 0.0;
            Stress(0,2) = 0.0;
            Stress(1,2) = 0.0;
            Stress(1,0) = 0.0;
            Stress(2,0) = 0.0;
            Stress(2,1) = 0.0;
        }

        Mat3 ShearStress = Stress - (I1/3.0)*I;
        real Dnorm2 = DoubleDot(D, D);
        if (Dnorm2 < 0.0 && Dnorm2 > -tol) Dnorm2 = 0.0;
        real normD = sqrt(Dnorm2);
        real trD = D.trace();
        real trSD = DoubleDot(Stress, D);
        // Match the CPU/SPH path: guarded Stress is used for invariants, but
        // the actual stress state S multiplies the hypoplastic rate terms.
        Mat3 StressForRate = S;
        // D and W are rates, so S_rate is a stress rate; dt is applied once below.
        Mat3 S_rate = (W*S - S*W)
                    + (c1hp*I1)*D
                    + (c2hp*trD + c3hp*trSD/I1)*StressForRate
                    + (c4hp*normD)*(StressForRate + ShearStress);

        Mat3 S_trial = S + dt*S_rate;
        S_trial = 0.5*(S_trial + S_trial.transpose());

        real mean_stress = S_trial.trace()/3.0;
        if (mean_stress > Cohesion)
        {
            S_trial = S_trial - (mean_stress - Cohesion)*I;
        }

        // Keep the output strain consistent with the CPU/SPH bookkeeping.
        E = E + dt*D;
        S = S_trial;
        Pressure = fmax(0.0,-S.trace()/3.0);
    }

    __host__ __device__ void NewtonianFluidModel(real dt)
    {
        // Weakly-compressible Newtonian water: tau = 2*mu*dev(D),
        // p from Tait/Cole EOS, and S = tau - pI in this MPM sign convention.
        const real tol = 1.0e-12;
        Mat3 I;
        Mat3 D = 0.5*(L + L.transpose());
        const real trD = D(0,0) + D(1,1) + D(2,2);   // compute once
        //Mat3 E_trial = E + dt*D;
        real rho0 = (V0 > tol) ? m/V0 : 0.0;
        real rho  = (V  > tol) ? m/V : rho0;
        real p = EOSCU(EOSType,rho,rho0,Cs,EOSexp);
        p += VNRLPressureCU(rho,trD,FluidDx,Cs,AVLin,AVQuad);
        // Keep the EOS gauge pressure signed; expansion needs a restoring pressure.
        Pressure = p;

        Mat3 Ddev = D - (trD/3.0)*I;
        S = 2.0*Mu*Ddev - Pressure*I;

        // Strain is not a fluid state variable; keep it as an output diagnostic.
        // E = E_trial;
    }

};

__host__ void UploadParticle (ParticleCU & PCu, const Particle & Par)
{
    PCu.Nnodes= 0;
    PCu.Tag   = Par.Tag   ;         
    PCu.V0    = Par.V0    ;        
    PCu.V     = Par.V     ;        
    PCu.Vt    = Par.Vt    ;        
    PCu.m     = Par.m     ;       
    PCu.x0 .x = Par.x0 (0);         
    PCu.x0 .y = Par.x0 (1);         
    PCu.x0 .z = Par.x0 (2);         
    PCu.x  .x = Par.x  (0);         
    PCu.x  .y = Par.x  (1);         
    PCu.x  .z = Par.x  (2);         
    PCu.xb .x = Par.xb (0);         
    PCu.xb .y = Par.xb (1);         
    PCu.xb .z = Par.xb (2);         
    PCu.v  .x = Par.v  (0);         
    PCu.v  .y = Par.v  (1);         
    PCu.v  .z = Par.v  (2);         
    PCu.vf .x = Par.vf (0);         
    PCu.vf .y = Par.vf (1);         
    PCu.vf .z = Par.vf (2);         
    PCu.b  .x = Par.b  (0);         
    PCu.b  .y = Par.b  (1);         
    PCu.b  .z = Par.b  (2);         
    PCu.h  .x = Par.h  (0);         
    PCu.h  .y = Par.h  (1);         
    PCu.h  .z = Par.h  (2);         
    PCu.S     = Par.S     ;
    PCu.E     = Par.E     ;    
    PCu.L     = Par.L     ;
    PCu.F     = Par.F     ;         
    PCu.Dx0.x = Par.Dx0(0);         
    PCu.Dx0.y = Par.Dx0(1);         
    PCu.Dx0.z = Par.Dx0(2);         
    PCu.Dx .x = Par.Dx (0);        
    PCu.Dx .y = Par.Dx (1);        
    PCu.Dx .z = Par.Dx (2);        
    PCu.vxf   = Par.vxf   ;
    PCu.vyf   = Par.vyf   ;
    PCu.vzf   = Par.vzf   ;

    PCu.Constitutive = Par.Constitutive;
    PCu.G     = Par.G     ;
    PCu.K     = Par.K     ;
    PCu.Cohesion     = Par.Cohesion     ;
    PCu.Friction = Par.Friction;
    PCu.Dilation = Par.Dilation;
    // Angles are already radians in the main code; precompute trig for kernels.
    PCu.sin_phi = sin(PCu.Friction);
    PCu.cos_phi = cos(PCu.Friction);
    PCu.sin_psi = sin(PCu.Dilation);
    PCu.cos_psi = cos(PCu.Dilation);

    // Upload Drucker-Prager specific data
    PCu.Fail     = Par.Fail;
    PCu.ct       = Par.ct;
    PCu.Sigmab   = Par.Sigmab;
    PCu.Strainb  = Par.Strainb;
    PCu.FilStress = Par.FilStress;
    PCu.c1hp = Par.c1hp;
    PCu.c2hp = Par.c2hp;
    PCu.c3hp = Par.c3hp;
    PCu.c4hp = Par.c4hp;
    PCu.HypoPressureMin = Par.HypoPressureMin;
    PCu.Mu              = Par.Mu;
    PCu.Cs              = Par.Cs;
    PCu.EOSType         = Par.EOSType;
    PCu.EOSexp          = Par.EOSexp;
    PCu.AVLin           = Par.AVLin;
    PCu.AVQuad          = Par.AVQuad;
    PCu.FluidDx         = Par.FluidDx;
    PCu.Pressure        = Par.Pressure;
}

__host__ void DnloadParticle (const ParticleCU & PCu, Particle & Par)
{
    Par.V      = PCu.V;
    Par.Vt     = PCu.Vt;
    Par.x  (0) = PCu.x  .x;         
    Par.x  (1) = PCu.x  .y;         
    Par.x  (2) = PCu.x  .z;         
    Par.v  (0) = PCu.v  .x;         
    Par.v  (1) = PCu.v  .y;         
    Par.v  (2) = PCu.v  .z;         
    copyMat3(Par.S,PCu.S) ;
    copyMat3(Par.E,PCu.E) ;
    copyMat3(Par.L,PCu.L) ;
    copyMat3(Par.F,PCu.F) ;
    Par.Dx(0) = PCu.Dx.x;
    Par.Dx(1) = PCu.Dx.y;
    Par.Dx(2) = PCu.Dx.z;

    // Download Drucker-Prager specific data
    Par.Fail    = PCu.Fail;
    Par.ct      = PCu.ct;
    copyMat3(Par.Sigmab, PCu.Sigmab);
    copyMat3(Par.Strainb, PCu.Strainb);
    copyMat3(Par.FilStress, PCu.FilStress);
    Par.c1hp = PCu.c1hp;
    Par.c2hp = PCu.c2hp;
    Par.c3hp = PCu.c3hp;
    Par.c4hp = PCu.c4hp;
    Par.HypoPressureMin = PCu.HypoPressureMin;
    Par.Mu              = PCu.Mu;
    Par.Cs              = PCu.Cs;
    Par.EOSType         = PCu.EOSType;
    Par.EOSexp          = PCu.EOSexp;
    Par.AVLin           = PCu.AVLin;
    Par.AVQuad          = PCu.AVQuad;
    Par.FluidDx         = PCu.FluidDx;
    Par.Pressure        = PCu.Pressure;
}

#endif

}
#endif //MECHSYS_MPM_PARTICLE_H
