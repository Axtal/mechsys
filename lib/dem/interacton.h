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

#ifndef MECHSYS_DEM_INTERACTON_H
#define MECHSYS_DEM_INTERACTON_H

// Std Lib
#include <math.h>
//#include <map>
#include <unordered_map>
#include <vector>
#include <utility>

// MechSys
#include <mechsys/dem/particle.h>

namespace DEM
{

// typedefs
typedef std::unordered_map<size_t,Vec3_t> FrictionMap_t;
//typedef std::map<std::pair<int,int>,Vec3_t> FrictionMap_t;
typedef Array<std::pair<int,int> > ListContacts_t;

class Interacton   //General class for interactons
{
public:

    // Constructor and destructor
    Interacton () {}           ///< Default constructor
    virtual ~Interacton () {}  ///< Destructor

    // Methods
    virtual bool UpdateContacts   (double alpha, Vec3_t const & Per = OrthoSys::O, size_t const iter=0) =0;    ///< Update contacts by verlet algorithm
    virtual bool CalcForce        (double dt = 0.0, Vec3_t const & Per = OrthoSys::O, size_t const iter=0, size_t const contactlaw=0) =0; ///< Calculates the contact force between particles
    virtual void UpdateParameters (size_t contactlaw=0) =0;                ///< Update the parameters

    // Data
    Particle     * P1;        ///< First particle
    Particle     * P2;        ///< Second particle
    size_t         I1;        ///< Index of the first particle
    size_t         I2;        ///< Index of the second particle
    Vec3_t         F1;        ///< Provisional force  for particle 1
    Vec3_t         F2;        ///< Provisional force  for particle 2
    Vec3_t         T1;        ///< Provisional torque for particle 1
    Vec3_t         T2;        ///< Provisional torque for particle 2
};

class CInteracton: public Interacton // Interacton for collision
{
public:
    // Constructor
    CInteracton (Particle * Pt1, Particle * Pt2, size_t contaclaw=0); ///< Constructor requires pointers to both particles
    CInteracton () {};

    // Methods
    virtual bool UpdateContacts   (double alpha, Vec3_t const & Per = OrthoSys::O, size_t const iter=0);    ///< Update contacts by verlet algorithm
    virtual bool CalcForce        (double dt = 0.0, Vec3_t const & Per = OrthoSys::O, size_t const iter=0, size_t const contaclaw = 0); ///< Calculates the contact force between particles
    virtual void UpdateParameters (size_t contactlaw=0);              ///< Update the parameters

    // Data
    bool           First;     ///< Is it the first collision?
    bool           BothFree;  ///< Are both particles free?
    double         Kn;        ///< Normal stiffness
    double         Kt;        ///< Tengential stiffness
    double         Gn;        ///< Normal viscous coefficient
    double         Gt;        ///< Tangential viscous coefficient
    double         Mu;        ///< Microscpic coefficient of friction
    double         Epot;      ///< Potential elastic energy
    double         dEvis;     ///< Energy dissipated in viscosity at time step
    double         dEfric;    ///< Energy dissipated by friction at time step
    size_t         Nc;        ///< Number of contacts
    size_t         Nsc;       ///< Number of sliding contacts
    size_t         Nr;        ///< Number of rolling contacts (only for spheres)
    Vec3_t         Fn;        ///< Normal force between elements
    Vec3_t         Fnet;      ///< Net normal force
    Vec3_t         Ftnet;     ///< Net tangential force
    //Vec3_t         Xc;        ///< Net Position of the contact
    Vec3_t         Branch;    ///< Branch vector
    ListContacts_t Lee;       ///< List of edge-edge contacts 
    ListContacts_t Lvf;       ///< List of vertex-face contacts 
    ListContacts_t Lfv;       ///< List of face-vertex contacts
    ListContacts_t Lvt;       ///< List of vertex-torus contacts
    ListContacts_t Ltv;       ///< List of torus-vertex contacts
    ListContacts_t Lvc;       ///< List of vertex-cylinder contacts
    ListContacts_t Lcv;       ///< List of cylinder-vertex contacts
    FrictionMap_t  Fdee;      ///< Static friction displacement for pair of edges
    FrictionMap_t  Fdvf;      ///< Static friction displacement for pair of vertex-face
    FrictionMap_t  Fdfv;      ///< Static friction displacement for pair of face-vertex
    FrictionMap_t  Fdvt;      ///< Static friction displacement for pair of vertex-torus
    FrictionMap_t  Fdtv;      ///< Static friction displacement for pair of torus-vertex
    FrictionMap_t  Fdvc;      ///< Static friction displacement for pair of vertex-cylinder
    FrictionMap_t  Fdcv;      ///< Static friction displacement for pair of cylinder-vertex
protected:
    template<typename FeatureA_T, typename FeatureB_T>
    bool _update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, ListContacts_t & L, double dt, Vec3_t const & Per = OrthoSys::O, size_t const iter=0);
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_contacts        (FeatureA_T & A, FeatureB_T & B, ListContacts_t & L, double alpha, Vec3_t const & Per = OrthoSys::O, size_t const iter=0);
};

class CInteractonSphere: public CInteracton // Collision interacton for spheres
{
public:
    // Methods 
    CInteractonSphere (Particle * Pt1, Particle * Pt2, size_t contaclaw=0); ///< Constructor requires pointers to both particles
    bool UpdateContacts (double alpha, Vec3_t const & Per = OrthoSys::O, size_t const iter=0);                 ///< Update contacts by verlet algorithm
    bool CalcForce (double dt = 0.0, Vec3_t const & Per = OrthoSys::O, size_t const iter=0, size_t contactlaw=0);    ///< Calculates the contact force between particles
    void UpdateParameters (size_t contactlaw=0);                           ///< Update the parameters

    // Data
    //ListContacts_t Lvv;                            ///< List of edge-edge contacts 
    //FrictionMap_t  Fdvv;      ///< Static friction displacement for pair of edges
    Vec3_t         Fdvv;                           ///< Static Friction displacement for the vertex vertex pair
    Vec3_t         Fdr;                            ///< Rolling displacement 
    double         beta;                           ///< Rolling stiffness coefficient
    double         eta;                            ///< Plastic moment coefficient
    double         Bn;                             ///< Elastic normal constant for the cohesion
    double         Bt;                             ///< Elastic tangential constant for the cohesion
    double         eps;                            ///< Maximun strain before fracture
    bool           cohesion;                       ///< A flag to determine if cohesion must be used for spheres or not

};

class BInteracton: public Interacton // Interacton for cohesion
{
public:
    BInteracton (Particle * Pt1, Particle * Pt2, size_t Fi1, size_t Fi2); ///< Constructor requires pointers to both particles and alos the indixes of the size they share

    // Methods
    bool UpdateContacts (double alpha, Vec3_t const & Per = OrthoSys::O, size_t const iter=0);    ///< Update contacts by verlet algorithm
    bool CalcForce      (double dt = 0.0, Vec3_t const & Per = OrthoSys::O, size_t const iter=0, size_t const contactlaw=0); ///< Calculates the contact force between particles
    void UpdateParameters (size_t contactlaw=0);              ///< Update the parameters in case they change

    // Data
    size_t IF1;                            ///< Index of the shared face for particle 1
    size_t IF2;                            ///< Index of the shared face for particle 2
    double Area;                           ///< Area of the shared side
    double Bn;                             ///< Elastic normal constant for the cohesion
    double Bt;                             ///< Elastic tangential constant for the cohesion
    double Bm;                             ///< Elastic for the cohesion torque
    double Gn;                             ///< Dissipative normal constant
    double Gt;                             ///< Dissipative tangential constant
    double L0;                             ///< Equilibrium distance
    double An;                             ///< Angular displacement
    double eps;                            ///< Maximun strain before fracture
    double eta;                            ///< ratio between normal and tangential breaking threshold
    bool   valid;                          ///< Check if the bound has not been broken
    double s1,t1;                          ///< Planar coordinates for face F1
    double s2,t2;                          ///< Planar coordinates for face F2
    Vec3_t Fnet;                           ///< Net Normal force excerted by the interacton
    Vec3_t Ftnet;                          ///< Net Tangential force excerted by the interacton
    Vec3_t xnet;                           ///< Position where the force is applied

};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// Collision interacton

inline CInteracton::CInteracton (Particle * Pt1, Particle * Pt2, size_t contaclaw)
{
    First           = true;
    BothFree        = Pt1->IsFree()&&Pt2->IsFree();
    P1              = Pt1;
    P2              = Pt2;
    I1              = P1->Index;
    I2              = P2->Index;
    //double r1       = pow(P1->Props.V,1.0/3.0);
    //double r2       = pow(P2->Props.V,1.0/3.0);
    //Kn              = (r1+r2)*ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
    //Kt              = (r1+r2)*ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
    if (contaclaw==0)
    {
        Kn              = ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
        Kt              = ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
    }
    else if(contaclaw==1)
    {
        if (Pt1->Verts.Size()==1)
        {
            Kn              = ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn)*Pt1->Props.R;
            Kt              = ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt)*Pt1->Props.R;
        } 
        else
        {
            Kn              = ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
            Kt              = ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
        }
    }
    double me       = ReducedValue(Pt1->Props.m ,Pt2->Props.m );
    Gn              = 2*ReducedValue(Pt1->Props.Gn,Pt2->Props.Gn);
    Gt              = 2*ReducedValue(Pt1->Props.Gt,Pt2->Props.Gt);
    if (Gn < -0.001)
    {
        if (fabs(Gn)>1.0) throw new Fatal("CInteracton the restitution coefficient is greater than 1");
        Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
        Gt = 0.0;
    }
    Gn *= me;
    Gt *= me;
    //Mu              = 2*ReducedValue(Pt1->Props.Mu,Pt2->Props.Mu);
    if (Pt1->Props.Mu>1.0e-12&&Pt2->Props.Mu>1.0e-12)
    {
        if (!Pt1->IsFree())      Mu = Pt1->Props.Mu;
        else if (!Pt2->IsFree()) Mu = Pt2->Props.Mu;
        else                     Mu = std::max(Pt1->Props.Mu,Pt2->Props.Mu);
    }
    else 
    {
        Mu          = 0.0;
    }
    //if (Pt2->Tag<-1) std::cout << Mu << " " << Pt1->Tag << " " << Pt2->Tag << std::endl;
    CalcForce(0.0);

}

inline bool CInteracton::UpdateContacts (double alpha, Vec3_t const & Per, size_t const iter)
{
    Lee.Resize(0);
    Lvf.Resize(0);
    Lfv.Resize(0);
    Lvt.Resize(0);
    Ltv.Resize(0);
    Lvc.Resize(0);
    Lcv.Resize(0);
    Vec3_t Pert = OrthoSys::O;
    if (BothFree) Pert = Per;
    if (Distance(P1->x,P2->x,Pert)<=P1->Dmax+P2->Dmax+2*alpha)
    {
        _update_contacts (P1->Edges     ,P2->Edges     ,Lee,alpha,Pert);
        _update_contacts (P1->Verts     ,P2->Faces     ,Lvf,alpha,Pert);
        _update_contacts (P1->Faces     ,P2->Verts     ,Lfv,alpha,Pert);
        _update_contacts (P1->Verts     ,P2->Tori      ,Lvt,alpha,Pert);
        _update_contacts (P1->Tori      ,P2->Verts     ,Ltv,alpha,Pert);
        _update_contacts (P1->Verts     ,P2->Cylinders ,Lvc,alpha,Pert);
        _update_contacts (P1->Cylinders ,P2->Verts     ,Lcv,alpha,Pert);
        if (Lee.Size()>0||Lvf.Size()>0||Lfv.Size()>0||Ltv.Size()>0||Lvt.Size()>0||Lcv.Size()>0||Lvc.Size()>0) return true;
        //if (Lee.Size()>0||Lvf.Size()>0||Lfv.Size()>0) return true;
        else return false;
    }
    else return false;
}

inline bool CInteracton::CalcForce (double dt, Vec3_t const & Per, size_t const iter, size_t const contactlaw)
{
    bool   overlap = false;
    Epot   = 0.0;
    dEvis  = 0.0;
    dEfric = 0.0;
    Nc     = 0;
    Nsc    = 0;
    Nr     = 0;
    Fnet   = OrthoSys::O;
    Ftnet  = OrthoSys::O;
    //Xc     = OrthoSys::O;
    F1     = OrthoSys::O;
    F2     = OrthoSys::O;
    T1     = OrthoSys::O;
    T2     = OrthoSys::O;
    Vec3_t Pert;
    if (BothFree) // I had to do this due to the two definitions for BranchVec function
    {
        BranchVec(P2->x,P1->x,Branch,Per);
        Pert = Per;
    }
    else
    {
        Branch = P1->x-P2->x;
        Pert = OrthoSys::O;
    }
    if (norm(Branch) > P1->Dmax + P2->Dmax) return false;
    //std::cout << "Lee " << Lee.Size() << std::endl;
    //std::cout << "Lvf " << Lvf.Size() << std::endl;
    //std::cout << "Lfv " << Lfv.Size() << std::endl;
    //std::cout << Pert << std::endl;
    //std::cout << I1 << std::endl;
    //std::cout << I2 << std::endl;
    //std::cout << Lvf[0].first << std::endl;
    //std::cout << Lvf[0].second << std::endl;
    if (_update_disp_calc_force (P1->Edges     ,P2->Edges     ,Fdee,Lee,dt,Pert,iter)) overlap = true;
    if (_update_disp_calc_force (P1->Verts     ,P2->Faces     ,Fdvf,Lvf,dt,Pert,iter)) overlap = true;
    if (_update_disp_calc_force (P1->Faces     ,P2->Verts     ,Fdfv,Lfv,dt,Pert,iter)) overlap = true;
    if (_update_disp_calc_force (P1->Verts     ,P2->Tori      ,Fdvt,Lvt,dt,Per)) overlap = true;
    if (_update_disp_calc_force (P1->Tori      ,P2->Verts     ,Fdtv,Ltv,dt,Per)) overlap = true;
    if (_update_disp_calc_force (P1->Verts     ,P2->Cylinders ,Fdvc,Lvc,dt,Per)) overlap = true;
    if (_update_disp_calc_force (P1->Cylinders ,P2->Verts     ,Fdcv,Lcv,dt,Per)) overlap = true;
    //
    return overlap;
}

inline void CInteracton::UpdateParameters (size_t contactlaw)
{
    Kn              = 2*ReducedValue(P1->Props.Kn,P2->Props.Kn);
    Kt              = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
    double me       = ReducedValue(P1->Props.m ,P2->Props.m );
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn);
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt);
    if (Gn < -0.001)
    {
        if (fabs(Gn)>1.0) throw new Fatal("CInteracton the restitution coefficient is greater than 1");
        Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
        Gt = 0.0;
    }
    Gn *= me;
    Gt *= me;
    if (P1->Props.Mu>1.0e-12&&P2->Props.Mu>1.0e-12)
    {
        if (!P1->IsFree())      Mu = P1->Props.Mu;
        else if (!P2->IsFree()) Mu = P2->Props.Mu;
        else                    Mu = std::max(P1->Props.Mu,P2->Props.Mu);
    }
    else 
    {
        Mu          = 0.0;
    }
}

template<typename FeatureA_T, typename FeatureB_T>
inline bool CInteracton::_update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, ListContacts_t & L, double dt, Vec3_t const & Per, size_t const iter)
{
    // update
    for (size_t k=0; k<L.Size(); ++k)
    {
        size_t i = L[k].first;
        size_t j = L[k].second;
        if (!Overlap((*A[i]), (*B[j]),P1->Props.R,P2->Props.R,Per)) continue;
        Vec3_t xi, xf, s;
        Distance ((*A[i]), (*B[j]), xi, xf, s, Per);
        double dist  = norm(s);
        double delta = P1->Props.R + P2->Props.R - dist;
        //if (First&&delta>0.0)
        //{
            //std::cout << "Initial overlapping detected "  << std::endl;
            //std::cout << "overlap between         " << I1                << " and " << I2        <<  std::endl; 
            //std::cout << "Overlap                 " << delta             <<  std::endl; 
            //std::cout << "Particle's tags         " << P1->Tag           << " and " << P2->Tag   <<  std::endl; 
            //std::cout << "Memory address          " << P1                << " and " << P2        <<  std::endl; 
            //std::cout << "Memory address pair     " << this              << std::endl; 
            //std::cout << "Per Position particle 1 " << xi                << std::endl;
            //std::cout << "Per Position particle 2 " << xf                << std::endl;
            //std::cout << "Position particle 1     " << P1->x             << std::endl;
            //std::cout << "Position particle 2     " << P2->x             << std::endl;
            //std::cout << "Velocity particle 1     " << P1->v             << std::endl;
            //std::cout << "Velocity particle 2     " << P2->v             << std::endl;
            //std::cout << "Ang Velocity particle 1 " << P1->w             << std::endl;
            //std::cout << "Ang Velocity particle 2 " << P2->w             << std::endl;
            //std::cout << "Mass particle 1         " << P1->Props.m       << std::endl;
            //std::cout << "Mass particle 2         " << P2->Props.m       << std::endl;
            //std::cout << "Diameter particle 1     " << P1->Dmax          << std::endl;
            //std::cout << "Diameter particle 2     " << P2->Dmax          << std::endl;
            //std::cout << "Sradius particle 1      " << P1->Props.R       << std::endl;
            //std::cout << "Sradius particle 2      " << P2->Props.R       << std::endl;
            //std::cout << "Number of faces  1      " << P1->Faces.Size()  << std::endl;
            //std::cout << "Number of faces  2      " << P2->Faces.Size()  << std::endl;
            //return true;
        //}
        First = false;
        if (delta>0)
        {
#ifdef USE_CHECK_OVERLAP
            //if (delta > 0.4*std::min(P1->Props.R,P2->Props.R))
            if (delta > 0.4*(P1->Props.R+P2->Props.R))
            {
                std::cout << std::endl; 
                std::cout << "Maximun overlap between " << P1->Index         << " and " << P2->Index <<  std::endl; 
                std::cout << "Is the first collision? " << First             << std::endl; 
                std::cout << "Overlap                 " << delta             <<  std::endl; 
                std::cout << "Particle's tags         " << P1->Tag           << " and " << P2->Tag   <<  std::endl; 
                std::cout << "Memory address          " << P1                << " and " << P2        <<  std::endl; 
                std::cout << "Position particle 1     " << P1->x             << std::endl;
                std::cout << "Position particle 2     " << P2->x             << std::endl;
                std::cout << "Velocity particle 1     " << P1->v             << std::endl;
                std::cout << "Velocity particle 2     " << P2->v             << std::endl;
                std::cout << "Ang Velocity particle 1 " << P1->w             << std::endl;
                std::cout << "Ang Velocity particle 2 " << P2->w             << std::endl;
                std::cout << "Mass particle 1         " << P1->Props.m       << std::endl;
                std::cout << "Mass particle 2         " << P2->Props.m       << std::endl;
                std::cout << "Diameter particle 1     " << P1->Dmax          << std::endl;
                std::cout << "Diameter particle 2     " << P2->Dmax          << std::endl;
                std::cout << "Sradius particle 1      " << P1->Props.R       << std::endl;
                std::cout << "Sradius particle 2      " << P2->Props.R       << std::endl;
                std::cout << "Number of faces  1      " << P1->Faces.Size()  << std::endl;
                std::cout << "Number of faces  2      " << P2->Faces.Size()  << std::endl;
                P1->Tag = 10000;
                P2->Tag = 10000;
                return true;
                //throw new Fatal("Interacton::_update_disp_calc_force: Maximun overlap detected between particles %d(%d) and %d(%d)",P1->Index,P1->Tag,P2->Index,P2->Tag);
            }
#endif
            // Count a contact
            Nc++;

            // update force
            Vec3_t n = s/dist;
            double d = (P1->Props.R*P1->Props.R-P2->Props.R*P2->Props.R+dist*dist)/(2*dist);
            Vec3_t x1c = xi+n*d;
            Vec3_t x2c = xf-n*(dist-d);
            //Xc += x;
            Vec3_t t1,t2,x1,x2;
            Rotation(P1->w,P1->Q,t1);
            Rotation(P2->w,P2->Q,t2);
            x1 = x1c - P1->x;
            x2 = x2c - P2->x;
            Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
            Vec3_t vt = vrel - dot(n,vrel)*n;
            Fn = Kn*delta*n;
            Fnet += Fn;
            //std::pair<int,int> p;
            //p = std::make_pair(i,j);
            size_t p = HashFunction(i,j);
            if (FMap.count(p)==0) FMap[p] = OrthoSys::O;
            FMap[p] += vt*dt;
            FMap[p] -= dot(FMap[p],n)*n;
            Vec3_t tan = FMap[p];
            if (norm(tan)>0.0) tan/=norm(tan);
            if (norm(FMap[p])>Mu*norm(Fn)/Kt)
            {
                // Count a sliding contact
                Nsc++;
                FMap[p] = Mu*norm(Fn)/Kt*tan;
                dEfric += Kt*dot(FMap[p],vt)*dt;
            }
            Ftnet += Kt*FMap[p];
            Vec3_t F = Fn + Kt*FMap[p] + Gn*dot(n,vrel)*n + Gt*vt;



            if (dot(F,n)<0) F-=dot(F,n)*n;
            //P1->F    += -F;
            //P2->F    +=  F;
            //P1->Comp += norm(Fn);
            //P2->Comp += norm(Fn);
            F1   += -F;
            F2   +=  F;
            dEvis += (Gn*dot(vrel-vt,vrel-vt)+Gt*dot(vt,vt))*dt;
            // torque
            Vec3_t T, Tt;
            Tt = cross (x1,F);
            Quaternion_t q;
            Conjugate (P1->Q,q);
            Rotation  (Tt,q,T);
            //P1->T -= T;
            T1 -= T;
            Tt = cross (x2,F);
            Conjugate (P2->Q,q);
            Rotation  (Tt,q,T);
            //P2->T += T;
            T2 += T;
            //Transfering the branch vector information
            //Vec3_t nor = n;
            //Vec3_t dif = P2->x - P1->x;
            //if (P1->IsFree()&&P2->IsFree())
            //{
                //for (size_t m=0;m<3;m++)
                //{
                    //for (size_t n=0;n<3;n++)
                    //{
                        //P1->M(m,n)  -= F(m)*x1(n);
                        //P2->M(m,n)  += F(m)*x2(n);
                        //P1->B(m,n)  += nor(m)*nor(n);
                        //P2->B(m,n)  += nor(m)*nor(n);
                        //this->B(m,n) = nor(m)*nor(n);
                    //}
                //}
            //}
    //std::cout << P1->Index << " " << P2->Index << std::endl;

            // potential energy
            Epot += 0.5*Kn*delta*delta+0.5*Kt*dot(FMap[p],FMap[p]);
        }
    }
    return false;
}

template<typename FeatureA_T, typename FeatureB_T>
inline void CInteracton::_update_contacts (FeatureA_T & A, FeatureB_T & B, ListContacts_t & L, double alpha, Vec3_t const & Per, size_t const iter)
{
    for (size_t i=0; i<A.Size(); ++i)
    for (size_t j=0; j<B.Size(); ++j)
    {
        if (!Overlap ((*A[i]), (*B[j]),P1->Props.R + alpha,P2->Props.R + alpha,Per)) continue;
        if (Distance ((*A[i]), (*B[j]),Per)<=P1->Props.R+P2->Props.R+2*alpha)
        {
            std::pair<int,int> p;
            p = std::make_pair(i,j);
            L.Push(p);
        }
    }
}

//Collision interacton for spheres

inline CInteractonSphere::CInteractonSphere (Particle * Pt1, Particle * Pt2, size_t contactlaw)
{
    First           = true;
    BothFree        = Pt1->IsFree()&&Pt2->IsFree();
    P1   = Pt1;
    P2   = Pt2;
    I1 = P1->Index;
    I2 = P2->Index;
    if (I1==I2) throw new Fatal("CInteractonSphere indexes of the particles are the same");
    double me       = ReducedValue(Pt1->Props.m ,Pt2->Props.m );
    Gn              = 2*ReducedValue(Pt1->Props.Gn,Pt2->Props.Gn);
    Gt              = 2*ReducedValue(Pt1->Props.Gt,Pt2->Props.Gt);
    if (contactlaw==0) 
    {
        Kn   = 2*ReducedValue(Pt1->Props.Kn,Pt2->Props.Kn);
        Kt   = 2*ReducedValue(Pt1->Props.Kt,Pt2->Props.Kt);
        if (Gn < 0.0)
        {
            if (fabs(Gn)>1.0) throw new Fatal("CInteractonSphere the restitution coefficient is greater than 1");
            Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
            Gt = 0.0;
        }
        Gn *= me;
        Gt *= me;
    }

    else if (contactlaw==1)
    {
        double Re  = ReducedValue(Pt1->Props.R,Pt2->Props.R);
        double nu1 = 0.5*Pt1->Props.Kn/Pt1->Props.Kt-1.0;
        double nu2 = 0.5*Pt2->Props.Kn/Pt2->Props.Kt-1.0;
        double Ye  = 1.0/((1.0-nu1*nu1)/Pt1->Props.Kn + (1.0-nu2*nu2)/Pt2->Props.Kn);
        double Ge  = 1.0/(2.0*(2.0-nu1)*(1.0+nu1)/Pt1->Props.Kn + 2.0*(2.0-nu2)*(1.0+nu2)/Pt2->Props.Kn);
        Kn  = 4.0/3.0*sqrt(Re)*Ye;
        Kt  = 8.0*sqrt(Re)*Ge;

        if (Gn < 0.0)
        {
            double b = sqrt(pow(log(-Gn),2.0)/(M_PI*M_PI+pow(log(-Gn),2.0)));
            Gn = 2.0*sqrt(5.0/6.0)*b*sqrt(me*2.0*Ye*sqrt(Re));
            Gt = 2.0*sqrt(5.0/6.0)*b*sqrt(me*8.0*Ge*sqrt(Re));
        }
    }

    //Mu   = 2*ReducedValue(Pt1->Props.Mu,Pt2->Props.Mu);
    if (Pt1->Props.Mu>1.0e-12&&Pt2->Props.Mu>1.0e-12)
    {
        if (!Pt1->IsFree())      Mu = Pt1->Props.Mu;
        else if (!Pt2->IsFree()) Mu = Pt2->Props.Mu;
        else                     Mu = std::max(Pt1->Props.Mu,Pt2->Props.Mu);
    }
    else 
    {
        Mu          = 0.0;
    }
    beta = 2*ReducedValue(Pt1->Props.Beta,Pt2->Props.Beta);
    eta  = 2*ReducedValue(Pt1->Props.Eta,Pt2->Props.Eta);
    Bn   = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn);
    Bt   = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt);
    eps  = 2*ReducedValue(P1->Props.eps,P2->Props.eps);

    if (eps<0.0)
    {
        cohesion = true;
        if (Bn<Kn) throw new Fatal("CInteractonSphere::CInteractonSphere if cohesion is applied, Bn must be larger than Kn");
    }
    else
    {
        cohesion = false;
    }
    
    Nc = 0;
    Nsc = 0;

    Epot = 0.0;
    Fdr  = OrthoSys::O;
    Fdvv = OrthoSys::O;

    //std::pair<int,int> p;
    //p = std::make_pair(0,0);
    //Lvv.Push(p);

    CalcForce(0.0);
}

inline bool CInteractonSphere::CalcForce(double dt, Vec3_t const & Per, size_t const iter, size_t contactlaw)
{
    Epot   = 0.0;
    dEvis  = 0.0;
    dEfric = 0.0;
    Nc     = 0;
    Nsc    = 0;
    Nr     = 0;
    Fnet   = OrthoSys::O;
    Ftnet  = OrthoSys::O;
    //Xc     = OrthoSys::O;
    F1     = OrthoSys::O;
    F2     = OrthoSys::O;
    T1     = OrthoSys::O;
    T2     = OrthoSys::O;


    Vec3_t xi = P1->x;
    Vec3_t xf = P2->x;
    if (BothFree) BranchVec(xf,xi,Branch,Per);
    else Branch = xi-xf;

    double dist = norm(Branch);
    double delta = P1->Props.R + P2->Props.R - dist;
    //if (First&&delta>0.0)
    //{
        //std::cout << "Initial overlapping detected "  << std::endl;
        //std::cout << "overlap between         " << I1                << " and " << I2        <<  std::endl; 
        //std::cout << "Overlap                 " << delta             <<  std::endl; 
        //std::cout << "Particle's tags         " << P1->Tag           << " and " << P2->Tag   <<  std::endl; 
        //std::cout << "Memory address          " << P1                << " and " << P2        <<  std::endl; 
        //std::cout << "Memory address pair     " << this              << std::endl; 
        //std::cout << "Per Position particle 1 " << xi                << std::endl;
        //std::cout << "Per Position particle 2 " << xf                << std::endl;
        //std::cout << "Position particle 1     " << P1->x             << std::endl;
        //std::cout << "Position particle 2     " << P2->x             << std::endl;
        //std::cout << "Velocity particle 1     " << P1->v             << std::endl;
        //std::cout << "Velocity particle 2     " << P2->v             << std::endl;
        //std::cout << "Ang Velocity particle 1 " << P1->w             << std::endl;
        //std::cout << "Ang Velocity particle 2 " << P2->w             << std::endl;
        //std::cout << "Mass particle 1         " << P1->Props.m       << std::endl;
        //std::cout << "Mass particle 2         " << P2->Props.m       << std::endl;
        //std::cout << "Diameter particle 1     " << P1->Dmax          << std::endl;
        //std::cout << "Diameter particle 2     " << P2->Dmax          << std::endl;
        //std::cout << "Sradius particle 1      " << P1->Props.R       << std::endl;
        //std::cout << "Sradius particle 2      " << P2->Props.R       << std::endl;
        //std::cout << "Number of faces  1      " << P1->Faces.Size()  << std::endl;
        //std::cout << "Number of faces  2      " << P2->Faces.Size()  << std::endl;
        //return true;
    //}

    First = false;
    if (delta>0.0)
    {
#ifdef USE_CHECK_OVERLAP
        if (delta > 0.4*(P1->Props.R+P2->Props.R))
        {
            std::cout << std::endl; 
            std::cout << "Maximun overlap between " << I1                << " and " << I2        <<  std::endl; 
            std::cout << "Is the first collision? " << First             << std::endl; 
            std::cout << "Overlap                 " << delta             <<  std::endl; 
            std::cout << "Particle's tags         " << P1->Tag           << " and " << P2->Tag   <<  std::endl; 
            std::cout << "Memory address          " << P1                << " and " << P2        <<  std::endl; 
            std::cout << "Memory address pair     " << this              << std::endl; 
            std::cout << "Per Position particle 1 " << xi                << std::endl;
            std::cout << "Per Position particle 2 " << xf                << std::endl;
            std::cout << "Position particle 1     " << P1->x             << std::endl;
            std::cout << "Position particle 2     " << P2->x             << std::endl;
            std::cout << "Velocity particle 1     " << P1->v             << std::endl;
            std::cout << "Velocity particle 2     " << P2->v             << std::endl;
            std::cout << "Ang Velocity particle 1 " << P1->w             << std::endl;
            std::cout << "Ang Velocity particle 2 " << P2->w             << std::endl;
            std::cout << "Mass particle 1         " << P1->Props.m       << std::endl;
            std::cout << "Mass particle 2         " << P2->Props.m       << std::endl;
            std::cout << "Diameter particle 1     " << P1->Dmax          << std::endl;
            std::cout << "Diameter particle 2     " << P2->Dmax          << std::endl;
            std::cout << "Sradius particle 1      " << P1->Props.R       << std::endl;
            std::cout << "Sradius particle 2      " << P2->Props.R       << std::endl;
            std::cout << "Number of faces  1      " << P1->Faces.Size()  << std::endl;
            std::cout << "Number of faces  2      " << P2->Faces.Size()  << std::endl;
            P1->Tag = 10000;
            P2->Tag = 20000;
            return true;
        }
#endif
        //Count a contact
        Nc++;

        //update force
        Vec3_t n = -Branch/dist;
        double d = (P1->Props.R*P1->Props.R-P2->Props.R*P2->Props.R+dist*dist)/(2*dist);
        Vec3_t x1c = xi+n*d;
        Vec3_t x2c = xf-n*(dist-d);
        //Xc += x;
        Vec3_t t1,t2,x1,x2;
        Rotation(P1->w,P1->Q,t1);
        Rotation(P2->w,P2->Q,t2);
        x1 = x1c- xi;
        x2 = x2c- xf;
        Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
        Vec3_t vt = vrel - dot(n,vrel)*n;
        Vec3_t F  = OrthoSys::O;
        Vec3_t Ft = OrthoSys::O;
        if (contactlaw==0)
        {
            Fn  = Kn*delta*n;
            
            Fnet += Fn;
            Fdvv += vt*dt;
            Fdvv -= dot(Fdvv,n)*n;
            Vec3_t tan = Fdvv;
            if (norm(tan)>0.0) tan/=norm(tan);
            if (norm(Fdvv)>Mu*norm(Fn)/Kt)
            {
                 //Count a sliding contact
                Nsc++;
                Fdvv = Mu*norm(Fn)/Kt*tan;
                dEfric += Kt*dot(Fdvv,vt)*dt;
            }
            Ftnet += Kt*Fdvv;

            //Calculating the rolling resistance torque
            double Kr = beta*Kt;
            Vec3_t Vr = P1->Props.R*P2->Props.R*cross(Vec3_t(t1 - t2),n)/(P1->Props.R+P2->Props.R);
            Fdr += Kr*Vr*dt;
            Fdr -= dot(Fdr,n)*n;
            
            tan = Fdr;
            if (norm(tan)>0.0) tan/=norm(tan);
            if (norm(Fdr)>eta*Mu*norm(Fn))
            {
                Fdr = eta*Mu*norm(Fn)*tan;
                Nr++;
            }

            Ft = -Fdr;
            
            F = Fn + Kt*Fdvv + Gn*dot(n,vrel)*n + Gt*vt;
            dEvis += (Gn*dot(vrel-vt,vrel-vt)+Gt*dot(vt,vt))*dt;
        }
        else if (contactlaw==1)
        {
            Fn  = Kn*sqrt(delta)*delta*n;

            Fnet += Fn;
            Fdvv += vt*dt;
            Fdvv -= dot(Fdvv,n)*n;
            Vec3_t tan = Fdvv;
            if (norm(tan)>0.0) tan/=norm(tan);
            if (norm(Fdvv)>Mu*norm(Fn)/(Kt*sqrt(delta)))
            {
                 //Count a sliding contact
                Nsc++;
                Fdvv = Mu*norm(Fn)/(Kt*sqrt(delta))*tan;
                dEfric += Kt*sqrt(delta)*dot(Fdvv,vt)*dt;
            }
            Ftnet += Kt*sqrt(delta)*Fdvv;

            //Calculating the rolling resistance torque
            double Kr = beta*Kt*sqrt(delta);
            Vec3_t Vr = P1->Props.R*P2->Props.R*cross(Vec3_t(t1 - t2),n)/(P1->Props.R+P2->Props.R);
            Fdr += Kr*Vr*dt;
            Fdr -= dot(Fdr,n)*n;
            
            tan = Fdr;
            if (norm(tan)>0.0) tan/=norm(tan);
            if (norm(Fdr)>eta*Mu*norm(Fn))
            {
                Fdr = eta*Mu*norm(Fn)*tan;
                Nr++;
            }

            Ft = -Fdr;
            
            F = Fn + Kt*sqrt(delta)*Fdvv + Gn*sqrt(sqrt(delta))*dot(n,vrel)*n + Gt*sqrt(sqrt(delta))*vt;
            dEvis += (Gn*sqrt(sqrt(delta))*dot(vrel-vt,vrel-vt)+Gt*sqrt(sqrt(delta))*dot(vt,vt))*dt;
        }
        

        F1   += -F;
        F2   +=  F;
        //torque
        Vec3_t T, Tt;
        Tt = cross (x1,F) - P1->Props.R*cross(n,Ft);
        Quaternion_t q;
        Conjugate (P1->Q,q);
        Rotation  (Tt,q,T);
        T1 -= T;
        Tt = cross (x2,F) - P2->Props.R*cross(n,Ft);
        Conjugate (P2->Q,q);
        Rotation  (Tt,q,T);
        T2 += T;

        Fn        = Fdr;
    }

    //If there is at least a contact, increase the coordination number of the particles
    //if (Nc>0) 
    //{
        //P1->Cn++;
        //P2->Cn++;
        //if (!P1->IsFree()) P2->Bdry = true;
        //if (!P2->IsFree()) P1->Bdry = true;
    //}
    //return overlap;
    return false;
}

inline bool CInteractonSphere::UpdateContacts (double alpha, Vec3_t const & Per, size_t const iter)
{
    if (Distance(P1->x,P2->x,Per)<=P1->Dmax+P2->Dmax+2*alpha) return true;
    else return false;
}

inline void CInteractonSphere::UpdateParameters (size_t contactlaw)
{
    Kn   = 2*ReducedValue(P1->Props.Kn,P2->Props.Kn);
    Kt   = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
    double me       = ReducedValue(P1->Props.m ,P2->Props.m );
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn);
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt);
    if (contactlaw==0) 
    {
        Kn   = 2*ReducedValue(P1->Props.Kn,P2->Props.Kn);
        Kt   = 2*ReducedValue(P1->Props.Kt,P2->Props.Kt);
        if (Gn < 0.0)
        {
            if (fabs(Gn)>1.0) throw new Fatal("CInteractonSphere the restitution coefficient is greater than 1");
            Gn = 2.0*sqrt((pow(log(-Gn),2.0)*(Kn/me))/(M_PI*M_PI+pow(log(-Gn),2.0)));
            Gt = 0.0;
        }
        Gn *= me;
        Gt *= me;
    }

    else if (contactlaw==1)
    {
        double Re  = ReducedValue(P1->Props.R,P2->Props.R);
        double nu1 = 0.5*P1->Props.Kn/P1->Props.Kt-1.0;
        double nu2 = 0.5*P2->Props.Kn/P2->Props.Kt-1.0;
        double Ye  = 1.0/((1.0-nu1*nu1)/P1->Props.Kn + (1.0-nu2*nu2)/P2->Props.Kn);
        double Ge  = 1.0/(2.0*(2.0-nu1)*(1.0+nu1)/P1->Props.Kn + 2.0*(2.0-nu2)*(1.0+nu2)/P2->Props.Kn);
        Kn  = 4.0/3.0*sqrt(Re)*Ye;
        Kt  = 8.0*sqrt(Re)*Ge;

        if (Gn < 0.0)
        {
            double b = sqrt(pow(log(-Gn),2.0)/(M_PI*M_PI+pow(log(-Gn),2.0)));
            Gn = 2.0*sqrt(5.0/6.0)*b*sqrt(me*2.0*Ye*sqrt(Re));
            Gt = 2.0*sqrt(5.0/6.0)*b*sqrt(me*8.0*Ge*sqrt(Re));
        }
    }

    if (P1->Props.Mu>1.0e-12&&P2->Props.Mu>1.0e-12)
    {
        if (!P1->IsFree())      Mu = P1->Props.Mu;
        else if (!P2->IsFree()) Mu = P2->Props.Mu;
        else                    Mu = std::max(P1->Props.Mu,P2->Props.Mu);
    }
    else 
    {
        Mu          = 0.0;
    }
    beta = 2*ReducedValue(P1->Props.Beta,P2->Props.Beta);
    eta  = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);
    Bn   = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn);
    Bt   = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt);
    eps  = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
}

//Cohesion interacton

inline BInteracton::BInteracton (Particle * Pt1, Particle * Pt2, size_t Fi1, size_t Fi2)
{
    P1              = Pt1;
    P2              = Pt2;
    IF1             = Fi1;
    if (Fi2>=P2->Faces.Size())
    {
        IF2 = P2->Faces.Size()-1;
        Area            = P1->Faces[IF1]->Area();
    }
    else
    {
        IF2 = Fi2;
        Area            = 0.5*(P1->Faces[IF1]->Area()+P2->Faces[IF2]->Area());
    }
    I1              = P1->Index;
    I2              = P2->Index;
    Bn              = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn)*Area;
    Bt              = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt)*Area;
    Bm              = 2*ReducedValue(P1->Props.Bm,P2->Props.Bm)*Area;
    double me       = ReducedValue(Pt1->Props.m ,Pt2->Props.m );
    Gn              = 2*ReducedValue(Pt1->Props.Gn,Pt2->Props.Gn)*me;
    Gt              = 2*ReducedValue(Pt1->Props.Gt,Pt2->Props.Gt)*me;
    eps             = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
    eta             = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);

    Vec3_t n1,n2;
    P1->Faces[IF1]->Normal(n1);
    P2->Faces[IF2]->Normal(n2);
    Vec3_t c1,c2;
    An         = 0.0;
    valid      = true;
    P1->Faces[IF1]->Centroid(c1);
    P2->Faces[IF2]->Centroid(c2);
    L0         = dot(n1,c2-c1);
    if (Fi2>=P2->Faces.Size()) c2 = c1;
    Vec3_t V   = 0.5*(c1+c2);

    Face * F   = P1->Faces[IF1];
    double a   = dot(*F->Edges[0]->X0-V , F->Edges[0]->dL);
    double b   = dot(F->Edges[0]->dL    , F->Edges[0]->dL);
    double c   = dot(F->Edges[0]->dL    , F->Edges[1]->dL);
    double d   = dot(*F->Edges[0]->X0-V , F->Edges[1]->dL);
    double f   = dot(F->Edges[1]->dL    , F->Edges[1]->dL);
    s1         = (c*d-a*f)/(b*f-c*c);
    t1         = (a*c-b*d)/(b*f-c*c);
    
    Vec3_t p1  = -s1*F->Edges[0]->dL - t1*F->Edges[1]->dL;

    F          = P2->Faces[IF2];
    a          = dot(*F->Edges[0]->X0-V , F->Edges[0]->dL);
    b          = dot(F->Edges[0]->dL    , F->Edges[0]->dL);
    c          = dot(F->Edges[0]->dL    , F->Edges[1]->dL);
    d          = dot(*F->Edges[0]->X0-V , F->Edges[1]->dL);
    f          = dot(F->Edges[1]->dL    , F->Edges[1]->dL);
    s2         = (c*d-a*f)/(b*f-c*c);
    t2         = (a*c-b*d)/(b*f-c*c);

    Vec3_t p2  = -s2*F->Edges[0]->dL - t2*F->Edges[1]->dL;
    
    double arg = dot(p1,p2)/(norm(p1)*norm(p2));
    if (arg> 1.0) arg= 1.0;
    if (arg<-1.0) arg=-1.0;
    An         = acos(arg);
    if (dot(cross(p1,p2),n1)<0) An = -An;
}

inline bool BInteracton::UpdateContacts (double alpha, Vec3_t const & Per, size_t const iter)
{
    return valid;
}

inline bool BInteracton::CalcForce(double dt, Vec3_t const & Per, size_t const iter, size_t const contactlaw)
{
    F1 = OrthoSys::O;
    F2 = OrthoSys::O;
    T1 = OrthoSys::O;
    T2 = OrthoSys::O;
    if (valid)
    {
        // Calculate the normal vector and centroid of the contact face
        Vec3_t n1,n2,n;
        P1->Faces[IF1]->Normal(n1);
        P2->Faces[IF2]->Normal(n2);
        n = 0.5*(n1-n2);
        n/=norm(n);

        Face * F    = P1->Faces[IF1];
        Vec3_t pro1 = *F->Edges[0]->X0 + s1*F->Edges[0]->dL + t1*F->Edges[1]->dL;
        Vec3_t p1   = -s1*F->Edges[0]->dL - t1*F->Edges[1]->dL;
        F           = P2->Faces[IF2];
        Vec3_t pro2 = *F->Edges[0]->X0 + s2*F->Edges[0]->dL + t2*F->Edges[1]->dL;
        Vec3_t p2   = -s2*F->Edges[0]->dL - t2*F->Edges[1]->dL;

        // force
        double delta = (dot(pro2-pro1,n)-L0)/L0;
        if (delta<0.0) delta = 0.0;

        xnet        = 0.5*(pro1+pro2);
        Vec3_t t1,t2,x1,x2;
        Rotation(P1->w,P1->Q,t1);
        Rotation(P2->w,P2->Q,t2);
        x1 = xnet - P1->x;
        x2 = xnet - P2->x;
        Vec3_t td   = pro2-pro1-dot(pro2-pro1,n)*n;
        //Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
        //Vec3_t vt   = vrel - dot(n,vrel)*n;
        //Fnet        = -Bn*(delta)*n - Gn*dot(n,vrel)*n;
        //Ftnet       = -Bt*td/L0 - Gt*vt;
        Fnet        = -Bn*(delta)*n;
        Ftnet       = -Bt*td/L0;

        //Adding forces and torques
        //P1->F -= Fnet;
        //P2->F += Fnet;
        Vec3_t FT = Fnet + Ftnet;
        F1 -= FT;
        F2 += FT;
        Vec3_t T, Tt;
        Tt = cross (x1,FT);
        Quaternion_t q1;
        Conjugate (P1->Q,q1);
        Rotation  (Tt,q1,T);
        //P1->T -= T;
        T1 -= T;
        Tt = cross (x2,FT);
        Quaternion_t q2;
        Conjugate (P2->Q,q2);
        Rotation  (Tt,q2,T);
        //P2->T += T;
        T2 += T;

        //Torque
        double arg = dot(p1,p2)/(norm(p1)*norm(p2));
        if (arg> 1.0) arg= 1.0;
        if (arg<-1.0) arg=-1.0;
        double Ant = acos(arg);
        if (dot(cross(p1,p2),n)<0) Ant = -Ant;
        Tt         = -Bm*(Ant-An)*n/L0;
        Rotation  (Tt,q1,T);
        //P1->T -= T;
        T1 -= T;
        Rotation  (Tt,q2,T);
        //P2->T += T;
        T2 += T;

        //Breaking point
        if ((fabs(delta)/(eta*eps)+norm(td)/(L0*eps)+0.0*fabs(Ant-An)*L0>1.0)&&(eps>=0.0))
        //double en = fabs(delta)/(eta*eps);
        //double et = norm(td)/(L0*eps);
        //if ((en*en + et*et > 1.0)&&(eps>=0.0))
        {
            valid = false;
        }
    }
    return false;
}

inline void BInteracton::UpdateParameters (size_t contactlaw)
{
    Bn              = 2*ReducedValue(P1->Props.Bn,P2->Props.Bn)*Area;
    Bt              = 2*ReducedValue(P1->Props.Bt,P2->Props.Bt)*Area;
    Bm              = 2*ReducedValue(P1->Props.Bm,P2->Props.Bm)*Area;
    double me       = ReducedValue(P1->Props.m ,P2->Props.m );
    Gn              = 2*ReducedValue(P1->Props.Gn,P2->Props.Gn)*me;
    Gt              = 2*ReducedValue(P1->Props.Gt,P2->Props.Gt)*me;
    eps             = 2*ReducedValue(P1->Props.eps,P2->Props.eps);
    eta             = 2*ReducedValue(P1->Props.Eta,P2->Props.Eta);
    I1              = P1->Index;
    I2              = P2->Index;
}

#ifdef USE_CUDA
//////////////////////////////// CUDA IMPLEMENTATION /////////////////////////////
struct InteractonCU
{
    bool           BothFree;  ///< Are both particles free?
    real3          Fnf;       ///< Prescribed normal force
    real3          Ftf;       ///< Prescribed tangential force
    //real           Dmax1;     ///< Circumscribed diameter of particle 1
    //real           Dmax2;     ///< Circumscribed diameter of particle 2
    real           Kn;        ///< Normal stiffness
    real           Kt;        ///< Tengential stiffness
    real           Gn;        ///< Normal viscous coefficient
    real           Gt;        ///< Tangential viscous coefficient
    real           Mu;        ///< Microscpic coefficient of friction
    real           Beta;      ///< Rolling stiffness coefficient
    real           Eta;       ///< Plastic moment coefficient        
};

struct ComInteractonCU
{
    size_t         I1;        ///< Index of the first particle
    size_t         I2;        ///< Index of the second particle
    real3          Fnnet;     ///< Normal force between particles
    real3          Ftnet;     ///< Net tangential force
};

struct DynInteractonCU
{
    real           Dmax1;     ///< Circumscribed diameter of feature 1
    real           Dmax2;     ///< Circumscribed diameter of feature 2
    size_t         Idx;       ///< index of common interacton
    size_t         IF1;       ///< Index of the first feature
    size_t         IF2;       ///< Index of the second feature
    real3          F;         ///< Provisional force  
    real3          Fn;        ///< Normal force between elements
    real3          Ft;        ///< Net tangential force
    real3          Fr;        ///< Rolling resistance tangential force
};

#endif

}
#endif //  MECHSYS_DEM_INTERACTON_H
