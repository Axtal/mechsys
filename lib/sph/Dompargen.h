/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2017 Maziar Gholami                                    *
 * Copyright (C) 2020 Mario Trujillo                                    *
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

#ifndef MECHSYS_SPH_DOMAIN_PARGEN_H
#define MECHSYS_SPH_DOMAIN_PARGEN_H

// Particle generation

// Sihgle particle addition

inline void Domain::AddSphere (int Tag,Vec3_t const & X, double R, double rho)
{
    // vertices
    Array<Vec3_t> V(1);
    V[0] = X;

    // edges
    Array<Array <int> > E(0); // no edges

    // faces
    Array<Array <int> > F(0); // no faces

    // add particle
    DEMParticles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    Quaternion_t q;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    DEMParticles[DEMParticles.Size()-1]->Q          = q;
    DEMParticles[DEMParticles.Size()-1]->Props.V    = (4.0/3.0)*M_PI*R*R*R;
    DEMParticles[DEMParticles.Size()-1]->Props.m    = rho*(4.0/3.0)*M_PI*R*R*R;
    DEMParticles[DEMParticles.Size()-1]->I          = (2.0/5.0)*DEMParticles[DEMParticles.Size()-1]->Props.m*R*R;
    DEMParticles[DEMParticles.Size()-1]->x          = X;
    DEMParticles[DEMParticles.Size()-1]->Ekin       = 0.0;
    DEMParticles[DEMParticles.Size()-1]->Erot       = 0.0;
    DEMParticles[DEMParticles.Size()-1]->Dmax       = R;
    DEMParticles[DEMParticles.Size()-1]->PropsReady = true;
    DEMParticles[DEMParticles.Size()-1]->Index      = DEMParticles.Size()-1;

}

inline void Domain::AddSegment (int Tag, const Vec3_t & X0, const Vec3_t & X1, double R, double rho)
{
    Vec3_t X  = 0.5*(X0+X1);
    Vec3_t dX = X1-X0;

    AddPlane(Tag,X,R,norm(X1-X0),2.0*R,rho,M_PI/2.0,&OrthoSys::e0);

    double Angle = -acos(dot(dX,OrthoSys::e0)/norm(dX));
    Quaternion_t q;
    NormalizeRotation (Angle,OrthoSys::e2,q);

    DEMParticles[DEMParticles.Size()-1]->Rotate(q,X);
}

inline void Domain::AddPlane (int Tag, const Vec3_t & X, double R, double Lx, double Ly, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(4);
    double lx = Lx/2.0, ly = Ly/2.0;
    V[0] = -lx, -ly, 0.0;
    V[1] =  lx, -ly, 0.0;
    V[2] =  lx,  ly, 0.0;
    V[3] = -lx,  ly, 0.0;

    // edges
    Array<Array <int> > E(4);
    for (size_t i=0; i<4; ++i) E[i].Resize(2);
    E[ 0] = 0, 1;
    E[ 1] = 1, 2;
    E[ 2] = 2, 3;
    E[ 3] = 3, 0;

    // faces
    Array<Array <int> > F(1);
    F[0].Resize(4);
    F[0] = 0, 3, 2, 1;

    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = 0.;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    DEMParticles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    DEMParticles[DEMParticles.Size()-1]->Q          = q;
    DEMParticles[DEMParticles.Size()-1]->Props.V    = Lx*Ly*2*R;
    DEMParticles[DEMParticles.Size()-1]->Props.m    = rho*Lx*Ly*2*R;
    DEMParticles[DEMParticles.Size()-1]->I          = (1.0/12.0)*(Ly*Ly+4*R*R),(1.0/12.0)*(Lx*Lx+4*R*R),(1.0/12.0)*(Lx*Lx+Ly*Ly);
    DEMParticles[DEMParticles.Size()-1]->I         *= DEMParticles[DEMParticles.Size()-1]->Props.m;
    DEMParticles[DEMParticles.Size()-1]->x          = X;
    DEMParticles[DEMParticles.Size()-1]->Ekin       = 0.0;
    DEMParticles[DEMParticles.Size()-1]->Erot       = 0.0;
    DEMParticles[DEMParticles.Size()-1]->Dmax       = sqrt(Lx*Lx+Ly*Ly)+R;
    DEMParticles[DEMParticles.Size()-1]->PropsReady = true;
    DEMParticles[DEMParticles.Size()-1]->Index      = DEMParticles.Size()-1;
    // clean up
    if (!ThereisanAxis) delete Axis;
}
#endif
