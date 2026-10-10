/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                  *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// Save/Load round trip on the test_01 scene (the sliding block).
//
//   dom:  run K steps, Save, then continue in memory for M more steps
//   ld :  Load the same file into a FRESH domain, continue for the same M steps
//
// and the two are then compared field by field.  Run with one thread the whole thing is
// serial and deterministic, so the comparison must be exact -- anything that differs is
// something Save or Load failed to carry.  (On the GPU the contact forces accumulate with
// atomicAdd, so the two will not agree to the last bit; use Nproc=1 on the CPU path.)
//
// Usage:  test_saveload [Nproc] [K] [M]

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

#include <cstdio>
#include <cmath>
#include <algorithm>

using Util::PI;
using DEM::Domain;
using DEM::Particle;

int main(int argc, char ** argv) try
{
    size_t Nproc = (argc>1)?atoi(argv[1]):1;     // 1 = serial and deterministic
    size_t K     = (argc>2)?atoi(argv[2]):300;   // steps before saving
    size_t M     = (argc>3)?atoi(argv[3]):300;   // steps after loading
    double dt    = 1.0e-5;

    printf("\n=== test_saveload ===\n");
    printf("  Nproc=%zu   save after %zu steps (t=%g), then %zu more (t=%g)\n",
           Nproc,K,K*dt,M,(K+M)*dt);

    // ---------------------------------------------------------------------------------
    // The scene: test_01's sliding block, plus a pair of stacked spheres away from the
    // cube.  The spheres exercise the SPHERE friction maps, CInteractonSphere::Fdvv
    // (tangential displacement) and Fdr (rolling).  Two are needed because the library
    // only builds a CInteractonSphere when BOTH particles are single-vertex spheres; a
    // sphere against the cube or the floor gives an ordinary CInteracton and would fill
    // Fdee/Fdvf/Fdfv instead.  The upper one slides and rolls on the lower one.
    //
    // Heights matter: AddPlane puts the floor's top surface at z=+0.1 (note the cube's
    // bottom sits at 0.2), so a sphere of radius 0.3 rests with its centre at 0.40.
    // Starting it lower embeds it in the floor and the Kn=1e8 contact fires it away.
    // ---------------------------------------------------------------------------------
    Domain dom;
    dom.AddCube  (-1,Vec3_t(-5.0,0.0,0.7),0.1,1.0,1.0,PI/2.0,&OrthoSys::e0);  // the block
    dom.AddPlane (-2,OrthoSys::O,0.1,100.0,100.0,1.0);                        // the floor
    dom.AddSphere(-3,Vec3_t(2.0,0.0,0.40),0.3,1.0);                           // lower sphere
    dom.AddSphere(-4,Vec3_t(2.0,0.0,1.00),0.3,1.0);                           // upper sphere

    dom.GetParticle(-2)->FixVeloc();                        // the floor is fixed
    dom.GetParticle(-1)->v = Vec3_t(10.0,0.0,0.0);          // the block slides
    dom.GetParticle(-4)->v = Vec3_t( 1.0,0.0,0.0);          // the top sphere slides ...
    dom.GetParticle(-4)->w = Vec3_t( 0.0,5.0,0.0);          // ... and rolls

    dom.Initialize(dt);
    for (size_t i=0;i<dom.Particles.Size();i++)              // gravity on the free ones
        if (dom.Particles[i]->IsFree())
            dom.Particles[i]->Ff = Vec3_t(0.0,0.0,-dom.Particles[i]->Props.m*9.8);

    Dict B;
    B.Set(-1,"Gt Kn Kt",0.0,1.0e8,5.0e7);
    B.Set(-2,"Gt Kn Kt",0.0,1.0e8,5.0e7);
    B.Set(-3,"Gt Kn Kt",0.0,1.0e8,5.0e7);
    B.Set(-4,"Gt Kn Kt",0.0,1.0e8,5.0e7);
    dom.SetProps(B);

    dom.FricCoeff[std::pair<int,int>(-1,-2)] = 0.1;          // per-tag friction overrides
    dom.FricCoeff[std::pair<int,int>(-3,-4)] = 0.1;

    // ---------------------------------------------------------------------------------
    // Run, save, and load into a fresh domain.  ld must be left EMPTY: Load appends its
    // particles to whatever is already there.
    // ---------------------------------------------------------------------------------
    dom.Solve(K*dt,dt,K*dt,NULL,NULL,NULL,false,Nproc);
    dom.Save ("test_saveload");

    Domain ld;
    ld.Load("test_saveload");

    // ---- the two must agree exactly going in ----------------------------------------
    double dx=0.0,dv=0.0,dxb=0.0,dw=0.0,dq=0.0;
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        Particle * a = dom.Particles[i];
        Particle * b = ld .Particles[i];
        for (int k=0;k<3;k++)
        {
            dx  = std::max(dx ,fabs(a->x (k)-b->x (k)));
            dv  = std::max(dv ,fabs(a->v (k)-b->v (k)));
            dxb = std::max(dxb,fabs(a->xb(k)-b->xb(k)));
            dw  = std::max(dw ,fabs(a->w (k)-b->w (k)));
        }
        for (int k=0;k<4;k++) dq = std::max(dq,fabs(a->Q(k)-b->Q(k)));
    }
    bool same = (dx==0.0&&dv==0.0&&dxb==0.0&&dw==0.0&&dq==0.0);
    printf("\n  at the checkpoint   max|dx|=%.3e max|dv|=%.3e max|dxb|=%.3e max|dw|=%.3e max|dQ|=%.3e   %s\n",
           dx,dv,dxb,dw,dq,same?"BIT IDENTICAL":"*** DIFFERS ***");
    bool ok = same;

    // ---------------------------------------------------------------------------------
    // Continue both, then compare again.
    // ---------------------------------------------------------------------------------
    dom.Solve((K+M)*dt,dt,(K+M)*dt,NULL,NULL,NULL,false,Nproc);
    ld .Solve((K+M)*dt,dt,(K+M)*dt,NULL,NULL,NULL,false,Nproc);

    dx=dv=dxb=dw=dq=0.0;
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        Particle * a = dom.Particles[i];
        Particle * b = ld .Particles[i];
        for (int k=0;k<3;k++)
        {
            dx  = std::max(dx ,fabs(a->x (k)-b->x (k)));
            dv  = std::max(dv ,fabs(a->v (k)-b->v (k)));
            dxb = std::max(dxb,fabs(a->xb(k)-b->xb(k)));
            dw  = std::max(dw ,fabs(a->w (k)-b->w (k)));
        }
        for (int k=0;k<4;k++) dq = std::max(dq,fabs(a->Q(k)-b->Q(k)));
    }
    same = (dx==0.0&&dv==0.0&&dxb==0.0&&dw==0.0&&dq==0.0);
    printf("  after %zu more   max|dx|=%.3e max|dv|=%.3e max|dxb|=%.3e max|dw|=%.3e max|dQ|=%.3e   %s\n",
           M,dx,dv,dxb,dw,dq,same?"BIT IDENTICAL":"*** DIFFERS ***");
    ok = ok && same;

    // ---- how much sphere friction history is there, and is it the same on both? ------
    Domain * two[2] = {&dom,&ld};
    const char * nm[2] = {"dom","ld "};
    for (int q=0;q<2;q++)
    {
        size_t nvv=0,nrr=0;
        double mvv=0.0,mrr=0.0;
        for (auto it=two[q]->PairtoCInt.begin();it!=two[q]->PairtoCInt.end();++it)
        {
            DEM::CInteracton * Ci = it->second;
            if (two[q]->Particles[Ci->I1]->Verts.Size()!=1||two[q]->Particles[Ci->I2]->Verts.Size()!=1) continue;
            DEM::CInteractonSphere * Cs = dynamic_cast<DEM::CInteractonSphere *>(Ci);
            if (Cs==NULL) continue;
            if (norm(Cs->Fdvv)>0.0) { nvv++; mvv=std::max(mvv,norm(Cs->Fdvv)); }
            if (norm(Cs->Fdr )>0.0) { nrr++; mrr=std::max(mrr,norm(Cs->Fdr )); }
        }
        printf("  %s sphere history   Fdvv non-zero=%zu max=%.4e   Fdr non-zero=%zu max=%.4e\n",
               nm[q],nvv,mvv,nrr,mrr);
    }

    printf("\n  RESULT: %s\n", ok?"PASS - the reloaded domain is bit-for-bit the same"
                                   :"FAIL - the reloaded domain differs");
    return ok?0:1;
}
MECHSYS_CATCH
