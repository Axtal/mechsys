/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2016 Sergio Galindo                                    *
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

/////////////////////////////LBM OpenCL implementation////////////////////

typedef struct lbm_aux
{
    size_t     Nl;          ///< Number of Lattices
    size_t     Nneigh;      ///< Number of Neighbors
    size_t     NCPairs;     ///< Number of cell pairs
    size_t     Nx;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Ny;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Nz;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Op[27];      ///< Array with the opposite directions for bounce back calculation
    double3    C[27];       ///< Collection of discrete velocity vectors
    double     EEk[27];     ///< Dyadic product of discrete velocities for LES calculation
    double     W[27];       ///< Collection of discrete weights
    double     Tau[2];      ///< Collection of characteristic collision times
    double     G[2];        ///< Collection of cohesive constants for multiphase simulation
    double     Gs[2];       ///< Collection of cohesive constants for multiphase simulation
    double     Rhoref[2];   ///< Collection of cohesive constants for multiphase simulation
    double     Psi[2];      ///< Collection of cohesive constants for multiphase simulation
    double     Gmix;        ///< Repulsion constant for multicomponent simulation
    double     Cs;          ///< Lattice speed
    double     Sc;          ///< Smagorinsky constant
    
} d_lbm_aux;

double FeqFluid (size_t const k, double const rho, double3 const vel, global const struct lbm_aux * lbmaux)
{
    double VdotV = dot(vel,vel);
    double VdotC = dot(vel,lbmaux[0].C[k]);
    double Cs    = lbmaux[0].Cs;
    return lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
    //double Feq1temp = VdotC/Cs;
    //double Feq2temp = 0.5*VdotC*VdotC/(Cs*Cs) - 0.5*VdotV/(Cs*Cs);
    //double Feq3temp = (VdotC/(6.0*Cs))*VdotC*VdotC/(Cs*Cs) - (3.0*VdotC/(6.0*Cs))*VdotV/(Cs*Cs); 
    //return lbmaux[0].W[k]*rho*(1.0 + Feq1temp + Feq2temp + Feq3temp);
}

void Initialize(size_t iv, global double * F, global double * Rho, global double3 * Vel, double const rho, double3 const vel, global const struct lbm_aux * lbmaux)
{
    Rho[iv] = 0.0;
    Vel[iv] = (double3)(0.0,0.0,0.0);
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        F[iv*lbmaux[0].Nneigh + k] = FeqFluid(k,rho,vel,lbmaux);
        Rho[iv] += F[iv*lbmaux[0].Nneigh + k];
        Vel[iv] += F[iv*lbmaux[0].Nneigh + k]*lbmaux[0].C[k];
    }
    Vel[iv] *= lbmaux[0].Cs/Rho[iv];
}

void kernel CheckUpLoad (global struct lbm_aux * lbmaux)
{
    //printf("Nl          %lu \n",  lbmaux[0].Nl     );
    //printf("Nneigh      %lu \n",  lbmaux[0].Nneigh );
    //printf("NCP         %lu \n",  lbmaux[0].NCPairs);
    //printf("Dim      %d %lu \n",0,lbmaux[0].Nx );
    //printf("Dim      %d %lu \n",1,lbmaux[0].Ny );
    //printf("Dim      %d %lu \n",2,lbmaux[0].Nz );
    //printf("Sc          %f \n"  ,lbmaux[0].Sc );
    //printf("Tau_0        %f \n"  ,lbmaux[0].Tau[0]);
    //for (size_t i=0;i < lbmaux[0].Nneigh;i++)
    //{
        //printf("C      %d %f %f %f \n",i,lbmaux[0].C[i].x,lbmaux[0].C[i].y,lbmaux[0].C[i].z);
    //}
    //for (size_t i=0;i < lbmaux[0].Nneigh;i++)
    //{
        //printf("EEk    %d %f       \n",i,lbmaux[0].EEk[i]);
        //printf("Op    %lu %lu       \n",i,lbmaux[0].Op[i]);
    //}
    //double Feq = FeqFluid(3,1.0,(double3)(0.2,0.0,0.0),lbmaux);
    //printf(" %f \n",Feq);
}

void kernel ApplyForcesSC(global const bool * IsSolid, global double3* BForce, global const double * Rho, global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = ic/(lbmaux[0].Nx*lbmaux[0].Ny);

    for (size_t k=1;k<lbmaux[0].Nneigh;k++)
    {
        size_t inx = (size_t)((long)icx + (long)lbmaux[0].C[k].x + (long)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((long)icy + (long)lbmaux[0].C[k].y + (long)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((long)icz + (long)lbmaux[0].C[k].z + (long)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;
        
        double psic = 1.0;
        double psin = 1.0;
        double G    = lbmaux[0].G[0];
        if (!IsSolid[ic]) psic = lbmaux[0].Psi[0]*exp(-lbmaux[0].Rhoref[0]/Rho[ic]);
        if (!IsSolid[in]) psin = lbmaux[0].Psi[0]*exp(-lbmaux[0].Rhoref[0]/Rho[in]);
        else              G    = lbmaux[0].Gs[0];
        
        BForce[ic] += -G*lbmaux[0].W[k]*psic*psin*lbmaux[0].C[k];
    }
}
                                                                                                              
void kernel ApplyForcesSCMP(global const bool * IsSolid0 , global const bool * IsSolid1 ,
                            global       double3* BForce0, global       double3* BForce1,
                            global const double * Rho0   , global const double * Rho1   , 
                            global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = ic/(lbmaux[0].Nx*lbmaux[0].Ny);

    for (size_t k=1;k<lbmaux[0].Nneigh;k++)
    {
        size_t inx = (size_t)((long)icx + (long)lbmaux[0].C[k].x + (long)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((long)icy + (long)lbmaux[0].C[k].y + (long)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((long)icz + (long)lbmaux[0].C[k].z + (long)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;
        
        
        double psic = 0.0;
        double psin = 0.0;
        double G    = lbmaux[0].G[0];
        if (!IsSolid0[ic]) psic = lbmaux[0].Psi[0]*exp(-lbmaux[0].Rhoref[0]/Rho0[ic]);
        if (!IsSolid0[in]) psin = lbmaux[0].Psi[0]*exp(-lbmaux[0].Rhoref[0]/Rho0[in]);
        else               G    = lbmaux[0].Gs[0];
        
        BForce0[ic] += -G*lbmaux[0].W[k]*psic*psin*lbmaux[0].C[k];

        psic        = 0.0;
        psin        = 0.0;
        G           = lbmaux[0].G[1];
        if (!IsSolid1[ic]) psic = lbmaux[0].Psi[1]*exp(-lbmaux[0].Rhoref[1]/Rho1[ic]);
        if (!IsSolid1[in]) psin = lbmaux[0].Psi[1]*exp(-lbmaux[0].Rhoref[1]/Rho1[in]);
        else               G    = lbmaux[0].Gs[1];
        
        BForce1[ic] += -G*lbmaux[0].W[k]*psic*psin*lbmaux[0].C[k];

        psic        = 1.0;
        psin        = 1.0;
        G           = lbmaux[0].Gmix;
        if (!IsSolid0[ic]) psic = Rho0[ic];
        if (!IsSolid1[in]) psin = Rho1[in];
        else               G    = lbmaux[0].Gs[0];
        
        BForce0[ic] += -G*lbmaux[0].W[k]*psic*psin*lbmaux[0].C[k];

        psic        = 1.0;
        psin        = 1.0;
        G           = lbmaux[0].Gmix;
        if (!IsSolid1[ic]) psic = Rho1[ic];
        if (!IsSolid0[in]) psin = Rho0[in];
        else               G    = lbmaux[0].Gs[1];
        
        BForce1[ic] += -G*lbmaux[0].W[k]*psic*psin*lbmaux[0].C[k];
    }
    
}

void kernel ApplyForcesMP(global const bool * IsSolid0 , global const bool * IsSolid1 ,
                          global       double3* BForce0, global       double3* BForce1,
                          global const double * Rho0   , global const double * Rho1   , 
                          global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = ic/(lbmaux[0].Nx*lbmaux[0].Ny);

    
    for (size_t k=1;k<lbmaux[0].Nneigh;k++)
    {
        size_t inx = (size_t)((long)icx + (long)lbmaux[0].C[k].x + (long)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((long)icy + (long)lbmaux[0].C[k].y + (long)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((long)icz + (long)lbmaux[0].C[k].z + (long)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;
        
        double psic = 1.0;
        double psin = 1.0;
        double G    = lbmaux[0].Gmix;
        if (!IsSolid0[ic]) psic = Rho0[ic];
        if (!IsSolid1[in]) psin = Rho1[in];
        else               G    = lbmaux[0].Gs[0];
        
        BForce0[ic] += -G*lbmaux[0].W[k]*psic*psin*lbmaux[0].C[k];

        psic        = 1.0;
        psin        = 1.0;
        G           = lbmaux[0].Gmix;
        if (!IsSolid1[ic]) psic = Rho1[ic];
        if (!IsSolid0[in]) psin = Rho0[in];
        else               G    = lbmaux[0].Gs[1];
        
        BForce1[ic] += -G*lbmaux[0].W[k]*psic*psin*lbmaux[0].C[k];
    }
}

void kernel CollideSC    (global const bool * IsSolid, global double * F, global double * Ftemp, global const double3* BForce, global const double3* Vel, global const double * Rho, global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    if (!IsSolid[ic])
    {
        double3 vel  = Vel[ic]+lbmaux[0].Tau[0]*BForce[ic]/Rho[ic];
        //double3 vel  = Vel[ic]+BForce[ic]/Rho[ic];
        //double3 vel  = Vel[ic]+0.5*BForce[ic]/Rho[ic];
        //double3 vel  = Vel[ic];
        double rho   = Rho[ic];
        //double VdotV = dot(vel,vel);
        //double Cs    = lbmaux[0].Cs;
        double tau   = lbmaux[0].Tau[0];
        double NonEq[27];
        double Feq  [27];
        double Q = 0.0;
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            //double VdotC = dot(vel,lbmaux[0].C[k]);
            //Feq  [k]     = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
            Feq[k]       = FeqFluid(k,rho,vel,lbmaux);
            NonEq[k]     = F[ic*lbmaux[0].Nneigh + k] - Feq[k];
            Q           += NonEq[k]*NonEq[k]*lbmaux[0].EEk[k];
        }
        Q = sqrt(2.0*Q);
        tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*lbmaux[0].Sc/rho));
        

        bool valid = true;
        double alpha = 1.0;
        while (valid)
        {
            valid = false;
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                //double Fk = 3.0*(1 - 0.5/lbmaux[0].Tau[0])*dot(BForce[ic],Cs*lbmaux[0].C[k]-vel)*Feq[k]/(rho*Cs*Cs);
                //Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k] - alpha*(NonEq[k]/tau - Fk);
                Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k] - alpha*(NonEq[k]/tau);
                if (Ftemp[ic*lbmaux[0].Nneigh + k]<-1.0e-12)
                {
                    //double temp = F[ic*lbmaux[0].Nneigh + k]/(NonEq[k]/tau - Fk);
                    double temp = tau*F[ic*lbmaux[0].Nneigh + k]/(NonEq[k]);
                    if (temp<alpha) alpha = temp;
                    valid = true;
                }
            }
        }
    }
    else
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
        }
    }
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k]; 
    }
}

void kernel CollideMP    (global const bool * IsSolid0 , global const bool * IsSolid1 ,
                          global double * F0           , global double * F1           , 
                          global double * Ftemp0       , global double * Ftemp1       , 
                          global const double3* BForce0, global const double3* BForce1,
                          global const double3* Vel0   , global const double3* Vel1   ,
                          global const double * Rho0   , global const double * Rho1   , 
                          global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    double3 Vmix = Rho0[ic]*Vel0[ic]/lbmaux[0].Tau[0] + Rho1[ic]*Vel1[ic]/lbmaux[0].Tau[1];
    Vmix = Vmix/(Rho0[ic]/lbmaux[0].Tau[0] + Rho1[ic]/lbmaux[0].Tau[1]);
    if (!IsSolid0[ic])
    {

        double  rho   = Rho0[ic];
        double3 vel   = Vmix + lbmaux[0].Tau[0]*BForce0[ic]/rho;
        double  Cs    = lbmaux[0].Cs; 
        double  VdotV = dot(vel,vel);
        double  tau   = lbmaux[0].Tau[0];
        bool    valid = true;
        double  alpha = 1.0;
        while (valid)
        {
            valid = false;
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                double VdotC = dot(vel,lbmaux[0].C[k]);
                double NonEq = F0[ic*lbmaux[0].Nneigh + k] - lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                Ftemp0[ic*lbmaux[0].Nneigh + k] = F0[ic*lbmaux[0].Nneigh + k] - alpha*(NonEq/tau);
                if (Ftemp0[ic*lbmaux[0].Nneigh + k]<-1.0e-12)
                {
                    double temp = tau*F0[ic*lbmaux[0].Nneigh + k]/(NonEq);
                    if (temp<alpha) alpha = temp;
                    valid = true;
                }
            }
        }
    }
    else
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Ftemp0[ic*lbmaux[0].Nneigh + k] = F0[ic*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
        }
    }
    if (!IsSolid1[ic])
    {
        double  rho  = Rho1[ic];
        double3 vel  = Vmix + lbmaux[0].Tau[1]*BForce1[ic]/rho;
        double  Cs    = lbmaux[0].Cs; 
        double  VdotV = dot(vel,vel);
        double  tau   = lbmaux[0].Tau[1];
        bool    valid = true;
        double  alpha = 1.0;
        while (valid)
        {
            valid = false;
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                double VdotC = dot(vel,lbmaux[0].C[k]);
                double NonEq = F1[ic*lbmaux[0].Nneigh + k] - lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                Ftemp1[ic*lbmaux[0].Nneigh + k] = F1[ic*lbmaux[0].Nneigh + k] - alpha*(NonEq/tau);
                if (Ftemp1[ic*lbmaux[0].Nneigh + k]<-1.0e-12)
                {
                    double temp = tau*F1[ic*lbmaux[0].Nneigh + k]/(NonEq);
                    if (temp<alpha) alpha = temp;
                    valid = true;
                }
            }
        }
    }
    else
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Ftemp1[ic*lbmaux[0].Nneigh + k] = F1[ic*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
        }
    }
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        F0[ic*lbmaux[0].Nneigh + k] = Ftemp0[ic*lbmaux[0].Nneigh + k]; 
        F1[ic*lbmaux[0].Nneigh + k] = Ftemp1[ic*lbmaux[0].Nneigh + k]; 
    }
}

void kernel CollideAD    (global const bool * IsSolid0 , global const bool * IsSolid1 ,
                          global double * F0           , global double * F1           , 
                          global double * Ftemp0       , global double * Ftemp1       , 
                          global const double3* BForce0, global const double3* BForce1,
                          global const double3* Vel0   , global const double3* Vel1   ,
                          global const double * Rho0   , global const double * Rho1   , 
                          global const struct lbm_aux * lbmaux)
{
     
}


void kernel Stream1    (global double * F, global double * Ftemp, global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = ic/(lbmaux[0].Nx*lbmaux[0].Ny);

    for (size_t k=1;k<lbmaux[0].Nneigh;k++)
    {
        size_t inx = (size_t)((long)icx + (long)lbmaux[0].C[k].x + (long)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((long)icy + (long)lbmaux[0].C[k].y + (long)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((long)icz + (long)lbmaux[0].C[k].z + (long)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;
        Ftemp[in*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k];
    }
}

void kernel Stream2     (global const bool * IsSolid, global double * F, global double * Ftemp, global double3* BForce, global double3* Vel, global double * Rho, global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k];
    }
    Rho   [ic] = 0.0;
    Vel   [ic] = (double3)(0.0,0.0,0.0);
    BForce[ic] = (double3)(0.0,0.0,0.0);
    if (!IsSolid[ic])
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Rho[ic] += F[ic*lbmaux[0].Nneigh + k];
            Vel[ic] += F[ic*lbmaux[0].Nneigh + k]*lbmaux[0].C[k];
            //if (ic==0) printf("k: %lu %f %f %f %f \n",k,Rho[ic],Vel[ic].x,Vel[ic].y,Vel[ic].z);
        }
        Vel[ic] *= lbmaux[0].Cs/Rho[ic];
    }
}
