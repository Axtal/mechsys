/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
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


#ifndef MECHSYS_LBMDEM_DOMAIN_H
#define MECHSYS_LBMDEM_DOMAIN_H

// Mechsys
#include <mechsys/flbm/Domain.h>
#include <mechsys/dem/domain.h>
#ifdef USE_CUDA
#include <mechsys/lbmdem/lbmdem.cuh>
#endif

#include <chrono>

namespace LBMDEM
{

struct ParticleCellPair
{
    iVec3_t ICell;        ///< Index of the cell
    size_t IPar;          ///< Index of the particle
    Array<size_t> IGeo;   ///< Array of index of the geometric feature
};

struct MtData;

class Domain
{
public:

    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    // Constructors
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    double                    nu, ///< Viscosity of the fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    Domain (LBMethod  Method,     ///< Type of array, for example D2Q9
    double                nu,     ///< Viscosity of the fluid
    char const *     DEMfile,     ///< A DEM save file 
    double                dx,     ///< Space spacing
    double                dt,     ///< Time step
    Vec3_t  xmin=OrthoSys::O,     ///< minimun point of the domain to be rendered
    Vec3_t  xmax=OrthoSys::O      ///< maximun point of the domain to be rendered
    );    
    //Methods
    void Reset();                    ///< Reset LBM grid
    void ImprintLattice();           ///< Imprint the DEM particles into the LBM grid
    void ResetParCell();             ///< Reset the information of particle cell contacts
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);            ///< Solve the Domain dynamics
    
#ifdef USE_HDF5    
    void WriteXDMF         (char const * FileKey);                                      ///< Save a xdmf file for visualization
#endif

    //Methods for CUDA 
    #ifdef USE_CUDA
    void UpLoadDevice (size_t Nc, bool first);                                ///< Upload the buffers into the coprocessor device
    void DnLoadDevice (size_t Nc);                                            ///< Download the buffers from the coprocessor device
    #endif

    //Variables
    double           dt;                ///< Time Step
    double           dx;                ///< Grid size
    double        Alpha;                ///< Verlet distance
    size_t        Nproc;                ///< Number of processors
    size_t      idx_out;                ///< The discrete time step for output
    size_t         iter;                ///< Iteration counter
    String      FileKey;                ///< File Key for output files
    void *     UserData;                ///< User Data
    double         Time;                ///< Simulation time variable
    bool       Finished;                ///< Boolen flag to signal the end of simulation
    bool      PeriodicX;                ///< Flag to signal periodic boundary conditions in X direction
    bool      PeriodicY;                ///< Flag to signal periodic boundary conditions in Y direction
    bool      PeriodicZ;                ///< Flag to signal periodic boundary conditions in Z direction
    FLBM::Domain LBMDOM;                ///< The LBM domain
    DEM ::Domain DEMDOM;                ///< The DEM domain
    MtData *        MTD;                ///< Multithread data
    Array <ParticleCellPair> ParCellPairs; ///< Pairs of cells and particles

#ifdef USE_CUDA
    size_t                                 Nthread=256;       ///< Number of CUDA threads
    thrust::device_vector<real>            bOmeis;            ///< Buffer with the DEM collision information
    thrust::device_vector<real>            bGamma;            ///< Buffer with the volume fraction information
    thrust::device_vector<real>            bGammaf;           ///< Buffer with the volume fraction information
    thrust::device_vector<ParCellPairCU>   bPaCe;             ///< Buffer with the information of Particle cell pairs
    thrust::device_vector<size_t>          bPaCeF;            ///< Buffer with the geometric features assigned to a Particle cell pair
    thrust::device_vector<size_t>          bPaCeV;            ///< Buffer with the Particle cell pair information for spheres
    lbmdem_aux                             lbmdemaux;         ///< Auxiliary data for lbmdem simulations                                                          
    real                                 * pOmeis;            ///< Pointer to Buffer with the DEM collision information
    real                                 * pGamma;            ///< Pointer to Buffer with the volume fraction information
    real                                 * pGammaf;           ///< Pointer to Buffer with the prescribed volume fraction information 
    ParCellPairCU                        * pPaCe;             ///< Pointer to particle cell pair
    size_t                               * pPaCeF;            ///< Pointer to particle cell pair
    size_t                               * pPaCeV;            ///< Pointer to particle cell pair
    lbmdem_aux                           * plbmdemaux;        ///< pointer to auxiliary data
                                                        
#endif
};

struct MtData
{
    size_t                        ProcRank; ///< Rank of the thread
    size_t                          N_Proc; ///< Total number of threads
    LBMDEM::Domain *                   Dom; ///< Pointer to the lbm domain
    double                             Dmx; ///< Maximun displacement
    double                              dt; ///< Time step
    //Array <Vec3_t>                    FLBM; ///< Array of contributions to the force  for each core
    //Array <Vec3_t>                    TLBM; ///< Array of contributions to the torque for each core
    Array<ParticleCellPair>            LPC; ///< A temporal array of possible particle cell contacts
    Array<ParticleCellPair>           LPCP; ///< Same array for periodic boundary information
};

inline Domain::Domain(LBMethod TheMethod, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Nproc  = 1;
    idx_out= 0;
    iter   = 0;
    dt = Thedt;
    dx = Thedx;
    Time = 0.0;
    Alpha = 0.05;
    PeriodicX= false;
    PeriodicY= false;
    PeriodicZ= false;
    LBMDOM = FLBM::Domain(TheMethod, Thenu, TheNdim, Thedx, Thedt);
    DEMDOM = DEM ::Domain();

    LBMDOM.Omeis = new double *** [TheNdim(0)];
    LBMDOM.Gamma = new double **  [TheNdim(0)];
    LBMDOM.Gammaf= new double **  [TheNdim(0)];
    for (size_t ix=0; ix< TheNdim(0); ix++)
    {
        LBMDOM.Omeis [ix] = new double ** [TheNdim(1)];
        LBMDOM.Gamma [ix] = new double *  [TheNdim(1)];
        LBMDOM.Gammaf[ix] = new double *  [TheNdim(1)];
        for (size_t iy=0; iy< TheNdim(1); iy++)
        {
            LBMDOM.Omeis [ix][iy] = new double * [TheNdim(2)];
            LBMDOM.Gamma [ix][iy] = new double   [TheNdim(2)];
            LBMDOM.Gammaf[ix][iy] = new double   [TheNdim(2)];
            for (size_t iz=0; iz< TheNdim(2); iz++)
            {
                //LBMDOM.Gammaf[ix][iy][iz] = 0.0;
                LBMDOM.Omeis [ix][iy][iz] = new double [LBMDOM.Nneigh];
            }
        }
    }
}

inline Domain::Domain(LBMethod TheMethod, double Thenu, char const * DEMfile, double Thedx, double Thedt, Vec3_t xmin, Vec3_t xmax)
{
    Nproc  = 1;
    idx_out= 0;
    iter   = 0;
    dt = Thedt;
    dx = Thedx;
    Time = 0.0;
    Alpha = 0.05;
    PeriodicX= false;
    PeriodicY= false;
    PeriodicZ= false;
    DEMDOM = DEM ::Domain();
    DEMDOM.Load(DEMfile);
    Vec3_t Xmax,Xmin;
    DEMDOM.BoundingBoxAll(Xmin,Xmax);
    //iVec3_t shift(1,1,1);
    if (fabs(DEMDOM.Xmax-DEMDOM.Xmin)>1.0e-12) 
    {
        PeriodicX = true;
        //shift(0)   = 0;
    }
    else if(fabs(xmax(0)-xmin(0))>1.0e-12)
    {
        Xmin(0) = xmin(0);
        Xmax(0) = xmax(0);
    }

    if (fabs(DEMDOM.Ymax-DEMDOM.Ymin)>1.0e-12) 
    {
        PeriodicY = true;
        //shift(1)   = 0;
    }
    else if (fabs(xmax(1)-xmin(1))>1.0e-12)
    {
        Xmin(1) = xmin(1);
        Xmax(1) = xmax(1);
    }

    if (fabs(DEMDOM.Zmax-DEMDOM.Zmin)>1.0e-12) 
    {
        PeriodicZ = true;
        //shift(2)   = 0;
    }
    else if (fabs(xmax(2)-xmin(2))>1.0e-12)
    {
        Xmin(2) = xmin(2);
        Xmax(2) = xmax(2);
    }
    
    Vec3_t Transport(-Xmin);
    for (size_t i=0; i<DEMDOM.Particles.Size(); i++) DEMDOM.Particles[i]->Translate(Transport);

    //iVec3_t TheNdim = (Xmax-Xmin)/dx+shift;
    iVec3_t TheNdim;
    TheNdim(0) = round((Xmax(0)-Xmin(0))/dx);
    TheNdim(1) = round((Xmax(1)-Xmin(1))/dx);
    TheNdim(2) = round((Xmax(2)-Xmin(2))/dx);


    LBMDOM = FLBM::Domain(TheMethod, Thenu, TheNdim, Thedx, Thedt);

    LBMDOM.Omeis = new double *** [TheNdim(0)];
    LBMDOM.Gamma = new double **  [TheNdim(0)];
    LBMDOM.Gammaf= new double **  [TheNdim(0)];
    for (size_t ix=0; ix< TheNdim(0); ix++)
    {
        LBMDOM.Omeis [ix] = new double ** [TheNdim(1)];
        LBMDOM.Gamma [ix] = new double *  [TheNdim(1)];
        LBMDOM.Gammaf[ix] = new double *  [TheNdim(1)];
        for (size_t iy=0; iy< TheNdim(1); iy++)
        {
            LBMDOM.Omeis [ix][iy] = new double * [TheNdim(2)];
            LBMDOM.Gamma [ix][iy] = new double   [TheNdim(2)];
            LBMDOM.Gammaf[ix][iy] = new double   [TheNdim(2)];
            for (size_t iz=0; iz< TheNdim(2); iz++)
            {
                //LBMDOM.Gammaf[ix][iy][iz] = 0.0;
                LBMDOM.Omeis [ix][iy][iz] = new double [LBMDOM.Nneigh];
            }
        }
    }

}

void Domain::ImprintLattice()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i = 0;i<ParCellPairs.Size();i++)
    {
        DEM::Particle  * Pa   = DEMDOM.Particles[ParCellPairs[i].IPar];
        Vec3_t C = dx*ParCellPairs[i].ICell;
        size_t n = ParCellPairs[i].ICell(0);
        size_t m = ParCellPairs[i].ICell(1);
        size_t l = ParCellPairs[i].ICell(2);
        Vec3_t  Xtemp,Xs,Xstemp,S,St;
        double len = 12.0*dx,minl = Pa->Dmax;
        Vec3_t  Pert = DEMDOM.Per;
        if (!Pa->IsFree()) Pert = OrthoSys::O;
        if (DEM::Distance(C,Pa->x,Pert)>Pa->Dmax) continue;
        Vec3_t Nor = OrthoSys::O;
        if (ParCellPairs[i].IGeo.Size()>0) 
        {
            if (Pa->Faces.Size()>0)
            {
                DEM::Distance(C,*Pa->Faces[ParCellPairs[i].IGeo[0]],Xtemp,Xs,S,Pert);
                minl = norm(S);
                Nor = Pa->Faces[ParCellPairs[i].IGeo[0]]->Nor;
                for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                {
                    DEM::Distance(C,*Pa->Faces[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp,St,Pert);
                    if (norm(St) < minl)
                    {
                        minl = norm(St);
                        Xs   = Xstemp;
                        Nor = Pa->Faces[ParCellPairs[i].IGeo[j]]->Nor;
                        S   = St;
                    }
                }
            }
            else if (Pa->Edges.Size()>0)
            {
                DEM::Distance(C,*Pa->Edges[ParCellPairs[i].IGeo[0]],Xtemp,Xs,S,Pert);
                minl = norm(S);
                for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                {
                    DEM::Distance(C,*Pa->Edges[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp,St,Pert);
                    if (norm(St) < minl)
                    {
                        minl = norm(St);
                        Xs   = Xstemp;
                        S    = St;
                    }
                }
            }
            else if (Pa->Verts.Size()>0)
            {
                DEM::Distance(C,*Pa->Verts[ParCellPairs[i].IGeo[0]],Xtemp,Xs,S,Pert);
                minl = norm(S);
                for (size_t j=1;j<ParCellPairs[i].IGeo.Size();j++)
                {
                    DEM::Distance(C,*Pa->Verts[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp,St,Pert);
                    if (norm(St) < minl)
                    {
                        minl = norm(Xtemp-Xstemp);
                        Xs   = Xstemp;
                        S    = St;
                    }
                }
            }
            double dotpro = -dot(S,Nor);
            if (dotpro>0.0||fabs(dotpro)<0.95*minl||Pa->Faces.Size()<4||!Pa->Closed) 
            {
                Vec3_t Xst = S + C;
                len = DEM::SphereCube(Xst,C,Pa->Props.R,dx);
            }
        }
        if (fabs(len)<1.0e-12) continue;
        double Tau = LBMDOM.Tau[0];
        double gamma  = len/(12.0*dx);
        if (gamma<LBMDOM.Gamma[n][m][l]) continue;
        LBMDOM.Gamma[n][m][l] = gamma;
        //Vec3_t B      = C - Pa->x;
        Vec3_t B;
        DEM::BranchVec(Pa->x,C,B,Pert);
        Vec3_t tmp;
        Rotation(Pa->w,Pa->Q,tmp);
        Vec3_t VelP   = Pa->v + cross(tmp,B);
        double rho = LBMDOM.Rho[0][n][m][l];
        double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
        //double Bn  = gamma;
        size_t ncells = LBMDOM.Nneigh;
        Vec3_t Flbm = OrthoSys::O;
        for (size_t k=0;k<ncells;k++)
        {
            double Fvpp     = LBMDOM.Feq(LBMDOM.Op[k],rho,VelP);
            double Fvp      = LBMDOM.Feq(k           ,rho,VelP);
            double Omega    = LBMDOM.F[0][n][m][l][LBMDOM.Op[k]] - Fvpp - (LBMDOM.F[0][n][m][l][k] - Fvp);
            //cell->Omeis[k] += Omega;
            //cell->Omeis[k] += gamma*Omega;
            LBMDOM.Omeis[n][m][l][k] = Omega;
            Flbm += -Bn*Omega*LBMDOM.Cs*LBMDOM.Cs*dx*dx*LBMDOM.C[k];
        }
        Vec3_t Tlbm,Tt;
        Tt =           cross(B,Flbm);
        Quaternion_t q;
        Conjugate    (Pa->Q,q);
        Rotation     (Tt,q,Tlbm);
        //std::cout << "1" << std::endl;
        omp_set_lock      (&Pa->lck);
        Pa->F          += Flbm;
        Pa->Flbm       += Flbm;
        Pa->T          += Tlbm;
        omp_unset_lock    (&Pa->lck);
        //MTD[omp_get_thread_num()].FLBM[Pa->Index] += Flbm;
        //MTD[omp_get_thread_num()].TLBM[Pa->Index] += Tlbm;
    }
}

void Domain::ResetParCell()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].LPC.Resize(0);
        MTD[i].LPCP.Resize(0);
    }

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip<DEMDOM.Particles.Size(); ip++)
    {
        DEM::Particle * Pa = DEMDOM.Particles[ip];
        if (Pa->Bdry) continue;
        bool free = Pa->IsFree();
        int nbmin = floor(free&&PeriodicX ? (Pa->MinX()-2.0*Alpha-dx)/dx : std::max(0.0,double(Pa->MinX()-2.0*Alpha-dx)/dx));
        int mbmin = floor(free&&PeriodicY ? (Pa->MinY()-2.0*Alpha-dx)/dx : std::max(0.0,double(Pa->MinY()-2.0*Alpha-dx)/dx));
        int lbmin = floor(free&&PeriodicZ ? (Pa->MinZ()-2.0*Alpha-dx)/dx : std::max(0.0,double(Pa->MinZ()-2.0*Alpha-dx)/dx));
        int nbmax = ceil (free&&PeriodicX ? (Pa->MaxX()+2.0*Alpha+dx)/dx : std::min(double(LBMDOM.Ndim(0)-1),double(Pa->MaxX()+2.0*Alpha+dx)/dx));
        int mbmax = ceil (free&&PeriodicY ? (Pa->MaxY()+2.0*Alpha+dx)/dx : std::min(double(LBMDOM.Ndim(1)-1),double(Pa->MaxY()+2.0*Alpha+dx)/dx));
        int lbmax = ceil (free&&PeriodicZ ? (Pa->MaxZ()+2.0*Alpha+dx)/dx : std::min(double(LBMDOM.Ndim(2)-1),double(Pa->MaxZ()+2.0*Alpha+dx)/dx));
        for (int nb = nbmin;nb<= nbmax;nb++)
        for (int mb = mbmin;mb<= mbmax;mb++)
        for (int lb = lbmin;lb<= lbmax;lb++)
        {
            size_t n     = (nb+LBMDOM.Ndim(0))%LBMDOM.Ndim(0);
            size_t m     = (mb+LBMDOM.Ndim(1))%LBMDOM.Ndim(1);
            size_t l     = (lb+LBMDOM.Ndim(2))%LBMDOM.Ndim(2);
            double x     = dx*n;
            double y     = dx*m;
            double z     = dx*l;
            Vec3_t  C(x,y,z);
            Vec3_t  Pert = DEMDOM.Per;
            if (!free) Pert = OrthoSys::O;
            if ((DEM::Distance(C,Pa->x,Pert)>2.0*Alpha+2.0*dx+Pa->Dmax)||LBMDOM.IsSolid[0][n][m][l]) continue;
            ParticleCellPair NewPCP;
            NewPCP.IPar = ip;
            NewPCP.ICell= iVec3_t(n,m,l);
            bool valid = true;
            Vec3_t Nor = OrthoSys::O;
            if (Pa->Faces.Size()>0)
            {
                double minl = 2.0*Pa->Dmax;
                Vec3_t B,Xs = OrthoSys::O;
                for (size_t j=0;j<Pa->Faces.Size();j++)
                {
                    Vec3_t Xstemp,Xtemp,S;
                    DEM::Distance(C,*Pa->Faces[j],Xtemp,Xstemp,S,Pert);
                    double dist = norm(S);
                    if (dist<minl)
                    {
                        Nor  = Pa->Faces[j]->Nor;
                        minl = dist;
                        Xs   = Xstemp;
                        B    = -S;
                    }
                    if (dist<2.0*Alpha+Pa->Props.R)
                    {
                        if (Pa->Faces[j]->Area()<2.0*M_PI*Pa->Props.R*Pa->Props.R)
                        {
                            continue;
                        }
                        NewPCP.IGeo.Push(j);
                    }
                }
                double dotpro = dot(B,Nor);
                if ((dotpro>0.0||fabs(dotpro)<0.95*minl)&&Pa->Faces.Size()>3&&minl>2.0*Alpha+Pa->Props.R)       valid = false;
                else if (minl>2.0*Alpha+Pa->Props.R&&(Pa->Faces.Size()<4||!Pa->Closed))                         valid = false;
                else if (Pa->Faces.Size()>3&&Pa->Closed&&NewPCP.IGeo.Size()==0&&!Pa->IsInsideFaceOnly(C,Pert))  valid = false;
            }
            else if (Pa->Edges.Size()>0)
            {
                for (size_t j=0;j<Pa->Edges.Size();j++)
                {
                    if (DEM::Distance(C,*Pa->Edges[j],Pert)<2.0*Alpha+Pa->Props.R) 
                    {
                        NewPCP.IGeo.Push(j);
                    }
                    else valid = false;
                }
            }
            else if (Pa->Verts.Size()>0)
            {
                for (size_t j=0;j<Pa->Verts.Size();j++)
                {
                    if (DEM::Distance(C,*Pa->Verts[j],Pert)<2.0*Alpha+Pa->Props.R)
                    {
                        NewPCP.IGeo.Push(j);
                    }
                    else valid = false;
                }
            }
            if (valid) MTD[omp_get_thread_num()].LPC.Push(NewPCP);
        }
    }
    ParCellPairs.Resize(0);
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LPC.Size();j++)
        {
            ParCellPairs.Push(MTD[i].LPC[j]);
        }
    }
}

void Domain::Reset()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ix=0; ix < LBMDOM.Ndim(0); ix++)
    for (size_t iy=0; iy < LBMDOM.Ndim(1); iy++)
    for (size_t iz=0; iz < LBMDOM.Ndim(2); iz++)
    {
        LBMDOM.Gamma[ix][iy][iz] = (double) LBMDOM.IsSolid[0][ix][iy][iz];
        for (size_t k=0; k<LBMDOM.Nneigh; k++)
        {
            LBMDOM.Omeis[ix][iy][iz][k] = 0.0;
        }
    }

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < DEMDOM.Particles.Size(); ip++)
    {
        DEMDOM.Particles[ip]->F    = DEMDOM.Particles[ip]->Ff;
        DEMDOM.Particles[ip]->Flbm = OrthoSys::O;
        DEMDOM.Particles[ip]->T    = DEMDOM.Particles[ip]->Tf;
        //for (size_t i=0;i<Nproc;i++)
        //{
            //MTD[i].FLBM[ip] = OrthoSys::O;
            //MTD[i].TLBM[ip] = OrthoSys::O;
        //}
    }

}

void Domain::WriteXDMF(char const * FileKey)
{
}

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{
    FileKey.Printf("%s",TheFileKey);
    idx_out = 0;
    Finished = false;
    //Assigning the value for the time step to the domain variable
    DEMDOM.Dt = dt;
    Nproc        = TheNproc;
    DEMDOM.Nproc = TheNproc;
    LBMDOM.Nproc = TheNproc;
    DEMDOM.Initialize (dt);

    MTD = new LBMDEM::MtData[Nproc];
    DEMDOM.MTD = new DEM::MtData[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = this;
        MTD[i].Dmx      = 0.0;
        //MTD[i].FLBM.Resize(DEMDOM.Particles.Size());
        //MTD[i].TLBM.Resize(DEMDOM.Particles.Size());
        DEMDOM.MTD[i].N_Proc   = Nproc;
        DEMDOM.MTD[i].ProcRank = i;
        DEMDOM.MTD[i].Dom      = &DEMDOM;
        DEMDOM.MTD[i].Dmx      = 0.0;
    }

    DEMDOM.Vs = 0.0;
    DEMDOM.Ms = 0.0;
    DEMDOM.MaxDmax =  0.0;
    DEMDOM.Alpha   = Alpha;
    double MaxKn   =  0.0;
    double MaxBn   =  0.0;
    double MinDmax = -1.0;
    double MinMass = -1.0;
    for (size_t i=0; i<DEMDOM.Particles.Size(); i++) 
    { 
        if (DEMDOM.Particles[i]->IsFree())
        {
            DEMDOM.Vs += DEMDOM.Particles[i]->Props.V;
            DEMDOM.Ms += DEMDOM.Particles[i]->Props.m;
            if (DEMDOM.Particles[i]->Dmax     > DEMDOM.MaxDmax) DEMDOM.MaxDmax = DEMDOM.Particles[i]->Dmax;
            if (DEMDOM.Particles[i]->Props.Kn > MaxKn         ) MaxKn   = DEMDOM.Particles[i]->Props.Kn;
            if (DEMDOM.Particles[i]->Dmax     < MinDmax||(MinDmax<0.0)) MinDmax = DEMDOM.Particles[i]->Dmax;
            if (DEMDOM.Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = DEMDOM.Particles[i]->Props.m;
            DEMDOM.FreePar.Push(i);
        }
        else DEMDOM.NoFreePar.Push(i);
        DEMDOM.Particles[i]->Index = i;
    }
    for (size_t i=0; i<DEMDOM.BInteractons.Size(); i++)
    {
        double pbn = std::max(DEMDOM.BInteractons[i]->Bn/DEMDOM.BInteractons[i]->L0,DEMDOM.BInteractons[i]->Bt/DEMDOM.BInteractons[i]->L0);
        if (pbn > MaxBn) MaxBn = pbn;
    }

    if (PeriodicX) 
    {
        DEMDOM.Xmin = 0.0;
        DEMDOM.Xmax = dx*LBMDOM.Ndim(0);
    }
    if (PeriodicY) 
    {
        DEMDOM.Ymin = 0.0;
        DEMDOM.Ymax = dx*LBMDOM.Ndim(1);
    }
    if (PeriodicZ) 
    {
        DEMDOM.Zmin = 0.0;
        DEMDOM.Zmax = dx*LBMDOM.Ndim(2);
    }

    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1                                  , TERM_RST);
    printf("%s  Number of cores                           =  %zd%s\n"      ,TERM_CLR2, Nproc                                 , TERM_RST);
    printf("%s  Time step                                 =  %g%s\n"       ,TERM_CLR2, dt                                                         , TERM_RST);
    printf("%s  Grid size                                 =  %g%s\n"       ,TERM_CLR2, dx                                                         , TERM_RST);
    printf("%s  C parameter                               =  %g%s\n"       ,TERM_CLR2, dx/dt                                                      , TERM_RST);
    printf("%s  Number of LBM cells                       =  %zd%s\n"      ,TERM_CLR2, LBMDOM.Ncells                                              , TERM_RST);
    for (size_t i=0;i<LBMDOM.Nl;i++)
    {
    printf("%s  Tau of Lattice %zd                          =  %g%s\n"        ,TERM_CLR2, i, LBMDOM.Tau[i]                                          , TERM_RST);
    }
    printf("%s  Total mass   of free particles            =  %g%s\n"        ,TERM_CLR4, DEMDOM.Ms                                                 , TERM_RST);
    printf("%s  Total volume of free particles            =  %g%s\n"        ,TERM_CLR4, DEMDOM.Vs                                                 , TERM_RST);
    printf("%s  Total number of particles                 =  %zd%s\n"       ,TERM_CLR2, DEMDOM.Particles.Size()                                   , TERM_RST);
    printf("%s  Time step                                 =  %g%s\n"        ,TERM_CLR2, dt                                                        , TERM_RST);
    printf("%s  Simulated Time                            =  %g%s\n"        ,TERM_CLR2, Tf                                                        , TERM_RST);
    printf("%s  Verlet distance                           =  %g%s\n"        ,TERM_CLR2, Alpha                                                     , TERM_RST);
    printf("%s  Suggested Time Step                       =  %g%s\n"        ,TERM_CLR5, 0.1*sqrt(MinMass/(MaxKn+MaxBn))                           , TERM_RST);
    printf("%s  Suggested Verlet distance                 =  %g or  %g%s\n" ,TERM_CLR5, 0.5*MinDmax, 0.25*(MinDmax + DEMDOM.MaxDmax)              , TERM_RST);
    if (PeriodicX)
    printf("%s  Periodic Boundary conditions in X between =  %g and %g%s\n" ,TERM_CLR5, DEMDOM.Xmin, DEMDOM.Xmax                                  , TERM_RST);
    if (PeriodicY)
    printf("%s  Periodic Boundary conditions in Y between =  %g and %g%s\n" ,TERM_CLR5, DEMDOM.Ymin, DEMDOM.Ymax                                  , TERM_RST);
    if (PeriodicZ)
    printf("%s  Periodic Boundary conditions in Z between =  %g and %g%s\n" ,TERM_CLR5, DEMDOM.Zmin, DEMDOM.Zmax                                  , TERM_RST);
    
    if (Alpha>DEMDOM.Beta*DEMDOM.MaxDmax&&DEMDOM.MaxDmax>1.0e-12)
    {
        Alpha = DEMDOM.Beta*DEMDOM.MaxDmax;
        DEMDOM.Alpha = Alpha;
        printf("%s  Verlet distance changed to                =  %g%s\n"   ,TERM_CLR2, Alpha                                    , TERM_RST);
    }
    fflush(stdout); 

    DEMDOM.Per = Vec3_t(DEMDOM.Xmax-DEMDOM.Xmin,DEMDOM.Ymax-DEMDOM.Ymin,DEMDOM.Zmax-DEMDOM.Zmin);

    Reset();
    ResetParCell();
    ImprintLattice();
    DEMDOM.UpdateContacts();
#ifdef USE_CUDA
    DEMDOM.UpLoadDevice(Nproc,true);
    LBMDOM.UpLoadDevice(Nproc);
    UpLoadDevice(Nproc,true);
    cudaDeviceProp prop;
    cudaGetDeviceProperties (&prop,0);
    std::cout 
        << TERM_CLR2 
        << "  Using GPU:                                =  " << prop.name              << TERM_RST << std::endl
        << "  Using Number of CUDA threads:             =  " << Nthread                << TERM_RST << std::endl
        << "  Number of vertices                        =  " << DEMDOM.demaux.nverts   << TERM_RST << std::endl
        << "  Number of edges                           =  " << DEMDOM.demaux.nedges/2 << TERM_RST << std::endl
        << "  Number of faces                           =  " << DEMDOM.demaux.nfaces/2 << TERM_RST << std::endl;
#endif

    size_t iter_b = iter;
    size_t iter_t = 0;
    size_t numup  = 0;

    double tout = Time;
    while (Time<Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time>=tout)
        {
#ifdef USE_CUDA
            DEMDOM.DnLoadDevice(Nproc,true);
            LBMDOM.DnLoadDevice(Nproc);
            DnLoadDevice(Nproc);
#endif
            if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            if (TheFileKey!=NULL)
            {
                String fndem,fbdem;
                fndem.Printf    ("%s_dem_%04d", TheFileKey, idx_out);
                fbdem.Printf    ("%s_dem_bf_%04d", TheFileKey, idx_out);
                String fnlbm;
                fnlbm.Printf    ("%s_lbm_%04d", TheFileKey, idx_out);
                DEMDOM.WriteXDMF    (fndem.CStr());
                DEMDOM.WriteBF      (fbdem.CStr());
                LBMDOM.WriteXDMF_DEM(fnlbm.CStr());
            }
            idx_out++;
            tout += dtOut;
        }
#ifdef USE_CUDA
        //auto start = std::chrono::high_resolution_clock::now();

        cudaReset<<<LBMDOM.Nl*LBMDOM.Ncells/Nthread+1,Nthread>>>(LBMDOM.pIsSolid,pGammaf,pGamma,pOmeis,LBMDOM.plbmaux);

        DEM::Reset<<<(DEMDOM.demaux.nparts+DEMDOM.demaux.ncoint)/Nthread+1,Nthread>>>(DEMDOM.pParticlesCU,DEMDOM.pDynParticlesCU,DEMDOM.pComInteractons,DEMDOM.pdemaux);
        DEM::CalcForceVV<<<DEMDOM.demaux.nvvint/Nthread+1,Nthread>>>(DEMDOM.pInteractons,DEMDOM.pComInteractons, DEMDOM.pDynInteractonsVV, DEMDOM.pParticlesCU, DEMDOM.pDynParticlesCU, DEMDOM.pdemaux);
        DEM::CalcForceEE<<<DEMDOM.demaux.neeint/Nthread+1,Nthread>>>(DEMDOM.pEdgesCU,DEMDOM.pVertsCU,DEMDOM.pInteractons, DEMDOM.pComInteractons, DEMDOM.pDynInteractonsEE, DEMDOM.pParticlesCU, DEMDOM.pDynParticlesCU, DEMDOM.pdemaux);
        DEM::CalcForceVF<<<DEMDOM.demaux.nvfint/Nthread+1,Nthread>>>(DEMDOM.pFacesCU,DEMDOM.pFacidCU,DEMDOM.pVertsCU,DEMDOM.pInteractons, DEMDOM.pComInteractons, DEMDOM.pDynInteractonsVF, DEMDOM.pParticlesCU, DEMDOM.pDynParticlesCU, DEMDOM.pdemaux);
        DEM::CalcForceFV<<<DEMDOM.demaux.nfvint/Nthread+1,Nthread>>>(DEMDOM.pFacesCU,DEMDOM.pFacidCU,DEMDOM.pVertsCU,DEMDOM.pInteractons, DEMDOM.pComInteractons, DEMDOM.pDynInteractonsFV, DEMDOM.pParticlesCU, DEMDOM.pDynParticlesCU, DEMDOM.pdemaux);
        //cudaDeviceSynchronize();
        
        //auto stop  = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        //std::cout << "DEM Calculation of forces " << duration.count() << std::endl;

        //start = std::chrono::high_resolution_clock::now();

        cudaImprintLatticeVC<<<lbmdemaux.nvc/Nthread+1,Nthread>>>(pPaCeV,DEMDOM.pParticlesCU,DEMDOM.pDynParticlesCU,LBMDOM.pRho,pGamma,pOmeis,LBMDOM.pF,DEMDOM.pdemaux,LBMDOM.plbmaux,plbmdemaux);
        cudaImprintLatticeFC<<<lbmdemaux.nfc/Nthread+1,Nthread>>>(pPaCe,pPaCeF,DEMDOM.pFacesCU,DEMDOM.pFacidCU,DEMDOM.pVertsCU,DEMDOM.pParticlesCU,DEMDOM.pDynParticlesCU,LBMDOM.pRho,pGamma,pOmeis,LBMDOM.pF,DEMDOM.pdemaux,LBMDOM.plbmaux,plbmdemaux);

        //cudaDeviceSynchronize();
        //stop  = std::chrono::high_resolution_clock::now();

        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        //std::cout << "LBM-DEM Imprint           " << duration.count() << std::endl;

        //start = std::chrono::high_resolution_clock::now();

        DEM::Translate<<<DEMDOM.demaux.nparts/Nthread+1,Nthread>>>(DEMDOM.pVertsCU, DEMDOM.pParticlesCU, DEMDOM.pDynParticlesCU, DEMDOM.pdemaux);
        DEM::Rotate   <<<DEMDOM.demaux.nparts/Nthread+1,Nthread>>>(DEMDOM.pVertsCU, DEMDOM.pParticlesCU, DEMDOM.pDynParticlesCU, DEMDOM.pdemaux);
        DEM::MaxD     <<<DEMDOM.demaux.nverts/Nthread+1,Nthread>>>(DEMDOM.pVertsCU, DEMDOM.pVertsoCU, DEMDOM.pMaxDCU, DEMDOM.pdemaux);

        real maxdis = 0.0;
        thrust::device_vector<real>::iterator it=thrust::max_element(DEMDOM.bMaxDCU.begin(),DEMDOM.bMaxDCU.end());
        maxdis = *it;

        DEMDOM.demaux.iter++;
        DEMDOM.demaux.Time += dt;

        if (maxdis>Alpha)
        {
            numup++;
            iter_t+= iter - iter_b;
            iter_b = iter;
            DEMDOM.UpdateContactsDevice();
            ResetParCell();
            UpLoadDevice(Nproc,false);
        }

        //cudaDeviceSynchronize();
        //stop  = std::chrono::high_resolution_clock::now();
//
        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
//
        //std::cout << "DEM Final                " << duration.count() << std::endl;

        //start = std::chrono::high_resolution_clock::now();

        FLBM::cudaCollideSCDEM<<<LBMDOM.Nl*LBMDOM.Ncells/Nthread+1,Nthread>>>(LBMDOM.pIsSolid,LBMDOM.pF,LBMDOM.pFtemp,LBMDOM.pBForce,LBMDOM.pVel,LBMDOM.pRho,pGamma,pOmeis,LBMDOM.plbmaux);
        real * tmp = LBMDOM.pF;
        LBMDOM.pF = LBMDOM.pFtemp;
        LBMDOM.pFtemp = tmp;
        FLBM::cudaStream1<<<LBMDOM.Nl*LBMDOM.Ncells/Nthread+1,Nthread>>>(LBMDOM.pF,LBMDOM.pFtemp,LBMDOM.plbmaux);
        tmp = LBMDOM.pF;
        LBMDOM.pF = LBMDOM.pFtemp;
        LBMDOM.pFtemp = tmp;
        //FLBM::cudaStream2<<<LBMDOM.Nl*LBMDOM.Ncells/Nthread+1,Nthread>>>(LBMDOM.pIsSolid,LBMDOM.pF,LBMDOM.pFtemp,LBMDOM.pBForce,LBMDOM.pVel,LBMDOM.pRho,LBMDOM.plbmaux);
        cudaStream2<<<LBMDOM.Nl*LBMDOM.Ncells/Nthread+1,Nthread>>>(LBMDOM.pIsSolid,pGamma,pOmeis,LBMDOM.pF,LBMDOM.pFtemp,LBMDOM.pBForce,LBMDOM.pVel,LBMDOM.pRho,LBMDOM.plbmaux);

        //cudaDeviceSynchronize();
        //stop  = std::chrono::high_resolution_clock::now();

        //duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        //std::cout << "LBM Final                " << duration.count() << std::endl;
#else
        //std::cout << "1" <<std::endl;
        //Initialize domain
        Reset();


        //std::cout << "2" <<std::endl;
        //Calculate DEM forces
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<DEMDOM.Interactons.Size(); i++)
        {
		    if (DEMDOM.Interactons[i]->CalcForce(dt,DEMDOM.Per,iter))
            {
                String f_error(FileKey+"_error");
                DEMDOM.Save     (f_error.CStr());
                DEMDOM.WriteXDMF(f_error.CStr());
                std::cout << "Maximun overlap detected between particles at time " << Time << std::endl;
                std::cout << "Iteration number                                   " << iter << std::endl;
                sleep(1);
                throw new Fatal("Maximun overlap detected between particles");
            }
            omp_set_lock  (&DEMDOM.Interactons[i]->P1->lck);
            DEMDOM.Interactons[i]->P1->F += DEMDOM.Interactons[i]->F1;
            DEMDOM.Interactons[i]->P1->T += DEMDOM.Interactons[i]->T1;
            omp_unset_lock(&DEMDOM.Interactons[i]->P1->lck);
            omp_set_lock  (&DEMDOM.Interactons[i]->P2->lck);
            DEMDOM.Interactons[i]->P2->F += DEMDOM.Interactons[i]->F2;
            DEMDOM.Interactons[i]->P2->T += DEMDOM.Interactons[i]->T2;
            omp_unset_lock(&DEMDOM.Interactons[i]->P2->lck);
        }
        
        //std::cout << "3" <<std::endl;
        //Imprint the DEM into the LBM lattice to calcualte forces
        ImprintLattice();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Nproc;i++)
        {
            MTD[i].Dmx = 0.0;
        }

        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<DEMDOM.Particles.Size(); i++)
        {
            //std::cout << "1" << std::endl;
		    DEMDOM.Particles[i]->Translate(dt);
            //std::cout << "2" << std::endl;
		    DEMDOM.Particles[i]->Rotate(dt);
            //std::cout << "3" << std::endl;
            if (DEMDOM.Particles[i]->MaxDisplacement()>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = DEMDOM.Particles[i]->MaxDisplacement();
        }

        double maxdis = 0.0;
        for (size_t i=0;i<Nproc;i++)
        {
            if (maxdis<MTD[i].Dmx) maxdis = MTD[i].Dmx;
        }

        if (maxdis>Alpha)
        {
            DEMDOM.UpdateContacts();
            ResetParCell();
        }

        LBMDOM.CollideSCDEM();
        LBMDOM.StreamSC();
#endif
        
        Time += dt; 
        DEMDOM.iter++;
        iter++;
    }
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    printf("%s  Total number of iterations                                = %zd%s\n",TERM_CLR4,iter, TERM_RST);
    if (numup>0)
    {
    printf("%s  Average number of iterations between contact list updates = %zd%s\n",TERM_CLR4,iter_t/numup, TERM_RST);
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

#ifdef USE_CUDA
inline void Domain::UpLoadDevice(size_t Nc, bool first)
{
    if (first)
    {
        
        thrust::host_vector<real>  hOmeis         (LBMDOM.Ncells*LBMDOM.Nneigh);
        thrust::host_vector<real>  hGamma         (LBMDOM.Ncells); 
        thrust::host_vector<real>  hGammaf        (LBMDOM.Ncells); 
    
        for (size_t nz=0;nz<LBMDOM.Ndim(2);nz++)
        {
            #pragma omp parallel for schedule(static) num_threads(Nc)
            for (size_t ny=0;ny<LBMDOM.Ndim(1);ny++)
            {
                for (size_t nx=0;nx<LBMDOM.Ndim(0);nx++)
                {
                    size_t Nm = nx + ny*LBMDOM.Ndim(0) + nz*LBMDOM.Ndim(1)*LBMDOM.Ndim(0);
                    hGamma [Nm]        = LBMDOM.Gamma [nx][ny][nz];
                    hGammaf[Nm]        = LBMDOM.Gammaf[nx][ny][nz];
                    for (size_t nn=0;nn<LBMDOM.Nneigh;nn++)
                    {
                        size_t Nn = nn + nx*LBMDOM.Nneigh + ny*LBMDOM.Ndim(0)*LBMDOM.Nneigh + nz*LBMDOM.Ndim(1)*LBMDOM.Ndim(0)*LBMDOM.Nneigh; 
                        hOmeis[Nn] = LBMDOM.Omeis[nx][ny][nz][nn];
                    }
                }
            }
        }

        bOmeis = hOmeis;
        bGamma = hGamma;
        bGammaf= hGammaf;

        pOmeis = thrust::raw_pointer_cast(bOmeis .data());
        pGamma = thrust::raw_pointer_cast(bGamma .data());
        pGammaf= thrust::raw_pointer_cast(bGammaf.data());
    }
    
    size_t   idvc = 0;
    size_t   idfc = 0;
    size_t   icfc = 0;
    size_t * Ivc = new size_t[2*ParCellPairs.Size()];
    size_t * Ifc = new size_t[2*ParCellPairs.Size()];
    size_t * Icc = new size_t[2*ParCellPairs.Size()];
    for(size_t npc=0;npc<ParCellPairs.Size();npc++)
    {
        Ivc[2*npc  ] = idvc; 
        Ifc[2*npc  ] = idfc; 
        Icc[2*npc  ] = icfc; 
        if(DEMDOM.Particles[ParCellPairs[npc].IPar]->Verts.Size()==1)
        {
            idvc++;
        }
        else
        {
            idfc ++;
            icfc += ParCellPairs[npc].IGeo.Size();
        }
        Ivc[2*npc+1] = idvc; 
        Ifc[2*npc+1] = idfc; 
        Icc[2*npc+1] = icfc; 
    }
    lbmdemaux.nvc = idvc;
    lbmdemaux.nfc = idfc;
    thrust::host_vector<ParCellPairCU>   hPaCe (  lbmdemaux.nfc);
    thrust::host_vector<size_t>          hPaCeF(           icfc);
    thrust::host_vector<size_t>          hPaCeV(2*lbmdemaux.nvc);

    #pragma omp parallel for schedule(static) num_threads(Nc)
    for(size_t npc=0;npc<ParCellPairs.Size();npc++)
    {
        if(DEMDOM.Particles[ParCellPairs[npc].IPar]->Verts.Size()==1)
        {
            hPaCeV[2*Ivc[2*npc]  ] = FLBM::Pt2idx(ParCellPairs[npc].ICell,LBMDOM.Ndim);
            hPaCeV[2*Ivc[2*npc]+1] = ParCellPairs[npc].IPar;
        }
        else
        {
            ParCellPairCU PCU;
            PCU.Ic = FLBM::Pt2idx(ParCellPairs[npc].ICell,LBMDOM.Ndim);
            PCU.Ip = ParCellPairs[npc].IPar;
            PCU.Nfi = Icc[2*npc  ];
            PCU.Nff = Icc[2*npc+1];
            for (size_t i=Icc[2*npc];i<Icc[2*npc+1];i++) 
            {
                hPaCeF[i] = ParCellPairs[npc].IGeo[i-Icc[2*npc]]+DEMDOM.Particles[ParCellPairs[npc].IPar]->Nfi;
            }
            hPaCe[Ifc[2*npc]] = PCU;
        }
    }
    delete [] Ivc;
    delete [] Ifc;
    delete [] Icc;

    bPaCe .resize(hPaCe .size());
    bPaCeF.resize(hPaCeF.size());
    bPaCeV.resize(hPaCeV.size());

    thrust::copy(hPaCe .begin(),hPaCe .end(),bPaCe .begin());
    thrust::copy(hPaCeF.begin(),hPaCeF.end(),bPaCeF.begin());
    thrust::copy(hPaCeV.begin(),hPaCeV.end(),bPaCeV.begin());

    pPaCe  = thrust::raw_pointer_cast(bPaCe .data());
    pPaCeF = thrust::raw_pointer_cast(bPaCeF.data());
    pPaCeV = thrust::raw_pointer_cast(bPaCeV.data());

    if (!first) cudaFree  (plbmdemaux);
    cudaMalloc(&plbmdemaux, sizeof(lbmdem_aux));
    cudaMemcpy(plbmdemaux, &lbmdemaux, sizeof(lbmdem_aux), cudaMemcpyHostToDevice);
}

inline void Domain::DnLoadDevice(size_t Nc)
{
    thrust::host_vector<real>  hGamma  (LBMDOM.Ncells);

    hGamma = bGamma;
    
    for (size_t nz=0;nz<LBMDOM.Ndim(2);nz++)
    {
        #pragma omp parallel for schedule(static) num_threads(Nc)
        for (size_t ny=0;ny<LBMDOM.Ndim(1);ny++)
        {
            for (size_t nx=0;nx<LBMDOM.Ndim(0);nx++)
            {
                size_t Nm = nx + ny*LBMDOM.Ndim(0) + nz*LBMDOM.Ndim(1)*LBMDOM.Ndim(0);
                LBMDOM.Gamma[nx][ny][nz]    = hGamma[Nm];
            }
        }
    }
}
#endif
}
#endif //MECHSYS_LBMDEM_DOMAIN_H
