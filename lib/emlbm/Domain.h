/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2014 Sergio Galindo                                    *
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


#ifndef MECHSYS_EMLBM_DOMAIN_H
#define MECHSYS_EMLBM_DOMAIN_H

// STD
#include <map>
#include <vector>
#include <utility>
#include <set>

// MechSys
#include <mechsys/emlbm/Lattice.h>

using std::set;
using std::map;
using std::pair;
using std::make_pair;

namespace EMLBM
{

class Domain
{
public:
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructors
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    Array<double>         nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Special constructor with only one component, the parameters are the same as above
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    double                nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step
    
    //Methods
#ifdef USE_HDF5
    void WriteXDMF         (char const * FileKey);  ///< Write the domain data in xdmf file
#endif

    void Initialize     (double dt=0.0);                                                                                              ///< Set the particles to a initial state and asign the possible insteractions
    void Collide        (size_t n = 0, size_t Np = 1);                                                                                ///< Apply the interaction forces and the collision operator
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                                ///< Solve the Domain dynamics

#ifdef USE_THREAD
    Array<pair<size_t, size_t> >                ListPosPairs;         ///< List of all possible particles pairs
#endif
    //Data
    bool                                         Initialized;         ///< System (particles and interactons) initialized ?
    bool                                              PrtVec;         ///< Print Vector data into the xdmf-h5 files
    bool                                            Finished;         ///< Has the simulation finished
    String                                           FileKey;         ///< File Key for output files
    Array <Lattice>                                      Lat;         ///< Fluid Lattices
    double                                              Time;         ///< Time of the simulation
    double                                                dt;         ///< Timestep
    void *                                          UserData;         ///< User Data
    size_t                                           idx_out;         ///< The discrete time step
    size_t                                              Step;         ///< The space step to reduce the size of the h5 file for visualization
    size_t                                             Nproc;         ///< Number of cores used for the simulation
};

#ifdef USE_THREAD
struct MtData
{
    size_t                  ProcRank; ///< Rank of the thread
    size_t                    N_Proc; ///< Total number of threads
    EMLBM::Domain *              Dom; ///< Pointer to the lbm domain
    double                        dt; ///< Time step
};


void * GlobalCollide (void * Data)
{
    EMLBM::MtData & dat = (*static_cast<EMLBM::MtData *>(Data));
    dat.Dom->Collide(dat.ProcRank, dat.N_Proc);
    return NULL;
}

void * GlobalStream1 (void * Data)
{
    EMLBM::MtData & dat = (*static_cast<EMLBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].Stream1(dat.ProcRank, dat.N_Proc);
    }
    return NULL;
}

void * GlobalStream2 (void * Data)
{
    EMLBM::MtData & dat = (*static_cast<EMLBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].Stream2(dat.ProcRank, dat.N_Proc);
    }
    return NULL;
}

void * GlobalCalcField (void * Data)
{
    EMLBM::MtData & dat = (*static_cast<EMLBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].CalcField(dat.ProcRank, dat.N_Proc);
    }
    return NULL;
}

#endif

inline Domain::Domain(LBMethod Method, Array<double> Tau, iVec3_t Ndim, double Thedx, double Thedt)
{
    Initialized = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (Tau.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one Lattice please");
    for (size_t i=0;i<Tau.Size();i++)
    {
        Lat.Push(Lattice(Method,Tau[i],Ndim,Thedx,Thedt));
    }
    Time   = 0.0;
    dt     = Thedt;
    Step   = 1;
    PrtVec = true;
    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Size()*Lat[0].Ncells,TERM_RST);
}

inline Domain::Domain(LBMethod Method, double Tau, iVec3_t Ndim, double dx, double Thedt)
{
    Initialized = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    Lat.Push(Lattice(Method,Tau,Ndim,dx,Thedt));
    Time   = 0.0;
    dt     = Thedt;
    Step   = 1;
    PrtVec = true;

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Size()*Lat[0].Ncells,TERM_RST);
}

#ifdef USE_HDF5

inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Lat[0].Ndim[0]/Step;
    size_t  Ny = Lat[0].Ndim[1]/Step;
    size_t  Nz = Lat[0].Ndim[2]/Step;
    for (size_t j=0;j<Lat.Size();j++)
    {
        // Creating data sets
        float * Eps       = new float[  Nx*Ny*Nz];
        float * Char      = new float[  Nx*Ny*Nz];
        float * Phi       = new float[  Nx*Ny*Nz];
        float * Cur       = new float[3*Nx*Ny*Nz];
        float * Avec      = new float[3*Nx*Ny*Nz];
        float * Bvec      = new float[3*Nx*Ny*Nz];
        float * Evec      = new float[3*Nx*Ny*Nz];

        size_t i=0;
        for (size_t m=0;m<Lat[0].Ndim(2);m+=Step)
        for (size_t l=0;l<Lat[0].Ndim(1);l+=Step)
        for (size_t n=0;n<Lat[0].Ndim(0);n+=Step)
        {
            double eps    = 0.0;
            double phi    = 0.0;
            double cha    = 0.0;
            Vec3_t avec   = OrthoSys::O;
            Vec3_t cur    = OrthoSys::O;
            Vec3_t bvec   = OrthoSys::O;
            Vec3_t evec   = OrthoSys::O;

            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                eps      += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Eps;
                phi      += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->A[0];
                avec(0)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->A[1];
                avec(1)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->A[2];
                avec(2)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->A[3];
                cha      += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->J[0];
                cur (0)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->J[1];
                cur (1)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->J[2];
                cur (2)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->J[3];
                bvec(0)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->B[0];
                bvec(1)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->B[1];
                bvec(2)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->B[2];
                evec(0)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->E[0];
                evec(1)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->E[1];
                evec(2)  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->E[2];
            }
            eps  /= Step*Step*Step;
            cha  /= Step*Step*Step;
            phi  /= Step*Step*Step;
            cur  /= Step*Step*Step;
            avec /= Step*Step*Step;
            bvec /= Step*Step*Step;
            evec /= Step*Step*Step;
            Eps [i]      = (float) eps;
            Phi [i]      = (float) phi;
            Avec[3*i  ]  = (float) avec(0);
            Avec[3*i+1]  = (float) avec(1);
            Avec[3*i+2]  = (float) avec(2);
            Char[i]      = (float) cha;
            Cur [3*i  ]  = (float) cur (0);
            Cur [3*i+1]  = (float) cur (1);
            Cur [3*i+2]  = (float) cur (2);
            Bvec[3*i  ]  = (float) bvec(0);
            Bvec[3*i+1]  = (float) bvec(1);
            Bvec[3*i+2]  = (float) bvec(2);
            Evec[3*i  ]  = (float) evec(0);
            Evec[3*i+1]  = (float) evec(1);
            Evec[3*i+2]  = (float) evec(2);
            i++;
        }


        //Write the data
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("ScalPot_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Phi );
        dsname.Printf("Charge_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Char);
        dsname.Printf("Epsilon_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Eps );
        if (PrtVec)
        {
            dims[0] = 3*Nx*Ny*Nz;
            dsname.Printf("VecPot_%d",j);
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Avec    );
            dsname.Printf("Current_%d",j);
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Cur     );
            dsname.Printf("MagField_%d",j);
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Bvec    );
            dsname.Printf("ElecField_%d",j);
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Evec    );
        }
        dims[0] = 1;
        int N[1];
        N[0] = Nx;
        dsname.Printf("Nx");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Ny;
        dsname.Printf("Ny");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Nz;
        dsname.Printf("Nz");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);

        delete [] Eps     ;
        delete [] Char    ;
        delete [] Phi     ;
        delete [] Cur     ;
        delete [] Avec    ;
        delete [] Bvec    ;
        delete [] Evec    ;
    }


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"EMLBM_Mesh\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*Lat[0].dx << " " << Step*Lat[0].dx  << " " << Step*Lat[0].dx  << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    for (size_t j=0;j<Lat.Size();j++)
    {
    oss << "     <Attribute Name=\"Epsilon_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Epsilon_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Charge_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Charge_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ScalPot_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ScalPot_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    if (PrtVec)
    {
    oss << "     <Attribute Name=\"Current_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Current_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"VecPot_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/VecPot_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"MagField_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/MagField_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ElecField_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ElecField_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    }
    oss << "   </Grid>\n";
    }
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

#endif

void Domain::Collide (size_t n, size_t Np)
{
	size_t Ni = Lat[0].Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Lat[0].Ncells : Fn = (n+1)*Ni;
#ifdef USE_OMP
    In = 0;
    Fn = Lat[0].Ncells;
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            for (size_t k=0;k<c->Nneigh;k++)
            {
                c->Htemp[k] = c->H[k] - (c->H[k] - c->Heq(k))/Lat[j].Tau;
                for (size_t mu=0;mu<4;mu++)
                {
                    c->Ftemp[mu][k] = c->F[mu][k] - (c->F[mu][k] - c->Feq(mu,k))/Lat[j].Tau;
                    c->Gtemp[mu][k] = c->G[mu][k] - (c->G[mu][k] - c->Geq(mu,k))/Lat[j].Tau;
                    //if (c->Index[0]==Lat[0].Ndim[0]/2&&c->Index[1]==Lat[0].Ndim[1]/2&&c->Index[2]==Lat[0].Ndim[2]/2) std::cout << c->Feq(mu,k) << " " << c->Geq(mu,k) << " " << c->A[mu] << " " << c->Sig[mu] << " " << mu << " " << k << std::endl;
                    //if (c->Index[0]==15&&c->Index[1]==15&&c->Index[2]==15) std::cout << c->Feq(mu,k) << " " << c->Geq(mu,k) << " " << c->A[mu] << " " << c->Sig[mu] << " " << c->G[mu][k] << " " << mu << " " << k << std::endl;
                    //if (c->Index[0]==15&&c->Index[1]==15&&c->Index[2]==15) std::cout << c->Ftemp[mu][k] << " " << c->F[mu][k] << " " << mu << " " << k << std::endl;
                }
            }
            for (size_t k=0;k<c->Nneigh;k++)
            {
                c->H[k] = c->Htemp[k];
                for (size_t mu=0;mu<4;mu++)
                {
                    c->F[mu][k] = c->Ftemp[mu][k];
                    c->G[mu][k] = c->Gtemp[mu][k];
                    //if (c->Index[0]==Lat[0].Ndim[0]/2&&c->Index[1]==Lat[0].Ndim[1]/2&&c->Index[2]==Lat[0].Ndim[2]/2) std::cout << c->Feq(mu,k) << " " << c->Geq(mu,k) << " " << c->A[mu] << " " << c->Sig[mu] << " " << mu << " " << k << std::endl;
                    //if (c->Index[0]==15&&c->Index[1]==15&&c->Index[2]==15) std::cout << c->Feq(mu,k) << " " << c->Geq(mu,k) << " " << c->A[mu] << " " << c->Sig[mu] << " " << c->G[mu][k] << " " << mu << " " << k << std::endl;
                    //if (c->Index[0]==15&&c->Index[1]==15&&c->Index[2]==15) std::cout << c->Ftemp[mu][k] << " " << c->F[mu][k] << " " << mu << " " << k << std::endl;
                }
            }
        }
    }   
}

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{

    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Finished = false;

    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                   , TERM_RST);
    for (size_t i=0;i<Lat.Size();i++)
    {
    printf("%s  Tau of Lattice %zd                 =  %g%s\n"       ,TERM_CLR2, i, Lat[i].Tau                        , TERM_RST);
    }

    Nproc = TheNproc;

    for (size_t j=0;j<Lat.Size();j++)
    {
        for (size_t i=0;i<Lat[j].Ncells;i++)
        {
            Lat[j].Cells[i]->Initialize();
            Lat[j].Cells[i]->CalcProp();
        }
    }


#ifdef USE_THREAD
    EMLBM::MtData MTD[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = this;
        MTD[i].dt       = Lat[0].dt;
    }
    pthread_t thrs[Nproc];   
    
#else

#endif
    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time >= tout)
        {
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    WriteXDMF(fn.CStr());
                    #else
                    //WriteVTK (fn.CStr());
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }


#ifdef USE_THREAD
        //GlobalCollide
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalCollide, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //GlobalStream1
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalStream1, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //GlobalStream2
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalStream2, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //GlobalCalcField
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalCalcField, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
#elif USE_OMP 
        Collide(1,Nproc);
        for (size_t i=0;i<Lat.Size();i++)
        {
            Lat[i].Stream1  (1,Nproc);
            Lat[i].Stream2  (1,Nproc);
            Lat[i].CalcField(1,Nproc);
        }

#endif

        Time += dt;
    }
    // last output
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}
}


#endif

