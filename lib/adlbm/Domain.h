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


#ifndef MECHSYS_ADLBM_DOMAIN_H
#define MECHSYS_ADLBM_DOMAIN_H

// STD
#include <map>
#include <vector>
#include <utility>
#include <set>

// MechSys
#include <mechsys/adlbm/Lattice.h>

using std::set;
using std::map;
using std::pair;
using std::make_pair;

namespace ADLBM
{

class Domain
{
public:
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructors
    Domain (
    LBMethod              Method,
    double                Thenu,
    double                Thedif,
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Methods
#ifdef USE_HDF5
    void WriteXDMF         (char const * FileKey);  ///< Write the domain data in xdmf file
#endif

    void Initialize     (double dt=0.0);                                                                                              ///< Set the particles to a initial state and asign the possible insteractions
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                                ///< Solve the Domain dynamics

    //Data
    bool                                         Initialized;         ///< System (particles and interactons) initialized ?
    bool                                              PrtVec;         ///< Print Vector data into the xdmf-h5 files
    bool                                            Finished;         ///< Has the simulation finished
    String                                           FileKey;         ///< File Key for output files
    String                                          Concname;         ///< String to identify the solute
    String                                          Fluxname;         ///< String to identify the solute flux
    Lattice                                              Lat;         ///< Fluid Lattices
    double                                              Time;         ///< Time of the simulation
    double                                                dt;         ///< Timestep
    void *                                          UserData;         ///< User Data
    size_t                                           idx_out;         ///< The discrete time step
    size_t                                              Step;         ///< The space step to reduce the size of the h5 file for visualization
    size_t                                             Nproc;         ///< Number of cores used for the simulation
};

inline Domain::Domain(LBMethod Method, double Thenu, double Thedif, iVec3_t Ndim, double Thedx, double Thedt)
{
    Initialized = false;
    Util::Stopwatch stopwatch;
    if (Ndim(2) >1&&(Method==D2Q9||Method==D2Q5))  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (Ndim(2)==1&&(Method==D3Q15||Method==D3Q19)) throw new Fatal("LBM::Domain: Ndim(2) is 1. Either change the method to D2Q9 or increase the z-dimension");
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    Lat = Lattice(Method, Thenu, Thedif, Ndim,Thedx,Thedt);
    Time   = 0.0;
    dt     = Thedt;
    Step   = 1;
    PrtVec = true;
    Concname= "Concentration";
    Fluxname= "Massflux";
    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Ncells,TERM_RST);

}


#ifdef USE_HDF5

inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Lat.Ndim[0]/Step;
    size_t  Ny = Lat.Ndim[1]/Step;
    size_t  Nz = Lat.Ndim[2]/Step;
    // Creating data sets
    float * Rho       = new float[  Nx*Ny*Nz];
    float * Temp      = new float[  Nx*Ny*Nz];
    float * Gamma     = new float[  Nx*Ny*Nz];
    float * Dif       = new float[  Nx*Ny*Nz];
    float * Vel       = new float[3*Nx*Ny*Nz];
    float * Flux      = new float[3*Nx*Ny*Nz];

    size_t i=0;
    for (size_t m=0;m<Lat.Ndim(2);m+=Step)
    for (size_t l=0;l<Lat.Ndim(1);l+=Step)
    for (size_t n=0;n<Lat.Ndim(0);n+=Step)
    {
        double rho    = 0.0;
        double temp   = 0.0;
        double gamma  = 0.0;
        double dif    = 0.0;
        Vec3_t vel    = OrthoSys::O;
        Vec3_t flux   = OrthoSys::O;

        for (size_t ni=0;ni<Step;ni++)
        for (size_t li=0;li<Step;li++)
        for (size_t mi=0;mi<Step;mi++)
        {
            rho      += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Rho;
            temp     += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Temp;
            gamma    += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->IsSolid;
            dif      += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Dif;
            vel (0)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel [0];
            vel (1)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel [1];
            vel (2)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel [2];
            flux(0)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Flux[0];
            flux(1)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Flux[1];
            flux(2)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Flux[2];
        }
        rho   /= Step*Step*Step;
        temp  /= Step*Step*Step;
        gamma /= Step*Step*Step;
        dif   /= Step*Step*Step;
        vel   /= Step*Step*Step;
        flux  /= Step*Step*Step;
        Rho    [i]      = (float) rho;
        Temp   [i]      = (float) temp;
        Gamma  [i]      = (float) gamma;
        Dif    [i]      = (float) dif;
        Vel    [3*i  ]  = (float) vel (0);
        Vel    [3*i+1]  = (float) vel (1);
        Vel    [3*i+2]  = (float) vel (2);
        Flux   [3*i  ]  = (float) flux(0);
        Flux   [3*i+1]  = (float) flux(1);
        Flux   [3*i+2]  = (float) flux(2);
        i++;
    } 
    //Write the data
    hsize_t dims[1];
    dims[0] = Nx*Ny*Nz;
    String dsname;
    dsname.Printf("Density");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Rho );
    dsname.Printf(Concname.CStr());
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Temp);
    dsname.Printf("Gamma");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Gamma);
    dsname.Printf("Diffusion");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Dif);
    if (PrtVec)
    {
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vel );
        dsname.Printf(Fluxname.CStr());
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Flux);
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

    delete [] Rho     ;
    delete [] Temp    ;
    delete [] Gamma   ;
    delete [] Dif     ;
    delete [] Vel     ;
    delete [] Flux    ;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"ADLBM_Mesh\" GridType=\"Uniform\">\n";
    if (Lat.Ndim[2]==1)
    {
    oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Ny << " " << Nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> " << Step*Lat.dx  << " " << Step*Lat.dx  << "\n";
    size_t t = Nz; Nz = Nx; Nx=t;
    }
    else
    {
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*Lat.dx << " " << Step*Lat.dx  << " " << Step*Lat.dx  << "\n";
    }
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Density" << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Density" << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"" << Concname.CStr() << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/"<< Concname.CStr() << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Gamma" << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Gamma" << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Diffusion" << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Diffusion" << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    if (PrtVec)
    {
    oss << "     <Attribute Name=\"Velocity" << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity" << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"" << Fluxname.CStr() << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/" << Fluxname.CStr() << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    }
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

#endif

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{

    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Finished = false;

    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"     ,TERM_CLR2, dt                                     , TERM_RST);
    printf("%s  Tau of Lattice                   =  %g%s\n"     ,TERM_CLR2, Lat.Tau                                , TERM_RST);
    printf("%s  C   of Lattice                   =  %g%s\n"     ,TERM_CLR2, Lat.dx/Lat.dt                          , TERM_RST);

    Nproc = TheNproc;

    for (size_t i=0;i<Lat.Ncells;i++)
    {
        Lat.Cells[i]->Tauc = 3.0*Lat.Cells[i]->Dif*Lat.dt/(Lat.dx*Lat.dx) + 0.5;
    }



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
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }


#ifdef USE_OMP 
        Lat.Collide  (Nproc);
        Lat.Stream1  (Nproc);
        Lat.Stream2  (Nproc);
        Lat.CalcProps(Nproc);
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

