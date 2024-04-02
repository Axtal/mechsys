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

#ifndef MECHSYS_LBM_LATTICE_H
#define MECHSYS_LBM_LATTICE_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// MechSys
#include <mechsys/lbm/Cell.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>

class Lattice
{
public:
    //typedefs
    typedef void (*ptFun_t) (Lattice & Lat, void * UserData);

    //Constructors
    Lattice () {};            //Default
    Lattice (LBMethod Method, double nu, iVec3_t Ndim, double dx, double dt);

    //Methods
    void Solve(double Tf, double dtOut, ptFun_t ptSetup=NULL, ptFun_t ptReport=NULL,
               char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);   ///< Solve the LBM equation in time
    void Stream(size_t n = 0, size_t Np = 1);                                       ///< Stream the velocity distributions
    void Stream1(size_t n = 0, size_t Np = 1);                                      ///< Stream the velocity distributions
    void Stream2(size_t n = 0, size_t Np = 1);                                      ///< Stream the velocity distributions
    void SolidDisk(Vec3_t const & X, double R);                                     ///< Add a solid fixed disk in space
    void SetZeroGamma(size_t n=0, size_t Np = 1);                                   ///< Set an initial value for the fluid/solid ratio of each cell
    double Psi(double);                                                             ///< Interaction potential
    double SolidFraction();                                                         ///< Solid fraction
    void ApplyForce();                                                              ///< Apply molecular forces
    void Collide();                                                                 ///< apply the collision operator
    void CollideAlt();                                                              ///< apply the collision operator
    void BounceBack(size_t n = 0, size_t Np = 1);                                   ///< apply interaction with solids
    void WriteVTK(char const * FileKey);                                            ///< Write the state in a VTK file
#ifdef USE_HDF5
    void WriteXDMF(char const * FileKey);                                           ///< Write the state in a XDMF file
#endif
    Cell * GetCell(iVec3_t const & v);                                              ///< Get pointer to cell at v


     

    //Data
    size_t                                    idx_out;          // The discrete time step
    size_t                                    Ncells;           // Number of cells per lattice
    double                                    Time;             // The current time
    double                                    G;                // Interaction strength
    double                                    Gs;               // Interaction strength with solids
    double                                    Nu;               // Real viscosity
    iVec3_t                                   Ndim;             // Integer dimension of the domain
    double                                    dx;               // grid space
    double                                    dt;               // time step
    double                                    Tau;              // Relaxation time
    double                                    Rhoref;           // Values for th intermolecular force
    double                                    Psiref;           // 
    //Array<Cell *>                             Cells;            // Array of pointer cells
    Cell                                   ** Cells;            // Array of pointer cells
    //Array<std::pair<Cell *, Cell*> >          CellPairs;        // Array of pairs of cells for interaction
    void *                                    UserData;         // User Data
};

inline Lattice::Lattice(LBMethod TheMethod, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Nu   = Thenu;
    Ndim = TheNdim;
    dx   = Thedx;
    dt   = Thedt;
    Tau  = 3.0*Nu*dt/(dx*dx) + 0.5;
    Rhoref = 200.0;
    Psiref = 4.0;
    G      = 0.0;
    Gs     = 0.0;

    //Cells.Resize(Ndim[0]*Ndim[1]*Ndim[2]);
    Cells = new Cell * [Ndim[0]*Ndim[1]*Ndim[2]];
    Ncells = Ndim[0]*Ndim[1]*Ndim[2];
    size_t n = 0;
    for (size_t k=0;k<Ndim[2];k++)
    for (size_t j=0;j<Ndim[1];j++)
    for (size_t i=0;i<Ndim[0];i++)
    {
        //Cells[n] =  new Cell(n,TheMethod,iVec3_t(i,j,k),Ndim,dx/dt,Tau);
        Cells[n] = new Cell(n,TheMethod,iVec3_t(i,j,k),Ndim,dx/dt,Tau);
        n++;
    } 
    //for (size_t i=0;i<Ndim[0]*Ndim[1]*Ndim[2];i++)
    //{
        //Cell * c = Cells[i];
        //for (size_t j=1;j<c->Nneigh;j++)
        //{
            //Cell * nb     = Cells[c->Neighs[j]];
            //if (nb->ID>c->ID)
            //{
                //std::pair<Cell *, Cell*> p;
                //p = std::make_pair(c,nb);
                //CellPairs.Push(p);
            //}
        //}
    //}
}

inline void Lattice::Stream(size_t n, size_t Np)
{
	size_t Ni = Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Ncells : Fn = (n+1)*Ni;
    // Assign temporal distributions
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule(static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    for (size_t j=1;j<Cells[i]->Nneigh;j++)
    {
        Cells[Cells[i]->Neighs[j]]->Ftemp[j] = Cells[i]->F[j];
    }

    //Swap the distribution values
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule(static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        double * Ftemp   = Cells[i]->F;
        Cells[i]->F      = Cells[i]->Ftemp;
        Cells[i]->Ftemp  = Ftemp;
        //for (size_t j=1;j<Cells[i]->Nneigh;j++)
        //{
            //Cells[i]->F[j] = Cells[i]->Ftemp[j];
        //}
        Cells[i]->Rho = Cells[i]->VelDen(Cells[i]->Vel);
    }
}

inline void Lattice::Stream1(size_t n, size_t Np)
{
	size_t Ni = Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Ncells : Fn = (n+1)*Ni;
    // Assign temporal distributions
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule(static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    for (size_t j=1;j<Cells[i]->Nneigh;j++)
    {
        //std::cout << i << " " << j << " " << Cells[i]->Neighs[j] << " " <<  Cells[i]->IsSolid << std::endl;
        Cells[Cells[i]->Neighs[j]]->Ftemp[j] = Cells[i]->F[j];
    }
}

inline void Lattice::Stream2(size_t n, size_t Np)
{
	size_t Ni = Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Ncells : Fn = (n+1)*Ni;
    //Swap the distribution values
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule(static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        double * Ftemp   = Cells[i]->F;
        Cells[i]->F      = Cells[i]->Ftemp;
        Cells[i]->Ftemp  = Ftemp;
        //for (size_t j=1;j<Cells[i]->Nneigh;j++)
        //{
            //Cells[i]->F[j] = Cells[i]->Ftemp[j];
        //}
        Cells[i]->Rho = Cells[i]->VelDen(Cells[i]->Vel);
    }
}

inline void Lattice::SolidDisk(Vec3_t const & X, double R)
{
    for (size_t n=std::max(0.0,double(X(0)-R-dx)/dx);n<=std::min(double(Ndim(0)-1),double(X(0)+R+dx)/dx);n++)
    for (size_t m=std::max(0.0,double(X(1)-R-dx)/dx);m<=std::min(double(Ndim(1)-1),double(X(1)+R+dx)/dx);m++)
    for (size_t l=std::max(0.0,double(X(2)-R-dx)/dx);l<=std::min(double(Ndim(2)-1),double(X(2)+R+dx)/dx);l++)
    {
        Cell  * cell = GetCell(iVec3_t(n,m,l));
        double x     = dx*cell->Index(0);
        double y     = dx*cell->Index(1);
        double z     = dx*cell->Index(2);
        Vec3_t  C(x,y,z);
        if (norm(C-X)<R) cell->IsSolid = true;
    }
}

inline void Lattice::SetZeroGamma(size_t n, size_t Np)
{
	size_t Ni = Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Ncells : Fn = (n+1)*Ni;
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule(static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        Cells[i]->Gammap = Cells[i]->Gamma;
        Cells[i]->Gamma  = 0.0;
        Cells[i]->BForce = Cells[i]->BForcef;
        size_t ncells = Cells[i]->Nneigh;
        for (size_t j=0;j<ncells;j++)
        {
            Cells[i]->Omeis[j] = 0.0;
        }
    }
}

inline double Lattice::Psi(double rho)
{
    return Psiref*exp(-Rhoref/rho);
    //double a = 0.0015;
    //double b = 1.0/3000;
    //double RT= 1.0/3.0;
    //double p=RT*(rho/(1-b*rho))-a*rho*rho/(1+rho*b);
    //if (RT*rho-p>0.0) return  sqrt( 2.0*(p - RT*rho)/(RT*G));
    //else              return -sqrt(-2.0*(p - RT*rho)/(RT*G));
}

inline double Lattice::SolidFraction()
{
    double Sf = 0.0;
    for (size_t i=0; i<Ncells; i++)
    {
        if (Cells[i]->IsSolid||Cells[i]->Gamma>0.0) Sf+=1.0;
    }
    return Sf/(Ncells);
}

inline void Lattice::ApplyForce()
{
    //for (size_t i=0;i<CellPairs.Size();i++)
    //{
        //Cell * c  = CellPairs[i].first;
        //Cell * nb = CellPairs[i].second;
        //if (fabs(c->Gamma-1.0)+fabs(nb->Gamma-1.0)<1.0e-12) continue;
        //double psi    = Psi(c->Rho/);
        //double nb_psi = Psi(nb->Rho/);
        //double C      = G;
        //if (nb->IsSolid)
        //{
            //C       = Gs;
            //nb_psi  = 1.0;
        //}
        //if (c->IsSolid)
        //{
//
        //}
    //}

    for (size_t i=0;i<Ncells;i++)
    {
        Cell * c = Cells[i];
        double psi = Psi(c->Rho);
        if (fabs(c->Gamma-1.0)<1.0e-12) continue;
        for (size_t j=1;j<c->Nneigh;j++)
        {
            Cell * nb     = Cells[c->Neighs[j]];
            double nb_psi = Psi(nb->Rho);
            double C      = G;
            if (nb->Gamma>0.0||nb->IsSolid)
            {
                nb_psi = 1.0;
                C      = Gs;
            }
            c->BForce    += -C*psi*c->W[j]*nb_psi*c->C[j];
        }
    }
}

inline void Lattice::Collide()
{
    double ome = 1.0/Tau;
    for (size_t i=0;i<Ncells    ;i++)
    {
        Cell * c = Cells[i];
        if (c->IsSolid) continue;
        //if (fabs(c->Gamma-1.0)<1.0e-12) continue;
        Vec3_t V;
        double rho = c->VelDen(V);
        double Bn  = (c->Gamma*(Tau-0.5))/((1.0-c->Gamma)+(Tau-0.5));
        for (size_t j=0;j<c->Nneigh;j++)
        {
            double Feqn = c->Feq(j,       V,rho);
            //double Fvp  = c->Feq(j,c->VelP,rho);
            c->F[j] = c->F[j] - (1 - Bn)*ome*(c->F[j] - Feqn) + Bn*c->Omeis[j];
        }
    }
}

inline void Lattice::CollideAlt()
{
    double ome = 1.0/Tau;
    for (size_t i=0;i<Ncells    ;i++)
    {
        Cell * c = Cells[i];
        if (c->IsSolid) continue;
        if (fabs(c->Gamma-1.0)<1.0e-12&&fabs(G)>1.0e-12) continue;
        Vec3_t V;
        double rho = c->VelDen(V);
        Vec3_t DV  = V + c->BForce*Tau/rho;
        //Vec3_t DV  = V + c->BForce*dt/rho;
        double Bn  = (c->Gamma*(Tau-0.5))/((1.0-c->Gamma)+(Tau-0.5));
        bool valid  = true;
        //double omel = ome;
        //double omet = ome;
        double alphal = 1.0;
        double alphat = 1.0;
        size_t num  = 0;
        while (valid)
        {
            valid = false;
            alphal  = alphat;
            for (size_t j=0;j<c->Nneigh;j++)
            {
                //double Feqn  = c->Feq(j,       V,rho);
                double FDeqn = c->Feq(j,      DV,rho);
                //double Fvp  = c->Feq(j, c->VelP,rho);

                //First method owen
                //c->F[j] = c->F[j] - (1 - Bn)*ome*(c->F[j] - Feqn) + Bn*c->Omeis[j] + dt*dt*(1 - Bn)*c->W[j]*dot(c->BForce,c->C[j])/(K*dx);

                //Second method LBM EOS
                //c->Ftemp[j] = c->F[j] - alphal*((1 - Bn)*(ome*(c->F[j] - Feqn) - FDeqn + Feqn) - Bn*c->Omeis[j]);
                
                //Third method sukop
                //c->Ftemp[j] = c->F[j] - (1 - Bn)*omel*(c->F[j] - FDeqn) + Bn*c->Omeis[j];
                c->Ftemp[j] = c->F[j] - alphal*((1 - Bn)*ome*(c->F[j] - FDeqn) - Bn*c->Omeis[j]);

                if (c->Ftemp[j]<0.0&&num<1)
                {
                    double temp = fabs(c->F[j]/((1 - Bn)*ome*(c->F[j] - FDeqn) - Bn*c->Omeis[j]));
                    if (temp<alphat) alphat = temp;
                    valid = true;
                    //std::cout << i << " " << c->F[j] << " " << FDeqn << " " << c->Omeis[j] << " " << Bn << " " << alphat << " " << temp << std::endl;
                }
            }
            num++;
            if (num>2) 
            {
                throw new Fatal("Lattice::Collide: Redefine your time step, the current value ensures unstability");
            }
        }
        for (size_t j=0;j<c->Nneigh;j++)
        {
            c->F[j] = fabs(c->Ftemp[j]);
        }
    }
}

inline void Lattice::BounceBack(size_t n, size_t Np)
{
	size_t Ni = Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Ncells : Fn = (n+1)*Ni;
#ifdef USE_OMP
    In = 0;
    Fn = Ncells;
    #pragma omp parallel for schedule(static) num_threads(Np)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        if (!Cells[i]->IsSolid) continue;
        for (size_t j = 0;j<Cells[i]->Nneigh;j++) Cells[i]->Ftemp[j] = Cells[i]->F[j];
        for (size_t j = 0;j<Cells[i]->Nneigh;j++) Cells[i]->F[j]     = Cells[i]->Ftemp[Cells[i]->Op[j]];
    }
}

#ifdef USE_HDF5
inline void Lattice::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Creating data sets
    float *Density   = new float[Ncells];
    float *Gamma     = new float[Ncells];
    float *Velocity  = new float[Ncells];
    float *MassFlux  = new float[Ncells];


    for (size_t i=0;i<Ncells;i++)
    {
        double rho;
        Vec3_t vel;
        rho = Cells[i]->VelDen(vel);
        Density  [i] = (float) rho;
        Gamma    [i] = (float) Cells[i]->IsSolid? 1.0:Cells[i]->Gamma;
        Velocity [i] = (float) norm(vel);
        MassFlux [i] = (float) rho*norm(vel);
    }

    //Write the data
    hsize_t dims[1];
    dims[0] = Ncells;
    H5LTmake_dataset_float(file_id,"Density" ,1,dims,Density );
    H5LTmake_dataset_float(file_id,"Gamma"   ,1,dims,Gamma   );
    H5LTmake_dataset_float(file_id,"Velocity",1,dims,Velocity);
    H5LTmake_dataset_float(file_id,"MassFlux",1,dims,MassFlux);

    delete [] Density ;
    delete [] Gamma   ;
    delete [] Velocity;
    delete [] MassFlux;

    //Closing the file
    H5Fclose(file_id);

	// Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Ndim(1) << " " << Ndim(0) << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Density\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Gamma\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"MassFlux\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/MassFlux\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
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

inline void Lattice::WriteVTK(char const * FileKey)
{
	// Header
	std::ostringstream oss;
	oss << "# vtk DataFile Version 2.0\n";
	oss << "TimeStep = " << idx_out << "\n";
	oss << "ASCII\n";
	oss << "DATASET STRUCTURED_POINTS\n";
	oss << "DIMENSIONS " << Ndim[0] << " " << Ndim[1] << " " << Ndim[2] << "\n";
	oss << "ORIGIN "     << 0       << " " << 0       << " " << 0       << "\n";
	oss << "SPACING "    << 1       << " " << 1       << " " << 1       << "\n";
	oss << "POINT_DATA " << Ndim[0]*Ndim[1]*Ndim[2]   << "\n";

	// Solid cells
	oss << "SCALARS Geom float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<Ncells; ++i)
	{
		if (Cells[i]->IsSolid) oss << "1.0\n";
		else                   oss << "0.0\n";
	}

	// Density field
	oss << "SCALARS Density float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<Ncells; i++)
		oss << Cells[i]->Rho << "\n";

	// Density field
	oss << "SCALARS Gamma float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<Ncells; i++)
		oss << Cells[i]->Gamma << "\n";

	oss << "VECTORS Velocity float\n";
	for (size_t i=0; i<Ncells; ++i)
	{
		Vec3_t v;  Cells[i]->Velocity(v);
		oss << v(0) << " " << v(1) << " " << v(2) << "\n";
	}
	oss << "VECTORS Mass_flux float\n";
	for (size_t i=0; i<Ncells; ++i)
	{
		Vec3_t v; Cells[i]->Velocity(v);
		v *= Cells[i]->Rho;
		oss << v(0) << " " << v(1) << " " << v(2) << "\n";
	}

    String fn(FileKey);
    fn.append(".vtk");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline Cell * Lattice::GetCell(iVec3_t const & v)
{
    return Cells[v[0] + v[1]*Ndim[0] + v[2]*Ndim[0]*Ndim[1]];
}

inline void Lattice::Solve(double Tf, double dtOut, ptFun_t ptSetup, ptFun_t ptReport,
                           char const * FileKey, bool RenderVideo, size_t Nproc)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    idx_out     = 0;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        Collide();
        BounceBack();
        Stream();
        
        if (Time >= tout)
        {
            if (FileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%08d", FileKey, idx_out);
                WriteVTK     (fn.CStr());
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }

        Time += dt;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}
#endif
