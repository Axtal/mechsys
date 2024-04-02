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


#ifndef MECHSYS_MPM_DOMAIN_H
#define MECHSYS_MPM_DOMAIN_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// Std lib
#ifdef USE_OMP
#include <omp.h>
#endif

//STD
#include<iostream>

// boost => to read json files
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

// Mechsys
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/mpm/Node.h>
#include <mechsys/mpm/shape_functions.h>
#include <mechsys/mpm/Corner.h>
#include <mechsys/dem/domain.h>
#include <mechsys/mesh/unstructured.h>


namespace MPM
{

struct BdryFac // A structure for the boundary facets data
{
    Mesh::Cell *  Ce; //Given boundary cell
    size_t Facets[3]; //Boundary facets
    Vec3_t       Nor; //Normal vector pointing out
    Vec3_t      Xmin; //Bounding box
    Vec3_t      Xmax; //Bounding box
};

inline size_t Pt2idx(iVec3_t & iv, iVec3_t & Dim) // Calculates the index of the cell at coordinates iv for a cubic lattice of dimensions Dim
{
    return iv(0) + iv(1)*Dim(0) + iv(2)*Dim(0)*Dim(1);
}

inline void   idx2Pt(size_t n, iVec3_t & iv, iVec3_t & Dim) // Calculates the coordinates from the index
{
    iv(0) = n%Dim(0);
    iv(1) = (n/Dim(0))%(Dim(1));
    iv(2) = n/(Dim(0)*Dim(1));
}

class Domain
{
public:
    //Constructor
    Domain() 
    {
        Nproc = 1; 
        Gn = 0.001;
        BodyMesh = NULL;
        Time = 0.0;
    };
    //Domain(MPMethod Method, ///< Type of material to solve
    //iVec3_t Ndim            ///< Vector of integers with the dimensions
    //);

    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Data
    void *               UserData;              ///< User Data
    size_t                  Nproc;              ///< Number of procesors for OMP computation
    size_t                idx_out;              ///< The discrete time step for output
    String                FileKey;              ///< File Key for output files
    double                     Dx;              ///< Grid size in physical units
    double                     Gn;              ///< Dissipative constant 
    double                     Dt;              ///< Time step in physical units
    double                   Time;              ///< Total simulation time
    double                   Mmin;              ///< The minimun mass of particle
    Vec3_t                   Xlow;              ///< Lower bound of cubic grid
    Vec3_t                  Xhigh;              ///< Upper bound of cubic grid
    Node *                  Nodes;              ///< Array of Nodes
    iVec3_t                  Ndim;              ///< Array with grid dimensions    
    size_t                 Nnodes;              ///< Total number of nodes
    Array<Particle *>   Particles;              ///< Array of particles
    Array<size_t >         VNodes;              ///< Array of valid nodes
    Array<Corner *>       Corners;              ///< Array of corners of mesh
    
    //Mesh Data
    Mesh::Unstructured * BodyMesh;              ///< Unstructured mesh for detailed geometries
    Array<BdryFac>           BFac;              ///< Array with boundary data


    //Methods
    void AddFromJson (int Tag, char const * Filename, double rho, double scale, size_t nx); ///< Loads a Jason mesh into a set of MPM particles
    void AddFromJsonMesh(int Tag, char const * Filename, Vec3_t const & dis/*dispalcement*/, Vec3_t const & xr, Vec3_t const & axis, double ang, double rho, double scale, size_t nx); ///< Loads a Jason mesh into a set of MPM particles, but also build an utility mesh
    void AddFromOBJMesh(int Tag, char const * Filename, Vec3_t const & dis/*dispalcement*/, Vec3_t const & xr, Vec3_t const & axis, double ang, double rho, double scale, size_t nx); ///< Loads a OBJ mesh into a set of MPM particles, but also build an utility mesh
    void AddRectangularBeam(int Tag, Vec3_t const & Xmin, Vec3_t const & Xmax, double rho, size_t nx);   ///< Creates a solid rectangular beam of particles between Xmin and Xmax with a number of divisions nx per x lenght and density rho
    void AddRectangularBeamMesh(int Tag, Vec3_t const & Xmin, Vec3_t const & Xmax, double rho, size_t nx);   ///< Creates a solid rectangular beam of particles between Xmin and Xmax with a number of divisions nx per x lenght and density rho from a tet mesh
    void ResizeDomain(Vec3_t const & Xmin, Vec3_t const & Xmax, double deltax);    ///< Creates a mesh between Xmin and Xmax limits with deltax as the grid size
    double ResizeDomainMesh(Vec3_t const & Xmin, Vec3_t const & Xmax, double beta);   ///< Creates a mesh between Xmin and Xmax limits with beta the ratio between the largest particle and the cell size
    void NodePosition(size_t in, Vec3_t & Xn);                                      ///< Gives a posiiton of a node in the coordinate system
    void BoundingBox(Vec3_t & minX, Vec3_t & maxX);                                 ///< Give dimensions of the bounding box enclosing the particles
    void ParticleToNode();                                                          ///< Imprinting particle information into the mesh
    void NodeToParticle();                                                          ///< Transfering back node information into particles
    void OneStepUSF();                                                             ///< One time step in the USF formulation
    void OneStepCPI();                                                             ///< One time step in the CPI formulation
    void UpdateMesh();                                                             ///< Update the mesh geometry if it is present
    void BoundaryMesh();                                                             ///< Find the boundary triangles if it is present
    void BoundaryMeshUpdate();                                                      ///< Update the mesh boundanry if it is present
    void PosMesh(Vec3_t & pos);                                             ///< Position the center of the bounding box at a given place
    void Solve(double Tf, double dt, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);            ///< Solve the Domain dynamics

#ifdef USE_HDF5    
    void WriteXDMF         (char const * FileKey);                                      ///< Save a xdmf file for visualization
#endif
};

inline void Domain::AddFromJson (int Tag, char const * Filename, double rho, double scale, size_t nx)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Introducing Mesh --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    DEM::Domain demdom;
    demdom.AddFromJson(-1, Filename, 0.0, 1.0, 1.0);
    Vec3_t Xmin,Xmax;
    demdom.BoundingBox(Xmin,Xmax);
    double deltax = (Xmax(0)-Xmin(0))/nx;
    size_t Nx = nx;
    size_t Ny = (Xmax(1)-Xmin(1))/deltax;
    size_t Nz = (Xmax(2)-Xmin(2))/deltax;
    size_t ip = 0;
    for (size_t ix=0;ix<Nx;ix++)
    for (size_t iy=0;iy<Ny;iy++)
    for (size_t iz=0;iz<Nz;iz++)
    {
        Vec3_t X0 = Xmin + deltax*Vec3_t(0.5+ix,0.5+iy,0.5+iz);
        if (demdom.Particles[0]->IsInside(X0)&&demdom.Particles[0]->IsInsideAlt(X0))
        {
            Particles.Push(new Particle(Tag, X0, OrthoSys::O, rho*deltax*deltax*deltax, deltax*deltax*deltax));
            ip++;
        }
    }
    printf("%s  Num of material points   = %zd%s\n",TERM_CLR2,ip,TERM_RST);
}

inline void Domain::AddFromJsonMesh (int Tag, char const * Filename, Vec3_t const & dis/*displacement*/, Vec3_t const & xr, Vec3_t const & axis, double ang, double rho, double scale, size_t nx)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Introducing body from Mesh --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    Array<Vec3_t>       V; //Array of vertices
    Array<Array <int> > F; //Array of facets
    try {
        boost::property_tree::ptree pt;
        boost::property_tree::read_json(Filename, pt);
        BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("verts")) {
            Vec3_t coords;
            int i = 0;
            BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("c")) {
                coords[i] = scale * boost::lexical_cast<double>(b.second.data());
                i++;
            }
            V.Push(coords);
        }
        //BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("edges")) {
            //Array<int> vids(2);
            //int i = 0;
            //BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("verts")) {
                //vids[i] = boost::lexical_cast<int>(b.second.data());
                //i++;
            //}
            //E.Push(vids);
        //}
        BOOST_FOREACH(boost::property_tree::ptree::value_type & a, pt.get_child("faces")) {
            Array<int>     vids;
            BOOST_FOREACH(boost::property_tree::ptree::value_type & b, a.second.get_child("verts")) {
                int vid = boost::lexical_cast<int>(b.second.data());
                vids .Push(vid);
            }
            F.Push(vids);
        }
        printf("[1;32mMPM::Domain.h AddFromJson: finished[0m\n");
    } catch (std::exception & e) {
        throw new Fatal("MPM::domain.h: AddFromJson failed:\n\t%s", e.what());
    }

    Quaternion_t Q;
    NormalizeRotation(ang, axis, Q);

    for (size_t iv = 0; iv < V.Size(); iv++)
    {
        Vec3_t xt = V[iv]-xr;
        Rotation(xt,Q,V[iv]);
        V[iv] += xr+dis;
    }

    // approximated mass center
    // Vec3_t mc (0.,0.,0.);

    // for (size_t iv=0;iv<V.Size();iv++)  mc += V[iv];
    // mc /= V.Size();

    // std::cout << "mc: " << mc << std::endl;

    // for (size_t iv=0;iv<V.Size();iv++)  V[iv] -= mc;

    BodyMesh = new Mesh::Unstructured(3);
    BodyMesh->Set(V.Size()/*number of points to define the goemetry*/,F.Size()/*number of faces to define it*/,1/*number of regions*/,0/*number of holes*/);
    double xmax = V[0](0);
    double xmin = V[0](0);
    double ymax = V[0](1);
    double ymin = V[0](1);
    double zmax = V[0](2);
    double zmin = V[0](2);
    for (size_t iv=1;iv<V.Size();iv++)
    {
        if (xmax < V[iv](0)) xmax = V[iv](0);
        if (xmin > V[iv](0)) xmin = V[iv](0);
        if (ymax < V[iv](1)) ymax = V[iv](1);
        if (ymin > V[iv](1)) ymin = V[iv](1);
        if (zmax < V[iv](2)) zmax = V[iv](2);
        if (zmin > V[iv](2)) zmin = V[iv](2);
    }
    double Lx = std::max(xmax - xmin,std::max(ymax-ymin,zmax-zmin));
    double Vmax = Lx*Lx*Lx/(nx*nx*nx); //maximun volume of cells, controls resolution
    Vec3_t mid(0.5*(xmax+xmin),0.5*(ymax+ymin),0.5*(zmax+zmin));
    BodyMesh->SetReg(0,0,Vmax,mid(0)/*x*/,mid(1)/*y*/,mid(2)/*z*/); // The values of L must be apoint inside the region to be meshed

    for (size_t iv=0;iv<V.Size();iv++)
    {
        BodyMesh->SetPnt(iv/*index of the point*/,-1/*tag of the point*/,V[iv](0),V[iv](1),V[iv](2)); 
    }

    for (size_t ic=0;ic<F.Size();ic++)
    {
        BodyMesh->SetFac(ic, -2, F[ic]/*array of indexes of the points defined before*/);
    }
    BodyMesh->Generate();

    for (size_t iv=0; iv < BodyMesh->Verts.Size(); iv++)
    {   
         //BodyMesh->Verts[iv]->C(0) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));
         //BodyMesh->Verts[iv]->C(1) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));
         //BodyMesh->Verts[iv]->C(2) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));
         //BodyMesh->Verts[iv]->C += Vec3_t(0.5*Lx,0.5*Ly,0.5*Lz) + Xmin;

         Corners.Push(new Corner(&BodyMesh->Verts[iv]->C));
    }

    for (size_t i=0; i<BodyMesh->Cells.Size(); ++i)
    {
        Vec3_t x0 = BodyMesh->Cells[i]->V[0]->C;
        Vec3_t x1 = BodyMesh->Cells[i]->V[1]->C;
        Vec3_t x2 = BodyMesh->Cells[i]->V[2]->C;
        Vec3_t x3 = BodyMesh->Cells[i]->V[3]->C;
        Vec3_t C  = 0.25*(x0+x1+x2+x3);
        
        Vec3_t a  = x1 - x0;
        Vec3_t b  = x2 - x0;
        Vec3_t c  = x3 - x0;

        double V  = fabs(dot(a,cross(b,c)))/6.0; 

        Particles.Push(new MPM::Particle(Tag, C, OrthoSys::O, V*rho, V));

        Particles[Particles.Size()-1]->Vcorner[0] = &BodyMesh->Cells[i]->V[0]->C;
        Particles[Particles.Size()-1]->Vcorner[1] = &BodyMesh->Cells[i]->V[1]->C;
        Particles[Particles.Size()-1]->Vcorner[2] = &BodyMesh->Cells[i]->V[2]->C;
        Particles[Particles.Size()-1]->Vcorner[3] = &BodyMesh->Cells[i]->V[3]->C;

        size_t i0 = BodyMesh->Cells[i]->V[0]->ID;
        size_t i1 = BodyMesh->Cells[i]->V[1]->ID;
        size_t i2 = BodyMesh->Cells[i]->V[2]->ID;
        size_t i3 = BodyMesh->Cells[i]->V[3]->ID;

        Particles[Particles.Size()-1]->Corners[0] = i0;
        Particles[Particles.Size()-1]->Corners[1] = i1;
        Particles[Particles.Size()-1]->Corners[2] = i2;
        Particles[Particles.Size()-1]->Corners[3] = i3;

        Particles[Particles.Size()-1]->CalcVolCPI(0.0);

        Vec3_t v0 = BodyMesh->Verts[BodyMesh->Cells[i]->V[0]->ID]->C;
        Vec3_t v1 = BodyMesh->Verts[BodyMesh->Cells[i]->V[1]->ID]->C;
        Vec3_t v2 = BodyMesh->Verts[BodyMesh->Cells[i]->V[2]->ID]->C;
        Vec3_t v3 = BodyMesh->Verts[BodyMesh->Cells[i]->V[3]->ID]->C;
        BodyMesh->Cells[i]->Xmin(0) = std::min(std::min(v0(0),v1(0)),std::min(v2(0),v3(0)));
        BodyMesh->Cells[i]->Xmin(1) = std::min(std::min(v0(1),v1(1)),std::min(v2(1),v3(1)));
        BodyMesh->Cells[i]->Xmin(2) = std::min(std::min(v0(2),v1(2)),std::min(v2(2),v3(2)));
        BodyMesh->Cells[i]->Xmax(0) = std::max(std::max(v0(0),v1(0)),std::max(v2(0),v3(0)));
        BodyMesh->Cells[i]->Xmax(1) = std::max(std::max(v0(1),v1(1)),std::max(v2(1),v3(1)));
        BodyMesh->Cells[i]->Xmax(2) = std::max(std::max(v0(2),v1(2)),std::max(v2(2),v3(2)));
        
    }

    printf("%s  Num of material points   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
    printf("%s  Num of corner   points   = %zd%s\n",TERM_CLR2,Corners.Size()  ,TERM_RST);
}

inline void Domain::AddFromOBJMesh (int Tag, char const * Filename, Vec3_t const & dis/*displacement*/, Vec3_t const & xr, Vec3_t const & axis, double ang, double rho, double scale, size_t nx)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Introducing body from Mesh --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    Array<Vec3_t>       V; //Array of vertices
    Array<Array <int> > F; //Array of facets

    
    String fn(Filename);
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    ifstream fi(fn.CStr());
    
    std::string line;
    while(getline(fi,line))
    {
        std::istringstream ss( line );
        char c;
        double x,y,z;
	  	if (line[0]=='v' && line[1]==' ')
	  	{
	  		ss >> c >> x >> y >> z;
	        Vec3_t v(x,y,z);
            v *= scale;
            V.Push(v);
	  	}
	  	else if (line[0]=='f')
	  	{
            std::vector<std::string> ls;
            std::string lse;
			while (ss>>lse)	ls.push_back(lse);
            Array<int> face;
	        for (size_t m=1; m<ls.size(); ++m)
	        {
	        	 face.Push(stoi(ls[m])-1);
	        }
            F.Push(face);
        }
    }

    Quaternion_t Q;
    NormalizeRotation(ang, axis, Q);

    for (size_t iv = 0; iv < V.Size(); iv++)
    {
        Vec3_t xt = V[iv]-xr;
        Rotation(xt,Q,V[iv]);
        V[iv] += xr+dis;
    }

    // approximated mass center
    Vec3_t mc (0.,0.,0.);

    for (size_t iv=0;iv<V.Size();iv++)  mc += V[iv];
    mc /= V.Size();

    // std::cout << "mc: " << mc << std::endl;

    for (size_t iv=0;iv<V.Size();iv++)  V[iv] -= mc;

    BodyMesh = new Mesh::Unstructured(3);
    BodyMesh->Set(V.Size()/*number of points to define the goemetry*/,F.Size()/*number of faces to define it*/,1/*number of regions*/,0/*number of holes*/);
    double xmax = V[0](0);
    double xmin = V[0](0);
    double ymax = V[0](1);
    double ymin = V[0](1);
    double zmax = V[0](2);
    double zmin = V[0](2);
    for (size_t iv=1;iv<V.Size();iv++)
    {
        if (xmax < V[iv](0)) xmax = V[iv](0);
        if (xmin > V[iv](0)) xmin = V[iv](0);
        if (ymax < V[iv](1)) ymax = V[iv](1);
        if (ymin > V[iv](1)) ymin = V[iv](1);
        if (zmax < V[iv](2)) zmax = V[iv](2);
        if (zmin > V[iv](2)) zmin = V[iv](2);
    }
    double Lx = std::max(xmax - xmin,std::max(ymax-ymin,zmax-zmin));
    double Vmax = Lx*Lx*Lx/(nx*nx*nx); //maximun volume of cells, controls resolution
    Vec3_t mid(0.5*(xmax+xmin),0.5*(ymax+ymin),0.5*(zmax+zmin));
    BodyMesh->SetReg(0,0,Vmax,mid(0)/*x*/,mid(1)/*y*/,mid(2)/*z*/); // The values of L must be apoint inside the region to be meshed

    for (size_t iv=0;iv<V.Size();iv++)
    {
        BodyMesh->SetPnt(iv/*index of the point*/,-1/*tag of the point*/,V[iv](0),V[iv](1),V[iv](2)); 
    }

    for (size_t ic=0;ic<F.Size();ic++)
    {
        BodyMesh->SetFac(ic, -2, F[ic]/*array of indexes of the points defined before*/);
    }
    BodyMesh->Generate();

    for (size_t iv=0; iv < BodyMesh->Verts.Size(); iv++)
    {   
         //BodyMesh->Verts[iv]->C(0) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));
         //BodyMesh->Verts[iv]->C(1) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));
         //BodyMesh->Verts[iv]->C(2) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));
         //BodyMesh->Verts[iv]->C += Vec3_t(0.5*Lx,0.5*Ly,0.5*Lz) + Xmin;
        BodyMesh->Verts[iv]->C += mc;
        Corners.Push(new Corner(&BodyMesh->Verts[iv]->C));
    }

    for (size_t i=0; i<BodyMesh->Cells.Size(); ++i)
    {
        Vec3_t x0 = BodyMesh->Cells[i]->V[0]->C;
        Vec3_t x1 = BodyMesh->Cells[i]->V[1]->C;
        Vec3_t x2 = BodyMesh->Cells[i]->V[2]->C;
        Vec3_t x3 = BodyMesh->Cells[i]->V[3]->C;
        Vec3_t C  = 0.25*(x0+x1+x2+x3);
        
        Vec3_t a  = x1 - x0;
        Vec3_t b  = x2 - x0;
        Vec3_t c  = x3 - x0;

        double V  = fabs(dot(a,cross(b,c)))/6.0; 

        Particles.Push(new MPM::Particle(Tag, C, OrthoSys::O, V*rho, V));

        Particles[Particles.Size()-1]->Vcorner[0] = &BodyMesh->Cells[i]->V[0]->C;
        Particles[Particles.Size()-1]->Vcorner[1] = &BodyMesh->Cells[i]->V[1]->C;
        Particles[Particles.Size()-1]->Vcorner[2] = &BodyMesh->Cells[i]->V[2]->C;
        Particles[Particles.Size()-1]->Vcorner[3] = &BodyMesh->Cells[i]->V[3]->C;

        size_t i0 = BodyMesh->Cells[i]->V[0]->ID;
        size_t i1 = BodyMesh->Cells[i]->V[1]->ID;
        size_t i2 = BodyMesh->Cells[i]->V[2]->ID;
        size_t i3 = BodyMesh->Cells[i]->V[3]->ID;

        Particles[Particles.Size()-1]->Corners[0] = i0;
        Particles[Particles.Size()-1]->Corners[1] = i1;
        Particles[Particles.Size()-1]->Corners[2] = i2;
        Particles[Particles.Size()-1]->Corners[3] = i3;

        Particles[Particles.Size()-1]->CalcVolCPI(0.0);

        Vec3_t v0 = BodyMesh->Verts[BodyMesh->Cells[i]->V[0]->ID]->C;
        Vec3_t v1 = BodyMesh->Verts[BodyMesh->Cells[i]->V[1]->ID]->C;
        Vec3_t v2 = BodyMesh->Verts[BodyMesh->Cells[i]->V[2]->ID]->C;
        Vec3_t v3 = BodyMesh->Verts[BodyMesh->Cells[i]->V[3]->ID]->C;
        BodyMesh->Cells[i]->Xmin(0) = std::min(std::min(v0(0),v1(0)),std::min(v2(0),v3(0)));
        BodyMesh->Cells[i]->Xmin(1) = std::min(std::min(v0(1),v1(1)),std::min(v2(1),v3(1)));
        BodyMesh->Cells[i]->Xmin(2) = std::min(std::min(v0(2),v1(2)),std::min(v2(2),v3(2)));
        BodyMesh->Cells[i]->Xmax(0) = std::max(std::max(v0(0),v1(0)),std::max(v2(0),v3(0)));
        BodyMesh->Cells[i]->Xmax(1) = std::max(std::max(v0(1),v1(1)),std::max(v2(1),v3(1)));
        BodyMesh->Cells[i]->Xmax(2) = std::max(std::max(v0(2),v1(2)),std::max(v2(2),v3(2)));
        
    }

    printf("%s  Num of material points   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
    printf("%s  Num of corner   points   = %zd%s\n",TERM_CLR2,Corners.Size()  ,TERM_RST);
}

inline void Domain::AddRectangularBeam(int Tag, Vec3_t const & Xmin, Vec3_t const & Xmax, double rho, size_t nx)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Introducing rectangular Beam --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    double deltax = (Xmax(0)-Xmin(0))/nx;
    size_t Nx = nx;
    size_t Ny = (Xmax(1)-Xmin(1))/deltax;
    size_t Nz = (Xmax(2)-Xmin(2))/deltax;
    size_t Np = Nx*Ny*Nz;
    //double Vol = (Xmax(0)-Xmin(0))*(Xmax(1)-Xmin(1))*(Xmax(2)-Xmin(2));
    double Vol = Np*deltax*deltax*deltax;
    double Mp = rho*Vol/Np;
    size_t ip = 0;
    for (size_t ix=0;ix<Nx;ix++)
    for (size_t iy=0;iy<Ny;iy++)
    for (size_t iz=0;iz<Nz;iz++)
    {
        Vec3_t X0 = Xmin + deltax*Vec3_t(0.5+ix,0.5+iy,0.5+iz);
        Particles.Push(new Particle(Tag, X0, OrthoSys::O, Mp, Vol/Np));
        ip++;
    }
    printf("%s  Num of material points   = %zd%s\n",TERM_CLR2,Np,TERM_RST);
}

inline void Domain::AddRectangularBeamMesh(int Tag, Vec3_t const & Xmin, Vec3_t const & Xmax, double rho, size_t nx)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Introducing rectangular Beam --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    BodyMesh = new Mesh::Unstructured(3);
    BodyMesh->Set(8/*number of points to define the goemetry*/,6/*number of faces to define it*/,1/*number of regions*/,0/*number of holes*/);
    double Lx = Xmax(0)-Xmin(0);
    double Ly = Xmax(1)-Xmin(1);
    double Lz = Xmax(2)-Xmin(2);
    double L  = std::max(Lx,std::max(Ly,Lz));
    double Vmax = L*L*L/(nx*nx*nx); //maximun volume of cells, controls resolution
    BodyMesh->SetReg(0,0,Vmax,0.5*Lx/*x*/,0.5*Ly/*y*/,0.5*Lz/*z*/); // The values of L must be apoint inside the region to be meshed
    Vec3_t x00 = Vec3_t(-0.5*Lx,-0.5*Ly,-0.5*Lz);
    Vec3_t x01 = Vec3_t( 0.5*Lx,-0.5*Ly,-0.5*Lz);
    Vec3_t x02 = Vec3_t( 0.5*Lx,-0.5*Ly, 0.5*Lz);
    Vec3_t x03 = Vec3_t(-0.5*Lx,-0.5*Ly, 0.5*Lz);
    Vec3_t x04 = Vec3_t(-0.5*Lx, 0.5*Ly,-0.5*Lz);
    Vec3_t x05 = Vec3_t( 0.5*Lx, 0.5*Ly,-0.5*Lz);
    Vec3_t x06 = Vec3_t( 0.5*Lx, 0.5*Ly, 0.5*Lz);
    Vec3_t x07 = Vec3_t(-0.5*Lx, 0.5*Ly, 0.5*Lz);
    BodyMesh->SetPnt( 0/*index of the point*/,-1/*tag of the point*/,x00(0),x00(1),x00(2)); 
    BodyMesh->SetPnt( 1/*index of the point*/,-1/*tag of the point*/,x01(0),x01(1),x01(2)); 
    BodyMesh->SetPnt( 2/*index of the point*/,-1/*tag of the point*/,x02(0),x02(1),x02(2)); 
    BodyMesh->SetPnt( 3/*index of the point*/,-1/*tag of the point*/,x03(0),x03(1),x03(2)); 
    BodyMesh->SetPnt( 4/*index of the point*/,-1/*tag of the point*/,x04(0),x04(1),x04(2)); 
    BodyMesh->SetPnt( 5/*index of the point*/,-1/*tag of the point*/,x05(0),x05(1),x05(2)); 
    BodyMesh->SetPnt( 6/*index of the point*/,-1/*tag of the point*/,x06(0),x06(1),x06(2)); 
    BodyMesh->SetPnt( 7/*index of the point*/,-1/*tag of the point*/,x07(0),x07(1),x07(2)); 
    BodyMesh->SetFac( 0, -2, Array<int>(0,1,2,3)/*array of indexes of the points defined before*/);
    BodyMesh->SetFac( 1, -2, Array<int>(4,7,6,5)/*array of indexes of the points defined before*/);
    BodyMesh->SetFac( 2, -2, Array<int>(0,3,7,4)/*array of indexes of the points defined before*/);
    BodyMesh->SetFac( 3, -2, Array<int>(1,5,6,2)/*array of indexes of the points defined before*/);
    BodyMesh->SetFac( 4, -2, Array<int>(0,4,5,1)/*array of indexes of the points defined before*/);
    BodyMesh->SetFac( 5, -2, Array<int>(3,2,6,7)/*array of indexes of the points defined before*/);
    BodyMesh->Generate();

    for (size_t iv=0; iv < BodyMesh->Verts.Size(); iv++)
    {   
         BodyMesh->Verts[iv]->C(0) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));
         BodyMesh->Verts[iv]->C(1) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));
         BodyMesh->Verts[iv]->C(2) *= (1.0+1.0e-4*(rand()/RAND_MAX-0.5));

         BodyMesh->Verts[iv]->C += Vec3_t(0.5*Lx,0.5*Ly,0.5*Lz) + Xmin;

         Corners.Push(new Corner(&BodyMesh->Verts[iv]->C));
    }

    for (size_t i=0; i<BodyMesh->Cells.Size(); ++i)
    {
        Vec3_t x0 = BodyMesh->Cells[i]->V[0]->C;
        Vec3_t x1 = BodyMesh->Cells[i]->V[1]->C;
        Vec3_t x2 = BodyMesh->Cells[i]->V[2]->C;
        Vec3_t x3 = BodyMesh->Cells[i]->V[3]->C;
        Vec3_t C  = 0.25*(x0+x1+x2+x3);
        
        Vec3_t a  = x1 - x0;
        Vec3_t b  = x2 - x0;
        Vec3_t c  = x3 - x0;

        double V  = fabs(dot(a,cross(b,c)))/6.0; 

        Particles.Push(new MPM::Particle(Tag, C, OrthoSys::O, V*rho, V));

        Particles[Particles.Size()-1]->Vcorner[0] = &BodyMesh->Cells[i]->V[0]->C;
        Particles[Particles.Size()-1]->Vcorner[1] = &BodyMesh->Cells[i]->V[1]->C;
        Particles[Particles.Size()-1]->Vcorner[2] = &BodyMesh->Cells[i]->V[2]->C;
        Particles[Particles.Size()-1]->Vcorner[3] = &BodyMesh->Cells[i]->V[3]->C;

        size_t i0 = BodyMesh->Cells[i]->V[0]->ID;
        size_t i1 = BodyMesh->Cells[i]->V[1]->ID;
        size_t i2 = BodyMesh->Cells[i]->V[2]->ID;
        size_t i3 = BodyMesh->Cells[i]->V[3]->ID;

        Particles[Particles.Size()-1]->Corners[0] = i0;
        Particles[Particles.Size()-1]->Corners[1] = i1;
        Particles[Particles.Size()-1]->Corners[2] = i2;
        Particles[Particles.Size()-1]->Corners[3] = i3;

        Particles[Particles.Size()-1]->CalcVolCPI(0.0);

        Vec3_t v0 = BodyMesh->Verts[BodyMesh->Cells[i]->V[0]->ID]->C;
        Vec3_t v1 = BodyMesh->Verts[BodyMesh->Cells[i]->V[1]->ID]->C;
        Vec3_t v2 = BodyMesh->Verts[BodyMesh->Cells[i]->V[2]->ID]->C;
        Vec3_t v3 = BodyMesh->Verts[BodyMesh->Cells[i]->V[3]->ID]->C;
        BodyMesh->Cells[i]->Xmin(0) = std::min(std::min(v0(0),v1(0)),std::min(v2(0),v3(0)));
        BodyMesh->Cells[i]->Xmin(1) = std::min(std::min(v0(1),v1(1)),std::min(v2(1),v3(1)));
        BodyMesh->Cells[i]->Xmin(2) = std::min(std::min(v0(2),v1(2)),std::min(v2(2),v3(2)));
        BodyMesh->Cells[i]->Xmax(0) = std::max(std::max(v0(0),v1(0)),std::max(v2(0),v3(0)));
        BodyMesh->Cells[i]->Xmax(1) = std::max(std::max(v0(1),v1(1)),std::max(v2(1),v3(1)));
        BodyMesh->Cells[i]->Xmax(2) = std::max(std::max(v0(2),v1(2)),std::max(v2(2),v3(2)));
        
    }

    printf("%s  Num of material points   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
    printf("%s  Num of corner   points   = %zd%s\n",TERM_CLR2,Corners.Size()  ,TERM_RST);
}

inline void Domain::ResizeDomain(Vec3_t const & Xmin, Vec3_t const & Xmax, double deltax)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing MPM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    double Lx = (Xmax(0)-Xmin(0));
    double Ly = (Xmax(1)-Xmin(1));
    double Lz = (Xmax(2)-Xmin(2));
    if (Lx<deltax||Ly<deltax||Lz<deltax)
    {
        throw new Fatal("Recheck mesh dimensions and grid size since no cells can be produced");
    }
    size_t Nx = Lx/deltax;
    size_t Ny = Ly/deltax;
    size_t Nz = Lz/deltax;
    Ndim      = iVec3_t(Nx,Ny,Nz);
    Nnodes    = Nx*Ny*Nz;
    Dx        = deltax;
    Nodes     = new Node [Nnodes];
    for (size_t iz=0;iz<Nz;iz++)
    for (size_t iy=0;iy<Ny;iy++)
    for (size_t ix=0;ix<Nx;ix++)
    {
        iVec3_t Pt(ix,iy,iz);
        size_t idx = Pt2idx(Pt,Ndim);
        Nodes[idx].X   = Pt;
        Nodes[idx].Idx = idx;
    }

    //Assigning domain limits
    Xlow  = Xmin;
    Xhigh = Xmax;
   
    printf("%s  Num of nodes   = %zd%s\n",TERM_CLR2,Nx*Ny*Nz,TERM_RST);
}

inline double Domain::ResizeDomainMesh(Vec3_t const & Xmin, Vec3_t const & Xmax, double beta)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing MPM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    double Lx = (Xmax(0)-Xmin(0));
    double Ly = (Xmax(1)-Xmin(1));
    double Lz = (Xmax(2)-Xmin(2));

    if (BodyMesh==NULL)
    {
        throw new Fatal("MPM::ResizeDomainMesh: this function can only be used for the CPDI version of the code");
    }

    double deltax = 0.0;
    double mindx  = BodyMesh->Cells[0]->Xmax(0) - BodyMesh->Cells[0]->Xmin(0);

	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        double dx = BodyMesh->Cells[ip]->Xmax(0) - BodyMesh->Cells[ip]->Xmin(0);
        double dy = BodyMesh->Cells[ip]->Xmax(1) - BodyMesh->Cells[ip]->Xmin(1);
        double dz = BodyMesh->Cells[ip]->Xmax(2) - BodyMesh->Cells[ip]->Xmin(2);
        double dmax = std::max(dx,std::max(dy,dz));
        double dmin = std::min(dx,std::min(dy,dz));
        if (deltax < beta*dmax) deltax = beta*dmax;
        if (mindx  > dmin) mindx = dmin;
    }

    if (Lx<deltax||Ly<deltax||Lz<deltax)
    {
        throw new Fatal("Recheck mesh dimensions and grid size since no cells can be produced");
    }
    size_t Nx = Lx/deltax;
    size_t Ny = Ly/deltax;
    size_t Nz = Lz/deltax;
    Ndim      = iVec3_t(Nx,Ny,Nz);
    Nnodes    = Nx*Ny*Nz;
    Dx        = deltax;
    Nodes     = new Node [Nnodes];
    for (size_t iz=0;iz<Nz;iz++)
    for (size_t iy=0;iy<Ny;iy++)
    for (size_t ix=0;ix<Nx;ix++)
    {
        iVec3_t Pt(ix,iy,iz);
        size_t idx = Pt2idx(Pt,Ndim);
        Nodes[idx].X   = Pt;
        Nodes[idx].Idx = idx;
    }

    //Assigning domain limits
    Xlow  = Xmin;
    Xhigh = Xmax;
   
    printf("%s  Num of nodes     = %zd%s\n",TERM_CLR2,Nx*Ny*Nz,TERM_RST);
    printf("%s  Nx = %zd Ny = %zd, Nz = %zd %s\n",TERM_CLR2,Nx,Ny,Nz,TERM_RST);
    printf("%s  Grid spacing     = %g%s\n" ,TERM_CLR2,deltax,TERM_RST);
    printf("%s  Smallest element size     = %g%s\n" ,TERM_CLR2,mindx,TERM_RST);
    printf("%s  Grid lower bound = %g %g %g %s \n",TERM_CLR2,Xmin(0),Xmin(1),Xmin(2),TERM_RST);
    printf("%s  Grid upper bound = %g %g %g %s \n",TERM_CLR2,Xmax(0),Xmax(1),Xmax(2),TERM_RST);

    return mindx;

}

inline void Domain::NodePosition(size_t in, Vec3_t & Xn)
{
    iVec3_t Pt;
    idx2Pt(in,Pt,Ndim);
    Xn = Vec3_t(Dx*Pt(0),Dx*Pt(1),Dx*Pt(2));
    Xn += Xlow;
}

inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    if (Particles.Size()==0) throw new Fatal("MPM::Domain::BoundingBox: There are no particles to build the bounding box");
    
    if (BodyMesh==NULL)
    {
        minX = Vec3_t(Particles[0]->x(0)-Particles[0]->Dx(0), Particles[0]->x(1)-Particles[0]->Dx(1), Particles[0]->x(2)-Particles[0]->Dx(2));
        maxX = Vec3_t(Particles[0]->x(0)+Particles[0]->Dx(0), Particles[0]->x(1)+Particles[0]->Dx(1), Particles[0]->x(2)+Particles[0]->Dx(2));
        for (size_t i=1; i<Particles.Size(); i++)
        {
            if (minX(0)>Particles[i]->x(0)-Particles[i]->Dx(0)) minX(0) = Particles[i]->x(0)-Particles[i]->Dx(0);
            if (maxX(0)<Particles[i]->x(0)+Particles[i]->Dx(0)) maxX(0) = Particles[i]->x(0)+Particles[i]->Dx(0);
            if (minX(1)>Particles[i]->x(1)-Particles[i]->Dx(1)) minX(1) = Particles[i]->x(1)-Particles[i]->Dx(1);
            if (maxX(1)<Particles[i]->x(1)+Particles[i]->Dx(1)) maxX(1) = Particles[i]->x(1)+Particles[i]->Dx(1);
            if (minX(2)>Particles[i]->x(2)-Particles[i]->Dx(2)) minX(2) = Particles[i]->x(2)-Particles[i]->Dx(2);
            if (maxX(2)<Particles[i]->x(2)+Particles[i]->Dx(2)) maxX(2) = Particles[i]->x(2)+Particles[i]->Dx(2);
        }
    }
    else
    {
        minX = *Corners[0]->x;
        maxX = *Corners[0]->x;
        for (size_t ic=1;ic<Corners.Size();ic++)
        {
            if (minX(0)> (*Corners[ic]->x)(0)) minX(0) = (*Corners[ic]->x)(0);
            if (maxX(0)< (*Corners[ic]->x)(0)) maxX(0) = (*Corners[ic]->x)(0);
            if (minX(1)> (*Corners[ic]->x)(1)) minX(1) = (*Corners[ic]->x)(1);
            if (maxX(1)< (*Corners[ic]->x)(1)) maxX(1) = (*Corners[ic]->x)(1);
            if (minX(2)> (*Corners[ic]->x)(2)) minX(2) = (*Corners[ic]->x)(2);
            if (maxX(2)< (*Corners[ic]->x)(2)) maxX(2) = (*Corners[ic]->x)(2);
        }
    }
}

#ifdef USE_HDF5    
inline void Domain::WriteXDMF (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Storing center of mass data
    float * Posvec  = new float[3*Particles.Size()];
    float * Velvec  = new float[3*Particles.Size()];
    int   * Tag     = new int  [  Particles.Size()];
    float * Sigma	= new float[6*Particles.Size()];
    float * Strain	= new float[6*Particles.Size()];

    for (size_t i=0;i<Particles.Size();i++)
    {

        Posvec[3*i  ] = float(Particles[i]->x(0));
        Posvec[3*i+1] = float(Particles[i]->x(1));
        Posvec[3*i+2] = float(Particles[i]->x(2));
        Velvec[3*i  ] = float(Particles[i]->v(0));
        Velvec[3*i+1] = float(Particles[i]->v(1));
        Velvec[3*i+2] = float(Particles[i]->v(2));
        Sigma [6*i  ] = float(Particles[i]->S(0,0));
        Sigma [6*i+1] = float(Particles[i]->S(0,1));
        Sigma [6*i+2] = float(Particles[i]->S(0,2));
        Sigma [6*i+3] = float(Particles[i]->S(1,1));
        Sigma [6*i+4] = float(Particles[i]->S(1,2));
        Sigma [6*i+5] = float(Particles[i]->S(2,2));
        Strain[6*i  ] = float(Particles[i]->E(0,0));
        Strain[6*i+1] = float(Particles[i]->E(0,1));
        Strain[6*i+2] = float(Particles[i]->E(0,2));
        Strain[6*i+3] = float(Particles[i]->E(1,1));
        Strain[6*i+4] = float(Particles[i]->E(1,2));
        Strain[6*i+5] = float(Particles[i]->E(2,2));
        Tag   [i]     = int  (Particles[i]->Tag);  
    }

    hsize_t dims[1];
    dims[0] = 3*Particles.Size();
    String dsname;
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("PVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dims[0] = Particles.Size();
    dsname.Printf("PTag");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tag   );
    dims[0] = 6*Particles.Size();
    dsname.Printf("Sigma");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma);
    dsname.Printf("Strain");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Strain);


    delete [] Posvec;
    delete [] Velvec;
    delete [] Tag;
    delete [] Sigma;
    delete [] Strain;

    if (BodyMesh != NULL)
    {
        int   * Con     = new int  [4*BodyMesh->Cells.Size()];
        float * Verts   = new float[3*BodyMesh->Verts.Size()];
        int   * TagM    = new int  [  BodyMesh->Cells.Size()];
        float * SigmaM	= new float[6*BodyMesh->Cells.Size()];
        float * StrainM	= new float[6*BodyMesh->Cells.Size()];
        
        //Saving the vertices of the tetrahedron
        for (size_t iv=0;iv<BodyMesh->Verts.Size();iv++)
        {
            Verts[3*iv+0] = BodyMesh->Verts[iv]->C(0);
            Verts[3*iv+1] = BodyMesh->Verts[iv]->C(1);
            Verts[3*iv+2] = BodyMesh->Verts[iv]->C(2);
        }

        //Saving the cell information
        for (size_t ic=0;ic<BodyMesh->Cells.Size();ic++)
        {
            TagM[ic] = Particles[ic]->Tag;
            SigmaM [6*ic  ] = float(Particles[ic]->S(0,0));
            SigmaM [6*ic+1] = float(Particles[ic]->S(0,1));
            SigmaM [6*ic+2] = float(Particles[ic]->S(0,2));
            SigmaM [6*ic+3] = float(Particles[ic]->S(1,1));
            SigmaM [6*ic+4] = float(Particles[ic]->S(1,2));
            SigmaM [6*ic+5] = float(Particles[ic]->S(2,2));
            StrainM[6*ic  ] = float(Particles[ic]->E(0,0));
            StrainM[6*ic+1] = float(Particles[ic]->E(0,1));
            StrainM[6*ic+2] = float(Particles[ic]->E(0,2));
            StrainM[6*ic+3] = float(Particles[ic]->E(1,1));
            StrainM[6*ic+4] = float(Particles[ic]->E(1,2));
            StrainM[6*ic+5] = float(Particles[ic]->E(2,2));
            Con    [4*ic  ] = BodyMesh->Cells[ic]->V[0]->ID;
            Con    [4*ic+1] = BodyMesh->Cells[ic]->V[1]->ID;
            Con    [4*ic+2] = BodyMesh->Cells[ic]->V[2]->ID;
            Con    [4*ic+3] = BodyMesh->Cells[ic]->V[3]->ID;
        }
        dims[0] = 3*BodyMesh->Verts.Size();
        dsname.Printf("Vertex");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Verts);
        dims[0] = 4*BodyMesh->Cells.Size();
        dsname.Printf("Con");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Con);
        dims[0] = BodyMesh->Cells.Size();
        dsname.Printf("MTag");
        H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,TagM);
        dims[0] = 6*BodyMesh->Cells.Size();
        dsname.Printf("MSigma");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,SigmaM);
        dsname.Printf("MStrain");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,StrainM);

        delete [] Con    ;
        delete [] Verts  ; 
        delete [] TagM   ; 
        delete [] SigmaM ; 
        delete [] StrainM; 
    }


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
    

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"MPM_points\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Sigma\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Sigma \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Strain\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Strain \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    if (BodyMesh != NULL)
    {
    oss << "   <Grid Name=\"MPM_mesh\">\n";
    oss << "     <Topology TopologyType=\"Tetrahedron\" NumberOfElements=\"" << BodyMesh->Cells.Size() << "\">\n";
    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << BodyMesh->Cells.Size() << " 4\">\n";
    oss << "        " << fn.CStr() <<":/Con \n";
    oss << "       </DataItem>\n";
    oss << "     </Topology>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << BodyMesh->Verts.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Vertex \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << BodyMesh->Cells.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/MTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    //oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    //oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    //oss << "        " << fn.CStr() <<":/MVelocity\n";
    //oss << "       </DataItem>\n";
    //oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Sigma\" AttributeType=\"Tensor6\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << BodyMesh->Cells.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/MSigma \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Strain\" AttributeType=\"Tensor6\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << BodyMesh->Cells.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/MStrain \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
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

inline void Domain::ParticleToNode()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].Reset();
    }
    VNodes.Resize(0);

    //Transfering info from particle to node
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        iVec3_t nlow,nhigh;
        Vec3_t Xp = Particles[ip]->x;
        Vec3_t Lp = Particles[ip]->Dx;
        Mat3_t Vsp = -Particles[ip]->V*Particles[ip]->S;
        //Vec3_t Fex = Particles[ip]->b + Particles[ip]->h;
        Vec3_t Fex = Particles[ip]->b;

        //nlow (0) = (size_t(trunc((Xp(0)-Lp(0)-Xlow(0))/Dx)-1+Ndim(0)))%Ndim(0);
        //nlow (1) = (size_t(trunc((Xp(1)-Lp(1)-Xlow(1))/Dx)-1+Ndim(1)))%Ndim(1);
        //nlow (2) = (size_t(trunc((Xp(2)-Lp(2)-Xlow(2))/Dx)-1+Ndim(2)))%Ndim(2);
        //nhigh(0) = (size_t(ceil ((Xp(0)+Lp(0)-Xlow(0))/Dx)+1+Ndim(0)))%Ndim(0);
        //nhigh(1) = (size_t(ceil ((Xp(1)+Lp(1)-Xlow(1))/Dx)+1+Ndim(1)))%Ndim(1);
        //nhigh(2) = (size_t(ceil ((Xp(2)+Lp(2)-Xlow(2))/Dx)+1+Ndim(2)))%Ndim(2);

        nlow (0) = (size_t)std::max((trunc((Xp(0)-0.5*Lp(0)-Xlow(0))/Dx)-1),0.0);
        nlow (1) = (size_t)std::max((trunc((Xp(1)-0.5*Lp(1)-Xlow(1))/Dx)-1),0.0);
        nlow (2) = (size_t)std::max((trunc((Xp(2)-0.5*Lp(2)-Xlow(2))/Dx)-1),0.0);
        nhigh(0) = (size_t)std::min((ceil ((Xp(0)+0.5*Lp(0)-Xlow(0))/Dx)+1),(double)Ndim(0)-1.0);
        nhigh(1) = (size_t)std::min((ceil ((Xp(1)+0.5*Lp(1)-Xlow(1))/Dx)+1),(double)Ndim(1)-1.0);
        nhigh(2) = (size_t)std::min((ceil ((Xp(2)+0.5*Lp(2)-Xlow(2))/Dx)+1),(double)Ndim(2)-1.0);

        for (size_t nx=nlow(0); nx<=nhigh(0); nx++)
        for (size_t ny=nlow(1); ny<=nhigh(1); ny++)
        for (size_t nz=nlow(2); nz<=nhigh(2); nz++)
        {
            double Sf=0.0; //Shape function value
            Vec3_t Gs=OrthoSys::O; //Shape function gradient
            //Vec3_t Xn((nx+0.5)*Dx,(ny+0.5)*Dx,(nz+0.5)*Dx);
            Vec3_t Xn(nx*Dx,ny*Dx,nz*Dx);
            Xn += Xlow;
            iVec3_t Pt(nx,ny,nz);
            GIMP3D(Xp,Xn,Dx,Lp,Sf,Gs);
            if (Sf>0.0)
            {
                size_t nc = Pt2idx(Pt,Ndim);
                Particles[ip]->Nodes.Push(nc);
                Particles[ip]->Sfunc.Push(Sf);
                Particles[ip]->GSf  .Push(Gs);
                double ms = Sf*Particles[ip]->m;
                Vec3_t VspGs;
                Mult(Vsp,Gs,VspGs);
                Vec3_t dF = Sf*Fex + VspGs;
                omp_set_lock(&Nodes[nc].lck);
                Nodes[nc].Mass += ms;
                Nodes[nc].Mn   += ms*Particles[ip]->v + Dt*dF;
                Nodes[nc].Fn   += dF;
                //Nodes[nc].Sn    = Nodes[nc].Sn + ms*Particles[ip]->S;
                Nodes[nc].valid = true;
                omp_unset_lock(&Nodes[nc].lck);
                if (isnan(norm(Particles[ip]->x))||isnan(norm(Particles[ip]->v))||isnan(norm(Fex)))
                {
                    std::cout << "nan value detected" << std::endl;
                    std::cout << "ip =" << ip <<std::endl;
                    std::cout << "M =" << ms               <<std::endl;
                    std::cout << "v =" << Particles[ip]->v <<std::endl;
                    std::cout << "x =" << Particles[ip]->x <<std::endl;
                    std::cout << "dF=" << dF  <<std::endl;
                    std::cout << "fh =" << Particles[ip]->h <<std::endl;
                    std::cout << "Mn=" << Nodes[nc].Mn <<std::endl;
                    std::cout << "sf=" << Sf <<std::endl;
                    std::cout << "Time =" << Time <<std::endl;
                    throw new Fatal("Domain::ParticletoNode nan value found");
                }
            }
        }
        if (Particles[ip]->Nodes.Size()<2)
        {
            std::cout << "Particle does not have enough nodes" << std::endl;
            std::cout << "Time " << Time << std::endl;
            std::cout << "Ip   " << ip   << std::endl;
            std::cout << "Position  " << Particles[ip]->x << std::endl;
            std::cout << "Dx  " << Particles[ip]->Dx << std::endl;
            std::cout << "Dx0  " << Particles[ip]->Dx0 << std::endl;
            std::cout << "F  " << Particles[ip]->F << std::endl;
            std::cout << "Free? =" << Particles[ip]->IsFree() <<std::endl;
            throw new Fatal("MPM:Domain:ParticleToNode No mesh cells intersected by particles");
        }
    }

    //Producing an array of only valid nodes
    for (size_t in=0; in<Nnodes; in++)
    {
        if (Nodes[in].valid) VNodes.Push(Nodes[in].Idx);
    }
}

inline void Domain::NodeToParticle()
{
    //Solving the equations at the nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].UpdateNode(Gn,Dt,Mmin);
        Nodes[in].FixNode();
    }

    //Transfering info from node to particle
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        set_to_zero(Particles[ip]->L);
        //Particles[ip]->x0 = Particles[ip]->x;
        for (size_t nn=0; nn<Particles[ip]->Nodes.Size(); nn++)
        {
            size_t nc = Particles[ip]->Nodes[nn];
            double sf = Particles[ip]->Sfunc[nn];
            Vec3_t gs = Particles[ip]->GSf[nn];
            Vec3_t an = Nodes[nc].Fn/Nodes[nc].Mass;
            Particles[ip]->v += sf*an*Dt;
            Particles[ip]->x += sf*Nodes[nc].Vn*Dt;
            if (isnan(norm(Particles[ip]->x))||isnan(norm(Particles[ip]->v)))
            {
                std::cout << "nan value detected" << std::endl;
                std::cout << "ip =" << ip <<std::endl;
                std::cout << "Free? =" << Particles[ip]->IsFree() <<std::endl;
                std::cout << "v =" << Particles[ip]->v <<std::endl;
                std::cout << "x =" << Particles[ip]->x <<std::endl;
                std::cout << "fh =" << Particles[ip]->h <<std::endl;
                std::cout << "Vn =" << Nodes[nc].Vn <<std::endl;
                std::cout << "Mn =" << Nodes[nc].Mn <<std::endl;
                std::cout << "Ma =" << Nodes[nc].Mass <<std::endl;
                std::cout << "an =" << an <<std::endl;
                std::cout << "sf =" << sf <<std::endl;
                std::cout << "Time =" << Time <<std::endl;
                throw new Fatal("Domain::NodetoParticle nan value found");
            }
            Mat3_t Gv;
            //Dyad(Nodes[nc].Vn,gs,Gv);
            Dyad(gs,Nodes[nc].Vn,Gv);
            Particles[ip]->L = Particles[ip]->L + Gv;
        }
        Particles[ip]->CalcVol(Dt);
        Particles[ip]->CalcStress(Dt);
        Particles[ip]->Reset(Dt);
    }

}

inline void Domain::OneStepUSF()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].Reset();
    }
    VNodes.Resize(0);
    
    //Transfering info from particle to node
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        iVec3_t nlow,nhigh;
        Vec3_t Xp = Particles[ip]->x;
        Vec3_t Lp = Particles[ip]->Dx;

        //nlow (0) = (size_t(trunc((Xp(0)-Lp(0)-Xlow(0))/Dx)-1+Ndim(0)))%Ndim(0);
        //nlow (1) = (size_t(trunc((Xp(1)-Lp(1)-Xlow(1))/Dx)-1+Ndim(1)))%Ndim(1);
        //nlow (2) = (size_t(trunc((Xp(2)-Lp(2)-Xlow(2))/Dx)-1+Ndim(2)))%Ndim(2);
        //nhigh(0) = (size_t(ceil ((Xp(0)+Lp(0)-Xlow(0))/Dx)+1+Ndim(0)))%Ndim(0);
        //nhigh(1) = (size_t(ceil ((Xp(1)+Lp(1)-Xlow(1))/Dx)+1+Ndim(1)))%Ndim(1);
        //nhigh(2) = (size_t(ceil ((Xp(2)+Lp(2)-Xlow(2))/Dx)+1+Ndim(2)))%Ndim(2);

        nlow (0) = (size_t)std::max((trunc((Xp(0)-0.5*Lp(0)-Xlow(0))/Dx)-1),0.0);
        nlow (1) = (size_t)std::max((trunc((Xp(1)-0.5*Lp(1)-Xlow(1))/Dx)-1),0.0);
        nlow (2) = (size_t)std::max((trunc((Xp(2)-0.5*Lp(2)-Xlow(2))/Dx)-1),0.0);
        nhigh(0) = (size_t)std::min((ceil ((Xp(0)+0.5*Lp(0)-Xlow(0))/Dx)+1),(double)Ndim(0)-1.0);
        nhigh(1) = (size_t)std::min((ceil ((Xp(1)+0.5*Lp(1)-Xlow(1))/Dx)+1),(double)Ndim(1)-1.0);
        nhigh(2) = (size_t)std::min((ceil ((Xp(2)+0.5*Lp(2)-Xlow(2))/Dx)+1),(double)Ndim(2)-1.0);

        for (size_t nx=nlow(0); nx<=nhigh(0); nx++)
        for (size_t ny=nlow(1); ny<=nhigh(1); ny++)
        for (size_t nz=nlow(2); nz<=nhigh(2); nz++)
        {
            double Sf=0.0; //Shape function value
            Vec3_t Gs=OrthoSys::O; //Shape function gradient
            //Vec3_t Xn((nx+0.5)*Dx,(ny+0.5)*Dx,(nz+0.5)*Dx);
            Vec3_t Xn(nx*Dx,ny*Dx,nz*Dx);
            Xn += Xlow;
            iVec3_t Pt(nx,ny,nz);
            GIMP3D(Xp,Xn,Dx,Lp,Sf,Gs);
            if (Sf>0.0)
            {
                size_t nc = Pt2idx(Pt,Ndim);
                Particles[ip]->Nodes.Push(nc);
                Particles[ip]->Sfunc.Push(Sf);
                Particles[ip]->GSf  .Push(Gs);
                double ms = Sf*Particles[ip]->m;
                omp_set_lock(&Nodes[nc].lck);
                Nodes[nc].Mass += ms;
                Nodes[nc].Mn   += ms*Particles[ip]->v;
                //Nodes[nc].Sn    = Nodes[nc].Sn + ms*Particles[ip]->S;
                Nodes[nc].valid = true;
                omp_unset_lock(&Nodes[nc].lck);
            }
        }
        if (Particles[ip]->Nodes.Size()<2)
        {
            std::cout << "Particle does not have enough nodes" << std::endl;
            std::cout << "Time " << Time << std::endl;
            std::cout << "Ip   " << ip   << std::endl;
            std::cout << "Position  " << Particles[ip]->x << std::endl;
            std::cout << "Dx  " << Particles[ip]->Dx << std::endl;
            std::cout << "Dx0  " << Particles[ip]->Dx0 << std::endl;
            std::cout << "F  " << Particles[ip]->F << std::endl;
            std::cout << "Free? =" << Particles[ip]->IsFree() <<std::endl;
            throw new Fatal("MPM:Domain::OneStepUSF No mesh cells intersected by particles");
        }
    }

    //Producing an array of only valid nodes
    for (size_t in=0; in<Nnodes; in++)
    {
        if (Nodes[in].valid) VNodes.Push(Nodes[in].Idx);
    }

    //Fixing Dirichlet nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].UpdateNode(0.0,0.0,Mmin);
        Nodes[in].FixNode();
    }

    //Transfering info from node to particle
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        set_to_zero(Particles[ip]->L);
        for (size_t nn=0; nn<Particles[ip]->Nodes.Size(); nn++)
        {
            size_t nc = Particles[ip]->Nodes[nn];
            Vec3_t gs = Particles[ip]->GSf[nn];
            Mat3_t Gv;
            Dyad(gs,Nodes[nc].Vn,Gv);
            Particles[ip]->L = Particles[ip]->L + Gv;
        }
        
        Particles[ip]->CalcVol(Dt);
        Particles[ip]->CalcStress(Dt);
        
        Mat3_t Vsp = -Particles[ip]->V*Particles[ip]->S;
        Vec3_t Fex = Particles[ip]->b + Particles[ip]->h;
        for (size_t nn=0; nn<Particles[ip]->Nodes.Size(); nn++)
        {
            size_t nc = Particles[ip]->Nodes[nn];
            double Sf = Particles[ip]->Sfunc[nn];
            Vec3_t Gs = Particles[ip]->GSf[nn];
            Vec3_t VspGs;
            Mult(Vsp,Gs,VspGs);
            Vec3_t dF = Sf*Fex + VspGs;
            omp_set_lock(&Nodes[nc].lck);            
            Nodes[nc].Fn   += dF;
            Nodes[nc].Mn   += Dt*dF;
            omp_unset_lock(&Nodes[nc].lck);
        }
    }

    //Solve equations at the nodes    
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].UpdateNode(Gn,Dt,Mmin);
        Nodes[in].FixNode();
    }

    //Updating the positions of the MPM particles
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        for (size_t nn=0; nn<Particles[ip]->Nodes.Size(); nn++)
        {
            size_t nc = Particles[ip]->Nodes[nn];
            double sf = Particles[ip]->Sfunc[nn];
            Vec3_t gs = Particles[ip]->GSf[nn];
            Vec3_t an = Nodes[nc].Fn/Nodes[nc].Mass;
            Particles[ip]->v += sf*an*Dt;
            Particles[ip]->x += sf*Nodes[nc].Vn*Dt;
        }
        Particles[ip]->Reset(Dt);
    }
}

inline void Domain::OneStepCPI()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].Reset();
    }
    VNodes.Resize(0);
    
    //std::cout << "1" << std::endl;
    //Getting information from the corners
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ic=0;ic<Corners.Size();ic++)
    {
        iVec3_t nlow,nhigh;
        Vec3_t Xp  = *Corners[ic]->x;
        Vec3_t Fex = Corners[ic]->h;

        //nlow (0) = (size_t)std::max((trunc((Xp(0)-Xlow(0))/Dx)-1),0.0);
        //nlow (1) = (size_t)std::max((trunc((Xp(1)-Xlow(1))/Dx)-1),0.0);
        //nlow (2) = (size_t)std::max((trunc((Xp(2)-Xlow(2))/Dx)-1),0.0);
        //nhigh(0) = (size_t)std::min((ceil ((Xp(0)-Xlow(0))/Dx)+1),(double)Ndim(0)-1.0);
        //nhigh(1) = (size_t)std::min((ceil ((Xp(1)-Xlow(1))/Dx)+1),(double)Ndim(1)-1.0);
        //nhigh(2) = (size_t)std::min((ceil ((Xp(2)-Xlow(2))/Dx)+1),(double)Ndim(2)-1.0);

        nlow (0) = (size_t) trunc((Xp(0)-Xlow(0))/Dx);
        nlow (1) = (size_t) trunc((Xp(1)-Xlow(1))/Dx);
        nlow (2) = (size_t) trunc((Xp(2)-Xlow(2))/Dx);
        nhigh(0) = nlow(0) + 1;
        nhigh(1) = nlow(1) + 1;
        nhigh(2) = nlow(2) + 1;


        if (nhigh(0)>Ndim(0)||nhigh(1)>Ndim(1)||nhigh(2)>Ndim(2))
        {
            std::cout << "nhigh out of bounds" << std::endl;
            std::cout << "Ic " << ic << std::endl;
            std::cout << "Xp " << Xp << std::endl;
            std::cout << "Fex " << Fex << std::endl;
            std::cout << "nhigh " << nhigh << std::endl;
            std::cout << "Ndim " << Ndim << std::endl;
            throw new Fatal("Domain::OneStepCPI nhigh out of bounds");
        }

        size_t ncs = 0;
        for (size_t nx=nlow(0); nx<=nhigh(0); nx++)
        for (size_t ny=nlow(1); ny<=nhigh(1); ny++)
        for (size_t nz=nlow(2); nz<=nhigh(2); nz++)
        {
            //std::cout << Xp << nlow << nhigh << Fex << " " << nx << " " << ny << " " << nz << std::endl;
            Vec3_t Xn(nx*Dx,ny*Dx,nz*Dx);
            Xn += Xlow;
            iVec3_t Pt(nx,ny,nz);
            double Sf = Shape3D(Xp,Xn,Dx);
            size_t nc = Pt2idx(Pt,Ndim);
            Corners[ic]->Nodes[ncs] = nc;
            Corners[ic]->Sfunc[ncs] = Sf;
            Vec3_t dF = Sf*Fex;
            //std::cout << Xp << nlow << nhigh << " " << nx << " " << ny << " " << nz << std::endl;
            omp_set_lock(&Nodes[nc].lck);
            Nodes[nc].Mn   += Dt*dF;
            Nodes[nc].Fn   += dF;
            Nodes[nc].valid = true;
            omp_unset_lock(&Nodes[nc].lck);
            //std::cout << Xp << nlow << nhigh << " " << nx << " " << ny << " " << nz << std::endl;
            if (isnan(norm(Fex)))
            {
                std::cout << "nan value detected at corner" << std::endl;
                std::cout << "dF=" << dF  <<std::endl;
                std::cout << "sf=" << Sf <<std::endl;
                std::cout << "Time =" << Time <<std::endl;
                throw new Fatal("Domain::OneStepCPI nan value found");
            }
            ncs ++;
        }
    }

    //std::cout << "2" << std::endl;
    //Computing shape functions at particles
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        double Vp  = Particles[ip]->Vt;
        Mat3_t Vsp = -Particles[ip]->V*Particles[ip]->S;
        Vec3_t Fex = Particles[ip]->b;
        //Go through each corner
        for (size_t ic=0;ic<4;ic++)
        {
            size_t ncor = Particles[ip]->Corners[ic];
            for (size_t ics=0;ics<8;ics++)
            {
                size_t nc = Corners[ncor]->Nodes[ics];
                double sf = Corners[ncor]->Sfunc[ics];
                double ms = sf*Particles[ip]->m;
                Vec3_t VspGs;
                Vec3_t Gs   = sf*Particles[ip]->Gradvec[ic]/(6.0*Vp);
                Mult(Vsp,Gs,VspGs);
                Vec3_t dF = 0.25*sf*Fex + VspGs;
                omp_set_lock(&Nodes[nc].lck);

                Nodes[nc].Mass += 0.25*ms;
                Nodes[nc].Mn   += 0.25*ms*Particles[ip]->v + Dt*dF;
                Nodes[nc].Fn   += dF;

                //Nodes[nc].Sn    = Nodes[nc].Sn + ms*Particles[ip]->S;
                omp_unset_lock(&Nodes[nc].lck);
                if (isnan(norm(Particles[ip]->x))||isnan(norm(Particles[ip]->v))||isnan(norm(Fex)))
                {
                    std::cout << "nan value detected" << std::endl;
                    std::cout << "ip =" << ip <<std::endl;
                    std::cout << "M =" << ms               <<std::endl;
                    std::cout << "v =" << Particles[ip]->v <<std::endl;
                    std::cout << "x =" << Particles[ip]->x <<std::endl;
                    std::cout << "dF=" << dF  <<std::endl;
                    std::cout << "Mn=" << Nodes[nc].Mn <<std::endl;
                    std::cout << "sf=" << sf <<std::endl;
                    std::cout << "Time =" << Time <<std::endl;
                    throw new Fatal("Domain::OneStepCPI nan value found");
                }
            }
        }
    }

    //std::cout << "3" << std::endl;
    //Producing an array of only valid nodes
    for (size_t in=0; in<Nnodes; in++)
    {
        if (Nodes[in].valid) VNodes.Push(Nodes[in].Idx);
    }

    //std::cout << "4" << std::endl;
    //Solving the equations at the nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].UpdateNode(Gn,Dt,Mmin);
        Nodes[in].FixNode();
    }

    //std::cout << "5" << std::endl;
    //Moving corners    
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ic=0;ic<Corners.Size();ic++)
    {
        if (Corners[ic]->IsFree())
        {
            for (size_t ics=0;ics<8;ics++)
            {
                size_t nc = Corners[ic]->Nodes[ics];
                double sf = Corners[ic]->Sfunc[ics];
                Vec3_t an = Nodes[nc].Fn/Nodes[nc].Mass;
                *Corners[ic]->x += sf*Nodes[nc].Vn*Dt;
                //Corners[ic]->v  += sf*an*Dt;
            }
        }
        else *Corners[ic]->x += Corners[ic]->vf*Dt;
        Corners[ic]->Reset(Dt);
    }

    //std::cout << "6" << std::endl;
    //Transfering info from node to particle
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        double Vp  = Particles[ip]->Vt;
        set_to_zero(Particles[ip]->L);
        //Go through each corner
        for (size_t ic=0;ic<4;ic++)
        {
            size_t ncor = Particles[ip]->Corners[ic];
            for (size_t ics=0;ics<8;ics++)
            {
                size_t nc = Corners[ncor]->Nodes[ics];
                double sf = Corners[ncor]->Sfunc[ics];
                Vec3_t Gs = sf*Particles[ip]->Gradvec[ic]/(6.0*Vp);
                Vec3_t an = Nodes[nc].Fn/Nodes[nc].Mass;
                Particles[ip]->v += 0.25*sf*an*Dt;
                //Particles[ip]->x += 0.25*sf*Nodes[nc].Vn*Dt;
                if (isnan(norm(Particles[ip]->x))||isnan(norm(Particles[ip]->v)))
                {
                    std::cout << "nan value detected" << std::endl;
                    std::cout << "ip =" << ip <<std::endl;
                    std::cout << "Free? =" << Particles[ip]->IsFree() <<std::endl;
                    std::cout << "v =" << Particles[ip]->v <<std::endl;
                    std::cout << "x =" << Particles[ip]->x <<std::endl;
                    std::cout << "fh =" << Particles[ip]->h <<std::endl;
                    std::cout << "Vn =" << Nodes[nc].Vn <<std::endl;
                    std::cout << "Mn =" << Nodes[nc].Mn <<std::endl;
                    std::cout << "Ma =" << Nodes[nc].Mass <<std::endl;
                    std::cout << "an =" << an <<std::endl;
                    std::cout << "sf =" << sf <<std::endl;
                    std::cout << "Time =" << Time <<std::endl;
                    throw new Fatal("Domain::NodetoParticle nan value found");
                }
                Mat3_t Gv;
                //Dyad(Nodes[nc].Vn,Gs,Gv);
                Dyad(Gs,Nodes[nc].Vn,Gv);
                Particles[ip]->L = Particles[ip]->L + Gv;
            }
        }
        Particles[ip]->CalcVolCPI(Dt);
        Particles[ip]->CalcStress(Dt);
        Particles[ip]->Reset(Dt);
    }

    //std::cout << "6" << std::endl;
    //Resetting corners    
	//#pragma omp parallel for schedule(static) num_threads(Nproc)
    //for (size_t ic=0;ic<Corners.Size();ic++)
    //{
        //Corners[ic]->Reset(Dt);
    //}
    //std::cout << "7" << std::endl;
}                       

inline void Domain::BoundaryMesh()
{   
    Array<Array <int> > Neigh (BodyMesh->Cells.Size());
    Array<Array <int> > FNeigh(BodyMesh->Cells.Size());
    BodyMesh->FindNeigh();
    for (size_t ic;ic<BodyMesh->Cells.Size();ic++)
    {
        for (Mesh::Neighs_t::const_iterator p=BodyMesh->Cells[ic]->Neighs.begin(); p!=BodyMesh->Cells[ic]->Neighs.end(); ++p)
        {
            Neigh[ic].Push(p->second.second->ID);
            FNeigh[ic].Push(p->second.first);
        }
        if (Neigh[ic].Size()<4) 
        {
            if (!FNeigh[ic].Has(0))
            {
                BdryFac Bf;
                Bf.Ce = BodyMesh->Cells[ic];
                Bf.Facets[0] = 0;
                Bf.Facets[1] = 3;
                Bf.Facets[2] = 2;
                size_t i0 = Bf.Ce->V[Bf.Facets[0]]->ID;
                size_t i1 = Bf.Ce->V[Bf.Facets[1]]->ID;
                size_t i2 = Bf.Ce->V[Bf.Facets[2]]->ID;
                size_t io = Bf.Ce->V[1]->ID;
                Vec3_t a  = BodyMesh->Verts[i1]->C-BodyMesh->Verts[i0]->C;
                Vec3_t b  = BodyMesh->Verts[i2]->C-BodyMesh->Verts[i1]->C;
                Vec3_t c  = BodyMesh->Verts[io]->C-BodyMesh->Verts[i0]->C;
                if (dot(c,cross(a,b))>0.0)
                {
                    Bf.Facets[0] = 2;
                    Bf.Facets[2] = 0;
                }
                BFac.Push(Bf);
            } 
            if (!FNeigh[ic].Has(1))
            {
                BdryFac Bf;
                Bf.Ce = BodyMesh->Cells[ic];
                Bf.Facets[0] = 0;
                Bf.Facets[1] = 1;
                Bf.Facets[2] = 3;
                size_t i0 = Bf.Ce->V[Bf.Facets[0]]->ID;
                size_t i1 = Bf.Ce->V[Bf.Facets[1]]->ID;
                size_t i2 = Bf.Ce->V[Bf.Facets[2]]->ID;
                size_t io = Bf.Ce->V[2]->ID;
                Vec3_t a  = BodyMesh->Verts[i1]->C-BodyMesh->Verts[i0]->C;
                Vec3_t b  = BodyMesh->Verts[i2]->C-BodyMesh->Verts[i1]->C;
                Vec3_t c  = BodyMesh->Verts[io]->C-BodyMesh->Verts[i0]->C;
                if (dot(c,cross(a,b))>0.0)
                {
                    Bf.Facets[0] = 3;
                    Bf.Facets[2] = 0;
                }
                BFac.Push(Bf);
            } 
            if (!FNeigh[ic].Has(2))
            {
                BdryFac Bf;
                Bf.Ce = BodyMesh->Cells[ic];
                Bf.Facets[0] = 0;
                Bf.Facets[1] = 2;
                Bf.Facets[2] = 1;
                size_t i0 = Bf.Ce->V[Bf.Facets[0]]->ID;
                size_t i1 = Bf.Ce->V[Bf.Facets[1]]->ID;
                size_t i2 = Bf.Ce->V[Bf.Facets[2]]->ID;
                size_t io = Bf.Ce->V[3]->ID;
                Vec3_t a  = BodyMesh->Verts[i1]->C-BodyMesh->Verts[i0]->C;
                Vec3_t b  = BodyMesh->Verts[i2]->C-BodyMesh->Verts[i1]->C;
                Vec3_t c  = BodyMesh->Verts[io]->C-BodyMesh->Verts[i0]->C;
                if (dot(c,cross(a,b))>0.0)
                {
                    Bf.Facets[0] = 1;
                    Bf.Facets[2] = 0;
                }
                BFac.Push(Bf);
            } 
            if (!FNeigh[ic].Has(3))
            {
                BdryFac Bf;
                Bf.Ce = BodyMesh->Cells[ic];
                Bf.Facets[0] = 1;
                Bf.Facets[1] = 2;
                Bf.Facets[2] = 3;
                size_t i0 = Bf.Ce->V[Bf.Facets[0]]->ID;
                size_t i1 = Bf.Ce->V[Bf.Facets[1]]->ID;
                size_t i2 = Bf.Ce->V[Bf.Facets[2]]->ID;
                size_t io = Bf.Ce->V[0]->ID;
                Vec3_t a  = BodyMesh->Verts[i1]->C-BodyMesh->Verts[i0]->C;
                Vec3_t b  = BodyMesh->Verts[i2]->C-BodyMesh->Verts[i1]->C;
                Vec3_t c  = BodyMesh->Verts[io]->C-BodyMesh->Verts[i0]->C;
                if (dot(c,cross(a,b))>0.0)
                {
                    Bf.Facets[0] = 3;
                    Bf.Facets[2] = 1;
                }
                BFac.Push(Bf);
            } 
        }
    }
    BoundaryMeshUpdate();
}

inline void Domain::BoundaryMeshUpdate()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t is=0;is<BFac.Size();is++)
    {
        size_t i0 = BFac[is].Ce->V[BFac[is].Facets[0]]->ID;
        size_t i1 = BFac[is].Ce->V[BFac[is].Facets[1]]->ID;
        size_t i2 = BFac[is].Ce->V[BFac[is].Facets[2]]->ID;
        Vec3_t x0 = BodyMesh->Verts[i0]->C;
        Vec3_t x1 = BodyMesh->Verts[i1]->C;
        Vec3_t x2 = BodyMesh->Verts[i2]->C;
        Vec3_t a  = x1 - x0;
        Vec3_t b  = x2 - x1;
        BFac[is].Nor  = cross(a,b);
        BFac[is].Nor /= norm(BFac[is].Nor);
        BFac[is].Xmin(0) = std::min(x0(0),std::min(x1(0),x2(0)));
        BFac[is].Xmin(1) = std::min(x0(1),std::min(x1(1),x2(1)));
        BFac[is].Xmin(2) = std::min(x0(2),std::min(x1(2),x2(2)));
        BFac[is].Xmax(0) = std::max(x0(0),std::max(x1(0),x2(0)));
        BFac[is].Xmax(1) = std::max(x0(1),std::max(x1(1),x2(1)));
        BFac[is].Xmax(2) = std::max(x0(2),std::max(x1(2),x2(2)));
        
        //std::cout << is << " " << BodyMesh->Verts[i0]->C << " " << BFac[is].Nor << " " << BFac[is].Xmin << " " << BFac[is].Xmax << std::endl;
        //std::cout << is << " " << BodyMesh->Verts[i1]->C << " " << BFac[is].Nor << " " << BFac[is].Xmin << " " << BFac[is].Xmax << std::endl;
        //std::cout << is << " " << BodyMesh->Verts[i2]->C << " " << BFac[is].Nor << " " << BFac[is].Xmin << " " << BFac[is].Xmax << std::endl;
    }
}

inline void Domain::UpdateMesh()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ic=0;ic<BodyMesh->Cells.Size();ic++)
    {
        Vec3_t v0 = BodyMesh->Verts[BodyMesh->Cells[ic]->V[0]->ID]->C;
        Vec3_t v1 = BodyMesh->Verts[BodyMesh->Cells[ic]->V[1]->ID]->C;
        Vec3_t v2 = BodyMesh->Verts[BodyMesh->Cells[ic]->V[2]->ID]->C;
        Vec3_t v3 = BodyMesh->Verts[BodyMesh->Cells[ic]->V[3]->ID]->C;
        BodyMesh->Cells[ic]->Xmin(0) = std::min(std::min(v0(0),v1(0)),std::min(v2(0),v3(0)));
        BodyMesh->Cells[ic]->Xmin(1) = std::min(std::min(v0(1),v1(1)),std::min(v2(1),v3(1)));
        BodyMesh->Cells[ic]->Xmin(2) = std::min(std::min(v0(2),v1(2)),std::min(v2(2),v3(2)));
        BodyMesh->Cells[ic]->Xmax(0) = std::max(std::max(v0(0),v1(0)),std::max(v2(0),v3(0)));
        BodyMesh->Cells[ic]->Xmax(1) = std::max(std::max(v0(1),v1(1)),std::max(v2(1),v3(1)));
        BodyMesh->Cells[ic]->Xmax(2) = std::max(std::max(v0(2),v1(2)),std::max(v2(2),v3(2)));
    }

    BoundaryMeshUpdate();
}

inline void Domain::PosMesh(Vec3_t & pos)
{
    if (BodyMesh==NULL) throw new Fatal("MPM::PosMesh can only be used if a body mesh is defined");
    Vec3_t Xmin,Xmax;
    BoundingBox(Xmin,Xmax);
    Vec3_t trans = pos - 0.5*(Xmax+Xmin);
    for (size_t ic=0;ic<Corners.Size();ic++)
    {
        *Corners[ic]->x += trans;
        Corners[ic]->xb += trans;
    }

    for (size_t ip=0;ip<Particles.Size();ip++)
    {
        Particles[ip]->x  += trans;
        Particles[ip]->x0 += trans;
    }

    UpdateMesh();
}

inline void Domain::Solve(double Tf, double dt, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{

    FileKey.Printf("%s",TheFileKey);
    Nproc = TheNproc;

    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1    , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                    , TERM_RST);
    
    idx_out     = 0;
    double tout = Time;
    Dt          = dt;
    //Initializing particles
    Mmin = Particles[0]->m;
    for (size_t ip=0;ip<Particles.Size();ip++)
    {
        Particles[ip]->xb = Particles[ip]->x - Particles[ip]->v*Dt;
        if (Mmin > Particles[ip]->m) Mmin = Particles[ip]->m;
    }
    Mmin *= 1.0e-12;
    if (BodyMesh!=NULL)
    {
        BoundaryMesh();
    }
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
        if (BodyMesh!=NULL)
        {
            OneStepCPI();
            UpdateMesh();
        }
        else
        {
            ParticleToNode();
            NodeToParticle();
        //OneStepUSF();
        }
        Time += dt;
        //std::cout << Time << std::endl;
    }
}
}
#endif //MECHSYS_MPM_DOMAIN_H
