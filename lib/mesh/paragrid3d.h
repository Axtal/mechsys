/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

#ifndef MECHSYS_PARAGRID3D_H
#define MECHSYS_PARAGRID3D_H

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/mesh/structured.h>

// enum
enum CellType { Inner_t=-1, BryIn_t=-2, BryOut_t=-3, Outer_t=-4 };

class ParaGrid3D
{
public:
    // enums
    enum MetisAlgo { Recursive_t=0, Kway_t=1, VKway_t=2 };

    // Constructor & Destructor
    ParaGrid3D (Array<int> const & N, Array<double> & L, char const * FileKey=NULL, MetisAlgo Algo=Recursive_t);

    // Methods
    int  FindCell  (Vec3_t const & X) const;                       ///< Get the Id of the cell containing the point with coordinates equal to X
    void FindCells (Vec3_t const & X, double R, int Ids[8]) const; ///< Get the Ids of the cells crossed by the cube centered at X

    // Data (read-only)
    Array<int>          N;             ///< Number of cells along each axis
    Array<double>       L;             ///< Size of the grid (xmin,xmax, ymin,ymax, zmin,zmax)
    Vec3_t              D;             ///< Dimension of each cell (dx,dy,dz)
    Array<int>          NeighProcs;    ///< Neighbour processors/partitions of this processor
    Array<CellType>     Id2Type;       ///< Id of cell to type: size==num cells
    Array<int>          Id2Proc;       ///< Id of cell to partition/processor number: size==num cells
    Array< Array<int> > Id2NeighProcs; ///< Id of cell to neighbour processors/partitions of cell: size == num cells. 
                                       ///<   Only valid for those cells belonging to this processor.
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline ParaGrid3D::ParaGrid3D (Array<int> const & TheN, Array<double> & TheL, char const * FileKey, MetisAlgo Algo)
    : N(TheN), L(TheL)
{
    // cell size
    D = (L[1]-L[0])/static_cast<double>(N[0]),
        (L[3]-L[2])/static_cast<double>(N[1]),
        (L[5]-L[4])/static_cast<double>(N[2]);

    // proc id and number of cells
    int my_id  = MPI::COMM_WORLD.Get_rank();
    int nprocs = MPI::COMM_WORLD.Get_size();
    int ncells = N[0]*N[1]*N[2];

    // find neighbours for METIS
    Array<int> Xadj(0, /*justone*/true);
    Array<int> Adjn;
    for (int k=0; k<N[2]; ++k)
    for (int j=0; j<N[1]; ++j)
    for (int i=0; i<N[0]; ++i)
    {
        int back   = (i-1) +  j   *N[0] +  k   *N[0]*N[1];
        int front  = (i+1) +  j   *N[0] +  k   *N[0]*N[1];
        int left   =  i    + (j-1)*N[0] +  k   *N[0]*N[1];
        int right  =  i    + (j+1)*N[0] +  k   *N[0]*N[1];
        int bottom =  i    +  j   *N[0] + (k-1)*N[0]*N[1];
        int top    =  i    +  j   *N[0] + (k+1)*N[0]*N[1];
        Array<int> neighs; // not including the opposite-diagonal ones (for better quality domain decomposition)
        if (i>0)      neighs.Push (back);
        if (i<N[0]-1) neighs.Push (front);
        if (j>0)      neighs.Push (left);
        if (j<N[1]-1) neighs.Push (right);
        if (k>0)      neighs.Push (bottom);
        if (k<N[2]-1) neighs.Push (top);
        Xadj.Push (Xadj.Last()+neighs.Size());
        for (size_t p=0; p<neighs.Size(); ++p) Adjn.Push (neighs[p]);
    }

    // domain decomposition
    Id2Proc.Resize    (ncells);
    Id2Proc.SetValues (0);
#if defined(HAS_PARMETIS) && defined(HAS_MPI)
    if (nprocs>1)
    {
        int wgtflag    = 0; // No weights
        int numflag    = 0; // zero numbering
        int options[5] = {0,0,0,0,0};
        if (Algo==Recursive_t)
        {
            int edgecut;
            METIS_PartGraphRecursive (&ncells, Xadj.GetPtr(), Adjn.GetPtr(), NULL, NULL, &wgtflag, &numflag, &nprocs, options, &edgecut, Id2Proc.GetPtr());
            if (my_id==0) printf("METIS: Recursive: Edgecut = %d\n",edgecut);
        }
        else if (Algo==Kway_t)
        {
            int edgecut;
            METIS_PartGraphKway (&ncells, Xadj.GetPtr(), Adjn.GetPtr(), NULL, NULL, &wgtflag, &numflag, &nprocs, options, &edgecut, Id2Proc.GetPtr());
            if (my_id==0) printf("METIS: Kway: Edgecut = %d\n",edgecut);
        }
        else if (Algo==VKway_t)
        {
            int volume;
            METIS_PartGraphVKway (&ncells, Xadj.GetPtr(), Adjn.GetPtr(), NULL, NULL, &wgtflag, &numflag, &nprocs, options, &volume, Id2Proc.GetPtr());
            if (my_id==0) printf("METIS: VKway: Volume = %d\n",volume);
        }
        else throw new Fatal("ParaGrid3D::ParaGrid3D: MetisAlgo==%d is invalid",Algo);
#else
        throw new Fatal("ParaGrid3D::ParaGrid3D: This method requires ParMETIS and MPI (if Part array is not provided)");
#endif
    }

    // find types and neighbour processors of all cells
    Id2Type      .Resize (ncells);
    Id2NeighProcs.Resize (ncells);
    Array<int> bryouts; // boundary cells outside this processor
    for (int k=0; k<N[2]; ++k)
    for (int j=0; j<N[1]; ++j)
    for (int i=0; i<N[0]; ++i)
    {
        int n = i + j*N[0] + k*N[0]*N[1];
        Id2Type[n] = Outer_t;
        if (Id2Proc[n]==my_id) // cell is in this processor
        {
            Id2Type[n] = Inner_t;
            for (int r=i-1; r<i+2; ++r)
            for (int s=j-1; s<j+2; ++s)
            for (int t=k-1; t<k+2; ++t) // look for all 26 neighbours of cell # n
            {
                if (!(r==i && s==j && t==k))
                {
                    int m = r + s*N[0] + t*N[0]*N[1];
                    if (r>=0 && r<N[0] &&
                        s>=0 && s<N[1] &&
                        t>=0 && t<N[2]) // cell # m is inside limits
                    {
                        if (Id2Proc[m]!=my_id) // cell # m is in another processor
                        {
                            Id2Type[n] = BryIn_t;
                            bryouts.XPush (m);
                            Id2NeighProcs[n].XPush (Id2Proc[m]);
                        }
                    }
                }
            }
        }
    }
    for (size_t i=0; i<bryouts.Size(); ++i)
    {
        int     m  = bryouts[i];
        Id2Type[m] = BryOut_t;
        NeighProcs.XPush (Id2Proc[m]);
    }

    // write VTU
    if (FileKey!=NULL)
    {
        std::ostringstream oss;
        size_t nn = (N[0]+1)*(N[1]+1)*(N[2]+1); // number of Nodes
        size_t nc = ncells;                     // number of Cells

        // constants
        size_t          nimax = 40;        // number of integers in a line
        size_t          nfmax =  6;        // number of floats in a line
        Util::NumStream nsflo = Util::_8s; // number format for floats

        // header
        oss << "<?xml version=\"1.0\"?>\n";
        oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        oss << "  <UnstructuredGrid>\n";
        oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << nc << "\">\n";

        // nodes: coordinates
        oss << "      <Points>\n";
        oss << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        size_t k = 0; oss << "        ";
        int nx = N[0]+1;
        int ny = N[1]+1;
        //int nz = N[2]+1;
        for (size_t n=0; n<nn; ++n)
        {
            int K =  n / (nx*ny);
            int J = (n % (nx*ny)) / nx;
            int I = (n % (nx*ny)) % nx;
            oss << "  " << nsflo << L[0] + I*D[0] << " ";
            oss <<         nsflo << L[2] + J*D[1] << " ";
            oss <<         nsflo << L[4] + K*D[2];
            k++;
            VTU_NEWLINE (n,k,nn,nfmax/3-1,oss);
        }
        oss << "        </DataArray>\n";
        oss << "      </Points>\n";

        // elements: connectivity, offsets, types
        oss << "      <Cells>\n";
        oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        nx = N[0];
        ny = N[1];
        //nz = N[2];
        for (size_t n=0; n<nc; ++n)
        {
            int K =  n / (nx*ny);
            int J = (n % (nx*ny)) / nx;
            int I = (n % (nx*ny)) % nx;
            int C[8][3] = {{I  ,  J  ,  K},
                           {I+1,  J  ,  K},
                           {I+1,  J+1,  K},
                           {I  ,  J+1,  K},
                           {I  ,  J  ,  K+1},
                           {I+1,  J  ,  K+1},
                           {I+1,  J+1,  K+1},
                           {I  ,  J+1,  K+1}};
            oss << "  ";
            for (size_t j=0; j<8; ++j) oss << C[j][0] + C[j][1]*(N[0]+1) + C[j][2]*(N[0]+1)*(N[1]+1) << " ";
            k++;
            VTU_NEWLINE (n,k,nc,nimax/8,oss);
        }
        oss << "        </DataArray>\n";
        oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        size_t offset = 0;
        for (size_t i=0; i<nc; ++i)
        {
            offset += 8;
            oss << (k==0?"  ":" ") << offset;
            k++;
            VTU_NEWLINE (i,k,nc,nimax,oss);
        }
        oss << "        </DataArray>\n";
        oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<nc; ++i)
        {
            oss << (k==0?"  ":" ") << 12; // VTK_HEXAHEDRON
            k++;
            VTU_NEWLINE (i,k,nc,nimax,oss);
        }
        oss << "        </DataArray>\n";
        oss << "      </Cells>\n";

        // data -- elements
        oss << "      <CellData Scalars=\"TheScalars\">\n";
        oss << "        <DataArray type=\"Int32\" Name=\"" << "celltype" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<nc; ++i)
        {
            oss << (k==0?"  ":" ") << Id2Type[i];
            k++;
            VTU_NEWLINE (i,k,nc,nimax,oss);
        }
        oss << "        </DataArray>\n";
        oss << "        <DataArray type=\"Int32\" Name=\"" << "proc" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        int max_neigh_procs = 0;
        for (size_t i=0; i<nc; ++i)
        {
            if (Id2NeighProcs[i].Size()>(size_t)max_neigh_procs) max_neigh_procs = Id2NeighProcs[i].Size();
            oss << (k==0?"  ":" ") << Id2Proc[i];
            k++;
            VTU_NEWLINE (i,k,nc,nimax,oss);
        }
        oss << "        </DataArray>\n";
        if (max_neigh_procs>0)
        {
            oss << "        <DataArray type=\"Int32\" Name=\"" << "neighproc" << "\" NumberOfComponents=\""<<max_neigh_procs<<"\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t i=0; i<nc; ++i)
            {
                oss << (k==0?"  ":" ");
                for (int j=0; j<max_neigh_procs; ++j)
                {
                    if (j>static_cast<int>(Id2NeighProcs[i].Size()-1)) oss << -1 << " ";
                    else oss << Id2NeighProcs[i][j] << " ";
                }
                k++;
                VTU_NEWLINE (i,k,nc,nimax,oss);
            }
            oss << "        </DataArray>\n";
        }
        oss << "      </CellData>\n";

        // Bottom
        oss << "    </Piece>\n";
        oss << "  </UnstructuredGrid>\n";
        oss << "</VTKFile>" << std::endl;

        // Write to file
        String fn(FileKey); fn.append(".vtu");
        std::ofstream of(fn.CStr(), std::ios::out);
        of << oss.str();
        of.close();
    }
}

inline int ParaGrid3D::FindCell (Vec3_t const & X) const
{
    // cell
    int I = static_cast<int>((X(0)-L[0])/D(0));
    int J = static_cast<int>((X(1)-L[2])/D(1));
    int K = static_cast<int>((X(2)-L[4])/D(2));
    int n = I + J*N[0] + K*N[0]*N[1];

    // check
    if (n<0) throw new Fatal("ParaGrid3D::FindCellId:: Point X=(%g,%g,%g) is outside grid: N=(%d,%d,%d), L=[x:(%g,%g) y:(%g,%g) z:(%g,%g)]", X(0),X(1),X(2), N[0],N[1],N[2], L[0],L[1],L[2],L[3],L[4],L[5]);
    return n;
}

inline void ParaGrid3D::FindCells (Vec3_t const & X, double R, int Ids[8]) const
{
    // coordinates of corners
    double xa=X(0)-R;   double xb=X(0)+R;
    double ya=X(1)-R;   double yb=X(1)+R;
    double za=X(2)-R;   double zb=X(2)+R;

    // check
    if (xa<L[0] || xb>L[1] ||
        ya<L[2] || yb>L[3] ||
        za<L[4] || zb>L[5]) throw new Fatal("ParaGrid3D::FindCellsIds: Corner of cube centered at (%g,%g,%g) with inner radius=%g is outside limits of domain (N=(%d,%d,%d), L=[x:(%g,%g) y:(%g,%g) z:(%g,%g)])", X(0),X(1),X(2),R, N[0],N[1],N[2], L[0],L[1],L[2],L[3],L[4],L[5]);

    // corners
    double C[8][3] = {{xa,  ya,  za},
                      {xb,  ya,  za},
                      {xa,  yb,  za},
                      {xb,  yb,  za},
                      {xa,  ya,  zb},
                      {xb,  ya,  zb},
                      {xa,  yb,  zb},
                      {xb,  yb,  zb}};

    // find cells ids
    for (size_t i=0; i<8; ++i)
    {
        // cell touched by corner
        int I = static_cast<int>((C[i][0]-L[0])/D(0));
        int J = static_cast<int>((C[i][1]-L[2])/D(1));
        int K = static_cast<int>((C[i][2]-L[4])/D(2));
        Ids[i] = I + J*N[0] + K*N[0]*N[1];

        // check
        if (Ids[i]<0) throw new Fatal("ParaGrid3D::FindCellsIds:: __internal_error: cell id (%d) cannot be negative",Ids[i]);
    }
}

#endif // MECHSYS_PARAGRID3D_H
