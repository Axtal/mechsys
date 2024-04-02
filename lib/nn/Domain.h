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


#ifndef MECHSYS_NN_DOMAIN_H
#define MECHSYS_NN_DOMAIN_H

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
#include<random>

// Mechsys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/numstreams.h>

namespace NN
{

inline size_t Pt2idx(iVec3_t iv, iVec3_t & Dim) // Calculates the index of the cell at coordinates iv for a cubic lattice of dimensions Dim
{
    return iv(0) + iv(1)*Dim(0) + iv(2)*Dim(0)*Dim(1);
}

inline void   idx2Pt(size_t n, iVec3_t & iv, iVec3_t & Dim) // Calculates the coordinates from the index
{
    iv(0) = n%Dim(0);
    iv(1) = (n/Dim(0))%(Dim(1));
    iv(2) = n/(Dim(0)*Dim(1));
}

inline double sigmoid (double x)
{
    return 1.0/(1.0 + exp(-x));
}

inline double dsigmoid(double x)
{
    double s = sigmoid(x);
    return s*(1.0-s);
}

inline size_t Training(char const * FileKey, Array<Array <double> > & I, Array<Array <double> > & O)
{
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    size_t val[1];
    H5LTread_dataset(file_id,"/NI",H5T_NATIVE_ULONG,val);
    size_t NI = val[0];

    H5LTread_dataset(file_id,"/NO",H5T_NATIVE_ULONG,val);
    size_t NO = val[0];

    H5LTread_dataset(file_id,"/NT",H5T_NATIVE_ULONG,val);
    size_t NT = val[0];

    I.Resize(NT);
    O.Resize(NT);

    double Iu[NT*NI];
    H5LTread_dataset_double(file_id,"/I",Iu);
    size_t nu=0;
    for (size_t nt=0;nt<NT;nt++)
    {
        I[nt].Resize(NI);
        for (size_t ni=0;ni<NI;ni++)
        {
            I[nt][ni] = Iu[nu];
            //std::cout << Iu[nu] << " ";
            nu++;
        }
        //std::cout << std::endl;
    }

    double Ou[NT*NI];
    H5LTread_dataset_double(file_id,"/O",Ou);
    nu=0;
    for (size_t nt=0;nt<NT;nt++)
    {
        O[nt].Resize(NO);
        for (size_t no=0;no<NO;no++)
        {
            O[nt][no] = Ou[nu];
            //std::cout << Ou[nu] << " ";
            nu++;
        }
        //std::cout << std::endl;
    }

    
    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);

    return NT;
}

class Domain
{
public:

    //Constructors
    Domain(char const * FileKey);        //Domain loading prexisting NN
    Domain(Array<Array <double> > & In,  //number of inputs
           Array<Array <double> > & Out, //number of outputs 
           size_t Nhidden,               //number of hidden layers
           size_t *Nn                    //array with the number of neurons per layer
           );
    Domain(Array<Array <double> > & In,  //number of inputs
           Array<Array <double> > & Out, //number of outputs 
           size_t Nn                     //number of neurons for the single layer
           );
    Domain(Array<Array <double> > & In,  //number of inputs
           Array<Array <double> > & Out, //number of outputs 
           size_t Nh,                    //Number of hidden layers
           size_t Nn                     //number of neurons per hidden layer
           );

    //Data
    Array<Array<Array <double> > >  W;                       ///< Array  of Weights
    Array<Array <double> >          B;                       ///< Array  of Biases
    //Array<Array <double> >          A;                       ///< Array  of Amplitudes
    Array<Array <double> >          I;                       ///< Array  of Inputs
    Array<Array <double> >          O;                       ///< Array  of Outputs
    Array<Array<Array <double> > > OT;                       ///< Array of outputs for training
    Array<Array<Array <double> > >  E;                       ///< Array  of Errors
    Array <size_t>                  N;                       ///< Number of neurons per hidden layer
    double                      Alpha;                       ///< Learning rate
    double                     Gerror;                       ///< Global error
    size_t                         NI;                       ///< Number of Inputs
    size_t                         NO;                       ///< Number of Outputs
    size_t                         NH;                       ///< Number of hidden layers
    size_t                         NT;                       ///< Number of training sets
    //size_t                      Nproc;                       ///< Number of computing cores


    //Methods
    void Initialize (size_t seed);                          ///< Initialize the weights
    void Forward (size_t Nproc);                            ///< Propagate outputs from training
    void BackErr (size_t Nproc);                            ///< Back propagate errors from training
    void Train   (size_t epochs, size_t Nproc);             ///< Train the NN
    void Predict (Array<double> & Input, Array<double> & Predict, size_t Nproc);   ///<Predict outcome with trained NN
    void TrainDat(char const * FileKey);                    ///< Load training set
    

    //Writing Methods
    void Save         (char const * FileKey);   ///< Save the NN
    void Load         (char const * FileKey);   ///< Load the NN

    omp_lock_t                   lck;                         ///< to protect variables in multithreading
};

inline Domain::Domain(char const * FileKey)
{
    omp_init_lock(&lck);
    Alpha = 0.3;
    Load(FileKey);
}

inline Domain::Domain(Array<Array <double> > & In, Array<Array <double> > & Out, size_t Nhidden, size_t *Nn)
{
    omp_init_lock(&lck);
    Alpha = 0.3;
    NT = In.Size();
    OT.Resize(NT);
    E .Resize(NT);
    NI = In [0].Size();
    NO = Out[0].Size();
    NH = Nhidden;
    I = In;
    O = Out;
    W.Resize(Nhidden+1);
    B.Resize(Nhidden+1);
    //A.Resize(Nhidden+1);
    N.Resize(Nhidden+1);
    for (size_t nh=0;nh<=Nhidden;nh++)
    {
        if (Nhidden==0)
        {
            W[nh].Resize(NI);
            for (size_t nl=0;nl<NI;nl++)
            {
                W[nh][nl].Resize(NO);
            }
            B[nh].Resize(NO);
            //A[nh].Resize(NO);
        }
        else if (nh==0)
        {
            W[nh].Resize(NI);
            for (size_t nl=0;nl<NI;nl++)
            {
                W[nh][nl].Resize(Nn[0]);
            }
            N[nh] = Nn[nh];
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
        else if (nh==Nhidden)
        {
            W[nh].Resize(Nn[Nhidden-1]);
            for (size_t nl=0;nl<Nn[Nhidden-1];nl++)
            {
                W[nh][nl].Resize(NO);
            }
            N[nh] = NO;
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
        else
        {
            W[nh].Resize(Nn[nh-1]);
            for (size_t nl=0;nl<Nn[nh-1];nl++)
            {
                W[nh][nl].Resize(Nn[nh]);
            }
            N[nh] = Nn[nh];
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
    }
    for (size_t nt=0;nt<NT;nt++)
    {
        OT[nt].Resize(Nhidden+1);
        E [nt].Resize(Nhidden+1);
        for (size_t nh=0;nh<OT[nt].Size();nh++)
        {
            OT[nt][nh].Resize(N[nh]);
            E [nt][nh].Resize(N[nh]);
        }
    }
}

inline Domain::Domain(Array<Array <double> > & In, Array<Array <double> > & Out, size_t NN)
{
    omp_init_lock(&lck);
    Alpha = 0.3;
    size_t Nhidden = 1;
    size_t Nn[1];
    Nn[0] = NN;
    NT = In.Size();
    OT.Resize(NT);
    E .Resize(NT);
    NI = In [0].Size();
    NO = Out[0].Size();
    NH    = Nhidden;
    I = In;
    O = Out;
    W.Resize(Nhidden+1);
    B.Resize(Nhidden+1);
    //A.Resize(Nhidden+1);
    N.Resize(Nhidden+1);
    for (size_t nh=0;nh<=Nhidden;nh++)
    {
        if (Nhidden==0)
        {
            W[nh].Resize(NI);
            for (size_t nl=0;nl<NI;nl++)
            {
                W[nh][nl].Resize(NO);
            }
            B[nh].Resize(NO);
            //A[nh].Resize(NO);
        }
        else if (nh==0)
        {
            W[nh].Resize(NI);
            for (size_t nl=0;nl<NI;nl++)
            {
                W[nh][nl].Resize(Nn[0]);
            }
            N[nh] = Nn[nh];
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
        else if (nh==Nhidden)
        {
            W[nh].Resize(Nn[Nhidden-1]);
            for (size_t nl=0;nl<Nn[Nhidden-1];nl++)
            {
                W[nh][nl].Resize(NO);
            }
            N[nh] = NO;
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
        else
        {
            W[nh].Resize(Nn[nh-1]);
            for (size_t nl=0;nl<Nn[nh-1];nl++)
            {
                W[nh][nl].Resize(Nn[nh]);
            }
            N[nh] = Nn[nh];
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
    }
    for (size_t nt=0;nt<NT;nt++)
    {
        OT[nt].Resize(Nhidden+1);
        E [nt].Resize(Nhidden+1);
        for (size_t nh=0;nh<OT[nt].Size();nh++)
        {
            OT[nt][nh].Resize(N[nh]);
            E [nt][nh].Resize(N[nh]);
        }
    }
}

inline Domain::Domain(Array<Array <double> > & In, Array<Array <double> > & Out, size_t Nhidden, size_t NN)
{
    omp_init_lock(&lck);
    Alpha = 0.3;
    size_t Nn[Nhidden];
    for (size_t nh=0;nh<Nhidden;nh++) Nn[nh] = NN;
    NT = In.Size();
    OT.Resize(NT);
    E .Resize(NT);
    NI = In [0].Size();
    NO = Out[0].Size();
    NH    = Nhidden;
    I = In;
    O = Out;
    W.Resize(Nhidden+1);
    B.Resize(Nhidden+1);
    //A.Resize(Nhidden+1);
    N.Resize(Nhidden+1);
    for (size_t nh=0;nh<=Nhidden;nh++)
    {
        if (Nhidden==0)
        {
            W[nh].Resize(NI);
            for (size_t nl=0;nl<NI;nl++)
            {
                W[nh][nl].Resize(NO);
            }
            B[nh].Resize(NO);
            //A[nh].Resize(NO);
        }
        else if (nh==0)
        {
            W[nh].Resize(NI);
            for (size_t nl=0;nl<NI;nl++)
            {
                W[nh][nl].Resize(Nn[0]);
            }
            N[nh] = Nn[nh];
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
        else if (nh==Nhidden)
        {
            W[nh].Resize(Nn[Nhidden-1]);
            for (size_t nl=0;nl<Nn[Nhidden-1];nl++)
            {
                W[nh][nl].Resize(NO);
            }
            N[nh] = NO;
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
        else
        {
            W[nh].Resize(Nn[nh-1]);
            for (size_t nl=0;nl<Nn[nh-1];nl++)
            {
                W[nh][nl].Resize(Nn[nh]);
            }
            N[nh] = Nn[nh];
            B[nh].Resize(N[nh]);
            //A[nh].Resize(N[nh]);
        }
    }
    for (size_t nt=0;nt<NT;nt++)
    {
        OT[nt].Resize(Nhidden+1);
        E [nt].Resize(Nhidden+1);
        for (size_t nh=0;nh<OT[nt].Size();nh++)
        {
            OT[nt][nh].Resize(N[nh]);
            E [nt][nh].Resize(N[nh]);
        }
    }
}

inline void Domain::Initialize(size_t seed = 0)
{
    //srand(seed);
    for (size_t nh=0;nh<W.Size();nh++)
    {
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0.0,1.0/sqrt(W[nh].Size()));    
        for (size_t ni=0;ni<W[nh].Size();ni++)
        {
            for (size_t no=0;no<W[nh][ni].Size();no++)
            {
                //double r = double(rand())/RAND_MAX;
                //W[nh][ni][no] = (2.0*r-1.0)/sqrt(W[nh].Size());
                W[nh][ni][no] = distribution(generator);
            }
        }
        for (size_t nn=0;nn<N[nh];nn++)
        {
            B[nh][nn] = 0.0;
            //A[nh][nn] = 1.0;
        }
    }
}

inline void Domain::BackErr (size_t Nproc)
{
    //Back propagating errors
    Gerror = 0.0;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nt=0;nt<NT;nt++)
    {
        for (int nh=W.Size()-1;nh>=0;nh--)
        {
            if (nh==W.Size()-1)
            {
                for (size_t nn=0;nn<N[nh];nn++)
                {
                    //E [nt][nh][nn] = (O[nt][nn]-OT[nt][nh][nn])*OT[nt][nh][nn]*(1.0-OT[nt][nh][nn])/NT; //Square cost
                    E [nt][nh][nn] = (O[nt][nn]-OT[nt][nh][nn])/NT; //Cross entropy cost
                    omp_set_lock(&lck);
                    Gerror += (O[nt][nn]-OT[nt][nh][nn])*(O[nt][nn]-OT[nt][nh][nn]);
                    omp_unset_lock(&lck);
                }
            }
            else
            {
                ParTMult(W[nh+1],E [nt][nh+1],E [nt][nh],1);
                for (size_t nn=0;nn<N[nh];nn++)
                {
                    E[nt][nh][nn] *= OT[nt][nh][nn]*(1.0-OT[nt][nh][nn]);
                }
            }
        }
    }
    //Calculating deltas and updating weights and biases
    for (size_t nh=0;nh<=NH;nh++)
    {
        if (nh>0)
        {
            for (size_t nl=0;nl<W[nh].Size();nl++)
            {
                #pragma omp parallel for schedule(static) num_threads(Nproc)
                for (size_t nn=0;nn<W[nh][nl].Size();nn++)
                {
                    double delta = 0.0;
                    for (size_t nt=0;nt<NT;nt++)
                    {
                        delta += E[nt][nh][nn]*OT[nt][nh-1][nl];
                        //delta += Alpha*E[nt][nh][nn]*OT[nt][nh][nn]*(1.0-OT[nt][nh][nn])*OT[nt][nh-1][nl];
                        //delta += Alpha*E[nt][nh][nn]*OT[nt][nh][nn]*(1.0-OT[nt][nh][nn]/A[nh][nn])*OT[nt][nh-1][nl];
                    }
                    W[nh][nl][nn] += Alpha*delta;
                }
            }
        }
        else
        {
            for (size_t nl=0;nl<W[nh].Size();nl++)
            {
                #pragma omp parallel for schedule(static) num_threads(Nproc)
                for (size_t nn=0;nn<W[nh][nl].Size();nn++)
                {
                    double delta = 0.0;
                    for (size_t nt=0;nt<NT;nt++)
                    {
                        delta += E[nt][nh][nn]*I[nt][nl];
                        //delta += Alpha*E[nt][nh][nn]*OT[nt][nh][nn]*(1.0-OT[nt][nh][nn]/A[nh][nn])*I[nt][nl];
                    }
                    W[nh][nl][nn] += Alpha*delta;
                    //std::cout << W[nh][nl][nn] << " ";
                }
            }
            //std::cout << std::endl;
        }
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t nn=0;nn<B[nh].Size();nn++)
        {
            //double deltaa = 0.0;
            double deltab = 0.0;
            for (size_t nt=0;nt<NT;nt++)
            {
                //deltaa += Alpha*E[nt][nh][nn]*OT[nt][nh][nn]/A[nh][nn];
                //deltab += Alpha*E[nt][nh][nn]*OT[nt][nh][nn]*(1.0-OT[nt][nh][nn]/A[nh][nn]);
                deltab += E[nt][nh][nn];
            }
            B[nh][nn] += Alpha*deltab;
            //A[nh][nn] += deltaa;
            //std::cout << A[nh][nn] << " ";
        }
        //std::cout << std::endl;
    }
}

inline void Domain::Forward (size_t Nproc)
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nt=0;nt<NT;nt++)
    {
        for (size_t nh=0;nh<=NH;nh++)
        {
            if (nh==0)
            {
                ParMult(W[nh],I[nt],OT[nt][nh],1);             
            }
            else
            {
                ParMult(W[nh],OT[nt][nh-1],OT[nt][nh],1);             
            }
            for (size_t nn=0;nn<N[nh];nn++)
            {
                //OT[nt][nh][nn] = A[nh][nn]*sigmoid(OT[nt][nh][nn]+B[nh][nn]);
                OT[nt][nh][nn] = sigmoid(OT[nt][nh][nn]+B[nh][nn]);
            }
        }
    }
}

inline void Domain::Train   (size_t epochs, size_t Nproc)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Training the NN ---------------------------------------------------------------------%s\n", TERM_CLR1, TERM_RST);
    printf("%s      Number of inputs            =  %zd%s\n", TERM_CLR4, NI   , TERM_RST);
    printf("%s      Number of outputs           =  %zd%s\n", TERM_CLR4, NO   , TERM_RST);
    printf("%s      Number of hidden layers     =  %zd%s\n", TERM_CLR4, NH   , TERM_RST);
    printf("%s      Number of neurons per layer =  %zd%s\n", TERM_CLR4, N[0] , TERM_RST);


    for (size_t ne=0;ne<epochs;ne++)
    {
        //std::cout << "1" << std::endl;
        Forward(Nproc);
        //std::cout << "2" << std::endl;
        BackErr(Nproc);
        //std::cout << Gerror << std::endl;
    }
    printf("%s      Final global error          =  %g%s\n", TERM_CLR4, Gerror, TERM_RST);
}

inline void Domain::Predict (Array<double> & Input, Array<double> & Predict, size_t Nproc)
{
    for (size_t nh=0;nh<=NH;nh++)
    {
        if (nh==0)
        {
            ParMult(W[nh],Input,OT[0][nh],Nproc);             
        }
        else
        {
            ParMult(W[nh],OT[0][nh-1],OT[0][nh],Nproc);
        }
        for (size_t nn=0;nn<N[nh];nn++)
        {
            OT[0][nh][nn] = sigmoid(OT[0][nh][nn]+B[nh][nn]);
        }
    }
    Predict = OT[0][NH];
}

inline void Domain::Save (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".hdf5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //Writing data to hdf5 file
    
    hsize_t dims[1];
    size_t val[1];
    String dsname;

    dims[0] = 1;
    val[0] = NI;
    dsname.Printf("NI");
    H5LTmake_dataset(file_id,dsname.CStr(),1,dims,H5T_NATIVE_ULONG,val);
    
    dims[0] = 1;
    val[0] = NO;
    dsname.Printf("NO");
    H5LTmake_dataset(file_id,dsname.CStr(),1,dims,H5T_NATIVE_ULONG,val);
    
    dims[0] = 1;
    val[0] = NH;
    dsname.Printf("NH");
    H5LTmake_dataset(file_id,dsname.CStr(),1,dims,H5T_NATIVE_ULONG,val);

    size_t Nu[NH+1];
    hsize_t dimsw[2];
    for (size_t nh=0;nh<=NH;nh++)
    {
        dsname.Printf("W_%d",nh);
        if (nh==0)
        {
            dimsw[0] = NI;
            dimsw[1] = N[nh];
        }
        else
        {
            dimsw[0] = N[nh-1];
            dimsw[1] = N[nh];
        }
        Nu[nh]   = N[nh];
        double Wu[dimsw[0]*dimsw[1]];
        size_t nu = 0;
        for (size_t ni=0;ni<dimsw[0];ni++)
        for (size_t nj=0;nj<dimsw[1];nj++)
        {
            Wu[nu] = W[nh][ni][nj];
            nu++;
        }
        H5LTmake_dataset_double(file_id,dsname.CStr(),2,dimsw,Wu);

        dsname.Printf("B_%d",nh);
        dims[0] = N[nh];
        double Bu[dims[0]];
        for (size_t nn=0;nn<N[nh];nn++)
        {
            Bu[nn] = B[nh][nn];
        }
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Bu);
    }

    dims[0] = NH+1;
    dsname.Printf("NN");
    H5LTmake_dataset(file_id,dsname.CStr(),1,dims,H5T_NATIVE_ULONG,Nu);

    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);  
}

inline void Domain::Load (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    size_t val[1];
    H5LTread_dataset(file_id,"/NI",H5T_NATIVE_ULONG,val);
    NI = val[0];

    H5LTread_dataset(file_id,"/NO",H5T_NATIVE_ULONG,val);
    NO = val[0];

    H5LTread_dataset(file_id,"/NH",H5T_NATIVE_ULONG,val);
    NH = val[0];

    W.Resize(NH+1);
    B.Resize(NH+1);
    N.Resize(NH+1);
    OT.Resize(1);
    OT[0].Resize(NH+1);


    size_t valn[NH+1];
    H5LTread_dataset(file_id,"/NN",H5T_NATIVE_ULONG,valn);
    
    for (size_t nh=0;nh<=NH;nh++)
    {
        N[nh] = valn[nh];

        String dsname;
        dsname.Printf("W_%d",nh);
        size_t NL,NN;
        if (nh==0)
        {
            NL = NI;
            NN = N[nh];
        }
        else
        {
            NL = N[nh-1];
            NN = N[nh];
        }
        double Wu[NL*NN];
        H5LTread_dataset_double(file_id,dsname.CStr(),Wu);
        
        size_t nu=0;
        W[nh].Resize(NL);
        for (size_t nl=0;nl<NL;nl++)
        {
            W[nh][nl].Resize(NN); 
            for (size_t nn=0;nn<NN;nn++)
            {
                W[nh][nl][nn] = Wu[nu];
                nu++;
            }
        }

        dsname.Printf("B_%d",nh);
        double Bu[N[nh]];
        H5LTread_dataset_double(file_id,dsname.CStr(),Bu);
        B[nh].Resize(N[nh]);
        
        for (size_t nn=0;nn<N[nh];nn++)
        {
            B[nh][nn] = Bu[nn];
        }

        OT[0][nh].Resize(N[nh]);
    }

    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);

}

inline void Domain::TrainDat (char const * FileKey)
{
    NT = NN::Training(FileKey,I,O);
    OT.Resize(NT);
    E .Resize(NT);
    for (size_t nt=0;nt<NT;nt++)
    {
        OT[nt].Resize(NH+1);
        E [nt].Resize(NH+1);
        for (size_t nh=0;nh<=NH;nh++)
        {
            OT[nt][nh].Resize(N[nh]);
            E [nt][nh].Resize(N[nh]);
        }
    }
}

}

#endif // MECHSYS_NN_DOMAIN_H
