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

/////////////////////// Test 01 fitting quadratic polynomial

// MechSys
#include <mechsys/nn/Domain.h>

int main(int argc, char **argv) try
{
    size_t Nproc = 1;
    if (argc>1) Nproc = atoi(argv[1]);

    size_t ni = 1;  //Number of inputs
    size_t no = 1;  //Number of outputs
    size_t nh = 1;  //Number of hidden layers
    size_t nn = 20; //Number of neurons per hidden layer
    size_t nt = 100; //Number of training sets

    //Preparing input and output data for training sets
    Array<Array <double> > In;
    Array<Array <double> > Out;
    In .Resize(nt);
    Out.Resize(nt);
    double dx = 1.0/nt;
    for (size_t i=0;i<nt;i++) //For each training set
    {
        In [i].Resize(ni);
        Out[i].Resize(no);
        for (size_t j=0;j<ni;j++)
        {
            In [i][j] = i*dx;
            //Out[i][j] = i*dx;
            Out[i][j] = i*dx*i*dx;
            //std::cout << In[i][j] << " " << Out[i][j] << std::endl;
        }
        //for (size_t j=0;j<no;j++)
        //{
            //Out[i][j] = i*i*dx*dx;
        //}
    }

    NN::Domain dom(In,Out,nh,nn);
    dom.Initialize();
    dom.Alpha = 5.0; //learning rate
    dom.Train(40000/*epochs*/,Nproc/*number of cores*/);
    
    //Try the trained NN
    In[0][0] = 0.4;
    dom.Predict(In[0],Out[0],1);
    std::cout << Out[0][0] << std::endl;

    //dom.Forward(1);
    //dom.BackErr(1);
    //std::cout << dom.E [0][1][0] << std::endl;
    //std::cout << dom.W [0][0][0] << std::endl;
    //std::cout << dom.W [0][0][1] << std::endl;
    //std::cout << dom.E [0][0][0] << std::endl;
    //std::cout << dom.E [0][0][1] << std::endl;
    //std::cout << dom.W [1][0][0] << std::endl;
    //std::cout << dom.W [1][1][0] << std::endl;

    //Saving the NN
    dom.Save("tnn01");
    return 0;
}
MECHSYS_CATCH

