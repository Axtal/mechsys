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

/////////////////////// Test 03 Loading properly formated training data 

// MechSys
#include <mechsys/nn/Domain.h>

int main(int argc, char **argv) try
{
    size_t Nproc = 1;
    if (argc>1) Nproc = atoi(argv[1]);

    //Input and output arrays
    Array<Array <double> > In;
    Array<Array <double> > Out;
    NN::Training("training",In,Out);

    size_t nh = 1;  //Number of hidden layers
    size_t nn = 50; //Number of neurons per hidden layer
    NN::Domain dom(In,Out,nh,nn);
    dom.Initialize();
    dom.Alpha = 0.5; //learning rate
    dom.Train(400000/*epochs*/,Nproc/*number of cores*/);
    dom.Save("tnn03");

    return 0;
}
MECHSYS_CATCH

