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

/////////////////////// Test 04 Loading NN and trianig data (run after tnn03)

// MechSys
#include <mechsys/nn/Domain.h>

int main(int argc, char **argv) try
{
    size_t Nproc = 1;
    if (argc>1) Nproc = atoi(argv[1]);

    NN::Domain dom("tnn03");
    dom.TrainDat("training");
    dom.Alpha = 0.2; //learning rate
    dom.Train(4000/*epochs*/,Nproc/*number of cores*/);
    dom.Save("tnn04");

    return 0;
}
MECHSYS_CATCH

