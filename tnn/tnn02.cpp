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

/////////////////////// Test 02 Loading NN 

// MechSys
#include <mechsys/nn/Domain.h>

int main(int argc, char **argv) try
{
    double input = 0.6;
    if (argc>1) input = atof(argv[1]);
    NN::Domain dom("square");
    //Try the trained NN
    Array<double> I(1);
    Array<double> O(1);
    I[0] = input;
    dom.Predict(I,O,1);
    std::cout << O[0] << std::endl;
    dom.Save("tnn02");
    return 0;
}
MECHSYS_CATCH

