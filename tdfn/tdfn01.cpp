/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
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

// MechSys

#include <mechsys/dfn/Loop.h>
#include <chrono>

int main(int argc, char **argv)
try
{
	auto start = std::chrono::steady_clock::now();
	srand((unsigned)time(NULL));
	double random_seed = random_double(0, 1);
	double vari = 0;
	String orientation_distribution;
	String percolation_direction;
	double array13[7]; //mean dip direction (dd), mean dip angle, Fisher constant, min dd, max dd, min dipangle, max dip angle
	DFN::Loop loop_1;
	
	
	std::ifstream oii("tdfn01.inp",std::ios::in);
	if(!oii)
	{
		std::cout<<"Please define a input.inp!\n";
		exit(0);
	}
	oii >> vari; 						oii.ignore(300,'\n');
	oii >> loop_1.Nproc;				oii.ignore(300,'\n');
	oii >> loop_1.times;				oii.ignore(300,'\n');
	oii >> loop_1.nt;					oii.ignore(300,'\n');
	oii >> loop_1.nk;					oii.ignore(300,'\n');
	oii >> loop_1.nv_MC_TIMES; 			oii.ignore(300,'\n');
	oii >> loop_1.nx; 					oii.ignore(300,'\n');
	oii >> loop_1.RatioLR;				oii.ignore(300,'\n');
	oii >> loop_1.R_a;					oii.ignore(300,'\n');
	oii >> loop_1.R_low;				oii.ignore(300,'\n');
	oii >> loop_1.R_up;					oii.ignore(300,'\n');
	oii >> orientation_distribution;	oii.ignore(300,'\n');
	oii >> percolation_direction;		oii.ignore(300,'\n');
	oii >> array13[0];					oii.ignore(300,'\n');
	oii >> array13[1];					oii.ignore(300,'\n');
	oii >> array13[2];					oii.ignore(300,'\n');
	oii >> array13[3];					oii.ignore(300,'\n');
	oii >> array13[4];					oii.ignore(300,'\n');
	oii >> array13[5];					oii.ignore(300,'\n');
	oii >> array13[6];					oii.ignore(300,'\n');
	oii.close();					
	loop_1.Loop_create_DFNs(vari, random_seed, orientation_distribution, percolation_direction, array13);
	
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double, std::micro> elapsed = end - start; // std::micro time (us)
	std::cout << "Running time: " << (double)(elapsed.count() * 1.0) * (0.000001) << "s" << std::endl;
	return 0;
}
MECHSYS_CATCH
