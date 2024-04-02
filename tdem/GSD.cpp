/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2015 To Huu Duc                                        *
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

//Packing algorithm, needs mysoil.gsd as input

#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <utility>
#include <boost/date_time/posix_time/posix_time.hpp> // for time to milliseconds

// Voro++
#include "voro++.hh"

// MechSys

#include <mechsys/dem/domain.h>
#include <mechsys/dem/gsd.h>
// using namespace
using namespace std;
using namespace DEM;
// this namespace


int main(int argc, char **argv)
{
	srand(time(NULL));
//	int s=rand()%100;
//	cout << s<<endl;
	DEM::Soil sand;
	DEM::SequentialPackingData packinfo;
	string path;
	if (argc<2)
	{
		path = "./1/parameters.gsd";
	}
	else
	{
		path = argv[1];
	}
	sand.ReadGsd(packinfo, path);
	cout << "Start Time: "<<sand.Now() << endl;
	sand.CompleteGsdData(packinfo);
	//sand.FindMinimumSpecimenSize(packinfo,8);
	sand.PrepareList(packinfo);
	cout << "Number of particles: "<<packinfo.gsd.numberparticles<<endl;
	// sand.PrintOut(packinfo);

	sand.SequentialPacking(packinfo);
	sand.DeleteUnusedParticles(packinfo);
	sand.ConstrictionSizeDistribution(packinfo);
	//packinfo.specimen.AddCylinder(-1000,Vec3_t(0.,0.,-packinfo.boundaryfactors[1]/2),packinfo.boundaryfactors[0]/2,Vec3_t(0.,0.,packinfo.boundaryfactors[1]/2),packinfo.boundaryfactors[0]/2,0.1,2.65);
	sand.SaveDomain(packinfo.specimen,packinfo.gsd.Soilname,1);
	sand.TextOut(packinfo);
	sand.SaveTetrahedraMesh(packinfo);
	cout << "End time: "<<sand.Now()<<endl;
//	sand.ReadDemParameters(packinfo,"sergio");
//	sand.DropDown(packinfo);
//	sand.SaveDomain(packinfo.specimen,packinfo.gsd.Soilname,1);
}
