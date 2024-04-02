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


#include <vector>

#include <mechsys/dem/domain.h>
using std::cout;
using std::endl;
using DEM::Domain;

int main( void )
{
    Domain dom;
    //dom.AddVoroPack (-1, 0.1, 10.0, 10.0, 10.0, 3, 3, 3, 1.0, true, 1200, 1.0,1.0);
    //dom.AddCube(-1,Vec3_t(0.0,0.0,-3.0),0.1,1.0,3.0,0.0,&OrthoSys::e0);
    //dom.AddCube(-2,Vec3_t(0.0,0.0, 3.0),0.1,1.0,3.0,0.0,&OrthoSys::e0);
    //dom.AddTetra(-1,Vec3_t(0.0,0.0,-3.0),0.1,1.0,3.0,0.0,&OrthoSys::e0);
    //dom.AddTetra(-2,Vec3_t(0.0,0.0, 3.0),0.1,1.0,3.0,0.0,&OrthoSys::e0);
    //dom.AddRecBox(-1,Vec3_t(0.0,0.0,-3.0),Vec3_t(1.0,2.0,4.0),0.1,3.0,0.0,&OrthoSys::e0);
    //dom.AddRecBox(-2,Vec3_t(0.0,0.0, 3.0),Vec3_t(1.0,2.0,4.0),0.1,3.0,0.0,&OrthoSys::e0);
    //dom.GetParticle(-2)->Erode(0.1);
    dom.AddCylinder(-1,Vec3_t(-1.0,0.0,0.0),0.5,Vec3_t(1.0,0.0,0.0),0.2,0.1,1.0);
    
    dom.Initialize();
    dom.Save("domainwrite");
    dom.WriteXDMF("domainwrite");
    
    return 0;
}
