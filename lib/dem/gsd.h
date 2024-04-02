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



// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <utility>
#include <boost/date_time/posix_time/posix_time.hpp> // for time to milliseconds
#include <string.h>
#include <stdio.h>

// Voro++
#include "voro++.hh"

// MechSys

#include <mechsys/dem/domain.h>
// using namespace
using namespace std;
// this namespace
namespace DEM
{
struct DemParameters
{
	double Kn;
	double Kt;
	double Gn;
	double Gt;
	double Gv;
	double Mu;
	double dt;
	double tf;
	double dtout;
	double roundratio;
	double Alpha;
	double Kinematicenergy;
	int numberprocessors;
	int numbershapes;
	int numberintervals;
	int visualization;
	int starttag;
	string FileKey;
};
struct SequentialFace
{
	bool boundaryface;
	int points[3];
	int point;
	int faceuse;
	double constrictionsize;		// radius not diameter
	Vec3_t normal;
	Vec3_t constrictioncentre;
};
struct SequentialTetrahedron
{
	int points[4];
	int faces[4];
	double poreradius;
	Vec3_t porecentre;
};
struct File
{
	string filein;
	string fileout;
	string parameterfile;
};
struct Interval
{
	bool ability;
	int firstparticle;														// particle number
	int lastparticle;
	int usingparticle;
	double begin;															// fraction
	double begindiameter;
	double diameter;														//mean diameter;
	double end;
	double enddiameter;
	double fraction;
	double mass;
	double numberparticles;
	double specialgravity;
	double volumeparticles;
	double shaperatios[10]; 	// 0-sphere 1- cube 2- tetrahedron 3-rectangular box 9-user defined
	int shapestates[10][4];	// 0 - numberparticles, 1- begin ID, 2 - end ID, 3 - current using
};
struct Gsd
{
	string Soilname;
	string arrangement; 	// particle arrangement 0 - customized, 1- layer, 2 - discrete
	string boundaryfile;
	double density;
	double mass;
	double porosity;
	double specialgravity;	// special gravity 0 - provided by intervals
	double volume;
	int maxfraction; 		// max percentage per intervals
	int numberboundaries;
	int numberintervals;
	int numberinters;
	int numberparticles;
	int numbershapes;
	double shapem[10];	//ratio between mass of shape and basic size
	double rectangularboxratios[3];
	vector <Interval> intervals;		// readin intervals
	vector <Interval> inters;			// calculating intervals
};
struct Csd
{
	vector <double> size;
	vector <double> fraction;
};
struct SequentialPackingData
{
	vector <SequentialFace> faces;
	vector <SequentialTetrahedron> tetrahedra;
	vector <bool> particleuses;
	vector <int> boundary;
	vector <int> threadnumbers;
	Vec3_t localroot;
	Vec3_t localsystem[3];
	Vec3_t position;
	bool check;
	bool checkboundary;
	bool checkoverlap;
	bool checkradius;
	bool facereuse;
	bool particlereuse;
	double approximation;
	double boundaryfactors[10];
	double boundarysizes[10];
	int facereusenumber;
	int numberopenfaces;
	int numberfaces;
	int numberprocessors;
	int numberunusedparticles;
	int overlappingpoint;
	int randomness;
	int usingface;
	int usinginter;
	int usingparticle;
	int temproraryparticle;
	string boundarytype;
	DEM::Csd csd;
	DEM::Gsd gsd;
	DEM::Domain specimen;
	DemParameters para;
};

// main classes
class Soil
{
public:
    // typedefs
    typedef void (*ptFun_t) (Soil & soil, void * UserData);

    // Constructor
    Soil(void * UserData=NULL);

    // Destructor
    ~Soil();
    // List of functions
    double ConstrictionSize(DEM::SequentialPackingData&packinfo, int p0, int p1, int p2);
    double PlaneConstriction(DEM::SequentialPackingData&packinfo, int face, int p3);
    double ShapeMass(DEM::SequentialPackingData&packinfo, int shape, double size, double density);
    int GetInterval(DEM::SequentialPackingData&packinfo, int p0);
    string Now();
    void AddParticle( DEM::SequentialPackingData&packinfo, int inter, int shape, double size, double density, Vec3_t position, double roundness=0.005);
    void BasicTetrahedron (DEM::SequentialPackingData&packinfo, int p[4]);
    void CheckBoundary(DEM::SequentialPackingData&packinfo, int p3);
    void CheckBoundaryConstriction(DEM::SequentialPackingData&packinfo, double size, Vec3_t position);
    void CheckOverlap(DEM::SequentialPackingData&packinfo, int p3);
    void CheckOverlapConstriction(DEM::SequentialPackingData&packinfo, double size, Vec3_t position);
    void CloseSequentialFace(DEM::SequentialPackingData&packinfo, int face);
    void CompleteGsdData(DEM::SequentialPackingData&packinfo);
    void CompleteInterval(DEM::SequentialPackingData&packinfo);
    void ConstrictionPlane(DEM::SequentialPackingData&packinfo, int face, int p3);
    void ConstrictionSizeDistribution(DEM::SequentialPackingData&packinfo);
    void CreateSequentialFace(DEM::SequentialPackingData&packinfo, int p0, int p1, int p2,int p3, bool sort=true);
    void CreateSequentialTetrahedron(DEM::SequentialPackingData&packinfo, int face, int p3);
    void DeleteUnusedParticles(DEM::SequentialPackingData&packinfo);
    void DrawBoundary(DEM::SequentialPackingData&packinfo);
    void DropDown(DEM::SequentialPackingData&packinfo);
    void EstablishLocalSystem(DEM::SequentialPackingData&packinfo, int face);
    void FindMinimumSpecimenSize(DEM::SequentialPackingData&packinfo, int coarsest=1);
    void FrozenTag(DEM::SequentialPackingData&packinfo, int tag);
    void MoveOn(DEM::SequentialPackingData&packinfo, int p3,int face);
    void ParticleInfo(DEM::SequentialPackingData&packinfo,int p0);
    void PlaceSequentialParticle(DEM::SequentialPackingData&packinfo,int face, int p3);							// for layer arrangement
    void PrepareList(DEM::SequentialPackingData&packinfo, bool randomness=false);
    void PrintOut(DEM::SequentialPackingData&packinfo);
	void PutSequentialParticle(DEM::SequentialPackingData&packinfo, int face, int p3);
	void ReadDemParameters(DEM::SequentialPackingData&packinfo, string FileKey="sand");
    void ReadGsd(DEM::SequentialPackingData&packinfo,string FileKey="sand");
    void SaveDomain(DEM::Domain&specimen,string filename, int outputtype);
    void SaveTetrahedraMesh(DEM::SequentialPackingData&packinfo);
    void SequentialPacking(DEM::SequentialPackingData&packinfo);
    void TextConstriction(DEM::SequentialPackingData&packinfo, int type=0);
    void TextOut(DEM::SequentialPackingData&packinfo);
    void TrySequentialParticle(DEM::SequentialPackingData&packinfo, int face, int p3);
    void UseParticle(DEM::SequentialPackingData&packinfo, int p3);
};

// Constructor & Destructor
inline Soil::Soil(void * UD)
{}

inline Soil::~Soil()
{}

inline double Soil::ConstrictionSize(DEM::SequentialPackingData&packinfo,int p0, int p1, int p2)
{
	Vec3_t x = packinfo.specimen.Particles[p1]->x-packinfo.specimen.Particles[p0]->x;
	double x2 = norm(x);																										// get coordinates of particles in local coordinate system
	x = x/x2;
	double x3 = dot(packinfo.specimen.Particles[p2]->x - packinfo.specimen.Particles[p0]->x, x);
	Vec3_t y = (packinfo.specimen.Particles[p2]->x-packinfo.specimen.Particles[p0]->x-x3*x);
	double y3 = norm(y);
	y=y/y3;
	double a1 = (packinfo.specimen.Particles[p0]->Dmax-packinfo.specimen.Particles[p1]->Dmax)/x2;
	double a2 = (packinfo.specimen.Particles[p0]->Dmax -packinfo.specimen.Particles[p2]->Dmax-a1*x3)/y3;
	double b1 = (pow(x2,2.0)+pow(packinfo.specimen.Particles[p0]->Dmax,2.0)-pow(packinfo.specimen.Particles[p1]->Dmax,2.0))/2.0/x2;
	double b2 = (pow(x3,2.0)+pow(y3,2.0)+pow(packinfo.specimen.Particles[p0]->Dmax,2.0)-pow(packinfo.specimen.Particles[p2]->Dmax,2.0)-2*b1*x3)/2/y3;
	double a = pow(a1,2.0)+pow(a2,2.0)-1;
	double b = a1*b1+a2*b2-packinfo.specimen.Particles[p0]->Dmax;
	double c = pow(b1, 2.0) +pow(b2, 2.0) - pow(packinfo.specimen.Particles[p0]->Dmax, 2.0);
	if (b*b-a*c >0.)
		{
			double constriction = (-b-pow(b*b-a*c,0.5))/a;
			if (constriction > packinfo.approximation)
			{
				packinfo.position = (a1*constriction+b1)*x+(a2*constriction+b2)*y+packinfo.specimen.Particles[p0]->x;
				return (constriction);
			}
			else
			{
				cout << "constriction problems \n";
				packinfo.check =false;
				return (0.0);
			}
		}
	else
		{
			cout << packinfo.numberfaces<<" no constriction \n";
			return (0.0);
		}
}

inline double Soil::PlaneConstriction(DEM::SequentialPackingData&packinfo, int face, int p3)
{
	// p0, p1, p3 -> new constriction, p0,p1,p2 -> plane on which centre of the constriction
	EstablishLocalSystem(packinfo,face);
	double planeconstriction;
	double x2 = norm(packinfo.specimen.Particles[packinfo.faces[face].points[1]]->x -packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x);
	double x3 = dot(packinfo.specimen.Particles[p3]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.localsystem[0]);
	double y3 = dot(packinfo.specimen.Particles[p3]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x, packinfo.localsystem[1]);
	double z3 = dot(packinfo.specimen.Particles[p3]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x, packinfo.localsystem[2]);
	double a1 = (packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax-packinfo.specimen.Particles[packinfo.faces[face].points[1]]->Dmax)/x2;
	double a2 = (pow(packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax,2.0)+pow(x2,2.0)-pow(packinfo.specimen.Particles[packinfo.faces[face].points[1]]->Dmax,2.0))/2/x2;
	double b1 = (packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax-packinfo.specimen.Particles[p3]->Dmax-a1*x3)/y3;
	double b2 = (pow(packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax,2.0)+pow(x3,2.0)+pow(y3,2.0)+pow(z3,2.0)-pow(packinfo.specimen.Particles[p3]->Dmax,2.0)-2*a2*x3)/2/y3;
	double a = pow(a1,2.0)+pow(b1,2.0)-1;
	double b = a1*a2+b1*b2-packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax;
	double c = pow(a2,2.0)+pow(b2,2.0)-pow(packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax,2.0);
	if (pow(b,2.0)-a*c>0)
	{
		planeconstriction=(-b-pow(pow(b,2.0)-a*c,0.5))/a;
		if (planeconstriction > packinfo.approximation)
		{
			packinfo.position = (a1*planeconstriction+a2)*packinfo.localsystem[0]+(b1*planeconstriction+b2)*packinfo.localsystem[1]+packinfo.localroot;
			return (planeconstriction);
		}
		else
		{
			planeconstriction=(-b+pow(pow(b,2.0)-a*c,0.5))/a;
			if (planeconstriction > packinfo.approximation)
			{
				cout<<"Wrong sign \n";
				packinfo.position = (a1*planeconstriction+a2)*packinfo.localsystem[0]+(b1*planeconstriction+b2)*packinfo.localsystem[1]+packinfo.localroot;
				return (planeconstriction);
			}
			else
			{
				packinfo.faces[face].boundaryface=true; 	// both results smaller than 0.
				cout << "Plane constriction error";
				return(0.0);
			}
		}
	}
	else
	{
		packinfo.faces[face].boundaryface=true; // p3 is behind a particle
		cout<<"Plane constriction error \n";
		return(0.0);
	}
}

inline double Soil::ShapeMass(DEM::SequentialPackingData&packinfo, int shape, double size, double density)
{
	double mass;
	switch (shape)
	{
	case 0:
		mass = density * acos(0.) * pow(size,3.) / 3.;
		break;
	case 1:
		mass = density * pow(size,3.);
		break;
	case 2:
		mass = density * pow(size,3.) / 3.;
		break;
	case 3:
		mass = density * pow(pow(2.,0.5)*size/(1+packinfo.gsd.rectangularboxratios[1]),3.) *packinfo.gsd.rectangularboxratios[1] * packinfo.gsd.rectangularboxratios[2];
		break;
	default :
		break;
	}
	return (mass);
}

inline int Soil::GetInterval(DEM::SequentialPackingData&packinfo, int p0)
{
	int inter = floor(abs(packinfo.specimen.Particles[p0]->Tag)/packinfo.gsd.numbershapes);
	return (inter);
}

inline string Soil::Now()
{
	const boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();	    // Get the time offset in current day
	const boost::posix_time::time_duration td = now.time_of_day();
	const long hours        = td.hours();
	const long minutes      = td.minutes();
	const long seconds      = td.seconds();
	const long milliseconds = td.total_milliseconds()-((hours * 3600 + minutes * 60 + seconds) * 1000);
	char buf[40];
	sprintf(buf, "%02ld:%02ld:%02ld.%03ld",hours, minutes, seconds, milliseconds);
	return buf;
}
// Procedures
inline void Soil::AddParticle( DEM::SequentialPackingData&packinfo, int inter, int shape, double size, double density, Vec3_t position, double roundness)
{
	switch (shape)
	{
	case 0:
		packinfo.specimen.AddSphere(packinfo.gsd.numbershapes*inter,position,size/2,density);
		break;
	case 1:
		packinfo.specimen.AddCube(-packinfo.gsd.numbershapes*inter-shape,position,roundness,size,density);
		packinfo.specimen.Particles[packinfo.specimen.Particles.Size()-1]->Erode(roundness);
		break;
	case 2:
		packinfo.specimen.AddTetra(-packinfo.gsd.numbershapes*inter-shape,position,roundness,pow(2.,0.5)*size,density);
		packinfo.specimen.Particles[packinfo.specimen.Particles.Size()-1]->Erode(roundness);
		break;
	case 3:
		double smallsize;
		smallsize= size/(1+packinfo.gsd.rectangularboxratios[1]);
		packinfo.specimen.AddRecBox(-packinfo.gsd.numbershapes*inter-shape,position,Vec3_t (smallsize,smallsize*packinfo.gsd.rectangularboxratios[1],smallsize*packinfo.gsd.rectangularboxratios[2]),roundness,density);
		packinfo.specimen.Particles[packinfo.specimen.Particles.Size()-1]->Erode(roundness);
		break;
	default:
		break;
	}
}

inline void Soil::BasicTetrahedron(DEM::SequentialPackingData&packinfo, int p[4])
{
	packinfo.numberunusedparticles = packinfo.specimen.Particles.Size();
	packinfo.numberopenfaces=0;
	cout << "Creating basic tetrahedron \n";
	double r[4];
	for (int i=0;i<4;i++)
	{
		r[i]=packinfo.specimen.Particles[p[i]]->Dmax;
	}
	packinfo.specimen.Particles[p[0]]->Position(Vec3_t(0.,0.,0.));
	packinfo.specimen.Particles[p[1]]->Position(Vec3_t(r[0]+r[1],0.,0.));
	double x = (pow(r[0]+r[1],2)+pow(r[0]+r[2],2)-pow(r[1]+r[2],2))/2/(r[0]+r[1]);
	packinfo.specimen.Particles[p[2]]->Position(Vec3_t(x, pow(pow(r[0]+r[2],2)-pow(x,2),0.5),0.));
	CreateSequentialFace(packinfo,p[0],p[1],p[2],p[3],true);
	PutSequentialParticle(packinfo,0,packinfo.usingparticle-3);
	CreateSequentialTetrahedron(packinfo, 0,packinfo.usingparticle-3);
	packinfo.faces[0].faceuse=1;																											// mark used side
	packinfo.usingface=0;
	packinfo.numberfaces =4;				// or 1 ?
	packinfo.numberopenfaces =4;			// or 1 ?
	for (int i=0; i<4; i++)																							//mark used particles
		{
			UseParticle(packinfo, p[i]);
		}
}

inline void Soil::CheckBoundary( DEM::SequentialPackingData&packinfo, int p3)
{
	packinfo.checkboundary=true;
	if (packinfo.boundarytype =="File")
	{
		for (int i=0; i<packinfo.boundary.size();i++)
		{

		}
	}
	else if (packinfo.boundarytype=="Sphere")
	{
		if (norm(packinfo.specimen.Particles[p3]->x)+ packinfo.specimen.Particles[p3]->Dmax >packinfo.boundaryfactors[0]/2)
			{
				packinfo.checkboundary=false;
			}
	}
	else if (packinfo.boundarytype=="Cube")
	{
		for (int i = 0; i < 3; i++)
			{
				if (abs(packinfo.specimen.Particles[p3]->x(i))+packinfo.specimen.Particles[p3]->Dmax > packinfo.boundaryfactors[i]/2)
					{
						packinfo.checkboundary= false;
					}
			}
	}
	else if (packinfo.boundarytype== "Cylinder")
	{
		if (abs(packinfo.specimen.Particles[p3]->x(2))+packinfo.specimen.Particles[p3]->Dmax > packinfo.boundaryfactors[1]/2)
			{
				packinfo.checkboundary=false;
			}
		else if (pow(pow(packinfo.specimen.Particles[p3]->x(0),2.)+pow(packinfo.specimen.Particles[p3]->x(1),2.),0.5)+packinfo.specimen.Particles[p3]->Dmax > packinfo.boundaryfactors[0]/2)
			{
				packinfo.checkboundary=false;
			}
	}
	else if (packinfo.boundarytype=="Box")
	{
		for (int i = 0; i < 3; i++)
			{
				if (abs(packinfo.specimen.Particles[p3]->x(i))+packinfo.specimen.Particles[p3]->Dmax > packinfo.boundaryfactors[i])
					{
						packinfo.checkboundary= false;
					}
			}
	}
}

inline void Soil::CheckBoundaryConstriction(DEM::SequentialPackingData&packinfo, double radius, Vec3_t position)
{
	packinfo.checkboundary=true;
	if (packinfo.boundarytype =="File")
	{
		for (int i=0; i<packinfo.boundary.size();i++)
		{

		}
	}
	else if (packinfo.boundarytype=="Sphere")
	{
		if (norm(position)+ radius >packinfo.boundaryfactors[0]/2)
			{
				packinfo.checkboundary=false;
			}
	}
	else if (packinfo.boundarytype=="Cube")
	{
		for (int i = 0; i < 3; i++)
			{
				if (abs(position(i)) + radius > packinfo.boundaryfactors[i]/2)
					{
						packinfo.checkboundary= false;
					}
			}
	}
	else if (packinfo.boundarytype== "Cylinder")
	{
		if (abs(position(2))+ radius > packinfo.boundaryfactors[1]/2)
			{
				packinfo.checkboundary=false;
			}
		else if (pow(pow(position(0),2.)+pow(position(1),2.),0.5)+radius> packinfo.boundaryfactors[0]/2)
			{
				packinfo.checkboundary=false;
			}
	}
	else if (packinfo.boundarytype=="Box")
	{
		for (int i = 0; i < 3; i++)
			{
				if (abs(position(i))+ radius > packinfo.boundaryfactors[i])
					{
						packinfo.checkboundary= false;
					}
			}
	}
}

inline void Soil::CheckOverlap( DEM::SequentialPackingData&packinfo, int p3)
{
//	bool checkoverlap =true;
//	for (int i =0; i<packinfo.numberprocessors; i++)
//	{
//		pthread_t()
//	}
//	return (checkoverlap);
	packinfo.checkoverlap=true;
	double overlapdistance = -packinfo.approximation;
	double distance;
	for (int i=0;i<packinfo.specimen.Particles.Size();i++)
	{
		if (packinfo.particleuses[i]and(!(i==p3)))
		{
			distance = norm(packinfo.specimen.Particles[p3]->x-packinfo.specimen.Particles[i]->x)-packinfo.specimen.Particles[p3]->Dmax-packinfo.specimen.Particles[i]->Dmax;
			if (distance <overlapdistance)
			{
				packinfo.overlappingpoint = i;
				overlapdistance =distance;
				packinfo.checkoverlap=false;
			}
		}
	}
}

inline void Soil::CheckOverlapConstriction( DEM::SequentialPackingData&packinfo, double radius, Vec3_t position)
{
//	bool checkoverlap =true;
//	for (int i =0; i<packinfo.numberprocessors; i++)
//	{
//		pthread_t()
//	}
//	return (checkoverlap);
	packinfo.checkoverlap=true;
	double overlapdistance = -packinfo.approximation;
	double distance;
	for (int i=0;i<packinfo.specimen.Particles.Size();i++)
	{
		if (packinfo.particleuses[i])
		{
			distance = norm(position-packinfo.specimen.Particles[i]->x)-radius-packinfo.specimen.Particles[i]->Dmax;
			if (distance <overlapdistance)
			{
				packinfo.overlappingpoint = i;
				overlapdistance =distance;
				packinfo.checkoverlap=false;
			}
		}
	}
}

inline void Soil::CloseSequentialFace(DEM::SequentialPackingData&packinfo, int face)
{
	packinfo.faces[face].faceuse =0;
	packinfo.numberopenfaces -=1;
}

inline void Soil::CompleteGsdData(DEM::SequentialPackingData&packinfo)
{
    // complete data
    if (packinfo.gsd.intervals[0].begin >0)
    {
    	packinfo.gsd.intervals.insert(packinfo.gsd.intervals.begin(),packinfo.gsd.intervals[0]);
    	packinfo.gsd.intervals[0].begin =0;
    }
    if (packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].begin==100)
    {
    	packinfo.gsd.intervals[packinfo.gsd.numberintervals-2].end=100;
    	packinfo.gsd.intervals[packinfo.gsd.numberintervals-2].enddiameter=packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].begindiameter;
    	packinfo.gsd.intervals.pop_back();
    	packinfo.gsd.numberintervals--;
    }
    else
    {
    	packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].end=100;
    	packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].enddiameter=packinfo.gsd.intervals[packinfo.gsd.numberintervals-1].begindiameter;
    }
    for (int i=0; i<packinfo.gsd.numberintervals-1;i++)
    {
    	packinfo.gsd.intervals[i].end = packinfo.gsd.intervals[i+1].begin;
    	packinfo.gsd.intervals[i].enddiameter=packinfo.gsd.intervals[i+1].begindiameter;
    }
    packinfo.gsd.numberinters =0;
    double count1=0;
    DEM::Interval in;
	for (int i=0; i<packinfo.gsd.numberintervals; i++)
	{
		if ((packinfo.gsd.intervals[i].end -packinfo.gsd.intervals[i].begin)>packinfo.gsd.maxfraction)
		{
			count1 = packinfo.gsd.intervals[i].begin;
			while ((count1 +packinfo.gsd.maxfraction)<packinfo.gsd.intervals[i].end)
			{
				in.begin =count1;
				if (packinfo.gsd.numberinters>0)
				{
					in.begindiameter=packinfo.gsd.inters[packinfo.gsd.numberinters-1].enddiameter;
				}
				else
				{
					in.begindiameter = packinfo.gsd.intervals[0].begindiameter;
				}
				count1 +=packinfo.gsd.maxfraction;
				in.end=count1;
				in.enddiameter= pow(10,log10(packinfo.gsd.intervals[i].begindiameter)+(in.end-packinfo.gsd.intervals[i].begin)/(packinfo.gsd.intervals[i].end-packinfo.gsd.intervals[i].begin)*(log10(packinfo.gsd.intervals[i].enddiameter)-log10(packinfo.gsd.intervals[i].begindiameter)));
				for (int j=0; j<packinfo.gsd.numbershapes;j++)
				{
					in.shaperatios[j]=packinfo.gsd.intervals[i].shaperatios[j];
				}
				packinfo.gsd.inters.push_back(in);
				packinfo.gsd.numberinters++;
			}
			if (count1 <packinfo.gsd.intervals[i].end)
			{
				in.begin=count1;
				in.begindiameter=packinfo.gsd.inters[packinfo.gsd.numberinters-1].enddiameter;
				count1 = packinfo.gsd.intervals[i].end;
				in.end=count1;
				in.enddiameter=packinfo.gsd.intervals[i].enddiameter;
				for (int j=0; j<packinfo.gsd.numbershapes;j++)
				{
					in.shaperatios[j]=packinfo.gsd.intervals[i].shaperatios[j];
				}
				packinfo.gsd.inters.push_back(in);
				packinfo.gsd.numberinters++;
			}

		}
		else
		{
			count1=packinfo.gsd.intervals[i].end;
			packinfo.gsd.inters.push_back(packinfo.gsd.intervals[i]);
			packinfo.gsd.numberinters ++;
		}
	}
	if (packinfo.boundarytype=="File")
	{
		// must calculate volume here
	}
	else if (packinfo.boundarytype=="Sphere")
	{
		packinfo.gsd.volume=acos(0.)/3*pow(packinfo.boundarysizes[0],3.);
	}
	else if (packinfo.boundarytype=="Cube")
	{
		packinfo.gsd.volume=pow(packinfo.boundarysizes[0],3.);
	}
	else if (packinfo.boundarytype=="Cylinder")
	{
		packinfo.gsd.volume=acos(0)/2.*pow(packinfo.boundarysizes[0],2.)*packinfo.boundarysizes[1];
	}
	else if (packinfo.boundarytype=="Box")
	{
		packinfo.gsd.volume=packinfo.boundarysizes[0]*packinfo.boundarysizes[1]*packinfo.boundarysizes[2];
	}
	if (packinfo.gsd.specialgravity==0)
	{
		double vol=0.;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			vol+=(packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin)/packinfo.gsd.inters[i].specialgravity;
		}
		packinfo.gsd.mass=packinfo.gsd.volume/(vol/100/(1-packinfo.gsd.porosity));
		packinfo.gsd.density = packinfo.gsd.mass/packinfo.gsd.volume;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			packinfo.gsd.inters[i].mass= (packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin)/100*packinfo.gsd.mass;
			packinfo.gsd.inters[i].volumeparticles=packinfo.gsd.inters[i].mass/packinfo.gsd.inters[i].specialgravity;
		}
	}
	else
	{
		packinfo.gsd.mass=packinfo.gsd.volume*packinfo.gsd.specialgravity*(1-packinfo.gsd.porosity);
		packinfo.gsd.density=packinfo.gsd.mass/packinfo.gsd.volume;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			packinfo.gsd.inters[i].mass=packinfo.gsd.mass/100*(packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin);
			packinfo.gsd.inters[i].specialgravity= packinfo.gsd.specialgravity;
			packinfo.gsd.inters[i].volumeparticles=packinfo.gsd.inters[i].mass/packinfo.gsd.specialgravity;
		}
	}
	for (int i=0;i<10;i++)
	{
		packinfo.boundaryfactors[i]*=packinfo.boundarysizes[i];
	}
	for (int i=0; i<3;i++)
	{
		packinfo.gsd.rectangularboxratios[i]/=packinfo.gsd.rectangularboxratios[0];
	}
	// Textout(packinfo);						//test completion process
}

inline void Soil::CompleteInterval(DEM::SequentialPackingData&packinfo)
{
    // complete data
	packinfo.gsd.intervals[0].end=packinfo.gsd.intervals[0].begin;
	packinfo.gsd.intervals[0].begin=0;
    for (int i=1; i<packinfo.gsd.numberintervals;i++)
    {
    	packinfo.gsd.intervals[i].end = packinfo.gsd.intervals[i].begin;
    	packinfo.gsd.intervals[i].begin=packinfo.gsd.intervals[i-1].end;
    	packinfo.gsd.intervals[i].enddiameter=packinfo.gsd.intervals[i].begindiameter;
    }
    packinfo.gsd.numberinters =0;
    double count1=0;
    DEM::Interval in;
	for (int i=0; i<packinfo.gsd.numberintervals; i++)
	{
		packinfo.gsd.inters.push_back(packinfo.gsd.intervals[i]);
		packinfo.gsd.numberinters ++;
	}
	if (packinfo.boundarytype=="File")
	{
		// must calculate volume here
	}
	else if (packinfo.boundarytype=="Sphere")
	{
		packinfo.gsd.volume=acos(0.)/3*pow(packinfo.boundarysizes[0],3.);
	}
	else if (packinfo.boundarytype=="Cube")
	{
		packinfo.gsd.volume=pow(packinfo.boundarysizes[0],3.);
	}
	else if (packinfo.boundarytype=="Cylinder")
	{
		packinfo.gsd.volume=acos(0)/2.*pow(packinfo.boundarysizes[0],2.)*packinfo.boundarysizes[1];
	}
	else if (packinfo.boundarytype=="Box")
	{
		packinfo.gsd.volume=packinfo.boundarysizes[0]*packinfo.boundarysizes[1]*packinfo.boundarysizes[2];
	}
	if (packinfo.gsd.specialgravity==0)
	{
		double vol=0.;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			vol+=(packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin)/packinfo.gsd.inters[i].specialgravity;
		}
		packinfo.gsd.mass=packinfo.gsd.volume/(vol/100/(1-packinfo.gsd.porosity));
		packinfo.gsd.density = packinfo.gsd.mass/packinfo.gsd.volume;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			packinfo.gsd.inters[i].mass= (packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin)/100*packinfo.gsd.mass;
			packinfo.gsd.inters[i].volumeparticles=packinfo.gsd.inters[i].mass/packinfo.gsd.inters[i].specialgravity;
		}
	}
	else
	{
		packinfo.gsd.mass=packinfo.gsd.volume*packinfo.gsd.specialgravity*(1-packinfo.gsd.porosity);
		packinfo.gsd.density=packinfo.gsd.mass/packinfo.gsd.volume;
		for (int i=0;i<packinfo.gsd.numberinters;i++)
		{
			packinfo.gsd.inters[i].mass=packinfo.gsd.mass/100*(packinfo.gsd.inters[i].end-packinfo.gsd.inters[i].begin);
			packinfo.gsd.inters[i].specialgravity= packinfo.gsd.specialgravity;
			packinfo.gsd.inters[i].volumeparticles=packinfo.gsd.inters[i].mass/packinfo.gsd.specialgravity;
		}
	}
	for (int i=0;i<10;i++)
	{
		packinfo.boundaryfactors[i]*=packinfo.boundarysizes[i];
	}
	for (int i=0; i<3;i++)
	{
		packinfo.gsd.rectangularboxratios[i]/=packinfo.gsd.rectangularboxratios[0];
	}
	// Textout(packinfo);						//test completion process
}

//inline void Soil::ConstrictionPlane(DEM::SequentialPackingData&packinfo, int face, int p3)
//{
//	// p0, p1, p3 -> new constriction, p0,p1,p2 -> plane on which centre of the constriction
//	EstablishLocalSystem(packinfo,face);
//	double x2 = norm(packinfo.specimen.Particles[packinfo.faces[face].points[1]]->x -packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x);
//	double x3 = dot(packinfo.specimen.Particles[p3]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.localsystem[0]);
//	double y3 = dot(packinfo.specimen.Particles[p3]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x, packinfo.localsystem[1]);
//	double z3 = dot(packinfo.specimen.Particles[p3]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x, packinfo.localsystem[2]);
//	double a1 = (packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax-packinfo.specimen.Particles[packinfo.faces[face].points[1]]->Dmax)/x2;
//	double a2 = (pow(packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax,2.0)+pow(x2,2.0)-pow(packinfo.specimen.Particles[packinfo.faces[face].points[1]]->Dmax,2.0))/2/x2;
//	double b1 = (packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax-packinfo.specimen.Particles[p3]->Dmax-a1*x3)/y3;
//	double b2 = (pow(packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax,2.0)+pow(x3,2.0)+pow(y3,2.0)+pow(z3,2.0)-pow(packinfo.specimen.Particles[p3]->Dmax,2.0)-2*a2*x3)/2/y3;
//	double a = pow(a1,2.0)+pow(b1,2.0)-1;
//	double b = a1*a2+b1*b2-packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax;
//	double c = pow(a2,2.0)+pow(b2,2.0)-pow(packinfo.specimen.Particles[packinfo.faces[face].points[0]]->Dmax,2.0);
//	if (pow(b,2.0)-a*c>0)
//	{
//		packinfo.faces[face].constrictionsize=(-b-pow(pow(b,2.0)-a*c,0.5))/a;
//		double r =(-b+pow(pow(b,2.0)-a*c,0.5))/a;
//		if ((r>0) and (r<packinfo.faces[face].constrictionsize))
//		{
//			cout<<"Wrong sign \n";
//		}
//	}
//	else
//	{
//		packinfo.faces[face].boundaryface=true; // p3 is behind a particle
//		cout<<"Plane constriction error \n";
//	}
//}

inline void Soil::ConstrictionSizeDistribution(DEM::SequentialPackingData&packinfo)
{
	bool check;
	int countboundary=0;
	int countoverlap=0;
	int countmultipleoverlap=0;
	for (int i=0; i<packinfo.numberfaces; i++)
	{
		if (packinfo.faces[i].constrictionsize>0)
		{
			check =true;
			CheckBoundaryConstriction(packinfo,packinfo.faces[i].constrictionsize,packinfo.faces[i].constrictioncentre);
			if (packinfo.checkboundary==false)
			{
				packinfo.faces[i].boundaryface=true;
				countboundary+=1;
				check=false;
				continue;
			}
			for (int j=0; j<packinfo.gsd.numberparticles;j++)		// cannot use checkoverlapconstriction because the constriction is resized continuously
			{
				if (norm(packinfo.faces[i].constrictioncentre-packinfo.specimen.Particles[j]->x)-(packinfo.faces[i].constrictionsize+packinfo.specimen.Particles[j]->Dmax)<-packinfo.approximation)
				{
					if (check==true)
					{
						countoverlap+=1;
						check=false;
					}
					else
					{
						countmultipleoverlap+=1;
					}
					packinfo.faces[i].constrictionsize=PlaneConstriction(packinfo,i,j);
					if (packinfo.faces[i].boundaryface==true)
					{
						break;
					}
					else
					{
						packinfo.faces[i].constrictioncentre=packinfo.position;
					}
				}
			}
		}
	}
	cout << "Boundary faces: "<< countboundary<<"\n";
	cout << "Overlap faces: "<< countoverlap<<"\n";
	cout << "Multiple overlap: "<<countmultipleoverlap<<"\n";
	TextConstriction(packinfo);
}

inline void Soil::CreateSequentialFace( DEM::SequentialPackingData&packinfo, int p0, int p1, int p2,int p3, bool sort)
{
	int p[3]={p0,p1,p2};
	if (sort)
	{
		int count;
		for (int i =0; i<2;i++)
			{
				for (int j=1; j<3; j++)
				{
					if (packinfo.specimen.Particles[p[i]]->Dmax< packinfo.specimen.Particles[p[j]]->Dmax)
						{
							count= p[i];
							p[i]=p[j];
							p[j]=count;
						}
				}
			}
	}
	DEM::SequentialFace facetemporary;

	Array<Vec3_t> V(3);
	for (int i=0;i<3;i++)
	{
		facetemporary.points[i]=p[i];
		V[i]= packinfo.specimen.Particles[p[i]]->x;
	}

	DEM::Face facedemtemprorary(V);															// add face in library
	facedemtemprorary.Normal(facetemporary.normal);
	facetemporary.normal /= norm(facetemporary.normal);
	facetemporary.point = p3;
	facetemporary.constrictionsize = ConstrictionSize(packinfo, facetemporary.points[0],facetemporary.points[1],facetemporary.points[2]);
	facetemporary.constrictioncentre = packinfo.position;
	facetemporary.boundaryface=false;
	double distance = dot(packinfo.specimen.Particles[p3]->x - packinfo.specimen.Particles[facetemporary.points[0]]->x, facetemporary.normal);
	if (distance>0)
	{
		facetemporary.faceuse = 1;
	}
	else
	{
		facetemporary.faceuse =-1;
	}
	if (packinfo.particlereuse)
	{
		packinfo.facereuse =false;
		for (int i=0; i< packinfo.numberfaces; i++)
		{
			if (facetemporary.points[0]== packinfo.faces[i].points[0])
			{
				if ((facetemporary.points[1]== packinfo.faces[i].points[1])and(facetemporary.points[2]== packinfo.faces[i].points[2]))
				{
					if (facetemporary.faceuse*packinfo.faces[i].faceuse <0)
					{
						CloseSequentialFace(packinfo,i);
					}
					else
					{
						packinfo.check =true;
					}
					packinfo.facereuse=true;
					packinfo.facereusenumber = i;
				}
			}
		}
	}
	if (packinfo.facereuse)
	{
		packinfo.facereuse = false;
	}
	else
	{
		packinfo.faces.push_back(facetemporary);												// add user-defined face
		packinfo.numberfaces +=1;
		packinfo.numberopenfaces +=1;
	}
}

inline void Soil::CreateSequentialTetrahedron(DEM::SequentialPackingData&packinfo, int face, int p3)
{
	DEM::SequentialTetrahedron temprorarytetrahedron;
	temprorarytetrahedron.points[0]=packinfo.faces[face].points[0];
	temprorarytetrahedron.points[1]=packinfo.faces[face].points[1];
	temprorarytetrahedron.points[2]=packinfo.faces[face].points[2];
	temprorarytetrahedron.points[3]=p3;
	temprorarytetrahedron.faces[0]= face;
	CreateSequentialFace(packinfo, packinfo.faces[face].points[0], packinfo.faces[face].points[1], p3, packinfo.faces[face].points[2]);
	if (packinfo.facereuse)
		{
			temprorarytetrahedron.faces[1]= packinfo.facereusenumber;
			packinfo.facereuse=false;
		}
	else
		{
			temprorarytetrahedron.faces[1]= packinfo.numberfaces-1;
		}
	CreateSequentialFace(packinfo,packinfo.faces[face].points[0], packinfo.faces[face].points[2], p3, packinfo.faces[face].points[1]);
	if (packinfo.facereuse)
		{
			temprorarytetrahedron.faces[2]= packinfo.facereusenumber;
			packinfo.facereuse=false;
		}
	else
		{
			temprorarytetrahedron.faces[2]= packinfo.numberfaces-1;
		}
	CreateSequentialFace( packinfo,packinfo.faces[face].points[1], packinfo.faces[face].points[2], p3, packinfo.faces[face].points[0]);
	if (packinfo.facereuse)
		{
			temprorarytetrahedron.faces[3]= packinfo.facereusenumber;
			packinfo.facereuse=false;
		}
	else
		{
			temprorarytetrahedron.faces[3]= packinfo.numberfaces-1;
		}
	packinfo.tetrahedra.push_back(temprorarytetrahedron);
}

inline void Soil::DeleteUnusedParticles(DEM::SequentialPackingData&packinfo)
{
	bool check=false;
	for (int i=0;i<packinfo.gsd.numberparticles;i++)
	{
		if (packinfo.particleuses[i]==0)
		{
			packinfo.specimen.Particles[i]->Tag =-1000;
			check =true;
		}
	}
	if (check)
	{
	    printf("\n%s--- Deleting unused particles --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
		Array <int> delpar;
		delpar.Push(-1000);
		packinfo.specimen.DelParticles(delpar);
	}
}

inline void Soil::DrawBoundary( DEM::SequentialPackingData&packinfo)
{
	// Build container
	if (packinfo.boundarytype=="File")
	{
	}
	else if (packinfo.boundarytype=="Sphere")
	{
		packinfo.specimen.AddSphere(-1000,Vec3_t(0.,0.,0.),packinfo.boundaryfactors[0],packinfo.gsd.specialgravity);
	}
	else if (packinfo.boundarytype=="Cube")
	{
		packinfo.specimen.AddCube(-1001,Vec3_t(0.,0.,0.),packinfo.gsd.inters[0].begindiameter/2,packinfo.boundaryfactors[0],packinfo.gsd.specialgravity);
	}
	else if (packinfo.boundarytype=="Cylinder")
	{
		packinfo.specimen.AddCylinder(-1002,Vec3_t(0.,0.,-packinfo.boundaryfactors[1]/2),packinfo.boundaryfactors[0],Vec3_t(0.,0.,packinfo.boundaryfactors[1]/2),packinfo.boundaryfactors[0],packinfo.gsd.inters[0].begindiameter/2,packinfo.gsd.specialgravity);
		packinfo.specimen.AddPlane(-1003,Vec3_t(0.,0.,packinfo.boundaryfactors[1]/2),packinfo.gsd.inters[0].begindiameter/2,1.2*packinfo.boundaryfactors[0],1.2*packinfo.boundaryfactors[0],packinfo.gsd.density,1.);
		packinfo.specimen.AddPlane(-1004,Vec3_t(0.,0.,-packinfo.boundaryfactors[1]/2),packinfo.gsd.inters[0].begindiameter/2,1.2*packinfo.boundaryfactors[0],1.2*packinfo.boundaryfactors[0],packinfo.gsd.density,1.);
	}
	else if (packinfo.boundarytype=="Box")
	{
		packinfo.specimen.GenBox(-1005, packinfo.boundaryfactors[0],packinfo.boundaryfactors[1],packinfo.boundaryfactors[2],packinfo.gsd.inters[0].begindiameter/2,1.2,false);
	}
	for (int i =0; i<packinfo.specimen.Particles.Size(); i++)
	{
		packinfo.boundary.push_back(i);
	}
}

inline void Soil::DropDown(DEM::SequentialPackingData&packinfo)
{
	// should have a cylinder boundary
	if (packinfo.para.starttag<0)
		{
			packinfo.specimen.GenBoundingBox(packinfo.para.starttag,0.01,1.2,false);
			for (int i=packinfo.para.starttag;i>packinfo.para.starttag-6;i--)
			    {
			        packinfo.specimen.GetParticle(i)->FixVeloc();
			    }
		}
    Dict P;
    for (int i=0;i<packinfo.para.numberintervals;i++)
    {
    	P.Set(/*Tag*/i*packinfo.para.numbershapes,"Kn Kt Gn Gt Gv Mu",packinfo.para.Kn,packinfo.para.Kt,packinfo.para.Gn,packinfo.para.Gt,packinfo.para.Gv, packinfo.para.Mu);
        for (int j=1; j<packinfo.para.numbershapes; j++)
        	{
        		P.Set(/*Tag*/-i*packinfo.para.numbershapes-j,"Kn Kt Gn Gt Gv Mu",packinfo.para.Kn,packinfo.para.Kt,packinfo.para.Gn,packinfo.para.Gt, packinfo.para.Gv, packinfo.para.Mu);
        	}
    }
    packinfo.specimen.SetProps(P);
    Vec3_t g(0.0,0.0,-9.8);
    //packinfo.specimen.Alpha = packinfo.para.Alpha;
    for (size_t i=0;i<packinfo.specimen.Particles.Size();i++)
    {
    	packinfo.specimen.Particles[i]->Ff = packinfo.specimen.Particles[i]->Props.m*g;
    	//posbefore.Push(particles.Particles[i]->x);
    }
    packinfo.specimen.Solve  (/*tf*/packinfo.para.tf, /*dt*/packinfo.para.dt, /*dtOut*/packinfo.para.dtout, NULL, NULL, /*filekey*/packinfo.para.FileKey.c_str(),/*Visit visualization*/packinfo.para.visualization,/*N_proc*/packinfo.para.numberprocessors, /*kinematic energy*/packinfo.para.Kinematicenergy);
}

inline void Soil::EstablishLocalSystem( DEM::SequentialPackingData&packinfo, int face)
{
	packinfo.localsystem[0]= packinfo.specimen.Particles[packinfo.faces[face].points[1]]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x;
	packinfo.localsystem[0]= packinfo.localsystem[0]/norm(packinfo.localsystem[0]);							// e1 vector
	packinfo.localsystem[2]= packinfo.faces[face].normal;													// e3 vector
	packinfo.localsystem[1]= -cross(packinfo.localsystem[0],packinfo.localsystem[2]);
	packinfo.localsystem[1]= packinfo.localsystem[1]/norm(packinfo.localsystem[1]);							// may be not in need
	packinfo.localroot = packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x;
}

inline void Soil::FindMinimumSpecimenSize(DEM::SequentialPackingData&packinfo, int coarsest)
{

	double coarsestdiameter = pow(10,(log10(packinfo.gsd.inters[packinfo.gsd.numberinters-1].begindiameter)+log10(packinfo.gsd.inters[packinfo.gsd.numberinters-1].enddiameter))/2);
	// should be aware of different shape
	double coarsestshapemass = ShapeMass(packinfo, 0, coarsestdiameter,packinfo.gsd.inters[packinfo.gsd.numberinters-1].specialgravity);; // mass of one sphere
	double scale = (coarsest* coarsestshapemass*(1+packinfo.approximation))/(packinfo.gsd.inters[packinfo.gsd.numberinters-1].end-packinfo.gsd.inters[packinfo.gsd.numberinters-1].begin)*100/packinfo.gsd.mass;
	packinfo.gsd.volume *= scale;
	for (int i=0;i<packinfo.gsd.numberboundaries;i++)						// should be aware of extended size
	{
		packinfo.boundarysizes[i]*=pow(scale,1/3.0);
		packinfo.boundaryfactors[i]*=pow(scale,1/3.0);
	}
	packinfo.gsd.mass *=scale;
	for (int i=0; i<packinfo.gsd.numberinters;i++)
	{
		packinfo.gsd.inters[i].mass*=scale;
	}
}

inline void Soil::FrozenTag(DEM::SequentialPackingData&packinfo, int tag=-1000)
{
	for (int i=0;i<packinfo.specimen.Particles.Size();i++)
	{
		if (packinfo.specimen.Particles[i]->Tag <= tag)
		{
			packinfo.specimen.Particles[i]->FixVeloc();
		}
	}
}

inline void Soil::MoveOn(DEM::SequentialPackingData&packinfo, int p3, int face)
{
	Vec3_t distance = packinfo.specimen.Particles[p3]->x-packinfo.specimen.Particles[packinfo.overlappingpoint]->x;
	Vec3_t projection = dot(distance,packinfo.faces[face].normal);
	double move = sqrt(pow(packinfo.specimen.Particles[p3]->Dmax+packinfo.specimen.Particles[packinfo.overlappingpoint]->Dmax,2.0)-pow(norm(distance),2.0)+pow(norm(projection),2.0))-norm(projection);
	Vec3_t finalposition = packinfo.specimen.Particles[p3]->x-packinfo.faces[face].faceuse*move*packinfo.faces[face].normal;
	packinfo.specimen.Particles[p3]->Position(finalposition);
}

inline void Soil::ParticleInfo(DEM::SequentialPackingData&packinfo,int p0)
{
	cout<< "Particle: "<<p0<<" Dmax: "<<packinfo.specimen.Particles[p0]->Dmax<<"\n";
	cout<< " 	x: "<<packinfo.specimen.Particles[p0]->x<< "\n";
}

inline void Soil::PlaceSequentialParticle(DEM::SequentialPackingData&packinfo,int face,int p3)
{
	packinfo.temproraryparticle=p3;
	PutSequentialParticle(packinfo,face,packinfo.temproraryparticle);
	CheckOverlap(packinfo, packinfo.temproraryparticle);
	CheckBoundary(packinfo,packinfo.temproraryparticle);
	if ((!packinfo.checkoverlap)and (packinfo.checkboundary)and(packinfo.checkradius))			// if overlap only
	{
		double distanceoverlap=0.;
		double distancepoint = dot(packinfo.specimen.Particles[packinfo.faces[face].point]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.faces[face].normal);
		while ((!packinfo.checkoverlap)and(packinfo.checkboundary))
		{
			distanceoverlap=dot(packinfo.specimen.Particles[packinfo.overlappingpoint]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.faces[face].normal);
			if (distanceoverlap*distancepoint <0)
			{
				packinfo.temproraryparticle =packinfo.overlappingpoint;
				packinfo.particlereuse =true;
				packinfo.checkoverlap=true;
			}
			else
			{
				MoveOn(packinfo,packinfo.temproraryparticle,face);
				CheckOverlap(packinfo,packinfo.temproraryparticle);
				CheckBoundary(packinfo,packinfo.temproraryparticle);
			}
		}
	}
}

inline void Soil::PrintOut(DEM::SequentialPackingData&packinfo)
{
	cout << packinfo.gsd.numberinters<<" "<<packinfo.gsd.numberparticles<< " "<<packinfo.numberunusedparticles<<"\n";
	cout << packinfo.gsd.mass<<" "<<packinfo.gsd.volume<<" "<<packinfo.boundarysizes[0]<<"\n";
	for (int i=0;i<packinfo.gsd.numberinters;i++)
	{
		cout << " "<<i<<" "<<packinfo.gsd.inters[i].ability<<" "<<packinfo.gsd.inters[i].numberparticles<<"\n";
	}
}

inline void Soil::PrepareList(DEM::SequentialPackingData&packinfo, bool randomness)
{
    printf("\n%s--- Preparing List of particles --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
	double shapem;
	double shape;
	double mass;
	double passingmass=0;
	packinfo.gsd.numberparticles =0;
	if (randomness)
	{
		// cout << "random\n";	// check whether random packing
		for (int i=packinfo.gsd.numberinters-1; i>-1;i--)
		{
			packinfo.gsd.inters[i].numberparticles=0;
			packinfo.gsd.inters[i].mass += passingmass;
		}
	}
	else
	{
		for (int i =packinfo.gsd.numberinters-1; i>-1;i--)
		{
			packinfo.gsd.inters[i].numberparticles=0;
			packinfo.gsd.inters[i].diameter=pow(10,(log10(packinfo.gsd.inters[i].begindiameter)+log10(packinfo.gsd.inters[i].enddiameter))/2);
			packinfo.gsd.inters[i].mass += passingmass;
			shape=0;
			for (int j=0; j<packinfo.gsd.numbershapes;j++)
			{
				shape+=packinfo.gsd.inters[i].shaperatios[j];
			}
			if (shape <=0)
			{
				throw new Fatal("Error: No proportion for particles in interval from <%s> % to <%d>%", packinfo.gsd.inters[i].begin, packinfo.gsd.inters[i].end);
			}
			else
			{
				for (int j=0;j<packinfo.gsd.numbershapes;j++)
				{
					packinfo.gsd.inters[i].shaperatios[j] /= shape;
				}
			}
			mass=0.;
			for (int j=0;j<packinfo.gsd.numbershapes;j++)
			{
				if (packinfo.gsd.inters[i].shaperatios[j]>0)
				{

					shapem= ShapeMass(packinfo, j, packinfo.gsd.inters[i].diameter,packinfo.gsd.inters[i].specialgravity);
					packinfo.gsd.inters[i].shapestates[j][0] = floor(packinfo.gsd.inters[i].shaperatios[j]* packinfo.gsd.inters[i].mass/shapem);
					packinfo.gsd.inters[i].numberparticles += packinfo.gsd.inters[i].shapestates[j][0];
					mass += packinfo.gsd.inters[i].shapestates[j][0]*shapem;
				}
				else
				{
					packinfo.gsd.inters[i].shapestates[j][0]=0;
				}
			}
			passingmass = packinfo.gsd.inters[i].mass -mass;
			packinfo.gsd.inters[i].mass = mass;
			if (packinfo.gsd.inters[i].numberparticles>0)
			{
				packinfo.gsd.inters[i].ability=true;
			}
			else
			{
				packinfo.gsd.inters[i].ability=false;
			}
		}
		// Add particles
		size_t count =0;
		for (int i=0; i<packinfo.gsd.numberinters;i++)
		{
			if (packinfo.gsd.inters[i].ability)
			{
				packinfo.gsd.inters[i].firstparticle=count;
				for (int j=0; j<packinfo.gsd.numbershapes; j++)
				{
					if (packinfo.gsd.inters[i].shapestates[j][0]>0)
					{
						packinfo.gsd.inters[i].shapestates[j][1]=count;
						packinfo.gsd.inters[i].shapestates[j][3]=count;
						for (int k=0; k<packinfo.gsd.inters[i].shapestates[j][0];k++)
						{
							AddParticle(packinfo, i, j, packinfo.gsd.inters[i].diameter,packinfo.gsd.inters[i].specialgravity,OrthoSys::O,packinfo.gsd.inters[0].begindiameter/2);
							packinfo.particleuses.push_back(false);
							count+=1;
						}
						packinfo.gsd.inters[i].shapestates[j][2]=count-1;
					}
					else
					{
						for (int k=2; k<4;k++)
						{
							packinfo.gsd.inters[i].shapestates[j][k]=-1;
						}
					}
				}
				packinfo.gsd.inters[i].lastparticle=count-1;
				packinfo.gsd.inters[i].usingparticle=packinfo.gsd.inters[i].lastparticle;
			}
		}
		packinfo.gsd.numberparticles = packinfo.specimen.Particles.Size();
	}
	//PrintOut(packinfo);			//Test list prepare process
}

inline void Soil::PutSequentialParticle(DEM::SequentialPackingData&packinfo, int face, int p3)
{
	EstablishLocalSystem(packinfo,face);															// calculate local coordinates system
	packinfo.checkradius =false;
	double r[4];
	for (int i =0;i<3;i++)
		{
			r[i] = packinfo.specimen.Particles[packinfo.faces[face].points[i]]->Dmax;
		}
	r[3]=packinfo.specimen.Particles[p3]->Dmax;
	double x2 = norm(packinfo.specimen.Particles[packinfo.faces[face].points[1]]->x -packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x);
	double a1 = (r[0] -r[1])/x2;
	double b1 = (pow(x2,2.)+pow(r[0],2.)-pow(r[1],2.))/2./x2;
	double x3 = dot(packinfo.specimen.Particles[packinfo.faces[face].points[2]]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.localsystem[0]);
	double y3 = dot(packinfo.specimen.Particles[packinfo.faces[face].points[2]]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x, packinfo.localsystem[1]);
	double a2 = (r[0]-r[2]-a1*x3)/y3;
	double b2 = (pow(x3,2.)+pow(y3,2.)+pow(r[0],2.)-pow(r[2],2.)-2.*b1*x3)/2./y3;
	double x4 = a1*r[3]+b1;
	double y4 = a2*r[3]+b2;
	double z4 = pow(r[0]+r[3],2)-pow(a1*r[3]+b1,2)-pow(a2*r[3]+b2,2);
	if (z4 > 0)
		{
			z4=-packinfo.faces[face].faceuse*pow(z4,0.5);
			Vec3_t position(x4,y4,z4);
			Vec3_t finalposition =  packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x;
			for (int i=0; i<3; i++)
				{
					for (int j=0; j<3; j++)
						{
							finalposition(i)+= position(j)*packinfo.localsystem[j](i);
						}
				}
			packinfo.specimen.Particles[p3]->Position(finalposition);
			packinfo.checkradius =true;
		}
}

inline void Soil::ReadDemParameters(DEM::SequentialPackingData&packinfo, string FileKey)
{
	string filename = FileKey;
	filename.append(".par");
	if (!Util::FileExists(filename))
	{
		string f="../src/";
		f.append(filename.c_str());
		filename =f;
		if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.c_str());
	}
    ifstream datain;
    datain.open(filename.c_str());
    printf("\n%s--- Loading DEM parameters file %s --------------------------------------------%s\n",TERM_CLR1,filename.c_str(),TERM_RST);
	datain >> packinfo.para.Kn;					datain.ignore(200,'\n');
	datain >> packinfo.para.Kt;					datain.ignore(200,'\n');
	datain >> packinfo.para.Gn;					datain.ignore(200,'\n');
	datain >> packinfo.para.Gt;					datain.ignore(200,'\n');
	datain >> packinfo.para.Gv;					datain.ignore(200,'\n');
	datain >> packinfo.para.Mu;					datain.ignore(200,'\n');
	datain >> packinfo.para.dt;					datain.ignore(200,'\n');
	datain >> packinfo.para.tf;					datain.ignore(200,'\n');
	datain >> packinfo.para.dtout;				datain.ignore(200,'\n');
	datain >> packinfo.para.FileKey;			datain.ignore(200,'\n');
	datain >> packinfo.para.visualization;		datain.ignore(200,'\n');
	datain >> packinfo.para.numberprocessors;	datain.ignore(200,'\n');
	datain >> packinfo.para.roundratio;			datain.ignore(200,'\n');
	datain >> packinfo.para.Alpha;     			datain.ignore(200,'\n');
	datain >> packinfo.para.Kinematicenergy;  	datain.ignore(200,'\n');
	datain >> packinfo.para.starttag;  			datain.ignore(200,'\n');
	datain >> packinfo.para.numbershapes;		datain.ignore(200,'\n');
	datain >> packinfo.para.numberintervals; 	datain.ignore(200,'\n');
	datain.close();
}

inline void Soil::ReadGsd( DEM::SequentialPackingData&packinfo, string FileKey)
{
	string filename=FileKey;
	if ((filename.substr(filename.length()-4).compare("gsd")==0)or(filename.length()<4))
	{
		filename.append(".gsd");
	}
	if (!Util::FileExists(filename))
	{
		string f="../src/";
		f.append(filename.c_str());
		filename =f;
		if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.c_str());
	}
    ifstream datain;
    datain.open(filename.c_str());
    printf("\n%s--- Loading GSD file %s --------------------------------------------%s\n",TERM_CLR1,filename.c_str(),TERM_RST);
    datain >> packinfo.gsd.Soilname;						datain.ignore(200,'\n');
    datain >> packinfo.boundarytype;						datain.ignore(200,'\n');
    if (packinfo.boundarytype =="File")
    {
    	datain >>packinfo.gsd.boundaryfile;					datain.ignore(200,'\n');
    	packinfo.gsd.boundaryfile.append(".h5");
    	if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",packinfo.gsd.boundaryfile.c_str());
    	// how to check overlap with boundary object ?
    	packinfo.specimen.Load(packinfo.gsd.boundaryfile.c_str());
    }
    else
    {
        datain >> packinfo.gsd.numberboundaries;
    	for (int i=0;i<packinfo.gsd.numberboundaries;i++)
    	{
    	   	datain >> packinfo.boundarysizes[i];
    	}
    	for (int i=0; i<packinfo.gsd.numberboundaries;i++)
    	{
    		datain >> packinfo.boundaryfactors[i];
    	}
    }
    												datain.ignore(200,'\n');
    datain >> packinfo.gsd.arrangement;				datain.ignore(200,'\n');		// type of grain arrangement
//	datain >> packinfo.gsd.particleshape;			datain.ignore(200,'\n');		// type of particle shape
    datain >> packinfo.numberprocessors; 			datain.ignore(200,'\n');
    datain >> packinfo.gsd.specialgravity;			datain.ignore(200,'\n');
    datain >> packinfo.gsd.porosity;				datain.ignore(200,'\n');
    datain >> packinfo.gsd.maxfraction;				datain.ignore(200,'\n');		// Maximum fraction of a divided interval
    datain >> packinfo.gsd.numberintervals;
    datain >> packinfo.gsd.numbershapes;			datain.ignore(200,'\n');
    DEM::Interval inter;
    for (int i=0; i<packinfo.gsd.numberintervals; i++)
    {
    	datain >> inter.begin;
    	datain >> inter.begindiameter;
    	datain >> inter.specialgravity;
    	for (int j=0; j <packinfo.gsd.numbershapes; j++)
    	{
    		datain >> inter.shaperatios[j];
    	}
    												datain.ignore(200,'\n');
    	packinfo.gsd.intervals.push_back(inter);
    }
    datain >> packinfo.approximation;
    for (int i =0; i<3; i++)
    {
    	datain >> packinfo.gsd.rectangularboxratios[i];								// rectangular size ratio
    }
    												datain.ignore(200,'\n');
    datain.close();
    //cout << packinfo.boundarytype<<"\n"; //test reading process
}

inline void Soil::SaveDomain(DEM::Domain&specimen,string filename, int outputtype)
{
	if (outputtype==1)
	{
		specimen.WriteXDMF(filename.c_str());
	}
	specimen.Save(filename.c_str());
}

inline void Soil::SaveTetrahedraMesh(DEM::SequentialPackingData&packinfo)
{
	ofstream dataout;
	String fn;
	fn.Printf("%scentres.%s",packinfo.gsd.Soilname.c_str(),"m");
	dataout.open(fn.c_str());
	for (int i=0; i<packinfo.specimen.Particles.Size();i++)
	{
		for (int j=0; j<3; j++)
		{
			dataout<<packinfo.specimen.Particles[i]->x(j)<<"\n";
		}
	}
	dataout.close();
	fn.Printf("%stetra.%s",packinfo.gsd.Soilname.c_str(),"m");
	dataout.open(fn.c_str());
	for (int i=0; i<packinfo.tetrahedra.size();i++)
	{
		for (int j=0; j<4; j++)
		{
			dataout<< packinfo.tetrahedra[i].points[j]<<"\n";
		}
	}
	dataout.close();
}

inline void Soil::SequentialPacking(DEM::SequentialPackingData&packinfo)
{
//	for (int i=0; i<packinfo.numberprocessors; i++)
//	{
//		packinfo.threadnumbers.push_back(int(i/packinfo.numberprocessors*specimen.Particles.Size()));
//	}
//	packinfo.threadnumbers.push_back(specimen.Particles.Size());
	packinfo.particlereuse =false;
	packinfo.facereuse =false;
	packinfo.check=true;
	packinfo.usinginter = packinfo.gsd.numberinters-1;
	packinfo.usingparticle = packinfo.specimen.Particles.Size()-1;
	int p[4]= {packinfo.usingparticle, packinfo.usingparticle-1, packinfo.usingparticle-2, packinfo.usingparticle-3};
	BasicTetrahedron(packinfo,p);	// create initial tetrahedron
	packinfo.usingparticle -=4;
	// packing loop
    printf("\n%s--- Sequentially packing --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
	while ((packinfo.numberopenfaces >0)and(packinfo.numberunusedparticles >0))
	{
//		if (packinfo.numberunusedparticles%1000==0)
//		{
//			cout << packinfo.numberunusedparticles<<"\n";
//		}
		while (packinfo.particleuses[packinfo.usingparticle])
			{
				packinfo.usingparticle -=1;
				if (packinfo.usingparticle==0)
					{
						break;
					}
				}
		while (packinfo.faces[packinfo.usingface].faceuse == 0)
			{
				packinfo.usingface +=1;
				if (packinfo.usingface == packinfo.numberfaces)
					{
						break;
					}
			}
//		if (packinfo.randomness>0)
//		{
//			packinfo.temproraryparticle=
//		}
//		else
//		{
//			packinfo.temproraryparticle= packinfo.usingparticle;
//		}
		packinfo.temproraryparticle= packinfo.usingparticle;
		if (packinfo.gsd.arrangement=="Layer")
		{
			PlaceSequentialParticle(packinfo, packinfo.usingface,packinfo.temproraryparticle);
		}
		else if (packinfo.gsd.arrangement=="Discrete")
		{
			TrySequentialParticle( packinfo, packinfo.usingface, packinfo.temproraryparticle);
		}
//		if (packinfo.gsd.arrangement=="Layer")
//		{
//			PlaceSequentialParticle(packinfo, packinfo.usingface,packinfo.usingparticle);
//		}
//		else if (packinfo.gsd.arrangement=="Discrete")
//		{
//			TrySequentialParticle( packinfo, packinfo.usingface, packinfo.usingparticle);
//		}
		if (packinfo.checkradius and packinfo.checkoverlap and packinfo.checkboundary)
		{
			CreateSequentialTetrahedron(packinfo,packinfo.usingface,packinfo.temproraryparticle);
			if (packinfo.particleuses[packinfo.temproraryparticle])
			{
				packinfo.particlereuse= false;
			}
			else
			{
				UseParticle(packinfo,packinfo.temproraryparticle);
			}
		}
		CloseSequentialFace(packinfo, packinfo.usingface);
	}


}

inline void Soil::TextConstriction(DEM::SequentialPackingData&packinfo, int type)
{
	ofstream dataout;
	String fn;
	fn.Printf("%s.%s",packinfo.gsd.Soilname.c_str(),"csd");
	dataout.open(fn.c_str());
	for (int i=0; i<packinfo.numberfaces;i++)
	{
		if ((packinfo.faces[i].constrictionsize>packinfo.approximation)and(packinfo.faces[i].boundaryface==false)) // type can be used to work with more options
		{
			dataout<< 2*packinfo.faces[i].constrictionsize<<"\n"; // export in diameter
		}
	}
	dataout.close();
}

inline void Soil::TextOut(DEM::SequentialPackingData&packinfo)
{
	ofstream dataout;
	String fn;
	fn.Printf("%s.%s",packinfo.gsd.Soilname.c_str(),"out");
	dataout.open(fn.c_str());
	dataout<< packinfo.gsd.numberparticles<< " number of particles \n";
	dataout<< packinfo.numberunusedparticles<<" number of unused particles \n";
	dataout<< packinfo.numberfaces<<" number of faces \n";
	dataout<< packinfo.numberopenfaces<< " number of open faces \n";
	dataout<< packinfo.boundarytype<<" "<<packinfo.gsd.volume<<" "<<packinfo.boundarysizes[0]<<"\n";
	for (int i=0;i<packinfo.gsd.numberinters;i++)
	{
		dataout<< " "<<i<<" "<<packinfo.gsd.inters[i].diameter<<" "<<packinfo.gsd.inters[i].numberparticles<<" "<<packinfo.gsd.inters[i].mass<<" "<<packinfo.gsd.inters[i].firstparticle<< " "<<packinfo.gsd.inters[i].lastparticle<<"\n";
	}
	dataout<<"List of unused particles \n";
	for (int i=0; i < packinfo.gsd.numberinters;i++)
	{
		if (packinfo.gsd.inters[i].ability)
		{
			dataout<< " "<<i<<" : "<<packinfo.gsd.inters[i].usingparticle<< " "<<packinfo.gsd.inters[i].firstparticle<<"\n";
		}
	}
	dataout.close();

}

inline void Soil::TrySequentialParticle( DEM::SequentialPackingData&packinfo, int face, int p3)
{
	packinfo.temproraryparticle = p3;
	PutSequentialParticle(packinfo,face,packinfo.temproraryparticle);
	CheckOverlap( packinfo, packinfo.temproraryparticle);
	CheckBoundary(packinfo,packinfo.temproraryparticle);
	if ((!packinfo.checkradius)or (!packinfo.checkboundary) or (!packinfo.checkoverlap))
	{
		double distancepoint =0;
		double distanceoverlap =0;
		int interval = GetInterval(packinfo, packinfo.temproraryparticle);
		distancepoint = dot(packinfo.specimen.Particles[packinfo.faces[face].point]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.faces[face].normal);
		if (packinfo.specimen.Particles[packinfo.temproraryparticle]->Dmax >packinfo.faces[face].constrictionsize)					// if particle size is bigger than constriction size
		{
			while ((!packinfo.checkboundary)or(!packinfo.checkradius)or(!packinfo.checkoverlap))									// try to reduce size of particles
			{
				interval-=1;
				if ((interval<0)or(packinfo.gsd.inters[interval].diameter<packinfo.faces[face].constrictionsize))					// if random list, so should be usingparticle diameter
				{
					if (interval>0)								// if size < constriction size, must +1 to -1 later
					{
						interval+=1;
					}
					break;
				}
				if (packinfo.gsd.inters[interval].ability)
				{
					packinfo.temproraryparticle=packinfo.gsd.inters[interval].usingparticle;
					PutSequentialParticle(packinfo,face,packinfo.temproraryparticle);
					CheckOverlap(packinfo,packinfo.temproraryparticle);
					CheckBoundary(packinfo,packinfo.temproraryparticle);
					if ((packinfo.checkboundary)and(packinfo.checkradius)and(!packinfo.checkoverlap))								// if overlap only
					{
						distanceoverlap=dot(packinfo.specimen.Particles[packinfo.overlappingpoint]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.faces[face].normal);
						while ((distancepoint*distanceoverlap>0)and (!packinfo.checkoverlap)and(packinfo.checkboundary))
						{
							MoveOn(packinfo,packinfo.temproraryparticle,face);
							CheckOverlap(packinfo,packinfo.temproraryparticle);
							CheckBoundary(packinfo,packinfo.temproraryparticle);
							if (!packinfo.checkoverlap)
							{
								distanceoverlap = dot(packinfo.specimen.Particles[packinfo.overlappingpoint]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.faces[face].normal);
							}
						}
					}
				}
			}
		}
		if ((!packinfo.checkradius)or (!packinfo.checkboundary)or(!packinfo.checkoverlap))						// if particle size is smaller than constriction size
		{
			packinfo.checkradius =true;
			while ((!packinfo.checkboundary)or(!packinfo.checkoverlap))
			{
				interval -=1;
				if (interval <0)
				{
					break;
				}
				if (packinfo.gsd.inters[interval].ability)
				{
					packinfo.temproraryparticle = packinfo.gsd.inters[interval].usingparticle;
					packinfo.specimen.Particles[packinfo.temproraryparticle]->Position(packinfo.faces[face].constrictioncentre-packinfo.approximation*packinfo.faces[face].normal*packinfo.faces[face].faceuse); // Try to put at the centre of the constriction. Can put next to two other particles
					CheckBoundary(packinfo,packinfo.temproraryparticle);
					CheckOverlap(packinfo,packinfo.temproraryparticle);
					if ((packinfo.checkboundary)and(!packinfo.checkoverlap))								// if overlap only
					{
						distanceoverlap=dot(packinfo.specimen.Particles[packinfo.overlappingpoint]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.faces[face].normal);
						while ((distancepoint*distanceoverlap>0)and (!packinfo.checkoverlap)and(packinfo.checkboundary))
						{
							MoveOn(packinfo,packinfo.temproraryparticle,face);
							CheckOverlap(packinfo,packinfo.temproraryparticle);
							CheckBoundary(packinfo,packinfo.temproraryparticle);
							if (!packinfo.checkoverlap)
							{
								distanceoverlap = dot(packinfo.specimen.Particles[packinfo.overlappingpoint]->x - packinfo.specimen.Particles[packinfo.faces[face].points[0]]->x,packinfo.faces[face].normal);
							}
						}
					}
				}
			}
		}

	}
}

inline void Soil::UseParticle( DEM::SequentialPackingData&packinfo, int p3)
{
	if (packinfo.particleuses[p3])
	{
		cout << "Particle " << p3 << " has been used before \n";
		return;
	}
	packinfo.particleuses[p3]=true;
	packinfo.numberunusedparticles-=1;
	int inter =GetInterval(packinfo, p3);
	packinfo.gsd.inters[inter].usingparticle -=1;
	if (packinfo.gsd.inters[inter].usingparticle < packinfo.gsd.inters[inter].firstparticle )
		{
			packinfo.gsd.inters[inter].ability =false;
			if (inter == packinfo.usinginter)
				{
					packinfo.usinginter -=1;
				}
		}
}
}

