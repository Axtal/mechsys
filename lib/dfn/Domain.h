/************************************************************************
            Array<Vec3_t > Verts_temp;
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

#ifndef MECHSYS_DFN_DOMAIN_H
#define MECHSYS_DFN_DOMAIN_H

// Hdf5
#include <hdf5.h>
#include <hdf5_hl.h>

//STL
#include <vector>
#include <algorithm>
#include <utility>

// MechSys
#include <mechsys/dfn/Fracture.h>
#include <mechsys/util/tree.h>
//#include <mechsys/dfn/Graph_i.h>
#include <mechsys/dfn/Graph_x.h>

// Boost
//#include <boost/config.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/connected_components.hpp>

namespace DFN
{

	class Domain
	{
	public:
		//Data
		Array<Fracture> Fractures;													  ///< Array of Fractures, this array stores all generated fractures
		Array<size_t> Connections;														  ///< Array of Fracture connectivity, record the tag / ID of fractures that are connected
		Array<Array<size_t>> Listofclusters;											  ///< List of fractures per cluster
		Array<Array<size_t>> Percolation_cluster;									  ///< three orders, N dimensions; first order means x direction, second means y, third means z; alone a percolation direction, there might be zero or more clusters;
		std::map<std::pair<size_t, size_t>, std::pair<Vec3_t, Vec3_t>> Intersections; ///< Map of line intersection between pairs of fractures
		Vec6_t Model_domain;														  ///< Top, bottom, front, back, left, right
		double n_I;																	  ///< Average number of intersections per fracture
		double P30;
		double P32_total;
		double P32_connected;
		double Percolation_parameter;
		double Ratio_of_P32; ///< probability of a fracture belonging to a percolation cluster
		double Expected_area_times_perimeter;

		Array<Fracture> Surfaces;														///< model surface
		Array<size_t> Connections_S;													///< thos fractures intersect with surfaces
		std::map<std::pair<size_t, size_t>, std::pair<Vec3_t, Vec3_t>> Intersections_S; ///< Map of line intersection between fractures and surfaces

		// Methods
		void Model_set(const double model_size[6]);			///< define model domain
		void AddSquareFracture(size_t Tag, Fracture &c);	///< Adding a square fracture
		bool Intersect_A(Fracture F1, Fracture F2);			///< JUST Identify if fractures are connected, but no intersection points are returned
		void Modify_fracture_attributes_Zmax(Fracture &F2); ///< modify fracture F2 (its Nvertices, vertexes and area), because F2 intersects surface(s)
		void Modify_fracture_attributes_Zmin(Fracture &F2);
		void Modify_fracture_attributes_Ymin(Fracture &F2);
		void Modify_fracture_attributes_Ymax(Fracture &F2);
		void Modify_fracture_attributes_Xmin(Fracture &F2);
		void Modify_fracture_attributes_Xmax(Fracture &F2);
		void Extract_intersection_between_surfaces_and_Fractures(Fracture &F2);
		bool Intersect(Fracture F1, Fracture F2); ///< Function to check if two fractures intersect, and return intersection
		void Clusters();						  ///< Check the connectivity array to form the clusters
		void Create_whole_model(const size_t n, const double random_seed, const double model_size[6], const String str, const double array11[3][2], const double array12[4], const double array13[7]);
		size_t Identify_percolation_clusters(String str);
		void WriteFrac(char const *FileKey);											///< Writing the Fracture network in h5 format
		void PlotMatlab_DFN(char const *FileKey);										///< Plot to Matlab script
		void PlotMatlab_DFN_and_Intersection(char const *FileKey);						///< Plot DFN and intersections
		void PlotMatlab_ORI_SCATTER(char const *FileKey);								///< Plot rose diagram
		void PlotMatlab_Traces_on_Model_surfaces(char const *FileKey);					///< Plot traces on surface
		void PlotMatlab_DFN_Highlight_Cluster(char const *FileKey);						///< Plot DFN highlighted by cluster values
		void PLotMatlab_DFN_Cluster_along_a_direction(char const *FileKey, string str); ///< Plot percolation cluster spanning model along x, y or z axis
		void Connectivity_uniform_orientation(string str);								///< When orientation data have uniform distribution, determines percolation-related varibales
		void Connectivity_fisher_orientation(string str);								///< not finished
		void Average_number_of_intersections_per_fracture();							///< Average_number_of_intersections_per_fracture
		void PlotMatlab_Radius_and_Area_kstest(char const *FileKey);					///< Use kstest tests if two groups of data are both following the same distribution, but seems not work
		void PlotMatlab_Radius_and_Perimeter_kstest(char const *FileKey);
		void DataFile_Radius_AreaAndPerimeter(char const *FileKey);				///< outputs the data
		void Expected_area_and_perimeter(double x0, double x1, double alpha_g); ///<  <A*P>
	};

	inline void Domain::AddSquareFracture(size_t Tag, Fracture &c)
	{

		if (Model_domain(4) <= c.Center(0) &&
			Model_domain(5) >= c.Center(0) &&
			Model_domain(2) <= c.Center(1) &&
			Model_domain(3) >= c.Center(1) &&
			Model_domain(1) <= c.Center(2) &&
			Model_domain(0) >= c.Center(2))
		{

			bool y1 = Intersect_A(Surfaces[0], c);
			bool y2 = Intersect_A(Surfaces[1], c);
			bool y3 = Intersect_A(Surfaces[2], c);
			bool y4 = Intersect_A(Surfaces[3], c);
			bool y5 = Intersect_A(Surfaces[4], c);
			bool y6 = Intersect_A(Surfaces[5], c);

			if (y1 == 1 || y2 == 1 || y3 == 1 || y4 == 1 || y5 == 1 || y6 == 1)
			{
				if (y1 == 1)
				{
					c.If_intersect_surfaces(0) = 1;
					//Modify_fracture_attributes_Zmax(c);
				}
				if (y2 == 1)
				{
					c.If_intersect_surfaces(1) = 1;
					//Modify_fracture_attributes_Zmin(c);
				}
				if (y3 == 1)
				{
					c.If_intersect_surfaces(2) = 1;
					//Modify_fracture_attributes_Ymin(c);
				}
				if (y4 == 1)
				{
					c.If_intersect_surfaces(3) = 1;
					//Modify_fracture_attributes_Ymax(c);
				}
				if (y5 == 1)
				{
					c.If_intersect_surfaces(4) = 1;
					//Modify_fracture_attributes_Xmin(c);
				}
				if (y6 == 1)
				{
					c.If_intersect_surfaces(5) = 1;
					//Modify_fracture_attributes_Xmax(c);
				}
			}
			Fractures.Push(c);
			Fractures[Fractures.Size() - 1].Tag = Fractures.Size() - 1;
			return;
		}
		else
		{

			bool y1 = Intersect_A(Surfaces[0], c);
			bool y2 = Intersect_A(Surfaces[1], c);
			bool y3 = Intersect_A(Surfaces[2], c);
			bool y4 = Intersect_A(Surfaces[3], c);
			bool y5 = Intersect_A(Surfaces[4], c);
			bool y6 = Intersect_A(Surfaces[5], c);

			if (y1 == 1 || y2 == 1 || y3 == 1 || y4 == 1 || y5 == 1 || y6 == 1)
			{
				if (y1 == 1)
				{
					c.If_intersect_surfaces(0) = 1;
					//Modify_fracture_attributes_Zmax(c);
				}
				if (y2 == 1)
				{
					c.If_intersect_surfaces(1) = 1;
					//Modify_fracture_attributes_Zmin(c);
				}
				if (y3 == 1)
				{
					c.If_intersect_surfaces(2) = 1;
					//Modify_fracture_attributes_Ymin(c);
				}
				if (y4 == 1)
				{
					c.If_intersect_surfaces(3) = 1;
					//Modify_fracture_attributes_Ymax(c);
				}
				if (y5 == 1)
				{
					c.If_intersect_surfaces(4) = 1;
					//Modify_fracture_attributes_Xmin(c);
				}
				if (y6 == 1)
				{
					c.If_intersect_surfaces(5) = 1;
					//Modify_fracture_attributes_Xmax(c);
				}
				Fractures.Push(c);
				Fractures[Fractures.Size() - 1].Tag = Fractures.Size() - 1;
				return;
			}
		}
		return;
	}

	inline void Domain::Clusters()
	{
		/*boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS > G;

		for (size_t ic=0; ic< 0.5*Connections.Size(); ic++)
		{
			boost::add_edge(Connections[2*ic],Connections[2*ic+1],G);
		}
		
		std::vector< size_t > component(boost::num_vertices(G));
        size_t num_clusters = boost::connected_components(G, &component[0]);

		Listofclusters.Resize(num_clusters);

		for (size_t ic=0; ic < component.size(); ic++)
		{
			Fractures[ic].Clus = component[ic];
			Listofclusters[component[ic]].Push(ic);
		}*/
		//Program te routine to find the pairs of fractures that are connected

		//Util::Tree tree(Connections);
		//tree.GetClusters(Listofclusters);

		Graph_i GK(Fractures.Size(), Connections);
		GK.CreateGraph_i(Listofclusters);

		//Graph_x GK(Fractures.Size(), Connections, Listofclusters);
		//std::cout<<"----Connections----\n";
		//std::cout<<Connections<<std::endl;
		//std::cout<<"-------------------\n";

		//std::cout <<"Listofclusters.Size(): " << Listofclusters.Size() << std::endl;
		for (size_t i = 0; i < Listofclusters.Size(); i++)
		{
			//std::cout << "debug 1.2" << std::endl;
			//std::cout<<Listofclusters[i]<<std::endl;
			for (size_t j = 0; j < Listofclusters[i].Size(); j++)
			{
				//std::cout << "debug 1.3" << std::endl;
				Fractures[Listofclusters[i][j]].Clus = i;
			}
			//std::cout << "debug 1.4" << std::endl;
		}

		//std::cout << "debug 1.5" << std::endl;
	}

	inline bool Domain::Intersect(Fracture F1, Fracture F2)
	{

		Vec3_t A1;
		A1 = 0, 0, 0;

		Vec3_t B1;
		B1 = 0, 0, 0;

		Vec3_t dis_vec = F2.Center - F1.Center;
		double distance = pow(dot(dis_vec, dis_vec), 0.5);
		if (distance > F2.Radius + F1.Radius)
		{
			return false;
		}
		/* std::cout<<"debug\n"; */
		double pi = acos(-1);
		size_t e1, e2;
		Parallel_or_not(F1.Plane_parameter, F2.Plane_parameter, e1, e2);
		///std::cout<<"\ne1:"<<e1<<"; e2: " << e2 << std::endl;
		if (e1 == 1)
		{
			if (e2 == 1)
			{
				//two infinite plane are overlaped
				//in real DFN generation, the interval of each kind of input parameter is large enough, which seldom and even does not lead to two overlapped fractures
				//because it is a random process
				return false;
			}
			else
			{
				//two infinite plane are parallel
				///std::cout<<"The two fractrues are parallel but not overlapped!\n";
				return false;
			}
		}
		else if (e1 == 0)
		{
			Array<Vec3_t> Verts_1;
			Array<Vec3_t> Verts_2;
			Verts_1.resize(F1.Nvertices);
			Verts_2.resize(F2.Nvertices);

			Vec3_t temp1;
			Find_vector_2(F1.Normal_vector, temp1);

			double R_angle_temp1 = 0;
			double x_temp = F1.Dip_angle;
			R_angle_temp1 = -x_temp * pi / 180;
			Quaternion_t Q_axis_1;

			if (F1.Dip_angle > 0.0001)
			{

				NormalizeRotation(R_angle_temp1, temp1, Q_axis_1);
				for (size_t i = 0; i < F1.Nvertices; ++i)
				{
					Rotation(F1.Verts[i], Q_axis_1, Verts_1[i]);
				}
				for (size_t i = 0; i < F2.Nvertices; ++i)
				{
					Rotation(F2.Verts[i], Q_axis_1, Verts_2[i]);
				}
			}
			else
			{
				for (size_t i = 0; i < F1.Nvertices; ++i)
					Verts_1[i] = F1.Verts[i];
				for (size_t i = 0; i < F2.Nvertices; ++i)
					Verts_2[i] = F2.Verts[i];
			}

			///std::cout<<"\nafter rotation, first: \n"<<Verts_1<<std::endl<<"    second:\n"<<Verts_2<<"\n";

			//-------a piece of debug code----
			for (size_t i = 0; i < F1.Nvertices; i++)
			{
				size_t j = i + 1 - (size_t)((i + 1) / F1.Nvertices) * (i + 1);
				if (abs(Verts_1[i](2) - Verts_1[j](2)) > 0.001)
				{
					std::cout << "Error!!! The Z values of all vertexes of 1st fracture should be the same!\n";
					exit(0);
				}
			}
			//--------------------------------

			double MAX_Z = Find_max_z_value(Verts_2);
			double MIN_Z = Find_min_z_value(Verts_2);
			if (Verts_1[0](2) >= MIN_Z && Verts_1[0](2) <= MAX_Z)
			{

				///--------intersection line segment between horizontal plane and 2nd fracture
				Array<Vec3_t> Intersection_1;
				Intersection_between_2ndpolygon_and_infinite_plane(F2.Nvertices, Verts_1[0](2), Verts_2, Intersection_1);
				//std::cout << Intersection_1 << std::endl;

				///---------now extending the line segment----
				Array<Vec3_t> Intersection_infinite;
				Output_a_relatively_infinite_line(300, Intersection_1, Intersection_infinite);
				//cout<<"Infinite line: \n"<<Intersection_infinite[0]<<std::endl<<Intersection_infinite[1]<<"\n";

				///--------now, determine the coordinates of endpoints of intersection line segment between extending line and 1st fracture
				size_t numOfIntersectionPoint_1; ///which lies in [0,2];
				Array<Vec3_t> Intersection_2;
				//	std::cout<<"\nIntersection_infinite and 1st fracture:\n";
				Intersection_between_line_segment_and_polygon(numOfIntersectionPoint_1, Verts_1, Intersection_infinite, Intersection_2);
				if (numOfIntersectionPoint_1 > 2)
				{
					std::cout << "Error!!! There shoule not be more than two intersection points!\n";
					exit(0);
				}
				//std::cout<<"\nExtending line segment intersect 1st fracture "<< numOfIntersectionPoint_1<<" times"<<std::endl;

				///-------now, determine the intersection section between 1st and 2nd polygonal fractures
				if (numOfIntersectionPoint_1 == 0)
				{
					///std::cout<<"No intersection.\n";
					return false;
				}
				else if (numOfIntersectionPoint_1 == 1)
				{
					std::cout << "Error, infinite line must intersect 0 or 2 sides of the 1st fracture.\nIt is impossible just intersect 1 sides!!!\n";
					exit(0);
				}
				else
				{
					/// ------------ determine the intersection between intersection_1 and 1st fracture
					Array<Vec3_t> Intersection_3;
					size_t numOfIntersectionPoint_2; ///which also lies in [0,2]
					Intersection_between_line_segment_and_polygon(numOfIntersectionPoint_2, Verts_1, Intersection_1, Intersection_3);
					if (numOfIntersectionPoint_2 > 2)
					{
						std::cout << "Error!!! There shoule be no more than two intersection points!\n";
						exit(0);
					};
					//	std::cout<<"The Intersection_1 intersect the 1st fracture "<<numOfIntersectionPoint_2<<" times\n";
					///------------- now, determine the include angle between intersection_1 and x-axis
					double beta_1;
					double m_1 = Intersection_1[1](1) - Intersection_1[0](1);
					double l_1 = Intersection_1[1](0) - Intersection_1[0](0);
					if (m_1 < 0)
					{
						m_1 = -m_1;
						l_1 = -l_1;
					}
					beta_1 = atan2(m_1, l_1) * 180.0 / pi;
					if (beta_1 < 0)
						beta_1 = 360 + beta_1;
					//	std::cout<<"Intersection_1: \n"<<Intersection_1[0]<<"\n"<<Intersection_1[1]<<"\n"<<"beta_1 is: "<<beta_1<<std::endl;
					//	std::cout<<"Intersection_2: \n"<<Intersection_2[0]<<"\n"<<Intersection_2[1]<<"\n---------------------\n";
					///-----------

					///------------a piece of debuging code---
					double beta_2;
					double m_2 = Intersection_infinite[1](1) - Intersection_infinite[0](1);
					double l_2 = Intersection_infinite[1](0) - Intersection_infinite[0](0);
					if (m_2 < 0)
					{
						m_2 = -m_2;
						l_2 = -l_2;
					}
					beta_2 = atan2(m_2, l_2) * 180.0 / pi;
					if (beta_2 < 0)
						beta_2 = 360 + beta_2;
					//std::cout<<"beta_2 is: "<<beta_2<<std::endl;
					if (abs(beta_1 - beta_2) > 0.001)
					{
						std::cout << "The angle of infinit line and x-axis is incorrect!!!\n";
						exit(0);
					}
					///---------------------------------------

					///-------
					if (numOfIntersectionPoint_2 == 0)
					{
						//this means the intersection line (between 2nd fracture and horizontal plane) might be inside or outside the 1st fracture

						Array<Vec3_t> Intersection_4;
						Array<Vec3_t> Intersection_5;
						Intersection_4.resize(2);
						Intersection_5.resize(2);

						Vec3_t axis_z_2;
						axis_z_2 = 0, 0, 1;
						Quaternion_t Q_axis_z_2;
						NormalizeRotation(((180 - beta_1) * pi / 180.0), axis_z_2, Q_axis_z_2);
						Rotation(Intersection_1[0], Q_axis_z_2, Intersection_4[0]);
						Rotation(Intersection_1[1], Q_axis_z_2, Intersection_4[1]);
						Rotation(Intersection_2[0], Q_axis_z_2, Intersection_5[0]);
						Rotation(Intersection_2[1], Q_axis_z_2, Intersection_5[1]);
						///std::cout<<Intersection_4[0]<<"\n"<<Intersection_4[1]<<std::endl<<Intersection_5[0]<<std::endl<<Intersection_5[1]<<"\n";

						///--------------a piece error report   ----
						if (abs(Intersection_4[0](1) - Intersection_4[1](1)) > 0.001 || abs(Intersection_5[0](1) - Intersection_5[1](1)) > 0.001 || abs(Intersection_4[0](1) - Intersection_5[1](1)) > 0.001)
						{
							std::cout << "Error!!! The y values of the two intersections shoule be the same in this step!\n";
							exit(0);
						}
						///-----------------------------------------

						Array<Vec3_t> Intersection_6;
						Intersection_6.resize(2);
						Intersection_6[0](0) = 0.001;
						Intersection_6[1](0) = 0.002;

						Intersection_of_1D_intervals(Intersection_4, Intersection_5, Intersection_6);
						///------------- a piece of test code
						if (Intersection_6[0](0) == 0.001 && Intersection_6[1](0) == 0.002)
						{
							///std::cout<<"No Intersection, because the intersection between 2nd fracture and horizontal plane is totally outside the 1st fracture\n";
							return false;
						};
						//-------------

						///if the test code (above) do not work, it indicate the intersection is totally inside the 1st fracture

						///std::cout<<"\nIntersection is totally inside the 1st fracture.\n";
						Array<Vec3_t> Intersection_7; // true intersection points
						Intersection_7.resize(2);
						Intersection_1[0](2) = Verts_1[0](2);
						Intersection_1[1](2) = Verts_1[0](2);

						if (F1.Dip_angle > 0.0001)
						{
							Vec3_t axis_z_4;
							axis_z_4 = temp1;
							Quaternion_t Q_axis_z_4;
							NormalizeRotation(-R_angle_temp1, axis_z_4, Q_axis_z_4);
							Rotation(Intersection_1[0], Q_axis_z_4, Intersection_7[0]);
							Rotation(Intersection_1[1], Q_axis_z_4, Intersection_7[1]);
						}
						else
						{
							Intersection_7[0] = Intersection_1[0];
							Intersection_7[1] = Intersection_1[1];
						}
						A1 = Intersection_7[0];
						B1 = Intersection_7[1];
					}
					else if (numOfIntersectionPoint_2 == 1)
					{
						//this means the intersection line (between 2nd fracture and horizontal plane) intersect with one edge of fracture 1, and one end is inside the fracture 1
						//but we need to know which end is inside the 1st fracture, so which end is closer to the 1st fracture center, that end must be inside the fracture 1
						//also, we need to rotate the center of the 1st fracture, since all vertexes have been rotated
						///std::cout<<"\nOne end of the intersection is inside the 1st fracture.\n";
						Array<Vec3_t> Intersection_4;
						Array<Vec3_t> Intersection_5;
						Intersection_4.resize(2);
						Intersection_5.resize(2);

						Vec3_t axis_z_2;
						axis_z_2 = 0, 0, 1;
						Quaternion_t Q_axis_z_2;
						NormalizeRotation(((180 - beta_1) * pi / 180.0), axis_z_2, Q_axis_z_2);
						Rotation(Intersection_1[0], Q_axis_z_2, Intersection_4[0]);
						Rotation(Intersection_1[1], Q_axis_z_2, Intersection_4[1]);
						Rotation(Intersection_2[0], Q_axis_z_2, Intersection_5[0]);
						Rotation(Intersection_2[1], Q_axis_z_2, Intersection_5[1]);

						//std::cout<<Intersection_4[0]<<"\n"<<Intersection_4[1]<<std::endl<<Intersection_5[0]<<std::endl<<Intersection_5[1]<<"\n";
						///--------------a piece error report   ----
						if (abs(Intersection_4[0](1) - Intersection_4[1](1)) > 0.001 || abs(Intersection_5[0](1) - Intersection_5[1](1)) > 0.001 || abs(Intersection_4[0](1) - Intersection_5[1](1)) > 0.001)
						{
							std::cout << "Error!!! The y values of the two intersections shoule be the same in this step!\n";
							exit(0);
						}
						///-----------------------------------------
						Array<Vec3_t> Intersection_6;
						Intersection_6.resize(2);
						Intersection_6[0](0) = 0.001;
						Intersection_6[1](0) = 0.002;

						Intersection_of_1D_intervals(Intersection_4, Intersection_5, Intersection_6);

						///------------- a piece of test code
						if (Intersection_6[0](0) == 0.001 && Intersection_6[1](0) == 0.002)
						{
							///std::cout<<"No Intersection, because the intersection between 2nd fracture and horizontal plane is totally outside the 1st fracture\n";
							return false;
						};
						///----------------------------------

						//---------------------------
						Intersection_6[0](1) = Intersection_4[0](1);
						Intersection_6[1](1) = Intersection_4[0](1);
						///std::cout<<"y value: "<<Intersection_4[0](1)<<"\n";

						Intersection_6[0](2) = Verts_1[0](2);
						Intersection_6[1](2) = Verts_1[0](2);

						//	std::cout<<"Intersection_6: "<<Intersection_6[0]<<"\n"<<Intersection_6[1]<<"\n";

						Array<Vec3_t> Intersection_7;
						Intersection_7.resize(2);

						Vec3_t axis_z_3;
						axis_z_3 = 0, 0, 1;
						Quaternion_t Q_axis_z_3;
						NormalizeRotation(((180 + beta_1) * pi / 180.0), axis_z_3, Q_axis_z_3);

						Rotation(Intersection_6[0], Q_axis_z_3, Intersection_7[0]);
						Rotation(Intersection_6[1], Q_axis_z_3, Intersection_7[1]);

						//	std::cout<<"Intersection_7: "<<Intersection_7[0]<<"\n"<<Intersection_7[1]<<"\n";

						Array<Vec3_t> Intersection_8; /// true Intersection point
						Intersection_8.resize(2);

						if (F1.Dip_angle > 0.0001)
						{
							Vec3_t axis_z_4;
							axis_z_4 = temp1;
							Quaternion_t Q_axis_z_4;
							NormalizeRotation(-R_angle_temp1, axis_z_4, Q_axis_z_4);
							Rotation(Intersection_7[0], Q_axis_z_4, Intersection_8[0]);
							Rotation(Intersection_7[1], Q_axis_z_4, Intersection_8[1]);
						}
						else
						{
							Intersection_8[0] = Intersection_7[0];
							Intersection_8[1] = Intersection_7[1];
						}

						A1 = Intersection_8[0];
						B1 = Intersection_8[1];
					}
					else if (numOfIntersectionPoint_2 == 2)
					{
						//this means the intersection line (between 2nd fracture and horizontal plane) might be inside the 1st fracture
						///	std::cout<<"\nBoth two ends of intersection are on the perimeter (sides) of the 1st fracture.\n";
						Array<Vec3_t> Intersection_6;
						Intersection_6.resize(2);
						Array<Vec3_t> Intersection_7;
						Intersection_7.resize(2);

						Intersection_6[0] = Intersection_3[0];
						Intersection_6[1] = Intersection_3[1];
						Intersection_6[0](2) = Verts_1[0](2);
						Intersection_6[1](2) = Verts_1[0](2);

						if (F1.Dip_angle > 0.0001)
						{
							Vec3_t axis_z_4;
							axis_z_4 = temp1;
							Quaternion_t Q_axis_z_4;
							NormalizeRotation(-R_angle_temp1, axis_z_4, Q_axis_z_4);
							Rotation(Intersection_6[0], Q_axis_z_4, Intersection_7[0]);
							Rotation(Intersection_6[1], Q_axis_z_4, Intersection_7[1]);
						}
						else
						{
							Intersection_7[0] = Intersection_6[0];
							Intersection_7[1] = Intersection_6[1];
						}

						A1 = Intersection_7[0];
						B1 = Intersection_7[1];
					}
				}
			}
			else
			{
				return false;
			}
		};

		///std::cout<<"Intersection points are: \n"<<A1<<"\n"<<B1<<std::endl;

		if (A1(0) != 0 && A1(1) != 0 && A1(2) != 0 && B1(0) != 0 && B1(1) != 0 && B1(2) != 0)
		{
			//std::cout << "Found intersection between Fracture[" << F1.Tag << "] and Fracture[" << F2.Tag << "]!" << std::endl;
			///Clusters(F1,F2, A1, B1);

			Connections.Push(F1.Tag); //Fracture F1 and F2 are connected
			Connections.Push(F2.Tag);
			std::pair<size_t, size_t> p = std::make_pair(F1.Tag, F2.Tag); //Saving the intersection line from x1 to x2 into the Intersections map for the pair of fracture 0 and 1
			Vec3_t x1, x2;
			x1 = A1;
			x2 = B1;
			Intersections[p] = std::make_pair(x1, x2);
			return true;
		}
		else
			return false;
	}

	inline void Domain::Model_set(const double model_size[6])
	{

		double xmin = model_size[0];
		double xmax = model_size[1];
		double ymin = model_size[2];
		double ymax = model_size[3];
		double zmin = model_size[4];
		double zmax = model_size[5];
		/// top-----------
		size_t Tag_1 = 0;
		size_t Clus_1 = -1;

		Array<Vec3_t> Verts_1;
		Verts_1.resize(4);
		Verts_1[0] = xmin, ymin, zmax;
		Verts_1[1] = xmax, ymin, zmax;
		Verts_1[2] = xmax, ymax, zmax;
		Verts_1[3] = xmin, ymax, zmax;

		Fracture Top(Tag_1, Clus_1, Verts_1);

		/// bottom-----------
		size_t Tag_2 = 1;
		size_t Clus_2 = -1;

		Array<Vec3_t> Verts_2;
		Verts_2.resize(4);
		Verts_2[0] = xmin, ymin, zmin;
		Verts_2[1] = xmax, ymin, zmin;
		Verts_2[2] = xmax, ymax, zmin;
		Verts_2[3] = xmin, ymax, zmin;

		Fracture Bottom(Tag_2, Clus_2, Verts_2);

		/// front-----------
		size_t Tag_3 = 2;
		size_t Clus_3 = -1;

		Array<Vec3_t> Verts_3;
		Verts_3.resize(4);
		Verts_3[0] = xmin, ymin, zmin;
		Verts_3[1] = xmax, ymin, zmin;
		Verts_3[2] = xmax, ymin, zmax;
		Verts_3[3] = xmin, ymin, zmax;

		Fracture Front(Tag_3, Clus_3, Verts_3);

		/// back-----------
		size_t Tag_4 = 3;
		size_t Clus_4 = -1;

		Array<Vec3_t> Verts_4;
		Verts_4.resize(4);
		Verts_4[0] = xmin, ymax, zmin;
		Verts_4[1] = xmax, ymax, zmin;
		Verts_4[2] = xmax, ymax, zmax;
		Verts_4[3] = xmin, ymax, zmax;

		Fracture Back(Tag_4, Clus_4, Verts_4);

		/// left-----------
		size_t Tag_5 = 4;
		size_t Clus_5 = -1;

		Array<Vec3_t> Verts_5;
		Verts_5.resize(4);
		Verts_5[0] = xmin, ymin, zmin;
		Verts_5[1] = xmin, ymax, zmin;
		Verts_5[2] = xmin, ymax, zmax;
		Verts_5[3] = xmin, ymin, zmax;

		Fracture Left(Tag_5, Clus_5, Verts_5);

		/// right-----------
		size_t Tag_6 = 5;
		size_t Clus_6 = -1;

		Array<Vec3_t> Verts_6;
		Verts_6.resize(4);
		Verts_6[0] = xmax, ymin, zmin;
		Verts_6[1] = xmax, ymax, zmin;
		Verts_6[2] = xmax, ymax, zmax;
		Verts_6[3] = xmax, ymin, zmax;

		Fracture Right(Tag_6, Clus_6, Verts_6);

		///--------push the six surfaces into Fractures, so now, remember, in Fractures, we really have (Size (of Fractures) minus six) fractures
		/// they are Fractures[0] to [5]
		Surfaces.Push(Top);
		Surfaces.Push(Bottom);
		Surfaces.Push(Front);
		Surfaces.Push(Back);
		Surfaces.Push(Left);
		Surfaces.Push(Right);

		Model_domain = zmax, zmin, ymin, ymax, xmin, xmax;
	};

	inline void Domain::Modify_fracture_attributes_Zmax(Fracture &F2) ///modify vertexes, area, Nvertices
	{
		//std::cout << "Zmax called\n";
		Array<Vec3_t> Verts_temp1;
		size_t nt = 0;
		double Surf = Surfaces[0].Verts[0](2);
		/// Fractures[0], Top, Zmax

		if (F2.Verts[0](2) <= Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{

				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if (F2.Verts[i](2) <= Surf && Surf < F2.Verts[ni](2))
				{

					Verts_temp1.Push(F2.Verts[i]);
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](2)) / n;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = Surf;

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](2) >= Surf && Surf >= F2.Verts[ni](2)))
				{
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](2)) / n;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = Surf;

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
			}
			//third
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				Verts_temp1.Push(F2.Verts[i]);
			}
		}
		else if (F2.Verts[0](2) > Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](2) >= Surf && Surf > F2.Verts[ni](2))
				{
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);
					double t = (Surf - F2.Verts[i](2)) / n;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = Surf;
					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](2) <= Surf && Surf <= F2.Verts[ni](2)))
				{
					Verts_temp1.Push(F2.Verts[i]);

					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](2)) / n;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = Surf;

					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
			}
		}

		///-----------------
		F2.Nvertices = Verts_temp1.Size();
		F2.Verts.resize(0);
		for (size_t i = 0; i < F2.Nvertices; ++i)
			F2.Verts.Push(Verts_temp1[i]);
		///-------------------Area
		///Heron's formula

		F2.Area = 0;
		for (size_t i = 0; i < F2.Nvertices - 2; ++i)
		{
			size_t j = i + 1;
			size_t k = i + 2 - (size_t)((i + 2) / F2.Nvertices) * (i + 2);
			double a, b, c, p;
			a = pow(dot((F2.Verts[0] - F2.Verts[j]), (F2.Verts[0] - F2.Verts[j])), 0.5);
			b = pow(dot((F2.Verts[j] - F2.Verts[k]), (F2.Verts[j] - F2.Verts[k])), 0.5);
			c = pow(dot((F2.Verts[k] - F2.Verts[0]), (F2.Verts[k] - F2.Verts[0])), 0.5);
			p = (a + b + c) / 2;

			double Area_1;
			if (a == 0 || b == 0 || c == 0)
				Area_1 = 0;
			else
				Area_1 = pow((p * (p - a) * (p - b) * (p - c)), 0.5);
			F2.Area = F2.Area + Area_1;
		}

		///--------------Perimeter
		F2.Perimeter = 0;
		for (size_t i = 0; i < F2.Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (F2.Verts.Size())) * (i + 1);
			double p = pow(dot((F2.Verts[i] - F2.Verts[j]), (F2.Verts[i] - F2.Verts[j])), 0.5);
			F2.Perimeter = F2.Perimeter + p;
		}
	}

	inline void Domain::Modify_fracture_attributes_Zmin(Fracture &F2) ///modify vertexes, area, Nvertices
	{
		//	std::cout << "Zmin called\n";
		Array<Vec3_t> Verts_temp1;
		size_t nt = 0;
		double Surf = Surfaces[1].Verts[0](2);
		/// Fractures[1], Bottom, Zmin

		if (F2.Verts[0](2) >= Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{

				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if (F2.Verts[i](2) >= Surf && Surf > F2.Verts[ni](2))
				{

					Verts_temp1.Push(F2.Verts[i]);
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](2)) / n;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = Surf;

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](2) <= Surf && Surf <= F2.Verts[ni](2)))
				{
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](2)) / n;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = Surf;

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
			}
			//third
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				Verts_temp1.Push(F2.Verts[i]);
			}
		}
		else if (F2.Verts[0](2) < Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](2) <= Surf && Surf < F2.Verts[ni](2))
				{
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);
					double t = (Surf - F2.Verts[i](2)) / n;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = Surf;
					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](2) >= Surf && Surf >= F2.Verts[ni](2)))
				{
					Verts_temp1.Push(F2.Verts[i]);

					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](2)) / n;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = Surf;

					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
			}
		}

		///-----------------
		F2.Nvertices = Verts_temp1.Size();
		F2.Verts.resize(0);
		for (size_t i = 0; i < F2.Nvertices; ++i)
			F2.Verts.Push(Verts_temp1[i]);
		///-------------------Area
		///Heron's formula

		F2.Area = 0;
		for (size_t i = 0; i < F2.Nvertices - 2; ++i)
		{
			size_t j = i + 1;
			size_t k = i + 2 - (size_t)((i + 2) / F2.Nvertices) * (i + 2);
			double a, b, c, p;
			a = pow(dot((F2.Verts[0] - F2.Verts[j]), (F2.Verts[0] - F2.Verts[j])), 0.5);
			b = pow(dot((F2.Verts[j] - F2.Verts[k]), (F2.Verts[j] - F2.Verts[k])), 0.5);
			c = pow(dot((F2.Verts[k] - F2.Verts[0]), (F2.Verts[k] - F2.Verts[0])), 0.5);
			p = (a + b + c) / 2;

			double Area_1;
			if (a == 0 || b == 0 || c == 0)
				Area_1 = 0;
			else
				Area_1 = pow((p * (p - a) * (p - b) * (p - c)), 0.5);
			F2.Area = F2.Area + Area_1;
		}

		F2.Perimeter = 0;
		for (size_t i = 0; i < F2.Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (F2.Verts.Size())) * (i + 1);
			double p = pow(dot((F2.Verts[i] - F2.Verts[j]), (F2.Verts[i] - F2.Verts[j])), 0.5);
			F2.Perimeter = F2.Perimeter + p;
		}
	}

	inline void Domain::Modify_fracture_attributes_Ymin(Fracture &F2) ///modify vertexes, area, Nvertices
	{
		//	std::cout << "Ymin called\n";
		Array<Vec3_t> Verts_temp1;
		size_t nt = 0;
		double Surf = Surfaces[2].Verts[0](1);
		/// Fractures[2], Front, Ymin

		if (F2.Verts[0](1) >= Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{

				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if (F2.Verts[i](1) >= Surf && Surf > F2.Verts[ni](1))
				{

					Verts_temp1.Push(F2.Verts[i]);
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](1)) / m;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = Surf;
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](1) <= Surf && Surf <= F2.Verts[ni](1)))
				{
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](1)) / m;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = Surf;
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
			}
			//third
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				Verts_temp1.Push(F2.Verts[i]);
			}
		}
		else if (F2.Verts[0](1) < Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](1) <= Surf && Surf < F2.Verts[ni](1))
				{
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);
					double t = (Surf - F2.Verts[i](1)) / m;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = Surf;
					temp_A(2) = t * n + F2.Verts[i](2);
					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](1) >= Surf && Surf >= F2.Verts[ni](1)))
				{
					Verts_temp1.Push(F2.Verts[i]);

					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](1)) / m;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = Surf;
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
			}
		}

		///-----------------
		F2.Nvertices = Verts_temp1.Size();
		F2.Verts.resize(0);
		for (size_t i = 0; i < F2.Nvertices; ++i)
			F2.Verts.Push(Verts_temp1[i]);
		///-------------------Area
		///Heron's formula

		F2.Area = 0;
		for (size_t i = 0; i < F2.Nvertices - 2; ++i)
		{
			size_t j = i + 1;
			size_t k = i + 2 - (size_t)((i + 2) / F2.Nvertices) * (i + 2);
			double a, b, c, p;
			a = pow(dot((F2.Verts[0] - F2.Verts[j]), (F2.Verts[0] - F2.Verts[j])), 0.5);
			b = pow(dot((F2.Verts[j] - F2.Verts[k]), (F2.Verts[j] - F2.Verts[k])), 0.5);
			c = pow(dot((F2.Verts[k] - F2.Verts[0]), (F2.Verts[k] - F2.Verts[0])), 0.5);
			p = (a + b + c) / 2;

			double Area_1;
			if (a == 0 || b == 0 || c == 0)
				Area_1 = 0;
			else
				Area_1 = pow((p * (p - a) * (p - b) * (p - c)), 0.5);
			F2.Area = F2.Area + Area_1;
		}
		F2.Perimeter = 0;
		for (size_t i = 0; i < F2.Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (F2.Verts.Size())) * (i + 1);
			double p = pow(dot((F2.Verts[i] - F2.Verts[j]), (F2.Verts[i] - F2.Verts[j])), 0.5);
			F2.Perimeter = F2.Perimeter + p;
		}
	}

	inline void Domain::Modify_fracture_attributes_Ymax(Fracture &F2) ///modify vertexes, area, Nvertices
	{
		//	std::cout << "Ymax called\n";
		Array<Vec3_t> Verts_temp1;
		size_t nt = 0;
		double Surf = Surfaces[3].Verts[0](1);
		/// Fractures[2], Back, Ymax

		if (F2.Verts[0](1) <= Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{

				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if (F2.Verts[i](1) <= Surf && Surf < F2.Verts[ni](1))
				{

					Verts_temp1.Push(F2.Verts[i]);
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](1)) / m;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = Surf;
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](1) >= Surf && Surf >= F2.Verts[ni](1)))
				{
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](1)) / m;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = Surf;
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
			}
			//third
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				Verts_temp1.Push(F2.Verts[i]);
			}
		}
		else if (F2.Verts[0](1) > Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](1) >= Surf && Surf > F2.Verts[ni](1))
				{
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);
					double t = (Surf - F2.Verts[i](1)) / m;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = Surf;
					temp_A(2) = t * n + F2.Verts[i](2);
					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](1) <= Surf && Surf <= F2.Verts[ni](1)))
				{
					Verts_temp1.Push(F2.Verts[i]);

					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](1)) / m;
					Vec3_t temp_A;
					temp_A(0) = t * l + F2.Verts[i](0);
					temp_A(1) = Surf;
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
			}
		}

		///-----------------
		F2.Nvertices = Verts_temp1.Size();
		F2.Verts.resize(0);
		for (size_t i = 0; i < F2.Nvertices; ++i)
			F2.Verts.Push(Verts_temp1[i]);
		///-------------------Area
		///Heron's formula

		F2.Area = 0;
		for (size_t i = 0; i < F2.Nvertices - 2; ++i)
		{
			size_t j = i + 1;
			size_t k = i + 2 - (size_t)((i + 2) / F2.Nvertices) * (i + 2);
			double a, b, c, p;
			a = pow(dot((F2.Verts[0] - F2.Verts[j]), (F2.Verts[0] - F2.Verts[j])), 0.5);
			b = pow(dot((F2.Verts[j] - F2.Verts[k]), (F2.Verts[j] - F2.Verts[k])), 0.5);
			c = pow(dot((F2.Verts[k] - F2.Verts[0]), (F2.Verts[k] - F2.Verts[0])), 0.5);
			p = (a + b + c) / 2;

			double Area_1;
			if (a == 0 || b == 0 || c == 0)
				Area_1 = 0;
			else
				Area_1 = pow((p * (p - a) * (p - b) * (p - c)), 0.5);
			F2.Area = F2.Area + Area_1;
		}

		F2.Perimeter = 0;
		for (size_t i = 0; i < F2.Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (F2.Verts.Size())) * (i + 1);
			double p = pow(dot((F2.Verts[i] - F2.Verts[j]), (F2.Verts[i] - F2.Verts[j])), 0.5);
			F2.Perimeter = F2.Perimeter + p;
		}
	}

	inline void Domain::Modify_fracture_attributes_Xmin(Fracture &F2) ///modify vertexes, area, Nvertices
	{
		//std::cout << "Xmin called\n";
		Array<Vec3_t> Verts_temp1;
		size_t nt = 0;
		double Surf = Surfaces[4].Verts[0](0);
		/// Fractures[4], Left, Xmin

		if (F2.Verts[0](0) >= Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{

				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if (F2.Verts[i](0) >= Surf && Surf > F2.Verts[ni](0))
				{

					Verts_temp1.Push(F2.Verts[i]);
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](0)) / l;
					Vec3_t temp_A;
					temp_A(0) = Surf;
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](0) <= Surf && Surf <= F2.Verts[ni](0)))
				{
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](0)) / l;
					Vec3_t temp_A;
					temp_A(0) = Surf;
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
			}
			//third
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				Verts_temp1.Push(F2.Verts[i]);
			}
		}
		else if (F2.Verts[0](0) < Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](0) <= Surf && Surf < F2.Verts[ni](0))
				{
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);
					double t = (Surf - F2.Verts[i](0)) / l;
					Vec3_t temp_A;
					temp_A(0) = Surf;
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = t * n + F2.Verts[i](2);
					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](0) >= Surf && Surf >= F2.Verts[ni](0)))
				{
					Verts_temp1.Push(F2.Verts[i]);

					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](0)) / l;
					Vec3_t temp_A;
					temp_A(0) = Surf;
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
			}
		}

		///-----------------
		F2.Nvertices = Verts_temp1.Size();
		F2.Verts.resize(0);
		for (size_t i = 0; i < F2.Nvertices; ++i)
			F2.Verts.Push(Verts_temp1[i]);
		///-------------------Area
		///Heron's formula

		F2.Area = 0;
		for (size_t i = 0; i < F2.Nvertices - 2; ++i)
		{
			size_t j = i + 1;
			size_t k = i + 2 - (size_t)((i + 2) / F2.Nvertices) * (i + 2);
			double a, b, c, p;
			a = pow(dot((F2.Verts[0] - F2.Verts[j]), (F2.Verts[0] - F2.Verts[j])), 0.5);
			b = pow(dot((F2.Verts[j] - F2.Verts[k]), (F2.Verts[j] - F2.Verts[k])), 0.5);
			c = pow(dot((F2.Verts[k] - F2.Verts[0]), (F2.Verts[k] - F2.Verts[0])), 0.5);
			p = (a + b + c) / 2;

			double Area_1;
			if (a == 0 || b == 0 || c == 0)
				Area_1 = 0;
			else
				Area_1 = pow((p * (p - a) * (p - b) * (p - c)), 0.5);
			F2.Area = F2.Area + Area_1;
		}

		F2.Perimeter = 0;
		for (size_t i = 0; i < F2.Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (F2.Verts.Size())) * (i + 1);
			double p = pow(dot((F2.Verts[i] - F2.Verts[j]), (F2.Verts[i] - F2.Verts[j])), 0.5);
			F2.Perimeter = F2.Perimeter + p;
		}
	}
	inline void Domain::Modify_fracture_attributes_Xmax(Fracture &F2) ///modify vertexes, area, Nvertices
	{
		//	std::cout << "Xmax called\n";
		Array<Vec3_t> Verts_temp1;
		size_t nt = 0;
		double Surf = Surfaces[5].Verts[0](0);
		/// Fractures[5], Right, Xmax

		if (F2.Verts[0](0) <= Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{

				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if (F2.Verts[i](0) <= Surf && Surf < F2.Verts[ni](0))
				{

					Verts_temp1.Push(F2.Verts[i]);
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](0)) / l;
					Vec3_t temp_A;
					temp_A(0) = Surf;
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](0) >= Surf && Surf >= F2.Verts[ni](0)))
				{
					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](0)) / l;
					Vec3_t temp_A;
					temp_A(0) = Surf;
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);
					nt = i + 1;
					break;
				}
			}
			//third
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				Verts_temp1.Push(F2.Verts[i]);
			}
		}
		else if (F2.Verts[0](0) > Surf)
		{
			// first
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](0) >= Surf && Surf > F2.Verts[ni](0))
				{
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);
					double t = (Surf - F2.Verts[i](0)) / l;
					Vec3_t temp_A;
					temp_A(0) = Surf;
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = t * n + F2.Verts[i](2);
					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				nt = F2.Nvertices;
			}
			//second
			for (size_t i = nt; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);

				if ((F2.Verts[i](0) <= Surf && Surf <= F2.Verts[ni](0)))
				{
					Verts_temp1.Push(F2.Verts[i]);

					//need push intersect point
					double l = F2.Verts[ni](0) - F2.Verts[i](0);
					double m = F2.Verts[ni](1) - F2.Verts[i](1);
					double n = F2.Verts[ni](2) - F2.Verts[i](2);

					double t = (Surf - F2.Verts[i](0)) / l;
					Vec3_t temp_A;
					temp_A(0) = Surf;
					temp_A(1) = t * m + F2.Verts[i](1);
					temp_A(2) = t * n + F2.Verts[i](2);

					Verts_temp1.Push(temp_A);

					nt = i + 1;
					break;
				}
				else
				{
					Verts_temp1.Push(F2.Verts[i]);
				}
			}
		}

		///-----------------
		F2.Nvertices = Verts_temp1.Size();
		F2.Verts.resize(0);
		for (size_t i = 0; i < F2.Nvertices; ++i)
			F2.Verts.Push(Verts_temp1[i]);
		///-------------------Area
		///Heron's formula

		F2.Area = 0;
		for (size_t i = 0; i < F2.Nvertices - 2; ++i)
		{
			size_t j = i + 1;
			size_t k = i + 2 - (size_t)((i + 2) / F2.Nvertices) * (i + 2);
			double a, b, c, p;
			a = pow(dot((F2.Verts[0] - F2.Verts[j]), (F2.Verts[0] - F2.Verts[j])), 0.5);
			b = pow(dot((F2.Verts[j] - F2.Verts[k]), (F2.Verts[j] - F2.Verts[k])), 0.5);
			c = pow(dot((F2.Verts[k] - F2.Verts[0]), (F2.Verts[k] - F2.Verts[0])), 0.5);
			p = (a + b + c) / 2;

			double Area_1;
			if (a == 0 || b == 0 || c == 0)
				Area_1 = 0;
			else
				Area_1 = pow((p * (p - a) * (p - b) * (p - c)), 0.5);
			F2.Area = F2.Area + Area_1;
		}
		F2.Perimeter = 0;
		for (size_t i = 0; i < F2.Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (F2.Verts.Size())) * (i + 1);
			double p = pow(dot((F2.Verts[i] - F2.Verts[j]), (F2.Verts[i] - F2.Verts[j])), 0.5);
			F2.Perimeter = F2.Perimeter + p;
		}
	}

	inline void Domain::Create_whole_model(const size_t n, const double random_seed, const double model_size[6], const String str, const double array11[3][2], const double array12[4], const double array13[7])
	{

		Random_function r1 = Random_function(random_seed);

		Model_set(model_size);

		for (size_t i = 0; i < n; ++i)
		{
			Fracture f(str, i, r1, array11, array12, array13);
			AddSquareFracture(i, f);
		};

		size_t nz = Fractures.Size();

		/*for (size_t i = 0; i < nz; ++i)
		{
			Extract_intersection_between_surfaces_and_Fractures(Fractures[i]);
		}*/

		if (nz == 0)
		{
			P32_total = 0;
			P32_connected = 0;
			P30 = 0;
			Percolation_parameter = 0;
			Ratio_of_P32 = 0;
			Expected_area_times_perimeter = 0;
			return;
		};

		for (size_t i = 0; i < nz - 1; ++i)
		{
			for (size_t j = i + 1; j < nz; ++j)
			{
				Intersect(Fractures[i], Fractures[j]);
			}
		}
//std::cout << "debug 1" << std::endl;
#pragma omp critical
		{
			Clusters();
		}

		//std::cout << "debug 2" << std::endl;
		//std::cout<<"Lower boundary of radius: "<<array12[2]<<", Upper boundary of radius: "<<array12[3]<<", Alpha: "<<array12[1]<<"\n";
		Expected_area_and_perimeter(array12[2], array12[3], array12[1]);
	};

	inline void Domain::WriteFrac(char const *FileKey)
	{
		/*
		//Geometric information
		size_t N_Faces = 0;
		size_t N_Verts = 0;

		for (size_t nf = 0; nf < Fractures.Size(); nf++)
		{
			N_Faces += Fractures[nf].Nvertices;
			N_Verts += Fractures[nf].Nvertices + 1;
		}

		if (N_Faces <= 0)
			throw new Fatal("DFN: no fractures to plot");

		//Geometric information
		float *Verts = new float[3 * N_Verts];
		size_t *FaceCon = new size_t[3 * N_Faces];

		size_t n_verts = 0;
		size_t n_faces = 0;
		size_t n_attrs = 0;

		//Atributes
		size_t *Tags = new size_t[N_Faces];
		size_t *Clus = new size_t[N_Faces];

		for (size_t nf = 0; nf < Fractures.Size(); nf++)
		{
			size_t n_reff = n_verts / 3;
			Vec3_t C = Fractures[nf].Center;
			Verts[n_verts++] = C(0);
			Verts[n_verts++] = C(1);
			Verts[n_verts++] = C(2);
			for (size_t nv = 0; nv < Fractures[nf].Nvertices; nv++)
			{
				Verts[n_verts++] = float(Fractures[nf].Verts[nv](0));
				Verts[n_verts++] = float(Fractures[nf].Verts[nv](1));
				Verts[n_verts++] = float(Fractures[nf].Verts[nv](2));
				FaceCon[n_faces++] = size_t(n_reff);
				FaceCon[n_faces++] = size_t(nv + n_reff + 1);
				FaceCon[n_faces++] = size_t((nv + 1) % Fractures[nf].Nvertices + n_reff + 1);

				//Writing the attributes
				Tags[n_attrs] = size_t(Fractures[nf].Tag);
				Clus[n_attrs] = size_t(Fractures[nf].Clus);
				n_attrs++;
			}
		}

		//Write the data
		String fn(FileKey);
		fn.append(".h5");
		hid_t file_id;
		file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		hsize_t dims[1];
		String dsname;
		dims[0] = 3 * N_Verts;
		dsname.Printf("Verts");
		H5LTmake_dataset_float(file_id, dsname.CStr(), 1, dims, Verts);
		dims[0] = 3 * N_Faces;
		dsname.Printf("FaceCon");
		H5LTmake_dataset_int(file_id, dsname.CStr(), 1, dims, FaceCon);
		dims[0] = N_Faces;
		dsname.Printf("Tag");
		H5LTmake_dataset_int(file_id, dsname.CStr(), 1, dims, Tags);
		dims[0] = N_Faces;
		dsname.Printf("Cluster");
		H5LTmake_dataset_int(file_id, dsname.CStr(), 1, dims, Clus);

		//Erasing the data
		delete[] Verts;
		delete[] FaceCon;
		delete[] Tags;
		delete[] Clus;

		//Closing the file
		H5Fflush(file_id, H5F_SCOPE_GLOBAL);
		H5Fclose(file_id);

		//Writing xmf file
		std::ostringstream oss;
		oss << "<?xml version=\"1.0\" ?>\n";
		oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
		oss << "<Xdmf Version=\"2.0\">\n";
		oss << " <Domain>\n";
		oss << "   <Grid Name=\"DFN_Fractures\">\n";
		oss << "     <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << N_Faces << "\">\n";
		oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << N_Faces << " 3\">\n";
		oss << "        " << fn.CStr() << ":/FaceCon \n";
		oss << "       </DataItem>\n";
		oss << "     </Topology>\n";
		oss << "     <Geometry GeometryType=\"XYZ\">\n";
		oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << N_Verts << " 3\" >\n";
		oss << "        " << fn.CStr() << ":/Verts \n";
		oss << "       </DataItem>\n";
		oss << "     </Geometry>\n";
		oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
		oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
		oss << "        " << fn.CStr() << ":/Tag \n";
		oss << "       </DataItem>\n";
		oss << "     </Attribute>\n";
		oss << "     <Attribute Name=\"Cluster\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
		oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
		oss << "        " << fn.CStr() << ":/Cluster \n";
		oss << "       </DataItem>\n";
		oss << "     </Attribute>\n";
		oss << "   </Grid>\n";
		oss << " </Domain>\n";
		oss << "</Xdmf>\n";

		fn = FileKey;
		fn.append(".xmf");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();*/
	}

	inline void Domain::PlotMatlab_ORI_SCATTER(char const *FileKey)
	{
		std::ostringstream oss;
		double pi = acos(-1);
		oss << "th = [";
		for (size_t i = 0; i < Fractures.Size(); ++i)
		{
			double DD = Fractures[i].Dip_direction;
			double alpha = 0;
			if (DD > 90)
				alpha = 450 - DD;
			else if (DD <= 90)
				alpha = 90 - DD;
			alpha = alpha * pi / 180.0;
			oss << alpha << " ";
		}
		oss << "];\nr = [";

		for (size_t i = 0; i < Fractures.Size(); ++i)
		{
			double DA = Fractures[i].Dip_angle;
			double beta = DA;
			beta = beta * pi / 180.0;
			oss << beta << " ";
		}
		oss << "];\npolarscatter(th,r,'filled');\nhold on;\nrlim([0 " << pi / 2 << "]);\n";
		oss << "hold on;\nrticks([" << pi / 12 << " " << 2 * pi / 12 << " " << 3 * pi / 12 << " " << 4 * pi / 12 << " " << 5 * pi / 12 << " " << 6 * pi / 12 << "]);\n";
		oss << "set(gca,'thetaticklabel',[]);\n";
		oss << "set(gca,'rticklabel',[]);";
		//Open Matlab script to plot
		String fn(FileKey);
		fn.append(".m");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	}

	inline void Domain::PlotMatlab_DFN(char const *FileKey)
	{
		//Writing data
		std::ostringstream oss;

		//Plotting the fractures
		for (size_t nf = 0; nf < Fractures.Size(); nf++)
		{
			size_t n_verts = Fractures[nf].Verts.Size();
			oss << "fill3([";
			for (size_t nv = 0; nv < n_verts + 1; ++nv)
			{
				size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
				oss << Fractures[nf].Verts[nv_1](0) << " ";
			}
			oss << "],[";
			for (size_t nv = 0; nv < n_verts + 1; ++nv)
			{
				size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
				oss << Fractures[nf].Verts[nv_1](1) << " ";
			}
			oss << "],[";
			for (size_t nv = 0; nv < n_verts + 1; ++nv)
			{
				size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
				oss << Fractures[nf].Verts[nv_1](2) << " ";
			}
			oss << "],[rand rand rand]);\ngrid on;\nhold on;\n";
		}

		//Plotting the model domain
		for (size_t i = 0; i < 6; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
				oss << "plot3(";
				oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
				oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
				oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
				oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
			}
		}
		double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
		double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
		double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
		oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

		//Open Matlab script to plot
		String fn(FileKey);
		fn.append(".m");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	}

	inline void Domain::PlotMatlab_DFN_and_Intersection(char const *FileKey)
	{
		//Writing data
		std::ostringstream oss;

		//Plotting the fractures
		for (size_t nf = 0; nf < Fractures.Size(); nf++)
		{
			size_t n_verts = Fractures[nf].Verts.Size();
			oss << "fill3([";
			for (size_t nv = 0; nv < n_verts + 1; ++nv)
			{
				size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
				oss << Fractures[nf].Verts[nv_1](0) << " ";
			}
			oss << "],[";
			for (size_t nv = 0; nv < n_verts + 1; ++nv)
			{
				size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
				oss << Fractures[nf].Verts[nv_1](1) << " ";
			}
			oss << "],[";
			for (size_t nv = 0; nv < n_verts + 1; ++nv)
			{
				size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
				oss << Fractures[nf].Verts[nv_1](2) << " ";
			}
			oss << "],[rand rand rand]);\ngrid on;\nhold on;\n";
		}

		//Plotting the intersections
		//std::cout<<Connections.Size()<<std::endl;
		for (size_t nc = 0; nc < Connections.Size() / 2; nc++)
		{
			size_t i1 = Connections[2 * nc];
			size_t i2 = Connections[2 * nc + 1];
			Vec3_t x1 = Intersections[std::make_pair(i1, i2)].first;
			Vec3_t x2 = Intersections[std::make_pair(i1, i2)].second;
			oss << "plot3(";
			oss << "\[" << x1(0) << " " << x2(0) << "\],";
			oss << "\[" << x1(1) << " " << x2(1) << "\],";
			oss << "\[" << x1(2) << " " << x2(2) << "\],";
			oss << "'color',[0 0 1],'Linewidth',3);\ngrid on;\nhold on;\n";
		}

		//Plotting the model domain
		for (size_t i = 0; i < 6; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
				oss << "plot3(";
				oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
				oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
				oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
				oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
			}
		}
		double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
		double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
		double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
		oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

		//Open Matlab script to plot
		String fn(FileKey);
		fn.append(".m");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	}

	inline bool Domain::Intersect_A(Fracture F1, Fracture F2)
	{

		Vec3_t A1;
		A1 = 0, 0, 0;

		Vec3_t B1;
		B1 = 0, 0, 0;

		Vec3_t dis_vec = F2.Center - F1.Center;
		double distance = pow(dot(dis_vec, dis_vec), 0.5);
		if (distance > F2.Radius + F1.Radius)
		{
			return false;
		}
		double pi = acos(-1);
		size_t e1, e2;
		Parallel_or_not(F1.Plane_parameter, F2.Plane_parameter, e1, e2);
		///std::cout<<"\ne1:"<<e1<<"; e2: " << e2 << std::endl;
		if (e1 == 1)
		{
			if (e2 == 1)
			{
				//two infinite plane are overlaped
				//in real DFN generation, the interval of each kind of input parameter is large enough, which seldom and even does not lead to two overlapped fractures
				//because it is a random process
				return false;
			}
			else
			{
				//two infinite plane are parallel
				///std::cout<<"The two fractrues are parallel but not overlapped!\n";
				return false;
			}
		}
		else if (e1 == 0)
		{
			Array<Vec3_t> Verts_1;
			Array<Vec3_t> Verts_2;
			Verts_1.resize(F1.Nvertices);
			Verts_2.resize(F2.Nvertices);

			Vec3_t temp1;
			Find_vector_2(F1.Normal_vector, temp1);
			//std::cout<<"Debug Tag\n"<<"temp1:\n"<<temp1<<"\n";
			double R_angle_temp1 = 0;
			double x_temp = F1.Dip_angle;
			R_angle_temp1 = -x_temp * pi / 180;
			Quaternion_t Q_axis_1;

			if (F1.Dip_angle > 0.0001)
			{

				NormalizeRotation(R_angle_temp1, temp1, Q_axis_1);
				for (size_t i = 0; i < F1.Nvertices; ++i)
				{
					Rotation(F1.Verts[i], Q_axis_1, Verts_1[i]);
				}
				for (size_t i = 0; i < F2.Nvertices; ++i)
				{
					Rotation(F2.Verts[i], Q_axis_1, Verts_2[i]);
				}
			}
			else
			{
				for (size_t i = 0; i < F1.Nvertices; ++i)
					Verts_1[i] = F1.Verts[i];
				for (size_t i = 0; i < F2.Nvertices; ++i)
					Verts_2[i] = F2.Verts[i];
			}

			///std::cout<<"\nafter rotation, first: \n"<<Verts_1<<std::endl<<"    second:\n"<<Verts_2<<"\n";

			//-------a piece of debug code----

			for (size_t i = 0; i < F1.Nvertices; i++)
			{
				size_t j = i + 1 - (size_t)((i + 1) / F1.Nvertices) * (i + 1);
				if (abs(Verts_1[i](2) - Verts_1[j](2)) > 0.001)
				{
					std::cout << "Error!!! The Z values of all vertexes of 1st fracture should be the same! (Intersect_A)\n";
					exit(0);
				}
			}
			//--------------------------------

			double MAX_Z = Find_max_z_value(Verts_2);
			double MIN_Z = Find_min_z_value(Verts_2);
			if (Verts_1[0](2) >= MIN_Z && Verts_1[0](2) <= MAX_Z)
			{
				///--------intersection line segment between horizontal plane and 2nd fracture
				Array<Vec3_t> Intersection_1;
				Intersection_between_2ndpolygon_and_infinite_plane(F2.Nvertices, Verts_1[0](2), Verts_2, Intersection_1);

				///---------now extending the line segment----
				Array<Vec3_t> Intersection_infinite;
				Output_a_relatively_infinite_line(300, Intersection_1, Intersection_infinite);
				//cout<<"Infinite line: \n"<<Intersection_infinite[0]<<std::endl<<Intersection_infinite[1]<<"\n";

				///--------now, determine the coordinates of endpoints of intersection line segment between extending line and 1st fracture
				size_t numOfIntersectionPoint_1; ///which lies in [0,2];
				Array<Vec3_t> Intersection_2;
				//	std::cout<<"\nIntersection_infinite and 1st fracture:\n";
				Intersection_between_line_segment_and_polygon(numOfIntersectionPoint_1, Verts_1, Intersection_infinite, Intersection_2);
				if (numOfIntersectionPoint_1 > 2)
				{
					std::cout << "Error!!! There shoule not be more than two intersection points! (Intersect_A)\n";
					exit(0);
				}
				//std::cout<<"\nExtending line segment intersect 1st fracture "<< numOfIntersectionPoint_1<<" times"<<std::endl;

				///-------now, determine the intersection section between 1st and 2nd polygonal fractures
				if (numOfIntersectionPoint_1 == 0)
				{
					///std::cout<<"No intersection.\n";
					return false;
				}
				else if (numOfIntersectionPoint_1 == 1)
				{
					std::cout << "Error, infinite line must intersect 0 or 2 sides of the 1st fracture. (Intersect_A)\nIt is impossible just intersect 1 sides!!! (Intersect_A)";
					exit(0);
				}
				else
				{
					/// ------------ determine the intersection between intersection_1 and 1st fracture
					Array<Vec3_t> Intersection_3;
					size_t numOfIntersectionPoint_2; ///which also lies in [0,2]
					Intersection_between_line_segment_and_polygon(numOfIntersectionPoint_2, Verts_1, Intersection_1, Intersection_3);
					if (numOfIntersectionPoint_2 > 2)
					{
						std::cout << "Error!!! There shoule be no more than two intersection points! (Intersect_A)\n";
						exit(0);
					};
					//	std::cout<<"The Intersection_1 intersect the 1st fracture "<<numOfIntersectionPoint_2<<" times\n";
					///------------- now, determine the include angle between intersection_1 and x-axis
					double beta_1;
					double m_1 = Intersection_1[1](1) - Intersection_1[0](1);
					double l_1 = Intersection_1[1](0) - Intersection_1[0](0);
					if (m_1 < 0)
					{
						m_1 = -m_1;
						l_1 = -l_1;
					}
					beta_1 = atan2(m_1, l_1) * 180.0 / pi;
					//	std::cout<<"Intersection_1: \n"<<Intersection_1[0]<<"\n"<<Intersection_1[1]<<"\n"<<"beta_1 is: "<<beta_1<<std::endl;
					//	std::cout<<"Intersection_2: \n"<<Intersection_2[0]<<"\n"<<Intersection_2[1]<<"\n---------------------\n";
					///-----------

					///------------a piece of debuging code---
					double beta_2;
					double m_2 = Intersection_infinite[1](1) - Intersection_infinite[0](1);
					double l_2 = Intersection_infinite[1](0) - Intersection_infinite[0](0);
					if (m_2 < 0)
					{
						m_2 = -m_2;
						l_2 = -l_2;
					}
					beta_2 = atan2(m_2, l_2) * 180.0 / pi;
					//std::cout<<"beta_1 is: "<<beta_1<<std::endl;
					//std::cout<<"beta_2 is: "<<beta_2<<std::endl;
					if (beta_2 < 0)
						beta_2 = 360 + beta_2;
					if (abs(beta_1 - beta_2) > 0.001)
					{
						if ((beta_1 == 0 && beta_2 == 180) || (beta_1 == 180 && beta_2 == 0))
						{
							//beta_1;
						}
						else
						{
							std::cout << "The angle of infinit line and x-axis is incorrect!!! (Intersect_A)\n";
							exit(0);
						}
					}
					///---------------------------------------

					///-------
					if (numOfIntersectionPoint_2 == 0)
					{
						//this means the intersection line (between 2nd fracture and horizontal plane) might be inside or outside the 1st fracture

						Array<Vec3_t> Intersection_4;
						Array<Vec3_t> Intersection_5;
						Intersection_4.resize(2);
						Intersection_5.resize(2);

						Vec3_t axis_z_2;
						axis_z_2 = 0, 0, 1;
						Quaternion_t Q_axis_z_2;
						NormalizeRotation(((180 - beta_1) * pi / 180.0), axis_z_2, Q_axis_z_2);
						Rotation(Intersection_1[0], Q_axis_z_2, Intersection_4[0]);
						Rotation(Intersection_1[1], Q_axis_z_2, Intersection_4[1]);
						Rotation(Intersection_2[0], Q_axis_z_2, Intersection_5[0]);
						Rotation(Intersection_2[1], Q_axis_z_2, Intersection_5[1]);
						///std::cout<<Intersection_4[0]<<"\n"<<Intersection_4[1]<<std::endl<<Intersection_5[0]<<std::endl<<Intersection_5[1]<<"\n";

						///--------------a piece error report   ----
						if (abs(Intersection_4[0](1) - Intersection_4[1](1)) > 0.001 || abs(Intersection_5[0](1) - Intersection_5[1](1)) > 0.001 || abs(Intersection_4[0](1) - Intersection_5[1](1)) > 0.001)
						{
							std::cout << "Error!!! The y values of the two intersections shoule be the same in this step! (Intersect_A)\n";
							exit(0);
						}
						///-----------------------------------------

						Array<Vec3_t> Intersection_6;
						Intersection_6.resize(2);
						Intersection_6[0](0) = 0.001;
						Intersection_6[1](0) = 0.002;

						Intersection_of_1D_intervals(Intersection_4, Intersection_5, Intersection_6);
						///------------- a piece of test code
						if (Intersection_6[0](0) == 0.001 && Intersection_6[1](0) == 0.002)
						{
							///std::cout<<"No Intersection, because the intersection between 2nd fracture and horizontal plane is totally outside the 1st fracture\n";
							return false;
						};
						//-------------

						///if the test code (above) do not work, it indicate the intersection is totally inside the 1st fracture

						///std::cout<<"\nIntersection is totally inside the 1st fracture.\n";
						Array<Vec3_t> Intersection_7; // true intersection points
						Intersection_7.resize(2);
						Intersection_1[0](2) = Verts_1[0](2);
						Intersection_1[1](2) = Verts_1[0](2);

						if (F1.Dip_angle > 0.0001)
						{
							Vec3_t axis_z_4;
							axis_z_4 = temp1;
							Quaternion_t Q_axis_z_4;
							NormalizeRotation(-R_angle_temp1, axis_z_4, Q_axis_z_4);
							Rotation(Intersection_1[0], Q_axis_z_4, Intersection_7[0]);
							Rotation(Intersection_1[1], Q_axis_z_4, Intersection_7[1]);
						}
						else
						{
							Intersection_7[0] = Intersection_1[0];
							Intersection_7[1] = Intersection_1[1];
						}
						A1 = Intersection_7[0];
						B1 = Intersection_7[1];
					}
					else if (numOfIntersectionPoint_2 == 1)
					{
						//this means the intersection line (between 2nd fracture and horizontal plane) intersect with one edge of fracture 1, and one end is inside the fracture 1
						//but we need to know which end is inside the 1st fracture, so which end is closer to the 1st fracture center, that end must be inside the fracture 1
						//also, we need to rotate the center of the 1st fracture, since all vertexes have been rotated
						///std::cout<<"\nOne end of the intersection is inside the 1st fracture.\n";
						Array<Vec3_t> Intersection_4;
						Array<Vec3_t> Intersection_5;
						Intersection_4.resize(2);
						Intersection_5.resize(2);

						Vec3_t axis_z_2;
						axis_z_2 = 0, 0, 1;
						Quaternion_t Q_axis_z_2;
						NormalizeRotation(((180 - beta_1) * pi / 180.0), axis_z_2, Q_axis_z_2);
						Rotation(Intersection_1[0], Q_axis_z_2, Intersection_4[0]);
						Rotation(Intersection_1[1], Q_axis_z_2, Intersection_4[1]);
						Rotation(Intersection_2[0], Q_axis_z_2, Intersection_5[0]);
						Rotation(Intersection_2[1], Q_axis_z_2, Intersection_5[1]);

						//std::cout<<Intersection_4[0]<<"\n"<<Intersection_4[1]<<std::endl<<Intersection_5[0]<<std::endl<<Intersection_5[1]<<"\n";
						///--------------a piece error report   ----
						if (abs(Intersection_4[0](1) - Intersection_4[1](1)) > 0.001 || abs(Intersection_5[0](1) - Intersection_5[1](1)) > 0.001 || abs(Intersection_4[0](1) - Intersection_5[1](1)) > 0.001)
						{
							std::cout << "Error!!! The y values of the two intersections shoule be the same in this step! (Intersect_A)\n";
							exit(0);
						}
						///-----------------------------------------
						Array<Vec3_t> Intersection_6;
						Intersection_6.resize(2);
						Intersection_6[0](0) = 0.001;
						Intersection_6[1](0) = 0.002;

						Intersection_of_1D_intervals(Intersection_4, Intersection_5, Intersection_6);

						///------------- a piece of test code
						if (Intersection_6[0](0) == 0.001 && Intersection_6[1](0) == 0.002)
						{
							///std::cout<<"No Intersection, because the intersection between 2nd fracture and horizontal plane is totally outside the 1st fracture\n";
							return false;
						};
						///----------------------------------

						//---------------------------
						Intersection_6[0](1) = Intersection_4[0](1);
						Intersection_6[1](1) = Intersection_4[0](1);
						///std::cout<<"y value: "<<Intersection_4[0](1)<<"\n";

						Intersection_6[0](2) = Verts_1[0](2);
						Intersection_6[1](2) = Verts_1[0](2);

						//	std::cout<<"Intersection_6: "<<Intersection_6[0]<<"\n"<<Intersection_6[1]<<"\n";

						Array<Vec3_t> Intersection_7;
						Intersection_7.resize(2);

						Vec3_t axis_z_3;
						axis_z_3 = 0, 0, 1;
						Quaternion_t Q_axis_z_3;
						NormalizeRotation(((180 + beta_1) * pi / 180.0), axis_z_3, Q_axis_z_3);

						Rotation(Intersection_6[0], Q_axis_z_3, Intersection_7[0]);
						Rotation(Intersection_6[1], Q_axis_z_3, Intersection_7[1]);

						//	std::cout<<"Intersection_7: "<<Intersection_7[0]<<"\n"<<Intersection_7[1]<<"\n";

						Array<Vec3_t> Intersection_8; /// true Intersection point
						Intersection_8.resize(2);

						if (F1.Dip_angle > 0.0001)
						{
							Vec3_t axis_z_4;
							axis_z_4 = temp1;
							Quaternion_t Q_axis_z_4;
							NormalizeRotation(-R_angle_temp1, axis_z_4, Q_axis_z_4);
							Rotation(Intersection_7[0], Q_axis_z_4, Intersection_8[0]);
							Rotation(Intersection_7[1], Q_axis_z_4, Intersection_8[1]);
						}
						else
						{
							Intersection_8[0] = Intersection_7[0];
							Intersection_8[1] = Intersection_7[1];
						}

						A1 = Intersection_8[0];
						B1 = Intersection_8[1];
					}
					else if (numOfIntersectionPoint_2 == 2)
					{
						//this means the intersection line (between 2nd fracture and horizontal plane) might be inside the 1st fracture
						///	std::cout<<"\nBoth two ends of intersection are on the perimeter (sides) of the 1st fracture.\n";
						Array<Vec3_t> Intersection_6;
						Intersection_6.resize(2);
						Array<Vec3_t> Intersection_7;
						Intersection_7.resize(2);

						Intersection_6[0] = Intersection_3[0];
						Intersection_6[1] = Intersection_3[1];
						Intersection_6[0](2) = Verts_1[0](2);
						Intersection_6[1](2) = Verts_1[0](2);

						if (F1.Dip_angle > 0.0001)
						{
							Vec3_t axis_z_4;
							axis_z_4 = temp1;
							Quaternion_t Q_axis_z_4;
							NormalizeRotation(-R_angle_temp1, axis_z_4, Q_axis_z_4);
							Rotation(Intersection_6[0], Q_axis_z_4, Intersection_7[0]);
							Rotation(Intersection_6[1], Q_axis_z_4, Intersection_7[1]);
						}
						else
						{
							Intersection_7[0] = Intersection_6[0];
							Intersection_7[1] = Intersection_6[1];
						}

						A1 = Intersection_7[0];
						B1 = Intersection_7[1];
					}
				}
			}
			else
			{
				return false;
			}
		};

		///std::cout<<"Intersection points are: \n"<<A1<<"\n"<<B1<<std::endl;

		if (A1(0) != 0 && A1(1) != 0 && A1(2) != 0 && B1(0) != 0 && B1(1) != 0 && B1(2) != 0)
		{
			return true;
		}
		else
			return false;
	}

	inline void Domain::Extract_intersection_between_surfaces_and_Fractures(Fracture &F2)
	{
		if (F2.If_intersect_surfaces(0) == 1)
		{
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](2) == Model_domain(0) && F2.Verts[ni](2) == Model_domain(0))
				{
					Connections_S.Push(Fractures[0].Tag); //Top surface and F2 are connected
					Connections_S.Push(F2.Tag);
					std::pair<size_t, size_t> p = std::make_pair(Fractures[0].Tag, F2.Tag); //Saving the intersection line from x1 to x2 into the Intersections map for the pair of fracture 0 and 1
					Vec3_t x1, x2;
					x1 = F2.Verts[i];
					x2 = F2.Verts[ni];
					Intersections_S[p] = std::make_pair(x1, x2);
					/* std::cout<<x1<<", "<<x2<<std::endl; */
				}
			}
		};

		if (F2.If_intersect_surfaces(1) == 1)
		{
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](2) == Model_domain(1) && F2.Verts[ni](2) == Model_domain(1))
				{
					Connections_S.Push(Fractures[1].Tag); //Bottom surface and F2 are connected
					Connections_S.Push(F2.Tag);
					std::pair<size_t, size_t> p = std::make_pair(Fractures[1].Tag, F2.Tag); //Saving the intersection line from x1 to x2 into the Intersections map for the pair of fracture 0 and 1
					Vec3_t x1, x2;
					x1 = F2.Verts[i];
					x2 = F2.Verts[ni];
					Intersections_S[p] = std::make_pair(x1, x2);
					/* std::cout<<x1<<", "<<x2<<std::endl; */
				}
			}
		};

		if (F2.If_intersect_surfaces(2) == 1)
		{
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](1) == Model_domain(2) && F2.Verts[ni](1) == Model_domain(2))
				{
					Connections_S.Push(Fractures[2].Tag); //Front surface and F2 are connected
					Connections_S.Push(F2.Tag);
					std::pair<size_t, size_t> p = std::make_pair(Fractures[2].Tag, F2.Tag); //Saving the intersection line from x1 to x2 into the Intersections map for the pair of fracture 0 and 1
																							//Writing data
					std::ostringstream oss;
					Vec3_t x1, x2;
					x1 = F2.Verts[i];
					x2 = F2.Verts[ni];
					Intersections_S[p] = std::make_pair(x1, x2);
					/* std::cout<<x1<<", "<<x2<<std::endl; */
				}
			}
		};

		if (F2.If_intersect_surfaces(3) == 1)
		{
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](1) == Model_domain(3) && F2.Verts[ni](1) == Model_domain(3))
				{
					Connections_S.Push(Fractures[3].Tag); //Back surface and F2 are connected
					Connections_S.Push(F2.Tag);
					std::pair<size_t, size_t> p = std::make_pair(Fractures[3].Tag, F2.Tag); //Saving the intersection line from x1 to x2 into the Intersections map for the pair of fracture 0 and 1
					Vec3_t x1, x2;
					x1 = F2.Verts[i];
					x2 = F2.Verts[ni];
					Intersections_S[p] = std::make_pair(x1, x2);
					/* std::cout<<x1<<", "<<x2<<std::endl; */
				}
			}
		};

		if (F2.If_intersect_surfaces(4) == 1)
		{
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](0) == Model_domain(4) && F2.Verts[ni](0) == Model_domain(4))
				{
					Connections_S.Push(Fractures[4].Tag); //Left surface and F2 are connected
					Connections_S.Push(F2.Tag);
					std::pair<size_t, size_t> p = std::make_pair(Fractures[4].Tag, F2.Tag); //Saving the intersection line from x1 to x2 into the Intersections map for the pair of fracture 0 and 1
					Vec3_t x1, x2;
					x1 = F2.Verts[i];
					x2 = F2.Verts[ni];
					Intersections_S[p] = std::make_pair(x1, x2);
					/* std::cout<<x1<<", "<<x2<<std::endl; */
				}
			}
		};

		if (F2.If_intersect_surfaces(5) == 1)
		{
			for (size_t i = 0; i < F2.Nvertices; ++i)
			{
				size_t ni = i + 1 - (size_t)((i + 1) / F2.Nvertices) * (i + 1);
				if (F2.Verts[i](0) == Model_domain(5) && F2.Verts[ni](0) == Model_domain(5))
				{
					Connections_S.Push(Fractures[5].Tag); //Right surface and F2 are connected
					Connections_S.Push(F2.Tag);
					std::pair<size_t, size_t> p = std::make_pair(Fractures[5].Tag, F2.Tag); //Saving the intersection line from x1 to x2 into the Intersections map for the pair of fracture 0 and 1
					Vec3_t x1, x2;
					x1 = F2.Verts[i];
					x2 = F2.Verts[ni];
					Intersections_S[p] = std::make_pair(x1, x2);
					/* std::cout<<x1<<", "<<x2<<std::endl; */
				}
			}
		};
	};

	inline void Domain::PlotMatlab_Traces_on_Model_surfaces(char const *FileKey)
	{
		//Writing data
		std::ostringstream oss;

		//plotting the intersections between model surfaces and fractures

		for (size_t nc = 0; nc < Connections_S.Size() / 2; nc++)
		{
			size_t i1 = Connections_S[2 * nc];
			size_t i2 = Connections_S[2 * nc + 1];
			Vec3_t x1 = Intersections_S[std::make_pair(i1, i2)].first;
			Vec3_t x2 = Intersections_S[std::make_pair(i1, i2)].second;
			oss << "plot3(";
			oss << "\[" << x1(0) << " " << x2(0) << "\],";
			oss << "\[" << x1(1) << " " << x2(1) << "\],";
			oss << "\[" << x1(2) << " " << x2(2) << "\],";
			oss << "'color',[0 0 1],'Linewidth',3);\ngrid on;\nhold on;\n";
		}

		//Plotting the model domain
		for (size_t i = 0; i < 6; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
				oss << "plot3(";
				oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
				oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
				oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
				oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
			}
		}
		double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
		double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
		double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
		oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

		//Open Matlab script to plot
		String fn(FileKey);
		fn.append(".m");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	}

	inline void Domain::PlotMatlab_DFN_Highlight_Cluster(char const *FileKey)
	{
		//Writing data
		std::ostringstream oss;

		//Plotting the fractures, they are distinguished by clus values
		for (size_t i = 0; i < Listofclusters.Size(); ++i)
		{
			double rand_1 = random_double(0, 1);
			double rand_2 = random_double(0, 1);
			double rand_3 = random_double(0, 1);
			for (size_t j = 0; j < Listofclusters[i].Size(); ++j)
			{
				size_t nf = Listofclusters[i][j];
				size_t n_verts = Fractures[nf].Verts.Size();
				oss << "fill3([";
				for (size_t nv = 0; nv < n_verts + 1; ++nv)
				{
					size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
					oss << Fractures[nf].Verts[nv_1](0) << " ";
				}
				oss << "],[";
				for (size_t nv = 0; nv < n_verts + 1; ++nv)
				{
					size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
					oss << Fractures[nf].Verts[nv_1](1) << " ";
				}
				oss << "],[";
				for (size_t nv = 0; nv < n_verts + 1; ++nv)
				{
					size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
					oss << Fractures[nf].Verts[nv_1](2) << " ";
				}
				oss << "],[" << rand_1 << " " << rand_2 << " " << rand_3 << "]);\ngrid on;\nhold on;\n";
			}
		}

		//Plotting the model domain
		for (size_t i = 0; i < 6; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
				oss << "plot3(";
				oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
				oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
				oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
				oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
			}
		}
		double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
		double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
		double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
		oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

		//Open Matlab script to plot
		String fn(FileKey);
		fn.append(".m");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	}
	inline size_t Domain::Identify_percolation_clusters(String str)
	{
		if (Fractures.Size() == 0)
			return 0;
		Percolation_cluster.Resize(3);
		for (size_t i = 0; i < Listofclusters.Size(); ++i) //Z direction
		{
			size_t q1 = 0, q2 = 0;
			for (size_t j = 0; j < Listofclusters[i].Size(); ++j)
			{
				if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(0) == 1)
					q1 = 1;
				if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(1) == 1)
					q2 = 1;
				if (q1 == 1 && q2 == 1)
				{
					Percolation_cluster[2].Push(i);
					break;
				}
			}
		}

		for (size_t i = 0; i < Listofclusters.Size(); ++i) //Y direction
		{
			size_t q1 = 0, q2 = 0;
			for (size_t j = 0; j < Listofclusters[i].Size(); ++j)
			{
				if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(2) == 1)
					q1 = 1;
				if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(3) == 1)
					q2 = 1;
				if (q1 == 1 && q2 == 1)
				{
					Percolation_cluster[1].Push(i);
					break;
				}
			}
		}

		for (size_t i = 0; i < Listofclusters.Size(); ++i) //X direction
		{
			size_t q1 = 0, q2 = 0;
			for (size_t j = 0; j < Listofclusters[i].Size(); ++j)
			{
				if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(4) == 1)
					q1 = 1;
				if (Fractures[Listofclusters[i][j]].If_intersect_surfaces(5) == 1)
					q2 = 1;
				if (q1 == 1 && q2 == 1)
				{
					Percolation_cluster[0].Push(i);
					break;
				}
			}
		}
		if (str == "x")
		{
			if (Percolation_cluster[0].Size() == 0)
			{
				return 0;
			}
			else
				return 1;
		}
		else if (str == "y")
		{
			if (Percolation_cluster[1].Size() == 0)
			{
				return 0;
			}
			else
				return 1;
		}
		else if (str == "z")
		{
			if (Percolation_cluster[2].Size() == 0)
			{
				return 0;
			}
			else
				return 1;
		}
		else
		{
			std::cout << "Please define a percolation direction with x, y, z\n";
			exit(0);
		}
	}

	inline void Domain::PLotMatlab_DFN_Cluster_along_a_direction(char const *FileKey, string str)
	{
		//Writing data
		std::ostringstream oss;

		//Plotting the fractures
		Array<size_t> temp;
		if (str == "x")
		{
			temp.resize(Percolation_cluster[0].Size());
			temp = Percolation_cluster[0];
		}
		else if (str == "y")
		{
			temp.resize(Percolation_cluster[1].Size());
			temp = Percolation_cluster[1];
		}
		else if (str == "z")
		{
			temp.resize(Percolation_cluster[2].Size());
			temp = Percolation_cluster[2];
		}
		else
		{
			std::cout << "please define the percolation direction with char 'x', 'y' or 'z'\n";
			exit(0);
		};
		for (size_t i = 0; i < temp.Size(); ++i)
		{
			double rand_1 = random_double(0, 1);
			double rand_2 = random_double(0, 1);
			double rand_3 = random_double(0, 1);
			for (size_t j = 0; j < Listofclusters[temp[i]].Size(); ++j)
			{
				size_t nf = Listofclusters[temp[i]][j];
				size_t n_verts = Fractures[nf].Verts.Size();
				oss << "fill3([";
				for (size_t nv = 0; nv < n_verts + 1; ++nv)
				{
					size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
					oss << Fractures[nf].Verts[nv_1](0) << " ";
				}
				oss << "],[";
				for (size_t nv = 0; nv < n_verts + 1; ++nv)
				{
					size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
					oss << Fractures[nf].Verts[nv_1](1) << " ";
				}
				oss << "],[";
				for (size_t nv = 0; nv < n_verts + 1; ++nv)
				{
					size_t nv_1 = nv - (size_t)(nv / n_verts) * n_verts;
					oss << Fractures[nf].Verts[nv_1](2) << " ";
				}
				oss << "],[" << rand_1 << " " << rand_2 << " " << rand_3 << "]);\ngrid on;\nhold on;\n";
			}
		}

		//Plotting the model domain
		for (size_t i = 0; i < 6; ++i)
		{
			for (size_t j = 0; j < 4; ++j)
			{
				size_t nj = j + 1 - (size_t)((j + 1) / 4) * (j + 1);
				oss << "plot3(";
				oss << "[" << Surfaces[i].Verts[j](0) << " " << Surfaces[i].Verts[nj](0) << "],";
				oss << "[" << Surfaces[i].Verts[j](1) << " " << Surfaces[i].Verts[nj](1) << "],";
				oss << "[" << Surfaces[i].Verts[j](2) << " " << Surfaces[i].Verts[nj](2) << "],";
				oss << "'color',[1 0 0],'Linewidth',3);\ngrid on;\nhold on;\n";
			}
		}
		double xmin_1 = Model_domain(4), xmax_1 = Model_domain(5);
		double ymin_1 = Model_domain(2), ymax_1 = Model_domain(3);
		double zmin_1 = Model_domain(1), zmax_1 = Model_domain(0);
		oss << "axis([" << xmin_1 << " " << xmax_1 << " " << ymin_1 << " " << ymax_1 << " " << zmin_1 << " " << zmax_1 << "])\nhold on;\nxlabel('x (m)');\nylabel('y (m)');\nzlabel('z (m)');\ntitle('DFN');\n";

		//Open Matlab script to plot
		String fn(FileKey);
		fn.append(".m");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	};

	inline void Domain::Connectivity_uniform_orientation(string str)
	{

		if (Fractures.Size() == 0)
		{
			return;
		}
		Array<size_t> temp;
		if (str == "x")
		{
			temp.resize(Percolation_cluster[0].Size());
			temp = Percolation_cluster[0];
		}
		else if (str == "y")
		{
			temp.resize(Percolation_cluster[1].Size());
			temp = Percolation_cluster[1];
		}
		else if (str == "z")
		{
			temp.resize(Percolation_cluster[2].Size());
			temp = Percolation_cluster[2];
		}
		else
		{
			std::cout << "please define the percolation direction with char 'x', 'y' or 'z'\n";
			exit(0);
		}
		P32_connected = 0;
		P32_total = 0;
		double Area_connected = 0;
		double Area_total = 0;
		double Model_volume = (Model_domain(0) - Model_domain(1)) * (Model_domain(3) - Model_domain(2)) * (Model_domain(5) - Model_domain(4));

		for (size_t i = 0; i < temp.Size(); ++i)
		{
			for (size_t j = 0; j < Listofclusters[temp[i]].Size(); ++j)
			{
				size_t nf = Listofclusters[temp[i]][j];
				Area_connected = Area_connected + Fractures[nf].Area;
			}
		}
		P32_connected = Area_connected / Model_volume;
		for (size_t i = 0; i < Fractures.Size(); ++i)
		{
			Area_total = Area_total + Fractures[i].Area;
		}
		P32_total = Area_total / Model_volume;

		Ratio_of_P32 = P32_connected / P32_total;

		P30 = Fractures.Size();

		Percolation_parameter = P30 * 0.5 * Expected_area_times_perimeter;
	};

	inline void Domain::Connectivity_fisher_orientation(string str){

	};

	inline void Domain::Average_number_of_intersections_per_fracture()
	{
		n_I = (Connections.Size() / 2) / Fractures.Size();
	};
	inline void Domain::PlotMatlab_Radius_and_Area_kstest(char const *FileKey)
	{
		std::ostringstream oss;
		oss << "x1=[\n";
		for (size_t i = 0; i < Fractures.Size(); ++i)
		{
			oss << Fractures[i].Radius << "\n";
		}
		oss << "];\n";

		oss << "x2=[\n";
		for (size_t i = 0; i < Fractures.Size(); ++i)
		{
			oss << Fractures[i].Area << "\n";
		}
		oss << "];\n";

		oss << "[h,p,k] = kstest2(x1,x2);\n";
		oss << "%if h = 1, means the two groups of data are not having similar distributions;\n";
		oss << "hold on;\n";
		oss << "nbins = 30;\n";
		oss << "subplot(2,1,1);\n";
		oss << "histogram(x1,30);\n";
		oss << "hold on;\n";
		oss << "subplot(2,1,2);\n";
		oss << "%histogram(x2);\n";
		oss << "histogram(x2,nbins);\n";

		String fn(FileKey);
		fn.append(".m");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	};

	inline void Domain::PlotMatlab_Radius_and_Perimeter_kstest(char const *FileKey)
	{
		std::ostringstream oss;
		oss << "x1=[\n";
		for (size_t i = 0; i < Fractures.Size(); ++i)
		{
			oss << Fractures[i].Radius << "\n";
		}
		oss << "];\n";

		oss << "x2=[\n";
		for (size_t i = 0; i < Fractures.Size(); ++i)
		{
			oss << Fractures[i].Perimeter << "\n";
		}
		oss << "];\n";

		oss << "[h,p,k] = kstest2(x1,x2);\n";
		oss << "%if h = 1, means the two groups of data are not having similar distributions;\n";
		oss << "hold on;\n";
		oss << "nbins = 30;\n";
		oss << "subplot(2,1,1);\n";
		oss << "histogram(x1,30);\n";
		oss << "hold on;\n";
		oss << "subplot(2,1,2);\n";
		oss << "%histogram(x2);\n";
		oss << "histogram(x2,nbins);\n";

		String fn(FileKey);
		fn.append(".m");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	};
	inline void Domain::DataFile_Radius_AreaAndPerimeter(char const *FileKey)
	{
		std::ostringstream oss;
		oss << "Radius\tArea\tPerimeter\n";
		for (size_t i = 0; i < Fractures.Size(); ++i)
		{
			if (i == Fractures.Size() - 1)
				oss << Fractures[i].Radius << "\t" << Fractures[i].Area << "\t" << Fractures[i].Perimeter;
			oss << Fractures[i].Radius << "\t" << Fractures[i].Area << "\t" << Fractures[i].Perimeter << "\n";
		};
		String fn(FileKey);
		fn.append(".txt");
		std::ofstream of(fn.CStr(), std::ios::out);
		of << oss.str();
		of.close();
	}

	inline void Domain::Expected_area_and_perimeter(double x0, double x1, double alpha_k)
	{
		double alpha_g = -alpha_k;
		//double C1 = (1-alpha_g)/((2-alpha_g)*(pow(x1,1-alpha_g)-pow(x0,1-alpha_g)));
		//double Expected_value_of_radius = C1*pow(x1,2-alpha_g) - C1*pow(x0,2-alpha_g);
		//std::cout<<"Expected_value_of_radius: "<<Expected_value_of_radius<<std::endl;

		double C2 = (1 - alpha_g) / (pow(x1, 1 - alpha_g) - pow(x0, 1 - alpha_g));
		Expected_area_times_perimeter = 8 * pow(2, 0.5) * C2 * pow(x1, 4 - alpha_g) / (4 - alpha_g) - 8 * pow(2, 0.5) * C2 * pow(x0, 4 - alpha_g) / (4 - alpha_g);
	}

}; //namespace DFN

#endif // MECHSYS_DEM_DOMAIN_H
