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

#ifndef MECHSYS_DFN_FRACTURE_H
#define MECHSYS_DFN_FRACTURE_H

// MechSys cross (0, 0, 0)

#include <mechsys/mesh/mesh.h>
#include <mechsys/linalg/quaternion.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/dfn/Random_function.h>
#include <mechsys/dfn/Math.h>

namespace DFN
{

	class Fracture
	{
	public:
		// Data
		size_t Tag;		  ///< Tag to classify fractures
		int Clus;		  ///< Tag of cluster
		size_t Nvertices; ///< Number of vertices
		double Radius;	  ///< radius of circle (controlling the generation of irregular polygonal fracture)
		double Dip_direction;
		double Dip_angle; ///< orientation
		double Area;
		double Perimeter;

		Array<Vec3_t> Verts; ///< Array of vertices in cyclic order
		Vec3_t Center;		 ///< the center of the fractures (circle-controlled)
		Vec3_t Normal_vector;
		Vec4_t Plane_parameter;		  ///< a,b,c,d
		Vec6_t If_intersect_surfaces; /// top, bottom, front, back, left, right
		// Constructor
		Fracture(String string, size_t Tag, Random_function &c,
				 const double array1[3][2], const double array2[4], const double array3[7]); ///< Fracture constructor as an array of vertices

		Fracture(size_t _Tag, int _Clus, Array<Vec3_t> _Verts); ///< used to represent model surface
	};

	inline Fracture::Fracture(String string, size_t T, Random_function &c, const double array1[3][2], const double array2[4], const double array3[7])
	{
		If_intersect_surfaces = 0, 0, 0, 0, 0, 0;
		Tag = T;
		Clus = -1;
		double pi = acos(-1);

		//--------------------randomly generated fracture center
		Center(0) = c.unifrm(array1[0][0], array1[0][1]);
		Center(1) = c.unifrm(array1[1][0], array1[1][1]);
		Center(2) = c.unifrm(array1[2][0], array1[2][1]);

		//--------------------radius of the circle
		//Radius = c.lognor(array2[0], array2[1], array2[2], array2[3]);
		Radius = c.powerlaw(array2[2], array2[3], array2[1]);
		//--------------------dip direction and dip angle
		if (string == "fisher")
		{

			Dip_direction = c.fisher_alpha(array3[0], array3[1], array3[2], array3[3], array3[4], array3[5], array3[6]);
			Dip_angle = c.fisher_beta(array3[0], array3[1], array3[2], array3[3], array3[4], array3[5], array3[6]);
		}
		else if (string == "uniform")
		{
			double l_1 = 0;
			double m_1 = 0;
			double n_1 = 0;

			while ((l_1 == 0 && m_1 == 0 && n_1 == 0) || ((isnan(l_1)) || (isnan(m_1)) || (isnan(n_1))))
			{
				l_1 = c.unifrm(-1, 1);
				m_1 = c.unifrm(-1, 1);
				n_1 = c.unifrm(-1, 1);
			}

			if (n_1 < 0)
			{
				l_1 = -l_1;
				m_1 = -m_1;
				n_1 = -n_1;
			}
			//std::cout<<" *** "<<l_1<<", "<<m_1<<", "<<n_1;

			double beta_tmp = acos(n_1) * 180.0 / pi;
			double alpha_tmp = atan2(m_1, l_1) * 180 / pi;

			if (alpha_tmp < 0)
				alpha_tmp = 360 + alpha_tmp;
			Dip_angle = beta_tmp;
			if (alpha_tmp <= 90)
				Dip_direction = 90 - alpha_tmp;
			else if (alpha_tmp > 90)
				Dip_direction = 450 - alpha_tmp;
		}
		else
		{
			std::cout << "Error! Please define orientation distribution!\n";
			exit(0);
		}

		///---------------------a piece of debuging code---
		if (Dip_direction > 360 || Dip_angle > 90 || Dip_direction < 0 || Dip_angle < 0)
		{
			std::cout << "Error!!! The orientation is incorrect!\n";
			exit(0);
		};

		//-------------------------------------------------

		//--------------------random number of vertexes
		Nvertices = 4; //random_integer(4, 7); //no more sides, becasue more sides, more likely close to a circle

		//--------------------normal vector--------
		Find_normal_vec(Dip_direction, Dip_angle, Normal_vector);

		//--------------------coordinates of all vertexes in order
		Verts.resize(Nvertices);

		double sub_angle1 = (-360.0 / Nvertices) * pi / 180.0;
		Vec3_t lower1;
		lower1 = 0.0, Radius, 0.0;
		Vec3_t upper1;
		Vec3_t axis_z;
		axis_z = 0.0, 0.0, 1.0;

		Quaternion_t Q_axis_z1;
		NormalizeRotation(sub_angle1, axis_z, Q_axis_z1);
		Rotation(lower1, Q_axis_z1, upper1);
		Vec3_t temp1;
		temp1(0) = random_double((size_t)(lower1(0)), (size_t)(upper1(0)));
		temp1(1) = pow((Radius * Radius - temp1(0) * temp1(0)), 0.5);
		temp1(2) = 0;
		Verts[0] = temp1;
		for (size_t i = 1; i < Nvertices; i++)
		{
			Vec3_t temp2;
			double sub_angle2 = (i) * (-360.0 / Nvertices) * pi / 180.0;
			Quaternion_t Q_axis_z2;

			NormalizeRotation(sub_angle2, axis_z, Q_axis_z2);
			Rotation(temp1, Q_axis_z2, temp2);
			Verts[i] = temp2;
		};

		Vec3_t temp3;
		Find_vector_2(Normal_vector, temp3);
		if (abs(temp3(0)) < 0.000001 && abs(temp3(1)) < 0.000001 && abs(temp3(2)) < 0.000001)
		{
			//Verts;
		}
		else
		{
			double R_angle_temp1 = 0;
			double x_temp = Dip_angle; ///it is better to create a new variable to represent the dip angle, because debuging shows direct use of 'Dip_angle' to calculate rotation angle leads wrong output
			R_angle_temp1 = x_temp * pi / 180;

			Quaternion_t Q_axis_1;

			NormalizeRotation(R_angle_temp1, temp3, Q_axis_1);

			for (size_t i = 0; i < Nvertices; i++)
			{
				Vec3_t temp4;
				Rotation(Verts[i], Q_axis_1, temp4);
				Verts[i] = temp4;
			};
		}
		for (size_t i = 0; i < Nvertices; i++)
		{
			Verts[i] = Verts[i] + Center;
		};
		///------------------a piece of debuging code------

		for (size_t i = 0; i < Nvertices; ++i)
		{
			Vec3_t temp_1 = Center - Verts[i];
			double temp_radius = dot(temp_1, temp_1);
			temp_radius = pow(temp_radius, 0.5);
			if (abs(temp_radius - Radius) > 0.001)
			{
				std::cout << "Error!!! Distance from center to a vertex is not equal to Radius";
				exit(0);
			}
		}
		///------------------------------------------------

		//--------------------plane equation parameters: a,b,c,d
		Plane_parameter(0) = Normal_vector(0);
		Plane_parameter(1) = Normal_vector(1);
		Plane_parameter(2) = Normal_vector(2);
		Plane_parameter(3) = -dot(Normal_vector, Center);

		///----------------------Area
		Area = 0;
		/*for (size_t i = 0; i < Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (Verts.Size())) * (i + 1);
			Vec3_t Perpendicular_foot = (Verts[i] + Verts[j]) / 2;
			double bottom_side = pow(dot((Verts[i] - Verts[j]), (Verts[i] - Verts[j])), 0.5);
			double height = pow(dot((Center - Perpendicular_foot), (Center - Perpendicular_foot)), 0.5);
			Area = Area + bottom_side * height * 0.5;
		}*/
		Area = pow(2 * Radius * Radius, 0.5);
		Area = Area * Area;
		///--------------Perimeter
		Perimeter = 0;
		Perimeter = pow(2 * Radius * Radius, 0.5);
		Perimeter = 4 * Perimeter;
		/*for (size_t i = 0; i < Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (Verts.Size())) * (i + 1);
			double p = pow(dot((Verts[i] - Verts[j]),(Verts[i] - Verts[j])),0.5);
			Perimeter = Perimeter+p;
		}*/
	};

	inline Fracture::Fracture(size_t _Tag, int _Clus, Array<Vec3_t> _Verts)
	{
		If_intersect_surfaces = 0, 0, 0, 0, 0, 0;
		Tag = _Tag;
		Clus = _Clus;
		double pi = acos(-1);
		Nvertices = _Verts.Size();
		Radius = pow(dot((_Verts[0] - _Verts[2]), (_Verts[0] - _Verts[2])), 0.5) * 0.5;
		double x1, y1, z1, x2, y2, z2, x3, y3, z3;
		double l, m, n, d;
		x1 = _Verts[0](0);
		y1 = _Verts[0](1);
		z1 = _Verts[0](2);

		x2 = _Verts[1](0);
		y2 = _Verts[1](1);
		z2 = _Verts[1](2);

		x3 = _Verts[2](0);
		y3 = _Verts[2](1);
		z3 = _Verts[2](2);

		l = (y3 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
		m = (x3 - x1) * (z2 - z1) - (x2 - x1) * (z3 - z1);
		n = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

		if (n < 0)
		{
			l = -l;
			m = -m;
			n = -n;
		};
		d = -(l * x1 + m * y1 + n * z1);
		Normal_vector = l, m, n;
		Plane_parameter = l, m, n, d;
		Verts.resize(Nvertices);
		for (size_t i = 0; i < Nvertices; ++i)
			Verts[i] = _Verts[i];

		double beta = acos(n / norm(Normal_vector)) * 180.0 / pi;
		double alpha = atan2(m, l) * 180.0 / pi;
		if (alpha < 0)
			alpha = 360 + alpha;

		Dip_angle = beta;
		if (alpha <= 90)
			Dip_direction = 90 - alpha;
		else if (alpha > 90)
			Dip_direction = 450 - alpha;
		Area = pow(Radius, 2) * 0.5 * 4;

		Center = (_Verts[0] + _Verts[2]) / 2;

		///--------------Perimeter
		Perimeter = 0;
		for (size_t i = 0; i < Verts.Size(); ++i)
		{
			size_t j = i + 1 - (size_t)((i + 1) / (Verts.Size())) * (i + 1);
			double p = pow(dot((Verts[i] - Verts[j]), (Verts[i] - Verts[j])), 0.5);
			Perimeter = Perimeter + p;
		}
		//std::cout<<Normal_vector<<std::endl;
		//std::cout<<Dip_direction<<", "<<Dip_angle<<std::endl;
	};

	/*
inline void Fracture::Centroid (Vec3_t & C)
{
    C = Vec3_t(0.0,0.0,0.0);
    for (size_t i=0; i<Nvertices; i++)
    {
        C += Verts[i];
    }
    C /= Nvertices;
}

inline void Fracture::Normal (Vec3_t & N)
{
    Vec3_t x1,x2;
    x1 = Verts[1] - Verts[0];
    x2 = Verts[2] - Verts[1];
    N  = cross(x1,x2);
    N /= norm(N);

}*/
} // namespace DFN

#endif // MECHSYS_DEM_FRACTURE_H
