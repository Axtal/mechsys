
#ifndef MATH_H
#define MATH_H

#include <iostream>
#include <ctime>
#include <cmath>
#include <mechsys/linalg/matvec.h>
#include <mechsys/dfn/Line.h>

#define Random(low, up) (rand() % (up - low + 1)) + low // used in function: random_integer

using namespace std;

size_t random_integer(size_t x, size_t y)
{
	return Random(x, y);
};

double random_double(size_t min, size_t max)//generate random numbers
{
	double m1 = (double)(rand() % 101) / 101;
	min++;
	double m2 = (double)((rand() % (max - min + 1)) + min);
	m2 = m2 - 1;
	return m1 + m2;
};

void Find_normal_vec(double dip_direction, double dip_angle, Vec3_t &a) /// based on orientation, find normal vector
{
	///spherical system firstly
	double alpha = 0, beta = 0;
	beta = dip_angle;
	if (dip_direction >= 90)
		alpha = 450 - dip_direction;
	else if (dip_direction <= 90)
		alpha = 90 - dip_direction;

	//------------------
	double pi = acos(-1);
	a(0) = sin(beta * pi / 180) * cos(alpha / 180.0 * pi);
	a(1) = sin(beta / 180.0 * pi) * sin(alpha / 180.0 * pi);
	a(2) = cos(beta / 180.0 * pi);
};

//the function below finds the vector that (1) is vertical to fracture normal vector; and (2) lies on the horizontal plane (z = 0)
void Find_vector_2(Vec3_t Normal_vector, Vec3_t &temp3) 
{
	temp3(2) = 0;
	if (Normal_vector(0) > 0)
	{
		if (Normal_vector(1) > 0)
		{
			temp3(0) = -Normal_vector(1);
			temp3(1) = Normal_vector(0);
		}
		else if (Normal_vector(1) < 0)
		{
			temp3(0) = Normal_vector(1);
			temp3(1) = -Normal_vector(0);
		}
		else if (Normal_vector(1) == 0)
		{
			temp3(0) = 0;
			temp3(1) = Normal_vector(0);
		}
	}
	else if (Normal_vector(0) < 0)
	{
		if (Normal_vector(1) < 0)
		{
			temp3(0) = -Normal_vector(1);
			temp3(1) = Normal_vector(0);
		}
		else if (Normal_vector(1) > 0)
		{
			temp3(0) = -Normal_vector(1);
			temp3(1) = Normal_vector(0);
		}
		else if (Normal_vector(1) == 0)
		{
			temp3(0) = 0;
			temp3(1) = Normal_vector(0);
		}
	}
	else if (Normal_vector(0) == 0)
	{
		if (Normal_vector(1) > 0)
		{
			temp3(0) = -Normal_vector(1);
			temp3(1) = 0;
		}
		else if (Normal_vector(1) < 0)
		{
			temp3(0) = -Normal_vector(1);
			temp3(1) = 0;
		}
		else if (Normal_vector(1) == 0)
		{
			temp3(0) = 0;
			temp3(1) = 0;
			return;
		}
	};
}

//the function below determine whether two infinite planes are parallel or not
void Parallel_or_not(Vec4_t plane_parameter1, Vec4_t plane_parameter2, size_t &index1, size_t &index2)
{
	double a1 = plane_parameter1(0);
	double b1 = plane_parameter1(1);
	double c1 = plane_parameter1(2);
	double d1 = plane_parameter1(3);
	double a2 = plane_parameter2(0);
	double b2 = plane_parameter2(1);
	double c2 = plane_parameter2(2);
	double d2 = plane_parameter2(3);

	index1 = 0; //when index1 = 1, means parallel; if index1 = 0, means two infinite planes intersect
	index2 = 0; //0 means parallel but not overlap, 1 means overlap

	if ((a1 == 0 && a2 != 0) || (a1 != 0 && a2 == 0))
	{
		index1 = 0;
		return;
	};
	if ((b1 == 0 && b2 != 0) || (b1 != 0 && b2 == 0))
	{
		index1 = 0;
		return;
	};
	if ((c1 == 0 && c2 != 0) || (c1 != 0 && c2 == 0))
	{
		index1 = 0;
		return;
	};

	if (a1 != 0 && b1 != 0 && c1 != 0 && a2 != 0 && b2 != 0 && c2 != 0)
	{
		double ratio_a = a1 / a2;
		double ratio_b = b1 / b2;
		double ratio_c = c1 / c2;

		if (abs(ratio_a - ratio_b) < 0.001 && abs(ratio_a - ratio_c) < 0.001 && abs(ratio_b - ratio_c) < 0.001)
		{
			index1 = 1; //parallel
			double ratio_d = d1 / d2;
			if (abs(ratio_d - ratio_a) < 0.001)
			{
				index2 = 1;
			}
			else
			{

				index2 = 0;
			}
		}
		else
		{
			index1 = 0;
		}
	}
	else
	{
		double array1[4] = {a1, b1, c1, d1};
		double array2[4] = {a2, b2, c2, d2};
		double array3[4] = {0, 0, 0, 0};

		size_t index_tmp = 0;
		for (size_t i = 0; i < 4; i++)
		{
			if (array1[i] == 0)
			{
				array3[i] = 1;
				if (i < 3)
					index_tmp++;
			};
		}

		if (index_tmp == 2)
		{
			//which means the two normal vectors perpendicular to xy, yz or xz plane
			index1 = 1;
			size_t k = 0;
			for (k = 0; k < 3; k++)
			{
				if (array3[k] != 1)
				{
					break;
				}
			}

			if (abs(d1 / array1[k] - d2 / array2[k]) < 0.001)
			{
				index2 = 1;
			}
			else
			{
				index2 = 0;
			}
		}
		else if (index_tmp == 1)
		{
			//only one of the three parameters (a, b and c) is zero
			size_t k1 = 0;
			size_t k2 = 0;

			for (k1 = 0; k1 < 3; k1++)
			{
				if (array1[k1] != 0)
				{
					k2 = k1 + 1;
					break;
				}
			};

			for (k2 = 0; k2 < 3; k2++)
			{
				if (array1[k2] != 0)
				{
					break;
				}
			};
			double ratio_k1 = array1[k1] / array2[k1];
			double ratio_k2 = array1[k2] / array2[k2];
			if (abs(ratio_k1 - ratio_k2) < 0.001)
			{
				index1 = 1;
				if ((d1 / array1[k1] - d2 / array2[k1]) < 0.001)
				{
					index2 = 1;
				}
				else
				{
					index2 = 0;
				}
			}
			else
			{
				index1 = 0;
			}
		}
	}
};

//the function below finds the maximum element
double Find_max_z_value(Array<Vec3_t> A)
{

	double Max = 0;
	for (size_t i = 0; i < A.size() - 1; i++)
	{
		if (i == 0)
			Max = max(A[0](2), A[1](2));
		else
			Max = max(Max, A[i + 1](2));
	};
	return Max;
};

//below finds minimum element
double Find_min_z_value(Array<Vec3_t> A)
{

	double Min = 0;
	for (size_t i = 0; i < A.size() - 1; i++)
	{
		if (i == 0)
			Min = min(A[0](2), A[1](2));
		else
			Min = min(Min, A[i + 1](2));
	}
	return Min;
};

//below finds the intersection between a fracture and a infinite plane (this plane must be horizontal) 
void Intersection_between_2ndpolygon_and_infinite_plane(size_t vertexes_of_polygon_, double z_zplane, Array<Vec3_t> array, Array<Vec3_t> &Intersection_1)
{

	double end_point1_x = 0;
	double end_point1_y = 0;
	double end_point2_x = 0;
	double end_point2_y = 0;

	size_t k = 0;
	for (size_t i = 0; i < vertexes_of_polygon_; i++)
	{
		if (i == vertexes_of_polygon_ - 1)
		{
			if ((array[i](2) <= z_zplane && z_zplane <= array[i + 1](2)) || (array[i](2) >= z_zplane && z_zplane >= array[i + 1](2)))
			{
				double t = (z_zplane - array[0](2)) / (array[i](2) - array[0](2));
				end_point1_x = array[0](0) + (array[i](0) - array[0](0)) * t;
				end_point1_y = array[0](1) + (array[i](1) - array[0](1)) * t;
				k = i + 1;

				break;
			}
		}
		else
		{
			if ((array[i](2) <= z_zplane && z_zplane <= array[i + 1](2)) || (array[i](2) >= z_zplane && z_zplane >= array[i + 1](2)))
			{
				double t = (z_zplane - array[i + 1](2)) / (array[i](2) - array[i + 1](2));
				end_point1_x = array[i + 1](0) + (array[i](0) - array[i + 1](0)) * t;
				end_point1_y = array[i + 1](1) + (array[i](1) - array[i + 1](1)) * t;
				k = i + 1;

				break;
			}
		}
	}

	for (size_t j = k; j < vertexes_of_polygon_; j++)
	{

		if (j == vertexes_of_polygon_ - 1)
		{

			if ((array[j](2) <= z_zplane && z_zplane <= array[0](2)) || (array[j](2) >= z_zplane && z_zplane >= array[0](2)))
			{
				double t = (z_zplane - array[0](2)) / (array[j][2] - array[0](2));

				end_point2_x = array[0](0) + (array[j](0) - array[0](0)) * t;
				end_point2_y = array[0](1) + (array[j](1) - array[0](1)) * t;

				j = vertexes_of_polygon_;
			}
		}
		else
		{
			if ((array[j](2) <= z_zplane && z_zplane <= array[j + 1](2)) || (array[j](2) >= z_zplane && z_zplane >= array[j + 1](2)))
			{
				double t = (z_zplane - array[j + 1](2)) / (array[j](2) - array[j + 1](2));

				end_point2_x = array[j + 1](0) + (array[j](0) - array[j + 1](0)) * t;
				end_point2_y = array[j + 1](1) + (array[j](1) - array[j + 1](1)) * t;

				j = vertexes_of_polygon_;
			}
		}
		if (j == vertexes_of_polygon_ - 1)
			break;
	}
	Intersection_1.resize(2);
	Intersection_1[0] = end_point1_x, end_point1_y, 0;
	Intersection_1[1] = end_point2_x, end_point2_y, 0;
}

//below finds a relatively infinte line, i.e., absolute values of coordinates of both two ends are very large
//also they are in horizontal plane, or say, 2D Cartesian system
void Output_a_relatively_infinite_line(double max, Array<Vec3_t> Intersection_1, Array<Vec3_t> &Intersection_infinite)
{
	double x0, y0, x1, y1, x2, y2, x3, y3;
	x0 = Intersection_1[0](0);
	y0 = Intersection_1[0](1);
	x1 = Intersection_1[1](0);
	y1 = Intersection_1[1](1);
	x2 = max;
	x3 = -max;
	if ((x1 - x0) != 0)
	{
		double t1 = (x2 - x0) / (x1 - x0);
		y2 = t1 * (y1 - y0) + y0;

		double t2 = (x3 - x0) / (x1 - x0);
		y3 = t2 * (y1 - y0) + y0;
	}
	else
	{
		y2 = x2;
		y3 = x3;
		x2 = x0;
		x3 = x0;
	}
	Intersection_infinite.resize(2);
	Intersection_infinite[0](0) = x2;
	Intersection_infinite[0](1) = y2;
	Intersection_infinite[1](0) = x3;
	Intersection_infinite[1](1) = y3;
}

//below determines whether two lines are intersect, also, in 2D Cartesian system
size_t is_intersect(Line myline1, Line myline2) //deterine if two lines intersect
{
	if (myline1.get_max_x() < myline2.get_min_x() ||
		myline2.get_max_x() < myline1.get_min_x() ||
		myline1.get_max_y() < myline2.get_min_y() ||
		myline2.get_max_y() < myline1.get_min_y())
		return 0; //0 means disconnected
	double res1 = (myline1.xa - myline1.xb) * (myline2.ya - myline1.yb) - (myline1.ya - myline1.yb) * (myline2.xa - myline1.xb);
	double res2 = (myline1.xa - myline1.xb) * (myline2.yb - myline1.yb) - (myline1.ya - myline1.yb) * (myline2.xb - myline1.xb);

	double res3 = (myline2.xa - myline2.xb) * (myline1.ya - myline2.yb) - (myline2.ya - myline2.yb) * (myline1.xa - myline2.xb);
	double res4 = (myline2.xa - myline2.xb) * (myline1.yb - myline2.yb) - (myline2.ya - myline2.yb) * (myline1.xb - myline2.xb);
	if (res1 * res2 <= 0 && res3 * res4 <= 0)
		return 1; //1 means connected
	else
		return 0;
};

//below finds the intersection between two line segments (if they have)
void intersection_between_two_line_segments(Array<Vec3_t> Intersection_infinite, Vec3_t A, Vec3_t B, double &intersection_x, double &intersection_y)
{
	double x0 = Intersection_infinite[0](0);
	double y0 = Intersection_infinite[0](1);
	double x1 = Intersection_infinite[1](0);
	double y1 = Intersection_infinite[1](1);
	double x2 = A(0);
	double y2 = A(1);
	double x3 = B(0);
	double y3 = B(1);

	if ((x1 - x0) != 0 && (x3 - x2) != 0)
	{
		double k1, k2;
		k1 = (y1 - y0) / (x1 - x0);
		k2 = (y3 - y2) / (x3 - x2);
		if (k1 != 0 && k2 != 0)
		{
			intersection_x = (y0 - y2 + k2 * x2 - k1 * x0) / (k2 - k1);
			intersection_y = k1 * (intersection_x - x0) + y0;
		}
		else if (k1 == 0 || k2 == 0)
		{
			if (k1 == 0)
			{
				intersection_y = y0;
				intersection_x = (intersection_y - y2) / ((y3 - y2) / (x3 - x2)) + x2;
			}
			else if (k2 == 0)
			{
				intersection_y = y2;
				intersection_x = (intersection_y - y0) / ((y1 - y0) / (x1 - x0)) + x0;
			}
		}
	}
	else
	{
		if ((x1 - x0) == 0)
		{
			intersection_x = x0;
			double k2 = (y3 - y2) / (x3 - x2);
			intersection_y = k2 * (intersection_x - x2) + y2;
		}
		else if ((x3 - x2) == 0)
		{
			intersection_x = x2;
			double k1 = (y1 - y0) / (x1 - x0);
			;
			intersection_y = k1 * (intersection_x - x0) + y0;
		}
	}
}

//below finds the intersection between two 1D intervals (if they have)
void Intersection_of_1D_intervals(Array<Vec3_t> Intersection_4, Array<Vec3_t> Intersection_5, Array<Vec3_t> &Intersection_6)
{ //first two are intersection points between horizontal plane and 2nd fracture
	double x1, x2, x3, x4, x5 = 0, x6 = 0;
	x1 = Intersection_4[0](0);
	x2 = Intersection_4[1](0);
	x3 = Intersection_5[0](0);
	x4 = Intersection_5[1](0);

	double low1, up1, low2, up2;
	if (x2 >= x1)
	{
		up1 = x2;
		low1 = x1;
	}
	else
	{
		up1 = x1;
		low1 = x2;
	}
	if (x4 >= x3)
	{
		up2 = x4;
		low2 = x3;
	}
	else
	{
		up2 = x3;
		low2 = x4;
	}
	if (low1 > up2 || low2 > up1)
	{
		//printf("NULL Intersection of two 1D intervals\n");
		return;
	}
	else if (up1 >= up2)
	{
		if (low1 > low2)
		{
			x5 = low1;
			x6 = up2;
		}
		else
		{
			x5 = low2;
			x6 = up2;
		}
	}
	else if (up2 > up1)
	{
		if (low2 > low1)
		{
			x5 = low2;
			x6 = up1;
		}
		else
		{
			x5 = low1;
			x6 = up1;
		}
	}
	Intersection_6[0](0) = x5;
	Intersection_6[0](1) = Intersection_4[0](1);
	Intersection_6[1](0) = x6;
	Intersection_6[1](1) = Intersection_4[0](1);
}

//below finds intersection between line segment and a polygon, also, in 2D Cartesian system
void Intersection_between_line_segment_and_polygon(size_t &numOfIntersectionPoint, Array<Vec3_t> Verts_1, Array<Vec3_t> Intersection_infinite, Array<Vec3_t> &Intersection_2) ///there are in the same horizontal plane
{

	double end_point3_x = 0; //points between infinite line and each edge of 1st fracture
	double end_point3_y = 0;
	double end_point4_x = 0;
	double end_point4_y = 0;

	size_t index_tmp2 = 0;
	size_t index_tmp3 = 0;
	for (size_t i = 0; i < Verts_1.size(); i++)
	{
		index_tmp3 = i;
		if (i != Verts_1.size() - 1)
		{
			Line l1(Intersection_infinite);
			Line l2(Verts_1[i], Verts_1[i + 1]);
			size_t r1 = is_intersect(l1, l2);
			if (r1 == 1)
			{
				index_tmp2++;
				intersection_between_two_line_segments(Intersection_infinite, Verts_1[i], Verts_1[i + 1], end_point3_x, end_point3_y);

				break;
			}
		}
		else
		{

			Line l1(Intersection_infinite);
			Line l2(Verts_1[i], Verts_1[0]);
			size_t r1 = is_intersect(l1, l2);
			if (r1 == 1)
			{
				index_tmp2++;
				intersection_between_two_line_segments(Intersection_infinite, Verts_1[i], Verts_1[0], end_point3_x, end_point3_y);

				break;
			}
		}
	}
	for (size_t i = index_tmp3 + 1; i < Verts_1.size(); i++)
	{
		if (i != Verts_1.size() - 1)
		{
			Line l1(Intersection_infinite);
			Line l2(Verts_1[i], Verts_1[i + 1]);
			size_t r1 = is_intersect(l1, l2);
			if (r1 == 1)
			{
				index_tmp2++;
				intersection_between_two_line_segments(Intersection_infinite, Verts_1[i], Verts_1[i + 1], end_point4_x, end_point4_y);
				break;
			}
		}
		else
		{
			Line l1(Intersection_infinite);
			Line l2(Verts_1[i], Verts_1[0]);
			size_t r1 = is_intersect(l1, l2);

			if (r1 == 1)
			{
				index_tmp2++;
				intersection_between_two_line_segments(Intersection_infinite, Verts_1[i], Verts_1[0], end_point4_x, end_point4_y);
				break;
			}
		}
	}
	Intersection_2.resize(2);
	Intersection_2[0] = end_point3_x, end_point3_y, 0;
	Intersection_2[1] = end_point4_x, end_point4_y, 0;
	numOfIntersectionPoint = index_tmp2;
}

//below outputs Matlab script showing the curve (percolation parameter vs. percolation probability)
void PLotMatlab_DFN_Connectivity(char const *FileKey, Array<double> X, Array<double> Y, double alpha)
{
	ostringstream oss;
	oss << "%The alpha (power law) is " << alpha << " %\n";
	oss << "clear;\n";

	oss << "x = [";
	for (size_t i = 0; i < X.Size(); ++i)
	{
		oss << X[i] << " ";
	}
	oss << "];\n";

	oss << "y = [";
	for (size_t i = 0; i < Y.Size(); ++i)
	{
		oss << Y[i] << " ";
	}
	oss << "];\nscatter(x,y,'k');\nxlabel('Percolation parameter');\nylabel('Percolation probability');\n";

	//Open Matlab script to plot
	String fn(FileKey);
	fn.append(".m");
	ofstream of(fn.CStr(), std::ios::out);
	of << oss.str();
	of.close();
};

//below outputs data, including model edge size, P32, P32_connected, P30, and percolation-related variables, etc
void Datafile_output(char const *FileKey, double L, Array<double> P32_total, Array<double> P32_connected, Array<double> P30, Array<double> Percolation_parameter, Array<double> Ratio_of_P32, Array<double> Percolation_probability)
{
	ostringstream oss;
	size_t n = P32_total.Size();
	oss << "Side_length"
		<< "\t"
		<< "P30[i]"
		<< "\t"
		<< "P32_connected[i]"
		<< "\t"
		<< "P32_total[i]"
		<< "\t"
		<< "Percolation_parameter[i]"
		<< "\t"
		<< "Ratio_of_P32[i]"
		<< "\t"
		<< "Percolation_probability[i]"
		<< "\n";
	for (size_t i = 0; i < n; ++i)
	{
		oss << L << "\t" << P30[i] << "\t" << P32_connected[i] << "\t" << P32_total[i] << "\t" << Percolation_parameter[i] << "\t" << Ratio_of_P32[i] << "\t" << Percolation_probability[i] << "\n";
	};

	String fn(FileKey);
	fn.append(".txt");
	ofstream of(fn.CStr(), std::ios::out);
	of << oss.str();
	of.close();
};

#endif
