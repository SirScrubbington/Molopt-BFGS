// Filename: Molecule.h
// Author: Scrubbs
// Date: 2018-8-17
// Description: Molecule.h header file

#ifndef I_v46cQspgT8KK6epE1iB7GAk686t4t
#define I_v46cQspgT8KK6epE1iB7GAk686t4t

#include "Matrix.h"
#include "Point.h"

#define NO_CROSSOVERS 0
#define ALLOW_CROSSOVERS 1

#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>

#ifndef D2R
	#define D2R
	const double DEG2RAD = M_PI/180.0f;
#endif // DEG2RAD

#ifndef R2D
	#define R2D
	const double RAD2DEG = 180.0f/M_PI;
#endif // RAD2DEG

#include <cstdio>

#include "include/Eigen/Core"
#include "include/LBFGS.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

using namespace LBFGSpp;

class Molecule
{
private:
public:
	// Data Structures
	
	Matrix<double> * distance;
		
	VectorXd alphas;
	VectorXd real;
	
	std::vector<Point> coords;
	
	int n;
	
	// Member Functions
	void random()
	{
		double valid = false;
		
		while(!valid)
		{
			for(int i = 0; i < n; i++)
			{
				alphas[i] = DEG2RAD * ((rand() % 180) - 180);
			}
			if (updateSystem() < 0.0f)
			{
				valid = true;
			}
		}
	}
	
	void getRealAngleList()
	{
		real[0] = alphas[0];
		
		for(int i=1;i<n;i++)
		{
			real[i] = real[i-1] + alphas[i];
			if(real[i] >= DEG2RAD*360)
				real[i] -= DEG2RAD*360;
		}
	}

	Point getCoords(double distance, double angle)
	{
		return {distance*cos(angle),distance*sin(angle)};
	}

	void getSystemCoordinates()
	{
		for(int i=2;i<n;i++)
		{
			//printf("%i\n",i);
			coords.at(i) = coords.at(i-1) + getCoords(1.0f,real[i-1]);
		}	
	}

	double getR(double x, double y)
	{
		return sqrt(pow(x,2) + pow(y,2));
	}

	void getEulerianDistanceMatrix()
	{
		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				if (i==j)
				{
					distance->set(i,j,0.0f);
				}
				else
				{
					Point d = (coords.at(i) - coords.at(j)).absolute();
					distance->set(i,j,getR(d.x,d.y));
				}
			}
		}
	}

	double getSystemEnergy()
	{
		double V = 0;
		
		for(int j=0;j<n;j++)
			for(int i=0;i<j;i++)
			{
				V = V + ((1.0f/pow(distance->get(i,j),12.0f))-(2.0f/pow(distance->get(i,j),6.0f)));
			}
		
		return V;
		
	}
		
	Molecule(int n)
	{
		this->n = n;
		
		distance = new Matrix<double>(n,n);
		
		alphas = VectorXd::Zero(n);
		real = VectorXd::Zero(n);
		
		coords = std::vector<Point>(n);
		
		coords.at(0) = {0.0f,0.0f};
		coords.at(1) = {1.0f,0.0f};
		
		updateSystem();

	}
	
	double updateSystem()
	{
		getRealAngleList();
		getSystemCoordinates();
		getEulerianDistanceMatrix();
		return getSystemEnergy();
	}
	
	double operator()(const VectorXd& x, VectorXd& grad)
	{
		double a, b, c, d;
		
		int g_n = grad.size();
		int x_n = x.size();
		
		//printf("yeet\n");
		alphas = x;
		
		auto cost = updateSystem();
		
		for(int m = 0; m < n; m++)
		{
			grad[m] = 0.0f;
			
			for(int i = 0; i < m - 1; i++)
			{
				for(int j = m + 1; j < n; j++)
				{
					a = 1.0f / pow(distance->get(i,j),14.0f);
					b = 1.0f / pow(distance->get(i,j),8.0f);
					
					c = ((coords[i].x - coords[j].x) * (coords[m].y - coords[i].y));
					d = ((coords[i].y - coords[j].y) * (coords[j].x - coords[m].x));
					
					grad[m] += (a - b) * (c + d);
				}
			}
			grad[m] *= -12.0f;
		}	
		
		return cost;
	}
};

#endif // I_v46cQspgT8KK6epE1iB7GAk686t4t
