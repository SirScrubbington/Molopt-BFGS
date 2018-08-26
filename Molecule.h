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

#include <cmath>

#ifndef DEG2RAD
	#define DEG2RAD
	const double D2R = (M_PI/180.0f);
#endif // DEG2RAD

#ifndef RAD2DEG
	#define RAD2DEG
	const double R2D = (180.0f/M_PI);
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
		
	VectorXd * alphas;
	VectorXd * real;
	VectorXd * gradients;
	
	std::vector<Point> * coords;
	
	int n;
	
	// Member Functions

	void printGradients()
	{
		for(int i=0;i<n;i++)
			printf("(%i , %f)\n",i,gradients[i]);
		printf("\n");
	}

	void printAlphas()
	{
		for(int i=0;i<n;i++)
			printf("%f\n",alphas[i]);
		printf("\n");
	}
		
	void printReal()
	{
		for(int i=0;i<n;i++)
			printf("%f\n",real[i]);
		printf("\n");
	}
	
	void printCoords()
	{
		for(int i=0;i<n;i++)
			printf("(%f %f)\n",coords[i].x,coords[i].y);
		printf("\n");
	}
	
	void printDistance()
	{
		for(int i=0;i<distance->y; i++)
		{
			for(int j=0;j<distance->x;j++)
			{
				printf("%f,",distance->get(i,j));
			}
			printf("\n");
		}
	}
	
	 * getRealAngleList()
	{
		real[0] = alphas[0];
		
		for(int i=1;i<n;i++)
		{
			real[i] = real[i-1] + alphas[i];
			if(real[i] >= DEG2RAD*360)
				real[i] -= DEG2RAD*360;
		}
		return real;
	}

	Point getCoords(double distance, double angle)
	{
		return {distance*cos(angle),distance*sin(angle)};
	}

	Point * getSystemCoordinates()
	{
		for(int i=2;i<n;i++)
		{
			//printf("%i\n",i);
			coords[i] = coords[i-1] + getCoords(1.0f,real[i-1]);
		}	
		
		return coords;
	}

	double getR(double x, double y)
	{
		return sqrt(pow(x,2) + pow(y,2));
	}

	Matrix<double> * getEulerianDistanceMatrix()
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
					Point d = (coords[i] - coords[j]).absolute();
					distance->set(i,j,getR(d.x,d.y));
				}
			}
		}
		return distance;
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
		
		alphas = (double*)malloc(sizeof(*alphas)*n);
		real = (double*)malloc(sizeof(*real)*n);
		gradients = (double*)malloc(sizeof(*gradients)*n); 
		
		coords = (Point*)malloc(sizeof(*coords)*n);
		
		for(int i=0;i<n;i++)
		{
			alphas[i]=0.0f;
			gradients[i]=0.0f;
			coords[i]={0.0f,0.0f};
		}
		
		coords[1]={1.0f,0.0f};
		
		real = getRealAngleList();
		
		coords = getSystemCoordinates();
		
		distance = getEulerianDistanceMatrix();
		
	}
	
	double updateSystem()
	{
		real = getRealAngleList();
		coords = getSystemCoordinates();
		distance = getEulerianDistanceMatrix();
		return getSystemEnergy();
	}
	
	// Check entire system for crossovers
	Point hasCrossovers()
	{
		
		Point * a, * b;
		
		for(int i = 0;i < n-1;i++)
		{
			for(int j = 0; j < n-1; j++)
			{
				if(checkCrossover(i,j))
				{
					return {(double)i,(double)j};
				}
			}
		}
		return {-1.0f,-1.0f};
	}
	
	// Check if crossover exists at index
	bool checkCrossover(int a, int b)
	{
		Point p1, p2, q1, q2;
		
		p1 = coords[a];
		p2 = coords[a+1];
		
		q1 = coords[b];
		q2 = coords[b+1];
		
		return (
			((q1.x-p1.x)*(p2.y-p1.y) - (q1.y-p1.y)*(p2.x-p1.x))
				* ((q2.x-p1.x)*(p2.y-p1.y) - (q2.y-p1.y)*(p2.x-p1.x)) < 0)
				&&
			(((p1.x-q1.x)*(q2.y-q1.y) - (p1.y-q1.y)*(q2.x-q1.x))
				* ((p2.x-q1.x)*(q2.y-q1.y) - (p2.y-q1.y)*(q2.x-q1.x)) < 0);
		
	}
	
	Molecule(const Molecule<double> &cpy)
	{
		
		n = cpy.n;
		
		distance = new Matrix<double>(n,n);
		
		alphas = (double*)malloc(sizeof(*alphas)*n);
		real = (double*)malloc(sizeof(*real)*n);
		gradients = (double*)malloc(sizeof(*gradients)*n); 
		
		coords = (Point*)malloc(sizeof(*coords)*n);
		
		for(int i = 0 ; i < n ; i++)
		{
			alphas[i] = cpy.alphas[i];
			real[i] = cpy.real[i];
			gradients[i] = cpy.gradients[i];
			coords[i] = cpy.coords[i];
			
			for(int j = 0 ; j < n ; j++)
			{
				distance->set(i,j,cpy.distance->get(i,j));
			}
		}
	}
		
	~Molecule()
	{
		free(distance);
		free(alphas);
		free(real);
		free(gradients);
		free(coords);
	}
};

#endif // I_v46cQspgT8KK6epE1iB7GAk686t4t
