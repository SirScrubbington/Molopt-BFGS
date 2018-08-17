// Filename: Molecule.h
// Author: Scrubbs
// Date: 2018-8-17
// Description: Molecule.h header file

#ifndef I_v46cQspgT8KK6epE1iB7GAk686t4t
#define I_v46cQspgT8KK6epE1iB7GAk686t4t

#include "Matrix.h"
#include "Point.h"

#define _USE_MATH_DEFINES

#include <cmath>

#ifndef DEG2RAD
	#define DEG2RAD (M_PI/180.0f)
#endif // DEG2RAD

#ifndef RAD2DEG
	#define RAD2DEG (180.0f/M_PI)
#endif // RAD2DEG

#include <cstdio>

template<typename moltyp> 

class Molecule
{
private:
public:
	// Data Structures
	Matrix<moltyp> * distance;
	
	moltyp * alphas;
	moltyp * real;
	moltyp * gradients; 
	
	Point * coords;
	
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
	
	moltyp * getRealAngleList()
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
			coords[i] = coords[i-1] + getCoords(1.0f,real[i-1]);
		}	
		
		return coords;
	}

	moltyp getR(moltyp x, moltyp y)
	{
		return sqrt(pow(x,2) + pow(y,2));
	}

	Matrix<moltyp> * getEulerianDistanceMatrix()
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

	moltyp getSystemEnergy()
	{
		moltyp V = 0;
		
		for(int j=0;j<n;j++)
			for(int i=0;i<j;i++)
			{
				V = V + ((1.0f/pow(distance->get(i,j),12.0f))-(2.0f/pow(distance->get(i,j),6.0f)));
			}
		
		return V;
		
	}
		
	Molecule(int n)
	{
		distance = new Matrix<moltyp>(n,n);
		
		alphas = (moltyp*)malloc(sizeof(*alphas)*n);
		real = (moltyp*)malloc(sizeof(*real)*n);
		gradients = (moltyp*)malloc(sizeof(*gradients)*n); 
		
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
