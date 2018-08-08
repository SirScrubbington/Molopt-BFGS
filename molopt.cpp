// Filename: molopt.cpp
// Author: Scrubbs
// Date: 2018-8-7
// Description: molopt implementation file

#ifndef I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
#define I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7

#include "molopt.h"

double * getGradientList(double * gradients, double * real, Point * coords, Matrix<double> * distance, int n)
{
	for(int m=0;m<n;m++)
	{
		gradients[m]=0.0f;
		for(int i=0;i<m;i++)
		{
			for(int j=m+1;j<n;j++)
			{
				gradients[m] += 
				((1.0f / pow(distance->get(i,j),14)) - (1.0f / pow(distance->get(i,j),8))) * 
				((coords[i].x-coords[j].x) * (coords[m].y - coords[i].y) + 
				 (coords[i].y-coords[j].y) * (coords[j].x - coords[m].x));
			}
		}
		gradients[m] *= -12.0f;
	}
	return gradients;
}

double * getRealAngleList(double * alphas, double * real, int n)
{
	real[0] = alphas[0];
	
	for(int i=1;i<n;i++)
		real[i] = real[i-1] + alphas[i];
	
	return real;
}

Point getCoords(double distance, double angle)
{
	return {distance*cos(DEG2RAD*angle),distance*sin(DEG2RAD*angle)};
}

Point * getSystemCoordinates(double * real, Point * coords, int n)
{
	for(int i=2;i<n;i++)
	{
		coords[i] = coords[i-1] + getCoords(1.0f,real[i-1]);
	}	
	
	return coords;
}

double getR(double x, double y)
{
	return sqrt(pow(x,2) + pow(y,2));
}

Matrix<double> * getEulerianDistanceMatrix(Point * coords, Matrix<double> * distance, int n)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			//i == j ? distance.set(i,j,0.0f) : distance.set(i,j,getR());
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

double getSystemEnergy(Point * coords, Matrix<double> * distance, int n)
{
	double V = 0;

	for(int i=0;i<n;i++)
	{
		for(int j=i;j<n;j++)
		{
			
			if (i==j)
				continue;
			
			//Point distance = (coords[i] - coords[j]).absolute();
			
			//double r = getR(distance.x,distance.y);
			
			V = V + ((1.0f/pow(distance->get(i,j),12.0f))-(2.0f/pow(distance->get(i,j),6.0f)));
		}
	}
	
	return V;
	
}

int main(int argc, char ** argv, char ** envp)
{
	int n=0;
	
	if (argc > 1)
		n = atoi(argv[1]);
	else
		printf("Usage: ./molopt[.exe/.out] natoms");
	
	if(!n)
		return 0;
	
	Matrix<double> * distance = new Matrix<double>(n,n);
	
	double * alphas = (double*)malloc(sizeof(*alphas)*n);
	
	double * real = (double*)malloc(sizeof(*real)*n);
	
	double * gradients = (double*)malloc(sizeof(*gradients)*n); 
	
	Point * coords = (Point*)malloc(sizeof(*coords)*n);
	
	for(int i=0;i<n;i++)
	{
		alphas[i]=0.0f;
		gradients[i]=0.0f;
		coords[i]={0.0f,0.0f};
	}
	
	coords[1]={1.0f,0.0f};
	
	//for(int i=0;i<n;i++)
	//	printf("%f\n",alphas[i]);
	//printf("\n");
	
	real = getRealAngleList(alphas,real,n);
	
	//for(int i=0;i<n;i++)
	//	printf("%f\n",real[i]);
	//printf("\n");
	
	coords = getSystemCoordinates(real,coords,n);
	
	//for(int i=0;i<n;i++)
	//	printf("(%f %f)\n",coords[i].x,coords[i].y);
	//printf("\n");
	
	distance = getEulerianDistanceMatrix(coords,distance,n);
	
	//for(int i=0;i<distance->y; i++)
	//{
	//	for(int j=0;j<distance->x;j++)
	//	{
	//		printf("%f,",distance->get(i,j));
	//	}
	//	printf("\n");
	//}
	
	double E = getSystemEnergy(coords,distance,n);
	
	gradients = getGradientList(gradients,real,coords,distance,n);
	
	for(int i=0;i<n;i++)
		printf("(%i , %f)\n",i,gradients[i]);
	printf("\n");
	
	printf("%f\n",E);
	
	
	{	

			real = getRealAngleList(alphas,real,n);

			coords = getSystemCoordinates(real,coords,n);

			distance = getEulerianDistanceMatrix(coords,distance,n);
			
			double E = getSystemEnergy(coords,distance,n);
			
			gradients = getGradientList(gradients,real,coords,distance,n);
	}
	return 0;
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
