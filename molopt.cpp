// Filename: molopt.cpp
// Author: Scrubbs
// Date: 2018-8-7
// Description: molopt implementation file

#ifndef I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
#define I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7

#include "molopt.h"

/*

lbfgsfloatval_t * getGradientList(lbfgsfloatval_t * gradients, lbfgsfloatval_t * real, Point * coords, Matrix<lbfgsfloatval_t> * distance, int n)
{
	for(int m=0;m<n;m++)
	{
		gradients[m]=0.0f;
		for(int i=0;i<m-1;i++)
		{
			for(int j=m+1;j<n;j++)
			{
				gradients[m-1] += 
					((1.0f / pow(distance->get(i,j),14.0f)) - (1.0f / pow(distance->get(i,j),8.0f))) * 
					((coords[i].x-coords[j].x) * (coords[m].y - coords[i].y) + 
					 (coords[i].y-coords[j].y) * (coords[j].x - coords[m].x));
			}
		}
		gradients[m] *= -12.0f;
	}
	return gradients;
}

lbfgsfloatval_t lbfgsGradientCallback
(
	void * instance, 
	const lbfgsfloatval_t *x,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
)
{
	Matrix<lbfgsfloatval_t> * distance = (Matrix<lbfgsfloatval_t>*)instance;
	//g = getGradientList(g,
}

int lbfgsfunc
(
	int n,
	lbfgsfloatval_t *x,
	lbfgsfloatval_t *ptrFx,
	lbfgs_evaluate_t gradientCallback,
	lbfgs_progress_t progressCallback,
	void * instance,
	lbfgs_parameter_t *param
)
{
	return lbfgs(n,x,ptrFx,gradientCallback,progressCallback,instance,param);
}

lbfgsfloatval_t updateSystem(lbfgsfloatval_t * alphas, lbfgsfloatval_t * real, Point * coords, Matrix<lbfgsfloatval_t> * distance, int n)
{
	real = getRealAngleList(alphas,real,n);
	coords = getSystemCoordinates(real,coords,n);
	distance = getEulerianDistanceMatrix(coords,distance,n);
	return getSystemEnergy(coords,distance,n);
}

bool accept(lbfgsfloatval_t E, lbfgsfloatval_t pE, lbfgsfloatval_t temp)
{	
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0,100);
	auto rand = std::bind(distribution,generator);
	
	if (isnan(E))
		return false;
	
	if (E <= pE)
		return true;
	
	lbfgsfloatval_t P = exp((-(E-pE))/(temp*10));
	+
	printf("%f, %f, %f, %f\n",E, pE, temp, P);
	
	if((P*100) < rand())
		return true;
	return false;
}

lbfgsfloatval_t moloptAnnealingBFGS(lbfgsfloatval_t * alphas, 
						   lbfgsfloatval_t * gradients, 
						   lbfgsfloatval_t * real, 
						   Point * coords,
						   Matrix<lbfgsfloatval_t> * distance,
						   int n, 
						   int tMax
)
{

	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(-180,180);
	
	auto randAngle = std::bind(distribution,generator);

	lbfgsfloatval_t E = updateSystem(alphas,real,coords,distance,n),pE,prev;

	lbfgsfloatval_t err = E - optimal[n-2];
	
	int k = 0;
	
	for(int t = 0; t < tMax; t++)
	{
		
		lbfgsfloatval_t temp = (lbfgsfloatval_t)t/(lbfgsfloatval_t)tMax;
		
		printf("%i, %f\n",t,temp);
		
		k++;
		
		for(int j = 0; j < k; j++)
		{
			
			printf("%i, %i\n",j,k);
			
			for(int i = 1;i < n; i++)
			{
				
				printf("%i\n",i);
				
				
				pE = E;
				
				
				
				prev = alphas[i];
				alphas[i]+=(DEG2RAD*randAngle());
				E = updateSystem(alphas,real,coords,distance,n);
				
				if accept(E,pE,temp)
				{
					
				}
				else
				{
					
				}
				
			}
		}
	}

	return E;
}

lbfgsfloatval_t moloptGeneticBFGS(lbfgsfloatval_t * alphas, 
						   lbfgsfloatval_t * gradients, 
						   lbfgsfloatval_t * real, 
						   Point * coords,
						   Matrix<lbfgsfloatval_t> * distance,
						   int n,
						   int populationSiz,
						   int keepSiz
)
{
	lbfgsfloatval_t ** gAlphas = (lbfgsfloatval_t**)malloc(sizeof(*gAlphas)*populationSiz);
	lbfgsfloatval_t ** gGradients = (lbfgsfloatval_t**)malloc(sizeof(*gGradients)*populationSiz);
	lbfgsfloatval_t ** gReal = (lbfgsfloatval_t**)malloc(sizeof(*gReal)*populationSiz);
	Point ** gCoords = (Point**)malloc(sizeof(*gCoords)*populationSiz);
	
	Matrix<lbfgsfloatval_t> ** gDistance = (Matrix<lbfgsfloatval_t>**)malloc(sizeof(*gDistance)*populationSiz);
	
	gAlphas[0]=alphas;
	gGradients[0]=gradients;
	gReal[0]=real;
	gCoords[0]=coords;
	gDistance[0]=distance;
	
	lbfgsfloatval_t pE = getSystemEnergy(coords,distance,n);
	
	return pE;
}

*/

int main(int argc, char ** argv, char ** envp)
{
	int n=0;
	
	if (argc > 1)
		n = atoi(argv[1]);
	else
		printf("Usage: ./molopt[.exe/.out] natoms");
	
	if(!n)
		return 0;
	
	Molecule<lbfgsfloatval_t> * solution = new Molecule<lbfgsfloatval_t>(n);
	
	std::ofstream fs;
	fs.open("alphas.txt");
	for(int i=0;i<n;i++)
	{
		//fs << coords[i].x << "," << coords[i].y << "\n";
	}
	fs.close();
	
	//printf("%f\n",optcost);
	
	return 0;
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
