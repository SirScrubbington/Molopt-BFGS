// Filename: molopt.cpp
// Author: Scrubbs
// Date: 2018-8-7
// Description: molopt implementation file

#ifndef I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
#define I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7

#include "molopt.h"

typedef double(*temperaturefunc)(double x);

lbfgsfloatval_t cfdd
(
	lbfgsfloatval_t fa,
	lbfgsfloatval_t fb,
	lbfgsfloatval_t h
)
{
	return ((fa) - (fb))/h;
}

lbfgsfloatval_t gradientCallbackCFDD
(
	void * instance,
	const lbfgsfloatval_t *x,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
)
{
	Molecule<lbfgsfloatval_t> * mol = (Molecule<lbfgsfloatval_t>*)instance;
	
	lbfgsfloatval_t a, b, fa, fb,temp,diffstep,cost;
	
	for(int i=0;i<n;i++)
	{
		//printf("%f, %f\n",x[i],mol->alphas[i]);
		mol->alphas[i] = x[i];
	}
	cost = mol->updateSystem();
	
	diffstep = 0.1;
	
	for(int i=1;i<n;i++)
	{
		temp = mol->alphas[i];
		
		a = mol->alphas[i] - diffstep;
		b = mol->alphas[i] + diffstep;
		
		mol->alphas[i] = a;
		
		fa = mol->updateSystem();
		
		mol->alphas[i] = b;
		
		fb = mol->updateSystem();
		
		mol->alphas[i] = temp;
		
		mol->gradients[i] = cfdd(fa,fb,diffstep);
		
		g[i] = mol->gradients[i];
		
		//printf("%f\n",g[i]);
		
	}
	
	//printf("\n");
	
	return cost;
}

lbfgsfloatval_t gradientCallback
(
	void * instance,
	const lbfgsfloatval_t *x,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
)
{
	Molecule<lbfgsfloatval_t> * mol = (Molecule<lbfgsfloatval_t>*)instance;
	
	for(int i=0;i<n;i++)
		mol->alphas[i] = x[i];
	
	auto cost = mol->updateSystem();
	
	auto coords = mol->coords;
	auto distance = mol->distance;
	auto gradients = mol->gradients;
	
	for(int m = 0; m < n; m++)
	{
		gradients[m] = 0.0f;
		
		for(int i = 0; i < m - 1;i++)
		{
			for(int j = m + 1; j < n; j++)
			{
				gradients[m] += 
				((1.0f / pow(distance->get(i,j),14.0f)) - (1.0f / pow(distance->get(i,j),8.0f))) * 
				(((coords[i].x-coords[j].x) * (coords[m].y - coords[i].y)) + 
				 ((coords[i].y-coords[j].y) * (coords[j].x - coords[m].x)));
			}
		}
		gradients[m] *= -12.0f;
	}
	
	//mol->printGradients();
	
	//g = gradients;

	for(int i=0;i<n;i++)
		g[i] = gradients[i];
	
	return cost;
}

int progressCallback
(
	void * instance,
	const lbfgsfloatval_t *x,
	const lbfgsfloatval_t *g,
	const lbfgsfloatval_t fx,
	const lbfgsfloatval_t xnorm,
	const lbfgsfloatval_t gnorm,
	const lbfgsfloatval_t step,
	int n,
	int k,
	int ls
)
{
	
}

lbfgsfloatval_t moloptAnnealingBFGS
(
	int n, 
	int niters,
	temperaturefunc tempfunc,
	Molecule<lbfgsfloatval_t> ** save
)
{
	
	std::mt19937 gen(time(0));
	
	std::uniform_int_distribution<int> dist(-180,180);
	
	std::uniform_real_distribution<double> acpt(0.0f,1.0f);
	
	auto randAngle = std::bind(dist,gen);
	auto accept = std::bind(acpt,gen);
	
	lbfgs_parameter_t param;
	lbfgs_parameter_init(&param);
	
	Molecule<lbfgsfloatval_t> * mol = new Molecule<lbfgsfloatval_t>(n);
	
	int k = 0;
	
	int nAccepted=0,nRejected=0;
	
	lbfgsfloatval_t E,pE,prev;
	
	for(int i = 0;i < niters;i++)
	{
		k++;
		
		auto T = tempfunc((double)niters/(double)i);
		
		//printf("%f\n",T);

		for(int j = 0;j < k;j++)
		{
			for(int a = 2; a < n; a++)
			{
				pE = mol->updateSystem();
				
				double tv = mol->alphas[i];
				mol->alphas[i] = (double)(DEG2RAD * randAngle());

				E = mol->updateSystem();
				
				// Annealing Algorithm
				
				if(E > pE)
				{
					auto delta = E - pE;
					
					if(accept() < exp(-delta/T))
					{
						// accept
						nAccepted++;
					}
					else
					{
						mol->alphas[i] = tv;
						nRejected++;
					}
				}
				else
				{
					// accept
					//nAccepted++;
				}
			}
		}
	}
	
	printf("Accepted: %i\nRejected: %i\n",nAccepted,nRejected);
	
	//printf("%f\n",mol->updateSystem());
	
	lbfgsfloatval_t * bfgsval = (lbfgsfloatval_t*)malloc(sizeof(lbfgsfloatval_t)*n);
	
	lbfgs(mol->n,mol->alphas,bfgsval,gradientCallback,progressCallback,mol,&param);
	
	*save = mol;
	
	return mol->updateSystem();
}

lbfgsfloatval_t moloptGeneticBFGS
(
	int n
)
{
	
}

double tempfunc(double x)
{
	return x * 0.001;
}

int main(int argc, char ** argv, char ** envp)
{
	#if     defined(USE_SSE) && defined(__SSE2__) && LBFGS_FLOAT == 64
		//printf("Use SSE2 optimization for 64bit double precision.\n");
		
	#elif   defined(USE_SSE) && defined(__SSE__) && LBFGS_FLOAT == 32
		//printf("Use SSE optimization for 32bit float precision.\n");
		
	#else
		//printf("No CPU specific optimization. \n");
	
	#endif
	
	int n=0;
	
	if (argc > 1)
		n = atoi(argv[1]);
	else
		printf("Usage: ./molopt[.exe/.out] natoms");
	
	if(!n)
		return 0;
	
	Molecule<lbfgsfloatval_t> * mol;
	
	//printf("%i\n",n);
	
	int nIters = 50;
	
	clock_t begin,end;
	
	begin = clock();
	
	lbfgsfloatval_t optcost = moloptAnnealingBFGS(n,nIters,tempfunc,&mol);

	end = clock();	
	
	std::ofstream fs;
	fs.open("alphas.txt");
	
	for(int i=0;i<n;i++)
	{
		fs << mol->coords[i].x << "," << mol->coords[i].y << "\n";
	}
	fs.close();
	double err = fabs((optimal[n-2]-optcost)/optimal[n-2])*100;
	double time = (double)(end - begin)/CLOCKS_PER_SEC;
	
	
	
	printf("Atoms\tOptimal\tFound\tError\tGen\tTime\n");
	printf("%i\t%2.2f\t%2.2f\t%2.2f\t%i\t%2.2f\n",n,optimal[n-2],optcost,err,nIters,time);
	return 0;
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
