// Filename: molopt.cpp
// Author: Scrubbs
// Date: 2018-8-7
// Description: molopt implementation file

#ifndef I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
#define I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7

#include "molopt.h"

typedef double(*temperaturefunc)(double x);

std::mt19937 * gen;

double randAngle()
{
	return DEG2RAD * (((*gen)() % 180) - 180);
}

double accept()
{
	return (double)(((*gen)()%100)/100.0f);
}

double tempfunc(double x)
{
	return x * 0.001;
}

double BeginOrientatedAnnealing
(
	int n, 
	int iters,
	temperaturefunc temp,
	Molecule ** molecule
)
{
	if (*molecule == NULL)
	{
		*molecule = new Molecule(n);
	}
	
	Molecule * mol = *molecule;

	double T,E,pE,prev;
	
	int k;
	
	for(int i = 0; i < n; i++)
	{
		k = 1;
		
		for(int j = 0; j < iters; j++)
		{
			T = temp((double)iters/(double)j);
			
			for(int a = 1; a < k; a++)
			{
				pE = mol->updateSystem();
				
				double tv = mol->alphas[i];
				
				mol->alphas[i] = (double)(DEG2RAD * randAngle());

				E = mol->updateSystem();

				if(E > pE)
				{
					double delta = (E - pE);
					
					if(E < 0.0f && (accept() < exp(-delta/T)))
					{
						//
					}
					else
					{
						mol->alphas[i] = tv;
					}
				}
			}
			
			k++;
			
		}
	}
	
	return mol->updateSystem();
}

double IterativeAnnealing
(
	int n, 
	int iters,
	temperaturefunc temp,
	Molecule ** molecule
)
{
	if (*molecule == NULL)
	{
		*molecule = new Molecule(n);
	}

	Molecule * mol = *molecule;
	
	double E,pE,prev;
	
	int k = 1;
	
	for(int i = 0; i < iters; i++)
	{
		
		auto T = temp((double)n/(double)i);

		for(int j = 0; j < k; j++)
		{
			for(int a = 1; a < n; a++)
			{
				pE = mol->updateSystem();
				
				double tv = mol->alphas[i];
					
				mol->alphas[i] = (double)(DEG2RAD * randAngle());

				E = mol->updateSystem();

				if(E > pE)
				{
					double delta = (E - pE);
						
					if(E < 0.0f && (accept() < exp(-delta/T)))
					{
						//
					}
					else
					{
						mol->alphas[i] = tv;
					}
				}
			}
		}
		k++;
	}
	return mol->updateSystem();
}

std::vector<int> getTraversal(int n)
{
	int a = 0, b = n-1;
	
	std::vector<int> indexes;
	
	while(a < b)
	{
		indexes.push_back(a);
		indexes.push_back(b);
		a++;
		b--;
	}
	if(a == b)
	{
		indexes.push_back(a);
	}
	
	return indexes;
}

double centreOrientatedAnnealing
(
	int n, 
	int iters,
	temperaturefunc temp,
	Molecule ** molecule
)
{
	if (*molecule == NULL)
	{
		*molecule = new Molecule(n);
	}

	Molecule * mol = *molecule;
	
	double E,pE,prev;
	
	int k = 1, index;
	
	auto indexes = getTraversal(n);
	
	
	for(int i = n-1; i >= 0; i--)
	{
		index = indexes.at(i);
		
		printf("%i\n",index);
		
		auto T = temp((double)n/(double)i);

		for(int j = 0; j < k; j++)
		{
			for(int a = 1; a < n; a++)
			{
				pE = mol->updateSystem();
				
				double tv = mol->alphas[index];
					
				mol->alphas[index] = (double)(DEG2RAD * randAngle());

				E = mol->updateSystem();

				if(E > pE)
				{
					double delta = (E - pE);
						
					if(E < 0.0f && (accept() < exp(-delta/T)))
					{
						//
					}
					else
					{
						mol->alphas[index] = tv;
					}
				}
			}
		}
		k++;
	}
	return mol->updateSystem();
}

double getTime(clock_t begin, clock_t end)
{
	return (double)(end-begin)/(CLOCKS_PER_SEC);
}

int main(int argc, char ** argv)
{
	int n = 0;
	
	if (argc > 1)
		n = atoi(argv[1]);
	else
		printf("Usage: ./molopt[.exe/.out] natoms");
	
	if(!n)
		return 0;
	
	gen = new std::mt19937(time(0));
	
	clock_t begin,end;
	
	// Constant testing data
	
	const int nIterations = 50; // Number of annealing passes
	
	LBFGSParam<double> param; // Default parameters for the LBFGS routine
	param.max_iterations = 50; // Maximum number of LBFGS iterations
	
	Molecule * aLinear,* aRandom,* bLinear,* bRandom,* cLinear,* cRandom;
	
	double alTime, arTime, blTime, brTime, clTime, crTime;
	
	double alCost, arCost, blCost, brCost, clCost, crCost;
	
	// loop over all atoms 'niters' times at once, then progress to next atom

	// linear configuration
	
	begin = clock();
	
	aLinear = new Molecule(n);
	
	end = clock();
	
	alTime = getTime(begin,end);
	
	// random configuration
	
	begin = clock();
	
	aRandom = new Molecule(n);
	aRandom->random();
	
	end = clock();
	
	arTime = getTime(begin,end);
	
	// loop over each molecules 'niters' times, once for each every iteration

	// linear configuration
	
	begin = clock();
	
	bLinear = new Molecule(n);
	
	end = clock();
	
	blTime = getTime(begin,end);
	
	// random configuration
		
	begin = clock();
	
	bRandom = new Molecule(n);
	bRandom->random();
	
	end = clock();
	
	brTime = getTime(begin,end);
		
	// loop over each molecule 'niters' times at once, starting at the center and expanding outward
	// (starting with n/2, -> n/2 + 1, n/2 -1, n/2 + 2, ...  n/2 - n/2, n/2 + n/2

	// linear configuration
	
	begin = clock();
	
	cLinear = new Molecule(n);
	
	clCost = centreOrientatedAnnealing(n,nIterations,tempfunc,&cLinear);
	
	end = clock();
	
	clTime = getTime(begin,end);
	
	// random configuration
	
	begin = clock();
	
	cRandom = new Molecule(n);
	cRandom->random();
	
	crCost = centreOrientatedAnnealing(n,nIterations,tempfunc,&cRandom);
	
	end = clock();
	
	crTime = getTime(begin,end);
	
	printf("%f, %f\n",clCost,crCost);
	
	return 0;
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
