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
	int n = (*gen)() % 180;
	return (double)n*DEG2RAD;
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
	
	for(int i = 0; i < n; i++)
	{
		T = temp((double)n/(double)i);
		
		for(int j = 0; j < iters; j++)
		{
			for(int a = 1; a < iters; a++)
			{
				pE = mol->updateSystem();
				
				double tv = mol->alphas[i];
				
				mol->alphas[i] = randAngle();
				
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
	}
	
	//printf("Done begin orientated annealing\n");
	
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
				
				double tv = mol->alphas[a];
					
				mol->alphas[a] = randAngle();

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
						mol->alphas[a] = tv;
					}
				}
			}
		}
		k++;
	}
	
	//printf("Done iterative orientated annealing\n");
	
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
		
		auto T = temp((double)n/(double)i);

		for(int j = 0; j < iters; j++)
		{
			for(int a = 1; a < iters; a++)
			{
				pE = mol->updateSystem();
				
				double tv = mol->alphas[index];

				mol->alphas[index] = randAngle();

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
	}
	
	//printf("Done centre orientated annealing\n");
	
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
	
	std::fstream fs;
	
	clock_t begin,end;
	
	// Constant testing data
	
	const int nIterations = 50; // Number of annealing passes
	
	LBFGSParam<double> param; // Default parameters for the LBFGS routine
	param.max_iterations = 50; // Maximum number of LBFGS iterations
	
	Molecule * aLinear,* aRandom,* bLinear,* bRandom,* cLinear,* cRandom;
	
	std::string directory = "data/";
	std::string subdir = std::to_string(n);
	std::string writepath = directory + subdir + "/";
	
	#if defined(_WIN32)
		_mkdir(writepath.c_str()); // can be used on Windows
	#else 
		mkdir(writepath.c_str(),0733); // can be used on non-Windows
	#endif
	
	std::string filenames[] = {"alData.csv","arData.csv","blData.csv","brData.csv","clData.csv","crData.csv"};
	for(int i = 0; i < 6; i++)
	{
		filenames[i] = writepath + filenames[i];
	}
	
	//double alTime, arTime, blTime, brTime, clTime, crTime;
	
	//double alCost, arCost, blCost, brCost, clCost, crCost;
	
	Molecule randmol(n);
	randmol.random();
	
	Molecule * molecules[6];
	double times[6] = {0.0f};
	double costs[6] = {0.0f};
	
	// loop over all atoms 'niters' times at once, then progress to next atom

	// linear configuration
	
	begin = clock();
	
	molecules[0] = new Molecule(n);
	
	costs[0] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[0]);
	
	end = clock();
	
	times[0] = getTime(begin,end);
	
	// random configuration
	
	begin = clock();

	molecules[1] = new Molecule(randmol);
	
	costs[1] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[1]);
	
	end = clock();
	
	times[1] = getTime(begin,end);
	
	// loop over each molecules 'niters' times, once for each every iteration

	// linear configuration
	
	begin = clock();
	
	molecules[2] = new Molecule(n);
	
	costs[2] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[2]);
	
	end = clock();
	
	times[2] = getTime(begin,end);
	
	// random configuration
		
	begin = clock();
	
	molecules[3] = new Molecule(randmol);
	
	costs[3] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[3]);
	
	end = clock();
	
	times[3] = getTime(begin,end);
		
	// loop over each molecule 'niters' times at once, starting at the center and expanding outward
	// (starting with n/2, -> n/2 + 1, n/2 -1, n/2 + 2, ...  n/2 - n/2, n/2 + n/2

	// linear configuration
	
	begin = clock();
	
	molecules[4] = new Molecule(n);
	
	costs[4] = centreOrientatedAnnealing(n,nIterations,tempfunc,&molecules[4]);
	
	end = clock();
	
	times[4] = getTime(begin,end);
	
	// random configuration
	
	begin = clock();
	
	molecules[5] = new Molecule(randmol);
	
	costs[5] = centreOrientatedAnnealing(n,nIterations,tempfunc,&molecules[5]);
	
	end = clock();
	
	times[5] = getTime(begin,end);
	
	for(int j = 0; j < 6; j++)
	{
		Molecule * mol = molecules[j];
		
		//fs.open();
		
		//fs.close();
		
		fs.open(filenames[j],std::fstream::app);
		
		if(fs)
		{
			for (int i=0;i<n;i++)
			{
				
				fs << mol->alphas[i];
				
				if(i<n-1)
				{
					fs << ",";
				}
				else
				{
					fs << "," << costs[j] << "\n";
				}
			}
		}
		
		fs.close();
	}
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
