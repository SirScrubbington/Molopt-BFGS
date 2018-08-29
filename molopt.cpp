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

double tempfuncanneal(double x)
{
	return x * 20.0f;
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
	
	for(int i = 1; i < n; i++)
	{
		T = temp((double)n/(double)i);
		
		for(int j = 0; j < iters; j++)
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
	
	//int k = 1;
	
	for(int i = 0; i < iters; i++)
	{
		
		auto T = temp((double)n/(double)i);

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
	
	LBFGSSolver<double> solver(param);
	
	std::string directory = "data/";
	std::string subdir = std::to_string(n);
	std::string writepath = directory + subdir + "/";
	
	makeDirectory(writepath.c_str());
	
	Molecule circmol(n), squaremol(n);
	
	double circlangle = DEG2RAD * 360.0f / n;
	double squarangle = DEG2RAD * 90.0f;
	
	for(int i = 0; i < n; i++)
	{
		circmol.alphas[i] = circlangle;
	}
	
	if(n / 3 > 0)
	{
		for(int i = (n/3); i < n; i += (n / 3))
		{
			if (i >= n) break;
			squaremol.alphas[i] = squarangle;
		}
	}
	
	Molecule * molecules[12];
	double times[12] = {0.0f};
	double costs[12] = {0.0f};
	
	// loop over all atoms 'niters' times at once, then progress to next atom

	// linear configuration
	
	begin = clock();
	
	molecules[0] = new Molecule(n);

	//costs[0] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[0]);
	
	end = clock();

	times[0] = getTime(begin,end);
	
	printf("Time 0: %f\n",times[0]);

	// circular configuration
	
	begin = clock();
	
	molecules[1] = new Molecule(circmol);
	molecules[1]->updateSystem();
	
	//costs[1] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[1]);
	
	end = clock();
	
	times[1] = getTime(begin,end);
	
	printf("Time 2: %f\n",times[1]);
	
	// square configuration
	
	begin = clock();

	molecules[2] = new Molecule(squaremol);
	molecules[2]->updateSystem();
	//costs[2] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[2]);
	
	end = clock();
	
	times[2] = getTime(begin,end);
	
	printf("Time 3: %f\n",times[2]);
	
	// loop over each molecules 'niters' times, once for each every iteration

	// linear configuration
	
	begin = clock();
	
	molecules[3] = new Molecule(n);
	molecules[3]->updateSystem();
	
	//costs[3] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[3]);
	
	end = clock();
	
	times[3] = getTime(begin,end);
	
	printf("Time 0: %f\n",times[3]);

	// circular configuration
	
	begin = clock();
	
	molecules[4] = new Molecule(circmol);
	molecules[4]->updateSystem();
	
	//costs[4] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[4]);
	
	end = clock();
	
	times[4] = getTime(begin,end);
	
	printf("Time 0: %f\n",times[4]);
	
	// square configuration
		
	begin = clock();
	
	molecules[5] = new Molecule(squaremol);
	molecules[5]->updateSystem();
	
	//costs[5] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[5]);
	
	end = clock();
	
	times[5] = getTime(begin,end);
	
	printf("Time 0: %f\n",times[5]);
	
	/*
	
	-- BFGS METHOD --
	
	*/
	
	begin = clock();
	
	molecules[6] = new Molecule(&(molecules[0]));
	
	costs[6]// = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[6]);
	
	end = clock();
	
	times[6] = getTime(begin,end);
	
	printf("Time 6: %f\n",times[6]);

	// circular configuration
	
	begin = clock();
	
	molecules[7] = new Molecule(&(molecules[1]));
	
	costs[7]// = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[7]);
	
	end = clock();
	
	times[7] = getTime(begin,end);
	
	printf("Time 7: %f\n",times[7]);
	
	// square configuration
	
	begin = clock();

	molecules[8] = new Molecule(&(molecules[2]));
	
	costs[8]// = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[8]);
	
	end = clock();
	
	times[8] = getTime(begin,end);
	
	printf("Time 8: %f\n",times[8]);
	
	// loop over each molecules 'niters' times, once for each every iteration

	// linear configuration
	
	begin = clock();
	
	molecules[9] = new Molecule(&(molecules[3]));
	
	costs[9];
	
	end = clock();
	
	times[9] = getTime(begin,end);
	
	printf("Time 9: %f\n",times[9]);

	// circular configuration
	
	begin = clock();
	
	molecules[10] = new Molecule(&(molecules[4]));
	
	costs[10];
	
	end = clock();
	
	times[10] = getTime(begin,end);
	
	printf("Time 0: %f\n",times[10]);
	
	// square configuration
		
	begin = clock();
	
	molecules[11] = new Molecule(&(molecules[5]));
	
	costs[11];
	
	end = clock();
	
	times[11] = getTime(begin,end);
	
	printf("Time 0: %f\n",times[11]);
	
	begin = clock();
	
	std::string dFolder = "data/",dSubfolder;
	std::string cFolder = "coords/",cSubfolder;
	
	dFolder = dFolder + std::to_string(n);
	cFolder = cFolder + std::to_string(n);
	
	makeDirectory(cFolder.c_str());
	makeDirectory(dFolder.c_str());
	
	cFolder = cFolder + "/";
	dFolder = dFolder + "/";
	
	for(int j = 0; j < 6; j++)
	{
		Molecule * mol = molecules[j];
		
		cSubfolder = cFolder + std::to_string(j) + ".csv";
		
		fs.open(cSubfolder,std::fstream::app);
		
		fs << j << "," << mol->n << ",";
		
		for(int i = 0; i < n; i++)
		{
			fs << mol->coords.at(i).x << ":" << mol->coords.at(i).y;
			
			if(i <n-1)
				fs << ",";
			
			else
				fs << "," << optimal[n-2] << "," << costs[j] << "," << times[j] << "\n";
		}
		
		fs.close();
		
		dSubfolder = dFolder + std::to_string(j) + ".csv";
		
		fs.open(dSubfolder,std::fstream::app);
		
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
					fs << "," << costs[j] << times[j] << "\n";
				}
			}
		}
		
		fs.close();
	}
	
	end = clock();
	
	printf("Time to write: %f seconds.\n",getTime(begin,end));
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
