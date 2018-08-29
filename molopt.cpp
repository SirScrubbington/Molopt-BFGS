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
	
	////printf("Done begin orientated annealing\n");
	
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

		for(int j = 0; j < iters; j++)
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
	}
	
	////printf("Done iterative orientated annealing\n");
	
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
	
	////printf("Done centre orientated annealing\n");
	
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
		//printf("Usage: ./molopt[.exe/.out] natoms");
	
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
	
	std::string filenames[] = 
	{
		"alData.csv",
		"arData.csv",
		"acData.csv",
		"asData.csv",
		
		"blData.csv",
		"brData.csv",
		"bcData.csv",
		"bsData.csv",
		
		"clData.csv",
		"crData.csv",
		"ccData.csv",
		"csData.csv"
	};
	
	for(int i = 0; i < 6; i++)
	{
		filenames[i] = writepath + filenames[i];
	}
	
	Molecule randmol(n), circmol(n), squaremol(n);
	
	randmol.random();
	
	double circlangle = DEG2RAD * 360.0f / (n+1);
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
	
	costs[0] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[0]);
	
	end = clock();
	
	times[0] = getTime(begin,end);
		
	// random configuration
	
	begin = clock();

	molecules[1] = new Molecule(randmol);
	
	costs[1] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[1]);
	
	end = clock();
	
	times[1] = getTime(begin,end);

	// circular configuration
	
	begin = clock();
	
	molecules[2] = new Molecule(circmol);
	
	costs[2] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[2]);
	
	end = clock();
	
	times[2] = getTime(begin,end);
	
	// square configuration
	
	begin = clock();

	molecules[3] = new Molecule(squaremol);
	
	costs[3] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[3]);
	
	end = clock();
	
	times[3] = getTime(begin,end);
	
	// loop over each molecules 'niters' times, once for each every iteration

	// linear configuration
	
	begin = clock();
	
	molecules[4] = new Molecule(n);
	
	costs[4] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[4]);
	
	end = clock();
	
	times[4] = getTime(begin,end);
		
	// random configuration
		
	begin = clock();
	
	molecules[5] = new Molecule(randmol);
	
	costs[5] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[5]);
	
	end = clock();
	
	times[5] = getTime(begin,end);

	// circular configuration
	
	begin = clock();
	
	molecules[6] = new Molecule(circmol);
	
	costs[6] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[6]);
	
	end = clock();
	
	times[6] = getTime(begin,end);
	
	// square configuration
		
	begin = clock();
	
	molecules[7] = new Molecule(squaremol);
	
	costs[7] = IterativeAnnealing(n,nIterations,tempfunc,&molecules[7]);
	
	end = clock();
	
	times[7] = getTime(begin,end);
	
	// loop over each molecule 'niters' times at once, starting at the center and expanding outward
	// (starting with n/2, -> n/2 + 1, n/2 -1, n/2 + 2, ...  n/2 - n/2, n/2 + n/2

	// linear configuration
	
	begin = clock();
	
	molecules[8] = new Molecule(n);
	
	costs[8];
	
	end = clock();
	
	times[8] = getTime(begin,end);
	
	// random configuration
	
	begin = clock();
	
	molecules[9] = new Molecule(randmol);
	
	costs[9];
	
	end = clock();
	
	times[9] = getTime(begin,end);
	
	// circular configuration
	
	begin = clock();
	
	molecules[10] = new Molecule(circmol);
	
	costs[10];
	
	end = clock();
	
	times[10] = getTime(begin,end);
	
	// square configuration
		
	begin = clock();
	
	molecules[11] = new Molecule(squaremol);
	
	costs[11];
	
	end = clock();
	
	times[11] = getTime(begin,end);
	
	std::string cfolder = "coords/",subfolder;
	cfolder = cfolder + std::to_string(n);
	
	makeDirectory(cfolder.c_str());
	
	cfolder = cfolder + "/";
	
	for(int j = 0; j < 12; j++)
	{
		subfolder = cfolder + std::to_string(j) + ".csv";
		
		Molecule * mol = molecules[j];
		
		fs.open(subfolder,std::fstream::app);
		
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
					fs << "," << costs[j] << times[j] << "\n";
				}
			}
		}
		
		fs.close();
	}
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
