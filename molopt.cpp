// Filename: molopt.cpp
// Author: Scrubbs
// Date: 2018-8-7
// Description: molopt implementation file

#ifndef I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
#define I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7

#include "molopt.h"

typedef double(*temperaturefunc)(double x);

double tempfunc(double x)
{
	return x * 0.001;
}

double moloptAnnealingBFGS
(
	int n, 
	int nIters,
	temperaturefunc tempfunc,
	Molecule ** save
)
{
	
	std::mt19937 gen(time(0));
	
	std::uniform_int_distribution<int> dist(-180,180);
	
	std::uniform_real_distribution<double> acpt(0.0f,1.0f);
	
	auto randAngle = std::bind(dist,gen);
	auto accept = std::bind(acpt,gen);
	
	lbfgs_parameter_t param;
	lbfgs_parameter_init(&param);
	
	Molecule * mol = new Molecule(n);
	
	int k = 0;
	
	int nAccepted=0,nRejected=0;
	
	double E,pE,prev;
	
	for(int i = 0;i < n;i++)
	{
		k++;
		
		auto T = tempfunc((double)n/(double)i);

		for(int j = 0;j < nIters;j++)
		{
			for(int a = 1; a < n; a++)
			{
				//printf("%i %i %i\n",i,j,a);
				pE = mol->updateSystem();
				
				double tv = mol->alphas[i];
				mol->alphas[i] = (double)(DEG2RAD * randAngle());

				E = mol->updateSystem();

				if(E > pE)
				{
					auto delta = E - pE;
					
					if(accept() < exp(-delta/T))
					{
						nAccepted++;
					}
					else
					{
						mol->alphas[i] = tv;
						nRejected++;
					}
				}
			}
		}
	}
	
	//printf("Accepted: %i\nRejected: %i\n",nAccepted,nRejected);

	*save = mol;
	
	return mol->updateSystem();
}

int main(int argc, char ** argv, char ** envp)
{
	int n = 0;
	
	if (argc > 1)
		n = atoi(argv[1]);
	else
		printf("Usage: ./molopt[.exe/.out] natoms");
	
	if(!n)
		return 0;
	
	Molecule * mol;

	int nIters = 50;
	
	clock_t begin,end;
	
	begin = clock();
	
	double optcost = moloptAnnealingBFGS(n,nIters,tempfunc,&mol);

	end = clock();	
	
	std::ofstream fs;
	fs.open("alphas.txt",std::fstream::out);
	
	for(int i=0;i<n;i++)
	{
		fs << mol->coords[i].x << "," << mol->coords[i].y << "\n";
	}
	fs.close();
	
	double err = fabs((optimal[n-2]-optcost)/optimal[n-2])*100;
	double time = (double)(end - begin)/CLOCKS_PER_SEC;

	std::string datafile = "data/";
	datafile.append(argv[1]);
	datafile.append(".csv");
	
	printf("%s\n",datafile.c_str());
	
	fs.open(datafile,std::fstream::app);
	
	if(fs)
	{
		fs << n << ",";
		
		for (int i=0;i<n;i++)
		{
			
			fs << mol->alphas[i];
			
			if(i<n-1)
			{
				fs << ",";
			}
			else
			{
				fs << "," << optcost << "\n";
			}
		}
	}
	
	printf("Atoms\tOptimal\tFound\tError\tGen\tTime\n");
	printf("%i\t%2.2f\t%2.2f\t%2.2f\t%i\t%2.2f\n",n,optimal[n-2],optcost,err,nIters,time);
	return 0;
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
