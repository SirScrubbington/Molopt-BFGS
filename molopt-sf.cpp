// Filename: molopt-sf.cpp
// Author: Scrubbs
// Date: 20nTests-8-7
// Description: This file was generated as a combination of all non-library files associated with the original molopt.cpp
// As a result, all of the original file headers have been kept in place so it is clear which section of the file belongs
// to which program. This file requires an /include directory, containing the header-only implementations of Eigen and LBFGSpp.
// Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
// LBFGSPP: https://github.com/yixuan/LBFGSpp


// Filename: molopt.cpp
// Author: Scrubbs
// Date: 20nTests-8-7
// Description: molopt implementation file

#ifndef I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7
#define I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7

// Filename: molopt.h
// Author: Scrubbs
// Date: 2018-8-7
// Description: molopt.h header file

#ifndef I_NfAv5xV18d4u1Tj758VJpX4i18fsL
#define I_NfAv5xV18d4u1Tj758VJpX4i18fsL

#include <string>
#include <random>
#include <fstream>
#include <iostream>
#include <functional>

#include <ctime>
#include <cfloat>
#include <cstdio>
#include <cstdlib>

#define _USE_MATH_DEFINES

#include <cmath>

#include "include/Eigen/Core"
#include "include/LBFGS.h"

// Filename: Point.h
// Author: Scrubbs
// Date: 2018-8-16
// Description: Point.h header file

#ifndef I_A1jUxXF381b5882Pp0Q1r2ps81P6v
#define I_A1jUxXF381b5882Pp0Q1r2ps81P6v

#include <cmath>

struct Point 
{
	double x,y;
	inline Point operator=(Point a)
	{
		x=a.x;
		y=a.y;
		return a;
	}
	inline Point operator+(Point a)
	{
		return {a.x+x,a.y+y};
	}
	inline Point operator-(Point a)
	{
		return {x-a.x,y-a.y};
	}
	inline bool operator==(Point a)
	{
		return (a.x==x && a.y==y);
	}
	inline Point absolute()
	{
		return {fabs(x),fabs(y)};
	}
};

#endif // I_A1jUxXF381b5882Pp0Q1r2ps81P6v

// Filename: Matrix.h
// Author: Scrubbs
// Date: 2018-8-16
// Description: Matrix.h header file

#ifndef I_HJE28Mf1pDDIHkjF0Dc418085tCF2
#define I_HJE28Mf1pDDIHkjF0Dc418085tCF2

#include <cstdio>
#include <stdexcept>

template <typename T>
class Matrix
{
public:
	int x,y;
	T * data;

	Matrix(int x, int y)
	{
		this->x = x;
		this->y = y;
		data = (T*)malloc(sizeof(*data)*x*y);
	}
	
	~Matrix()
	{
		free(data);
	}
	
	inline void set(int i, int j, T v)
	{
		if((i*x+j) > (x*y))
			throw std::invalid_argument("Attempting to access out of bounds index");
		
		data[i*x+j]=v;
	}
	
	inline T get(int i, int j)
	{
		if((i*x+j) > (x*y))
			throw std::invalid_argument("Attempting to access out of bounds index");
		
		return data[i*x+j];
	}
};

#endif // I_HJE28Mf1pDDIHkjF0Dc418085tCF2

// Filename: Molecule.h
// Author: Scrubbs
// Date: 2018-8-17
// Description: Molecule.h header file

#ifndef I_v46cQspgT8KK6epE1iB7GAk686t4t
#define I_v46cQspgT8KK6epE1iB7GAk686t4t


#define NO_CROSSOVERS 0
#define ALLOW_CROSSOVERS 1

#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>

#ifndef D2R
	#define D2R
	const double DEG2RAD = M_PI/180.0f;
#endif // DEG2RAD

#ifndef R2D
	#define R2D
	const double RAD2DEG = 180.0f/M_PI;
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
		
	VectorXd alphas;
	VectorXd real;
	
	std::vector<Point> coords;
	
	const double R2 = DEG2RAD*360;
	
	int n;
	
	void getRealAngleList()
	{
		real[0] = alphas[0];
		
		for(int i=1;i<n;i++)
		{
			real[i] = real[i-1] + alphas[i];
			if(real[i] >= R2)
				real[i] -= R2;
		}
	}

	Point getCoords(double distance, double angle)
	{
		return {distance*cos(angle),distance*sin(angle)};
	}

	void getSystemCoordinates()
	{
		for(int i=2;i<n;i++)
		{
			//printf("%i\n",i);
			coords.at(i) = coords.at(i-1) + getCoords(1.0f,real[i-1]);
		}	
	}

	double getR(double x, double y)
	{
		return sqrt(x*x + y*y);
	}

	void getEulerianDistanceMatrix()
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
					Point d = (coords.at(i) - coords.at(j)).absolute();
					distance->set(i,j,getR(d.x,d.y));
				}
			}
		}
	}

	double getSystemEnergy()
	{
		double V = 0;
		
		for(int j=0;j<n;j++)
			for(int i=0;i<j;i++)
			{
				V = V + ((1.0f/std::pow(distance->get(i,j),12.0f))-(2.0f/std::pow(distance->get(i,j),6.0f)));
			}
		
		return V;
		
	}
		
	Molecule(int n)
	{
		this->n = n;
		
		distance = new Matrix<double>(n,n);
		
		alphas = VectorXd::Zero(n);
		real = VectorXd::Zero(n);
		
		coords = std::vector<Point>(n);
		
		coords.at(0) = {0.0f,0.0f};
		coords.at(1) = {1.0f,0.0f};
		
		updateSystem();

	}
	
	double updateSystem()
	{
		getRealAngleList();
		getSystemCoordinates();
		getEulerianDistanceMatrix();
		return getSystemEnergy();
	}
	
	// calculate gradients for BFGS
	double operator()(const VectorXd& x, VectorXd& grad)
	{
		double a, b, c, d;

		alphas = x;
		
		auto cost = updateSystem();
		
		for(int m = 0; m < n; m++)
		{
			grad[m] = 0.0f;
			
			for(int i = 0; i < m - 1; i++)
			{
				for(int j = m + 1; j < n; j++)
				{
					a = 1.0f / std::pow(distance->get(i,j),14.0f);
					b = 1.0f / std::pow(distance->get(i,j),8.0f);
					
					c = (coords[i].x - coords[j].x) * (coords[m].y - coords[i].y);
					d = (coords[i].y - coords[j].y) * (coords[j].x - coords[m].x);
					
					grad[m] += (a - b) * (c + d);
				}
			}
			grad[m] *= -12.0f;
		}	
		
		return cost;
	}
};

#endif // I_v46cQspgT8KK6epE1iB7GAk686t4t


#include <string>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef D2R
	#define D2R
	const double DEG2RAD = M_PI/180.0f;
#endif // DEG2RAD

#ifndef R2D
	#define R2D
	const double RAD2DEG = 180.0f/M_PI;
#endif // RAD2DEG

void makeDirectory(const char * path)
{
#if defined(_WIN32)
	_mkdir(path); // can be used on Windows
#else 
	mkdir(path,0733); // can be used on non-Windows
#endif
}

#endif // I_NfAv5xV18d4u1Tj758VJpX4i18fsL


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

double bruteForceMolecule
(
	int n,
	int iters,
	Molecule ** molecule
)
{
	double bestAngle, curTemp, bestTemp;
	
	if (*molecule == NULL)
	{
		*molecule = new Molecule(n);
	}
	
	Molecule * mol = *molecule;
	
	bestTemp = mol->updateSystem();
	
	for(int i = 0; i < n; i++)
	{
		bestAngle = mol->alphas[i];
		
		for(int j = -180; j <= 180; j += 360/iters)
		{
			mol->alphas[i] = j;
			
			curTemp = mol->updateSystem();
			
			if(curTemp < bestTemp)
			{
				bestAngle = j;
				bestTemp = curTemp;
			}
			else
			{
				mol->alphas[i] = bestAngle;
			}
		}
	}
	return mol->updateSystem();
}

double IterativeAnnealing
(
	int n, 
	int iters,
	int passes,
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
	
	for(int i = 1; i < iters; i++)
	{
		T = temp((double)iters/(double)i);
		
		for(int j = 0; j < n; j++)
		{
			pE = mol->updateSystem();
				
			double tv = mol->alphas[j];
				
			mol->alphas[j] = randAngle();
				
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
					mol->alphas[j] = tv;
				}
			}
		}
	}

	return mol->updateSystem();
}

double getError(double target, double found)
{
	return (((target - found)/target)*100.0f);
}

double getTime(clock_t begin, clock_t end)
{
	return (double)(end-begin)/(CLOCKS_PER_SEC);
}

#define WRITE_ALPHAS 1
#define WRITE_COORDS 2
#define NO_PRINT 4

int main(int argc, char ** argv)
{
	int n = 0;
	
	short properties = 0, scorefile = 0;
	
	if (argc > 1)
	{
		n = atoi(argv[1]);
		
		for(int i = 2; i < argc; i++)
		{
			if(strcmp(argv[i],"-a") == 0) 
			{
				properties |= WRITE_ALPHAS;
			}
			else if(strcmp(argv[i],"-c") == 0)
			{
				properties |= WRITE_COORDS;
			}
			else if(strcmp(argv[i],"-q") == 0)
			{
				properties |= NO_PRINT;
			}
			else if(strcmp(argv[i],"-f") == 0 && argc > i+1)
			{
				scorefile = i + 1;
			}
		}
	}
	else
	{
		if(!(properties & NO_PRINT)) printf("Usage: ./molopt[.exe/.out] natoms [-a] [-c] [-np] [-f filename]");
	}
	
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	// 							   Algorithm Setup					  	 	  //
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////
	// 							Constant Testing Data						  //
	////////////////////////////////////////////////////////////////////////////
	
	const int nTests = 18; // the number of tests which will be performed
	
	const int nIterations = 50; // Number of annealing passes
	
	gen = new std::mt19937(time(0)); // Seed random number generator
	
	std::fstream fs; // Filestream for writing data
	
	clock_t begin,end; // clock structures for timing the execution times of each algorithm
	
	Molecule circmol(n), squaremol(n);
	
	std::vector<double> optimal;
	
	// if there is an existing file with best costs, read from it
	
	std::ifstream inFile;
	
	if (scorefile > 0)
	{
		inFile.open(argv[scorefile]);
	}
	else
	{
		inFile.open("found.csv");
	}
	
	if (inFile)
	{
		double outvar;
		
		while(inFile >> outvar)
		{
			optimal.push_back(outvar);
		}
	}
	
	// for any undefined best costs, assign that index of the best costs array to zero
	
	for(int i = optimal.size(); i < n-1; i++)
	{
		optimal.push_back(0.0f);
	}
	
	inFile.close();
	
	// configure the default angles for the square template
	double circlangle = DEG2RAD * 360.0f / n;
	
	for(int i = 0; i < n; i++)
	{
		circmol.alphas[i] = circlangle;
	}

	// configure the default angles for the circular template
	double squarangle = DEG2RAD * 90.0f;
	
	if(n / 3 > 0)
	{
		for(int i = (n/3); i < n; i += (n / 3))
		{
			if (i >= n) break;
			squaremol.alphas[i] = squarangle;
		}
	}
	
	Molecule * molecules[nTests];
	
	double times[nTests] = {0.0f};
	double costs[nTests] = {0.0f};
	int iters[nTests] = {0};
	
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	// 							Algorithm Application				  	 	  //
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	// loop over all atoms 'niters' times at once, then progress to next atom //
	////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////// linear configuration ///////////////////////////
	
	begin = clock();
	
	molecules[0] = new Molecule(n);

	costs[0] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[0]);
	
	end = clock();

	times[0] = getTime(begin,end);
	
	iters[0] = nIterations;

	////////////////////////// circular configuration ///////////////////////////
	
	begin = clock();
	
	molecules[1] = new Molecule(circmol);
	
	costs[1] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[1]);
	
	end = clock();
	
	times[1] = getTime(begin,end);
	
	iters[1] = nIterations;
	////////////////////////// square configuration ////////////////////////////
	
	begin = clock();

	molecules[2] = new Molecule(squaremol);
	
	costs[2] = BeginOrientatedAnnealing(n,nIterations,tempfunc,&molecules[2]);
	
	end = clock();
	
	times[2] = getTime(begin,end);
	
	iters[2] = nIterations;
	
	////////////////////////////////////////////////////////////////////////////
	// 				  annealing search atom optimisation method		  		  //
	////////////////////////////////////////////////////////////////////////////
	
	const int passes = 10; // number of repeats for repeated search algorithm
	
	/////////////////////////// linear configuration ///////////////////////////
	
	begin = clock();
	
	molecules[3] = new Molecule(n);
	
	costs[3] = IterativeAnnealing(n,nIterations,passes,tempfunc,&molecules[3]);
	
	end = clock();
	
	times[3] = getTime(begin,end);
	
	iters[3] = nIterations;

	////////////////////////// circular configuration ///////////////////////////
	
	begin = clock();
	
	molecules[4] = new Molecule(circmol);
	
	costs[4] = IterativeAnnealing(n,nIterations,passes,tempfunc,&molecules[4]);
	
	end = clock();
	
	times[4] = getTime(begin,end);
	
	iters[4] = nIterations;
		
	
	////////////////////////// square configuration ////////////////////////////
		
	begin = clock();
	
	molecules[5] = new Molecule(squaremol);
	
	costs[5] = IterativeAnnealing(n,nIterations,passes,tempfunc,&molecules[5]);
	
	end = clock();
	
	times[5] = getTime(begin,end);
	
	iters[5] = nIterations;
	
	
	////////////////////////////////////////////////////////////////////////////
	// 					brute force atom optimisation method		  	 	  //
	////////////////////////////////////////////////////////////////////////////
	
	const int angleIter = 50;
	
	/////////////////////////// linear configuration ///////////////////////////
	
	begin = clock();
	
	molecules[6] = new Molecule(n);
	
	costs[6] = bruteForceMolecule(n,angleIter,&molecules[6]);
	
	end = clock();
	
	times[6] = getTime(begin,end);
	
	iters[6] = angleIter;
	

	////////////////////////// circular configuration ///////////////////////////
	
	begin = clock();
	
	molecules[7] = new Molecule(circmol);
	
	costs[7] = bruteForceMolecule(n,angleIter,&molecules[7]);
	
	end = clock();
	
	times[7] = getTime(begin,end);
	
	iters[7] = angleIter;
	
	////////////////////////// square configuration ////////////////////////////
		
	begin = clock();
	
	molecules[8] = new Molecule(squaremol);
	
	costs[8] = bruteForceMolecule(n,angleIter,&molecules[8]);
	
	end = clock();
	
	times[8] = getTime(begin,end);
	
	iters[8] = angleIter;
	
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	// 								BFGS Optimizer					  	 	  //
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	
	LBFGSParam<double> param; // Default parameters for the LBFGS routine
	param.max_iterations = 50; // Maximum number of LBFGS iterations
	LBFGSSolver<double> solver(param); // LBFGS++ Solver object for applying the BFGS algorithm
	
	////////////////////////////////////////////////////////////////////////////
	// loop over all atoms 'niters' times at once, then progress to next atom //
	////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////// linear configuration ///////////////////////////
	
	begin = clock();
	
	molecules[9] = new Molecule(*molecules[0]);

	solver.minimize((*molecules[9]),molecules[9]->alphas,costs[9]);
	
	end = clock();

	times[9] = getTime(begin,end);
	
	iters[9] = param.max_iterations;

	////////////////////////// circular configuration ///////////////////////////
	
	begin = clock();
	
	molecules[10] = new Molecule(*molecules[1]);

	solver.minimize((*molecules[10]),molecules[10]->alphas,costs[10]);
	
	end = clock();
	
	times[10] = getTime(begin,end);
	
	iters[10] = param.max_iterations;

	////////////////////////// square configuration ////////////////////////////
	
	begin = clock();

	molecules[11] = new Molecule(*molecules[2]);
	
	solver.minimize((*molecules[11]),molecules[11]->alphas,costs[11]);
	
	end = clock();
	
	times[11] = getTime(begin,end);
	
	iters[11] = param.max_iterations;
	
	
	////////////////////////////////////////////////////////////////////////////
	// 					brute force atom optimisation method		  	 	  //
	////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////// linear configuration ///////////////////////////
	
	begin = clock();
	
	molecules[12] = new Molecule(*molecules[3]);
	
	solver.minimize((*molecules[12]),molecules[12]->alphas,costs[12]);
	
	end = clock();
	
	times[12] = getTime(begin,end);
	
	iters[12] = param.max_iterations;
	
	////////////////////////// circular configuration ///////////////////////////
	
	begin = clock();
	
	molecules[13] = new Molecule(*molecules[4]);
	
	solver.minimize((*molecules[13]),molecules[13]->alphas,costs[13]);
	
	end = clock();
	
	times[13] = getTime(begin,end);
	
	iters[13] = param.max_iterations;
	
	
	////////////////////////// square configuration ////////////////////////////
		
	begin = clock();
	
	molecules[14] = new Molecule(*molecules[5]);
	
	solver.minimize((*molecules[14]),molecules[14]->alphas,costs[14]);
	
	end = clock();
	
	times[14] = getTime(begin,end);
	
	iters[14] = param.max_iterations;
	
	////////////////////////////////////////////////////////////////////////////
	// 					brute force atom optimisation method		  	 	  //
	////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////// linear configuration ///////////////////////////
	
	begin = clock();
	
	molecules[15] = new Molecule(*molecules[6]);
	
	solver.minimize((*molecules[15]),molecules[15]->alphas,costs[15]);
	
	end = clock();
	
	times[15] = getTime(begin,end);
	
	iters[15] = param.max_iterations;
	

	////////////////////////// circular configuration ///////////////////////////
	
	begin = clock();
	
	molecules[16] = new Molecule(*molecules[7]);
	
	solver.minimize((*molecules[16]),molecules[16]->alphas,costs[16]);
	
	end = clock();
	
	times[16] = getTime(begin,end);
	
	iters[16] = param.max_iterations;
	
	
	////////////////////////// square configuration ////////////////////////////
		
	begin = clock();
	
	molecules[17] = new Molecule(*molecules[8]);
	
	solver.minimize((*molecules[17]),molecules[17]->alphas,costs[17]);
	
	end = clock();
	
	times[17] = getTime(begin,end);
	
	iters[17] = param.max_iterations;
		
	//if(!(properties & NO_PRINT)) printf("Time 17: %f\n",times[17]);
	
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	// 								Output Writing					  	 	  //
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////
	// 								Filepath Setup 							  //
	////////////////////////////////////////////////////////////////////////////
	
	// Setup data write directories
	
	for(int j = 0; j < nTests; j++)
	{
		if (costs[j] < optimal[n-2])
		{
			optimal[n-2] = costs[j];
		}
	}
	
	if (scorefile > 0)
	{
		fs.open(argv[scorefile],std::fstream::out);
	}
	else
	{
		fs.open("found.csv",std::fstream::out);
	}
	
	for(int i = 0; i < optimal.size(); i++)
	{
		fs << optimal[i] << "\n";
	}
	
	fs.close();
	
	if(!(properties & NO_PRINT)) 
	{
		printf("Atoms\tOptimal\tFound\tError\tGen\tTime\n");// begin table
		for(int i = 0; i < 18; i++)
		{
			printf("%i\t%2.2f\t%2.2f\t%2.2f\t%i\t%2.2f\n",n,optimal[n-2],costs[i],getError(optimal[n-2],costs[i]),iters[i],times[i]);	
		}
	}
	
	if(properties & WRITE_COORDS)
	{
		makeDirectory("coords");
		std::string cFolder = "coords/",cSubfolder;
		cFolder = cFolder + std::to_string(n);
		makeDirectory(cFolder.c_str());
		cFolder = cFolder + "/";
		
		for(int j = 0; j < nTests; j++)
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
					fs << "," << iters[j] << "," << optimal[n-2] << "," << costs[j] << "," << times[j] << "\n";
			}	
			fs.close();
		}
	}

	if(properties & WRITE_ALPHAS)
	{
		
		makeDirectory("alphas");
		std::string dFolder = "alphas/",dSubfolder;	
		dFolder = dFolder + std::to_string(n);
		makeDirectory(dFolder.c_str());
		dFolder = dFolder + "/";
		
		for(int j = 0; j < nTests; j++)
		{
			Molecule * mol = molecules[j];
			
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
						fs << "," << iters[j] << "," << costs[j] << "," << times[j] << "\n";
					}
				}
			}
			
			fs.close();
		}
	}
}

#endif // I_MY3sy2M43uDRW3krlK4Oy3Et3o1B7

