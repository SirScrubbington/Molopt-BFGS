void crossFunc(void * a, void * b, int n)
{
	
	lbfgsfloatval_t temp;
	
	Molecule<lbfgsfloatval_t> * molA = (Molecule<lbfgsfloatval_t> *)a;
	Molecule<lbfgsfloatval_t> * molB = (Molecule<lbfgsfloatval_t> *)b;
	
	if (molA->n < n)
		return;
	
	for(int i = 0; i < n; i++)
	{
		
		//printf("%i: A: %f, B: %f ->",i,molA->alphas[i],molB->alphas[i]);
		
		temp = molA->alphas[i];
		molA->alphas[i] = molB->alphas[i];
		molB->alphas[i] = temp;
		
		//printf("A: %f, B: %f\n",i,molA->alphas[i],molB->alphas[i]);
		
	}
}

/*
void mutFunc(void * arg, int n)
{
	Molecule<lbfgsfloatval_t> * mol = (Molecule<lbfgsfloatval_t> *)arg;
	
	for(int i=0;i<n;i++)
	{
		mol->alphas[rng()] = 0;
	}
}
*/

int sortFunc(void * a, void * b)
{
	Molecule<lbfgsfloatval_t> * molA = (Molecule<lbfgsfloatval_t> *)a;
	Molecule<lbfgsfloatval_t> * molB = (Molecule<lbfgsfloatval_t> *)b;
	
	lbfgsfloatval_t ca = molA->updateSystem();
	lbfgsfloatval_t cb = molB->updateSystem();
	
	return ca < cb;
}

lbfgsfloatval_t moloptGeneticBFGS
(
	int n,
	int popSiz,
	int swapSiz,
	int generations,
	int keep,
	double mutRate,
	//mutationFunc mutFunc,
	crossoverFunc crossFunc,
	sortingFunc sortFunc,
	Molecule<lbfgsfloatval_t> ** save
)
{
	std::mt19937 gen(time(0));	
	std::uniform_int_distribution<int> dist(-180,180);
	std::uniform_real_distribution<double> acpt(0.0f,1.0f);
	
	auto randAngle = std::bind(dist,gen);
	auto accept = std::bind(acpt,gen);
	
	bool hasCrossovers;
	
	if (keep % 2 != 0)
		keep+=1;
	
	std::vector<Molecule<lbfgsfloatval_t>*> population(popSiz);
	
	for(int i = 0 ; i < popSiz ; i++)
	{
		population[i] = new Molecule<lbfgsfloatval_t>(n);
	}

	for(int j = keep; j < popSiz; j++)
		{
			hasCrossovers = true;

			while(hasCrossovers)
			{
				for(int k = 0; k < n; k++)
				{
					population[j]->alphas[k] = randAngle();
				}
				if(population[j]->updateSystem() < 0)
				{
					hasCrossovers = false;
				}
			}
		}

	for(int i = 0; i < generations; i++)
	{
		std::sort(population.begin(),population.end(),sortFunc);
		
		/*
		if(i == generations-1)
		{

		}
		*/
		//break;
		
		for(int j = 0; j < (keep / 2); j++)
		{
			crossFunc(population[j],population[keep-j-1],swapSiz); 
		}
		
		for(int j = keep; j < popSiz; j++)
		{
			hasCrossovers = true;

			if(accept() < mutRate * ((double)generations-(double)i)/(double)generations)
			{
				moloptAnnealingBFGS(n,10,tempFunc,&population[j]);
			}
			else
			{	
				while(hasCrossovers)
				{
					for(int k = 0; k < n; k++)
					{
						population[j]->alphas[k] = randAngle();
					}
					if(population[j]->updateSystem() < 0)
					{
						hasCrossovers = false;
					}
				}
			}
		}
	}
	
	std::sort(population.begin(),population.end(),sortFunc);
	
	for(int j = 0; j < popSiz; j++)
	{
		printf("%f\n",population[j]->updateSystem());
	}
	
	for(int i=1;i<popSiz;i++)
		delete population[i];
	
	*save = population[0];
	
	return population[0]->updateSystem();
}