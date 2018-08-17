// Filename: Molecule.tpp
// Author: Scrubbs
// Date: 2018-8-17
// Description: Molecule implementation file

#ifndef I_B3G3Q56W4ssY2nQd25TE33mP5DPu0
#define I_B3G3Q56W4ssY2nQd25TE33mP5DPu0

moltyp * Molecule<moltyp>::getRealAngleList()
{
	real[0] = alphas[0];
	
	for(int i=1;i<n;i++)
	{
		real[i] = real[i-1] + alphas[i];
		if(real[i] >= DEG2RAD*360)
			real[i] -= DEG2RAD*360;
	}
	return real;
}

Point Molecule<moltyp>::getCoords()
{
	return {distance*cos(angle),distance*sin(angle)};
}

Point * Molecule<moltyp>::getSystemCoordinates()
{
	for(int i=2;i<n;i++)
	{
		coords[i] = coords[i-1] + getCoords(1.0f,real[i-1]);
	}	
	
	return coords;
}

moltyp Molecule<moltyp>::getR(moltyp x, moltyp y)
{
	return sqrt(pow(x,2) + pow(y,2));
}

Matrix<moltyp> * Molecule<moltyp>::getEulerianDistanceMatrix()
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
				Point d = (coords[i] - coords[j]).absolute();
				distance->set(i,j,getR(d.x,d.y));
			}
		}
	}
	return distance;
}

moltyp Molecule<moltyp>::getSystemEnergy()
{
	moltyp V = 0;
	
	for(int j=0;j<n;j++)
		for(int i=0;i<j;i++)
		{
			V = V + ((1.0f/pow(distance->get(i,j),12.0f))-(2.0f/pow(distance->get(i,j),6.0f)));
		}
	
	return V;
	
}
	
Molecule<moltyp>::Molecule(int n)
{
	distance = new Matrix<moltyp>(n,n);
	
	alphas = (moltyp*)malloc(sizeof(*alphas)*n);
	real = (moltyp*)malloc(sizeof(*real)*n);
	gradients = (moltyp*)malloc(sizeof(*gradients)*n); 
	
	coords = (Point*)malloc(sizeof(*coords)*n);
	
	for(int i=0;i<n;i++)
	{
		alphas[i]=0.0f;
		gradients[i]=0.0f;
		coords[i]={0.0f,0.0f};
	}
	
	coords[1]={1.0f,0.0f};
	
	real = getRealAngleList();

	coords = getSystemCoordinates();
	
	distance = getEulerianDistanceMatrix();
}
	
Molecule<moltyp>::~Molecule()
{
	free(distance);
	free(alphas);
	free(real);
	free(gradients);
	free(coords);
}

#endif // I_B3G3Q56W4ssY2nQd25TE33mP5DPu0
