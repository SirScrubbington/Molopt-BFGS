// Filename: nnet.h
// Author: Scrubbs
// Date: 2018-8-20
// Description: nnet.h header file

#ifndef I_YE5k7nwL6XDDn21iHmL1msf320e5v
#define I_YE5k7nwL6XDDn21iHmL1msf320e5v

#include <cmath>

template<typename nntyp>
class NeuralNetwork
{
private:
	nntyp * x;
	nntyp * w1;
	nntyp * w2;
	nntyp * y;
	nntyp * out;
	
	nntyp layer1;
	
	
	int n;
public:
	nntyp dot(nntyp * a, nntyp * b, int n)
	{
		nntyp dp = 0.0f;
		for(int i=0;i<n;i++)
		{
			dp += a[i] * b[i];
		}
		return dp;
	}
	
	nntyp sigmoid(nntyp t)
	{
		return t / (1.0f + fabs(x));
	}
	
	nntyp sigmoidD(nntyp t)
	{
		return sigmoid(t) * (1.0f - sigmoid(t));
	}
	
	NeuralNetwork(nntyp * x, nntyp * y, )
	{
		x = malloc(sizeof(nntyp)*n);
		w1 = malloc(sizeof(nntyp)*n);
		w2 = malloc(sizeof(nntyp)*n);
		y = malloc(sizeof(nntyp)*n);
		out = malloc(sizeof(nntyp)*n);
		int n;
	}
	
	void feedForward()
	{
		
	}
	
	void backPropogation()
	{
		
	}
};

#endif // I_YE5k7nwL6XDDn21iHmL1msf320e5v
