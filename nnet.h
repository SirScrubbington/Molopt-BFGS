// Filename: nnet.h
// Author: Scrubbs
// Date: 2018-8-20
// Description: nnet.h header file

#ifndef I_YE5k7nwL6XDDn21iHmL1msf320e5v
#define I_YE5k7nwL6XDDn21iHmL1msf320e5v

#include <functional>
#include <random>
#include <vector>

#include "Matrix.h"

template<typename nntyp>
class NeuralNetwork
{
private:
	
	int nX, nY, nH;
	
	std::vector<nntyp> * x;
	std::vector<nntyp> * h;
	std::vector<nntyp> * y;
	
	std::vector<nntyp> * xDelta;
	std::vector<nntyp> * hDelta;
	std::vector<nntyp> * yDelta;
	
	Matrix<nntyp> * w1;
	Matrix<nntyp> * w2;
		
public:

	double transfer(double activation)
	{
		return 1.0f / (1.0f + exp(-activation));
	}

	void forwardPropogate()
	{
		double activation, output;
		
		// input -> hidden
		
		std::vector<nntyp> newH();
		
		for(int i = 0; i < nH; i++)
		{
			
			activation = h[i];
			
			// hidden layer activation
			for(int j = 0; j < nX; j++)
			{
				activation += w1->get(i,j) * x.at(j);
			}
			
			// transfer function
			output = transfer(activation);
			
			//newH.push_back(output);
			h.at(i) = output;
		}
		
		// hidden -> output
		
		std::vector<nntyp> newY();
		
		for(int i = 0; i < nY; i++)
		{
			activation = y[i];
			
			// output layer activation
			for(int j = 0; j < nH; j++)
			{
				activation += w2->get(i,j) * h.at(j);
			}
			
			// transfer function
			output = transfer(activation);
			
			//newY.push_back(output);
			y.at(i) = output;
		}
	}
	
	double transferDerivative(double output)
	{
		return output * (1.0f - output);
	}
	
	void backpropY(std::vector<nntyp> * expected)
	{
		std::vector<nntyp> errors();
		
		// calculate errors
		for(int i = 0; i < nY; i++)
		{
			errors.push_back(expected.at(i) - y.at(i));
		}
		
		// calculate deltas
		for(int i = 0; i < nY; i++)
		{
			yDelta.at(i) = errors.at(i) * transferDerivative(y.at(i));
		}
	}
	
	void backpropH(std::vector<nntyp> * expected)
	{
		std::vector<nntyp> errors();
		
		double error;
		
		// calculate errors
		for(int i = 0; i < nH; i++)
		{
			error = 0.0f;
		
			for( int j = 0; j < nY; j++)
			{
				error += w2->get(i,j) * yDelta.at(j);
			}
			
			errors.append(error);
		}
		
		// calculate deltas
		for(int i = 0; i < nH; i++)
		{
			hDelta.at(i) = errors.at(i) * transferDerivative(h.at(i));
		}
	}
	
	void backpropX(std::vector<nntyp> * expected)
	{
		std::vector<nntyp> errors();
		
		double error;
		
		// calculate errors
		for(int i = 0; i < nX; i++)
		{
			for(int j = 0; j < nH; j++)
			{
				error += w1->get(i,j) * hDelta.at(j);
			}
			
			errors.append(error);
			
		}
		
		// calculate deltas
		for(int i = 0; i < nX; i++)
		{
			xDelta.at(i) = errors.at(i) * transferDerivative(x.at(i));
		}
	}
	
	void backPropogate(std::vector<nntyp> * expected)
	{
		
	}
	
	
	NeuralNetwork(int nX, int nY, int nH)
	{
		
		std::mt19937 gen(time(0));	
		std::uniform_real_distribution<double> weights(-1.0f,1.0f);
		
		auto wgt = std::bind(weights,gen);
		
		this->nX = nX;
		this->nY = nY;
		this->nH = nH;
		
		this->x = new std::vector<nntyp>(nX);
		this->h = new std::vector<nntyp>(nH);
		this->y = new std::vector<nntyp>(nY);
		
		this->xDelta = new std::vector<nntyp>(nX);
		this->hDelta = new std::vector<nntyp>(nH);
		this->yDelta = new std::vector<nntyp>(nY);
		
		this->w1 = new Matrix<nntyp>(nX,nH);
		this->w2 = new Matrix<nntyp>(nH,nY);
		
		for(int i = 0; i < nX; i++)
		{
			for(int j = 0; j < nH; j++)
			{
				w1->set(i,j,wgt());
			}
		}
		
		for(int i = 0; i < nH; i++)
		{
			for(int j = 0; j < nY; j++)
			{
				w2->set(i,j,wgt());
			}
		}
	}
	
	~NeuralNetwork()
	{
		delete x;
		delete y;
		delete h;
		
		delete w1;
		delete w2;
	}
};

#endif // I_YE5k7nwL6XDDn21iHmL1msf320e5v
