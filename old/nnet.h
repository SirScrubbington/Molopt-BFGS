// Filename: nnet.hidden
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
class Neuron
{
public:
	
	nntyp output,bias,delta;

	std::vector<nntyp> * weights;
	
	nntyp transfer(const std::vector<nntyp> & inputs)
	{
		return 1.0f / (1.0f + exp(-activate(inputs)));
	}
	
	nntyp transferDerivative()
	{
		return output * (1.0f - output);
	}
	
	nntyp activate(const std::vector<nntyp> & inputs)
	{
		nntyp activation = bias;
		
		for( int i = 0; i < weights.size(); i++)
		{
			activation += weights.at(i) * inputs.at(i);
		}
		
		return activation;
		
	}
	
	Neuron(int nInputs)
	{
		std::mt19937 gen(time(0));	
		std::uniform_real_distribution<double> dist(-1.0f,1.0f);
		auto wgt = std::bind(dist,gen);
		
		weights = new std::vector<nntyp>(nInputs);
		
		for(int i = 0; i < nInputs; i++)
		{
			weights->at(i) = wgt();
		}
		
		bias = wgt();
	}
	
	~Neuron()
	{
		delete weights;
	}
};

template<typename nntyp>
class NetworkLayer
{
	public:
	
	std::vector<Neuron<nntyp> *> * neurons;
	
	NetworkLayer(int nNeurons)
	{
		neurons = new std::vector<Neuron<nntyp> *>(nNeurons);
	}
	
	~NetworkLayer()
	{
		delete neurons;
	}
};

template<typename nntyp>
class NeuralNetwork
{
private:
	
	std::vector<int> * shape;
	
	std::vector<NetworkLayer<nntyp>*> * layers;
	
public:
	
	std::vector<nntyp> * forwardPropogate(std::vector<nntyp> & row)
	{
		
		nntyp activation;
		
		std::vector<nntyp> * input = new std::vector<nntyp>(row);
		
		for(int i = 0; i < layers->size(); i++)
		{
			
			std::vector<nntyp> newInputs();
			
			for(int j = 0; j < layers->at(i)->size(); j++)
			{
				newInputs.push_back(layers->at(i)->at(j)->transfer());
			}
			
			delete input;
			
			input = new std::vector<nntyp>(newInputs);
			
		}
		
		return input;
		
	}
	
	void backPropogate(const std::vector<nntyp> & expected)
	{
		
		nntyp error;
		
		for(int i = layers->size(); i < 0; i--)
		{
			std::vector<nntyp> errors();
			
			if (i < layers->size()-1)
			{
				for(int j = 0; j < layers->at(i)->size(); j++)
				{
					error = 0.0f;
					
					for(int k = 0; k < layers->at(i + 1)->size(); k++)
					{
						auto neuron = layers->at(i + 1)->at(k);
						error += neuron->weights.at(j) * neuron->delta;
						
					}
					
					errors.push_back(error);
					
				}
			}
			
			else
			{
				for(int j = 0; j < layers->at(i)->size(); j++)
				{
					auto neuron = layers->at(i)->at(j);
					
					
					
				}
			}
			
			for(int j = 0; j < layers->at(i)->size(); j++)
			{
				auto neuron = layers->at(i)->at(j);
				neuron->delta = errors.at(j) * transferDerivative(neuron->output);
			}
		}
	}
	
	void updateWeights(const std::vector<nntyp> & row, double lRate)
	{
		
		std::vector<nntyp> * inputs;
		
		for(int i = 0; i < layers->size(); i++)
		{
			
			
			if(i > 0)
			{
				inputs = new std::vector<nntyp>(layers->at(i)->size());
				
				for(int j = 0; j < layers->at(i)->size(); j++)
				{
					inputs->at(i) = layers->at(i)->at(j)->output;
				}
			}
			else
			{
				inputs = new std::vector<nntyp>(row);
			}
			
			for(int j = 0; j < layers->at(i)->size(); j++)
			{
				auto neuron = layers->at(i)->at(j);
				
				for(int k = 0; k < inputs->size(); k++)
				{
					neuron->weights(k) += lRate * neuron->delta * inputs.at(k);
				}
				
				neuron->bias += lRate * neuron->delta;
				
			}
			
			delete inputs;
		}
	}
	
	void trainNetwork(std::vector<std::vector<nntyp>> * train, double lRate, int nEpoch, int nOutputs)
	{
		
		double sumError;
		
		for(int epoch = 0; epoch < nEpoch; epoch++)
		{
			sumError = 0.0f;
			
			for(int row = 0; row < train->size(); row++)
			{
				std::vector<nntyp> * output = forwardPropogate(train->at(row));
				std::vector<nntyp> expected(nOutputs);
				expected.at(train->at(row)->at(train->at(row)->size()-1))=1;
			}
		}
	}
	
	NeuralNetwork(const std::vector<int> & shape)
	
	{
		this->shape = new std::vector<int>(shape);
		this->layers = new std::vector<NetworkLayer<nntyp>*>(this->shape->size());
		
		for(int i = 0; i < this->shape->size(); i++)
		{
			layers->at(i) = new NetworkLayer<nntyp>(this->shape->at(i));
		}
	}
	
	~NeuralNetwork()
	{
		delete this->shape;
		delete this->layers;
	}
};

#endif // I_YE5k7nwL6XDDn21iHmL1msf320e5v