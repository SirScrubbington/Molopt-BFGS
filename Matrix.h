// Filename: Matrix.h
// Author: Scrubbs
// Date: 2018-8-16
// Description: Matrix.h header file

#ifndef I_HJE28Mf1pDDIHkjF0Dc418085tCF2
#define I_HJE28Mf1pDDIHkjF0Dc418085tCF2

#include <cstdio>
#include <stdexcept>

template <class T>
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
	
	inline void set(int i, int j, double v)
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
