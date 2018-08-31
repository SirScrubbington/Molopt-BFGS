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
