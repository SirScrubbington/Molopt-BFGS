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

#include "Point.h"
#include "Matrix.h"
#include "Molecule.h"

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
