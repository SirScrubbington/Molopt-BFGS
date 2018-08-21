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

#include "bfgs/lbfgs.c"

#include "Point.h"
#include "Matrix.h"
#include "Molecule.h"

#ifndef DEG2RAD
	#define DEG2RAD (M_PI/180.0f)
#endif // DEG2RAD

#ifndef RAD2DEG
	#define RAD2DEG (180.0f/M_PI)
#endif // RAD2DEG

const lbfgsfloatval_t tolerance = pow(10,-6);

const lbfgsfloatval_t optimal[] = {
	-1.0f,
	-3.0f,
	-5.1f,
	-7.2f,
	-9.3f,
	-12.5f,
	-14.7f,
	-16.9f,
	-20.1f,
	-22.3f,
	-25.5f,
	-27.8f,
	-31.0f,
	-33.2f,
	-36.5f,
	-38.7f,
	-42.0f,
	-45.3f,
	-47.5f,
	-50.8f,
	-53.0f,
	-56.3f,
	-59.6f,
	-61.8f,
	-65.1f,
	-68.4f,
	-70.0f,
	-74.0f,
	-77.2f,
	-79.5f,
	-82.8f,
	-86.1f,
	-88.3f,
	-91.7f,
	-95.0f,
	-98.3f,
	-100.5f,
	-103.8f,
	-107.1f,
	-109.4f,
	-112.7f,
	-116.0f,
	-119.3f,
	-121.6f,
	-124.9f,
	-128.6f,
	-131.5f,
	-133.8f,
	-137.1f,
	-140.5f,
	-143.7f,
	-146.0f,
	-149.3f,
	-152.7f,
	-154.7f, // interpolated
	-157.4f, // interpolated
	-160.2f, // interpolated
	-163.0f, // interpolated
	-164.0f, // interpolated
	-166.1f
};

#endif // I_NfAv5xV18d4u1Tj758VJpX4i18fsL
