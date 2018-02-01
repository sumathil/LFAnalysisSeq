///// Basic Header File /////

#pragma once

#include <stdio.h>
#include <iostream>
#include <complex.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>

using namespace std;

#define num_buses 57
#define linefile 80

#define pi 3.14159265358979323846

typedef struct
{
	double real[linefile];
	double imag[linefile];
} rect;

typedef struct
{
	double real[num_buses][num_buses];
	double imag[num_buses][num_buses];
}rect_matrix;
