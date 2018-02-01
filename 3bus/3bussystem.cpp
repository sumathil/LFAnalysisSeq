#include <stdio.h>
#include <iostream>
#include <complex.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;
#define NUM_BUSES 3


typedef struct {
	float real[NUM_BUSES][NUM_BUSES];
	float imag[NUM_BUSES][NUM_BUSES];
} rect;


typedef struct {
	float mag[3][3];
	float ang[3][3];
} polar;

void mismatchpower(float *deltapower,polar *yp, float *busangle,float *busvolatge,float *apower);
void jacobian(float **jacobianmatrix, polar *yp, float *busangle, float *busvoltage);
void cg(float *output, float *deltapower, float **jacobianmatrix, float *tbestimated);
void matvec(float *outvec, float **inmat, float *invec);

int main()
{
	rect *y_rect = new rect;
	polar *y_polar = new polar;
	float *dpower = new float[NUM_BUSES];
	float *apower = new float[NUM_BUSES];
	float *busvolatage = new float[NUM_BUSES];
	float *busangle = new float[NUM_BUSES];
	float **jacobianm = new float *[NUM_BUSES];
	float *tbestimated_parameters = new float[NUM_BUSES];
	float *newestimates = new float[NUM_BUSES];
	float tol;
	float eps = 1e-4;
	for (unsigned int i = 0; i < NUM_BUSES; i++)
	{
		jacobianm[i] = new float[NUM_BUSES];
	}
	// Bus voltages
	busvolatage[0] = 1.05;
	busvolatage[1] = 1;
	busvolatage[2] = 1.04;

	//Actual power to the bus
	apower[0] = -4.0f;
	apower[1] = 2.0f;
	apower[2] = -2.5f;

	// Bus angles
	for (int i = 0; i < NUM_BUSES; i++)
	{
		busangle[i] = 0.0f;
		newestimates[i] = 0.0f;
	}
	tbestimated_parameters[0] = busangle[1];
	tbestimated_parameters[1] = busangle[2];
	tbestimated_parameters[2] = busvolatage[1];
	ifstream reactangle("rectangle.txt");
	for (int i = 0; i < NUM_BUSES; i++)
	{
		for (int j = 0; j < NUM_BUSES; j++)
		{
			if (reactangle.eof())
			{
				cout << "end of file reached" << endl;
				break;
			}
			reactangle >> y_rect->real[i][j] >> y_rect->imag[i][j];
		}
	}

	// Printing 
	/*for (int i = 0; i < NUM_BUSES; i++)
 	{
 			for (int j = 0; j < NUM_BUSES; j++)
 				{
 							cout << y_rect->real[i][j] << "+" << y_rect->imag[i][j] << "i" <<'\t';
 								}
 										cout << endl;
 											}
 												*/


	// Reading admittance matrix in r and theta
	ifstream polar("polar.txt");
	for (int i = 0; i < NUM_BUSES; i++)
	{
		for (int j = 0; j < NUM_BUSES; j++)
		{
			if (polar.eof())
			{
				cout << "end of file reached" << endl;
				break;
			}
			polar >> y_polar->mag[i][j] >> y_polar->ang[i][j];
		}
	}

	// printing the values
	/*for (int i = 0; i < NUM_BUSES; i++)
 	{
		for (int j = 0; j < NUM_BUSES; j++)
		{
			cout << y_polar->mag[i][j] << "/_" << y_polar->ang[i][j] <<'\t';
		}
		cout << endl;
	}*/

	for (int it = 0; it < 10; it++)
	{
	// Calculating the mismatch power
	mismatchpower(dpower, y_polar, busangle, busvolatage, apower);

		tol = 0.0f;
		for (int i = 0; i < NUM_BUSES; i++)
		{
			tol += dpower[i] * dpower[i];
		}
		tol = sqrt(tol);
		if (tol < eps)
		{
			cout << "number of iterations is " << it << endl;
			break;
		}
		jacobian(jacobianm, y_polar, busangle, busvolatage);

	
		cg(newestimates, dpower, jacobianm, tbestimated_parameters);
		for (int i = 0; i < NUM_BUSES; i++)
		{
			tbestimated_parameters[i] = tbestimated_parameters[i] + newestimates[i];
			cout << tbestimated_parameters[i] << endl;
		}
		busangle[1] = tbestimated_parameters[0];
		busangle[2] = tbestimated_parameters[1];
		busvolatage[1] = tbestimated_parameters[2];

	}

	delete[] dpower;
	delete[] apower;
	delete[] busvolatage;
	delete[] busangle;
	delete[] y_rect;
	delete[] y_polar;
	for (int i = 0; i < NUM_BUSES; i++)
	{

		delete[] jacobianm[i];
	}
	delete[] jacobianm;
	delete[] tbestimated_parameters;
    return 0;
}



void mismatchpower(float *deltapower,polar *yp, float *busangle, float *busvolatge,float *apower)
{
	float *powerflow = new float[NUM_BUSES];
	for (int j = 0; j < NUM_BUSES; j++)
	{
		powerflow[j] = 0.0f;
	}
	int k = 1;
	for (int j = 0; j < NUM_BUSES; j++)
	{
		powerflow[0] += yp->mag[k][j] * busvolatge[j] * busvolatge[1] * cos(yp->ang[k][j] + busangle[j] - busangle[1]);
		powerflow[1] += yp->mag[k+1][j] * busvolatge[j] * busvolatge[2] * cos(yp->ang[k+1][j] + busangle[j] - busangle[k+1]);
		powerflow[2] -= yp->mag[k][j] * busvolatge[j] * busvolatge[1] * sin(yp->ang[k][j] + busangle[j] - busangle[1]);
	}

	for (int i = 0; i < NUM_BUSES; i++)
	{
		deltapower[i] = apower[i] - powerflow[i];
	}
		delete[] powerflow;
}

void jacobian(float **jm, polar *yp, float *busangle, float *busvoltage)
{
	for (int i = 0; i < NUM_BUSES; i++)
	{
		for (int j = 0; j < NUM_BUSES; j++)
		{
			jm[i][j] = 0.0f;
		}
	}
	jm[0][0] = busvoltage[1] * busvoltage[0] * yp->mag[1][0]*sin(yp->ang[1][0] - busangle[1] + busangle[0]) + busvoltage[1]*busvoltage[2]*yp->mag[1][2]*sin(yp->ang[1][2]-busangle[1]+busangle[2]);
	jm[0][1] = -busvoltage[1] * busvoltage[2] * yp->mag[1][2] * sin(yp->ang[1][2] - busangle[1] + busangle[2]);
	jm[0][2] = busvoltage[0] * yp->mag[1][0] * cos(yp->ang[1][0] - busangle[1] + busangle[0]) + 2 * busvoltage[1] * yp->mag[1][1] * cos(yp->ang[1][1]) + busvoltage[2] * yp->mag[1][2] * cos(yp->ang[1][2] - busangle[1] + busangle[2]);
	jm[1][0] = -busvoltage[2] * busvoltage[1] * yp->mag[2][1] * sin(yp->ang[2][1] - busangle[2] + busangle[1]);
	jm[1][1] = busvoltage[2] * busvoltage[0] * yp->mag[2][0]*sin(yp->ang[2][0] - busangle[2] + busangle[0]) + busvoltage[2] * busvoltage[1] * yp->mag[2][1] * sin(yp->ang[2][1] - busangle[2] + busangle[1]);
	jm[1][2] = busvoltage[2] * yp->mag[2][1] * cos(yp->ang[2][1] - busangle[2] + busangle[1]);
	jm[2][0] = busvoltage[1] * busvoltage[0] * yp->mag[1][0] * cos(yp->ang[1][0] - busangle[1] + busangle[0]) + busvoltage[1] * busvoltage[2] * yp->mag[1][2] * cos(yp->ang[1][2] - busangle[1] + busangle[2]);
	jm[2][1] = -busvoltage[1] * busvoltage[2] * yp->mag[1][2] * cos(yp->ang[1][2] - busangle[1] + busangle[2]);
	jm[2][2] = -busvoltage[0] * yp->mag[1][0] * sin(yp->ang[1][0] - busangle[1] + busangle[0]) - 2 * busvoltage[1] * yp->mag[1][1] * sin(yp->ang[1][1]) - busvoltage[2] * yp->mag[1][2] * sin(yp->ang[1][2] - busangle[1] + busangle[2]);
}

void matvec(float *outvec, float **inmat, float *invec)
{
	for (int i = 0; i < NUM_BUSES; i++)
	{
		outvec[i] = 0.0f;
		for (int j = 0; j < NUM_BUSES; j++)
		{
			outvec[i] += inmat[i][j] * invec[j];
		}
	}
}


void cg(float *output, float *deltapower, float **jacobianmatrix, float *tbestimated)
{
	float epsilon = 1e-9;
	float *r = new float[NUM_BUSES];
	float alpha = 0.0f;
	float num = 0.0f;
	float *temp2 = new float[NUM_BUSES];
	float *p = new float[NUM_BUSES];
	float *x = new float[NUM_BUSES];
	float *temp = new float[NUM_BUSES];
	float den = 0.0f;
	float tol = 0.0f;
	float lambda;
	int iter = 100;

	matvec(temp, jacobianmatrix, tbestimated);
#pragma omp parallel for
	for (int i = 0; i < NUM_BUSES; i++)
	{
		r[i] = deltapower[i] - temp[i];
		p[i] = r[i];
		x[i] = tbestimated[i];
	}


	for (int k = 0; k < iter; k++)
	{
		num = 0.0f;
		den = 0.0f;
		matvec(temp2, jacobianmatrix, p);
#pragma omp parallel for
		for (int i = 0; i < NUM_BUSES; i++)
		{
			num += p[i] * r[i];
			den += p[i] * temp2[i];
		}
		alpha = num / den;
#pragma omp parallel for
		for (int i = 0; i < NUM_BUSES; i++)
		{
			output[i] = x[i] + alpha*p[i];

		}
		matvec(temp, jacobianmatrix, output);
		tol = 0.0f;
#pragma omp parallel for
		for (int i = 0; i < NUM_BUSES; i++)
		{
			r[i] = 0.0f;
			r[i] = deltapower[i] - temp[i];
			tol += r[i] * r[i];
		}
		tol = abs(tol);
		//cout << tol << endl;
		if (tol < epsilon)
		{
			cout << "Number of iterations to converge is " << k << endl;
			break;
			
		}
		den = 0.0f;
		num = 0.0f;
		matvec(temp2, jacobianmatrix, p);
#pragma omp parallel for
		for (int i = 0; i < NUM_BUSES; i++)
		{
			num -= r[i] * temp2[i];
			den += p[i] * temp2[i];
		}
		lambda = num / den;
#pragma omp parallel for
		for (int i = 0; i < NUM_BUSES; i++)
		{
			p[i] = r[i] + lambda*p[i];
			x[i] = output[i];
			//cout << p[i] << endl;
		}
		//cout << num << '\t' << den << endl;
		//cout << lambda << endl;
			}
	delete[] r;
	delete[] p;
	delete[] temp;
	delete[] temp2;
}



