#include "LF_function.h"


void bustype(int *PV, int *PQ, int *type_bus)
{
	int j = 0;
	int k = 0;
	for (int i = 0; i < num_buses; i++)
	{
		if (type_bus[i] == 2 || type_bus[i] == 1)
		{
			PV[j] = i;
			j++;
		}
		if (type_bus[i] == 3)
		{
			PQ[k] = i;
			k++;
		}
	}
}


void initializevector( double *vect, int n)
{
	for (int i=0;i<n;i++)
	{
		vect[i]=0.0f;
	}

}

void printvector(double *vect,int n)
{
	for(int i=0;i<n;i++)
	{
		cout<< vect[i]<<endl;
	}


}


void initializematrix( double **mat, int nrows, int ncols)
{
	for (int i=0;i<nrows;i++)
	{
		for(int j=0;j<ncols;j++)
		{
			mat[i][j]=0.0f;
		}
	}

}


void printmatrix( double **mat, int nrows, int ncols)
{
	for (int i=0;i<nrows;i++)
	{
		for(int j=0;j<ncols;j++)
		{
			cout<<mat[i][j]<<"\t";
		}
		cout<<endl;
	}

}

void Jacobian(double **J, double **J11, double **J12, double **J21, double **J22, int n)
{
	for(int i=0;i<n-1;i++)
	{
		for (int j=0;j<n-1;j++)
		{
			J[i][j] = J11[i][j];
		}
	}
	for(int i=0;i<n-1;i++)
	{
		for (int j=n-1;j<53;j++)
		{
			J[i][j] = J12[i][j-(n-1)];
		}
	}
	for(int i=n-1;i<53;i++)
	{
		for(int j=0;j<n-1;j++)
		{
			J[i][j] = J21[i-(n-1)][j];
		}

	}
	for(int i=n-1;i<53;i++)
	{
		for(int j=n-1;j<53;j++)
		{
			J[i][j] = J22[i-(n-1)][j-(n-1)];
		}

	}

}


/*void InverseGE(double **Jinv, double **J, int n)
{
	double **temp = new double *[106];
	double d;
	for (int tt=0;tt<106;tt++)
	{
		temp[tt]= new double[106];
	}
	initializematrix(temp,106,106);
	//printmatrix(J,53,53);
	//Assigning the matrix
	
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
		{
			temp[i][j]=J[i][j];
		}
	}
	//printmatrix(temp,53,53);
	// Gauss Elimination method	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2 * n; j++)
		{
			if (j == (i + n))
			{
				temp[i][j] = 1;
			}
		}
	}
	//printmatrix(temp,53,106);
	for (int i = n-1; i > 0; i--)
	{
		if (temp[i - 1][0] < temp[i][0])
		{
			for (int j = 0; j < n * 2; j++)
			{
				d = temp[i][j];
				temp[i][j] = temp[i - 1][j];
				temp[i - 1][j] = d;
			}
		}	
	}
	printmatrix(temp,53,106);
	/********** reducing to diagonal  matrix ***********/

/*	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2*n ; j++)
			if (j != i)
			{
				//cout<< temp[i][i] <<" is temp diagonal element" <<endl;
				d = temp[j][i] / temp[i][i];
				for (int k = 0; k < n * 2; k++)
					temp[j][k] -= temp[i][k] * d;
			}

	}
	//printmatrix(temp,53,106);

	for (int i = 0; i < n; i++)
	{
		d = temp[i][i];
		for (int j = 0; j < n * 2; j++)
			temp[i][j] = temp[i][j] / d;
	}
//	printmatrix(temp,106,106);

	for (int i = 0; i < n; i++)
	{
		for (int j = n ; j < n * 2; j++)
		{
			//cout << temp[i][j] << "    ";
			Jinv[i][j-n]=temp[i][j];
		}
		//cout << endl;
	}




	for(int tt=0;tt<106;tt++)
	{
		delete[] temp[tt];
	}

	delete[] temp;
}
*/

void InverseGE(double **Jinv, double **J, int n)
{
	double **temp = new double *[53];
	double **s = new double *[53];
	double *d = new double[53];
	double xtemp;
	double *y = new double[53];
	for (int tt = 0; tt<53; tt++)
	{
		temp[tt] = new double[53];
	}
	for (int tt = 0; tt<53; tt++)
	{
		s[tt] = new double[53];
	}
	initializematrix(temp, 53, 53);
	initializevector(y, 53);
	

	for (int i = 0; i<n; i++)
	{
		for (int j = 0; j<n; j++)
		{
			temp[i][j] = J[i][j];
		}
	}

	LU(temp, n);
	
	for (int m = 0; m < n; m++)
	{
		initializevector(d, 53);
		d[m] = 1.0;
		for (int i = 0;i < n; i++)
		{
			xtemp = 0.0;
			for (int j = 0; j <= i - 1; j++)
			{
				xtemp = xtemp + temp[i][j] * y[j];
			}
			y[i] = d[i] - xtemp;
		}
		for (int i = n - 1; i >= 0; i--)
		{
			xtemp = 0.0f;
			for (int j = i + 1; j < n; j++)
			{
				xtemp = xtemp + temp[i][j] * s[j][m];

			}
			s[i][m] = (y[i] - xtemp) / temp[i][i];
		}

	}
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Jinv[i][j] = s[i][j];
		}
	}

	for (int tt = 0; tt<53; tt++)
	{
		delete[] temp[tt];
	}

	delete[] temp;
	delete[] d;
	delete[] y;
	for (int tt = 0; tt<53; tt++)
	{
		delete[] s[tt];
	}

	delete[] s;
}


void LU(double **temp, int n)
{
	int i, j, k, m, an, am;
	double x;
	
	for (k = 0; k < n - 1; k++)
	{
		for (j = k + 1; j < n; j++)
		{
			x = temp[j][k] / temp[k][k];
			for (i = k; i < n; i++)
			{
				temp[j][i] = temp[j][i] - x*temp[k][i];
			}
			temp[j][k] = x;
		}
	}
	
}






void matvec(double *outvec, double **inmat, double *invec,int n)
{
        for (int i = 0; i < n; i++)
        {
                outvec[i] = 0.0f;
                for (int j = 0; j < n; j++)
                {
                        outvec[i] += inmat[i][j] * invec[j];
                }
        }
}

void pol2deg(double *out,double *in,int n)
{
 for(int i =0;i<n;i++)
        {
        out[i] = (in[i]*180)/pi;
        }

}


/*void loadflow(int num_buses,double *busvoltage,double *busangle,int Bmva)
{




}*/
