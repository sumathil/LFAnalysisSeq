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
		for (int j=n-1;j<22;j++)
		{
			J[i][j] = J12[i][j-13];
		}
	}
	for(int i=n-1;i<22;i++)
	{
		for(int j=0;j<n-1;j++)
		{
			J[i][j] = J21[i-13][j];
		}

	}
	for(int i=n-1;i<22;i++)
	{
		for(int j=n-1;j<22;j++)
		{
			J[i][j] = J22[i-13][j-13];
		}

	}

}


void InverseGE(double **Jinv, double **J, int n)
{
	double **temp = new double *[44];
	double d;
	for (int tt=0;tt<44;tt++)
	{
		temp[tt]= new double[44];
	}
	initializematrix(temp,44,44);
	//Assigning the matrix
	
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
		{
			temp[i][j]=J[i][j];
		}
	}
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
	for (int i = n-1; i > 0; i--)
	{
		if (temp[i - 1][0] < temp[i][0])
			for (int j = 0; j <= n * 2; j++)
			{
				d = temp[i][j];
				temp[i][j] = temp[i - 1][j];
				temp[i - 1][j] = d;
			}
	}
	/********** reducing to diagonal  matrix ***********/

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n * 2; j++)
			if (j != i)
			{
				d = temp[j][i] / temp[i][i];
				for (int k = 0; k < n * 2; k++)
					temp[j][k] -= temp[i][k] * d;
			}

	}


	for (int i = 0; i < n; i++)
	{
		d = temp[i][i];
		for (int j = 0; j < n * 2; j++)
			temp[i][j] = temp[i][j] / d;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = n ; j < n * 2; j++)
		{
			//cout << temp[i][j] << "    ";
			Jinv[i][j-n]=temp[i][j];
		}
		//cout << endl;
	}




	for(int tt=0;tt<44;tt++)
	{
		delete[] temp[tt];
	}

	delete[] temp;
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
