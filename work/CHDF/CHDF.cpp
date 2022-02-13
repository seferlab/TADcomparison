#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>

using namespace std;

int DomainCaller(int Length, int number, double** matrice, int* domain, int Size)
{
	double **Intra, **Extra, **T, **D, **R, **E_sum, **I_sum_square, **E_sum_square;
	
	Intra = new double*[Length];
	Extra = new double*[Length];
	D = new double*[Length];
	T = new double*[Length];
	R = new double*[Length];
	I_sum_square = new double*[Length];
	E_sum = new double*[Length];
	E_sum_square = new double*[Length];

	int i, j, k;
	for (i = 0; i<Length; i++)
	{
		
		Intra[i] = new double[Length];
		Extra[i] = new double[Length];
		T[i] = new double[Length];
		D[i] = new double[Length];
		R[i] = new double[Length];
		E_sum[i] = new double[Length];
		I_sum_square[i] = new double[Length];
		E_sum_square[i] = new double[Length];
	}


	for (i = 0; i<Length; i++)
	{
		for (j = 0; j<Length; j++)
		{
		
			Intra[i][j] = -1E100;
			Extra[i][j] = 0;
			T[i][j] = -1E100;
			D[i][j] = -1E100;
			R[i][j] = -1E100;
			E_sum[i][j] = -1E100;
			I_sum_square[i][j] = -1E100;
			E_sum_square[i][j] = -1E100;
		}
	}


	//estimate mu of R
	double sum;

	int Li = (Length / 4), Lj = 3 * Li;
	sum = 0;
	for (i = 0; i<Li; i++)
	{
		for (j = Lj + i; j<Length; j++)
		{
			sum = sum + matrice[i][j];
		}
	}
	double mu_hat = sum / ((Li + 1)*Li / 2);
	

	//Intra£¬sum of square error in a domain

	double mu, sum_square, N;

	T[0][0] = 0;
	D[0][0] = 0;
	E_sum[0][0] = 0;
	I_sum_square[0][0] = 0;
	Intra[0][0] = 0;

	for (k = 1; k<Length; k++)
	{
		T[0][k] = T[0][(k - 1)] + matrice[0][k];
		E_sum[0][k] = E_sum[0][(k - 1)] + pow(matrice[0][k], 2);
		D[k][k] = 0;
		I_sum_square[k][k] = 0;
		sum = 0; sum_square = 0;
		for (i = 0; i < k; i++)
		{
			sum = sum + matrice[i][k];
			sum_square = sum_square + pow(matrice[i][k], 2);
		}
		T[k][k] = T[(k - 1)][(k - 1)] + sum;
		E_sum[k][k] = E_sum[(k - 1)][(k - 1)] + sum_square;
		D[0][k] = T[k][k];
		I_sum_square[0][k] = E_sum[k][k];
		mu = D[0][k] / ((pow(k + 1, 2) - k - 1) / 2);
		Intra[0][k] = -(I_sum_square[0][k] - D[0][k] * mu);
		Intra[k][k] = 0;
	}


	//Extra£¬sum of square error outside a domain
	for (i = 1; i<(Length - 1); i++)
	{
		for (j = (i + 1); j<Length; j++)
		{
			sum = 0; sum_square = 0;
			for (k = 0; k <= i; k++)
			{
				sum = sum + matrice[k][j];
				sum_square = sum_square + pow(matrice[k][j], 2);
			}
			T[i][j] = T[i][(j - 1)] + sum;
			R[i][j] = T[(i - 1)][j] - T[(i - 1)][(i - 1)];
			D[i][j] = T[j][j] - T[(i - 1)][j];
			E_sum[i][j] = E_sum[i][(j - 1)] + sum_square;
			E_sum_square[i][j] = E_sum[(i - 1)][j] - E_sum[(i - 1)][(i - 1)];
			I_sum_square[i][j] = E_sum[j][j] - E_sum[(i - 1)][j];
			N = (j - i + 1)*(j - i) / 2;

			mu = D[i][j] / N;
			Intra[i][j] = -(I_sum_square[i][j] - D[i][j] * mu);

			N = i*(j - i + 1);
			Extra[i][j] = -(E_sum_square[i][j] + pow(mu_hat, 2)*N - 2 * mu_hat * R[i][j]);
		}
		sum = 0; sum_square = 0;
		for (k = 0; k <= i - 1; k++)
		{
			sum = sum + matrice[k][i];
			sum_square = sum_square + pow(matrice[k][i], 2);
		}
		N = i;
		Extra[i][i] = -(sum_square + pow(mu_hat, 2)*N - 2 * mu_hat * sum);

	}

	//release memory
	for (i = 0; i<Length; i++){
		delete[]E_sum[i];
		delete[]E_sum_square[i];
	}
	delete[]E_sum;
	delete[]E_sum_square;



	//dynamic planning
	double *vecteur, **I, **t_matrice;
	vecteur = new double[Length - 1];
	I = new double*[number];
	t_matrice = new double*[number - 1];
	for (i = 0; i<number - 1; i++)
	{
		*(I + i) = new double[Length];
		*(t_matrice + i) = new double[Length];
	}
	*(I + number - 1) = new double[Length];




	for (i = 0; i<number - 1; i++)
	{
		for (j = 0; j<Length; j++)
		{
			I[i][j] = -1E100;
			t_matrice[i][j] = -1;
		}
	}

	for (j = 0; j<Length; j++)
	{
		I[number - 1][j] = -1E100;
	}

	for (j = 0; j<Length - 1; j++)
	{
		vecteur[j] = -1E100;
		I[0][j] = Intra[0][j];
	}
	I[0][Length - 1] = Intra[0][Length - 1];


	int index, kk, l, u, uu, sec;
	double max, n_I_sum_square, n_D, neighbor, S;
	for (k = 1; k<number; k++)
	{
		kk = ((k+1) * Size) < (Length) ? ((k+1) * Size) : (Length);
    	for (l = k; l<kk; l++)
		{

			for (i = 0; i<(Length - 1); i++)
			{
				vecteur[i] = -1E100;
			}

			uu = (l - Size - 1) > 1 ? (l - Size - 1) : 1;
			for (u = uu; u <= l; u++)
			{
				if (k>1)
				{
					sec = (int)t_matrice[k - 2][u - 1] + 1;
					n_I_sum_square = I_sum_square[sec][l] - I_sum_square[sec][u - 1] - I_sum_square[u][l];
					n_D = D[sec][l] - D[sec][u - 1] - D[u][l];
					neighbor = -(n_I_sum_square - n_D*n_D / (l - u + 1) / (u - sec));

					if (l > u)
					{
						S = D[u][l] / ((l - u + 1)*(l - u) / 2) - n_D / (l - u + 1) / (u - sec);
					}
					else
					{
						S = 0 - n_D / (l - u + 1) / (u - sec);
					}

					vecteur[u - 1] = I[k - 1][u - 1] + Intra[u][l] + (Extra[sec][l] - Extra[sec][u - 1]) + neighbor + S*fabs(S)*(l - u + 1);

				}
				else
				{
					n_I_sum_square = I_sum_square[0][l] - I_sum_square[0][u - 1] - I_sum_square[u][l];
					n_D = D[0][l] - D[0][u - 1] - D[u][l];
					neighbor = -(n_I_sum_square - n_D*n_D / (l - u + 1) / (u - 0));
					if (l>u)
					{
						S = D[u][l] / ((l - u + 1)*(l - u) / 2) - n_D / (l - u + 1) / (u - 0);
					}
					else
					{
						S = 0 - n_D / (l - u + 1) / (u - 0);
					}


					vecteur[u - 1] = I[k - 1][u - 1] + Intra[u][l] + neighbor + S*fabs(S)*(l - u + 1);
				}



			}
			index = 0;
			max = vecteur[0];
			for (u = 0; u<Length - 1; u++)
			{
				if (vecteur[u]>max)
				{
					max = vecteur[u];
					index = u;
				}
			}
			I[k][l] = max;
			t_matrice[k - 1][l] = index;
		}
	}



	//get the best domains
	max = I[0][Length - 1];
	for (i = 0; i<number; i++)
	{
		if (I[i][Length - 1]>max)
		{
			index = i;
			max = I[i][Length - 1];
		}
	}
	for (i = 0; i<number; i++)
	{
		domain[i] = -1;
	}


	domain[index] = Length - 1;
	for (k = index - 1; k >= 0; k--)
	{
		domain[k] = (int)t_matrice[k][(int)domain[k + 1]];
	}


	//release memory
	for (i = 0; i<Length; i++){
		delete[]Intra[i];
		delete[]Extra[i];
		delete[]T[i];
		delete[]D[i];
		delete[]R[i];
		delete[]I_sum_square[i];
		
		if (i<number - 1)
		{
			delete[]t_matrice[i];
			delete[]I[i];
		}
	}


	delete[]Intra;
	delete[]Extra;
	delete[]t_matrice;
	delete[]I;
	delete[]vecteur;
	delete[]T;
	delete[]D;
	delete[]R;
	delete[]I_sum_square;

	return 1;
}


int main(int argc, char *argv[])
{
	char input[50], output[50];  //input and output file names
	strcpy(input,argv[1]);
	strcpy(output,argv[2]);

	//Length is the dimension of the Hi-C matrix. number is the maximum domain numbers. Size is the maximum domain size.
	int  Length,number,Size;  
	Length = atoi(argv[3]);
	number = atoi(argv[4]);
	Size = atoi(argv[5]);


	int *domain;
	double** matrice;
	matrice = new double*[Length];
	domain = new  int[number];
	for (int i = 0; i < Length; i++)
	{
		matrice[i] = new double[Length];
	}

	//read matrix
	fstream file1(input);
	for (int i = 0; i < Length; i++)
	{
		for (int j = 0; j < Length; j++)
		{
			file1 >> matrice[i][j];
		}
	}

	//call domains
	DomainCaller(Length, number, matrice, domain, Size);

	
	//output
	ofstream file2(output);
	file2 << 1 << "\t" << domain[0]+1 << endl;
	for (int i = 1; i < number; i++)
	{
		if (domain[i]>0)
		{
			file2 << domain[i - 1] + 2 << "\t" << domain[i] + 1 << endl;
		}
		
	}
	return 1;
}

