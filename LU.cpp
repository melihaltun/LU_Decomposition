/**
* @file LU.cpp
* @author Melih Altun @2015
**/

#include "LU.h"

/*subroutine for pivot function*/
void swap(double *a, double *b)
{
	double tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

/*Partial pivoting function to optain permutation matrix P and rotated version of input A*/
void pivot(double P[], double A[], int n)
{
	int i, j, k, max_j; 
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j)
				P[lin_index(i, j, n)] = 1;
			else
				P[lin_index(i, j, n)] = 0;
		}
	}
	for (i = 0; i < n-1; i++) {
		max_j = i;
		for (j = i + 1; j < n; j++) {
			if (fabs(A[lin_index(j, i, n)]) > fabs(A[lin_index(max_j, i, n)]))
				max_j = j;
		}
		if (max_j != i) {
			for (k = 0; k < n; k++) {
				swap(&P[lin_index(i, k, n)], &P[lin_index(max_j, k, n)]);
				swap(&A[lin_index(i, k, n)], &A[lin_index(max_j, k, n)]);
			}
		}
	}
}

/*LU decomposition without permutation matrix.
L and U are lower and upper triangular matrices obtained from input A 
L x U will yield some rotated version of the input.
This implementation will result in diagonal of U having ones as values.
The product of diagonal of L will give the determinant of the input.
Parameters: (outputs) lower triangular matrix, upper triangular matrix,
(inputs) input square matrix, size of input */
void LU(double L[], double U[], double A[], int n)
{
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (j < i)
				L[lin_index(j, i, n)] = 0;
			else {
				L[lin_index(j, i, n)] = A[lin_index(j, i, n)];
				for (k = 0; k < i; k++)
					L[lin_index(j, i, n)] = L[lin_index(j, i, n)] - L[lin_index(j, k, n)] * U[lin_index(k, i, n)];
			}
		}
		for (j = 0; j < n; j++) {
			if (j < i)
				U[lin_index(i, j, n)] = 0;
			else if (j == i)
				U[lin_index(i, j, n)] = 1;
			else {
				if (fabs(L[lin_index(i, i, n)]) < TOLR)
					U[lin_index(i, j, n)] = 0;
				else {
					U[lin_index(i, j, n)] = A[lin_index(i, j, n)] / L[lin_index(i, i, n)];
					for (k = 0; k < i; k++)
						U[lin_index(i, j, n)] = U[lin_index(i, j, n)] - ((L[lin_index(i, k, n)] * U[lin_index(k, j, n)]) / L[lin_index(i, i, n)]);
				}
			}
		}
	}
}

/*LU decompostion with permutation matrix 
L and U are lower and upper triangular matrices obtained from input A
P is the permutation matrix which satisfies: P x A = L x U
Parameters: (outputs) lower triangular matrix, upper triangular matrix, permutation matrix,
(inputs) input square matrix, size of input */
void LUP(double L[], double U[], double P[], double A[], int n)
{
	double *A2;
	A2 = new double[n*n];

	memcpy(A2, A, n*n*sizeof(double));
	pivot(P, A2, n);
	LU(L, U, A2, n);

	delete[] A2;
}

/*Calculates the determinant of a sqaure matrix using LU decomposition.
Parameters: (inputs) sqaure matrix A, size of A. Returns: determinant of A. */
double detLU(double A[], int n)
{
	int i;
	double det=1.0;
	double *L, *U;
	L = new double[n*n];
	U = new double[n*n];

	LU(L, U, A, n);
	for (i = 0; i < n; i++) {
		det *= L[lin_index(i, i, n)];
	}

	delete[] L;
	delete[] U;

	return det;
}