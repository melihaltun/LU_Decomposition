/**
* @file LU.h
* @author Melih Altun @2015
**/

#ifndef LU_DECOMP
#define LU_DECOMP

#include<math.h>
#include<string.h>

#if !defined(lin_index)
#define lin_index(i, j, numCol)  ( ((i)*(numCol))+(j) )   //2D to 1D array
#endif
#define TOLR 1E-12  //limit for relative accuracy of numbers. Anything less than that is practically zero.

void LU(double L[], double U[], double A[], int n);

void LUP(double L[], double U[], double P[], double A[], int n);

double detLU(double A[], int n);

#endif
