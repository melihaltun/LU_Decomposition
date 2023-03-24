# LU_decomposition
Lower upper matrix decomposition and determinant algorithm.

The method factorizes a square matrix into one upper(U) and one lower(L) triangular matrix. The method is useful in solving linear systems of equations, calculating the determininants and matrix inversions.

The source code provides three functions:

LU function outputs L and U matrices for a given square matrix A and size n.
L x U will yield some rotated form of A

LUP functions outputs L, U and P matrices for a given square matrix input A and size n.
P matrix is the permutation matrix that defines how the matrix A is rotated when L and U are calculated. The matrices satisfy the condition:
L x U = P x A

detLU function calculates determinant of a square matrix A using LU decomposition. Complexity is O(N^3) which beats the classical recursive solution with the time complexity O(N!)
