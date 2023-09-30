#pragma once

typedef struct
{
    int rows;
    int cols;
    double **data;
} Matrix;

/**
 * Initializes a matrix with the given number of rows and columns.
 *
 * @param rows The number of rows in the matrix.
 * @param cols The number of columns in the matrix.
 * @return A Matrix struct with the specified number of rows and columns and all elements initialized to 0.
 */
Matrix matInit(int rows, int cols);

/**
 * Transposes a given matrix.
 *
 * @param A The matrix to be transposed.
 * @return The transposed matrix.
 */
Matrix matTranspose(Matrix A);

/**
 * Multiplies a matrix by a scalar value.
 *
 * @param scalar The scalar value to multiply the matrix by.
 * @param A The matrix to be multiplied.
 * @return The resulting matrix after scalar multiplication.
 */
Matrix matScalarMultiply(double scalar, Matrix A);

/**
 * Calculates the sum of two matrices.
 *
 * @param A The first matrix to be added.
 * @param B The second matrix to be added.
 * @return The resulting matrix of the addition operation.
 */
Matrix matSum(Matrix A, Matrix B);

/**
 * Multiplies two matrices A and B and returns the resulting matrix.
 *
 * @param A The first matrix to be multiplied.
 * @param B The second matrix to be multiplied.
 * @return The resulting matrix of the multiplication.
 */
Matrix matMultiply(Matrix A, Matrix B);

/**
 * Calculates the inverse of a given matrix using Gauss-Jordan elimination method.
 *
 * @param A The matrix to calculate the inverse of.
 * @return The inverse of the given matrix.
 */
Matrix matInverse(Matrix A);

/**
 * Calculates the pivot matrix of a given matrix using Gauss-Jordan elimination method.
 *
 * @param A The matrix to calculate the pivot matrix of.
 * @return The pivot matrix of the given matrix.
 */
Matrix matComputePivot(Matrix A);

/**
 * Frees the memory allocated for the given matrix.
 *
 * @param A The matrix to be freed.
 */
void matFree(Matrix A);

/**
 * Prints the given matrix A.
 *
 * @param A The matrix to be printed.
 */
void matPrint(Matrix A);