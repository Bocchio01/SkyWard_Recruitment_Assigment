#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "matrices.h"

Matrix matInit(int rows, int cols)
{
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double **)malloc(rows * sizeof(double *));
    if (mat.data == NULL)
        exit(EXIT_FAILURE);

    for (int i = 0; i < rows; i++)
    {
        mat.data[i] = (double *)malloc(cols * sizeof(double));
        if (mat.data[i] == NULL)
            exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            mat.data[i][j] = 0;
        }
    }

    return mat;
}

Matrix matTranspose(Matrix A)
{
    Matrix result = matInit(A.cols, A.rows);

    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            result.data[j][i] = A.data[i][j];
        }
    }

    return result;
}

Matrix matScalarMultiply(double scalar, Matrix A)
{
    Matrix result = matInit(A.rows, A.cols);

    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            result.data[i][j] = scalar * A.data[i][j];
        }
    }

    return result;
}

Matrix matSum(Matrix A, Matrix B)
{
    assert(A.rows == B.rows && A.cols == B.cols);

    Matrix result = matInit(A.rows, A.cols);

    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            result.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }

    return result;
}

Matrix matMultiply(Matrix A, Matrix B)
{
    assert(A.cols == B.rows);

    Matrix result = matInit(A.rows, B.cols);

    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < B.cols; j++)
        {
            result.data[i][j] = 0.0;
            for (int k = 0; k < A.cols; k++)
            {
                result.data[i][j] += A.data[i][k] * B.data[k][j];
            }
        }
    }

    return result;
}

Matrix matInverse(Matrix A)
{
    assert(A.rows == A.cols);

    double temp;

    Matrix I = matInit(A.rows, A.cols);
    for (int i = 0; i < A.rows; i++)
    {
        I.data[i][i] = 1;
    }

    for (int k = 0; k < A.rows; k++)
    {
        temp = A.data[k][k];

        for (int j = 0; j < A.cols; j++)
        {
            A.data[k][j] /= temp;
            I.data[k][j] /= temp;
        }

        for (int i = 0; i < A.rows; i++)
        {
            if (i == k)
                continue;
            temp = A.data[i][k];
            for (int j = 0; j < A.cols; j++)
            {
                A.data[i][j] -= A.data[k][j] * temp;
                I.data[i][j] -= I.data[k][j] * temp;
            }
        }
    }

    return I;
}

void matPrint(Matrix A)
{
    printf("Matrix %dx%d:\n", A.rows, A.cols);
    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            printf("%.4f ", A.data[i][j]);
        }
        printf("\n");
    }
}