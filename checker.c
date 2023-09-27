#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "grader.h"

#define P0 101325   // Standard atmospheric pressure at sea level (Pa)
#define L -0.0065   // Temperature lapse rate (K/m)
#define T0 288.15   // Standard temperature at sea level (K)
#define g 9.80665   // Acceleration due to gravity (m/s²)
#define M 0.0289644 // Molar mass of Earth's air (kg/mol)
#define R 8.31446   // Universal gas constant (J/(mol·K))

#define dt 0.02 // 50 Hz, supposed knwon from sensor datasheet

typedef struct
{
    int rows;
    int cols;
    double **data;
} Matrix;

bool existPrevState;
Matrix state_prev; // State vector
Matrix A;          // State transition matrix
Matrix K;          // Kalman gain matrix
Matrix H;          // State to measurement map matrix

// Function to initialize a matrix
Matrix initializeMatrix(int rows, int cols)
{
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        mat.data[i] = (double *)malloc(cols * sizeof(double));
    }
    return mat;
}

// Function to sum two matrices
Matrix matrixSum(Matrix A, Matrix B)
{
    assert(A.rows == B.rows && A.cols == B.cols);

    Matrix result = initializeMatrix(A.rows, A.cols);

    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            result.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }

    return result;
}

// Function to subtract two matrices
Matrix matrixSubtract(Matrix A, Matrix B)
{
    assert(A.rows == B.rows && A.cols == B.cols);

    Matrix result = initializeMatrix(A.rows, A.cols);

    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            result.data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }

    return result;
}

// Function to multiply two matrices
Matrix matrixMultiply(Matrix A, Matrix B)
{
    assert(A.cols == B.rows);

    Matrix result = initializeMatrix(A.rows, B.cols);

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

// Function to calculate acceleration magnitude
double accelerationMagnitude(double acc_x, double acc_y, double acc_z)
{
    return sqrt(acc_x * acc_x + acc_y * acc_y + acc_z * acc_z);
}

// Function to calculate altitude
double altitude(double baro)
{
    return (T0 / L) * (1 - pow((baro / P0), ((R * L) / (g * M))));
}

void init()
{
    printf("Init!\n");
    existPrevState = false;

    state_prev = initializeMatrix(3, 1);

    A = initializeMatrix(3, 3);
    A.data[0][0] = 1;
    A.data[0][1] = dt;
    A.data[0][2] = 0.5 * dt * dt;
    A.data[1][0] = 0;
    A.data[1][1] = 1;
    A.data[1][2] = dt;
    A.data[2][0] = 0;
    A.data[2][1] = 0;
    A.data[2][2] = 1;

    K = initializeMatrix(3, 2); // Approximated as constant. Precomputed with MATLAB
    K.data[0][0] = 0.00962494879450118;
    K.data[0][1] = 0.00144078029023556;
    K.data[1][0] = 0.00399699353157653;
    K.data[1][1] = 0.0130627832358456;
    K.data[2][0] = 0.00144078029023556;
    K.data[2][1] = 0.0918483237143149;

    H = initializeMatrix(2, 3);
    H.data[0][0] = 1;
    H.data[0][1] = 0;
    H.data[0][2] = 0;
    H.data[1][0] = 0;
    H.data[1][1] = 0;
    H.data[1][2] = 1;
}

void update(float acc_x, float acc_y, float acc_z, float gyro_x, float gyro_y, float gyro_z, float baro)
{

    Matrix state_estimated = initializeMatrix(3, 1); // State vector estimated
    Matrix measurement = initializeMatrix(2, 1);     // Measurement vector
    measurement.data[0][0] = altitude(baro);
    measurement.data[1][0] = accelerationMagnitude(acc_x, acc_y, acc_z);

    printf("Update: %f,%f,%f,%f,%f,%f,%f\n", acc_x, acc_y, acc_z, gyro_x, gyro_y, gyro_z, baro);

    if (existPrevState == false)
    {
        state_prev.data[0][0] = measurement.data[0][0];
        state_prev.data[1][0] = 0;
        state_prev.data[2][0] = measurement.data[1][0];
        existPrevState = true;

        return;
    }

    state_estimated = matrixMultiply(A, state_prev);
    state_estimated = matrixSum(state_estimated, matrixMultiply(K, matrixSubtract(measurement, matrixMultiply(H, state_estimated))));

    printf("Estimated: %f,%f,%f\n", state_estimated.data[0][0], state_estimated.data[1][0], state_estimated.data[2][0]);

    state_prev = state_estimated;
}