#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "grader.h"

#define P0 101325    // Standard atmospheric pressure at sea level (Pa)
#define L -0.0065    // Temperature lapse rate (K/m)
#define T0 288.15    // Standard temperature at sea level (K)
#define g 9.80665    // Acceleration due to gravity (m/s²)
#define M 0.0289644  // Molar mass of Earth's air (kg/mol)
#define RGAS 8.31446 // Universal gas constant (J/(mol·K))

#define dt 0.02 // 50 Hz, supposed knwon from sensor datasheet

#define ACC_LIFTOFF_TRESHOLD 0.5 // Threshold for liftoff detection
#define ACC_LANDING_TRESHOLD 6   // Threshold for apogee detection

typedef enum
{
    STATE_INIT,
    STATE_PREFLIGHT,
    STATE_FLIGHT,
    STATE_DESCENT,
    STATE_LANDED
} state_t;

typedef struct
{
    int rows;
    int cols;
    double **data;
} Matrix;

state_t state_update;
Matrix state_prev; // State vector
Matrix A;          // State transition matrix
Matrix K;          // Kalman gain matrix
Matrix H;          // State to measurement map matrix
Matrix P;          // Covariance matrix
Matrix Q;          // Process noise covariance matrix
Matrix R;          // Measurement noise covariance matrix

// Function to initialize a matrix
Matrix matInit(int rows, int cols)
{
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
    {
        mat.data[i] = (double *)malloc(cols * sizeof(double));
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

// Function to transpose a matrix
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

// Function to multiply a matrix by a scalar
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

// Function to sum two matrices
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

// Function to multiply two matrices
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

// Function to perform matrix inversion (Gauss-Jordan elimination)
Matrix matInverse(Matrix A)
{
    assert(A.rows == A.cols);

    Matrix I = matInit(A.rows, A.cols);

    // Initialize I as the identity matrix
    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            I.data[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    int k, i, j;
    double temp;

    for (k = 0; k < A.rows; k++)
    {
        temp = A.data[k][k];

        for (j = 0; j < A.cols; j++)
        {
            A.data[k][j] /= temp;
            I.data[k][j] /= temp;
        }

        for (i = 0; i < A.rows; i++)
        {
            if (i == k)
                continue;
            temp = A.data[i][k];
            for (j = 0; j < A.cols; j++)
            {
                A.data[i][j] -= A.data[k][j] * temp;
                I.data[i][j] -= I.data[k][j] * temp;
            }
        }
    }

    return I;
}

// Function to calculate acceleration magnitude
double computeAccMagnitude(double acc_x, double acc_y, double acc_z)
{
    return sqrt(acc_x * acc_x + acc_y * acc_y + acc_z * acc_z);
}

// Function to calculate altitude
double computeAltitude(double baro)
{
    // Here the latter is the Taylor expansion of the former.
    // return -(RGAS * T0)/(g * M) * log(baro / P0);
    return (T0 / L) * (1 - pow((baro / P0), ((RGAS * L) / (g * M))));
}

// Function to print a matrix
void matPrint(Matrix A)
{
    printf("Matrix %dx%d:\n", A.rows, A.cols);
    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.cols; j++)
        {
            printf("%.2f ", A.data[i][j]);
        }
        printf("\n");
    }
}

void init()
{
    state_update = STATE_INIT;
    state_prev = matInit(3, 1);

    A = matInit(3, 3);
    A.data[0][0] = 1;
    A.data[0][1] = dt;
    A.data[0][2] = 0.5 * dt * dt;
    A.data[1][1] = 1;
    A.data[1][2] = dt;
    A.data[2][2] = 1;

    K = matInit(3, 2); // Approximated as constant

    Matrix I = matInit(3, 3);
    I.data[0][0] = 1;
    I.data[1][1] = 1;
    I.data[2][2] = 1;

    P = I;

    H = matInit(2, 3);
    H.data[0][0] = 1;
    H.data[1][2] = 1;

    R = matInit(2, 2);
    R.data[0][0] = 100;
    R.data[1][1] = 100;

    Q = matInit(3, 3);
    Q.data[2][2] = 1;

    for (int i = 0; i < 20; i++)
    {
        // printf("Iteration %d\n", i);
        // printf("P");
        // matPrint(P);
        // printf("matTranspose(H)");
        // matPrint(matTranspose(H));

        Matrix PH_t = matMultiply(P, matTranspose(H));
        // printf("PH_t");
        // matPrint(PH_t);

        Matrix HPH_t = matMultiply(H, PH_t);
        // printf("HPH_t");
        // matPrint(HPH_t);

        Matrix HPH_t_R = matSum(HPH_t, R);
        // printf("HPH_t_R");
        // matPrint(HPH_t_R);

        // printf("inverse HPH_t_R");
        // matPrint(matInverse(HPH_t_R));
        K = matMultiply(PH_t, matInverse(HPH_t_R));
        // printf("K");
        // matPrint(K);

        P = matMultiply(matSum(I, matScalarMultiply(-1, matMultiply(K, H))), P);
        P = matSum(matMultiply(A, matMultiply(P, matTranspose(A))), Q);
    }

    printf("K\t");
    matPrint(K);
}

void update(float acc_x, float acc_y, float acc_z, float gyro_x, float gyro_y, float gyro_z, float baro)
{

    Matrix state_estimated = matInit(3, 1); // State vector estimated
    Matrix measurement = matInit(2, 1);     // Measurement vector
    measurement.data[0][0] = computeAltitude(baro);
    measurement.data[1][0] = computeAccMagnitude(acc_x, acc_y, acc_z);

    if (state_update == STATE_INIT)
    {
        state_prev.data[0][0] = measurement.data[0][0];
        state_prev.data[1][0] = 0;
        state_prev.data[2][0] = measurement.data[1][0];
        state_update = STATE_PREFLIGHT;

        return;
    }

    state_estimated = matMultiply(A, state_prev);
    state_estimated = matSum(state_estimated, matMultiply(K, matSum(measurement, matScalarMultiply(-1, matMultiply(H, state_estimated)))));

    printf("Estimated: %f,%f,%f\n", state_estimated.data[0][0], state_estimated.data[1][0], state_estimated.data[2][0]);

    switch (state_update)
    {
    case STATE_PREFLIGHT:
        if (abs(state_estimated.data[2][0] - state_prev.data[2][0]) > ACC_LIFTOFF_TRESHOLD)
        {
            liftoff();
            state_update = STATE_FLIGHT;
        }
        break;

    case STATE_FLIGHT:
        if (state_estimated.data[0][0] < state_prev.data[0][0])
        {
            apogee();
            state_update = STATE_DESCENT;
        }
        break;

    case STATE_DESCENT:
        if (abs(state_estimated.data[2][0] - state_prev.data[2][0]) > ACC_LANDING_TRESHOLD)
        {
            landed();
            state_update = STATE_LANDED;
        }
        break;
    }

    state_prev = state_estimated;
}