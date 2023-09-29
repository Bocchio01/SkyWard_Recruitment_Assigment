/**
 * @file state_updater.c
 * @brief State updater module
 * @details This module implements the Kalman filter algorithm and performs the event detection.
 * @date 2023-09-29
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include "grader.h"
#include "matrices.h"
#include "state_updater.h"

// #define DEBUG

#define P0 101325    // Standard atmospheric pressure at sea level (Pa)
#define L -0.0065    // Temperature lapse rate (K/m)
#define T0 288.15    // Standard temperature at sea level (K)
#define g 9.80665    // Acceleration due to gravity (m/s^2)
#define M 0.0289644  // Molar mass of Earth's air (kg/mol)
#define RGAS 8.31446 // Universal gas constant (J/(mol·K))

#define dT 0.02 // 50 Hz, supposed knwon from sensor datasheet

#define THRESHOLD_ACC_LIFTOFF 0.5 // Threshold for liftoff detection (m/s^2)
#define THRESHOLD_ACC_LANDING 15  // Threshold for apogee detection (m/s^2)

#define VARIANCE_BARO 2 * 14.6 // Variance of the barometer supposed knwon from sensor datasheet
#define VARIANCE_ACC 2 * 15.7  // Variance of the accelerometer supposed knwon from sensor datasheet

typedef enum
{
    STATE_INIT,
    STATE_PREFLIGHT,
    STATE_FLIGHT,
    STATE_DESCENT,
    STATE_LANDED
} state_t;

typedef enum
{
    ALTITUDE,
    VELOCITY,
    ACCELERATION
} axis_t;

state_t current_state;     // Current state of flight
Matrix estimated_previous; // Estimated state vector at previous time step
Matrix A;                  // State transition matrix
Matrix K;                  // Kalman gain matrix
Matrix H;                  // State to measurement map matrix
Matrix P;                  // Covariance matrix
Matrix Q;                  // Process noise covariance matrix
Matrix R;                  // Measurement noise covariance matrix

void init()
{
    current_state = STATE_INIT;
    estimated_previous = matInit(3, 1);

    A = matInit(3, 3);
    A.data[0][0] = 1;
    A.data[0][1] = dT;
    A.data[0][2] = 0.5 * dT * dT;
    A.data[1][1] = 1;
    A.data[1][2] = dT;
    A.data[2][2] = 1;

    K = matInit(3, 2);

    H = matInit(2, 3);
    H.data[0][0] = 1;
    H.data[1][2] = 1;

    Matrix I = matInit(3, 3);
    I.data[0][0] = 1;
    I.data[1][1] = 1;
    I.data[2][2] = 1;

    P = I;

    Q = matInit(3, 3);
    Q.data[2][2] = 1;

    R = matInit(2, 2);
    R.data[0][0] = VARIANCE_BARO;
    R.data[1][1] = VARIANCE_ACC;

    // Compute Kalman gain matrix K iteratively
    for (int i = 0; i < 20; i++)
    {
        // K = P * H' * inv(H * P * H' + R);
        K = matMultiply(
            matMultiply(P, matTranspose(H)),
            matInverse(matSum(matMultiply(H, matMultiply(P, matTranspose(H))), R)));

        // P = (I - K * H) * P;
        P = matMultiply(matSum(I, matScalarMultiply(-1, matMultiply(K, H))), P);

        // P = A * P * A' + Q;
        P = matSum(matMultiply(A, matMultiply(P, matTranspose(A))), Q);
    }
}

void update(float acc_x, float acc_y, float acc_z, float gyro_x, float gyro_y, float gyro_z, float baro)
{

    Matrix estimated_current = matInit(3, 1); // Estimated state vector at current time step
    Matrix measurement = matInit(2, 1);       // Measurement vector
    measurement.data[0][0] = computeAltitude(baro);
    measurement.data[1][0] = computeAccMagnitude(acc_x, acc_y, acc_z);

    if (current_state == STATE_INIT)
    {
        estimated_previous.data[ALTITUDE][0] = measurement.data[0][0];
        estimated_previous.data[VELOCITY][0] = 0;
        estimated_previous.data[ACCELERATION][0] = measurement.data[1][0];
        current_state = STATE_PREFLIGHT;

        return;
    }

    estimated_current = matMultiply(A, estimated_previous);
    estimated_current = matSum(estimated_current, matMultiply(K, matSum(measurement, matScalarMultiply(-1, matMultiply(H, estimated_current)))));

    switch (current_state)
    {
    case STATE_PREFLIGHT:
        /**
            Liftoff is detected when the absolute difference between the current estimated acceleration and the gravity acceleration exceeds the liftoff acceleration threshold.
            Based on accelerometer data only.
        */
        if (fabs(estimated_current.data[ACCELERATION][0] - g) > THRESHOLD_ACC_LIFTOFF)
        {
            liftoff();
            current_state = STATE_FLIGHT;
        }
        break;

    case STATE_FLIGHT:
        /**
            Apogee is detected when the rocket starts to descend.
            Based on barometer data only.
        */
        if (estimated_current.data[ALTITUDE][0] < estimated_previous.data[ALTITUDE][0])
        {
            apogee();
            current_state = STATE_DESCENT;
        }
        break;

    case STATE_DESCENT:
        /**
            Landing is detected when the absolute difference between the current estimated acceleration and the gravity acceleration exceeds the landing acceleration threshold.
            Based on accelerometer data only.
        */
        if (fabs(estimated_current.data[ACCELERATION][0] - estimated_previous.data[ACCELERATION][0]) > THRESHOLD_ACC_LANDING)
        {
            landed();
            current_state = STATE_LANDED;
        }
        break;
    }

    estimated_previous = estimated_current;

#ifdef DEBUG
    printf("Estimated: %f,%f,%f\n",
           estimated_current.data[ALTITUDE][0],
           estimated_current.data[VELOCITY][0],
           estimated_current.data[ACCELERATION][0]);
#endif
}

double computeAltitude(double baro)
{
    // From the combination of the ideal gas law and the state equation for steady heavy fluids
    return -(RGAS * T0) / (g * M) * log(baro / P0);
}

double computeAccMagnitude(double acc_x, double acc_y, double acc_z)
{
    return sqrt(acc_x * acc_x + acc_y * acc_y + acc_z * acc_z);
}
