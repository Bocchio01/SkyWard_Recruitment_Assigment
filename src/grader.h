#pragma once

// This function should be called by your code when liftoff is detected
void liftoff();

// This function should be called by your code when the apogee is detected
void apogee();

// This function should be called by your code when landing is detected
void landed();

/**
 * Initializes the system.
 * Sets the initial state, matrices, and constants.
 */
void init();

/**
 * Updates the system state based on the given sensor data.
 * Detects liftoff, apogee, and landing.
 * Datas are processed using a Kalman filter with the gain matrix K approximated as constant.
 *
 * @param acc_x The x component of acceleration.
 * @param acc_y The y component of acceleration.
 * @param acc_z The z component of acceleration.
 * @param gyro_x The x component of gyroscope data.
 * @param gyro_y The y component of gyroscope data.
 * @param gyro_z The z component of gyroscope data.
 * @param baro The barometric pressure data.
 */
void update(float acc_x, float acc_y, float acc_z, float gyro_x, float gyro_y,
            float gyro_z, float baro);