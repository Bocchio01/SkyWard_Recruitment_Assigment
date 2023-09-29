#pragma once

/**
 * Computes the magnitude of the acceleration vector given its x, y, and z components.
 *
 * @param acc_x The x component of the acceleration vector.
 * @param acc_y The y component of the acceleration vector.
 * @param acc_z The z component of the acceleration vector.
 * @return The magnitude of the acceleration vector.
 */
double computeAccMagnitude(double acc_x, double acc_y, double acc_z);

/**
 * Computes the altitude based on the given barometric pressure.
 *
 * @param baro The barometric pressure in hPa.
 * @return The altitude in meters.
 */
double computeAltitude(double baro);