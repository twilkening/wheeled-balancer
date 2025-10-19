#pragma once

#include <stdint.h>
#include <LSM6.h>
#include "KalmanFilter.h"

// These variables must be defined in your sketch.
extern LSM6 imu;

const int32_t UPDATE_TIME_MS = 10;
const int16_t CALIBRATION_ITERATIONS = 1000;

// Call this in your setup() to initialize and calibrate the IMU.
void sensorsSetup();

// Call this in your loop() to read the sensors and update global variables.
void sensorsUpdate();

// apply Kalman filter and update angle and angleRate variables
void sensorsFilter();

// Returns true if the last update cycle was delayed to more than
// UPDATE_TIME_MS+1 milliseconds.  This could indicate
// computations being too long or interrupts that are delaying
// the loop.
bool sensorsUpdateDelayed();
