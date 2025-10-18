#pragma once

#include <stdint.h>
#include <LSM6.h>


// These variables must be defined in your sketch.
extern LSM6 imu;

int32_t UPDATE_TIME_MS = 10;

// Call this in your setup() to initialize and calibrate the IMU.
void readSensors();

// Returns true if the last update cycle was delayed to more than
// UPDATE_TIME_MS+1 milliseconds.  This could indicate
// computations being too long or interrupts that are delaying
// the loop.
bool balanceUpdateDelayed();
