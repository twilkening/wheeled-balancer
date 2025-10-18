#include <Wire.h>
#include "Balance.h"
#include "LQR_PIDController.h"

int32_t gYZero;
int32_t angle; // millidegrees
int32_t angleRate; // degrees/s
bool balanceUpdateDelayedStatus;

bool balanceUpdateDelayed()
{
  return balanceUpdateDelayedStatus;
}

void writeSensors()
{
  // Apply sensitivity gain to gyro readings: 35 mdps/LSB  (for FS = +/-1000 dps)
  angleRate = (imu.g.y - gYZero) / 29; // units: degrees/s

  // angle += angleRate * UPDATE_TIME_MS; // units: millidegrees
  
  // Calculate angle from accelerometer (reference for low frequencies)
  // Only update accel angle when robot is relatively still to avoid motion artifacts
  static int32_t accel_angle = 0;
  int32_t gyro_prediction = angle + angleRate * UPDATE_TIME_MS;

  if (abs(angleRate) < 50) { // Only trust accelerometer when angular velocity is low
    accel_angle = atan2(imu.a.z, imu.a.x) * 57296; // millidegrees
    
    // Integer complementary filter: (252*gyro_prediction + 4*accel_angle) / 256
    // Uses fixed-point arithmetic: alpha = 252/256 ≈ 0.984
    angle = ((gyro_prediction << 8) - (gyro_prediction << 2) + (accel_angle << 2)) >> 8;

  } else {
    angle = gyro_prediction;
  }

  
}

void readSensors()
{
  static uint16_t lastMillis;
  uint16_t ms = millis();
  static uint8_t count = 0;

  // Perform the balance updates at 100 Hz.
  if ((uint16_t)(ms - lastMillis) < UPDATE_TIME_MS) { return; }
  balanceUpdateDelayedStatus = ms - lastMillis > UPDATE_TIME_MS + 1;
  lastMillis = ms;
  
  imu.read();
  writeSensors();
  
}
