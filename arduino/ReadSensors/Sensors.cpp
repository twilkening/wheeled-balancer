#include <Wire.h>
#include "Sensors.h"
#include "KalmanFilter.h"

int32_t gYZero;
int32_t angle; // millidegrees
int32_t angleRate; // degrees/s
KalmanFilter gyroKalmanFilter;
KalmanFilter accKalmanFilter;
float gyroP[3][3];
float accP[3][3];
float P[3][3];
bool sensorsUpdateDelayedStatus;

bool sensorsUpdateDelayed()
{
  return sensorsUpdateDelayedStatus;
}

void sensorsSetup()
{
  // Initialize IMU.
  Wire.begin();
  if (!imu.init())
  {
    while(true)
    {
      Serial.println("Failed to detect and initialize IMU!");
      delay(200);
    }
  }
  imu.enableDefault();
  imu.writeReg(LSM6::CTRL2_G, 0b01011000); // 208 Hz, 1000 deg/s

  // Wait for IMU readings to stabilize.
  delay(1000);

  // Calibrate the gyro offset and initial angle
  // assume we are at rest
  int32_t totalGyro = 0;
  int32_t totalAcc = 0;
  for (int i = 0; i < CALIBRATION_ITERATIONS; i++)
  {
    imu.read();
    totalGyro += imu.g.y;

    // It's really calm, so use the accelerometer to measure the
    // robot's rest angle.  The atan2 function returns a result
    // in radians, so we multiply it by 180000/pi to convert it
    // to millidegrees.
    totalAcc = atan2(imu.a.z, imu.a.x) * 57296; // millidegrees
    
    delay(2);
  }
  gYZero = totalGyro / CALIBRATION_ITERATIONS;
  angle = totalAcc / CALIBRATION_ITERATIONS;
  angleRate = 0; // assume we are at rest

  // setup Kalman Filter matrices
  float dt = UPDATE_TIME_MS / 1000.0f; // convert ms to s
  float x_init_gyro[3] = { (float)angle, 0.0f, (float)gYZero }; // initial state: [angle, angleRate, gyroBias]
  float x_init_acc[3] = { (float)angle, 0.0f, (float)gYZero }; // initial state: [angle, angleRate, gyroBias]
  float R_gyro = 1e-1; // measurement noise variance
  float R_acc = 1; // measurement noise variance
  float Q[3] = { 1e-2, 1e-2, 1e-4 }; // process noise variances for [angle, angleRate, gyroBias]
  float P_init[3][3] =  { {1e-3, 0.0, 0.0},
                          {0.0, 1e-2, 0.0},
                          {0.0, 0.0, 1e-4} };
  float H_gyro[3] = { 0.0, 1.0, -1.0 }; // measurement matrix for gyro
  float H_acc[3] = { 1.0, 0.0, 0.0 }; // measurement matrix for accelerometer
  gyroKalmanFilter.initialize(dt, x_init_gyro, Q, R_gyro, P_init, H_gyro);
  accKalmanFilter.initialize(dt, x_init_acc, Q, R_acc, P_init, H_acc);

}


void sensorsFilterUpdate(uint16_t ms)
{
  // read from the gyro + accelerometer
  // -- Apply sensitivity gain to gyro readings: 35 mdps/LSB  (for FS = +/-1000 dps)
  float gyroRate = (imu.g.y - gYZero) / 29; // units: degrees/s
  float gyroAngle = angle + gyroRate * ms; // units: millidegress; difference equation integration
  float accAngle = (atan2(imu.a.z, imu.a.x) * 57296) / 1000.0f; // degrees

  // pass readings into the kalman filter objects
  gyroKalmanFilter.update(gyroRate);
  accKalmanFilter.update(accAngle);
  gyroKalmanFilter.predict();
  accKalmanFilter.predict();

  // get P and xhat estimates out of the kalman filter objects
  float* gyroState = gyroKalmanFilter.getState();
  float* accState = accKalmanFilter.getState();
  float* gyroP = gyroKalmanFilter.getCovariance();
  float* accP = accKalmanFilter.getCovariance();

  // print states for debugging using a single print call, at 20Hz
  static uint16_t lastPrintMillis = 0;
  // issues with printing floats on Arduino - workaround is to use 
  // if ((uint16_t)(ms - lastPrintMillis) > 50) {
  //   static char debugBuffer[175];  // Static buffer to avoid repeated allocation
  //   snprintf(debugBuffer, sizeof(debugBuffer),
  //           "Time(ms): %u, Raw: %.3f, %.3f, %.3f, Gyro State: %.3f, %.3f, %.3f, Acc State: %.3f, %.3f, %.3f\n",
  //           ms, accAngle, gyroRate, gYZero,
  //           gyroState[0], gyroState[1], gyroState[2],
  //           accState[0], accState[1], accState[2]);
  //   Serial.print(debugBuffer);
  //   lastPrintMillis = ms;
  // }
  if ((uint16_t)(ms - lastPrintMillis) > 50) {
    Serial.print("Time(ms): "); Serial.print(ms);
    Serial.print(", Raw: "); Serial.print(accAngle, 3);
    Serial.print(", "); Serial.print(gyroRate, 3);
    Serial.print(", "); Serial.print(gYZero);
    Serial.print(", Gyro State: "); Serial.print(gyroState[0], 3);
    Serial.print(", "); Serial.print(gyroState[1], 3);
    Serial.print(", "); Serial.print(gyroState[2], 3);
    Serial.print(", Acc State: "); Serial.print(accState[0], 3);
    Serial.print(", "); Serial.print(accState[1], 3);
    Serial.print(", "); Serial.println(accState[2], 3);
    lastPrintMillis = ms;
  }


  // TODO: implement sensor fusion of the two estimates - first just use gyro estimate so that I can see KF working
  // // sensor fusion of the estimates to get angle and angleRate
  // // may skip this step if using only gyro estimate due to computation time
  // // (especially when doing inverse calc of 3x3 matrix and checking for positive definiteness)
  // // P^-1 = (gyroP^-1 + accP^-1)
  // // xhat = P * (gyroP^-1 * gyroState + accP^-1 * accState)
  // float fusedAngle, fusedAngleRate, fusedGyroBias;
  // P[0][0] = 1.0f / gyroP[0] + 1.0f / accP[0];
  // P[1][1] = 1.0f / gyroP[4] + 1.0f / accP[4];
  // P[2][2] = 1.0f / gyroP[8] + 1.0f / accP[8]; 
  // float invDet = 1.0f / (P[0][0] * P[1][1] * P[2][2]); // diagonal matrix inverse
  // float P_inv[3][3];
  // for (int i = 0; i < 3; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     P_inv[i][j] = invDet * P[i][j];
  //   }
  // }
  // fusedAngle = P_inv[0][0] * (gyroState[0] / gyroP[0] +  accState[0] / accP[0]);
  // fusedAngleRate = P_inv[1][1] * (gyroState[1] / gyroP[4] +  accState[1] / accP[4]);
  // fusedGyroBias = P_inv[2][2] * (gyroState[2] / gyroP[8] +  accState[2] / accP[8]);

  // output to global variables
  angle = (int32_t)(gyroState[0] * 1000); // convert to millidegrees
  angleRate = (int32_t)(gyroState[1]); // degrees/s
  gYZero = (int32_t)(gyroState[2]);
}

void sensorsUpdate()
{
  // static variable to track time between updates
  static uint16_t lastMillis;
  uint16_t ms = millis();

  // check if it's time to update the sensors + control effort
  if ((uint16_t)(ms - lastMillis) < UPDATE_TIME_MS) { return; }
  // read sensors first
  imu.read(); // tells IMU to put new data in its registers

  // then apply control based on previous sensor readings and KF prediction
  // <--- apply control to TWIP here when ready, based on KF prediction of what state would be at this point --->
  // <--- reduces delay between filter estimation and control action --->

  // power on red LED if we were delayed
  sensorsUpdateDelayedStatus = ms - lastMillis > UPDATE_TIME_MS + 1;

  // finally, update the Kalman filter with the new sensor readings
  sensorsFilterUpdate(ms);

  lastMillis = ms;
    
}
