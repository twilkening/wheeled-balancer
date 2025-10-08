#include <Wire.h>
#include "Balance.h"
#include "LQR_PIDController.h"

int32_t gYZero;
int32_t angle; // millidegrees
int32_t angleRate; // degrees/s
int32_t distanceLeft;
int32_t speedLeft;
int32_t driveLeft;
int32_t distanceRight;
int32_t speedRight;
int32_t driveRight;
int16_t motorSpeed;
bool isBalancingStatus = false;
bool balanceUpdateDelayedStatus;

// LQR Controller - static allocation to avoid heap fragmentation
static LQR_PIDController angleLQR_PID(0, 0, 0, 0, 0, 0, 0, 0, 0); // Will be properly initialized in setup
bool lqr_pid_initialized = false;

bool isBalancing()
{
  return isBalancingStatus;
}

bool balanceUpdateDelayed()
{
  return balanceUpdateDelayedStatus;
}

void balanceSetup()
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

  // Calibrate the gyro.
  int32_t total = 0;
  for (int i = 0; i < CALIBRATION_ITERATIONS; i++)
  {
    imu.read();
    total += imu.g.y;
    delay(1);
  }

  gYZero = total / CALIBRATION_ITERATIONS;

  initializeLQR_PIDController();
}

// This function contains the core algorithm for balancing a
// Balboa 32U4 robot.
void balance()
{
  // Adjust toward angle=0 with timescale ~10s, to compensate for
  // gyro drift.  More advanced AHRS systems use the
  // accelerometer as a reference for finding the zero angle, but
  // this is a simpler technique: for a balancing robot, as long
  // as it is balancing, we know that the angle must be zero on
  // average, or we would fall over.
  // TODO - replace with accelerometer-based correction
  angle = angle * 999 / 1000;

  // Add safety check before using angleLQR_PID
  if (!lqr_pid_initialized) {
      return;
  }

  // Define reference values for LQR calculation
  int16_t ref_distance = 0, ref_speed = 0, ref_angle = 0, ref_angleRate = 0;
  int32_t distance = (distanceLeft + distanceRight) / 2;
  static int32_t speed_history[3] = {0, 0, 0};
  speed_history[2] = speed_history[1];
  speed_history[1] = speed_history[0];
  speed_history[0] = (speedLeft + speedRight) / 2;
  int32_t speed_filtered = (speed_history[0] + speed_history[1] + speed_history[2]) / 3;

  // Primary balance control using LQR
  int16_t angle_control = angleLQR_PID.calculate(ref_distance, ref_speed, ref_angle, ref_angleRate,
                                            &distance, &speed_filtered, &angle, &angleRate);

  // stall torque of motors is about 70 mN-m
  // PWM range is -400 to +400
  // => 400 / (70 mN-m * GEAR_RATIO) = 0.0449943 PWM cnt / mN-m
  // angle_control is in mN-m, so convert to PWM command 
  // (use int32_t to prevent overflow):
  motorSpeed = (int16_t)(((int32_t)angle_control * 45) / 1000);

  if (motorSpeed > MOTOR_SPEED_LIMIT)
  {
    motorSpeed = MOTOR_SPEED_LIMIT;
  }
  if (motorSpeed < -MOTOR_SPEED_LIMIT)
  {
    motorSpeed = -MOTOR_SPEED_LIMIT;
  }

  // Adjust for differences in the left and right distances; this
  // will prevent the robot from rotating as it rocks back and
  // forth due to differences in the motors, and it allows the
  // robot to perform controlled turns.
  int16_t distanceDiff = distanceLeft - distanceRight;

  motors.setSpeeds(
    motorSpeed + distanceDiff * DISTANCE_DIFF_RESPONSE / 100,
    motorSpeed - distanceDiff * DISTANCE_DIFF_RESPONSE / 100);

  // Serial.print("motorSpeed command: ");
  // Serial.println(motorSpeed);
}

void lyingDown()
{
  // Reset things so it doesn't go crazy.
  motorSpeed = 0;
  distanceLeft = 0;
  distanceRight = 0;
  motors.setSpeeds(0, 0);
  resetLQR_PIDControllers();

  if (angleRate > -2 && angleRate < 2)
  {
    // It's really calm, so use the accelerometer to measure the
    // robot's rest angle.  The atan2 function returns a result
    // in radians, so we multiply it by 180000/pi to convert it
    // to millidegrees.
    angle = atan2(imu.a.z, imu.a.x) * 57296;

    distanceLeft = 0;
    distanceRight = 0;
  }
}

void integrateGyro()
{
  // Apply sensitivity gain to gyro readings: 35 mdps/LSB  (for FS = +/-1000 dps)
  angleRate = (imu.g.y - gYZero) / 29; // units: degrees/s

  angle += angleRate * UPDATE_TIME_MS; // units: millidegrees
}

void integrateEncoders()
{
  static int16_t lastCountsLeft;
  int16_t countsLeft = encoders.getCountsLeft();
  speedLeft = (countsLeft - lastCountsLeft);
  distanceLeft += countsLeft - lastCountsLeft;
  lastCountsLeft = countsLeft;

  static int16_t lastCountsRight;
  int16_t countsRight = encoders.getCountsRight();
  speedRight = (countsRight - lastCountsRight);
  distanceRight += countsRight - lastCountsRight;
  lastCountsRight = countsRight;
}

void balanceDrive(int16_t leftSpeed, int16_t rightSpeed)
{
  driveLeft = leftSpeed;
  driveRight = rightSpeed;
}

// sets an offset to the distance and speed tracker to force the 
// balancer to think it has driven a certain amount, and thus return to "0"
void balanceDoDriveTicks()
{
  distanceLeft -= driveLeft;
  distanceRight -= driveRight;
  speedLeft -= driveLeft;
  speedRight -= driveRight;
}

void balanceResetEncoders()
{
  distanceLeft = 0;
  distanceRight = 0;
}

void balanceUpdateSensors()
{
  imu.read();
  integrateGyro();
  integrateEncoders();
}

void balanceUpdate()
{
  static uint16_t lastMillis;
  uint16_t ms = millis();
  static uint8_t count = 0;

  // Perform the balance updates at 100 Hz.
  if ((uint16_t)(ms - lastMillis) < UPDATE_TIME_MS) { return; }
  balanceUpdateDelayedStatus = ms - lastMillis > UPDATE_TIME_MS + 1;
  lastMillis = ms;

  balanceUpdateSensors();
  balanceDoDriveTicks();

  if (isBalancingStatus)
  {
    balance();

    // Stop trying to balance if we have been farther from
    // vertical than STOP_BALANCING_ANGLE for 5 counts.
    if (abs(angle) > STOP_BALANCING_ANGLE)
    {
      if (++count > 5)
      {
        isBalancingStatus = false;
        count = 0;
      }
    }
    else
    {
      count = 0;
    }
  }
  else
  {
    lyingDown();

    // Start trying to balance if we have been closer to
    // vertical than START_BALANCING_ANGLE for 5 counts.
    if (abs(angle) < START_BALANCING_ANGLE)
    {
      if (++count > 5)
      {
        isBalancingStatus = true;
        count = 0;
      }
    }
    else
    {
      count = 0;
    }
  }
}

// LQR Controller Functions
void initializeLQR_PIDController() {
  
  // Use int16_t with scaled integers for LQR gains (must fit in -32,767 to 32,767)
  const int16_t kx_scaled = 0;      // Position gain
  const int16_t kx_dot_scaled = (int16_t)(-1.4363 * 16.4913);   // Velocity gain [mN-m per cnts/10ms]
  const int16_t kth_scaled = (int16_t)(-58.7283 * 17.453);      // Angle gain [*uN-m* per millidegree]
  const int16_t kth_dot_scaled = (int16_t)(-1.6316 * 17.453);   // Angular velocity gain [mN-m per degree/s]
  const int16_t kp_scaled = static_cast<int16_t>(0.001 * 1200);     // Proportional gain [millidegrees per cnts]
  const int16_t ki_scaled = static_cast<int16_t>(0 * 1.2);      // Integral gain [millidegrees per (cnts*ms)] (set to 1 as minimum)
  const int32_t kd_scaled = static_cast<int32_t>(0.0001 * 1200000);   // Derivative gain [millidegrees per (cnts/ms)]
  const int32_t angle_integrator_limit = 1000; // millidegrees, <= 0 is no limit
  const int32_t angle_output_limit = 32000;   // <= 0 is no limit

  // Initialize the static LQR controller
  angleLQR_PID = LQR_PIDController(kx_scaled, kx_dot_scaled, kth_scaled, kth_dot_scaled,
                                    kp_scaled, ki_scaled, kd_scaled,
                                    angle_integrator_limit, angle_output_limit);
  lqr_pid_initialized = true;
}

void resetLQR_PIDControllers() {
  if (lqr_pid_initialized) angleLQR_PID.reset();
}
