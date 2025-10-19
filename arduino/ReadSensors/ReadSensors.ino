asm(".global _printf_float");
#include <Balboa32U4.h>
#include <Wire.h>
#include <LSM6.h>
#include "Sensors.h"

LSM6 imu;

void setup()
{
  // Uncomment these lines if your motors are reversed.
  // motors.flipLeftMotor(true);
  // motors.flipRightMotor(true);

  ledYellow(0);
  ledRed(1);
  sensorsSetup();
  ledRed(0);

  Serial.begin(57600);
}


void loop()
{
  sensorsUpdate();

  // Illuminate the red LED if the last full update was too slow.
  ledRed(sensorsUpdateDelayed());

}
