# Function Call Map for Two-Wheeled Inverted Pendulum Arduino Code

This document provides a comprehensive map of function interactions in the balancing robot codebase with integrated PID controller implementation.

## Main Program Flow (`BalancerPID.ino`)

### Setup Function
```
setup()
├── balanceSetup() → Initialize IMU, calibrate gyro, and setup PID controllers
└── Set LED states
```

### Main Loop
```
loop() (main control loop)
├── balanceUpdate() → Core balancing algorithm
├── isBalancing() → Check if robot is actively balancing
│   ├── TRUE (balancing):
│   │   ├── playSong() (optional) → Play melody
│   │   └── driveAround() (optional) → Autonomous driving
│   │       └── balanceDrive() → Set target speeds
│   └── FALSE (not balancing):
│       ├── buzzer.stopPlaying()
│       ├── balanceDrive(0, 0) → Reset speeds
│       └── Button handling:
│           ├── buttonA → standUp()
│           ├── buttonB → standUp() + enable driving
│           └── buttonC → standUp() + enable driving + song
├── balanceUpdateDelayed() → Check timing performance
└── LED feedback based on fallingAngleOffset
```

## Core Balancing Functions (`Balance.cpp`)

### Primary Control Loop
```
balanceUpdate() (called every 10ms from main loop)
├── Check timing (100 Hz update rate)
├── balanceUpdateSensors()
│   ├── imu.read() → Get IMU data
│   ├── integrateGyro() → Process gyro data
│   │   └── Updates: angleRate, angle
│   └── integrateEncoders() → Process encoder data
│       └── Updates: speedLeft, speedRight, distanceLeft, distanceRight
├── balanceDoDriveTicks() → Apply driving commands
└── State Management:
    ├── IF isBalancing == TRUE:
    │   ├── balance() → Execute PID control algorithm
    │   └── Check if angle > STOP_BALANCING_ANGLE (fall detection)
    └── IF isBalancing == FALSE:
        ├── lyingDown() → Reset state when robot is down
        └── Check if angle < START_BALANCING_ANGLE (pickup detection)
```

### Control Algorithm
```
balance() (core PID controller implementation)
├── Drift compensation: angle = angle * 999/1000
├── Safety check: Verify anglePID is initialized
├── PID angle control:
│   ├── anglePID->calculate(0, angle, UPDATE_TIME_MS) → Calculate control output
│   └── Scale and apply to motorSpeed: -= (angle_control/GEAR_RATIO) * 100/572
├── Apply motor speed limits
├── Calculate steering correction: distanceDiff = distanceLeft - distanceRight
└── motors.setSpeeds() → Apply final motor commands with steering adjustment
```

## Sensor Processing Functions

```
integrateGyro()
├── Apply sensitivity: angleRate = (imu.g.y - gYZero) / 29
└── Integrate: angle += angleRate * UPDATE_TIME_MS

integrateEncoders()
├── Calculate speed: speedLeft/Right = (currentCounts - lastCounts)
├── Integrate distance: distanceLeft/Right += countsDifference
└── Store current counts for next iteration

lyingDown() (when robot is horizontal)
├── Reset motor speeds and distances
├── IF angleRate is very small:
│   ├── Use accelerometer: angle = atan2(imu.a.z, imu.a.x) * 57296
│   └── Reset distance measurements
```

## PID Controller Functions (`PIDController.cpp`)

### PID Controller Class
```
PIDController::PIDController() (constructor)
├── Initialize gains: kp, ki, kd
├── Initialize state: e_prior=0, integral=0, first_run=true
├── Set limits: integral_limit, output_limit
└── Store current time: t_prior

PIDController::calculate() (main PID computation)
├── Calculate time step: dt = current_time - t_prior (or use provided dt_ms)
├── Calculate error: error = reference - measurement
├── Calculate derivative: derivative = (error - e_prior) / dt
├── Update integral with anti-windup:
│   ├── integral += error * dt
│   └── Clamp to ±integral_limit if enabled
├── Calculate output: output = kp*error + ki*integral + kd*derivative
├── Apply output saturation if enabled: clamp to ±output_limit
├── Update history: e_prior = error, t_prior = current_time
└── Return control output

PIDController::reset()
├── Clear error history: e_prior = 0
├── Clear integral: integral = 0
├── Reset timing: t_prior = millis(), first_run = true
```

### PID Integration Functions
```
initializeAnglePID() (called from balanceSetup())
├── Define scaled PID gains:
│   ├── kp_scaled = 15 (mN-m per millidegree)
│   ├── ki_scaled = 2 (mN-m per millidegree-ms)
│   └── kd_scaled = 2500 (mN-m per millidegree/ms)
├── Set limits:
│   ├── angle_integral_limit = 6000
│   └── angle_output_limit = 20000
└── Create PIDController instance: anglePID = new PIDController(...)

resetPIDControllers()
└── anglePID->reset() → Clear PID state when robot falls
```

## User Interface Functions

```
standUp() (kick-up sequence)
├── motors.setSpeeds(0, 0) → Stop motors
├── buzzer.play() → Audio feedback
├── LED indicators
├── motors.setSpeeds(-MOTOR_SPEED_LIMIT, -MOTOR_SPEED_LIMIT) → Reverse thrust
├── delay(400)
├── motors.setSpeeds(150, 150) → Forward thrust
├── Loop until angle < 60000 (detection of successful stand-up)
└── balanceResetEncoders() → Clear encoder history

driveAround() (autonomous driving pattern)
├── Calculate time-based speeds (8-second cycle)
├── Set leftSpeed and rightSpeed based on time
└── balanceDrive(leftSpeed, rightSpeed)

playSong()
├── Check if buzzer is already playing
└── buzzer.playFromProgramSpace(song)
```

## Data Flow

### Input Sources
- **IMU (LSM6)**: `imu.g.y` (gyro), `imu.a.x`, `imu.a.z` (accelerometer)
- **Encoders**: `encoders.getCountsLeft()`, `encoders.getCountsRight()`
- **Buttons**: `buttonA`, `buttonB`, `buttonC`

### Key Variables
- **`angle`**: Current tilt angle (millidegrees)
- **`angleRate`**: Angular velocity (degrees/s)
- **`motorSpeed`**: Base motor speed command
- **`distanceLeft/Right`**: Accumulated wheel distances
- **`speedLeft/Right`**: Current wheel speeds
- **`driveLeft/Right`**: User-commanded drive speeds
- **`anglePID`**: Pointer to angle PID controller instance
- **`gYZero`**: Gyro zero offset from calibration

### PID Controller Internal State
- **`e_prior`**: Previous error value for derivative calculation
- **`integral`**: Accumulated integral term with anti-windup
- **`t_prior`**: Previous time stamp for dt calculation
- **`first_run`**: Flag to handle first execution cycle

### Output
- **Motors**: `motors.setSpeeds(left, right)`
- **LEDs**: Status indicators
- **Buzzer**: Audio feedback

## Control Parameters

### PID Controller Gains (Angle Control)
- **`kp_scaled`** (15): Proportional gain - mN-m output per millidegree angle error
- **`ki_scaled`** (2): Integral gain - mN-m output per millidegree-ms accumulated error  
- **`kd_scaled`** (2500): Derivative gain - mN-m output per millidegree/ms angle rate error
- **`angle_integral_limit`** (6000): Anti-windup limit for integral accumulator
- **`angle_output_limit`** (20000): Saturation limit for PID output

### Timing and Thresholds
- **`UPDATE_TIME_MS`** (10): 100 Hz control loop
- **`START_BALANCING_ANGLE`** (45000): Begin balancing threshold
- **`STOP_BALANCING_ANGLE`** (70000): Stop balancing threshold
- **`MOTOR_SPEED_LIMIT`** (300): Maximum motor speed
- **`GEAR_RATIO`** (127): Mechanical gear ratio

## Architecture Summary

This implementation uses a classic inverted pendulum control architecture with a formal PID controller:

1. **Inner Loop**: PID-controlled angle stabilization using gyroscope feedback
   - Proportional term: Responds to current angle error
   - Integral term: Eliminates steady-state error with anti-windup protection
   - Derivative term: Provides damping against oscillations
2. **State Machine**: Manages transitions between balancing and non-balancing states
3. **User Interface**: Provides button controls and audio/visual feedback
4. **PID Controller**: Dedicated class with proper integral windup protection and output saturation

The control algorithm runs at 100 Hz with a dedicated PID controller that provides more predictable and tunable performance compared to the previous ad-hoc control law.

## File Structure

### Core Files
- **`BalancerPID.ino`**: Main Arduino sketch with setup() and loop() functions
- **`Balance.h`**: Header file with constants, function declarations, and PID controller interface
- **`Balance.cpp`**: Core balancing algorithm implementation with PID integration

### PID Controller Module  
- **`PIDController.h`**: PID controller class header with method declarations
- **`PIDController.cpp`**: PID controller implementation with anti-windup and saturation
- **`Function_Call_Map.md`**: This documentation file
