# Function Call Map for Two-Wheeled Inverted Pendulum Arduino Code

This document provides a comprehensive map of function interactions in the balancing robot codebase.

## Main Program Flow (`Balancer.ino`)

### Setup Function
```
setup()
├── balanceSetup() → Initialize IMU and calibrate gyro
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
balance() (core PID-like controller)
├── Drift compensation: angle = angle * 999/1000
├── Calculate risingAngleOffset = angleRate * ANGLE_RATE_RATIO + angle
├── PID calculation:
│   motorSpeed += (ANGLE_RESPONSE * risingAngleOffset +
│                  DISTANCE_RESPONSE * (distanceLeft + distanceRight) +
│                  SPEED_RESPONSE * (speedLeft + speedRight)) / 100 / GEAR_RATIO
├── Apply motor speed limits
├── Calculate steering correction: distanceDiff = distanceLeft - distanceRight
└── motors.setSpeeds() → Apply final motor commands
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

### Output
- **Motors**: `motors.setSpeeds(left, right)`
- **LEDs**: Status indicators
- **Buzzer**: Audio feedback

## Control Parameters

### PID-like Control Constants
- **`ANGLE_RESPONSE`** (11): Primary balance correction
- **`DISTANCE_RESPONSE`** (73): Position holding
- **`SPEED_RESPONSE`** (3300): Oscillation damping  
- **`DISTANCE_DIFF_RESPONSE`** (-50): Steering correction

### Timing and Thresholds
- **`UPDATE_TIME_MS`** (10): 100 Hz control loop
- **`START_BALANCING_ANGLE`** (45000): Begin balancing threshold
- **`STOP_BALANCING_ANGLE`** (70000): Stop balancing threshold
- **`MOTOR_SPEED_LIMIT`** (300): Maximum motor speed
- **`GEAR_RATIO`** (127): Mechanical gear ratio

## Architecture Summary

This implementation uses a classic inverted pendulum control architecture with cascaded loops:

1. **Inner Loop**: Stabilizes the tilt angle using gyroscope feedback
2. **Outer Loops**: Control position and handle user commands for autonomous behaviors
3. **State Machine**: Manages transitions between balancing and non-balancing states
4. **User Interface**: Provides button controls and audio/visual feedback

The control algorithm runs at 100 Hz and combines angle, position, and velocity feedback to maintain balance while allowing controlled movement and user interaction.