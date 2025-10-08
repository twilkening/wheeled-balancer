#include "LQR_PIDController.h"
#include <Arduino.h>
#include <math.h>
#include "MemoryFree.h"

// Base LQR Controller Implementation
LQR_PIDController::LQR_PIDController(int32_t kx_gain, int32_t kx_dot_gain, int32_t kth_gain, int32_t kth_dot_gain,
                                    int32_t kp_gain, int32_t ki_gain, int32_t kd_gain,
                                    int32_t integrator_limit_in, int32_t output_limit_in) 
    : kx(kx_gain), kx_dot(kx_dot_gain), kth(kth_gain), kth_dot(kth_dot_gain),
        kp(kp_gain), ki(ki_gain), kd(kd_gain), 
        integrator_limit(integrator_limit_in), output_limit(output_limit_in),
        last_error(0), integrator(0), t_prior(0), first_run(true) {
}

int32_t LQR_PIDController::calculate(
    int32_t reference0,
    int32_t reference1,
    int32_t reference2,
    int32_t reference3,
    int32_t measurement0,
    int32_t measurement1,
    int32_t measurement2,
    int32_t measurement3,
    uint32_t dt_ms
) {
    uint32_t current_time = millis();
    
    // Calculate time step
    int32_t dt;
    if (dt_ms > 0) {
        dt = dt_ms; // Use provided dt in milliseconds
    } else if (first_run) {
        dt = 10; // Default 10ms for first run
        first_run = false;
        t_prior = current_time;
        return 0; // No control action on first run
    } else {
        dt = (current_time - t_prior);
    }
    
    // Guard against zero or negative dt
    // especially important if millis() overflows
    if (dt <= 0) dt = 10; // Minimum 10ms

    int32_t offset = PID(reference0, measurement0, dt); // effective unit: millidegrees

    // Calculate control output
    // offset angle measurement by 3352 millidegrees to account for IMU mounting angle w.r.t. CoM location
    // **0.058508 radians * 180 deg / pi radians * 1000 millideg/deg = 3352 millidegrees
    int32_t output =  ( - kx        * (measurement0 - reference0)         // mNm
                        - kx_dot    * (measurement1 - reference1)         // mNm
                        - kth       * (measurement2 - 3352 - reference2 - offset) / 1000  // convert uNm to mNm
                        - kth_dot   * (measurement3 - reference3) )       // mNm
                        / 2; // divide by 2 bc two motors

    // Apply output saturation if enabled
    if (output_limit > 0) {
        if (output > output_limit) {
            output = output_limit;
        } else if (output < -output_limit) {
            output = -output_limit;
        }
    }

    // Update history
    t_prior = current_time;

    return output;
}

int32_t LQR_PIDController::PID(int32_t setpoint, int32_t measurement, uint32_t dt) {
    int32_t error = setpoint - measurement; // encoder counts

    // Proportional term [millidegrees/cnts]
    int32_t P_out = (kp * error); // millidegrees

    // Integral term [millidegree/(cnt*ms)]
    integrator += error * dt; // error in cnts * dt in ms = cnts*ms
    // Anti-windup: clamp integrator
    if (integrator_limit > 0) {
        if (integrator > integrator_limit) {
            integrator = integrator_limit;
        } else if (integrator < -integrator_limit) {
            integrator = -integrator_limit;
        }
    }
    int32_t I_out = (ki * integrator); // millidegrees

    // Derivative term [millidegrees/(cnt/ms)]
    int32_t derivative = dt > 0 ? (error - last_error) / dt : 0;  // cnts/ms
    int32_t D_out = (kd * derivative); // millidegrees

    // Save error for next derivative calculation
    last_error = error;

    // static uint32_t last_memory_print = 0;
    // if (millis() - last_memory_print > 1000) {  // Print every second
    //     Serial.print("Free memory: ");
    //     Serial.println(freeMemory());
    //     last_memory_print = millis();
    // }

    // Total output
    return P_out + I_out + D_out;
}

void LQR_PIDController::reset() {
    t_prior = millis();
    first_run = true;
    integrator = 0;
    last_error = 0;
}
