#include "PIDController.h"
#include <Arduino.h>
#include <math.h>

// Base PID Controller Implementation
PIDController::PIDController(int32_t kp_gain, int32_t ki_gain, int32_t kd_gain, 
                            int16_t integral_limit, int16_t output_limit) 
    : kp(kp_gain), ki(ki_gain), kd(kd_gain),
        e_prior(0), t_prior(0), integral(0),
        int_limit(integral_limit), output_limit(output_limit), first_run(true) {
}

int16_t PIDController::calculate(int32_t reference, int32_t measurement, uint32_t dt_ms) {
    uint32_t current_time = millis();
    
    // Calculate time step
    int32_t dt;
    if (dt_ms > 0) {
        dt = dt_ms; // Use provided dt in milliseconds
    } else if (first_run) {
        dt = 10; // Default 10ms for first run
        first_run = false;
        return 0; // No control action on first run
    } else {
        dt = (current_time - t_prior);
    }
    
    // Guard against zero or negative dt
    // especially important if millis() overflows
    if (dt <= 0) dt = 10; // Minimum 10ms

    // Calculate error
    int32_t error = reference - measurement; // millidegrees
    
    // Calculate derivative of error (with safe guard)
    int32_t derivative = (dt > 0) ? (error - e_prior) / dt : 0;
    
    // Update integral with anti-windup clamping
    integral += error * dt;
    if (int_limit > 0) {
        if (integral > int_limit) {
            integral = int_limit;
        } else if (integral < -int_limit) {
            integral = -int_limit;
        }
    }
    
    // Calculate control output
    int32_t output = kp * error + ki * integral + kd * derivative;
    
    // Apply output saturation if enabled
    if (output_limit > 0) {
        if (output > output_limit) {
            output = output_limit;
        } else if (output < -output_limit) {
            output = -output_limit;
        }
    }

    // Update history
    e_prior = error;
    t_prior = current_time;

    return (int16_t)output;
}

void PIDController::reset() {
    e_prior = 0;
    integral = 0;
    t_prior = millis();
    first_run = true;
}

void PIDController::setGains(int32_t kp_new, int32_t ki_new, int32_t kd_new) {
    kp = kp_new;
    ki = ki_new;
    kd = kd_new;
}

void PIDController::getGains(int32_t* kp_out, int32_t* ki_out, int32_t* kd_out) {
    if (kp_out) *kp_out = kp;
    if (ki_out) *ki_out = ki;
    if (kd_out) *kd_out = kd;
}

void PIDController::setIntegralLimit(int16_t limit) {
    int_limit = limit;
}

void PIDController::setOutputLimit(int16_t limit) {
    output_limit = limit;
}
