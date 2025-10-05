#pragma once

#include <stdint.h>

/**
 * PID Controller for Arduino-based balancing robot
 * Adapted from MuJoCo implementation for use with Balboa 32U4
 */
class PIDController {
private:
    // PID gains
    int32_t kp, ki, kd;
    
    // State variables
    int32_t e_prior;      // previous error
    uint32_t t_prior;    // previous time (milliseconds)
    int32_t integral;     // integral accumulator
    
    // Limits and settings
    int16_t int_limit;    // integral windup limit
    int16_t output_limit; // output saturation limit
    bool first_run;      // flag for first execution
    
public:
    /**
     * Constructor - initializes PID controller
     * @param kp_gain Proportional gain
     * @param ki_gain Integral gain  
     * @param kd_gain Derivative gain
     * @param integral_limit Anti-windup limit for integral term (0 = no limit)
     * @param output_limit Output saturation limit (0 = no limit)
     */
    PIDController(int32_t kp_gain, int32_t ki_gain, int32_t kd_gain, 
                    int16_t integral_limit = 0, int16_t output_limit = 0);
    
    /**
     * Calculate PID control output
     * @param reference Desired setpoint value
     * @param measurement Current measured value
     * @param dt_ms Time step in milliseconds (0 = auto-calculate from millis())
     * @return Control output
     */
    double calculate(int32_t reference, int32_t measurement, uint32_t dt_ms = 0);
    
    /**
     * Reset the PID controller (clears integral and derivative history)
     */
    void reset();
    
    /**
     * Update PID gains during runtime
     */
    void setGains(int32_t kp_new, int32_t ki_new, int32_t kd_new);
    
    /**
     * Get current PID gains
     */
    void getGains(int32_t* kp_out, int32_t* ki_out, int32_t* kd_out);
    
    /**
     * Set integral limit (anti-windup)
     */
    void setIntegralLimit(int16_t limit);
    
    /**
     * Set output saturation limit
     */
    void setOutputLimit(int16_t limit);
    
    /**
     * Get current integral value (for debugging/monitoring)
     */
    int32_t getIntegral() const { return integral; }
    
    /**
     * Get last error value (for debugging/monitoring) 
     */
    int32_t getLastError() const { return e_prior; }
};
