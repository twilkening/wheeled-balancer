#pragma once

#include <stdint.h>

/**
 * LQR Controller for Arduino-based balancing robot
 * Adapted from MuJoCo implementation for use with Balboa 32U4
 */
class LQR_PIDController {
private:
    // LQR gains - using int32_t for better performance (no casting needed)
    int32_t kx, kx_dot, kth, kth_dot; // State feedback gains
    // PID gains
    int32_t kp, ki, kd;

    // Limits and settings
    int32_t last_error; // Last error for derivative calculation
    int32_t integrator; // Integral term accumulator
    int32_t integrator_limit; // Integral windup limit
    uint32_t t_prior;      // time of prior update (ms)
    int32_t output_limit; // output saturation limit
    bool first_run;      // flag for first execution
    
public:
    /**
     * Constructor - initializes LQR controller
     * @param kx_gain Position gain
     * @param kx_dot_gain Velocity gain
     * @param kth_gain Angular gain
     * @param kth_dot_gain Angular velocity gain
     * @param output_limit Output saturation limit (0 = no limit)
     */
    LQR_PIDController(int32_t kx_gain, int32_t kx_dot_gain, int32_t kth_gain, int32_t kth_dot_gain,
                    int32_t kp_gain, int32_t ki_gain, int32_t kd_gain,
                    int32_t integrator_limit_in = 0, int32_t output_limit_in = 0);
    
    /**
     * Calculate LQR control output
     * @param reference0-3 Desired setpoint values for each state
     * @param measurement0-3 Current measured values for each state
     * @param dt_ms Time step in milliseconds (0 = auto-calculate from millis())
     * @return Control output
     */
    int32_t calculate(
        int32_t reference0,
        int32_t reference1,
        int32_t reference2,
        int32_t reference3,
        int32_t measurement0,
        int32_t measurement1,
        int32_t measurement2,
        int32_t measurement3,
        uint32_t dt_ms = 0
    );
    
    /**
     * Compute PID control output
     * @param setpoint Desired setpoint
     * @param measurement Current measured value
     * @param Kp Proportional gain
     * @param Ki Integral gain
     * @param Kd Derivative gain
     * @param integrator_limit Integral windup limit
     * @return PID control output
     */
    int32_t PID(int32_t setpoint, int32_t measurement, uint32_t dt);

    /**
     * Reset the LQR controller (clears integral and derivative history)
     */
    void reset();
    
};
