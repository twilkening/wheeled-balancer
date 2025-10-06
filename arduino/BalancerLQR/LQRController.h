#pragma once

#include <stdint.h>

/**
 * LQR Controller for Arduino-based balancing robot
 * Adapted from MuJoCo implementation for use with Balboa 32U4
 */
class LQRController {
private:
    // LQR gains
    // could make these const to free up memory if desired
    int32_t kx, kx_dot, kth, kth_dot; // State feedback gains
    
    // Limits and settings
    uint32_t t_prior;      // time of prior update (ms)
    int16_t output_limit; // output saturation limit
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
    LQRController(int16_t kx_gain, int16_t kx_dot_gain, int16_t kth_gain, int16_t kth_dot_gain,
                    int16_t output_limit = 0);
    
    /**
     * Calculate LQR control output
     * @param reference0-3 Desired setpoint values for each state
     * @param measurement0-3 Current measured values for each state
     * @param dt_ms Time step in milliseconds (0 = auto-calculate from millis())
     * @return Control output
     */
    int16_t calculate(
        int16_t reference0,
        int16_t reference1,
        int16_t reference2,
        int16_t reference3,
        int32_t* measurement0,
        int32_t* measurement1,
        int32_t* measurement2,
        int32_t* measurement3,
        uint32_t dt_ms = 0
    );
    
    /**
     * Reset the LQR controller (clears integral and derivative history)
     */
    void reset();
    
};
