#include "LQRController.h"
#include <Arduino.h>
#include <math.h>

// Base LQR Controller Implementation
LQRController::LQRController(int16_t kx_gain, int16_t kx_dot_gain, int16_t kth_gain, int16_t kth_dot_gain,
                            int16_t output_limit) 
    : kx(kx_gain), kx_dot(kx_dot_gain), kth(kth_gain), kth_dot(kth_dot_gain),
        t_prior(0), output_limit(output_limit), first_run(true) {
}

int16_t LQRController::calculate(
    int16_t reference0,
    int16_t reference1,
    int16_t reference2,
    int16_t reference3,
    int32_t* measurement0,
    int32_t* measurement1,
    int32_t* measurement2,
    int32_t* measurement3,
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
        return 0; // No control action on first run
    } else {
        dt = (current_time - t_prior);
    }
    
    // Guard against zero or negative dt
    // especially important if millis() overflows
    if (dt <= 0) dt = 10; // Minimum 10ms

    // Calculate control output
    int32_t output =  ( - (int32_t)kx        * (*measurement0 - reference0)         // mNm
                        - (int32_t)kx_dot    * (*measurement1 - reference1)         // mNm
                        - (int32_t)kth       * (*measurement2 - reference2) / 1000  // convert uNm to mNm
                        - (int32_t)kth_dot   * (*measurement3 - reference3) )       // mNm
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

    return (int16_t)output;
}

void LQRController::reset() {
    t_prior = millis();
    first_run = true;
}
