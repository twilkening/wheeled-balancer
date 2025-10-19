#pragma once

#include <stdint.h>

class KalmanFilter {
    public:
        KalmanFilter();
        // ~KalmanFilter(); // Destructor (if needed)

        void initialize(float dt, float x_init[3], float procNoise[3], float measNoise, float P_init[3][3], float H_init[3]);
        void update(float measurement);
        void predict();
        float* getState();
        float* getCovariance();
    private:
        float dt; // nominal time step
        float Q[3][3]; // Process noise variance
        float R; // Measurement noise variance
        float state[3]; // [angle, angleRate, gyroBias]
        float P[3][3]; // Error covariance matrix
        float F[3][3]; // State transition matrix
        float H[3]; // Measurement matrix
};
