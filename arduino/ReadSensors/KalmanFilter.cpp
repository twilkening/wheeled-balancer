#include "KalmanFilter.h"
#include <MatrixMath.h>
#include <string.h>  // For memcpy

// Kalman Filter constructor implementation
KalmanFilter::KalmanFilter() {
    dt = 0.01; // Default time step
    R = 1.0;  // Default measurement noise variance
    // Initialize matrices and vectors to zero
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Q[i][j] = 0.0;
            P[i][j] = 0.0;
            F[i][j] = 0.0;
        }
        state[i] = 0.0;
        H[i] = 0.0;
    }
}

// // Kalman Filter destructor implementation
// KalmanFilter::~KalmanFilter() {
// }

void KalmanFilter::initialize(
    float dt,
    float x_init[3],
    float procNoise[3],
    float measNoise,
    float P_init[3][3],
    float H_init[3])
{
    this->dt = dt; // seconds
    // Set process noise variance matrix Q, diagonal elements
    for (int i = 0; i < 3; i++) {
        Q[i][i] = procNoise[i];
    }
    R = measNoise;

    // Initialize state and error covariance
    state[0] = x_init[0]; // angle
    state[1] = x_init[1]; // angleRate
    state[2] = x_init[2]; // gyroBias
    P[0][0] = P_init[0][0];
    P[0][1] = P_init[0][1];
    P[0][2] = P_init[0][2];
    P[1][0] = P_init[1][0];
    P[1][1] = P_init[1][1];
    P[1][2] = P_init[1][2];
    P[2][0] = P_init[2][0];
    P[2][1] = P_init[2][1];
    P[2][2] = P_init[2][2];
    // Initialize measurement matrix H
    H[0] = H_init[0];
    H[1] = H_init[1];
    H[2] = H_init[2];
}

void KalmanFilter::update(float measurement) {
    // Measurement update
    float e = measurement - (H[0] * state[0] + H[1] * state[1] + H[2] * state[2]); // Innovation
    float P_Ht[3]; // P * Ht
    // Re = H * P * Ht + R
    float Re; // Innovation covariance
    P_Ht[0] = P[0][0] * H[0] + P[0][1] * H[1] + P[0][2] * H[2];
    P_Ht[1] = P[1][0] * H[0] + P[1][1] * H[1] + P[1][2] * H[2];
    P_Ht[2] = P[2][0] * H[0] + P[2][1] * H[1] + P[2][2] * H[2];
    Re = H[0] * P_Ht[0] + H[1] * P_Ht[1] + H[2] * P_Ht[2] + R;

    float K_fi[3]; // Kalman Filter gain; K_fi = P * Ht * Re^-1
    K_fi[0] = P_Ht[0] / Re;
    K_fi[1] = P_Ht[1] / Re;
    K_fi[2] = P_Ht[2] / Re;

    // Update state; xhat_(i|i) = xhat_(i|i-1) + K_fi * e
    state[0] += K_fi[0] * e;
    state[1] += K_fi[1] * e;
    state[2] += K_fi[2] * e;

    // Update error covariance; P_(i|i) = (I - K_fi * H) * P_(i|i-1)
    float I_KH[3][3] = {
        {1 -K_fi[0] * H[0],    -K_fi[0] * H[1],   -K_fi[0] * H[2]},
        {  -K_fi[1] * H[0],  1 -K_fi[1] * H[1],   -K_fi[1] * H[2]},
        {  -K_fi[2] * H[0],    -K_fi[2] * H[1], 1 -K_fi[2] * H[2]}
    };
    float newP[3][3];
    Matrix.Multiply((float*)I_KH, (float*)P, 3, 3, 3, (float*)newP);
    // Copy newP back to P using memcpy (more efficient than nested loops)
    memcpy(P, newP, sizeof(float) * 9);
}

void KalmanFilter::predict() {
    // Predict state
    state[0] += state[1] * dt;
    state[1] += 0.0; // No acceleration term in this simple model
    state[2] += 0.0; // No bias update in this simple model

    // Predict error covariance
    // P_(i+1|i) = F * P_(i|i) * Ft + Q
    // where,
    // F = [1, dt, 0
    //      0, 1,  0
    //      0, 0,  1]  
    // and,
    // F*P*F' =
    // [P0_0 + P1_0*T + T*(P0_1 + P1_1*T), P0_1 + P1_1*T, P0_2 + P1_2*T]
    // [                    P1_0 + P1_1*T,          P1_1,          P1_2]
    // [                    P2_0 + P2_1*T,          P2_1,          P2_2]
    P[0][0] = P[0][0] + dt * P[1][0] + dt * ( P[0][1] + dt * P[1][1] ) + Q[0][0];
    P[0][1] = P[0][1] + dt * P[1][1];
    P[0][2] = P[0][2] + dt * P[1][2];
    P[1][0] = P[1][0] + dt * P[1][1];
    P[1][1] = P[1][1] + Q[1][1];
    P[1][2] = P[1][2];
    P[2][0] = P[2][0] + dt * P[2][1];
    P[2][1] = P[2][1];
    P[2][2] = P[2][2] + Q[2][2];
}

float* KalmanFilter::getState() {
    return state;
}

float* KalmanFilter::getCovariance() {
    return &P[0][0];
}

