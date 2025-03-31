#define _USE_MATH_DEFINES
#include <cmath>
#include "equations.h"

const double g = 9.81;   // Gravity (m/s²)
const double m = 1.0;    // Pendulum mass (kg)
const double M = 5.0;    // Cart mass (kg)
const double L = 2.0;    // Pendulum length (m)
const double b = 0.1;    // Damping coefficient
const double u = 0.0;    // Control force (for now, zero)

State cart_pendulum_dynamics(double t, const State& y) {
    double x = y[0];      // Cart position
    double v = y[1];      // Cart velocity
    double theta = y[2];  // Pendulum angle
    double omega = y[3];  // Pendulum angular velocity

    double sinTheta = sin(theta);
    double cosTheta = cos(theta);
    double denom = (M + m * sinTheta * sinTheta);

    double dxdt = v;
    double dvdt = (u + m * sinTheta * (L * omega * omega + g * cosTheta) - b * v) / denom;
    double dthetadt = omega;
    double domegadt = (-u * cosTheta - m * L * omega * omega * sinTheta * cosTheta - (M + m) * g * sinTheta) / (L * denom);

    return {dxdt, dvdt, dthetadt, domegadt};
}
