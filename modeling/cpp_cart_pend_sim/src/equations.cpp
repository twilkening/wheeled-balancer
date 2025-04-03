#include <cmath> // For std::sqrt

// ODE: 
// dx0/dt = wn
// dx1/dt = -wc**2 * sin(theta) - 2*gamma*wn
// where wc = 1.0, gamma = 0.1, wn = 1.0
// and theta = x0, wn = x1
// This is a simple pendulum with damping
double f0(double t, double x0, double x1) {
    return x1; // return the derivatives as an array
}

double f1(double t, double x0, double x1) {
    double g = 9.81; // [m/s^2] gravitational acceleration
    double l = 1.0; // [m] length of the pendulum
    double wc = std::sqrt(g/l); // natural frequency
    double beta = 0.1; // damping ratio
    double m = 1.0; // [kg] mass of the pendulum bob
    double gamma = beta/(2 * m * l); // damping coefficient

    return  -wc * wc * sin(x0) - 2 * gamma * x1;
}