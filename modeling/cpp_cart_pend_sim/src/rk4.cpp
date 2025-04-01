#include <cmath>
#include <vector>
#include <array>
#include <iostream>

// RK4 class for solving ODEs using the Runge-Kutta 4th order method
// This class is designed to solve a system of two first-order ODEs
class rk4
{
    double xn[2]; // current states
    double k1[2], k2[2], k3[2], k4[2]; // RK4 coefficients
    const double h; // step size
    double t; // current time
    double t_end; // end time
    double (*f0)(double, double, double); // function pointer for the ODE0
    double (*f1)(double, double, double); // function pointer for the ODE1
public:
    rk4(double (*f0)(double, double, double), double (*f1)(double, double, double),
        double xinit[2], double t0, double t_end, double h)
        : f0(f0), f1(f1), t(t0), t_end(t_end), h(h)
    {
        // Initialize the states
        xn[0] = xinit[0];
        xn[1] = xinit[1];
        std::cout << "h: " << h << std::endl;
        std::cout << "t_end: " << t_end << std::endl;
        std::cout << "t: " << t << ", x0: " << xn[0] << ", x1: " << xn[1] << std::endl;
        std::cout << "f0: " << f0(t, xn[0], xn[1]) << std::endl;
        std::cout << "f1: " << f1(t, xn[0], xn[1]) << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "----------------------------------------" << std::endl;
    }
    void step(); // perform a single RK4 step
    bool isComplete() const; // check if the simulation is complete
    void run(); // run the RK4 solver
    std::array<double, 2> getState() const; // get the current state
};

void rk4::step()
{
    k1[0] = h * f0(t, xn[0], xn[1]);
    k1[1] = h * f1(t, xn[0], xn[1]);

    k2[0] = h * f0(t + h / 2, xn[0] + k1[0] / 2, xn[1] + k1[1] / 2);
    k2[1] = h * f1(t + h / 2, xn[0] + k1[0] / 2, xn[1] + k1[1] / 2);

    k3[0] = h * f0(t + h / 2, xn[0] + k2[0] / 2, xn[1] + k2[1] / 2);
    k3[1] = h * f1(t + h / 2, xn[0] + k2[0] / 2, xn[1] + k2[1] / 2);

    k4[0] = h * f0(t + h, xn[0] + k3[0], xn[1] + k3[1]);
    k4[1] = h * f1(t + h, xn[0] + k3[0], xn[1] + k3[1]);

    xn[0] += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
    xn[1] += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;

    t += h;
}

bool rk4::isComplete() const
{
    return t >= t_end;
}

void rk4::run()
{
    while (t < t_end) {
        step();
        std::cout << "t: " << t << ", x0: " << xn[0] << ", x1: " << xn[1] << std::endl;
        // You can also store the results in a vector or array if needed
    }
}

std::array<double, 2> rk4::getState() const
{
    return {xn[0], xn[1]};
}
