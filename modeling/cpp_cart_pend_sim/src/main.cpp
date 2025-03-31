/*In order to actually simulate the file... 
it's really not that hard to implement the RK4 is it?

equation... sin(x) 

let's see if we can do this!*/

#include "rk4.cpp"
#include <cmath>
#include <array>

std::array<double, 2> f(double t, const double x[2]) {
    // ODE: 
    // dx1/dt = wn
    // dx2/dt = -wc**2 * sin(theta) - 2*gamma*wn
    // where wc = 1.0, gamma = 0.1, wn = 1.0
    // and theta = x1, wn = x2
    // This is a simple pendulum with damping

    double wc = 1.0; // natural frequency
    double gamma = 0.1; // damping coefficient
    double f1 = x[1]; 
    double f2 = -wc * wc * sin(x[0]) - 2 * gamma * x[1];
    return {f1, f2}; // return the derivatives as an array
}

int main()
{
    // Define the ODE function
    // dx/dt = f(t, x)  where f is a function of t and x
    // For example, let's use the ODE dx/dt = sin(x)

    double t0 = 0.0; // initial time
    double t_end = 10.0; // end time
    double h = 0.01; // step size

    rk4 solver(f(), t0, t_end, h); // create RK4 solver instance
    solver.run(); // run the solver
    
    return 0;
}