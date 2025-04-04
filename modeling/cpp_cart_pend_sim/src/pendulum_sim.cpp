#include "rk4.cpp"
#include "pendulumODE.cpp"
#include <string>
#include <iostream>

int main()
{
    // Read simulation parameters from a file
    std::ifstream param_file("siminit.txt");
    if (!param_file.is_open()) {
        std::cerr << "Error opening siminit.txt" << std::endl;
        return 1;
    }
    double t0, t_end, h, theta0, omega0;

    // Use the helper function to read each parameter, line by line
    if (!readParameter(param_file, t0, "t0")) return 1;
    if (!readParameter(param_file, t_end, "t_end")) return 1;
    if (!readParameter(param_file, h, "step size h")) return 1;
    if (!readParameter(param_file, theta0, "theta")) return 1;
    if (!readParameter(param_file, omega0, "omega")) return 1;

    PendulumODE ode("pendparm.txt");
    ode.printParameters();
    auto f0 = [&ode](double t, double x0, double x1) { return ode.f0(t, x0, x1); };
    auto f1 = [&ode](double t, double x0, double x1) { return ode.f1(t, x0, x1); };
    
    const double pi = 3.141592653589793;
    double xinit[2] = {theta0, omega0}; // initial conditions [theta, wn]
    std::string fileout = "simout.txt"; // output file name

    rk4 solver(f0, f1, xinit, t0, t_end, h); // create RK4 solver instance
    solver.run(fileout); // run the solver
    
    return 0;
}