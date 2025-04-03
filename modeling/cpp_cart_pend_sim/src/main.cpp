#include "rk4.cpp"
#include "equations.cpp"
#include <string>

int main()
{
    // Define the ODE parameters
    double t0 = 0.0; // initial time
    double t_end = 10.0; // end time
    double h = 0.001; // step size
    const double pi = 3.141592653589793;
    double xinit[2] = {pi/2, 0}; // initial conditions [theta, wn]
    std::string fout = "simout.txt"; // output file name

    rk4 solver(f0, f1, xinit, t0, t_end, h); // create RK4 solver instance
    solver.run(fout); // run the solver
    
    return 0;
}