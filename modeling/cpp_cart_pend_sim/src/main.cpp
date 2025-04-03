#include "rk4.cpp"
#include "equations.cpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

// Helper function to read a parameter from a line
bool readParameter(std::ifstream& file, double& param, const std::string paramName) {
    std::string line;
    if (std::getline(file, line)) {
        std::stringstream ss(line.substr(0, line.find("//")));
        ss >> param;
        if (ss.fail()) {
            std::cerr << "Error reading " << paramName << std::endl;
            return false;
        }
        return true;
    } else {
        std::cerr << "Error reading " << paramName << std::endl;
        return false;
    }
}

int main()
{
    // Read ODE parameters from a file
    std::ifstream param_file("siminit.txt");
    if (!param_file.is_open()) {
        std::cerr << "Error opening siminit.txt" << std::endl;
        return 1;
    }
    double t0, t_end, h, theta, omega;

    // Use the helper function to read each parameter
    if (!readParameter(param_file, t0, "t0")) return 1;
    if (!readParameter(param_file, t_end, "t_end")) return 1;
    if (!readParameter(param_file, h, "step size h")) return 1;
    if (!readParameter(param_file, theta, "theta")) return 1;
    if (!readParameter(param_file, omega, "omega")) return 1;

    
    const double pi = 3.141592653589793;
    double xinit[2] = {theta, omega}; // initial conditions [theta, wn]
    std::string fout = "simout.txt"; // output file name

    rk4 solver(f0, f1, xinit, t0, t_end, h); // create RK4 solver instance
    solver.run(fout); // run the solver
    
    return 0;
}