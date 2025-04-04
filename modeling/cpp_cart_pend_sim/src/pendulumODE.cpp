#include <cmath> // For std::sqrt
#include "readParameter.cpp" // For readParameter function
#include <iostream> // For std::cout, std::endl
#include <fstream> // For std::ifstream
#include <string> // For std::string

// ODE to be implemented in this class:
// f0 = dx0/dt = wn
// f1 = dx1/dt = -wc**2 * sin(theta) - 2*gamma*wn
// where wc = 1.0, gamma = 0.1, wn = 1.0
// and theta = x0, wn = x1
// This is a simple pendulum with damping

class PendulumODE {
private:
    double g; // gravitational acceleration
    double l; // length of the pendulum
    double wc; // natural frequency
    double beta; // damping ratio
    double m; // mass of the pendulum bob
    double gamma; // damping coefficient
public:
    // Constructor to initialize the parameters
    PendulumODE(const std::string& filename) {
        std::ifstream param_file(filename);
        if (!param_file.is_open()) {
            std::cerr << "Error opening " << filename << std::endl;
            exit(1); // Exit the program if the file cannot be opened
        }

        if (!readParameter(param_file, g, "g")) exit(1);
        if (!readParameter(param_file, l, "l")) exit(1);
        if (!readParameter(param_file, beta, "beta")) exit(1);
        if (!readParameter(param_file, m, "m")) exit(1); // what's the difference between return 1 and exit(1)?
        wc = std::sqrt(g / l); // natural frequency
        gamma = beta / (2 * m * l); // damping coefficient
    }

    // ODE function f0
    double f0(double t, double x0, double x1) {
        return x1;
    }

    // ODE function f1
    double f1(double t, double x0, double x1) {
        return -wc * wc * sin(x0) - 2 * gamma * x1;
    }

    // Debugging: Print the parameters
    void printParameters() const {
        std::cout << "g: " << g << ", l: " << l << ", beta: " << beta << ", m: " << m << std::endl;
    }
};

// Example usage (if needed for testing)
/*
int main() {
    PendulumODE ode("pendparm.txt");
    ode.printParameters();

    double t = 0.0, x0 = 1.0, x1 = 0.0;
    std::cout << "f0: " << ode.f0(t, x0, x1) << std::endl;
    std::cout << "f1: " << ode.f1(t, x0, x1) << std::endl;

    return 0;
}
*/
