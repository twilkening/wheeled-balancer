#include <iostream>
#include <fstream>
#include "rk4.h"
#include "equations.h"

int main() {
    std::ofstream file("simulation_output.csv");
    file << "time,x,v,theta,omega\n";  // CSV Header

    double t = 0.0, dt = 0.01, t_end = 10.0;  // Time settings
    State y = {0.0, 0.0, M_PI / 6, 0.0};  // Initial state: x, v, theta, omega

    while (t <= t_end) {
        file << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "\n";
        y = rk4_step(cart_pendulum_dynamics, t, y, dt);
        t += dt;
    }

    file.close();
    std::cout << "Simulation complete! Output saved to simulation_output.csv\n";
    return 0;
}
