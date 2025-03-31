#include <cmath>
#include <vector>

class rk4
{
    double xn[2]; // current states
    double k1[2], k2[2], k3[2], k4[2]; // RK4 coefficients
    const double h; // step size
    double t; // current time
    double t_end; // end time
    double (*f)(double, double); // function pointer for the ODE
public:
    rk4(double (*f)(double, double), double xinit[2], double t0, double t_end, double h)
        : f(f), t(t0), t_end(t_end), h(h)
    {
        // Initialize the states
        xn[0] = xinit[0];
        xn[1] = xinit[1];
    }

    void step()
    {
        k1[0] = h * f(t, xn[0]);
        k1[1] = h * f(t, xn[1]);

        k2[0] = h * f(t + h / 2, xn[0] + k1[0] / 2);
        k2[1] = h * f(t + h / 2, xn[1] + k1[1] / 2);

        k3[0] = h * f(t + h / 2, xn[0] + k2[0] / 2);
        k3[1] = h * f(t + h / 2, xn[1] + k2[1] / 2);

        k4[0] = h * f(t + h, xn[0] + k3[0]);
        k4[1] = h * f(t + h, xn[1] + k3[1]);

        xn[0] += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
        xn[1] += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;

        t += h;
    }

    bool isComplete() const
    {
        return t >= t_end;
    }

    void run()
    {
        while (t < t_end) {
            step();
        }
    }

    std::vector<double> getState() const
    {
        return {xn[0], xn[1]};
    }

}