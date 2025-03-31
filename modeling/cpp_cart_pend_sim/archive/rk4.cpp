#include "rk4.h"

State rk4_step(const DerivativeFunc& derivs, double t, const State& y, double dt) {
    State k1 = derivs(t, y);
    State y1 = y;
    for (size_t i = 0; i < y.size(); ++i) y1[i] += dt * k1[i] / 2.0;

    State k2 = derivs(t + dt / 2.0, y1);
    State y2 = y;
    for (size_t i = 0; i < y.size(); ++i) y2[i] += dt * k2[i] / 2.0;

    State k3 = derivs(t + dt / 2.0, y2);
    State y3 = y;
    for (size_t i = 0; i < y.size(); ++i) y3[i] += dt * k3[i];

    State k4 = derivs(t + dt, y3);

    State y_next = y;
    for (size_t i = 0; i < y.size(); ++i)
        y_next[i] += dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;

    return y_next;
}
