#ifndef RK4_H
#define RK4_H

#include <vector>
#include <functional>

using State = std::vector<double>;
using DerivativeFunc = std::function<State(double, const State&)>;

State rk4_step(const DerivativeFunc& derivs, double t, const State& y, double dt);

#endif
