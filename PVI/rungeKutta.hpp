#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP
#include "pvi.hpp"

Estado rungeKutta4(FuncaoPVI F, Estado S0, double t0, double tf, double dt) {
    Estado S = S0;
    double t = t0;
    int n = calcularNumeroDePassos(t0, tf, dt);

    for (int i = 0; i < n; i++) {
        S = rk4Passo(F, S, t, dt);
        t += dt;
    }
    return S;
}
#endif