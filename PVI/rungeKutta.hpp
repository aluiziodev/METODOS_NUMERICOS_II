#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP
#include "pvi.hpp"

Estado rungeKutta4(FuncaoPVI F, Estado S0, double t0, double tf, double dt) {
    Estado S = S0;
    double t = t0;
    int n = (int)std::round((tf - t0) / dt);

    for (int i = 0; i < n; i++) {
        Estado F1 = F(S, t);
        Estado F2 = F(somar(S, escalar(dt/2.0, F1)), t + dt/2.0);
        Estado F3 = F(somar(S, escalar(dt/2.0, F2)), t + dt/2.0);
        Estado F4 = F(somar(S, escalar(dt, F3)), t + dt);

        Estado soma(S.size());
        for (size_t j = 0; j < S.size(); j++)
            soma[j] = (F1[j] + 2*F2[j] + 2*F3[j] + F4[j]);

        S = somar(S, escalar(dt/6.0, soma));
        t += dt;
    }
    return S;
}
#endif