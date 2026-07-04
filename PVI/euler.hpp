#ifndef EULER_HPP
#define EULER_HPP
#include "pvi.hpp"

// EULER EXPLÍCITO: S_{i+1} = S_i + dt * F(S_i, t_i)
Estado eulerExplicito(FuncaoPVI F, Estado S0, double t0, double tf, double dt) {
    Estado S = S0;
    double t = t0;
    int n = (int)std::round((tf - t0) / dt);
    for (int i = 0; i < n; i++) {
        S = somar(S, escalar(dt, F(S, t)));
        t += dt;
    }
    return S;
}

// EULER IMPLÍCITO: S_{i+1} = S_i + dt * F(S_{i+1}, t_{i+1})
// Para o PVI-1 (linear), a equação pode ser isolada algebricamente:
//   y_{i+1} = y_i / (1 - dt*(2/3))
double eulerImplicitoPVI1(double y0, double t0, double tf, double dt) {
    double y = y0, t = t0;
    int n = (int)std::round((tf - t0) / dt);
    for (int i = 0; i < n; i++) {
        y = y / (1.0 - dt * (2.0/3.0));
        t += dt;
    }
    return y;
}

// Para o PVI-2 (sistema linear 2x2), resolvido analiticamente por passo
// (mesma dedução da Aula #24, Seção 4.2.2)
Estado eulerImplicitoPVI2(Estado S0, double t0, double tf, double dt,
                           double g, double k, double m) {
    Estado S = S0;
    double t = t0;
    int n = (int)std::round((tf - t0) / dt);
    for (int i = 0; i < n; i++) {
        double v = S[0], y = S[1];
        double vNovo = (m / (m + k*dt)) * (v - g*dt);
        double yNovo = y + (m*dt / (m + k*dt)) * (v - g*dt);
        S = {vNovo, yNovo};
        t += dt;
    }
    return S;
}
#endif