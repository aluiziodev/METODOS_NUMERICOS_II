#ifndef PVI_HPP
#define PVI_HPP
#include <vector>
#include <functional>
#include <cmath>
#include <cstddef>
#include <stdexcept>

using Estado = std::vector<double>;
using FuncaoPVI = std::function<Estado(const Estado&, double)>;

inline Estado somar(const Estado &a, const Estado &b) {
    Estado r(a.size());
    for (size_t i = 0; i < a.size(); i++) r[i] = a[i] + b[i];
    return r;
}

inline Estado escalar(double c, const Estado &a) {
    Estado r(a.size());
    for (size_t i = 0; i < a.size(); i++) r[i] = c * a[i];
    return r;
}

inline int calcularNumeroDePassos(double t0, double tf, double dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument("dt deve ser positivo");
    }
    return std::max(0, static_cast<int>(std::llround((tf - t0) / dt)));
}

inline Estado rk4Passo(const FuncaoPVI &F, const Estado &S, double t, double dt) {
    Estado F1 = F(S, t);
    Estado F2 = F(somar(S, escalar(dt / 2.0, F1)), t + dt / 2.0);
    Estado F3 = F(somar(S, escalar(dt / 2.0, F2)), t + dt / 2.0);
    Estado F4 = F(somar(S, escalar(dt, F3)), t + dt);

    Estado soma(S.size());
    for (size_t j = 0; j < S.size(); j++) {
        soma[j] = F1[j] + 2.0 * F2[j] + 2.0 * F3[j] + F4[j];
    }

    return somar(S, escalar(dt / 6.0, soma));
}
#endif