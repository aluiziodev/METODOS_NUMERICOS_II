#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
 
#include "PVI/pvi.hpp"
 
namespace {
 
// Método de Runge-Kutta de terceira ordem, equações (47) a (52) da Aula #25.
//
// F1 = F(Si, ti)
// Smeio = Si + (dt/2) * F1                    -> estimativa em ti + dt/2   (48)
// F2 = F(Smeio, ti + dt/2)                     -> (49)
// Sbar = Si + dt * (-F1 + 2*F2)                -> estimativa em ti + dt    (50)
// F3 = F(Sbar, ti + dt)                        -> (51)
// Si+1 = Si + dt * (1/6 F1 + 4/6 F2 + 1/6 F3)  -> atualização final       (52)
template <typename Sistema>
Estado rk3Passo(Sistema sistema, const Estado &Si, double ti, double dt) {
    const Estado F1 = sistema(Si, ti);
 
    Estado Smeio(Si.size());
    for (std::size_t k = 0; k < Si.size(); ++k) {
        Smeio[k] = Si[k] + (dt / 2.0) * F1[k];
    }
    const Estado F2 = sistema(Smeio, ti + dt / 2.0);
 
    Estado Sbar(Si.size());
    for (std::size_t k = 0; k < Si.size(); ++k) {
        Sbar[k] = Si[k] + dt * (-F1[k] + 2.0 * F2[k]);
    }
    const Estado F3 = sistema(Sbar, ti + dt);
 
    Estado Sprox(Si.size());
    for (std::size_t k = 0; k < Si.size(); ++k) {
        Sprox[k] = Si[k] + dt * ((1.0 / 6.0) * F1[k] + (4.0 / 6.0) * F2[k] +
                                  (1.0 / 6.0) * F3[k]);
    }
    return Sprox;
}
 
struct ResultadoTrajetoria {
    double ymax = 0.0;
    double tmax = 0.0;
    double ttotal = 0.0;
    double vImpacto = 0.0;
    bool encontrouMax = false;
    bool encontrouImpacto = false;
};
 
// PVI-2 modificado (Tarefa 17, item 1):
// t0 = 0s, v0 = 5 m/s, y0 = 200 m, k = 0.25 kg/s, m = 2 kg.
Estado sistemaQuedaLivre(const Estado &S, double t) {
    (void)t;
    const double v = S[0];
    const double y = S[1];
    (void)y;
 
    const double g = 9.81;
    const double k = 0.25;
    const double m = 2.0;
 
    Estado dS(2);
    dS[0] = -g - (k / m) * v;
    dS[1] = v;
    return dS;
}
 
ResultadoTrajetoria simularParaquedista(double dt) {
    const double t0 = 0.0;
    const double tFinalMax = 100.0;
 
    Estado S{5.0, 200.0};
    double t = t0;
 
    ResultadoTrajetoria resultado;
    resultado.ymax = S[1];
    resultado.tmax = t;
 
    while (t < tFinalMax) {
        const Estado Sprev = S;
        const double tprev = t;
        S = rk3Passo(sistemaQuedaLivre, Sprev, tprev, dt);
        t += dt;
 
        const double vPrev = Sprev[0];
        const double yPrev = Sprev[1];
        const double vNext = S[0];
        const double yNext = S[1];
 
        if (!resultado.encontrouMax && vPrev > 0.0 && vNext <= 0.0) {
            const double alpha = (0.0 - vPrev) / (vNext - vPrev);
            resultado.tmax = tprev + alpha * dt;
            resultado.ymax = yPrev + alpha * (yNext - yPrev);
            resultado.encontrouMax = true;
        }
 
        if (!resultado.encontrouImpacto && yPrev > 0.0 && yNext <= 0.0) {
            const double alpha = (0.0 - yPrev) / (yNext - yPrev);
            resultado.ttotal = tprev + alpha * dt;
            resultado.vImpacto = vPrev + alpha * (vNext - vPrev);
            resultado.encontrouImpacto = true;
            break;
        }
    }
 
    if (!resultado.encontrouMax) {
        resultado.tmax = t;
        resultado.ymax = S[1];
    }
 
    if (!resultado.encontrouImpacto) {
        resultado.ttotal = t;
        resultado.vImpacto = S[0];
    }
 
    return resultado;
}
 
}  // namespace
 
int main() {
    std::cout << "Tarefa 17: Runge-Kutta de terceira ordem (eqs. 47-52)\n";
    std::cout << "PVI-2 modificado: v0=5 m/s, y0=200 m, k=0.25 kg/s, m=2 kg\n";
    std::cout << "======================================================\n";
    std::cout << std::left
              << std::setw(12) << "dt"
              << std::setw(16) << "ymax"
              << std::setw(16) << "tmax"
              << std::setw(16) << "ttotal"
              << std::setw(16) << "v_impacto" << '\n';
 
    // Valores de Δt da Tabela 1 (Aula #25).
    const std::vector<double> passos = {0.1, 0.01, 0.001, 0.0001};
 
    for (double dt : passos) {
        const ResultadoTrajetoria r = simularParaquedista(dt);
        std::cout << std::fixed << std::setprecision(6)
                  << std::setw(12) << dt
                  << std::setw(16) << r.ymax
                  << std::setw(16) << r.tmax
                  << std::setw(16) << r.ttotal
                  << std::setw(16) << r.vImpacto << '\n';
    }
 
    return 0;
}