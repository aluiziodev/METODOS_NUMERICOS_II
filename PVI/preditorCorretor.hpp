#ifndef PREDITOR_CORRETOR_HPP
#define PREDITOR_CORRETOR_HPP
#include "pvi.hpp"
#include "rungeKutta.hpp"
#include <vector>

// Fase de inicialização (RK4) + predição (Adams-Bashforth 4a ordem) +
// correção (Adams-Moulton 4a ordem)
Estado preditorCorretor4(FuncaoPVI F, Estado S0, double t0, double tf, double dt) {
    int n = (int)std::round((tf - t0) / dt);

    std::vector<Estado> S(n+1);
    std::vector<double> t(n+1);
    std::vector<Estado> Fs(n+1);

    S[0] = S0; t[0] = t0;

    // Inicialização: os 3 primeiros passos via Runge-Kutta 4
    for (int i = 0; i < 3 && i < n; i++) {
        Estado F1 = F(S[i], t[i]);
        Estado F2 = F(somar(S[i], escalar(dt/2.0, F1)), t[i] + dt/2.0);
        Estado F3 = F(somar(S[i], escalar(dt/2.0, F2)), t[i] + dt/2.0);
        Estado F4 = F(somar(S[i], escalar(dt, F3)), t[i] + dt);

        Estado soma(S[i].size());
        for (size_t j = 0; j < S[i].size(); j++)
            soma[j] = (F1[j] + 2*F2[j] + 2*F3[j] + F4[j]);

        S[i+1] = somar(S[i], escalar(dt/6.0, soma));
        t[i+1] = t[i] + dt;
    }

    for (int i = 0; i <= std::min(3, n); i++) Fs[i] = F(S[i], t[i]);

    // Fase de predição + correção (a partir do 4o ponto)
    for (int i = 3; i < n; i++) {
        Estado &Fi = Fs[i], &Fim1 = Fs[i-1], &Fim2 = Fs[i-2], &Fim3 = Fs[i-3];

        // Preditor (Adams-Bashforth 4a ordem)
        Estado soma(S[i].size());
        for (size_t j = 0; j < S[i].size(); j++)
            soma[j] = 55*Fi[j] - 59*Fim1[j] + 37*Fim2[j] - 9*Fim3[j];
        Estado Spred = somar(S[i], escalar(dt/24.0, soma));
        double tProx = t[i] + dt;
        Estado Fpred = F(Spred, tProx);

        // Corretor (Adams-Moulton 4a ordem)
        Estado somaCorr(S[i].size());
        for (size_t j = 0; j < S[i].size(); j++)
            somaCorr[j] = 9*Fpred[j] + 19*Fi[j] - 5*Fim1[j] + Fim2[j];
        Estado Scorr = somar(S[i], escalar(dt/24.0, somaCorr));

        S[i+1] = Scorr;
        t[i+1] = tProx;
        Fs[i+1] = F(Scorr, tProx);
    }

    return S[n];
}
#endif