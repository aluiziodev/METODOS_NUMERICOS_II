#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "PVI/pvi.hpp"
#include "PVI/rungeKutta.hpp"

// -----------------------------------------------------------------------
// Tarefa 18 (Aula #26): Método preditor-corretor de quarta ordem
//
// Parte 1 - o método (mesmas fórmulas de preditorCorretor.hpp, obtidas pelo
// mesmo processo usado para deduzir o preditor-corretor de 3a ordem: g(t) é
// um polinômio de interpolação passando pelas derivadas dos k+1 = 4 últimos
// estados; aqui reaproveitamos diretamente essas fórmulas):
//
//   Inicialização: S1, S2, S3 obtidos via Runge-Kutta de 4a ordem (T1-T8)
//
//   Preditor (Adams-Bashforth 4a ordem):
//     Sbar_{i+1} = S_i + dt/24 * (55 F_i - 59 F_{i-1} + 37 F_{i-2} - 9 F_{i-3})
//
//   Corretor (Adams-Moulton 4a ordem):
//     S_{i+1} = S_i + dt/24 * (9 F(Sbar_{i+1}) + 19 F_i - 5 F_{i-1} + F_{i-2})
//
// Parte 2 - repetir a Tarefa 17 (Aula #25): trajetória de um paraquedista,
// obtendo ymax, tmax, ttotal e a velocidade de impacto no mar, para os
// valores de Δt da Tabela 1 (0.1, 0.01, 0.001, 0.0001).
// -----------------------------------------------------------------------

namespace {

struct ResultadoTrajetoria {
    double ymax = 0.0;
    double tmax = 0.0;
    double ttotal = 0.0;
    double vImpacto = 0.0;
    bool encontrouMax = false;
    bool encontrouImpacto = false;
};

// PVI-2 modificado (mesmo problema da Tarefa 17):
// t0 = 0s, v0 = 5 m/s, y0 = 200 m, k = 0.25 kg/s, m = 2 kg.
Estado sistemaQuedaLivre(const Estado &S, double t) {
    (void)t;
    const double v = S[0];

    const double g = 9.81;
    const double k = 0.25;
    const double m = 2.0;

    Estado dS(2);
    dS[0] = -g - (k / m) * v;
    dS[1] = v;
    return dS;
}

// Verifica, entre dois estados consecutivos, se a altura máxima (v muda de
// sinal) ou o impacto no mar (y cruza zero) ocorreram, refinando o instante
// por interpolação linear. Retorna true se o impacto foi encontrado.
bool atualizarEventos(const Estado &Sprev, const Estado &Snext, double tprev,
                       double dt, ResultadoTrajetoria &resultado) {
    const double vPrev = Sprev[0];
    const double yPrev = Sprev[1];
    const double vNext = Snext[0];
    const double yNext = Snext[1];

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
        return true;
    }
    return false;
}

ResultadoTrajetoria simularParaquedista(double dt) {
    const double t0 = 0.0;
    const double tFinalMax = 200.0;

    std::vector<Estado> S;
    std::vector<double> t;
    std::vector<Estado> Fs;

    S.push_back(Estado{5.0, 200.0});
    t.push_back(t0);
    Fs.push_back(sistemaQuedaLivre(S[0], t[0]));

    ResultadoTrajetoria resultado;
    resultado.ymax = S[0][1];
    resultado.tmax = t[0];

    bool impacto = false;

    // Fase de inicialização: os 3 primeiros passos via Runge-Kutta de 4a
    // ordem (fórmulas T1-T8), como pedido no enunciado.
    for (int i = 0; i < 3 && !impacto; i++) {
        Estado Snext = rk4Passo(sistemaQuedaLivre, S[i], t[i], dt);
        double tNext = t[i] + dt;

        impacto = atualizarEventos(S[i], Snext, t[i], dt, resultado);

        S.push_back(Snext);
        t.push_back(tNext);
        Fs.push_back(sistemaQuedaLivre(Snext, tNext));
    }

    // Fase de predição + correção (Adams-Bashforth-Moulton de 4a ordem),
    // a partir do 4o ponto (i = 3), igual a preditorCorretor.hpp.
    int i = 3;
    while (!impacto && static_cast<int>(S.size()) > i && t[i] < tFinalMax) {
        const Estado &Fi = Fs[i], &Fim1 = Fs[i - 1], &Fim2 = Fs[i - 2],
                     &Fim3 = Fs[i - 3];

        // Preditor (Adams-Bashforth 4a ordem)
        Estado soma(S[i].size());
        for (std::size_t j = 0; j < S[i].size(); j++) {
            soma[j] = 55.0 * Fi[j] - 59.0 * Fim1[j] + 37.0 * Fim2[j] -
                      9.0 * Fim3[j];
        }
        Estado Spred = somar(S[i], escalar(dt / 24.0, soma));
        double tProx = t[i] + dt;
        Estado Fpred = sistemaQuedaLivre(Spred, tProx);

        // Corretor (Adams-Moulton 4a ordem)
        Estado somaCorr(S[i].size());
        for (std::size_t j = 0; j < S[i].size(); j++) {
            somaCorr[j] = 9.0 * Fpred[j] + 19.0 * Fi[j] - 5.0 * Fim1[j] +
                          Fim2[j];
        }
        Estado Scorr = somar(S[i], escalar(dt / 24.0, somaCorr));

        impacto = atualizarEventos(S[i], Scorr, t[i], dt, resultado);

        S.push_back(Scorr);
        t.push_back(tProx);
        Fs.push_back(sistemaQuedaLivre(Scorr, tProx));
        i++;
    }

    if (!resultado.encontrouMax) {
        resultado.tmax = t.back();
        resultado.ymax = S.back()[1];
    }

    if (!resultado.encontrouImpacto) {
        resultado.ttotal = t.back();
        resultado.vImpacto = S.back()[0];
    }

    return resultado;
}

}  // namespace

int main() {
    std::cout << "Tarefa 18: Preditor-Corretor de quarta ordem "
                 "(Adams-Bashforth-Moulton)\n";
    std::cout << "Inicializacao via Runge-Kutta de 4a ordem (eqs. T1-T8)\n";
    std::cout << "PVI-2 modificado: v0=5 m/s, y0=200 m, k=0.25 kg/s, m=2 kg\n";
    std::cout << "==============================================================\n";
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