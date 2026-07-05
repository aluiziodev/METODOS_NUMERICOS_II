#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "PVI/euler.hpp"
#include "PVI/preditorCorretor.hpp"
#include "PVI/rungeKutta.hpp"

using Estado = std::vector<double>;

namespace {
constexpr double g = 9.81;
constexpr double c = 1.0;
constexpr double m = 2.0;

double fPVI1(double t, double y) {
    (void)t;
    return (2.0 / 3.0) * y;
}

std::vector<double> fPVI2(double t, const std::vector<double> &S) {
    (void)t;
    const double v = S[0];
    const double y = S[1];
    (void)y;

    std::vector<double> dS(2);
    dS[0] = -g - (c / m) * v;
    dS[1] = v;
    return dS;
}

Estado fPVI1Wrapper(const Estado &S, double t) {
    return {fPVI1(t, S[0])};
}

Estado fPVI2Wrapper(const Estado &S, double t) {
    return fPVI2(t, S);
}

bool verificarErro(const std::string &nome, double aproximado, double exato, double tolerancia) {
    const double erro = std::abs(aproximado - exato);
    std::cout << std::fixed << std::setprecision(8)
              << nome << " -> aproximado=" << aproximado
              << ", exato=" << exato
              << ", erro=" << erro << '\n';
    return erro <= tolerancia;
}

bool verificarErroEstado(const std::string &nome, const Estado &aproximado, const Estado &exato, double tolerancia) {
    bool ok = true;
    for (std::size_t i = 0; i < aproximado.size(); ++i) {
        ok = ok && verificarErro(nome + "[" + std::to_string(i) + "]", aproximado[i], exato[i], tolerancia);
    }
    return ok;
}

bool testarPVI1() {
    std::cout << "\nCaso teste: PVI 1\n";
    std::cout << "-----------------\n";

    const double t0 = 0.0;
    const double tf = 1.0;
    const double dt = 0.1;
    const Estado y0{1.0};

    const Estado yEuler = eulerExplicito(fPVI1Wrapper, y0, t0, tf, dt);
    const double exatoEuler = std::exp((2.0 / 3.0) * tf);
    const bool okEuler = verificarErro("Euler explícito", yEuler[0], exatoEuler, 0.08);

    const Estado yRK4 = rungeKutta4(fPVI1Wrapper, y0, t0, tf, dt);
    const bool okRK4 = verificarErro("Runge-Kutta 4", yRK4[0], exatoEuler, 1e-4);

    const Estado yPC4 = preditorCorretor4(fPVI1Wrapper, y0, t0, tf, dt);
    const bool okPC4 = verificarErro("Preditor-Corretor 4", yPC4[0], exatoEuler, 1e-4);

    const bool ok = okEuler && okRK4 && okPC4;
    std::cout << "Resultado PVI 1: " << (ok ? "passou" : "falhou") << "\n";
    return ok;
}

bool testarPVI2() {
    std::cout << "\nCaso teste: PVI 2\n";
    std::cout << "-----------------\n";

    const double t0 = 0.0;
    const double tf = 1.0;
    const double dt = 0.01;
    const Estado estado0{0.0, 0.0};

    const Estado yEuler = eulerExplicito(fPVI2Wrapper, estado0, t0, tf, dt);
    const double vExato = -(m * g / c) * (1.0 - std::exp(-(c / m) * tf));
    const double yExato = -(m * g / c) * tf + (m * m * g / (c * c)) * (1.0 - std::exp(-(c / m) * tf));
    const Estado exato{vExato, yExato};
    const bool okEuler = verificarErroEstado("Euler explícito", yEuler, exato, 0.05);

    const Estado yRK4 = rungeKutta4(fPVI2Wrapper, estado0, t0, tf, dt);
    const bool okRK4 = verificarErroEstado("Runge-Kutta 4", yRK4, exato, 1e-3);

    const Estado yPC4 = preditorCorretor4(fPVI2Wrapper, estado0, t0, tf, dt);
    const bool okPC4 = verificarErroEstado("Preditor-Corretor 4", yPC4, exato, 1e-3);

    const bool ok = okEuler && okRK4 && okPC4;
    std::cout << "Resultado PVI 2: " << (ok ? "passou" : "falhou") << "\n";
    return ok;
}
}  // namespace

int main() {
    std::cout << "Teste dos métodos de PVI\n";
    std::cout << "========================\n";

    const bool okPVI1 = testarPVI1();
    const bool okPVI2 = testarPVI2();

    std::cout << "\nResumo: " << ((okPVI1 && okPVI2) ? "todos os testes passaram" : "alguns testes falharam") << "\n";
    return (okPVI1 && okPVI2) ? 0 : 1;
}
