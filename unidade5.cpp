/*
 * main.cpp
 * --------
 * Resolve, usando o header fdm.hpp (Método das Diferenças Finitas), o PVC:
 *
 *   u''(x) + 7 u'(x) - u(x) = 2,   x em [0, 2]
 *   u(0) = 10
 *   u(2) = 1
 *
 * com passo dx = 0.1, e imprime a tabela de valores aproximados
 * u(0.1), u(0.2), ..., u(1.9), incluindo também os valores de contorno
 * u(0) = 10 e u(2) = 1. Também grava o resultado em "solucao.csv".
 *
 * Compilar:   g++ -std=c++17 -O2 -Wall -Wextra -o solve_pvc main.cpp
 * Executar:   ./solve_pvc
 */

#include <iostream>
#include "PVC/metodoDifFinitas.hpp"

int main() {
    // Coeficientes de a(x) u'' + b(x) u' + c(x) u = f(x):
    fdm::BVP2 problema{
        [](double /*x*/) { return 1.0; },   // a(x): coeficiente de u''
        [](double /*x*/) { return 7.0; },   // b(x): coeficiente de u'
        [](double /*x*/) { return -1.0; },  // c(x): coeficiente de u
        [](double /*x*/) { return 2.0; }    // f(x): termo fonte
    };

    const double x0 = 0.0, xN = 2.0;
    const double u0 = 10.0, uN = 1.0;
    const double dx = 0.1;

    try {
        fdm::FiniteDifferenceSolver solver(problema);
        fdm::Result r = solver.solve(x0, xN, u0, uN, dx);
        std::cout << "Metodo das Diferencas Finitas (FDM) aplicado ao PVC:\n";
        std::cout << "Solucao aproximada de u''(x) + 7 u'(x) - u(x) = 2\n";
        std::cout << "u(0) = " << u0 << " , u(2) = " << uN << " , dx = " << dx << "\n";
        std::cout << "Numero de subintervalos: " << r.n << "\n\n";
        r.print();

    } catch (const std::exception &e) {
        std::cerr << "Erro ao resolver o PVC: " << e.what() << "\n";
        return 1;
    }

    return 0;
}