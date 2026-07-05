/*
 * main_fem.cpp
 * ------------
 * Resolve, usando o header fem.hpp (Método dos Elementos Finitos,
 * formulação de Galerkin com elementos lineares), o PVC:
 *
 *   u''(x) + 7 u'(x) - u(x) = 2,   x em [0, 2]
 *   u(0) = 10
 *   u(2) = 1
 *
 * com N = 20 elementos (h = 0.1), e imprime a tabela de valores
 * aproximados u(0.1), u(0.2), ..., u(1.9), incluindo os valores de
 * contorno u(0) = 10 e u(2) = 1. Grava também "solucao_fem.csv".
 *
 * Compilar:   g++ -std=c++17 -O2 -Wall -Wextra -o solve_pvc_fem main_fem.cpp
 * Executar:   ./solve_pvc_fem
 */

#include <iostream>
#include "PVC/metodoElemFinitos.hpp"

int main() {
    // Coeficientes de a(x) u'' + b(x) u' + c(x) u = f(x):
    fem::BVP2 problema{
        [](double /*x*/) { return 1.0; },   // a(x): coeficiente de u''
        [](double /*x*/) { return 7.0; },   // b(x): coeficiente de u'
        [](double /*x*/) { return -1.0; },  // c(x): coeficiente de u
        [](double /*x*/) { return 2.0; }    // f(x): termo fonte
    };

    const double x0 = 0.0, xN = 2.0;
    const double u0 = 10.0, uN = 1.0;
    const double dx = 0.1;

    try {
        fem::GalerkinFEMSolver solver(problema);
        fem::Result r = solver.solve(x0, xN, u0, uN, dx);
        std::cout << "Metodo dos Elementos Finitos (FEM) aplicado ao PVC:\n";
        std::cout << "Solucao aproximada de u''(x) + 7 u'(x) - u(x) = 2\n";
        std::cout << "u(0) = " << u0 << " , u(2) = " << uN << " , h = " << r.h << "\n";
        std::cout << "Numero de elementos: " << r.n << "\n\n";
        r.print();

    } catch (const std::exception &e) {
        std::cerr << "Erro ao resolver o PVC: " << e.what() << "\n";
        return 1;
    }

    return 0;
}