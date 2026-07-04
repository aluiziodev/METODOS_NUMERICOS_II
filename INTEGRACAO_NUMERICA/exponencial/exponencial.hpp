#ifndef EXPONENCIAL_HPP
#define EXPONENCIAL_HPP
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include "../../utils/function.hpp"

// QUADRATURA EXPONENCIAL SIMPLES (SE)
// Transformação: x = mid + half*tanh(t),  t in (-inf, inf)
// Boa para integrais com singularidade leve nas extremidades.

double intgExpSimples(Function &f, double a, double b, double h, int N) {
    double mid = (a + b) / 2.0;
    double half = (b - a) / 2.0;
    double sum = 0.0;
    for (int k = -N; k <= N; k++) {
        double t = k * h;
        double phi = tanh(t);
        double dphi = 1.0 / (cosh(t) * cosh(t));
        double x = mid + half * phi;
        sum += f.value(x) * dphi;
    }
    return half * h * sum;
}

// QUADRATURA EXPONENCIAL DUPLA (DE) - Takahasi-Mori
// Transformação: x = mid + half*tanh((pi/2)*sinh(t))
// Convergência muito mais rápida (dupla-exponencial), ideal para
// integrais com singularidade forte ou fraca nas extremidades.

double intgExpDupla(Function &f, double a, double b, double h, int N) {
    double mid = (a + b) / 2.0;
    double half = (b - a) / 2.0;
    double sum = 0.0;
    for (int k = -N; k <= N; k++) {
        double t = k * h;
        double s = sinh(t);
        double phi = tanh((M_PI / 2.0) * s);
        double dphi = (M_PI / 2.0) * cosh(t) / pow(cosh((M_PI / 2.0) * s), 2);
        double x = mid + half * phi;
        sum += f.value(x) * dphi;
    }
    return half * h * sum;
}

// Teste de convergência: mantém o truncamento T = N*h fixo e refina h,
// dobrando N a cada iteração (analogamente aos testes de Newton-Cotes).
void testExponencial(double (*metodo)(Function&, double, double, double, int),
                      Function &f, double a, double b, double T, double tol,
                      int maxIter = 20) {
    double h = T / 8.0;
    int N = (int)(T / h);
    double i_prev = metodo(f, a, b, h, N);
    double erro = 1.0;
    int iter = 0;
    std::cout << "Iteracao: " << iter << ", h: " << h << ", N: " << N
               << ", Integral: " << i_prev << ", Erro: N/A" << std::endl;
    while (erro > tol && iter < maxIter) {
        h /= 2.0;
        N = (int)(T / h);
        double i_cur = metodo(f, a, b, h, N);
        erro = fabs(i_cur - i_prev);
        i_prev = i_cur;
        iter++;
        std::cout << "Iteracao: " << iter << ", h: " << h << ", N: " << N
                   << ", Integral: " << i_prev << ", Erro: " << erro << std::endl;
    }
}
#endif