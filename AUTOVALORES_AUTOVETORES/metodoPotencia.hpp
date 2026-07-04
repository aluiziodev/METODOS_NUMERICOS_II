#ifndef METODO_POTENCIA_HPP
#define METODO_POTENCIA_HPP
#include "algebraLinear.hpp"
#include <iostream>

struct ResultadoAutovalor {
    double autovalor;
    Vetor autovetor;
    int iteracoes;
};

// MÉTODO DA POTÊNCIA REGULAR: converge para o autovalor de maior módulo
ResultadoAutovalor metodoPotencia(const Matriz &A, Vetor x0,
                                   double tol = 1e-8, int maxIter = 200) {
    Vetor x = normalizar(x0);
    double lambdaAntigo = 0.0;

    for (int k = 0; k < maxIter; k++) {
        Vetor y = matVec(A, x);
        x = normalizar(y);
        double lambda = produtoInterno(x, matVec(A, x)); // quociente de Rayleigh

        if (fabs(lambda - lambdaAntigo) < tol)
            return {lambda, x, k+1};
        lambdaAntigo = lambda;
    }
    return {lambdaAntigo, x, maxIter};
}

// MÉTODO DA POTÊNCIA INVERSA: converge para o autovalor de MENOR módulo
// (equivale à potência regular aplicada a A^-1)
ResultadoAutovalor metodoPotenciaInversa(const Matriz &A, Vetor x0,
                                          double tol = 1e-8, int maxIter = 200) {
    Vetor x = normalizar(x0);
    double lambdaAntigo = 0.0;

    for (int k = 0; k < maxIter; k++) {
        Vetor y = resolverSistema(A, x);   // resolve A*y = x  (equivale a y = A^-1 x)
        x = normalizar(y);
        double lambda = produtoInterno(x, matVec(A, x));

        if (fabs(lambda - lambdaAntigo) < tol)
            return {lambda, x, k+1};
        lambdaAntigo = lambda;
    }
    return {lambdaAntigo, x, maxIter};
}

// MÉTODO DA POTÊNCIA COM DESLOCAMENTO: converge para o autovalor
// mais PRÓXIMO do deslocamento (shift) escolhido
ResultadoAutovalor metodoPotenciaDeslocamento(const Matriz &A, Vetor x0,
                                               double shift,
                                               double tol = 1e-8, int maxIter = 200) {
    Matriz B = subtrairEscalarDiagonal(A, shift); // B = A - shift*I
    Vetor x = normalizar(x0);
    double lambdaAntigo = 0.0;

    for (int k = 0; k < maxIter; k++) {
        Vetor y = resolverSistema(B, x);
        x = normalizar(y);
        double lambda = produtoInterno(x, matVec(A, x)); // Rayleigh sobre A, não B

        if (fabs(lambda - lambdaAntigo) < tol)
            return {lambda, x, k+1};
        lambdaAntigo = lambda;
    }
    return {lambdaAntigo, x, maxIter};
}
#endif