#ifndef METODO_POTENCIA_HPP
#define METODO_POTENCIA_HPP
#include "../utils/algebraLinear.hpp"
#include <iostream>

struct ResultadoAutovalor {
    double autovalor;
    Vetor autovetor;
    int iteracoes;
};

ResultadoAutovalor metodoPotencia(const Matriz &A, Vetor x0, double tol = 1e-8, int maxIter = 200) {
    Vetor x = normalizar(x0);
    double lambdaAntigo = 0.0;

    for (int k = 0; k < maxIter; k++) {
        Vetor y = matVec(A, x);
        x = normalizar(y);
        double lambda = produtoInterno(x, matVec(A, x));

        if (fabs(lambda - lambdaAntigo) < tol)
            return {lambda, x, k+1};
        lambdaAntigo = lambda;
    }
    return {lambdaAntigo, x, maxIter};
}

ResultadoAutovalor metodoPotenciaInversa(const Matriz &A, Vetor x0, double tol = 1e-8, int maxIter = 200) {
    Vetor x = normalizar(x0);
    double lambdaAntigo = 0.0;

    for (int k = 0; k < maxIter; k++) {
        Vetor y = resolverSistema(A, x); 
        x = normalizar(y);
        double lambda = produtoInterno(x, matVec(A, x));

        if (fabs(lambda - lambdaAntigo) < tol)
            return {lambda, x, k+1};
        lambdaAntigo = lambda;
    }
    return {lambdaAntigo, x, maxIter};
}


ResultadoAutovalor metodoPotenciaDeslocamento(const Matriz &A, Vetor x0, double shift, double tol = 1e-8, int maxIter = 200) {
    Matriz B = subtrairEscalarDiagonal(A, shift); 
    Vetor x = normalizar(x0);
    double lambdaAntigo = 0.0;

    for (int k = 0; k < maxIter; k++) {
        Vetor y = resolverSistema(B, x);
        x = normalizar(y);
        double lambda = produtoInterno(x, matVec(A, x)); 

        if (fabs(lambda - lambdaAntigo) < tol)
            return {lambda, x, k+1};
        lambdaAntigo = lambda;
    }
    return {lambdaAntigo, x, maxIter};
}
#endif