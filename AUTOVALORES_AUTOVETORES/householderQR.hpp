#ifndef HOUSEHOLDER_QR_HPP
#define HOUSEHOLDER_QR_HPP
#include "../utils/algebraLinear.hpp"


struct ResultadoQR {
    Vetor autovalores;
    int iteracoes;
};



Matriz multiplicarMatrizes(const Matriz &A, const Matriz &B) {
    int n = A.size();
    Matriz C(n, Vetor(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

void imprimirMatriz(const Matriz &A) {
    for (const auto &linha : A) {
        for (double v : linha) printf("%10.6f ", v);
        printf("\n");
    }
}


Matriz householderTridiagonalizacao(Matriz A, Matriz &Hacumulada) {
    int n = A.size();
    Hacumulada = identidade(n);

    for (int k = 0; k < n-2; k++) {
        Vetor x(n-k-1);
        for (int i = k+1; i < n; i++) x[i-k-1] = A[i][k];

        double normaX = norma(x);
        double alpha = (x[0] >= 0) ? -normaX : normaX;

        Vetor v = x;
        v[0] -= alpha;
        v = normalizar(v);

        Matriz H = identidade(n);
        for (int i = 0; i < (int)v.size(); i++)
            for (int j = 0; j < (int)v.size(); j++)
                H[k+1+i][k+1+j] -= 2.0 * v[i] * v[j];

        Matriz HA(n, Vetor(n, 0.0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int l = 0; l < n; l++)
                    HA[i][j] += H[i][l] * A[l][j];

        Matriz novoA(n, Vetor(n, 0.0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int l = 0; l < n; l++)
                    novoA[i][j] += HA[i][l] * H[l][j];

        A = novoA;
        Hacumulada = multiplicarMatrizes(Hacumulada, H); // H = H1*H2*H3*...
    }
    return A;
}


void passoQRGivens(Matriz &T) {
    int n = T.size();
    Matriz Q = identidade(n);
    Matriz R = T;

    for (int i = 0; i < n-1; i++) {
        double a = R[i][i], b = R[i+1][i];
        double r = sqrt(a*a + b*b);
        double c = a/r, s = b/r;

        for (int j = 0; j < n; j++) {
            double r1 = R[i][j], r2 = R[i+1][j];
            R[i][j]   =  c*r1 + s*r2;
            R[i+1][j] = -s*r1 + c*r2;
        }
        for (int j = 0; j < n; j++) {
            double q1 = Q[j][i], q2 = Q[j][i+1];
            Q[j][i]   =  c*q1 + s*q2;
            Q[j][i+1] =  -s*q1 + c*q2;
            
        }
    }

    // T = R * Q
    Matriz novoT(n, Vetor(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int l = 0; l < n; l++)
                novoT[i][j] += R[i][l] * Q[l][j];
    T = novoT;
}



ResultadoQR iteracaoQR(Matriz T, double tol = 1e-10, int maxIter = 500) {
    int n = T.size();
    for (int k = 0; k < maxIter; k++) {
        passoQRGivens(T);

        double somaSubdiagonal = 0.0;
        for (int i = 0; i < n-1; i++) somaSubdiagonal += fabs(T[i+1][i]);

        if (somaSubdiagonal < tol) {
            Vetor autovalores(n);
            for (int i = 0; i < n; i++) autovalores[i] = T[i][i];
            return {autovalores, k+1};
        }
    }
    Vetor autovalores(n);
    for (int i = 0; i < n; i++) autovalores[i] = T[i][i];
    return {autovalores, maxIter};
}


ResultadoQR autovaloresSimetrica(const Matriz &A, double tol = 1e-10, int maxIter = 500) {
    Matriz Hdescartavel; // nao usada aqui, so para satisfazer a assinatura
    Matriz T = householderTridiagonalizacao(A, Hdescartavel);
    return iteracaoQR(T, tol, maxIter);
}
#endif