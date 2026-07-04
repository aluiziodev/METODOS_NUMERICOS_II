#ifndef HOUSEHOLDER_QR_HPP
#define HOUSEHOLDER_QR_HPP
#include "../utils/algebraLinear.hpp"

// CAIXA PRETA 1: reduz matriz simétrica A a forma tridiagonal T = H^T A H
Matriz householderTridiagonalizacao(Matriz A) {
    int n = A.size();
    for (int k = 0; k < n-2; k++) {
        Vetor x(n-k-1);
        for (int i = k+1; i < n; i++) x[i-k-1] = A[i][k];

        double normaX = norma(x);
        double alpha = (x[0] >= 0) ? -normaX : normaX;

        Vetor v = x;
        v[0] -= alpha;
        v = normalizar(v);

        // H = I - 2*v*v^T (aplicado no bloco inferior-direito)
        Matriz H = identidade(n);
        for (int i = 0; i < (int)v.size(); i++)
            for (int j = 0; j < (int)v.size(); j++)
                H[k+1+i][k+1+j] -= 2.0 * v[i] * v[j];

        // A = H * A * H
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
    }
    return A;
}

// Um passo de decomposição QR via rotações de Givens (eficiente p/ tridiagonal)
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
            // (Q acumulado transposto; sinais ajustados para Q = produto de G^T)
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

// CAIXA PRETA 2: itera QR sobre a matriz tridiagonal até convergir
// para uma matriz (quase) diagonal, cuja diagonal são os autovalores
struct ResultadoQR {
    Vetor autovalores;
    int iteracoes;
};

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

// Pipeline completo: Householder -> QR (as "duas caixas pretas em sequência")
ResultadoQR autovaloresSimetrica(const Matriz &A, double tol = 1e-10, int maxIter = 500) {
    Matriz T = householderTridiagonalizacao(A);
    return iteracaoQR(T, tol, maxIter);
}
#endif