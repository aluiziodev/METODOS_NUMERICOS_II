#ifndef ALGEBRA_LINEAR_HPP
#define ALGEBRA_LINEAR_HPP
#include <vector>
#include <cmath>
using Matriz = std::vector<std::vector<double>>;
using Vetor  = std::vector<double>;

Vetor matVec(const Matriz &A, const Vetor &x) {
    int n = A.size();
    Vetor y(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            y[i] += A[i][j] * x[j];
    return y;
}

double norma(const Vetor &v) {
    double s = 0.0;
    for (double vi : v) s += vi*vi;
    return sqrt(s);
}

Vetor normalizar(const Vetor &v) {
    double n = norma(v);
    Vetor r(v.size());
    for (size_t i = 0; i < v.size(); i++) r[i] = v[i]/n;
    return r;
}

double produtoInterno(const Vetor &a, const Vetor &b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); i++) s += a[i]*b[i];
    return s;
}

Matriz identidade(int n) {
    Matriz I(n, Vetor(n, 0.0));
    for (int i = 0; i < n; i++) I[i][i] = 1.0;
    return I;
}

Matriz subtrairEscalarDiagonal(const Matriz &A, double shift) {
    Matriz B = A;
    for (size_t i = 0; i < A.size(); i++) B[i][i] -= shift;
    return B;
}

// Resolve B*y = x por eliminação de Gauss com pivotamento parcial
Vetor resolverSistema(Matriz B, Vetor x) {
    int n = B.size();
    for (int i = 0; i < n; i++) B[i].push_back(x[i]); // matriz aumentada

    for (int k = 0; k < n; k++) {
        int piv = k;
        for (int i = k+1; i < n; i++)
            if (fabs(B[i][k]) > fabs(B[piv][k])) piv = i;
        std::swap(B[k], B[piv]);

        for (int i = k+1; i < n; i++) {
            double m = B[i][k] / B[k][k];
            for (int j = k; j <= n; j++)
                B[i][j] -= m * B[k][j];
        }
    }

    Vetor y(n);
    for (int i = n-1; i >= 0; i--) {
        double soma = B[i][n];
        for (int j = i+1; j < n; j++) soma -= B[i][j]*y[j];
        y[i] = soma / B[i][i];
    }
    return y;
}
#endif