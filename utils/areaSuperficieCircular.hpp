#ifndef AREA_SUPERFICIE_CIRCULAR_HPP
#define AREA_SUPERFICIE_CIRCULAR_HPP
#include <cmath>
#include "function2D.hpp"
 
// =============================================================================
// Aula #15 - Problema 2:
// Área da superfície z = f(x,y) sobre a região elíptica/circular
//      U = { (x,y) | x^2/a^2 + y^2/b^2 <= 1 }
//
// Reaproveita a classe Function2D (já usada em integralDupla.hpp) para
// f(x,y), calculando as derivadas parciais por diferenças finitas centrais
// (pois Function2D só expõe value(x,y), sem derivada analítica).
//
// PASSO 1) Mudança de variável 1 (igual ao caso da elipse, eq. 24):
//      x(alfa,beta) = a*alfa*cos(beta)
//      y(alfa,beta) = b*alfa*sin(beta)      , alfa em [0,1], beta em [0,2pi]
//      |J1| = a*b*alfa   =>   dA = a*b*alfa dalfa dbeta
//
// PASSO 2) Mudança de variável 2 (igual à solução do Problema 1, eq. 34):
//      alfa(xi)  = 1/2 + 1/2 * xi      , xi  em [-1,1]   =>  dalfa = (1/2) dxi
//      beta(eta) = pi  + pi  * eta     , eta em [-1,1]   =>  dbeta = pi   deta
//      |J2| = pi/2
//
// PASSO 3) Quadratura de Gauss-Legendre com 3 pontos em cada direção.
// =============================================================================
 
// --- derivadas parciais de f(x,y) por diferença finita central ---
double dfdx(Function2D &f, double x, double y, double h = 1e-6) {
    return (f.value(x + h, y) - f.value(x - h, y)) / (2.0 * h);
}
 
double dfdy(Function2D &f, double x, double y, double h = 1e-6) {
    return (f.value(x, y + h) - f.value(x, y - h)) / (2.0 * h);
}
 
// --- área da superfície z=f(x,y) sobre a região elíptica x^2/a^2+y^2/b^2<=1 ---
double areaSuperficieRegiaoEliptica(Function2D &f, double a, double b) {
 
    // Nós e pesos de Gauss-Legendre com 3 pontos em [-1,1]
    double nodes[3]   = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
    double weights[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
 
    double soma = 0.0;
 
    for (int i = 0; i < 3; i++) {
        double eta  = nodes[i];
        double beta = M_PI + M_PI * eta;                 // passo 2 (beta)
 
        double somaInterna = 0.0;
        for (int j = 0; j < 3; j++) {
            double xi    = nodes[j];
            double alfa  = 0.5 + 0.5 * xi;                // passo 2 (alfa)
 
            double x = a * alfa * cos(beta);              // passo 1
            double y = b * alfa * sin(beta);              // passo 1
 
            double fx = dfdx(f, x, y);
            double fy = dfdy(f, x, y);
 
            double integrando = sqrt(fx*fx + fy*fy + 1.0) * (a * b * alfa); // |J1|
 
            somaInterna += weights[j] * integrando;
        }
        soma += weights[i] * somaInterna;
    }
 
    double jacobiano2 = 0.5 * M_PI;   // |J2| = (1/2)*(pi)
    return soma * jacobiano2;
}
 
// Caso particular: região circular de raio R (a = b = R)
double areaSuperficieRegiaoCircular(Function2D &f, double R) {
    return areaSuperficieRegiaoEliptica(f, R, R);
}
 
#endif