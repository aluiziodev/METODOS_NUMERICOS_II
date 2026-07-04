#ifndef INTEGRAL_DUPLA_HPP
#define INTEGRAL_DUPLA_HPP
#include "../../utils/function2D.hpp"
#include "gaussLegendre.hpp"

// Integral dupla de f(x,y) sobre o retângulo [ax,bx] x [ay,by],
// usando Gauss-Legendre de 4 pontos em cada dimensão (aninhado).
double intgDuplaGaussLegendre(Function2D &f, double ax, double bx,
                               double ay, double by) {
    double sqrt65 = sqrt(6.0 / 5.0);
    double nodes[4] = {
        sqrt((3.0 - 2.0*sqrt65) / 7.0),
        sqrt((3.0 + 2.0*sqrt65) / 7.0),
       -sqrt((3.0 - 2.0*sqrt65) / 7.0),
       -sqrt((3.0 + 2.0*sqrt65) / 7.0)
    };
    double weights[4] = {
        (18.0 + sqrt(30.0)) / 36.0,
        (18.0 - sqrt(30.0)) / 36.0,
        (18.0 + sqrt(30.0)) / 36.0,
        (18.0 - sqrt(30.0)) / 36.0
    };

    double midX = (ax+bx)/2.0, halfX = (bx-ax)/2.0;
    double midY = (ay+by)/2.0, halfY = (by-ay)/2.0;

    double total = 0.0;
    for (int i = 0; i < 4; i++) {
        double y = midY + halfY * nodes[i];
        double innerSum = 0.0;
        for (int j = 0; j < 4; j++) {
            double x = midX + halfX * nodes[j];
            innerSum += weights[j] * f.value(x, y);
        }
        total += weights[i] * (halfX * innerSum);
    }
    return halfY * total;
}

// Área de uma região retangular (f(x,y) = 1)
double areaRetangulo(double ax, double bx, double ay, double by) {
    Function2D um("1");
    return intgDuplaGaussLegendre(um, ax, bx, ay, by);
}

// Volume sob a superfície z = f(x,y), sobre [ax,bx] x [ay,by]
double volumeSobSuperficie(Function2D &f, double ax, double bx,
                            double ay, double by) {
    return intgDuplaGaussLegendre(f, ax, bx, ay, by);
}
#endif