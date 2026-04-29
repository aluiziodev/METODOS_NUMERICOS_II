#ifndef HERMITE_HPP
#define HERMITE_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../../utils/function.hpp"

// GAUSS-HERMITE 2 PONTOS

double intgGaussHermite_2pontos(Function &f){
    double x_a1 = -1.0/sqrt(2.0);
    double x_a2 = 1.0/sqrt(2.0);
    double w_a = sqrt(M_PI)/2.0;
    return f.value(x_a1)*w_a + f.value(x_a2)*w_a;
}

// GAUSS-HERMITE 3 PONTOS

double intgGaussHermite_3pontos(Function &f){
    double x_a1 = -sqrt(3.0)/sqrt(2.0);
    double x_a2 = 0.0;
    double x_a3 = sqrt(3.0)/sqrt(2.0);
    double w_a = sqrt(M_PI)/6.0;
    double w_b = 2.0*sqrt(M_PI)/3.0;
    return f.value(x_a1)*w_a + f.value(x_a2)*w_b + f.value(x_a3)*w_a;
}

// GAUSS-HERMITE 4 PONTOS


double intgGaussHermite_4pontos(Function &f){
    double x_a1 = -sqrt(2.0*sqrt(6.0) + 6)/2.0;
    double x_a2 = -sqrt(-2.0*sqrt(6.0) + 6)/2.0;
    double x_a3 = sqrt(-2.0*sqrt(6.0) + 6)/2.0;
    double x_a4 = sqrt(2.0*sqrt(6.0) + 6)/2.0;;
    double w_a = sqrt(M_PI)*(3.0 - sqrt(6.0)) / 12;
    double w_b = sqrt(M_PI)*(3.0 + sqrt(6.0)) / 12;
    return f.value(x_a1)*w_a + f.value(x_a2)*w_b + f.value(x_a3)*w_b + f.value(x_a4)*w_a;
}


#endif