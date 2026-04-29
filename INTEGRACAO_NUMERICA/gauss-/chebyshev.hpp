#ifndef CHEBYSHEV_HPP
#define CHEBYSHEV_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../../utils/function.hpp"


// GAUSS-CHEBYSHEV 2 PONTOS

double intgGaussChebyshev_2pontos(Function &f){
    double x_a1 = -1.0/sqrt(2.0);
    double x_a2 = 1.0/sqrt(2.0);
    double w_a = M_PI/2.0;
    return f.value(x_a1)*w_a + f.value(x_a2)*w_a;
}

// GAUSS-CHEBYSHEV 3 PONTOS

double intgGaussChebyshev_3pontos(Function &f){
    double x_a1 = -sqrt(3.0)/2.0;
    double x_a2 = 0.0;
    double x_a3 = sqrt(3.0)/2.0;
    double w_a = M_PI/3.0;
    return f.value(x_a1)*w_a + f.value(x_a2)*w_a + f.value(x_a3)*w_a;
}

// GAUSS-CHEBYSHEV 4 PONTOS

double intgGaussChebyshev_4pontos(Function &f){
    double x_a1 = -sqrt(2.0 + sqrt(2.0))/2.0;
    double x_a2 = -sqrt(2.0 - sqrt(2.0))/2.0;
    double x_a3 = sqrt(2.0 - sqrt(2.0))/2.0;
    double x_a4 = sqrt(2.0 + sqrt(2.0))/2.0;
    double w_a = M_PI/4.0;
    return f.value(x_a1)*w_a + f.value(x_a2)*w_a + f.value(x_a3)*w_a + f.value(x_a4)*w_a;
}


#endif