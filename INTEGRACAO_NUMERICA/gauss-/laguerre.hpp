#ifndef LAGUERRE_HPP
#define LAGUERRE_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../../utils/function.hpp"

// GAUSS-LAGUERRE 2 PONTOS

double intgGaussLaguerre_2pontos(Function &f){
    double x_a1 = 2 - sqrt(2.0);
    double x_a2 = 2 + sqrt(2.0);
    double w_a1 = (1/4.0) * x_a2;
    double w_a2 = (1/4.0) * x_a1;
    return f.value(x_a1)*w_a1 + f.value(x_a2)*w_a2;
}

// GAUSS-LAGUERRE 3 PONTOS

double intgGaussLaguerre_3pontos(Function &f){
    double x_a1 = 0.4157745568;
    double x_a2 = 2.2942803603;
    double x_a3 = 6.2899450829;
    double w_a1 = 0.7110930099;
    double w_a2 = 0.2785177336;
    double w_a3 = 0.0103892565;
    return f.value(x_a1)*w_a1 + f.value(x_a2)*w_a2 + f.value(x_a3)*w_a3;
}

double intgGaussLaguerre_4pontos(Function &f){
    double x_a1 = 0.3225476896;
    double x_a2 = 1.7457611012;
    double x_a3 = 4.5366202969;
    double x_a4 = 9.3950709123;
    double w_a1 = 0.6031541043;
    double w_a2 = 0.3574186924;
    double w_a3 = 0.0388879085;
    double w_a4 = 0.0005392947;
    return f.value(x_a1)*w_a1 + f.value(x_a2)*w_a2 + f.value(x_a3)*w_a3 + f.value(x_a4)*w_a4;
}







#endif