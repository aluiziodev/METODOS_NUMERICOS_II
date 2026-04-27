#ifndef GAUSS_LEGENDRE_HPP
#define GAUSS_LEGENDRE_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../../utils/function.hpp"

//QUADRATURA DE GAUSS-LEGENDRE 2 PONTOS

double intgGaussLegendre_2pontos(Function &f, double a, double b){
    double x_a1 = (a+b)/2 - (b-a)/2 * sqrt(1.0/3.0);
    double x_a2 = (a+b)/2 + (b-a)/2 * sqrt(1.0/3.0);
    return ((b-a)/2) * (f.value(x_a1) + f.value(x_a2));
}

double intgGaussLegendre_3pontos(Function &f, double a, double b){
    double x_a1 = (a+b)/2 - (b-a)/2 * sqrt(3.0/5.0);
    double x_a2 = (a+b)/2;
    double x_a3 = (a+b)/2 + (b-a)/2 * sqrt(3.0/5.0);
    return ((b-a)/18) * (5*f.value(x_a1) + 8*f.value(x_a2) + 5*f.value(x_a3));
}


#endif
