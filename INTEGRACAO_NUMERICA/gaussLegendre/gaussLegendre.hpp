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


//QUADRATURA DE GAUSS-LEGENDRE 3 PONTOS

double intgGaussLegendre_3pontos(Function &f, double a, double b){
    double x_a1 = (a+b)/2 - (b-a)/2 * sqrt(3.0/5.0);
    double x_a2 = (a+b)/2;
    double x_a3 = (a+b)/2 + (b-a)/2 * sqrt(3.0/5.0);
    return ((b-a)/18) * (5*f.value(x_a1) + 8*f.value(x_a2) + 5*f.value(x_a3));
}


//QUADRATURA DE GAUSS-LEGENDRE 4 PONTOS

double intgGaussLegendre_4pontos(Function &f, double a, double b){
    double sqrt65 = sqrt(6.0 / 5.0);
    double x_a1 = sqrt((3.0 - 2.0*sqrt65) / 7.0);
    double x_a2 = sqrt((3.0 + 2.0*sqrt65) / 7.0);
    double x_a3 = -sqrt((3.0 - 2.0*sqrt65) / 7.0);
    double x_a4 = -sqrt((3.0 + 2.0*sqrt65) / 7.0);
    double w_a1 = (18.0 + sqrt(30.0)) / 36.0;
    double w_a2 = (18.0 - sqrt(30.0)) / 36.0;
    double mid  = (a + b) / 2.0;
    double half = (b - a) / 2.0;
    return half * (w_a1*f.value(mid+half*x_a1) + 
                w_a2*f.value(mid+half*x_a2) + 
                w_a1*f.value(mid+half*x_a3) + 
                w_a2*f.value(mid+half*x_a4));
}

double intgGaussLegendre(Function &f, double a, double b, int pontos=4){
    switch (pontos){
        case 2:
            return intgGaussLegendre_2pontos(f, a, b);
        case 3:
            return intgGaussLegendre_3pontos(f, a, b);
        case 4:
            return intgGaussLegendre_4pontos(f, a, b);
        default:
            cout << "Numero de pontos invalido. Usando 4 pontos por padrao." << endl;
            return intgGaussLegendre_4pontos(f, a, b);
    }
}


#endif
