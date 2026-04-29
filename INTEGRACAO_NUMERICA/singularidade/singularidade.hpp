#ifndef SINGULARIDADE_HPP
#define SINGULARIDADE_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../../utils/function.hpp"


// MUDANÇA DE VARIAVEL PARA INTEGRAL COM SINGULARIDADE

double x_expSimples(double s, double a, double b){
    return ((a+b)/2) + ((b-a)/2)*tanh(s);
}

double dxds_expSimples(double s, double a, double b){
    return ((b-a)/2)*(1/pow(cosh(s),2));
}

double f_expSimples(Function &f, double s, double a, double b){
    return f.value(x_expSimples(s, a, b)) * dxds_expSimples(s, a, b);
}


double x_expDupla(double s, double a, double b){
    return ((a+b)/2) + ((b-a)/2)*tanh((M_PI/2)*sinh(s));
}

double dxds_expDupla(double s, double a, double b){
    return ((b-a)/2)*(M_PI*cosh(s)/(2*pow(cosh((M_PI/2)*sinh(s)),2)));
}

double f_expDupla(Function &f, double s, double a, double b){
    return f.value(x_expDupla(s, a, b)) * dxds_expDupla(s, a, b);
}


double intgGaussLegendre_4pontos_expSimples(Function &f,
                                    double s0, double s1,
                                    double a, double b){
    double sqrt65 = sqrt(6.0 / 5.0);

    double x1 =  sqrt((3.0 - 2.0*sqrt65) / 7.0);
    double x2 =  sqrt((3.0 + 2.0*sqrt65) / 7.0);
    double x3 = -x1;
    double x4 = -x2;

    double w1 = (18.0 + sqrt(30.0)) / 36.0;
    double w2 = (18.0 - sqrt(30.0)) / 36.0;

    double mid  = (s0 + s1)/2.0;
    double half = (s1 - s0)/2.0;

    auto F = [&](double s){
        return f_expSimples(f, s, a, b);
    };

    return half * (
        w1*F(mid + half*x1) +
        w2*F(mid + half*x2) +
        w1*F(mid + half*x3) +
        w2*F(mid + half*x4)
    );
}





#endif