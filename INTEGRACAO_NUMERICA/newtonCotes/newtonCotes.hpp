#ifndef NEWTON_COTES_HPP
#define NEWTON_COTES_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../../utils/function.hpp"


//FORMULAS FECHADAS

double intgNewtonCotesFechada_grau1(Function &f, double a, double b){
    return (b-a)/2.0 * (f.value(a) + f.value(b));
}

double intgNewtonCotesFechada_grau2(Function &f, double a, double b){
    double h = (b-a)/2.0;
    double x0 = a;
    double x1 = a + h;
    double x2 = b;
    return (h/3.0) * (f.value(x0) + 4*f.value(x1) + f.value(x2));
}

double intgNewtonCotesFechada_grau3(Function &f, double a, double b){
    double h = (b-a)/3.0;
    double x0 = a;
    double x1 = a + h;
    double x2 = a + 2*h;
    double x3 = b;
    return (3*h/8.0) * (f.value(x0) + 3*f.value(x1) + 3*f.value(x2) + f.value(x3));
}


double intgNewtonCotesFechada_grau4(Function &f, double a, double b){
    double h = (b - a) / 4.0;

    double x0 = a;
    double x1 = a + h;
    double x2 = a + 2*h;
    double x3 = a + 3*h;
    double x4 = b;
    return (2*h/45.0) * (
        7*f.value(x0) +
        32*f.value(x1) +
        12*f.value(x2) +
        32*f.value(x3) +
        7*f.value(x4)
    );

}

double intgNewtonCotesFechada(Function &f, double a, double b,int grau=4){
    switch (grau){
        case 1:
            return intgNewtonCotesFechada_grau1(f, a, b);
        case 2:
            return intgNewtonCotesFechada_grau2(f, a, b);
        case 3:
            return intgNewtonCotesFechada_grau3(f, a, b);
        case 4:
            return intgNewtonCotesFechada_grau4(f, a, b);
        default:
            cout << "Grau inválido para a fórmula fechada." << endl;
            return -1.0;
    }
}


//FORMULAS ABERTAS

double intgNewtonCotesAberta_grau1(Function &f, double a, double b){
    double h = (b - a) / 3.0;
    double x0 = a + h;
    double x1 = a + 2*h;
    return (b-a)/2.0 * (f.value(x0) + f.value(x1));
}

double intgNewtonCotesAberta_grau2(Function &f, double a, double b){
    double h = (b - a) / 4.0;
    double x0 = a + h;
    double x1 = a + 2*h;
    double x2 = a + 3*h;
    return (4*h/3.0) * (2*f.value(x0) - f.value(x1) + 2*f.value(x2));
}

double intgNewtonCotesAberta_grau3(Function &f, double a, double b){
    double h = (b - a) / 5.0;
    double x0 = a + h;
    double x1 = a + 2*h;
    double x2 = a + 3*h;
    double x3 = a + 4*h;
    return (5*h/24.0) * (
        11*f.value(x0) + 
        f.value(x1) + 
        f.value(x2) + 
        11*f.value(x3)
    );
}


double intgNewtonCotesAberta_grau4(Function &f, double a, double b){
    double h = (b - a) / 6.0;

    double x0 = a + h;
    double x1 = a + 2*h;
    double x2 = a + 3*h;
    double x3 = a + 4*h;
    double x4 = a + 5*h;
    return (3*h/10.0) * (
        11*f.value(x0) -
        14*f.value(x1) +
        26*f.value(x2) -
        14*f.value(x3) +
        11*f.value(x4)
    );

}

double intgNewtonCotesAberta(Function &f, double a, double b,int grau=4){
    switch (grau){
        case 1:
            return intgNewtonCotesAberta_grau1(f, a, b);
        case 2:
            return intgNewtonCotesAberta_grau2(f, a, b);
        case 3:
            return intgNewtonCotesAberta_grau3(f, a, b);
        case 4:
            return intgNewtonCotesAberta_grau4(f, a, b);
        default:
            cout << "Grau inválido para a fórmula aberta." << endl;
            return -1.0;
    }
}



#endif