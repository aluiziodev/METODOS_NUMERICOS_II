#ifndef PARTICIONAMENTO_HPP
#define PARTICIONAMENTO_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "../../utils/function.hpp"
#include "newtonCotes.hpp"



double partitionNewtonCotesFechada(Function &f, double a, double b, int n, int grau){
    double h = (b - a) / n;
    double sum = 0.0;
    for(int i = 0; i < n; i+=grau){
        double x0 = a + i*h;
        double x1 = a + (i+grau)*h;
        sum += intgNewtonCotesFechada(f, x0, x1, grau);
    }
    return sum;
}


double partitionNewtonCotesAberta(Function &f, double a, double b, int n, int grau){
    double h = (b - a) / n;
    double sum = 0.0;
    for(int i = 0; i < n; i+=(grau+1)){
        double x0 = a + i*h;
        double x1 = a + (i+grau+1)*h;
        sum += intgNewtonCotesAberta(f, x0, x1, grau);
    }
    return sum;
}


void testFechada(double (*metodo)(Function&, double, double, int, int),
            Function &f, double a, double b, double tol, int grau=4, int maxIter = 10000) {
    int n = grau; 
    double i = metodo(f, a, b, n, grau);
    double erro = 1.0;
    int iter = 0;
    cout << "Teste da Regra de Newton-Cotes Fechada" << endl;
    while (erro > tol && iter < maxIter) {
        cout << "Iteracao: " << iter << ", n: " << n << ", Integral: " << fixed << 
        setprecision(6) << i << ", Erro: " << erro << endl;
        n *= 2;
        double j = metodo(f, a, b, n, grau);
        erro = fabs(j - i);
        i = j;
        iter++;
    }

    cout << "Resultado: " << i << endl;
}


void testAberta(double (*metodo)(Function&, double, double, int, int),
            Function &f, double a, double b, double tol, int grau=4, int maxIter = 10000) {
    int n = grau+1; 
    double i = metodo(f, a, b, n, grau);
    double erro = 1.0;
    int iter = 0;
    cout << "Teste da Regra de Newton-Cotes Aberta" << endl;
    while (erro > tol && iter < maxIter) {
        cout << "Iteracao: " << iter << ", n: " << n << ", Integral: " << fixed << 
        setprecision(6) << i << ", Erro: " << erro << endl;
        n *= 2;
        double j = metodo(f, a, b, n, grau);

        erro = fabs(j - i);
        i = j;
        iter++;
    }

    cout << "Resultado: " << i << endl;
}

#endif