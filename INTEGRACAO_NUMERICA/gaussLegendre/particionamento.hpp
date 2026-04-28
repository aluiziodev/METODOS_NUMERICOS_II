#ifndef PARTICIONAMENTO_HPP
#define PARTICIONAMENTO_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "../../utils/function.hpp"
#include "gaussLegendre.hpp"



double partitionGaussLegendre(Function &f, double a, double b, int n, int pontos=4){
    double h = (b - a) / n;
    double sum = 0.0;
    for(int i = 0; i < n; i++){
        double x0 = a + i*h;
        double x1 = a + (i+1)*h;
        sum += intgGaussLegendre(f, x0, x1, pontos);
    }
    return sum;
}

void testPartition(double (*metodo)(Function&, double, double, int, int),
            Function &f, double a, double b, double tol, int pontos=4, int maxIter = 10000) {
    int n = pontos; 
    double i = metodo(f, a, b, n, pontos);
    double erro = 1.0;
    int iter = 0;
    cout << "Teste da Regra de Gauss-Legendre" << endl;
    cout << "Iteracao: " << iter << ", n: " << n << ", Integral: " << fixed << 
            setprecision(6) << i << ", Erro: " << "N/A" << endl;
    while (erro > tol && iter < maxIter) {
        n *= 2;
        double j = metodo(f, a, b, n, pontos);
        erro = fabs(j - i);
        i = j;
        iter++;
        cout << "Iteracao: " << iter << ", n: " << n << ", Integral: " << fixed << 
        setprecision(6) << i << ", Erro: " << erro << endl;
        
        
    }

    cout << "Resultado: " << i << endl;
}



#endif