#ifndef PARTICIONAMENTO_HPP
#define PARTICIONAMENTO_HPP
#include <iomanip>
#include <vector>
#include "../../utils/function.hpp"
#include "singularidade.hpp"



double partitionSolucao2_expSimples(Function &f, double a, double b, double c, int n){
    double s0 = -c, s1 = c;
    double h = (s1 - s0) / n;
    double soma = 0.0;
    for(int i = 0; i < n; i++){
        double sa = s0 + i*h;
        double sb = s0 + (i+1)*h;
        soma += intgGaussLegendre_4pontos_expSimples(f, sa, sb, a, b);
    }
    return soma;
}

double partitionSolucao2_expDupla(Function &f, double a, double b, double c, int n){
    double s0 = -c, s1 = c;
    double h = (s1 - s0) / n;
    double soma = 0.0;
    for(int i = 0; i < n; i++){
        double sa = s0 + i*h;
        double sb = s0 + (i+1)*h;
        soma += intgGaussLegendre_4pontos_expDupla(f, sa, sb, a, b);
    }
    return soma;
}



double refinaEmN_expSimples(Function &f, double a, double b, double c,
                             double tol, int nInicial = 4, int maxIter = 25){
    int n = nInicial;
    double anterior = partitionSolucao2_expSimples(f, a, b, c, n);
    double erro = 1.0;
    int iter = 0;
    while(erro > tol && iter < maxIter){
        n *= 2;
        double atual = partitionSolucao2_expSimples(f, a, b, c, n);
        erro = fabs(atual - anterior);
        anterior = atual;
        iter++;
    }
    return anterior;
}

double refinaEmN_expDupla(Function &f, double a, double b, double c,
                           double tol, int nInicial = 4, int maxIter = 25){
    int n = nInicial;
    double anterior = partitionSolucao2_expDupla(f, a, b, c, n);
    double erro = 1.0;
    int iter = 0;
    while(erro > tol && iter < maxIter){
        n *= 2;
        double atual = partitionSolucao2_expDupla(f, a, b, c, n);
        erro = fabs(atual - anterior);
        anterior = atual;
        iter++;
    }
    return anterior;
}



void testSolucao2(Function &f, double a, double b, double tol, vector<double> &cs_simples, 
    vector<double> &cs_dupla, double simetria = 1.0){

    
    cout << "Exponencial SIMPLES" << endl;
    cout << fixed << setprecision(8);
    for(double c : cs_simples){
        double I = simetria * refinaEmN_expSimples(f, a, b, c, tol);
        cout << "c = " << c << "  |  I = " << I << endl;
    }

    cout << "Exponencial DUPLA" << endl;
    for(double c : cs_dupla){
        double I = simetria * refinaEmN_expDupla(f, a, b, c, tol);
        cout << "c = " << c << "  |  I = " << I << endl;
    }
}


#endif 