#ifndef PARTICIONAMENTO_HPP
#define PARTICIONAMENTO_HPP
#include <iomanip>
#include <vector>
#include "../../utils/function.hpp"
#include "singularidade.hpp"

// =============================================================================
// SOLUCAO 2 (Aula #14): transforma o problema original com singularidade em
// x=a e/ou x=b (limites FINITOS) numa integral sobre toda a reta em s, via
// exponencial simples ou dupla (ja implementadas em singularidade.hpp), e
// resolve essa integral truncando-a em [-c, +c] com uma Gauss-Legendre
// composta (particionada em n subintervalos de 4 pontos cada).
//
//      I = int_a^b f(x) dx = int_{-inf}^{+inf} fbar(s) ds
//                          ~= int_{-c}^{+c} fbar(s) ds
//                          ~= soma sobre as n particoes de GL-4pontos
// =============================================================================


// -----------------------------------------------------------------------
// Composicao (particionamento) da integral em s no intervalo [-c, c],
// usando n subintervalos e Gauss-Legendre de 4 pontos em cada um --
// exatamente a mesma ideia de partitionGaussLegendre (particionamento.hpp),
// so que aqui o integrando ja e f_expSimples/f_expDupla.
// -----------------------------------------------------------------------

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


// -----------------------------------------------------------------------
// Para um valor de c fixo, refina n (dobrando) ate a variacao entre duas
// iteracoes sucessivas ficar abaixo da tolerancia "tol" (igual filosofia
// do testPartition() do particionamento.hpp).
// -----------------------------------------------------------------------

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


// -----------------------------------------------------------------------
// Roda a Solucao 2 completa para uma lista de valores crescentes de c,
// refinando n em cada um deles, e imprime uma tabela igual a Tabela 3
// da Aula 14 (c x I), tanto para exponencial simples quanto para dupla.
//
// CUIDADO (mencionado na Aula #14): na exponencial dupla aparecem termos
// tipo exp(exp(s)), que estouram a representacao de ponto flutuante para
// valores de c relativamente pequenos (tipicamente c > ~4). Por isso a
// lista de valores de c testados para a dupla deve ser bem mais modesta
// que a da simples.
// -----------------------------------------------------------------------

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