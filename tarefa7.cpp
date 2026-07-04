#include <iostream>
#include <vector>
#include <math.h>
#include "utils/function.hpp"
#include "INTEGRACAO_NUMERICA/singularidade/particionamento.hpp"

using namespace std;

// -----------------------------------------------------------------------
// PROBLEMA 1:  I = int_{-1}^{1}  dx / cbrt(x^2)   = 6
//   f(x) = 1 / (x^2)^(1/3)
//
//   ATENCAO: a singularidade aqui esta em x = 0, que e um ponto INTERIOR
//   de [-1,1] (nao um dos extremos). As mudancas de variavel por
//   exponencial simples/dupla foram deduzidas na Aula #13/#14 para
//   singularidades NOS EXTREMOS a ou b -- se aplicadas diretamente em
//   [-1,1], a singularidade cai exatamente no centro do dominio em s
//   (s=0), e a convergencia da quadratura composta fica muito lenta
//   (a integral ainda converge, mas exige n absurdamente grande para
//   uma tolerancia de 1e-8, de forma parecida com o que a Aula #13
//   mostrou ao tentar Newton-Cotes Fechada num ponto de singularidade
//   nao removido).
//
//   Solucao: como f(x) = x^(-2/3) e uma funcao PAR, usa-se a simetria
//            I = int_{-1}^{1} f(x) dx = 2 * int_{0}^{1} f(x) dx
//   transformando o problema em uma integral em [0, 1] com singularidade
//   NO EXTREMO x = 0 -- exatamente o mesmo tipo de problema tratado nos
//   exemplos da Aula #13/#14 (lembra o exemplo 1/sqrt(x) em [0,1]).
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// PROBLEMA 2:  I = int_{-2}^{0}  dx / sqrt(4 - x^2)  = pi/2
//   f(x) = 1 / sqrt(4 - x^2)
//   Singularidade no limite inferior x = -2 (caso classico de singularidade
//   em ponto extremo do intervalo, igual aos exemplos da Aula #13/#14).
// -----------------------------------------------------------------------


int main(){

    double tol = 1e-8;

    // PROBLEMA 1    
    cout << "-----------------------------------------------------" << endl;
    cout << " PROBLEMA 1: I = int_{-1}^{1} dx / cbrt(x^2)  (exato = 6)" << endl;
    cout << " (usando simetria: I = 2 * int_{0}^{1} dx / cbrt(x^2))" << endl;
    cout << "-----------------------------------------------------" << endl;

    Function funcao1("1 / root(x*x, 3)");
    double a1 = 0.0, b1 = 1.0;  
    vector<double> cs_simples_1 = {2, 4, 6, 8, 10, 12};
    vector<double> cs_dupla_1   = {1.5, 2, 2.5, 3, 3.5};
    double simetria = 2.0; // f(x) = 1 / cbrt(x^2) e uma funcao PAR, entao I = 2 * int_{0}^{1} f(x) dx

    testSolucao2(funcao1, a1, b1, tol, cs_simples_1, cs_dupla_1, simetria);

    //  PROBLEMA 2 
    cout << "-----------------------------------------------------" << endl;
    cout << " PROBLEMA 2: I = int_{-2}^{0} dx / sqrt(4 - x*x)  (exato = pi/2 = "
         << (M_PI/2.0) << ")" << endl;
    cout << "-----------------------------------------------------" << endl;

    Function funcao2("1 / sqrt(4 - (x*x))");
    double a2 = -2.0, b2 = 0.0;

    vector<double> cs_simples_2 = {2, 4, 6, 8, 10, 12};
    vector<double> cs_dupla_2   = {1.5, 2, 2.5, 3, 3.5};

    testSolucao2(funcao2, a2, b2, tol, cs_simples_2, cs_dupla_2);

    return 0;
}
