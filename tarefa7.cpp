#include <iostream>
#include <vector>
#include <math.h>
#include "utils/function.hpp"
#include "INTEGRACAO_NUMERICA/singularidade/particionamento.hpp"

using namespace std;


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
