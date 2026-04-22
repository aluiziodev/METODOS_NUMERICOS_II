#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "INTEGRACAO_NUMERICA/newtonCotes/newtonCotes.hpp" 
#include "INTEGRACAO_NUMERICA/newtonCotes/particionamento.hpp" 
#include "utils/function.hpp"


using namespace std;

int main(){
    cout << "INTEGRACAO NUMERICA - NEWTON COTES" << endl;
    cout << "-------------------------------------------------------------\n\n\n\n" << endl;

    cout << "TESTE 1: f(x) = x^4, a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    Function f1("(sin(2*x) + 4*(x^2) + 3*x)^2");
    cout << "GRAU 4 - FECHADA" << endl;
    cout << "Valor da integral: " << endl;
    testFechada(partitionNewtonCotesFechada, f1, 0, 1, 1e-6, 3);
    cout << "GRAU 4 - ABERTA" << endl;
    cout << "Valor da integral: " << endl;
    testAberta(partitionNewtonCotesAberta, f1, 0, 1, 1e-6, 3);
    cout << "Valor exato: " << 17.8774 << endl << endl;

    /*cout << "TESTE 2: f(x) = x^2, a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    Function f2("x^2");
    cout << "GRAU 4 - FECHADA" << endl; 
    cout << "Valor da integral: " << endl;
    testFechada(partitionNewtonCotesFechada, f2, 0, 1, 1e-6);
    cout << "GRAU 4 - ABERTA" << endl;
    cout << "Valor da integral: " << endl;
    testAberta(partitionNewtonCotesAberta, f2, 0, 1, 1e-6);
    cout << "Valor exato: " << 0.333333 << endl << endl;

    cout << "TESTE 3: f(x) = exp(x), a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    Function f3("exp(x)");
    cout << "GRAU 4 - FECHADA" << endl;
    cout << "Valor da integral: " << endl;
    testFechada(partitionNewtonCotesFechada, f3, 0, 1, 1e-6);
    cout << "GRAU 4 - ABERTA" << endl;
    cout << "Valor da integral: " << endl;
    testAberta(partitionNewtonCotesAberta, f3, 0, 1, 1e-6);
    cout << "Valor exato: " << 1.71828 << endl << endl;

    cout << "TESTE 4: f(x) = exp(x) + x, a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    Function f4("exp(x) + x");
    cout << "GRAU 4 - FECHADA" << endl;
    cout << "Valor da integral: " << endl;
    testFechada(partitionNewtonCotesFechada, f4, 0, 1, 1e-6);
    cout << "GRAU 4 - ABERTA" << endl;
    cout << "Valor da integral: " << endl;
    testAberta(partitionNewtonCotesAberta, f4, 0, 1, 1e-6);
    cout << "Valor exato: " << 2.21828 << endl << endl;
    */

    return 0;

}


