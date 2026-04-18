#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "INTEGRACAO_NUMERICA/newtonCotes.hpp" 
#include "utils/function.hpp"


using namespace std;

int main(){
    cout << "INTEGRACAO NUMERICA - NEWTON COTES" << endl;
    cout << "-------------------------------------------------------------\n\n\n\n" << endl;

    cout << "TESTE 1: f(x) = x^4, a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    Function f1("x^4");
    cout << "GRAU 4 - FECHADA" << endl;
    cout << "Valor da integral: " << intgNewtonCotesFechada_grau4(f1, 0, 1) << endl;
    cout << "GRAU 4 - ABERTA" << endl;
    cout << "Valor da integral: " << intgNewtonCotesAberta_grau4(f1, 0, 1) << endl;
    cout << "Valor exato: " << 0.2 << endl << endl;

    cout << "TESTE 2: f(x) = x^2, a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    Function f2("x^2");
    cout << "GRAU 4 - FECHADA" << endl;
    cout << "Valor da integral: " << intgNewtonCotesFechada_grau4(f2, 0, 1) << endl;
    cout << "GRAU 4 - ABERTA" << endl;
    cout << "Valor da integral: " << intgNewtonCotesAberta_grau4(f2, 0, 1) << endl;
    cout << "Valor exato: " << 0.333333 << endl << endl;

    cout << "TESTE 3: f(x) = exp(x), a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    Function f3("exp(x)");
    cout << "GRAU 4 - FECHADA" << endl;
    cout << "Valor da integral: " << intgNewtonCotesFechada_grau4(f3, 0, 1) << endl;
    cout << "GRAU 4 - ABERTA" << endl;
    cout << "Valor da integral: " << intgNewtonCotesAberta_grau4(f3, 0, 1) << endl;
    cout << "Valor exato: " << 1.71828 << endl << endl;

    cout << "TESTE 4: f(x) = exp(x) + x, a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    Function f4("exp(x) + x");
    cout << "GRAU 4 - FECHADA" << endl;
    cout << "Valor da integral: " << intgNewtonCotesFechada_grau4(f4, 0, 1) << endl;
    cout << "GRAU 4 - ABERTA" << endl;
    cout << "Valor da integral: " << intgNewtonCotesAberta_grau4(f4, 0, 1) << endl;
    cout << "Valor exato: " << 2.21828 << endl << endl;


    return 0;

}


