#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "INTEGRACAO_NUMERICA/gaussLegendre/particionamento.hpp"
#include "utils/function.hpp"

using namespace std;

int main(){
    Function f1("(sin(2*x) + 4*(x^2) + 3*x)^2");
    double a = 0;
    double b = 1;
    cout << "Resultado 2 pontos: " << intgGaussLegendre_2pontos(f1, a, b) << endl;
    cout << "Resultado 3 pontos: " << intgGaussLegendre_3pontos(f1, a, b) << endl;
    cout << "Resultado 4 pontos: " << intgGaussLegendre_4pontos(f1, a, b) << endl;


    cout << "INTEGRACAO NUMERICA - GAUSS-LEGENDRE" << endl;
    cout << "-------------------------------------------------------------\n\n\n\n" << endl;
    cout << "TESTE 1: f(x) = (sin(2*x) + 4*(x^2) + 3*x)^2, a = 0, b = 1" << endl;
    cout << "-------------------------------------------------------------" << endl;
    cout << "2 PONTOS" << endl;
    cout << "Valor da integral: " << endl;
    testPartition(partitionGaussLegendre, f1, 0, 1, 0.000001, 2);
    cout << "3 PONTOS" << endl;
    cout << "Valor da integral: " << endl;
    testPartition(partitionGaussLegendre, f1, 0, 1, 0.000001, 3);
    cout << "4 PONTOS" << endl;
    cout << "Valor da integral: " << endl;
    testPartition(partitionGaussLegendre, f1, 0, 1, 0.000001, 4);
    return 0;
}


