#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "INTEGRACAO_NUMERICA/gaussLegendre/gaussLegendre.hpp" 
#include "utils/function.hpp"

using namespace std;

int main(){
    Function f1("(sin(2*x) + 4*(x^2) + 3*x)^2");
    double a = 0;
    double b = 1;
    cout << "Resultado 2 pontos: " << intgGaussLegendre_2pontos(f1, a, b) << endl;
    cout << "Resultado 3 pontos: " << intgGaussLegendre_3pontos(f1, a, b) << endl;
    return 0;
}