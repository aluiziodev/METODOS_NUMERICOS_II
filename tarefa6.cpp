#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "INTEGRACAO_NUMERICA/gauss-/chebyshev.hpp"
#include "INTEGRACAO_NUMERICA/gauss-/hermite.hpp"
#include "INTEGRACAO_NUMERICA/gauss-/laguerre.hpp"
#include "utils/function.hpp"




int main () {
    Function f1("x^6");
    cout << "f(x) = x^6" << endl;

    cout << "\n\nResultado Gauss-Hermite 2 pontos: " << intgGaussHermite_2pontos(f1) << endl;
    cout << "Resultado Gauss-Hermite 3 pontos: " << intgGaussHermite_3pontos(f1) << endl;
    cout << "Resultado Gauss-Hermite 4 pontos: " << intgGaussHermite_4pontos(f1) << endl;

    cout << "\n\nResultado Gauss-Laguerre 2 pontos: " << intgGaussLaguerre_2pontos(f1) << endl;
    cout << "Resultado Gauss-Laguerre 3 pontos: " << intgGaussLaguerre_3pontos(f1) << endl;
    cout << "Resultado Gauss-Laguerre 4 pontos: " << intgGaussLaguerre_4pontos(f1) << endl;

    cout << "Resultado Gauss-Chebyshev 2 pontos: " << intgGaussChebyshev_2pontos(f1) << endl;
    cout << "Resultado Gauss-Chebyshev 3 pontos: " << intgGaussChebyshev_3pontos(f1) << endl;
    cout << "Resultado Gauss-Chebyshev 4 pontos: " << intgGaussChebyshev_4pontos(f1) << endl;


    return 0;
}