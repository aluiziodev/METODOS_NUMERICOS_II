#include <cmath>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include "DERIVADA_NUMERICA/derivadaNum.hpp" 
#include "utils/function.hpp"

using namespace std;


int main() {
    double x = 1.0;
    double dlt = 0.5;
    Function f("sin(x) + x*x");
    double d_exact = cos(x) + 2*x;
    
    cout << fixed << setprecision(10);
    cout << "Funcao: f(x) = sin(x) + x^2" << endl;
    cout << "Derivada exata em x = " << x << " : " << d_exact << endl << endl;
    cout << "Iteracao | h | Forward | ErroF | Central | ErroC | Backward | ErroB " << endl;
    cout << "--------------------------------------------------------------------------" << endl;

    for (int i = 0; i < 3; i++) {
            cout << i + 1
             <<" | "  << dlt
             <<" | "  << dfNum1_FW(f, x, dlt)
             <<" | "  << fabs(dfNum1_FW(f, x, dlt) - d_exact)
             <<" | "  << dfNum1_CT(f, x, dlt)
             <<" | "  << fabs(dfNum1_CT(f, x, dlt) - d_exact)
             <<" | "  << dfNum1_BW(f, x, dlt)
             <<" | "  << fabs(dfNum1_BW(f, x, dlt) - d_exact) << endl;
        dlt /= 10.0;
    }

    return 0;
}
