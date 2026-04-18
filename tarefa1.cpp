#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "DERIVADA_NUMERICA/derivadaNum.hpp" 
#include "utils/function.hpp"




using namespace std;



int main(){
    Function f("sqrt(exp(3*x) + 4*x^2)");
    double x = 2.0;
    double dlt = 0.5;
    double error = 1.0;
    int k = 1;
    cout << "        dlt(k)    |    f(x)    |    f''(x)    |    e(x)    " << endl;
    cout << "-------------------------------------------------------------" << endl;
    while(error>=0.00001){
        double fx = f.value(x);
        double dfx2 = dfNum2_Ordem4_Ct(f, x, dlt);
        error = abs((dfx2 - dfNum2_Ordem4_Ct(f, x, dlt*2))/dfx2);
        cout << "k="<< k << "|     " << dlt << "    |    " << fx << "    |    " << dfx2 << "    |    " << error << endl;
        dlt = dlt / 2;
        k++;
    }

    return 0;
}