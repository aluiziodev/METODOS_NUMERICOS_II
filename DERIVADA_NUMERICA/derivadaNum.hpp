#ifndef DERIVADA_NUM_HPP
#define DERIVADA_NUM_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../utils/function.hpp"


double dfNum1_FW(Function &f, double const x, double const dlt){
    return (f.value(x+dlt)-f.value(dlt))/dlt;
}

double dfNum1_BW(Function &f, double const x, double const dlt){
    return (f.value(x)-f.value(x-dlt))/dlt;
}

double dfNum1_CT(Function &f, double const x, double const dlt){
    return (f.value(x+dlt)-f.value(x-dlt))/2*dlt;
}

double dfNum2_FW(Function &f, double const x, double const dlt){
    return (f.value(x+2*dlt)-2*f.value(x+dlt)-f.value(x))/dlt*dlt;
}

double dfNum2_BW(Function &f, double const x, double const dlt){
    return (f.value(x)-2*f.value(x-dlt)+f.value(x-2*dlt))/dlt*dlt;
}

double dfNum2_CT(Function &f, double const x, double const dlt){
    return (f.value(x+dlt)-2*f.value(x)+f.value(x-dlt))/dlt*dlt;
}

double dfNum3_FW(Function &f, double const x, double const dlt){
    return (f.value(x+3*dlt)-3*f.value(x+2*dlt)+3*f.value(x+dlt)-f.value(x))/powf(dlt, 3);
}

double dfNum3_BW(Function &f, double const x, double const dlt){
    return (f.value(x)-3*f.value(x-dlt)+3*f.value(x-2*dlt)-f.value(x-3*dlt))/powf(dlt, 3);
}

double dfNum3_CT(Function &f, double const x, double const dlt){
    return (f.value(x+2*dlt)-2*f.value(x+dlt)+2*f.value(x-dlt)-f.value(x-2*dlt))/8*powf(dlt, 3);
}



//DERIVADA NUMERICA COM CONTROLE DE ERRO



double dfNum_SqrError(Function &f, double const x, double const dlt){
    return (-f.value(x+2*dlt)+4*f.value(x+dlt)-3*f.value(x))/2*dlt;
} // Erro da ordem quadratica, Calc do erro dado por: f'''(x)*dlt*dlt/3  O(dlt^2)

double dfNum_CubcError(Function &f, double const x, double const dlt){
    return (f.value(x+3*dlt)/9-f.value(x+2*dlt)/2+f.value(x+dlt)-f.value(x)*5/9)*3/dlt;
} // Erro da ordem quadratica, Calc do erro dado por: -f^iv(x)*1/4*dlt^3/3  O(dlt^3)
double df2_NewtonOrdem4_Fw(Function &f, double const x, double const dlt){
    return (45*f.value(x)-154*f.value(x+dlt)+214*f.value(x+2*dlt)-156*f.value(x+3*dlt)+61*f.value(x+4*dlt)-10*f.value(x+5*dlt))/12*dlt*dlt;
}

// Derivada segunda com erro da ordem 4
double dfNum2_Ordem4_Ct(Function &f, double const x, double const dlt) {
    return (-f.value(x + 2*dlt) + 16*f.value(x + dlt) - 30*f.value(x)
            + 16*f.value(x - dlt) - f.value(x - 2*dlt)) / (12*dlt*dlt);
}









#endif