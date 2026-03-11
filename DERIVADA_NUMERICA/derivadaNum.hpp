#ifndef DERIVADA_NUM_HPP
#define DERIVADA_NUM_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>

float f(float const x){
    return sin(x)*powf(x, 3)*5;
}

float dfNum1_FW(float const x, float const dlt){
    return (f(x+dlt)-f(dlt))/dlt;
}

float dfNum1_BW(float const x, float const dlt){
    return (f(x)-f(x-dlt))/dlt;
}

float dfNum1_CT(float const x, float const dlt){
    return (f(x+dlt)-f(x-dlt))/2*dlt;
}

float dfNum2_FW(float const x, float const dlt){
    return (f(x+2*dlt)-2*f(x+dlt)-f(x))/dlt*dlt;
}

float dfNum2_BW(float const x, float const dlt){
    return (f(x)-2*f(x-dlt)+f(x-2*dlt))/dlt*dlt;
}

float dfNum2_CT(float const x, float const dlt){
    return (f(x+dlt)-2*f(x)+f(x-dlt))/dlt*dlt;
}

float dfNum3_FW(float const x, float const dlt){
    return (f(x+3*dlt)-3*f(x+2*dlt)+3*f(x+dlt)-f(x))/powf(dlt, 3);
}

float dfNum3_BW(float const x, float const dlt){
    return (f(x)-3*f(x-dlt)+3*f(x-2*dlt)-f(x-3*dlt))/powf(dlt, 3);
}

float dfNum3_CT(float const x, float const dlt){
    return (f(x+2*dlt)-2*f(x+dlt)+2*f(x-dlt)-f(x-2*dlt))/8*powf(dlt, 3);
}

//DERIVADA NUMERICA COM CONTROLE DE ERRO



float dfNum_SqrError(float const x, float const dlt){
    return (-f(x+2*dlt)+4*f(x+dlt)-3*f(x))/2*dlt;
} // Erro da ordem quadratica, Calc do erro dado por: f'''(x)*dlt*dlt/3  O(dlt^2)

float dfNum_CubcError(float const x, float const dlt){
    return (f(x+3*dlt)/9-f(x+2*dlt)/2+f(x+dlt)-f(x)*5/9)*3/dlt;
} // Erro da ordem quadratica, Calc do erro dado por: -f^iv(x)*1/4*dlt^3/3  O(dlt^3)


float df2_NewtonOrdem4_Fw(float const x, float const dlt){
    return (45*x-154*(x+dlt)+214*(x+2*dlt)-156*(x+3*dlt)+61*(x+4*dlt)-10*(x+5*dlt))/12*dlt*dlt;
}




/*float dfNumN_FW(float const x, float const dlt, int const n){
    if (n==0) return f(x);
    else return (dfNumN_FW(x+dlt, dlt, n-1)-dfNumN_FW(x, dlt, n-1))/dlt;
}

float dfNumN_BW(float const x, float const dlt, int const n){
    if (n==0) return f(x);
    else return (dfNumN_FW(x, dlt, n-1)-dfNumN_FW(x-dlt, dlt, n-1))/dlt;
}

float dfNumN_CT(float const x, float const dlt, int const n){
    if (n==0) return f(x);
    else return (dfNumN_FW(x+dlt, dlt, n-1)-dfNumN_FW(x-dlt, dlt, n-1))/2*dlt;
}*/










#endif