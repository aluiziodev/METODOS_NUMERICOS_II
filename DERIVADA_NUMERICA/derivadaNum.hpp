#ifndef DERIVADA_NUM_HPP
#define DERIVADA_NUM_HPP
#include <bits/stdc++.h>

float f(float const x){
    return sin(x)*pow(x, 3)*5;
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









#endif