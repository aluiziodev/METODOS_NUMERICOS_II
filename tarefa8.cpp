#include <iostream>
#include <iomanip>
#include "utils/function2D.hpp"
#include "utils/areaSuperficieCircular.hpp"
 
// Exemplo de uso: Aula 15, Problema 2
// Paraboloide hiperbolico z = f(x,y) = 0.2(x^2 - y^2)
// Regiao U = { (x,y) | x^2/1600 + y^2/1600 <= 1 }  ->  circulo de raio R = 40
 
int main(){
    Function2D f("0.2*(x^2 - y^2)");
 
    double R = 40.0; // raio: x^2/1600+y^2/1600<=1  =>  R = sqrt(1600) = 40
    double A = areaSuperficieRegiaoCircular(f, R);
 
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Area da superficie (Problema 2) = " << A << " m^2" << std::endl;
 
    // Caso mais geral: regiao eliptica x^2/a^2 + y^2/b^2 <= 1
    // double A_elipse = areaSuperficieRegiaoEliptica(f, a, b);
 
    return 0;
}