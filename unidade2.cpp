#include <iostream>
#include <math.h>
#include "utils/function.hpp"
#include "utils/function2D.hpp"
#include "INTEGRACAO_NUMERICA/exponencial/exponencial.hpp"
#include "INTEGRACAO_NUMERICA/gaussLegendre/integralDupla.hpp"

using namespace std;

int main() {

    // ============================================================
    // TESTE 1: Exponencial Simples e Dupla - função singular
    // f(x) = 1/sqrt(x), [0,1], integral exata = 2.0
    // ============================================================
    cout << "TESTE 1: f(x) = 1/sqrt(x), [0,1], exato = 2.0" << endl;
    cout << "----------------------------------------------" << endl;
    Function f1("1/sqrt(x)");

    cout << "\n-- Exponencial Simples (SE) --" << endl;
    testExponencial(intgExpSimples, f1, 0.0, 1.0, /*T=*/3.0, /*tol=*/1e-6);

    cout << "\n-- Exponencial Dupla (DE) --" << endl;
    testExponencial(intgExpDupla, f1, 0.0, 1.0, /*T=*/4.0, /*tol=*/1e-10);

    // ============================================================
    // TESTE 2: Exponencial Dupla - função suave, para comparar
    // f(x) = e^x, [0,1], integral exata = e - 1 ≈ 1.718281828
    // ============================================================
    cout << "\n\nTESTE 2: f(x) = e^x, [0,1], exato = " << (exp(1.0) - 1.0) << endl;
    cout << "----------------------------------------------" << endl;
    Function f2("exp(x)");
    testExponencial(intgExpDupla, f2, 0.0, 1.0, /*T=*/4.0, /*tol=*/1e-12);

    // ============================================================
    // TESTE 3: Integral dupla - volume sob z = x^2 + y^2
    // sobre [0,1] x [0,1], exato = 2/3 ≈ 0.666666...
    // ============================================================
    cout << "\n\nTESTE 3: Volume sob z = x^2 + y^2, [0,1]x[0,1]" << endl;
    cout << "----------------------------------------------" << endl;
    Function2D f3("x^2 + y^2");
    double volume = volumeSobSuperficie(f3, 0.0, 1.0, 0.0, 1.0);
    cout << "Volume calculado: " << volume << endl;
    cout << "Volume exato:     " << (2.0/3.0) << endl;
    cout << "Erro absoluto:    " << fabs(volume - 2.0/3.0) << endl;

    // ============================================================
    // TESTE 4: Área de um retângulo simples (f(x,y) = 1)
    // sobre [0,2] x [0,3] -> área exata = 6.0
    // ============================================================
    cout << "\n\nTESTE 4: Area do retangulo [0,2] x [0,3]" << endl;
    cout << "----------------------------------------------" << endl;
    double area = areaRetangulo(0.0, 2.0, 0.0, 3.0);
    cout << "Area calculada: " << area << endl;
    cout << "Area exata:     " << 6.0 << endl;

    return 0;
}