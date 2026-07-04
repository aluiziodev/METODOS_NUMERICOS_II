#include <iostream>
#include "utils/algebraLinear.hpp"
#include "AUTOVALORES_AUTOVETORES/metodoPotencia.hpp"
#include "AUTOVALORES_AUTOVETORES/householderQR.hpp"

using namespace std;

void imprimirVetor(const Vetor &v) {
    cout << "[ ";
    for (double x : v) cout << x << " ";
    cout << "]" << endl;
}

int main() {
    // Matriz simétrica de teste: autovalores exatos = 1.0, 2.26794919, 5.73205081
    Matriz A = {
        {4.0, 1.0, 2.0},
        {1.0, 2.0, 0.0},
        {2.0, 0.0, 3.0}
    };
    Vetor x0 = {1.0, 1.0, 1.0};

    cout << "TESTE 1: Metodo da Potencia Regular (maior autovalor)" << endl;
    auto r1 = metodoPotencia(A, x0);
    cout << "Autovalor: " << r1.autovalor << ", Iteracoes: " << r1.iteracoes << endl;
    cout << "Autovetor: "; imprimirVetor(r1.autovetor);

    cout << "\nTESTE 2: Metodo da Potencia Inversa (menor autovalor)" << endl;
    auto r2 = metodoPotenciaInversa(A, x0);
    cout << "Autovalor: " << r2.autovalor << ", Iteracoes: " << r2.iteracoes << endl;
    cout << "Autovetor: "; imprimirVetor(r2.autovetor);

    cout << "\nTESTE 3: Metodo da Potencia com Deslocamento (shift=2.0)" << endl;
    auto r3 = metodoPotenciaDeslocamento(A, x0, 2.0);
    cout << "Autovalor: " << r3.autovalor << ", Iteracoes: " << r3.iteracoes << endl;
    cout << "Autovetor: "; imprimirVetor(r3.autovetor);

    cout << "\nTESTE 4: Householder + QR (todos os autovalores)" << endl;
    auto r4 = autovaloresSimetrica(A);
    cout << "Autovalores: "; imprimirVetor(r4.autovalores);
    cout << "Iteracoes QR: " << r4.iteracoes << endl;

    return 0;
}