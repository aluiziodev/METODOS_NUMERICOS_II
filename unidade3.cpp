#include <iostream>
#include <iomanip>
#include "utils/algebraLinear.hpp"
#include "AUTOVALORES_AUTOVETORES/metodoPotencia.hpp"
#include "AUTOVALORES_AUTOVETORES/householderQR.hpp"

using namespace std;

void imprimirVetor(const Vetor &v) {
    cout << "[ ";
    for (double x : v) cout << fixed << setprecision(6) << x << " ";
    cout << "]";
}

void resolverMatriz(const string &nome, const Matriz &A, const Vetor &x0,
                     const vector<double> &shifts) {
    cout << "===================================================" << endl;
    cout << "MATRIZ " << nome << endl;
    cout << "===================================================" << endl;

    // Algoritmo 2.1a: Potência Regular -> maior autovalor
    auto r1 = metodoPotencia(A, x0);
    cout << "\n[Potencia Regular] Maior autovalor:" << endl;
    cout << "Autovalor: " << r1.autovalor << ", Iteracoes: " << r1.iteracoes << endl;
    cout << "Autovetor: "; imprimirVetor(r1.autovetor); cout << endl;

    // Algoritmo 2.1b: Potência Inversa -> menor autovalor (em modulo)
    auto r2 = metodoPotenciaInversa(A, x0);
    cout << "\n[Potencia Inversa] Menor autovalor:" << endl;
    cout << "Autovalor: " << r2.autovalor << ", Iteracoes: " << r2.iteracoes << endl;
    cout << "Autovetor: "; imprimirVetor(r2.autovetor); cout << endl;

    // Algoritmo 2.2: Potência com Deslocamento -> autovalores restantes
    int idx = 1;
    for (double s : shifts) {
        auto r3 = metodoPotenciaDeslocamento(A, x0, s);
        cout << "\n[Potencia com Deslocamento, shift=" << s << "] Autovalor "
             << idx++ << ":" << endl;
        cout << "Autovalor: " << r3.autovalor << ", Iteracoes: " << r3.iteracoes << endl;
        cout << "Autovetor: "; imprimirVetor(r3.autovetor); cout << endl;
    }
    cout << endl;
}

int main() {
    // ============================================================
    // MATRIZ A1 (3x3) - 1 autovalor intermediario (shift=2.0)
    // ============================================================
    Matriz A1 = {
        {5.0, 2.0, 1.0},
        {2.0, 3.0, 1.0},
        {1.0, 1.0, 2.0}
    };
    Vetor x0_A1 = {1.0, 1.0, 1.0};
    resolverMatriz("A1", A1, x0_A1, {2.0});

    // ============================================================
    // MATRIZ A2 (3x3) - 1 autovalor intermediario (shift=-10.0)
    // ============================================================
    Matriz A2 = {
        {-14.0, 1.0, -2.0},
        {  1.0,-1.0,  1.0},
        { -2.0, 1.0,-11.0}
    };
    Vetor x0_A2 = {1.0, 1.0, 1.0};
    resolverMatriz("A2", A2, x0_A2, {-10.0});

    // ============================================================
    // MATRIZ A3 (5x5) - 3 autovalores intermediarios (shifts=12,24,31)
    // ============================================================
    Matriz A3 = {
        {40.0, 8.0, 4.0, 2.0, 1.0},
        { 8.0,30.0,12.0, 6.0, 2.0},
        { 4.0,12.0,20.0, 1.0, 2.0},
        { 2.0, 6.0, 1.0,25.0, 4.0},
        { 1.0, 2.0, 2.0, 4.0, 5.0}
    };
    Vetor x0_A3 = {1.0, 1.0, 1.0, 1.0, 1.0};
    resolverMatriz("A3", A3, x0_A3, {12.0, 24.0, 31.0});

    return 0;
}