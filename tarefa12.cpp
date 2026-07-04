#include <iostream>
#include "utils/algebraLinear.hpp"
#include "AUTOVALORES_AUTOVETORES/householderQR.hpp"
#include "AUTOVALORES_AUTOVETORES/metodoPotencia.hpp"

using namespace std;

int main() {
    Matriz A = {
        {40.0, 8.0, 4.0, 2.0, 1.0},
        { 8.0,30.0,12.0, 6.0, 2.0},
        { 4.0,12.0,20.0, 1.0, 2.0},
        { 2.0, 6.0, 1.0,25.0, 4.0},
        { 1.0, 2.0, 2.0, 4.0, 5.0}
    };
    int n = A.size();
    Vetor x0(n, 1.0);

    // ---- PARTE 1: Householder -> A# e H ----
    Matriz H;
    Matriz Atil = householderTridiagonalizacao(A, H);

    cout << "Matriz tridiagonal A#:" << endl;
    imprimirMatriz(Atil);

    cout << "\nMatriz acumulada H:" << endl;
    imprimirMatriz(H);

    // ---- PARTE 3: Metodos da potencia sobre A# ----
    cout << "\n===== POTENCIA SOBRE A# =====" << endl;

    auto pRegular = metodoPotencia(Atil, x0);
    cout << "Regular (maior): " << pRegular.autovalor
         << " (" << pRegular.iteracoes << " iter.)" << endl;

    auto pInversa = metodoPotenciaInversa(Atil, x0);
    cout << "Inversa (menor): " << pInversa.autovalor
         << " (" << pInversa.iteracoes << " iter.)" << endl;

    vector<double> shifts = {12.0, 24.0, 31.0};
    vector<ResultadoAutovalor> intermediarios;
    for (double s : shifts) {
        auto r = metodoPotenciaDeslocamento(Atil, x0, s);
        intermediarios.push_back(r);
        cout << "Deslocamento (shift=" << s << "): " << r.autovalor
             << " (" << r.iteracoes << " iter.)" << endl;
    }

    // ---- PARTE 4: autovetores de A# -> autovetores de A (v_A = H * v_A#) ----
    cout << "\n===== AUTOVETORES DE A (via H * v) =====" << endl;

    Vetor vA_regular = matVec(H, pRegular.autovetor);
    cout << "Do maior autovalor: ";
    for (double v : vA_regular) cout << v << " ";
    cout << endl;

    Vetor vA_inversa = matVec(H, pInversa.autovetor);
    cout << "Do menor autovalor: ";
    for (double v : vA_inversa) cout << v << " ";
    cout << endl;

    for (size_t i = 0; i < intermediarios.size(); i++) {
        Vetor vA = matVec(H, intermediarios[i].autovetor);
        cout << "Shift=" << shifts[i] << ": ";
        for (double v : vA) cout << v << " ";
        cout << endl;
    }

    // ---- PARTE 5: autovalores de A (iguais aos de A#) ----
    cout << "\n===== AUTOVALORES DE A =====" << endl;
    cout << pRegular.autovalor << ", " << pInversa.autovalor << ", ";
    for (auto &r : intermediarios) cout << r.autovalor << ", ";
    cout << endl;

    // ---- Verificacao cruzada: pipeline Householder+QR (suas duas caixas pretas) ----
    cout << "\n===== VERIFICACAO: pipeline autovaloresSimetrica() =====" << endl;
    auto resultadoQR = autovaloresSimetrica(A);
    cout << "Autovalores (Householder+QR): ";
    for (double lambda : resultadoQR.autovalores) cout << lambda << " ";
    cout << "\nIteracoes QR: " << resultadoQR.iteracoes << endl;

    return 0;
}