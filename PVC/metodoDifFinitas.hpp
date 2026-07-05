/*
 * fdm.hpp
 * -------
 * Biblioteca (somente header, C++11+) que implementa o Método das
 * Diferenças Finitas (MDF) para Problemas de Valor de Contorno (PVC)
 * unidimensionais governados por uma Equação Diferencial Ordinária
 * linear de 2a ordem:
 *
 *      a(x) * u''(x) + b(x) * u'(x) + c(x) * u(x) = f(x),   x em [x0, xN]
 *      u(x0) = u0   (condição de contorno de Dirichlet à esquerda)
 *      u(xN) = uN   (condição de contorno de Dirichlet à direita)
 *
 * O domínio é dividido em N subintervalos iguais de comprimento h,
 * gerando N+1 nós (i = 0, 1, ..., N). Os valores u(x0) e u(xN) são
 * conhecidos (condições de contorno); as incógnitas são u_1, ..., u_{N-1}.
 *
 * A equação diferencial é discretizada usando diferenças centrais:
 *
 *   u''(x_i) ~= [u_{i-1} - 2 u_i + u_{i+1}] / h^2
 *   u'(x_i)  ~= [u_{i+1} - u_{i-1}] / (2 h)
 *
 * o que gera, para cada nó interno i, uma equação da forma
 *
 *   A_i * u_{i-1} + B_i * u_i + C_i * u_{i+1} = D_i
 *
 * com
 *   A_i = a(x_i)/h^2 - b(x_i)/(2h)
 *   B_i = -2 a(x_i)/h^2 + c(x_i)
 *   C_i = a(x_i)/h^2 + b(x_i)/(2h)
 *   D_i = f(x_i)
 *
 * O sistema tridiagonal resultante (N-1 equações, N-1 incógnitas) é
 * resolvido pelo Algoritmo de Thomas (eliminação de Gauss especializada
 * para matrizes tridiagonais), com custo O(N).
 *
 * Uso típico:
 *
 *   fdm::BVP2 problema {
 *       [](double){ return 1.0; },   // a(x)
 *       [](double){ return 7.0; },   // b(x)
 *       [](double){ return -1.0; },  // c(x)
 *       [](double){ return 2.0; }    // f(x)
 *   };
 *
 *   fdm::FiniteDifferenceSolver solver(problema);
 *   fdm::Result r = solver.solve(0.0, 2.0, 10.0, 1.0, 0.1);
 *   r.print();
 *
 * Compilação: header-only, basta #include "fdm.hpp" (padrão C++11 ou
 * superior, por causa de std::function e das lambdas).
 */

#ifndef METODO_DIF_FINITAS_HPP
#define METODO_DIF_FINITAS_HPP

#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <cstdio>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

namespace fdm {

/* Tipo de coeficiente/termo-fonte que pode variar com x.
 * Para coeficientes constantes, basta usar uma lambda que ignora x. */
using CoefFn = std::function<double(double)>;

/* ------------------------------------------------------------------ */
/* Descrição do PVC: a(x) u'' + b(x) u' + c(x) u = f(x)                */
/* ------------------------------------------------------------------ */
struct BVP2 {
    CoefFn a; /* coeficiente de u''(x) */
    CoefFn b; /* coeficiente de u'(x)  */
    CoefFn c; /* coeficiente de u(x)   */
    CoefFn f; /* termo fonte / lado direito */
};

/* ------------------------------------------------------------------ */
/* Resultado: nós da grade e solução aproximada                       */
/* ------------------------------------------------------------------ */
struct Result {
    std::vector<double> x; /* posições dos nós, tamanho n+1 */
    std::vector<double> u; /* solução aproximada, tamanho n+1 */
    double h = 0.0;        /* espaçamento efetivo da grade */
    int n = 0;             /* número de subintervalos */

    int size() const { return static_cast<int>(x.size()); }

    void print(std::ostream &os = std::cout) const {
        os << "   x        u(x)\n";
        os << "--------------------\n";
        os << std::fixed << std::setprecision(6);
        for (int i = 0; i <= n; ++i) {
            os << "x = " << std::setprecision(2) << x[i]
               << "   u = " << std::setprecision(6) << u[i] << "\n";
        }
    }

    /* Grava a solução em um arquivo CSV (colunas x,u). */
    void to_csv(const std::string &filename) const {
        std::ofstream fp(filename);
        if (!fp) throw std::runtime_error("Nao foi possivel abrir " + filename);
        fp << "x,u\n";
        fp << std::fixed << std::setprecision(6);
        for (int i = 0; i <= n; ++i) {
            fp << std::setprecision(4) << x[i] << "," << std::setprecision(6) << u[i] << "\n";
        }
    }
};

/* ------------------------------------------------------------------ */
/* Algoritmo de Thomas para sistemas tridiagonais                     */
/* ------------------------------------------------------------------ */
/*
 * Resolve o sistema tridiagonal:
 *
 *   sub[i]*x[i-1] + diag[i]*x[i] + sup[i]*x[i+1] = rhs[i],  i = 0..m-1
 *
 * onde sub[0] e sup[m-1] não são usados (não existem vizinhos fora
 * do domínio). Lança std::runtime_error se algum pivô se anular.
 */
inline std::vector<double> thomas_solve(const std::vector<double> &sub,
                                         const std::vector<double> &diag,
                                         const std::vector<double> &sup,
                                         const std::vector<double> &rhs)
{
    const int m = static_cast<int>(diag.size());
    if (m <= 0) throw std::invalid_argument("Sistema vazio.");

    std::vector<double> c_star(m), d_star(m), sol(m);

    double piv = diag[0];
    if (std::fabs(piv) < 1e-14)
        throw std::runtime_error("Pivo nulo no algoritmo de Thomas (i=0).");

    c_star[0] = sup[0] / piv;
    d_star[0] = rhs[0] / piv;

    for (int i = 1; i < m; ++i) {
        piv = diag[i] - sub[i] * c_star[i - 1];
        if (std::fabs(piv) < 1e-14)
            throw std::runtime_error("Pivo nulo no algoritmo de Thomas (i=" + std::to_string(i) + ").");
        c_star[i] = sup[i] / piv;
        d_star[i] = (rhs[i] - sub[i] * d_star[i - 1]) / piv;
    }

    sol[m - 1] = d_star[m - 1];
    for (int i = m - 2; i >= 0; --i) {
        sol[i] = d_star[i] - c_star[i] * sol[i + 1];
    }

    return sol;
}

/* ------------------------------------------------------------------ */
/* Solver do PVC por Diferenças Finitas                                */
/* ------------------------------------------------------------------ */
class FiniteDifferenceSolver {
public:
    explicit FiniteDifferenceSolver(BVP2 problem) : problem_(std::move(problem)) {}

    /*
     * Resolve o PVC no intervalo [x0, xN] com condições de contorno
     * u(x0) = u0 e u(xN) = uN, usando passo aproximadamente igual a dx
     * (o número de subintervalos é arredondado para garantir uma grade
     * uniforme com N = round((xN - x0) / dx) subintervalos).
     */
    Result solve(double x0, double xN, double u0, double uN, double dx) const {
        if (dx <= 0.0 || xN <= x0)
            throw std::invalid_argument("Intervalo ou passo invalidos.");

        const int n = static_cast<int>(std::lround((xN - x0) / dx));
        if (n < 2)
            throw std::invalid_argument("Numero de subintervalos insuficiente (use um dx menor).");

        const double h = (xN - x0) / n;

        Result r;
        r.n = n;
        r.h = h;
        r.x.resize(n + 1);
        r.u.resize(n + 1);
        for (int i = 0; i <= n; ++i) r.x[i] = x0 + i * h;

        const int m = n - 1; /* número de incógnitas: nós internos 1..n-1 */
        std::vector<double> sub(m), diag(m), sup(m), rhs(m);

        for (int k = 0; k < m; ++k) {
            const int i = k + 1;
            const double xi = r.x[i];
            const double a = problem_.a(xi);
            const double b = problem_.b(xi);
            const double c = problem_.c(xi);
            const double f = problem_.f(xi);

            const double A_i = a / (h * h) - b / (2.0 * h); /* coef. de u_{i-1} */
            const double B_i = -2.0 * a / (h * h) + c;      /* coef. de u_i     */
            const double C_i = a / (h * h) + b / (2.0 * h); /* coef. de u_{i+1} */
            double D_i = f;

            diag[k] = B_i;

            if (i - 1 == 0) {
                D_i -= A_i * u0;   /* valor conhecido vai para o lado direito */
                sub[k] = 0.0;
            } else {
                sub[k] = A_i;
            }

            if (i + 1 == n) {
                D_i -= C_i * uN;   /* valor conhecido vai para o lado direito */
                sup[k] = 0.0;
            } else {
                sup[k] = C_i;
            }

            rhs[k] = D_i;
        }

        const std::vector<double> sol = thomas_solve(sub, diag, sup, rhs);

        r.u.front() = u0;
        r.u.back()  = uN;
        for (int k = 0; k < m; ++k) r.u[k + 1] = sol[k];

        return r;
    }

private:
    BVP2 problem_;
};

} /* namespace fdm */

#endif /* FDM_HPP */