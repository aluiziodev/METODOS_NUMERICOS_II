/*
 * fem.hpp
 * -------
 * Biblioteca (somente header, C++11+) que implementa o Método dos
 * Elementos Finitos (MEF), na formulação de Galerkin com resíduos
 * ponderados, para Problemas de Valor de Contorno (PVC) unidimensionais
 * governados por uma Equação Diferencial Ordinária linear de 2a ordem:
 *
 *      a(x) * u''(x) + b(x) * u'(x) + c(x) * u(x) = f(x),   x em [x0, xN]
 *      u(x0) = u0   (condição de contorno de Dirichlet à esquerda)
 *      u(xN) = uN   (condição de contorno de Dirichlet à direita)
 *
 * ------------------------------------------------------------------
 * Formulação fraca (resíduos ponderados)
 * ------------------------------------------------------------------
 * Multiplicando a ED por uma função peso w(x) que se anula no contorno
 * (pois as condições são de Dirichlet) e integrando em [x0, xN]:
 *
 *   integral( a u'' w + b u' w + c u w ) dx = integral( f w ) dx
 *
 * Como as funções de interpolação de Lagrange lineares (funções-chapéu)
 * usadas em u_h(x) = sum_i phi_i(x) u_i não são duas vezes diferenciáveis
 * (são polinômios de 1o grau em cada elemento), integra-se por partes o
 * termo com a derivada segunda:
 *
 *   integral( a u' w ) dx = [a u' w] - integral( a u' w' ) dx
 *
 * O termo de contorno [a u' w] se anula porque w = 0 em x0 e em xN.
 * Assim, a formulação fraca fica
 *
 *   - integral( a u' w' ) dx + integral( b u' w ) dx
 *                            + integral( c u w  ) dx = integral( f w ) dx,
 *   para toda função peso w admissível.
 *
 * Substituindo u_h(x) = sum_k phi_k(x) u_k e tomando w = phi_j (método de
 * Galerkin: função peso no mesmo espaço da solução), obtém-se, para cada
 * nó j onde não há condição de Dirichlet, uma equação
 *
 *   sum_k [ -a (phi_k', phi_j') + b (phi_k', phi_j) + c (phi_k, phi_j) ] u_k
 *       = (f, phi_j)
 *
 * onde (.,.) denota a integral do produto no domínio. Isso gera o sistema
 * linear K u = F.
 *
 * ------------------------------------------------------------------
 * Discretização e montagem por elementos
 * ------------------------------------------------------------------
 * O domínio é dividido em N elementos (subintervalos) de comprimento h,
 * com N+1 nós. Em cada elemento el_e = [x_i, x_{i+1}], usam-se as duas
 * funções de interpolação lineares locais (N1 associada ao nó esquerdo,
 * N2 ao nó direito). Parametrizando o elemento em coordenadas naturais
 * (parametrização de Legendre) xi em [-1, 1]:
 *
 *   N1(xi) = (1 - xi)/2         N2(xi) = (1 + xi)/2
 *   dN1/dx = -1/h               dN2/dx =  1/h
 *   x(xi) = x_i + h*(xi+1)/2 ,  dx/dxi = h/2
 *
 * Como phi_i e phi_j só têm interseção não-nula dentro de um mesmo
 * elemento (ou quando i == j), a matriz global K é tridiagonal — cada
 * elemento contribui apenas para as 2x2 entradas associadas aos seus
 * dois nós. As integrais dos elementos são calculadas por quadratura
 * de Gauss-Legendre de 2 pontos (exata para os polinômios de grau <= 3
 * envolvidos quando a, b, c, f são constantes; aproximada, mas precisa,
 * quando esses coeficientes variam com x).
 *
 * Após a montagem, as condições de contorno de Dirichlet são impostas
 * eliminando as equações dos nós 0 e N e passando suas contribuições
 * (multiplicadas por u0 e uN) para o lado direito do sistema reduzido,
 * exatamente como descrito na Aula #28 do curso. O sistema tridiagonal
 * resultante, de tamanho N-1, é resolvido pelo Algoritmo de Thomas.
 *
 * Uso típico:
 *
 *   fem::BVP2 problema {
 *       [](double){ return 1.0; },   // a(x)
 *       [](double){ return 7.0; },   // b(x)
 *       [](double){ return -1.0; },  // c(x)
 *       [](double){ return 2.0; }    // f(x)
 *   };
 *
 *   fem::GalerkinFEMSolver solver(problema);
 *   fem::Result r = solver.solve(0.0, 2.0, 10.0, 1.0, 0.1);
 *   r.print();
 *
 * Compilação: header-only, basta #include "fem.hpp" (C++11 ou superior).
 */

#ifndef METODO_ELEM_FINITOS_HPP
#define METODO_ELEM_FINITOS_HPP

#include <vector>
#include <array>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace fem {

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
/* Resultado: nós da grade e solução aproximada nos nós                */
/* ------------------------------------------------------------------ */
struct Result {
    std::vector<double> x; /* posições dos nós, tamanho n+1 */
    std::vector<double> u; /* solução aproximada nos nós, tamanho n+1 */
    double h = 0.0;        /* comprimento de cada elemento (malha uniforme) */
    int n = 0;             /* número de elementos */

    int size() const { return static_cast<int>(x.size()); }

    void print(std::ostream &os = std::cout) const {
        os << "   x        u(x)\n";
        os << "--------------------\n";
        for (int i = 0; i <= n; ++i) {
            os << "x = " << std::fixed << std::setprecision(2) << x[i]
               << "   u = " << std::setprecision(6) << u[i] << "\n";
        }
    }

    void to_csv(const std::string &filename) const {
        std::ofstream fp(filename);
        if (!fp) throw std::runtime_error("Nao foi possivel abrir " + filename);
        fp << "x,u\n";
        for (int i = 0; i <= n; ++i) {
            fp << std::fixed << std::setprecision(4) << x[i] << ","
               << std::setprecision(6) << u[i] << "\n";
        }
    }
};

/* ------------------------------------------------------------------ */
/* Algoritmo de Thomas para sistemas tridiagonais                     */
/* ------------------------------------------------------------------ */
inline std::vector<double> thomas_solve(std::vector<double> sub,
                                         std::vector<double> diag,
                                         std::vector<double> sup,
                                         std::vector<double> rhs)
{
    const int m = static_cast<int>(diag.size());
    if (m <= 0) throw std::invalid_argument("Sistema vazio.");

    std::vector<double> sol(m);

    double piv = diag[0];
    if (std::fabs(piv) < 1e-14)
        throw std::runtime_error("Pivo nulo no algoritmo de Thomas (i=0).");

    sup[0] /= piv;
    rhs[0] /= piv;

    for (int i = 1; i < m; ++i) {
        piv = diag[i] - sub[i] * sup[i - 1];
        if (std::fabs(piv) < 1e-14)
            throw std::runtime_error("Pivo nulo no algoritmo de Thomas (i=" + std::to_string(i) + ").");
        sup[i] /= piv;
        rhs[i] = (rhs[i] - sub[i] * rhs[i - 1]) / piv;
    }

    sol[m - 1] = rhs[m - 1];
    for (int i = m - 2; i >= 0; --i) {
        sol[i] = rhs[i] - sup[i] * sol[i + 1];
    }

    return sol;
}

/* ------------------------------------------------------------------ */
/* Solver do PVC pelo Método dos Elementos Finitos (Galerkin, P1)       */
/* ------------------------------------------------------------------ */
class GalerkinFEMSolver {
public:
    explicit GalerkinFEMSolver(BVP2 problem) : problem_(std::move(problem)) {}

    /*
     * Resolve o PVC no intervalo [x0, xN] com condições de contorno
     * u(x0) = u0 e u(xN) = uN, usando N elementos lineares de
     * comprimento aproximadamente igual a dx (N = round((xN-x0)/dx),
     * garantindo malha uniforme).
     */
    Result solve(double x0, double xN, double u0, double uN, double dx) const {
        if (dx <= 0.0 || xN <= x0)
            throw std::invalid_argument("Intervalo ou passo invalidos.");

        const int n = static_cast<int>(std::lround((xN - x0) / dx)); /* nº de elementos */
        if (n < 2)
            throw std::invalid_argument("Numero de elementos insuficiente (use um dx menor).");

        const double h = (xN - x0) / n;

        Result r;
        r.n = n;
        r.h = h;
        r.x.resize(n + 1);
        r.u.resize(n + 1);
        for (int i = 0; i <= n; ++i) r.x[i] = x0 + i * h;

        /* Armazenamento tridiagonal da matriz global (nós 0..n).
         * diag[i]  = coeficiente de u_i na equação do nó i
         * sub[i]   = coeficiente de u_{i-1} na equação do nó i
         * sup[i]   = coeficiente de u_{i+1} na equação do nó i
         * rhs[i]   = lado direito da equação do nó i               */
        std::vector<double> diag(n + 1, 0.0), sub(n + 1, 0.0), sup(n + 1, 0.0), rhs(n + 1, 0.0);

        /* Quadratura de Gauss-Legendre de 2 pontos em [-1, 1] */
        const double gp = 1.0 / std::sqrt(3.0);
        const std::array<double, 2> xi_pts = { -gp, gp };
        const std::array<double, 2> w_pts  = { 1.0, 1.0 };

        /* Montagem elemento a elemento (assembly), como descrito na
         * "montagem alternativa" da Aula #28: cada elemento contribui
         * com uma matriz 2x2 e um vetor 2x1 para os nós que o compõem. */
        for (int e = 0; e < n; ++e) {
            const int i_left = e, i_right = e + 1;
            const double x_left = r.x[i_left];

            double Ke[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
            double Fe[2] = { 0.0, 0.0 };

            for (int q = 0; q < 2; ++q) {
                const double xi = xi_pts[q];
                const double w  = w_pts[q];

                const double N1 = (1.0 - xi) / 2.0;
                const double N2 = (1.0 + xi) / 2.0;
                const double dN1dx = -1.0 / h;
                const double dN2dx =  1.0 / h;
                const double N[2]    = { N1, N2 };
                const double dNdx[2] = { dN1dx, dN2dx };

                const double x_phys = x_left + h * (xi + 1.0) / 2.0;
                const double jac = h / 2.0; /* dx/dxi */

                const double a_val = problem_.a(x_phys);
                const double b_val = problem_.b(x_phys);
                const double c_val = problem_.c(x_phys);
                const double f_val = problem_.f(x_phys);

                const double weight = w * jac;

                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        Ke[j][k] += weight * ( -a_val * dNdx[k] * dNdx[j]
                                               + b_val * dNdx[k] * N[j]
                                               + c_val * N[k]    * N[j] );
                    }
                    Fe[j] += weight * f_val * N[j];
                }
            }

            /* Espalha (assembly) a contribuição do elemento na matriz global tridiagonal */
            diag[i_left]  += Ke[0][0];
            diag[i_right] += Ke[1][1];
            sup[i_left]   += Ke[0][1]; /* liga nó i_left  ao vizinho i_right */
            sub[i_right]  += Ke[1][0]; /* liga nó i_right ao vizinho i_left */
            rhs[i_left]   += Fe[0];
            rhs[i_right]  += Fe[1];
        }

        /* Imposição das condições de contorno de Dirichlet: as equações
         * dos nós 0 e n são descartadas; suas contribuições conhecidas
         * (u0 e uN) são passadas para o lado direito das equações dos
         * nós internos vizinhos (nós 1 e n-1). */
        const int m = n - 1; /* número de incógnitas (nós internos 1..n-1) */
        std::vector<double> A_sub(m), A_diag(m), A_sup(m), A_rhs(m);

        for (int k = 0; k < m; ++k) {
            const int i = k + 1;

            A_diag[k] = diag[i];
            A_rhs[k]  = rhs[i];

            if (i - 1 == 0) {
                A_rhs[k] -= sub[i] * u0;
                A_sub[k] = 0.0;
            } else {
                A_sub[k] = sub[i];
            }

            if (i + 1 == n) {
                A_rhs[k] -= sup[i] * uN;
                A_sup[k] = 0.0;
            } else {
                A_sup[k] = sup[i];
            }
        }

        const std::vector<double> sol = thomas_solve(A_sub, A_diag, A_sup, A_rhs);

        r.u.front() = u0;
        r.u.back()  = uN;
        for (int k = 0; k < m; ++k) r.u[k + 1] = sol[k];

        return r;
    }

private:
    BVP2 problem_;
};

} /* namespace fem */

#endif /* FEM_HPP */