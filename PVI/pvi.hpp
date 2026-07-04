#ifndef PVI_HPP
#define PVI_HPP
#include <vector>
#include <functional>
#include <cmath>
using Estado = std::vector<double>;
using FuncaoPVI = std::function<Estado(const Estado&, double)>;

Estado somar(const Estado &a, const Estado &b) {
    Estado r(a.size());
    for (size_t i = 0; i < a.size(); i++) r[i] = a[i] + b[i];
    return r;
}
Estado escalar(double c, const Estado &a) {
    Estado r(a.size());
    for (size_t i = 0; i < a.size(); i++) r[i] = c * a[i];
    return r;
}
#endif