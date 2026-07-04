#ifndef FUNCTION2D_HPP
#define FUNCTION2D_HPP
#include "../libs/exprtk.hpp"
using namespace exprtk;

class Function2D {
    std::string exp;
public:
    Function2D(std::string exp) : exp(exp) {}
    double value(double x, double y) {
        symbol_table<double> table;
        table.add_variable("x", x);
        table.add_variable("y", y);
        table.add_constants();
        expression<double> expr;
        expr.register_symbol_table(table);
        parser<double> parser;
        parser.compile(exp, expr);
        return expr.value();
    }
};
#endif