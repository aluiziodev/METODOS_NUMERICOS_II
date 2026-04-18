#ifndef FUNCTION_HPP
#define FUNCTION_HPP
#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
#include "../libs/exprtk.hpp"

using namespace std;
using namespace exprtk;


class Function{
    private:
    string exp;

    public:
    Function(string exp): exp(exp) {}

    double value(double x){
        symbol_table<double> table;
        table.add_variable("x", x);
        table.add_constants();

        expression<double> expr;
        expr.register_symbol_table(table);

        parser<double> parser;
        parser.compile(exp, expr);

        return expr.value();
    }

};


#endif 