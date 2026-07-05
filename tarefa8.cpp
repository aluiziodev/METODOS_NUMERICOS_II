#include <iostream>
#include <iomanip>
#include "utils/function2D.hpp"
#include "utils/areaSuperficieCircular.hpp"

using namespace std;
 
int main(){
    Function2D f("0.2*(x^2 - y^2)");
 
    double R = 40.0; // raio: x^2/1600+y^2/1600<=1  =>  R = sqrt(1600) = 40
    double A = areaSuperficieRegiaoCircular(f, R);
 
    cout << fixed << setprecision(2);
    cout << "Area da superficie (Problema 2) = " << A << " m^2" << endl;
 
 
    return 0;
}