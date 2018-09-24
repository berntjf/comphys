#include <cmath>
#include "head.h"

double f(double x) {
    return 100.0*exp(-10.0*x);
}

double u_anal(double x) {
    return 1. - (1. - exp(-10.))*x - exp(-10.*x);
}

double error_max(int n, arma::vec v, arma::vec u) {
    double error_max = log10(abs((v[1]-u[1])/u[1]));
    double current;
    for (int i = 1; i < n+1; i++) {
        current = log10(fabs((v[i]-u[i])/u[i]));
        if (current > error_max) {
            error_max = current;
        }
    }
    return error_max;
}
