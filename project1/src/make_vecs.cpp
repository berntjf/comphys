#include <armadillo>
#include "head.h"
using namespace arma;

double* make_abc(int n, float abc) {
    double *abc_i = new double[n+2];
    for (int i = 0; i < n + 1; i++)
        abc_i[i] = abc;
    return abc_i;
}

vec make_y(int n, double (*f)(double)) {
    vec y_i = zeros<vec>(n+2);
    double h = 1. / (n + 1);
    double h_sq = pow(h, 2.);
    for (int i = 0; i < n+1; i++)
        y_i[i] = h_sq * f(i*h);
    return y_i;
}

vec make_h(int n) {
    double h = 1. / (n + 1);
    vec h_i = zeros<vec>(n+2);
    for (int i = 0; i < n+2; i++) {
        h_i[i] = h*i;
    }
    return h_i;
}

vec make_u_anal(int n) {
    vec u = zeros<vec>(n+2);
    double h = 1. / (n + 1);
    for (int i=0; i<n+1; i++) {
        u[i] = u_anal(i*h);
    }
    return u;
}

mat make_A(int n) {
    mat A(n+2, n+2);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if ( (i == j+1) || (i+1 == j) ) {
                A(i,j) = -1.;
            } else if ( i == j ) {
                A(i,j) = 2.;
            } else {
                A(i,j) = 0.;
            }
    return A;
}
