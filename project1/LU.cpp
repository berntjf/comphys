#include "src/head.h"
#include "src/find_nums.cpp"
#include "src/make_vecs.cpp"
#include "src/timing.cpp"
using namespace arma;

int do_for_one_n(int n) {
    // Define numbers and arrays
    double h = 1. / (n + 1);
    vec y = make_y(n, f);
    double *a_i = make_abc(n, -1.);
    double *b_i = make_abc(n, 2.);
    double *c_i = make_abc(n, -1.);
    mat A = make_A(n);
    vec v;

    // Numerical calculation
    auto start = time_start();
    solve(v, A, y);
    cout << "LU time: "; time_finish(start);
    vec u = zeros<vec>(n+2);
    for (int i=0; i<n+1; i++) {
        u[i] = u_anal(i*h);
    }

    // Print
    cout << "Error when n=" << n << ": " << error_max(n, v, u) << endl;
    cout << "log(h): " << log10(h) << endl;
    cout << '\n';

    return 0;
}

int main() {
    for (int i = 1; i < 4; i++) {
        do_for_one_n(pow(10, i) );
    }

    return 0;
}
