#include "src/head.h"
#include "src/find_nums.cpp"
#include "src/make_vecs.cpp"
#include "src/timing.cpp"
#include "src/make_table_file.cpp"
using namespace arma;

vec rref_tridiagonal(int n, double* a, double* b, double* c, double h, vec f) {
    double *b_new = new double[n+2];
    double *f_new = new double[n+2];
    vec v = zeros<vec>(n+2);
    double a_div_bnew;

    b_new[1] = b[1];
    f_new[1] = f[1];
    for (int i=1; i < n; i++) {
        a_div_bnew = a[i] / b_new[i];
        b_new[i+1] = b[i+1] - a_div_bnew*c[i];
        f_new[i+1] = f[i+1] - a_div_bnew*f_new[i];
    }
	for (int i=n; i>0; i--) {
		v[i] = (f_new[i] - c[i]*v[i+1])/b_new[i];
	}

    return v;
}

vec rref_tridiagonal_3num(int n, double a, double b, double c, double h, vec f) {
    double *b_new = new double[n+2];
    double *f_new = new double[n+2];
    vec v = zeros<vec>(n+2);
    double a_div_bnew;

    b_new[1] = b;
    f_new[1] = f[1];
    for (int i=1; i < n; i++) {
        a_div_bnew = a / b_new[i];
        b_new[i+1] = b - a_div_bnew*c;
        f_new[i+1] = f[i+1] - a_div_bnew*f_new[i];
    }
	for (int i=n; i>0; i--) {
		v[i] = (f_new[i] - c*v[i+1])/b_new[i];
	}
    return v;
}

// int main(int argc, char *argv[]) {
int do_for_one_n(int n) {
    // Define numbers
    double h = 1. / (n + 1);

    // Numerical calculation
    double *a = make_abc(n, -1);
    double *b = make_abc(n, 2);
    double *c = make_abc(n, -1);
    vec y = make_y(n, f);
    auto start = time_start();
    vec v = rref_tridiagonal(n, a, b, c, h, y);
    cout << "Tridiagonal time: "; time_finish(start);
    start = time_start();
    vec v_3num = rref_tridiagonal_3num(n, -1., 2., -1., h, y);
    cout << "Tridiagonal 3num time: "; time_finish(start);


    // Analytical calculation
    vec u = make_u_anal(n);

    // Make file
    if (n < 10000) {
        vec *out_array = new vec[4];
        out_array[0] = make_h(n);
        out_array[1] = u;
        out_array[2] = v;
        out_array[3] = v_3num;
        make_table_file("tridiagonal_n=" + to_string(n), n+2, 4, out_array);
    }


    // Print
    cout << "Error when n=" << n << ": " << error_max(n, v, u) << endl;
    cout << "Error 3num when n=" << n << ": " << error_max(n, v_3num, u) << endl;
    cout << "log(h): " << log10(h) << endl;
    cout << '\n';

    return 0;
}

int main() {
    for (int i = 1; i < 8; i++) {
        do_for_one_n(pow(10, i) );
    }

    return 0;
}

// Tables:
// x | Analytix | Numerical
