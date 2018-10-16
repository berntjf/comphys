// Project 2b
#include "head.h"
#include "toeplitz.cpp"
using namespace arma;

/*
Tester det analytiske uttrykket mot armadillos egenverdi-funksjon.
*/
int test_lambda_anal() {
    int N = 50;
    double h = 1.0/N;
    double d = 2.0/pow(h, 2);
    double a = -1.0/pow(h, 2);
    mat A(N-1, N-1);
    for (int i = 0; i < N-1; i++) {
        for (int j = 0; j < N-1; j++) {
            if (i == j) {
                A(i,j) = d;
            } else if ((i+1 == j) || (i == j+1)) {
                A(i,j) = a;
            } else {
                A(i,j) = 0.0;
            }
        }
    }
    vec v = eig_sym(A);

    double computed;
    double expected;
    for (int j = 1; j < N; j++) {
        computed = v(j-1);
        expected = lambda_anal(j, N, a, d);
        if ((computed/expected)>1.1 || (expected/computed)>1.1) {
            cout << "Computed " << computed << ", expected " << expected << " for j=" << j << endl;
            return 1;
        }
    }

    return 0;
}

/*
Tester om find_kl_when_a2_max-funksjonen i toepliz gir riktig resultat.
*/
int test_find_kl_when_a2_max() {
    int N = 5;
    double** matrix = new double*[N];
    for (int i = 0; i < N; i++) {
        matrix[i] = new double[N];
        for (int j = 0; j < N; j++) {
            if      (i == 0)           matrix[i][j] =  2.0;
            else if ((i==2) && (j==3)) matrix[i][j] = -3.0;
            else                       matrix[i][j] =  0.0;
        }
    }
    int* computed = find_kl_when_a2_max(matrix, N, 0.1);
    int com_k = computed[0];
    int com_l = computed[1];
    int exp_k = 2;
    int exp_l = 3;
    if ((com_k != exp_k) || (com_l != exp_l)) {
        cout << "Computed " << com_k << " and " << com_l << ", expected " << exp_k << " and " << exp_l << endl;
        return 1;
    }
    return 0;
}

int test_jacobi() {
    int result = 0;
    double a, b, computed, expected;
    for (int N = 5; N < 20; N = N + 5) {
        a = -pow(N, 2.0);
        b = 2*pow(N, 2.0);
        double** A = make_buckling_beam_matrix(N);
        jacobi_buckling_beam(N, true, 200);
        for (int j = 1; j < N; j++) {
            computed = A[j-1][j-1];
            expected = lambda_anal(j, N, a, b);
            if (((computed/expected) > 1.01) || ((computed/expected) < 0.99)) {
                cout << "Computed " << A[j-1][j-1] << ", expected " << lambda_anal(j, N, a, b) << " as lambda_j for j=" << j << endl;
                result = 1;
            }
        }
        delete A;
    }
    return result;
}

int main() {
    test_lambda_anal();
    test_find_kl_when_a2_max();
    test_jacobi();
}
