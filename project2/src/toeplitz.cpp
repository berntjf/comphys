// Project 2b
#include "head.h"
#include <iomanip>

double lambda_anal(int j, int N, double a, double d) {
    return d + 2*a*cos(j*M_PI/(N+1));
}

/*
Finner det høyeste kvadratet av et ikke-diagonalt element i en kvadratisk matrise.
*/
int* find_kl_when_a2_max(double** A, int N, double epsilon) {
    double a2_max = 0.0;
    double current;
    int k = 0;
    int l = 0;
    for (int i = 0; i < N-1; ++i) {
        for (int j = 0; j < N-1; ++j) {
            if (i == j) continue;
            current = pow(A[i][j], 2.0);
            if (current > a2_max) {
                a2_max = current;
                k = i;
                l = j;
            }
        }
    }
    int* kl = new int[2];
    if (a2_max < epsilon) {
        kl[0] = -1;
        kl[1] = -1;
    } else {
        kl[0] = k;
        kl[1] = l;
    }
    return kl;
}

/*
Lager matrisen A.
*/
double** make_A_matrix(double* a, double* d, int N) {
    double** A = new double*[N+1];
    for (int i = 0; i < N+1; i++) {
        A[i] = new double[N+1];
    }
    for (int i = 0; i < N; i++) {
        A[i][i]   = d[i];
        A[i][i+1] = a[i];
        A[i+1][i] = a[i+1];
    }
    A[N][N] = d[N];
    return A;
}

double** make_buckling_beam_matrix(int N) {
    double* a = new double[N+1];
    double* d = new double[N+1];
    double N_sq = pow(N, 2.0);
    double a_value = -N_sq;
    double d_value = 2*N_sq;
    for (int i = 0; i < N+1; i++) {
        a[i] = a_value;
        d[i] = d_value;
    }
    double** A = make_A_matrix(a, d, N);
    delete[] a;
    delete[] d;
    return A;
}

double** make_quantum_dots_one_electron_matrix(int N) {
    double* a = new double[N+1];
    double* d = new double[N+1];
    double N_sq = pow(N, 2.0);
    double a_value = -N_sq;
    double d_value = 2*N_sq;
    for (int i = 0; i < N; i++) {
        a[i] = a_value;
        d[i] = d_value + pow(i+1, 2.0)/N_sq;
    }
    double** A = make_A_matrix(a, d, N);
    delete[] a;
    delete[] d;
    return A;
}

double** make_quantum_dots_two_electrons_matrix(int N, double omega_r) {
    double* a = new double[N+1];
    double* d = new double[N+1];
    double N_sq = pow(N, 2.0);
    double omega_r_sq = pow(omega_r, 2.0);
    double a_value = -N_sq;
    double d_value = 2*N_sq;
    for (int i = 0; i < N; i++) {
        a[i] = a_value;
        d[i] = d_value + omega_r_sq*pow(i+1, 2.0)/N_sq + N/(i+1.0);
    }
    double** A = make_A_matrix(a, d, N);
    delete[] a;
    delete[] d;
    return A;
}

void jacobi(double** A, int N, bool minimal_print, int max_iter) {
    double epsilon = 1e-8;
    int k, l;

    double** B = new double*[N+1];
    for (int i = 0; i < N+1; ++i) {
        B[i] = new double[N+1];
	}

    double tau, t, c, s;

	int count = 0;
    while (count < max_iter) {
        int* kl = find_kl_when_a2_max(A, N, epsilon);
        k = kl[0];
        l = kl[1];
        if (k == -1) {
            delete B;
            break;
        }
        tau = (A[l][l] - A[k][k]) / (2*A[k][l]);
        t = -tau + sqrt(1 + pow(tau, 2.0));
        c = 1/sqrt(1 + pow(t, 2.0));
        s = t*c;

        B[k][k] = A[k][k]*c*c - 2*A[k][l]*c*s + A[l][l]*s*s;
        B[l][l] = A[l][l]*c*c + 2*A[k][l]*c*s + A[k][k]*s*s;
        B[k][l] = (A[k][k] - A[l][l])*c*s + A[k][l]*(c*c - s*s);
        for (int i = 0; i < N+1; i++) {
            if ((i != k) && (i != l)) {
				B[i][i] = A[i][i];
				B[i][k] = A[i][k]*c - A[i][l]*s;
				B[i][l] = A[i][l]*c + A[i][k]*s;
			}
		}
        for (int i = 0; i < N-1; i++) {
            for (int j = 0; j < N-1; j++) {
                A[i][j] = B[i][j];
            }
        }
        count++;
    }

	if (!(minimal_print)) {
        cout << "Sim. transf. necessary when n=" << N << ": " << count << endl;
    }
    if ((N < 10) && !(minimal_print)) {
		for (int i = 0; i < N-1; i++) {
		    for (int j = 0; j < N-1; j++) {
		        cout << setw(12) << A[i][j] << " ";
		    }
		    cout << endl;
		}
	}
    delete[] A;
}

void jacobi_buckling_beam(int N, bool minimal_print, int max_iter) {
    jacobi(make_buckling_beam_matrix(N), N, minimal_print, max_iter);
}

void jacobi_quantum_dots_one_electron(int N, bool minimal_print, int max_iter) {
    jacobi(make_quantum_dots_one_electron_matrix(N), N, minimal_print, max_iter);
}

void jacobi_quantum_dots_two_electrons(int N, double omega_r, bool minimal_print, int max_iter) {
    jacobi(make_quantum_dots_two_electrons_matrix(N, omega_r), N, minimal_print, max_iter);
}
