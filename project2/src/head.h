#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <armadillo>
using namespace std;

double lambda_anal(int j, int N, double a, double d);
int* find_kl_when_a2_max(double** A, int N, double epsilon);
double** make_buckling_beam_matrix(int N);
void jacobi(double** A, int N, bool minimal_print, int max_iter);
void jacobi_buckling_beam(int N, bool minimal_print, int max_iter);
void jacobi_quantum_dots_one_electron(int N, bool minimal_print, int max_iter);
int test_lambda_anal();
int test_find_kl_when_a2_max();
double** make_quantum_dots_two_electrons_matrix(int N, double omega_r);
