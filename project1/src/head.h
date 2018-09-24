#include <iostream> // c-out etc.
#include <fstream> // // file handling
#include <string>
#include <cmath>
#include <armadillo>
using namespace std;

double* make_abc(int , float );
arma::mat make_A(int );
double f(double );
double u_anal(double );
double error_max(int , double* , double* );
arma::vec rref_tridiagonal(int , double* , double* , double* , double , double* );
arma::vec rref_tridiagonal_3num(int , double , double , double , double , double* );
auto time_start();
int time_finish(auto start);
