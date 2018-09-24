#include <iostream>
#include <armadillo>
using namespace arma;
using namespace std;

int main()
  {
   mat A = randu<mat>(5,5);
   vec b = randu<vec>(5);

  A.print("A =");
  b.print("b=");
  // solve Ax = b
  vec x = solve(A,b);
  // print x
  x.print("x=");
  // find LU decomp of A, if needed, P is the permutation matrix
  mat L, U;
  lu(L,U,A);
  // print l
  L.print(" L= ");
  // print U
  U.print(" U= ");
  //Check that A = LU
  (A-L*U).print("Test of LU decomposition");
    return 0;
  }
