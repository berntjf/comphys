// Project 2b
#include "head.h"
#include "toeplitz.cpp"

int main() {
    int max_iter = 200;
    jacobi_quantum_dots_two_electrons(8,  0.01, false,  max_iter);
	jacobi_quantum_dots_two_electrons(8,  0.50, false,  max_iter);
	jacobi_quantum_dots_two_electrons(8,  1.00, false,  max_iter);
	jacobi_quantum_dots_two_electrons(8,  5.00, false,  max_iter);

    return 0;
}
