// Project 2b
#include "head.h"
#include "toeplitz.cpp"

int main() {
    int max_iter = 200;
    jacobi_quantum_dots_one_electron(6,   false, max_iter);
	jacobi_quantum_dots_one_electron(11,  false,  max_iter);
	jacobi_quantum_dots_one_electron(101, true,  max_iter);

    return 0;
}
