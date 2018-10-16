// Project 2b
#include "head.h"
#include "toeplitz.cpp"

int main() {
    int max_iter = 200;
    jacobi_buckling_beam(6,   false, max_iter);
	jacobi_buckling_beam(11,  false,  max_iter);
	jacobi_buckling_beam(101, true,  max_iter);

    return 0;
}
