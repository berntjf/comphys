#include <string>
#include <stdarg.h>
#include "head.h"

int make_table_file(string filename, int length_of_arrays, int number_of_arrays, vec* a...) {
    string fileout = "out/" + filename + ".dat";
    ofstream ofile;
    ofile.open(fileout);
    for (int i = 0; i < length_of_arrays; i++ ) {
        for (int j = 0; j < number_of_arrays; j++) {
            ofile << a[j][i] << " ";
        }
        ofile << endl;
    }
    ofile.close();

    return 0;
}

/*
int main() {
    double* a = new double[4];
    double* b = new double[4];
    a[0] = 0;
    a[1] = 1;
    a[2] = 2;
    a[3] = 3;
    b[0] = 4;
    b[1] = 5;
    b[2] = 6;
    b[3] = 7;
    double** ab = new double*[2];
    ab[0] = a;
    ab[1] = b;
    make_table_file("TEST", 4, 2, ab);

    return 0;
}
*/
