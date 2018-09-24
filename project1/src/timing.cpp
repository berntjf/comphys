#include <chrono>
#include "head.h"


auto time_start() {
    auto start = chrono::high_resolution_clock::now();
    return start;
}

int time_finish(auto start) {
    auto finish = chrono::high_resolution_clock::now();
    cout << chrono::duration_cast<chrono::microseconds>(finish-start).count() << " mcs\n";

    return 0;
}
