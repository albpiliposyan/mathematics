#include <iostream>
#include "matrix.h"

// TODO C++20 modules instead of headers

int main() {

    // // matrix.print();
    Matrix<double> a(3,3);
    std::cout << "Enter 3x3 A:" << std::endl;
    a.initialize();

    Matrix<double> u = a.cholesky_decomposition();
    std::cout << std::endl;

    std::cout << "U:" << std::endl;
    u.print();
    std::cout << std::endl;

    return 0;
}
