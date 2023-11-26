#include <iostream>
#include "matrix.h"


int main() {

    const int n = 3;
    Matrix<double> A(n, n);
    std::cout << "Enter A matrix:" << std::endl;
    A.insert();

    std::cout << "\n";
	Matrix<double> A_inv = A.inverse();
    A_inv.print();

    return 0;
}

