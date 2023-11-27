#include <iostream>
#include "matrix.h"


int main() {

    const int n = 3;
    Matrix<double> A(n, n);
    std::cout << "Enter A matrix:" << std::endl;
    A.insert();


	Matrix<double> X_0(n, n);
	std::cout << "Enter X_0 approximation of A matrix:" << std::endl;
	X_0.insert();

	Matrix<double> A_inv = A.inverse_adjustment(X_0, 15);
	std::cout << "\n\n\nResult:\n";
	A_inv.print();


    // std::cout << "\n";
	// Matrix<double> A_inv = A.inverse_lu();
    // A_inv.print();

    return 0;
}

