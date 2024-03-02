#include <iostream>
#include "matrix.h"


int main() {
    std::cout << std::fixed;

    int n;
    std::cout << "Enter the size of the square matrix:\n";
    std::cin >> n;
    std::cout << "\n";

    Matrix<double> A(n, n);
    std::cout << "Enter A matrix:" << std::endl;
    A.insert();
    std::cout << "\n";

    std::vector<Matrix<double>> qr = A.qr_decomposition();
	std::cout << "Q matrix:" << std::endl;
    qr[0].print();
    std::cout << "\n";

	std::cout << "R matrix:" << std::endl;
    qr[1].print();
    std::cout << "\n";

    std::cout << "Check if A = QR:\n";
    (qr[0] * qr[1]).print();

	std::cout << "Check if Q columns are orthogonal:" << std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dp = dot_product(qr[0].get_matrix()[i], qr[0].get_matrix()[j]);
            std::cout << "(" << i << ", " << j << ") = " <<  dp << "\n";
        }
    }

    return 0;
}

