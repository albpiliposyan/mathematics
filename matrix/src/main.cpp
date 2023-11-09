#include <iostream>
#include "matrix.h"

// TODO C++20 modules instead of headers

int main() {


    const int n = 3;
    const double eps = pow(10, -14);

    Matrix<double> A(n, n);
    std::cout << "Enter A matrix:" << std::endl;
    A.initialize();

    std::vector<double> b(n);
    std::cout << "Enter b vector:" << std::endl;
    for (int i = 0; i < b.size(); ++i) {
        std::cin >> b[i];
    }

    std::vector<double> solution = A.jacobi_method(b, eps);

    std::cout << "\n\nSolution:" << std::endl;
    for (int i = 0; i < solution.size(); ++i) {
        std::cout << solution[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
