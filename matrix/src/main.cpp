#include <iostream>
#include "matrix.h"
    // TODO C++20 modules instead of headers

int main() {
    Matrix<double> matrix(4, 5);
    matrix.initialize();
    std::cout << std::endl;
    // std::cout << matrix.determinant() << std::endl;
    
    std::vector<double> solution = matrix.gaussEliminationMethod();
    
    if (solution.size() != 0) {
        std::cout << "The solutions to the system of linear equations are: " << std::endl;
        for (double& element : solution) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }

    // matrix.print();
    return 0;
}
