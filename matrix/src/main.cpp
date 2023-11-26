#include <iostream>
#include "matrix.h"


int main() {

    const int n = 3;
    Matrix<double> A(n, n);
    std::cout << "Enter A matrix:" << std::endl;
    A.insert();
    
    Matrix<double> A_inv = A.inverse_lu();
    std::cout << "\nA inverse is: \n";
    A_inv.print();

    return 0;
}
