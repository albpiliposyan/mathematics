#include <iostream>
#include "complex.h"


int main() {
    const Complex<double> a(1, 2),  b(2, 3);
    Complex<double> c = (a + b);
    std::cout << "1+2i + 2+3i = ";
    std::cout << c << std::endl;
    std::cout << "1+2i - 2+3i = ";
    c = a - b;
    std::cout << c << std::endl;
    std::cout << "1+2i * 2+3i = ";
    c = a * b;
    std::cout << c << std::endl;
    c = a / b;
    std::cout << "1+2i / 2+3i = ";
    std::cout << c << std::endl;
    c = a.conjugate();
    std::cout << c << std::endl;

    std::cout << "Norm b: " << b.magnitude() << std::endl;
    return 0;
}
