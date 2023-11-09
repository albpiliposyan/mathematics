#include <iostream>
#include "abelian_group.h"
#include "field.h"
#include <cassert>


int main() {

    // Checking the Group operations
    Group<int> a{1}, b{2};
    Group<int> c = a + b;

    // Checking the Ring operations
    Ring<int> d{3}, e{4};
    Ring<int> f = d * e;

    // Field operations
    std::cout << "Field operations:" << std::endl;
    Field<int> k(5), l(6);
    std::cout << "k = " << k << "\nl = " << l << std::endl;
    std::cout << "k + l = " << k + l << std::endl;
    std::cout << "k - l = " << k - l << std::endl;
    std::cout << "k * l = " << k * l << std::endl;
    std::cout << "k / l = " << k / l << std::endl;

    return 0;
}



