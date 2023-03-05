#include <iostream>
#include "abelian_group.h"
#include "field.h"
#include <cassert>


int main() {

    // Checking the Group operations
    Group<int> a{1}, b{2};
    Group<int> c = a + b;
    // associativity
    assert((a + b) + c == a + (b + c));
    // identity
    assert(a + Group<int>(0) == a);
    // inverse
    assert(a + (-a) == Group<int>(0));


    // Checking the Ring operations
    Ring<int> d{3}, e{4};
    Ring<int> f = d * e;
    // associativity
    assert((d * e) * f == d * (e * f));
    // identity
    assert(d * Ring<int>(1) == d);
    // commutativity
    assert(d * e == e * d);
    // inverse
    assert(d * (e + f) == d*e + d*f);
    assert((d + e) * f == d*f + e*f);
    

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



