#ifndef UTILS_H
#define UTILS_H

// GCD of two numbers
int gcd(int, int);
// input the equations in a vector<int, int>

std::vector<std::pair<int, int>> input_equations();
// print the equations

bool print_pairs_int(std::vector<std::pair<int, int>> );
// inverse of  a (mod m)
int modular_inverse(int, int);

// chinese remainder theorem algorithm
std::pair<int, int> chinese(std::vector<int>, std::vector<int>);


#endif
