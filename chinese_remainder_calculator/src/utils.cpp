#include <iostream>
#include <vector>

// GCD of two numbers
int gcd(int a, int b) {
	if (a == b)
		return a;
	if (a > b)
		return gcd(a - b, b);
	else
		return gcd(a, b - a);
}

// input the equations in a vector<int, int>
std::vector<std::pair<int, int>> input_equations() {
	int eq_number;
	std::cout << "Enter the number of equations: ";
	std::cin >> eq_number;
	std::vector<std::pair <int, int> > equations;
	int a, b;
	for (int i = 0; i < eq_number; ++i) {
		std::cin >> a >> b;
		a = a % b;						// a % b
		equations.push_back(std::make_pair(a, b));
	}

	for (int i = 0; i < equations.size(); ++i)
		for (int j = i + 1; j < equations.size(); ++j)
			if (gcd(equations[i].second, equations[j].second) != 1) {
				std::cerr << "The M numbers should be relatively prime." << std::endl;
				exit(1);
			}
	return equations;
}

// print the equations
bool print_pairs_int(std::vector<std::pair<int, int>> equations) {
	std::cout << std::endl;
	for (int i = 0; i < equations.size(); ++i) {
		std::cout << "x = " << equations[i].first<< " (mod " << equations[i].second << ")" << std::endl;
	}
	return false;
}

// inverse of  a (mod m)
int modular_inverse(int a, int m) {
	a = a % m;
	for (int i = 1; i < m; ++i)
		if (a*i % m == 1)
			return i;
	return -1;
}

// chinese remainder theorem algorithm
std::pair<int, int> chinese(std::vector<int> A, std::vector<int> M) {
	int size = A.size();
	std::vector<int> N,
	   	N_inverse;

	int N_prod = 1; 	// new modulo
	for (int i = 0; i < size; ++i)
		N_prod *= M[i];

	for (int i = 0; i < size; ++i) {
		float tmp = N_prod / M[i];
		N.push_back(tmp);
	}

	for (int i = 0; i < size; ++i)
		N_inverse.push_back(modular_inverse(N[i], M[i]));

	int x = 0;
	for (int i = 0; i < size; ++i)
		x += A[i] * N[i] * N_inverse[i];
	x = x % N_prod;

	return std::make_pair(x, N_prod);
}


