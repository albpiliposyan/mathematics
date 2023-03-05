// TUTORIAL: https://www.youtube.com/watch?v=zIFehsBHB8o&t=362s

#include <iostream>
#include <vector>
#include "utils.h"

int main() {

	std::vector<std::pair<int, int>> equations = input_equations();
	print_pairs_int(equations);

	std::vector<int> A, M;
	for (int i = 0; i < equations.size(); i++) {
		A.push_back(equations[i].first);
		M.push_back(equations[i].second);
	}

	std::pair<int, int> result = chinese(A, M);
	std::cout << "\nThe result: \n" <<
		"x = " << result.first << " (mod " << result.second << ")" << std::endl;

	//std::cout << modular_inverse(6, 13) << std::endl;

	return 0;
}
