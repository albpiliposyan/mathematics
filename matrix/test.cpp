#include <iostream>
#include "include/matrix.h"


void test_gauss(TwoDVector<double> input_vec, TwoDVector<double> test_vec) {
	Matrix<double> input(input_vec); Matrix<double> test(test_vec);

	Matrix<double> inverse = input.inverse_gauss();

	std::cout << ">>Inverse (Gauss)\ninput:" << std::endl;
	inverse.print();
	std::cout << "test:" << std::endl;
	test.print();
	std::cout << std::endl;
}


void test_lu(TwoDVector<double> input_vec, TwoDVector<double> test_vec) {
	Matrix<double> input(input_vec);
	Matrix<double> test(test_vec);


	Matrix<double> inverse = input.inverse_lu();

	std::cout << ">>Inverse (LU)\ninput:" << std::endl;
	inverse.print();
	std::cout << "test:" << std::endl;
	test.print();
	std::cout << std::endl;
}

void test_adjustment(TwoDVector<double> input_vec, TwoDVector<double> adj_vec, const int iter_count, TwoDVector<double> test_vec) {
	Matrix<double> input(input_vec);
	Matrix<double> adj(adj_vec);
	Matrix<double> test(test_vec);


	Matrix<double> inverse = input.inverse_adjustment(adj, iter_count);

	std::cout << ">>Inverse (Adjustment)\nviolated:" << std::endl;
	adj.print();
	std::cout << "inverse:" << std::endl;
	inverse.print();
	std::cout << "test:" << std::endl;
	test.print();
	std::cout << std::endl;
}

int main() {
	test_gauss({{1, 1, 1}, {1, 2, 4}, {1, 3, 9}},
			{{3, -3, 1}, {-2.5, 4, -1.5}, {0.5, -1, 0.5}});

	test_gauss({{1, 2, 3}, {1, 3, 5}, {1, -1, 0}},
			{{1.6666667, -1, 0.3333333}, {1.6666667, 1, 0.3333333}, {1.3333337, -1, 0.3333333}});

	// result = test_gauss({{4, 0, 0}, {0, 6, 0}, {0, 0, 1}},
	// 		{{3, -3, 1}, {-2.5, 4, -1.5}, {0.5, -1, 0.5}});

	test_lu({{1, 1, 1}, {1, 2, 4}, {1, 3, 9}},
			{{3, -3, 1}, {-2.5, 4, -1.5}, {0.5, -1, 0.5}});

	test_lu({{4, 0, 0}, {0, 6, 0}, {0, 0, 1}},
			{{0.25, 0, 0}, {0, 0.1666667, 0}, {0, 0, 1}});


	const int iter_count = 10;
	test_adjustment({{1, 1, 1}, {1, 2, 4}, {1, 3, 9}},
					{{3.001, -2.9, 1.2}, {-2.50001, 3.907, -1.5001}, {0.501, -1.0001, 0.52}},
					iter_count,
					{{3, -3, 1}, {-2.5, 4, -1.5}, {0.5, -1, 0.5}});

	test_adjustment({{4, 0, 0}, {0, 6, 0}, {0, 0, 1}},
					{{0.25, -0.008, 0.001}, {0.0001, 0.15, 0.0005}, {-0.004, 0.09, 1.05}},
					iter_count,
					{{0.25, 0, 0}, {0, 0.1666667, 0}, {0, 0, 1}});

    return 0;
}













