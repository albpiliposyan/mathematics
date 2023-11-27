#ifndef MATRIX_H
#define MATRIX_H


#include <numeric>
#include <vector>
#include <stdexcept>
#include <cmath>


template <class T>
concept Field = std::regular<T> && requires(T x, T y) {
    {x + y} -> std::convertible_to<T>;
    {x - y} -> std::convertible_to<T>;
    {x * y} -> std::convertible_to<T>;
    {x / y} -> std::convertible_to<T>;
};


template <Field T>
class Matrix;


template <Field T>
using TwoDVector = std::vector<std::vector<T>>;


template <Field U>
Matrix<U> operator+(Matrix<U> lhs, const Matrix<U>& rhs) {
    lhs += rhs;
    return lhs;
}


template <Field U>
 Matrix<U> operator-(Matrix<U> lhs, const Matrix<U>& rhs) {
    lhs -= rhs;
    return lhs;
}


template <Field U>
 Matrix<U> operator*(Matrix<U> lhs, const Matrix<U>& rhs) {
    lhs *= rhs;
    return lhs;
}



template <Field T>
class Matrix
{
public:
    Matrix() = default;
    Matrix(int rows, int cols) : m_rows(rows), m_cols(cols),  m_matrix(m_rows, std::vector<T>(m_cols)) {}
    Matrix(TwoDVector<T> mat) : m_rows(mat.size()), m_cols(mat[0].size()), m_matrix(mat) {}
    Matrix(const Matrix& other) : m_rows(other.m_rows), m_cols(other.m_cols), m_matrix(other.m_matrix) {}
    Matrix(const Matrix&& other) : m_rows(std::move(other.m_rows)), m_cols(std::move(other.m_cols)), m_matrix(std::move(other.m_matrix)) {}
public:
    constexpr bool insert();
    constexpr bool print() noexcept;
	constexpr Matrix<T> get_identity();
    constexpr T determinant() const;
    Matrix<T>& transpose() noexcept;
    Matrix<T>& swap_rows(const unsigned int& row1, const unsigned int& row2);
    Matrix<T>& swap_columns(const unsigned int& col1, const unsigned int& col2);
    Matrix<T>& multiply(const T& scalar);
    Matrix<T>& multiply_row(const unsigned int& row, const T& scalar);
    Matrix<T>& add_multiple_of_row(const unsigned int& row1,
                        const unsigned int& row2, const T& scalar);
    constexpr bool check_zero_row(const unsigned int row) const noexcept;
    constexpr bool check_zero_column(const unsigned int col) const noexcept;
    constexpr unsigned int number_of_zero_rows() const noexcept;
    constexpr unsigned int number_of_zero_columns() const noexcept;
    constexpr bool swap_non_zero_rows_to_top() noexcept;
    constexpr bool swap_non_zero_columns_to_left() noexcept;
    Matrix<T>& row_reduced_form();
    Matrix<T>& row_reduced_echelon_form();
    constexpr Matrix<T> inverse_lu();
	constexpr Matrix<T> inverse_gauss();
	constexpr Matrix<T> inverse_adjustment(Matrix<T> X_prev, const unsigned int num_of_iterations);
	constexpr Matrix<T> inverse();
    constexpr std::vector<T> gauss_elimination_method();
    Matrix<T> cholesky_decomposition();
    constexpr std::vector<Matrix<T>> crout_decomposition();
    constexpr std::vector<T> jacobi_method(std::vector<T>, const double);
    constexpr std::vector<T> lower_triang_equation(const std::vector<T>);
    constexpr std::vector<T> upper_triang_equation(const std::vector<T>);
    constexpr std::vector<T> system_of_equations_lu(const std::vector<T> b);
public:
    constexpr Matrix<T>& operator=(const Matrix<T>& other) noexcept;
    constexpr Matrix<T>& operator=(const Matrix<T>&& other) noexcept;
    constexpr Matrix<T>& operator+=(const Matrix<T>& rhs);
    constexpr Matrix<T>& operator-=(const Matrix<T>& rhs);
    constexpr Matrix<T>& operator*=(const Matrix<T>& rhs);
    constexpr std::vector<T>& operator[](int row) noexcept;
    const constexpr std::vector<T>& operator[](int row) const noexcept;
public:
    template <Field U>
    friend Matrix<U> operator+(Matrix<U> lhs, const Matrix<U>& rhs);
    template <Field U>
    friend Matrix<U> operator-(Matrix<U> lhs, const Matrix<U>& rhs);
    template <Field U>
    friend Matrix<U> operator*(Matrix<U> lhs, const Matrix<U>& rhs);
public:
    constexpr TwoDVector<T>& get_matrix() noexcept;
    constexpr const TwoDVector<T> get_matrix() const noexcept;
    constexpr bool set_matrix(TwoDVector<T> other) noexcept;
    constexpr unsigned int rows() const noexcept;
    constexpr unsigned int cols() const noexcept;
    constexpr bool set_rows(unsigned int row) noexcept;
    constexpr bool set_cols(unsigned int col) noexcept;
private:
    constexpr T cofactor(const unsigned int& row, const unsigned int& col) const;
    constexpr double vector_max_diff(const std::vector<T>, const std::vector<T>);
private:
    unsigned int m_rows;
    unsigned int m_cols;
    TwoDVector<T> m_matrix;
};


template <Field T>
constexpr bool Matrix<T>::insert() {
    for (unsigned int i = 0; i < rows(); ++i) {
        for (unsigned int j = 0; j < cols(); ++j) {
            std::cin >> get_matrix()[i][j];
        }
    }
    return true;
}


template <Field T>
constexpr bool Matrix<T>::print() noexcept {
    for (const std::vector<T>& row : get_matrix()) {
        for (const T& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
    return true;
}


template <Field T>
constexpr Matrix<T> Matrix<T>::get_identity() {
	Matrix<T> tmp(rows(), rows());
	for (unsigned int i = 0; i < rows(); ++i) {
		tmp.get_matrix()[i][i] = 1;
	}
	return tmp;
}

template <Field T>
constexpr T Matrix<T>::determinant() const {
    if (rows() != cols()) {
        throw std::invalid_argument("Matrix must be square");
    }
    T det = T();
    if (rows() == 1) {
        return get_matrix()[0][0];

    }
    if (rows() == 2) {
        return get_matrix()[0][0] * get_matrix()[1][1] -
            get_matrix()[0][1] * get_matrix()[1][0];
    }

    for (unsigned int i = 0; i < cols(); ++i) {
        det += get_matrix()[0][i] * cofactor(0, i);
    }

    return det;
}


template <Field T>
constexpr T Matrix<T>::cofactor(const unsigned int& row, const unsigned int& col) const {
    Matrix<T> submatrix(rows() - 1, cols() - 1);
    int subi = 0, subj = 0;
    for (unsigned int i = 0; i < rows(); ++i) {
        if (i == row) {
            continue;
        }
        subj = 0;
        for (unsigned int j = 0; j < cols(); ++j) {
            if (j == col) {
                continue;
            }
            submatrix[subi][subj] = get_matrix()[i][j];
            ++subj;
        }
        ++subi;
    }
    T det = submatrix.determinant();
    return ((row + col) % 2 == 0 ? 1 : -1) * det;
}


template <Field T>
Matrix<T>& Matrix<T>::transpose() noexcept {
    TwoDVector<T> transposed(cols(), std::vector<T>(rows()));
    for (unsigned int i = 0; i < rows(); ++i) {
        for (unsigned int j = 0; j < cols(); ++j) {
            transposed[j][i] = get_matrix()[i][j];
        }
    }
    set_matrix(std::move(transposed));
    return *this;
}


template <Field T>
Matrix<T>& Matrix<T>::swap_rows(const unsigned int& row1, const unsigned int& row2) {
    if (row1 < 0 || row1 >= rows() || row2 < 0 || row2 >= rows())
        throw std::out_of_range("Invalid row index");
    std::swap(get_matrix()[row1], get_matrix()[row2]);
    return *this;
}


template <Field T>
Matrix<T>& Matrix<T>::swap_columns(const unsigned int& col1, const unsigned int& col2) {
    if (col1 < 0 || col1 >= cols() || col2 < 0 || col2 >= cols())
        throw std::out_of_range("Invalid column index");
    for (unsigned int i = 0; i < rows(); ++i)
        std::swap(m_matrix[i][col1], m_matrix[i][col2]);
    return *this;
}


template <Field T>
Matrix<T>& Matrix<T>::multiply(const T& scalar) {
    if (scalar == 0) {
        throw std::invalid_argument("The scalar must be non-zero");
    }
    for (unsigned int i = 0; i < rows(); ++i) {
        for (unsigned int j = 0; j < cols(); ++j) {
            get_matrix()[i][j] *= scalar;
        }
    }
    return *this;
}


template <Field T>
Matrix<T>& Matrix<T>::multiply_row(const unsigned int& row, const T& scalar) {
    if (scalar == 0) {
        throw std::invalid_argument("The scalar must be non-zero");
    }
    if (row >= rows()) {
        throw std::invalid_argument("Invalid row index");
    }
    for (unsigned int i = 0; i < cols(); ++i) {
        get_matrix()[row][i] *= scalar;
    }
    return *this;
}


template <Field T>
Matrix<T>& Matrix<T>::add_multiple_of_row(const unsigned int& row1,
                    const unsigned int& row2, const T& scalar) {
    if (scalar == 0) {
        throw std::invalid_argument("The scalar must be non-zero");
    }
    if (row1 >= rows() || row2 >= rows()) {
        throw std::invalid_argument("Invalid row index");
    }
    for (unsigned int i = 0; i < cols(); ++i) {
        get_matrix()[row1][i] += scalar * get_matrix()[row2][i];
    }
    return *this;
}


template <Field T>
constexpr bool Matrix<T>::check_zero_row(const unsigned int row) const noexcept {
    for (unsigned int i = 0; i < cols(); ++i) {
        if (get_matrix()[row][i] != 0) {
            return false;
        }
    }
    return true;
}


template <Field T>
constexpr bool Matrix<T>::check_zero_column(const unsigned int col) const noexcept{
    for (unsigned int i = 0; i < rows(); ++i) {
        if (get_matrix()[i][col] != 0) {
            return false;
        }
    }
    return true;
}


template <Field T>
constexpr bool Matrix<T>::swap_non_zero_rows_to_top() noexcept {
    unsigned int row_index = rows() - 1;
    for (unsigned int i = 0; i <= row_index; ++i) {
        if (check_zero_row(i)) {
            swap_rows(i--, row_index--);
        }
    }
    return true;
}


template <Field T>
constexpr bool Matrix<T>::swap_non_zero_columns_to_left() noexcept {
    unsigned int col_index = 0;
    for (unsigned int i = cols() - 1; i > col_index; --i) {
        if (check_zero_column(i)) {
            swap_columns(i++, col_index++);
        }
    }
    return true;
}


template <Field T>
constexpr unsigned int Matrix<T>::number_of_zero_rows() const noexcept {
    unsigned int count = 0, row_zeroes;
    for (unsigned int i = 0; i < rows(); ++i) {
        row_zeroes = 0;
        for (unsigned int j = 0; j < cols(); ++j) {
            if (get_matrix()[i][j] == 0) {
                ++row_zeroes;
            }
        }
        if (row_zeroes == cols()) {
            ++count;
        }
    }
    return count;
}


template <Field T>
constexpr unsigned int Matrix<T>::number_of_zero_columns() const noexcept {
    unsigned int count = 0, col_zeroes;
    for (unsigned int i = 0; i < cols(); ++i) {
        col_zeroes = 0;
        for (unsigned int j = 0; j < rows(); ++j) {
            if (get_matrix()[j][i] == 0) {
                ++col_zeroes;
            }
        }
        if (col_zeroes == rows()) {
            ++count;
        }
    }
    return count;
}


template <Field T>
Matrix<T>& Matrix<T>::row_reduced_form() {
    unsigned int zero_rows_count = number_of_zero_rows();
    if (zero_rows_count == rows()) {
        return *this;
    }
    unsigned int zero_cols_count = number_of_zero_columns();
    if (zero_rows_count > 0) {
        swap_non_zero_rows_to_top();
    }
    if (zero_cols_count > 0) {
        swap_non_zero_columns_to_left();
    }

    for (unsigned int i = 0; i < rows() - zero_rows_count; ++i) {
        for (unsigned int j = i + 1; j < rows() - zero_rows_count; ++j) {
            if (std::abs(get_matrix()[i][zero_cols_count + i]) <
                                    std::abs(get_matrix()[j][zero_cols_count + i])) {
                swap_rows(i, j);
            }
        }
    }
    for (unsigned int i = 0; i < rows() - zero_rows_count - 1; ++i) {
        for (unsigned int j = i + 1; j < rows() - zero_rows_count; ++j) {
            if (get_matrix()[i][zero_cols_count + i] != 0) {
                add_multiple_of_row(j, i, (-1) * get_matrix()[j][zero_cols_count + i] /
                                        get_matrix()[i][zero_cols_count + i]);
            }
        }
    }
    swap_non_zero_rows_to_top();
    swap_non_zero_columns_to_left();
    return *this;
}

template <Field T>
Matrix<T>& Matrix<T>::row_reduced_echelon_form() {
	row_reduced_form();
    for (int i = 0; i < rows(); ++i) {
		multiply_row(i, 1 / get_matrix()[i][i]);
    }
    return *this;
}

template <Field T>
constexpr Matrix<T> Matrix<T>::inverse_gauss() {
	if (determinant() == 0) {
		throw std::invalid_argument("The matrix is not invertible.");
	}

	Matrix<T> tmp(rows(), 2 * rows());
	for (int i = 0; i < rows(); ++i) {
		for (int j = 0; j < cols(); ++j) {
			tmp.get_matrix()[i][j] = get_matrix()[i][j];
		}
	}
	for (int i = 0; i < tmp.rows(); ++i) {
		tmp.get_matrix()[i][tmp.rows() + i] = 1;
	}

	tmp.row_reduced_echelon_form();

	for (int i = 0; i < tmp.rows(); ++i) {
		for (int j = i; j < tmp.rows(); ++j) {
			if (i != j) {
				tmp.add_multiple_of_row(i, j, (-1) * tmp.get_matrix()[i][j]
						/ tmp.get_matrix()[i][i]);
			}
		}
	}
	Matrix<T> inverse(tmp.rows(), tmp.rows());
	for (int i = 0; i < inverse.rows(); ++i) {
		for (int j = 0; j < inverse.cols(); ++j) {
			inverse.get_matrix()[i][j] = tmp.get_matrix()[i][tmp.rows() + j];
		}
	}

	return inverse;
}


template <Field T>
constexpr Matrix<T> Matrix<T>::inverse() {
	return inverse_gauss();
}


template <Field T>
constexpr Matrix<T> Matrix<T>:: inverse_adjustment(Matrix<T> X_prev, const unsigned int num_of_iterations) {
	const unsigned int n = X_prev.rows();

	Matrix<T> F_prev(n, n);
	Matrix<T> X_k = X_prev;

	for (int k = 0; k < num_of_iterations; ++k) {
		X_prev = X_k;

		Matrix<T> prod = (*this) * X_prev;
		Matrix<T> identity = get_identity();
		;
		F_prev = identity - prod;

		F_prev += identity;

		X_k = X_prev * F_prev;

		std::cout << std::endl;
		X_k.print();
	}


	return X_k;
}


template <Field T>
constexpr std::vector<T> Matrix<T>::gauss_elimination_method() {
    Matrix tmp(*this);
    tmp.row_reduced_form();
    unsigned int zero_rows_count = tmp.number_of_zero_rows();
    unsigned int zero_cols_count = tmp.number_of_zero_columns();

    unsigned int last_row_check = tmp.rows() - zero_rows_count - 1;
    if (tmp.get_matrix()[last_row_check][tmp.cols() - 2] == 0) {
        std::cerr << "There is no solution for this system of equations." << std::endl;
        return std::vector<T>();
    }
    if (tmp.rows() - zero_rows_count + 1 != tmp.cols() - zero_cols_count) {
        std::cerr << "There are infinite number of solution to this system of equations" << std::endl;
        return std::vector<T>();
    }

    std::vector<T> solution;
    solution.resize(cols() - zero_cols_count - 1);

    for (unsigned int i = tmp.rows() - zero_rows_count - 1; i >= 0 && i <= tmp.rows(); --i) {
        solution[i] = tmp.get_matrix()[i][zero_cols_count + solution.size()];
        for (unsigned int j = i + 1; j < solution.size(); ++j) {
            if (i != j) {
                solution[i] -=  solution[j] * tmp.get_matrix()[i][zero_cols_count + j];
            }
        }
        solution[i] /= tmp.get_matrix()[i][zero_cols_count + i];
    }
    return solution;
}


template <Field T>
Matrix<T> Matrix<T>::cholesky_decomposition() {
    if (rows() != cols()) {
        throw std::invalid_argument("The matrix must be a square matrix for cholesky decomposition.");
    }
    const int n = rows();
    Matrix<T> U(n, n);

    U.get_matrix()[0][0] = sqrt(get_matrix()[0][0]);

    for (int i = 1; i < n; ++i) {
        U.get_matrix()[0][i] = get_matrix()[0][i] / U.get_matrix()[0][0];
    }

    double tmp_sum;
    for (int i = 1; i < n; ++i) {
        tmp_sum = 0;
        for (int k = 0; k < i; ++k) {
            tmp_sum += U.get_matrix()[k][i] * U.get_matrix()[k][i];
        }
        U.get_matrix()[i][i] = sqrt(get_matrix()[i][i] - tmp_sum);

        for (int j = i + 1; j < n; ++j) {
            tmp_sum = 0;
            for (int k = 0; k < i; ++k) {
                tmp_sum += U.get_matrix()[k][i] * U.get_matrix()[k][j];
            }
            U.get_matrix()[i][j] = (get_matrix()[i][j] - tmp_sum) / U.get_matrix()[i][i];
        }
    }

    return U;
}


template <Field T>
constexpr Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) noexcept {
    if (this != &other) {
        m_matrix = other.m_matrix;
    }
    return *this;
}


template <Field T>
constexpr Matrix<T>& Matrix<T>::operator=(const Matrix<T>&& other) noexcept {
    if (this != &other) {
        m_matrix = other.m_matrix;
    }
    return *this;
}


template <Field T>
constexpr Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) {
    if (rows() != rhs.rows() || cols() != rhs.cols())
        throw std::invalid_argument("Matrix sizes don't match");
    for (unsigned int i = 0; i < rows(); ++i) {
        for (unsigned int j = 0; j < cols(); ++j) {
            get_matrix()[i][j] += rhs[i][j];
        }
    }
    return *this;
}


template <Field T>
constexpr Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) {
    if (rows() != rhs.rows() || cols() != rhs.cols())
        throw std::invalid_argument("Matrix sizes don't match");
    for (unsigned int i = 0; i < rows(); ++i) {
        for (unsigned int j = 0; j < cols(); ++j) {
            get_matrix()[i][j] -= rhs[i][j];
        }
    }
    return *this;
}


template <Field T>
constexpr Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs) {
    if (rows() != rhs.cols())
        throw std::invalid_argument("Matrix sizes don't match");
    Matrix<T> result(rows(), rhs.cols());
    for (unsigned int i = 0; i < result.rows(); ++i)
        for (unsigned int j = 0; j < result.cols(); ++j)
            for (unsigned int k = 0; k < cols(); ++k)
                result[i][j] += get_matrix()[i][k] * rhs[k][j];
    // set_matrix(result);
    *this = result;
    return *this;
}


template <Field T>
constexpr std::vector<T>& Matrix<T>::operator[](int row) noexcept {
    return m_matrix[row];
}


template <Field T>
const constexpr std::vector<T>& Matrix<T>::operator[](int row) const noexcept {
    return m_matrix[row];
}


template <Field T>
constexpr TwoDVector<T>& Matrix<T>::get_matrix() noexcept {
    return m_matrix;
}


template <Field T>
constexpr const TwoDVector<T> Matrix<T>::get_matrix() const noexcept {
    return m_matrix;
}


template <Field T>
constexpr bool Matrix<T>::set_matrix(TwoDVector<T> other)  noexcept {
    m_matrix = std::move(other);
    return true;
}


template <Field T>
constexpr unsigned int Matrix<T>::rows() const noexcept {
    return m_rows;
}


template <Field T>
constexpr unsigned int Matrix<T>::cols() const noexcept {
    return m_cols;
}


template <Field T>
constexpr bool Matrix<T>::set_rows(unsigned int row) noexcept {
    get_matrix().resize(row);
    return true;
}


template <Field T>
constexpr bool Matrix<T>::set_cols(unsigned int col) noexcept {
    for (auto& row : get_matrix()) {
        row.resize(col);
    }
    return true;
}

template <Field T>
constexpr std::vector<Matrix<T>> Matrix<T>::crout_decomposition() {
    int n = rows();
    Matrix<T> L(n, n);
    Matrix<T> U(n, n);

    double sum = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += L[i][k] * U[k][j];
            }
            L[i][j] = get_matrix()[i][j] - sum;
        }
        U[i][i] = 1;
        for (int j = i + 1; j < n; ++j) {
            sum = 0;
            for (int k = 0; k < i; ++k) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = (get_matrix()[i][j] - sum) / L[i][i];
        }
    }
    std::vector<Matrix<T>> result;

    result.push_back(L);
    result.push_back(U);

    return result;
}


template <Field T>
constexpr double Matrix<T>::vector_max_diff(const std::vector<T> a, const std::vector<T> b) {
    int n = a.size();
    if (n != b.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    double max = -1;
    for (int i = 0; i < n; ++i) {
        double diff = std::abs(a[i] - b[i]);
        if (max < diff) {
            max = diff;
        }
    }
    return max;
}


template <Field T>
constexpr std::vector<T> Matrix<T>::jacobi_method(std::vector<T> b, const double eps) {
    int n = rows();
    std::vector<T> x(n, 1);
    std::vector<T> x_prev(n);

    int counter = 1;
    while (vector_max_diff(x, x_prev) >= eps) {
        std::cout << counter++  << "-th iteration"<< std::endl;
        for (int i = 0; i < x.size(); ++i) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;

        x_prev = x;
        for (int i = 0; i < n; ++i) {
            double sum_1 = 0;
            double sum_2 = 0;
            for (int j = 0; j < i; ++j) {
                sum_1 += get_matrix()[i][j] * x_prev[j];
            }
            sum_1 *= -1;
            for (int j = i + 1; j < n; ++j) {
                sum_2 += get_matrix()[i][j] * x_prev[j];
            }
            sum_2 *= -1;
            x[i] = (sum_1 + sum_2 + b[i]) / get_matrix()[i][i];
        }
    }
    return x;
}


template <Field T>
constexpr std::vector<T> Matrix<T>::lower_triang_equation(const std::vector<T> b) {
    const int n = b.size();
    std::vector<T> solution(n, 0);

    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += solution[j] * get_matrix()[i][j];
        }
        solution[i] = (b[i] - sum) / get_matrix()[i][i];
    }
    return solution;
}




template <Field T>
constexpr std::vector<T> Matrix<T>::upper_triang_equation(const std::vector<T> b) {
    const int n = b.size();
    std::vector<T> solution(n, 0);

    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = n; j > i; --j) {
            sum += solution[j] * get_matrix()[i][j];
        }
        solution[i] = (b[i] - sum) / get_matrix()[i][i];
    }
    return solution;
}

template <Field T>
constexpr std::vector<T> Matrix<T>::system_of_equations_lu(const std::vector<T> b) {
    std::vector<Matrix<T>> LU = crout_decomposition();
    std::vector<T> y = LU[0].lower_triang_equation(b);
    std::vector<T> x = LU[1].upper_triang_equation(y);

    return x;
}


template <Field T>
constexpr Matrix<T> Matrix<T>::inverse_lu() {
    if (determinant() == 0) {
        throw std::invalid_argument("The matrix is not invertible");
    }
    const int n = rows();
    Matrix<T> A_inv_t(n, n);

    std::vector<T> kronecker_vec(n, 0);
    kronecker_vec[0] = 1;
    for (int i = 0; i < n; ++i) {
        A_inv_t[i] = system_of_equations_lu(kronecker_vec);
        if (i + 1 < n) {
            kronecker_vec[i] = 0;
            kronecker_vec[i + 1] = 1;
        }
    }
    Matrix<T> A_inv = A_inv_t.transpose();

    return A_inv;
}



#endif // MATRIX_H
