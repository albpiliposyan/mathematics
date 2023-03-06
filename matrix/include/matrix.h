#pragma once

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
    constexpr bool initialize();
    constexpr bool print() noexcept;
    constexpr T determinant() const;
    constexpr Matrix<T>& transpose() noexcept;
    constexpr Matrix<T>& swapRows(const unsigned int& row1, const unsigned int& row2);
    constexpr Matrix<T>& swapColumns(const unsigned int& col1, const unsigned int& col2);
    constexpr Matrix<T>& multiply(const T& scalar);
    constexpr Matrix<T>& multiplyRow(const unsigned int& row, const T& scalar);
    constexpr Matrix<T>& addMultipleOfRow(const unsigned int& row1,
                        const unsigned int& row2, const T& scalar);
    constexpr bool checkZeroRow(const unsigned int row) const noexcept;
    constexpr bool checkZeroColumn(const unsigned int col) const noexcept;
    constexpr unsigned int numberOfZeroRows() const noexcept;
    constexpr unsigned int numberOfZeroColumns() const noexcept;
    constexpr bool swapNonZeroRowsToTop() noexcept;
    constexpr bool swapNonZeroColumnsToLeft() noexcept;
    constexpr Matrix<T>& reduceToRowEchelonForm();
    constexpr std::vector<T> gaussEliminationMethod();
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
    constexpr const TwoDVector<T>& get_matrix() const noexcept;
    constexpr bool set_matrix(TwoDVector<T> other) noexcept;
    constexpr unsigned int rows() const noexcept;
    constexpr unsigned int cols() const noexcept;
    constexpr bool set_rows(unsigned int row) noexcept;
    constexpr bool set_cols(unsigned int col) noexcept;
private:
    constexpr T cofactor(const unsigned int& row, const unsigned int& col) const;
private:
    unsigned int m_rows;
    unsigned int m_cols;
    TwoDVector<T> m_matrix;
};



template <Field T>
constexpr bool Matrix<T>::initialize() {
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
    return true;
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
constexpr Matrix<T>& Matrix<T>::transpose() noexcept {
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
constexpr Matrix<T>& Matrix<T>::swapRows(const unsigned int& row1, const unsigned int& row2) {
    if (row1 < 0 || row1 >= rows() || row2 < 0 || row2 >= rows())
        throw std::out_of_range("Invalid row index");
    std::swap(get_matrix()[row1], get_matrix()[row2]);
    return *this;
}

template <Field T>
constexpr Matrix<T>& Matrix<T>::swapColumns(const unsigned int& col1, const unsigned int& col2) {
    if (col1 < 0 || col1 >= cols() || col2 < 0 || col2 >= cols())
        throw std::out_of_range("Invalid column index");
    for (unsigned int i = 0; i < rows(); ++i)
        std::swap(m_matrix[i][col1], m_matrix[i][col2]);
    return *this;
}

template <Field T>
constexpr Matrix<T>& Matrix<T>::multiply(const T& scalar) {
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
constexpr Matrix<T>& Matrix<T>::multiplyRow(const unsigned int& row, const T& scalar) {
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
constexpr Matrix<T>& Matrix<T>::addMultipleOfRow(const unsigned int& row1, 
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
constexpr bool Matrix<T>::checkZeroRow(const unsigned int row) const noexcept {
    for (unsigned int i = 0; i < cols(); ++i) {
        if (get_matrix()[row][i] != 0) {
            return false;
        }
    }
    return true;
}

template <Field T>
constexpr bool Matrix<T>::checkZeroColumn(const unsigned int col) const noexcept{
    for (unsigned int i = 0; i < rows(); ++i) {
        if (get_matrix()[i][col] != 0) {
            return false;
        }
    }
    return true;
}

template <Field T>
constexpr bool Matrix<T>::swapNonZeroRowsToTop() noexcept {
    unsigned int row_index = rows() - 1;
    for (unsigned int i = 0; i <= row_index; ++i) {
        if (checkZeroRow(i)) {
            swapRows(i--, row_index--);
        }
    }
    return true;
}

template <Field T>
constexpr bool Matrix<T>::swapNonZeroColumnsToLeft() noexcept {
    unsigned int col_index = 0;
    for (unsigned int i = cols() - 1; i > col_index; --i) {
        if (checkZeroColumn(i)) {
            swapColumns(i++, col_index++);
        }
    }
    return true;
}

template <Field T>
constexpr unsigned int Matrix<T>::numberOfZeroRows() const noexcept {
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
constexpr unsigned int Matrix<T>::numberOfZeroColumns() const noexcept {
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
constexpr Matrix<T>& Matrix<T>::reduceToRowEchelonForm() {
    unsigned int zero_rows_count = numberOfZeroRows();
    if (zero_rows_count == rows()) {
        return *this;
    }
    unsigned int zero_cols_count = numberOfZeroColumns();
    if (zero_rows_count > 0) {
        swapNonZeroRowsToTop();
    }
    if (zero_cols_count > 0) {
        swapNonZeroColumnsToLeft();
    }

    for (unsigned int i = 0; i < rows() - zero_rows_count; ++i) {
        for (unsigned int j = i + 1; j < rows() - zero_rows_count; ++j) {
            if (std::abs(get_matrix()[i][zero_cols_count + i]) < 
                                    std::abs(get_matrix()[j][zero_cols_count + i])) {
                swapRows(i, j);
            }
        }
    }
    for (unsigned int i = 0; i < rows() - zero_rows_count - 1; ++i) {
        for (unsigned int j = i + 1; j < rows() - zero_rows_count; ++j) {
            if (get_matrix()[i][zero_cols_count + i] != 0) {
                addMultipleOfRow(j, i, (-1) * get_matrix()[j][zero_cols_count + i] / 
                                        get_matrix()[i][zero_cols_count + i]);
            }
        }
    }
    swapNonZeroRowsToTop();
    swapNonZeroColumnsToLeft();
    // print();
    return *this;
}

template <Field T>
constexpr std::vector<T> Matrix<T>::gaussEliminationMethod() {
    Matrix tmp(*this);
    tmp.reduceToRowEchelonForm();
    unsigned int zero_rows_count = tmp.numberOfZeroRows();
    unsigned int zero_cols_count = tmp.numberOfZeroColumns();    
    
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
    TwoDVector<T> result(rows(), std::vector<T>(rhs.cols(), 0));
    for (unsigned int i = 0; i < result.size(); ++i)
        for (unsigned int j = 0; j < result[0].size(); ++j)
            for (unsigned int k = 0; k < cols(); ++k)
                result[i][j] += get_matrix()[i][k] * rhs[k][j];
    set_matrix(result);
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
constexpr const TwoDVector<T>& Matrix<T>::get_matrix() const noexcept {
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
    this->get_matrix().resize(row);
    return true;
}

template <Field T>
constexpr bool Matrix<T>::set_cols(unsigned int col) noexcept {
    for (auto& row : this->get_matrix()) {
        row.resize(col);
    }
    return true;
}
