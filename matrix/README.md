# Matrix Class

This is a C++ implementation of a Matrix class, providing a range of functions for performing operations on matrices.


## Requirements

The Matrix class requires a C++ compiler that supports C++20 or later.

## Usage

To use the Matrix class, simply include the Matrix.hpp header file in your source code and create a Matrix object with a given number of rows and columns.

```cpp
#include "Matrix.h"

// create a 3x3 matrix
Matrix<double> mat(3, 3); // creates a 3x3 matrix of type double

// fill the matrix with values
mat[0] = {1, 2, 3};
mat[1] = {4, 5, 6};
mat[2] = {7, 8, 9};

// find the determinant
std::cout << "Determiant: " << mat.determinant() << std::endl;

// reduce to row echelon form
mat.reduceToRowEchelonForm();
mat.print();
```

You can then use the various member functions to perform matrix operations, such as multiplication, addition, finding the determinant, reduction to row echelon form and gaussian elimination.



## Functionality

The Matrix class provides the following functionality:

 - Initialization of a matrix object with a given number of rows and columns
 - Setting and getting the value of individual elements in the matrix
 - Printing the matrix to standard output
 - Calculating the determinant of the matrix
 - Transposing the matrix
 - Swapping rows or columns of the matrix
 - Multiplying the matrix by a scalar or a row by a scalar
 - Adding a multiple of one row to another row
 - Checking if a row or column is all zeros
 - Swapping all non-zero rows to the top and all non-zero columns to the left
 - Reducing the matrix to row echelon form
 - Solving a system of linear equations using Gauss elimination method
 - Overloading of the assignment, addition, subtraction, and multiplication operators
 - Getting and setting the matrix object as a two-dimensional vector


## Work in Progress

**Note:** This class is currently a work in progress. Some functions may not be fully implemented, and others may undergo changes as development continues.
