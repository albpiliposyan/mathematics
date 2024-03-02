#ifndef HELPERS_H
#define HELPERS_H

#include <iostream>
#include <cmath>
#include <vector>


template <typename T>
T dot_product(const std::vector<T> a, const std::vector<T> b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors sizes don't match");
    }
    T result = 0;
    const int n = a.size();
    for (int i = 0; i < n; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

template <typename T>
T vector_norm(const std::vector<T> vec) {
    T sum = dot_product(vec, vec);
    return std::sqrt(sum);
}

template <typename T>
double vector_max_diff(const std::vector<T> a, const std::vector<T> b) {
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



#endif //HELPERS_H
