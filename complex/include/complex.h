#ifndef COMPLEX_H
#define COMPLEX_H

#include <stdexcept>
#include <utility>
#include <cmath>

template <typename T>
concept ComplexNumber = std::regular<T> && requires(T a, T b) {
    { a + b } -> std::same_as<T>;
    { a - b } -> std::same_as<T>;
};



template <ComplexNumber T>
class Complex;


template <ComplexNumber U>
Complex<U> operator+(Complex<U> lhs, const Complex<U>& rhs) {
    lhs += rhs;
    return lhs;
}
template <ComplexNumber U>
Complex<U> operator-(Complex<U> lhs, const Complex<U>& rhs) {
    lhs -= rhs;
    return lhs;
}
template <ComplexNumber U>
Complex<U> operator*(Complex<U> lhs, const Complex<U>& rhs) {
    lhs *= rhs;
    return lhs;
}
template <ComplexNumber U>
Complex<U> operator/(Complex<U> lhs, const Complex<U>& rhs) {
    lhs /= rhs;
    return lhs;
}
template <ComplexNumber U>
std::istream& operator>>(std::istream& in, Complex<U>& c) {
    U real, imag;
    in >> real >> imag;
    c.set_real(real);
    c.set_imag(imag);
    return in;
}
template <ComplexNumber U>
std::ostream& operator<<(std::ostream& out, const Complex<U>& c) {
    if (c.real() != 0)
        std::cout << c.real();
    if (c.real() != 0 && c.imag() > 0)
        std::cout << "+";
    else if (c.imag() == 0)
        return out;
    std::cout << c.imag() << "i";
    return out;
}


template <ComplexNumber T>
class Complex
{
public:
    Complex() = default;
    Complex(T x = T(),T y = T()) : m_real(x), m_imag(y) {}
    Complex(const Complex& other) : m_real(other.m_real), m_imag(other.m_imag) {}
    Complex(Complex&& other) noexcept :
            m_real(std::move(other.m_real)), m_imag(std::move(other.m_imag)) {
        other.set_real(0);
        other.set_imag(0);
    }
    ~Complex() = default;
public:
    constexpr Complex<T>& operator+=(const Complex<T>& rhs);
    constexpr Complex<T>& operator-=(const Complex<T>& rhs);
    constexpr Complex<T>& operator*=(const Complex<T>& rhs);
    constexpr Complex<T>& operator/=(const Complex<T>& rhs);
    Complex<T>& operator=(const Complex<T>& other);
    bool operator==(const Complex<T>& rhs);
    template <ComplexNumber U>
    friend Complex<T> operator+(Complex<T> lhs, const Complex<T>& rhs);
    template <ComplexNumber U>
    friend Complex<T> operator-(Complex<T> lhs, const Complex<T>& rhs);
    template <ComplexNumber U>
    friend Complex<T> operator*(Complex<T> lhs, const Complex<T>& rhs);
    template <ComplexNumber U>
    friend Complex<T> operator/(Complex<T> lhs, const Complex<T>& rhs);
    template <ComplexNumber U>
    friend std::istream& operator>>(std::istream& in, Complex<T>& c);
    template <ComplexNumber U>
    friend std::ostream& operator<<(std::ostream& out, const Complex<T>& c);
public:
    Complex<T> conjugate() const;
    double magnitude() const;
public:
    T real() const;
    T imag() const;
    bool set_real(T x);
    bool set_imag(T y);
private:
    T m_real;
    T m_imag;
};

#endif

template <ComplexNumber T>
constexpr Complex<T>& Complex<T>::operator+=(const Complex<T>& rhs) {
    this->set_real(this->real() + rhs.real());
    this->set_imag(this->imag() + rhs.imag());
    return *this;
}

template <ComplexNumber T>
constexpr Complex<T>& Complex<T>::operator-=(const Complex<T>& rhs) {
    this->set_real(this->real() - rhs.real());
    this->set_imag(this->imag() - rhs.imag());
    return *this;
}

template <ComplexNumber T>
constexpr Complex<T>& Complex<T>::operator*=(const Complex<T>& rhs) {
    T this_real = this->real(), rhs_real = rhs.real();
    this->set_real(this->real() * rhs.real() - this->imag() * rhs.imag());
    this->set_imag(this_real * rhs.imag() + this->imag() * rhs_real);
    return *this;
}

template <ComplexNumber T>
constexpr Complex<T>& Complex<T>::operator/=(const Complex<T>& rhs) {
    if (rhs.real() == 0 && rhs.imag() == 0)
        throw std::logic_error("Division by zero.");
    T this_real = this->real(), rhs_real = rhs.real();
    this->set_real((this->real() * rhs.real() + rhs.imag() * this->imag())
            / (std::pow(rhs.real(), 2) + std::pow(rhs.imag(), 2)));
    this->set_imag((this->imag() * rhs_real - this_real * rhs.imag())
            / (std::pow(rhs_real, 2) + std::pow(rhs.imag(), 2)));
    return *this;
}

template <ComplexNumber T>
Complex<T>& Complex<T>::operator=(const Complex<T>& other) {
    this->set_real(other.real());
    this->set_imag(other.imag());
    return *this;
}

template <ComplexNumber T>
bool Complex<T>::operator==(const Complex<T>& rhs) {
    return this->real() == rhs.real() && this->imag() == rhs.imag();
}

template <ComplexNumber T>
Complex<T> Complex<T>::conjugate() const {
    return Complex<T>(this->real(), (-1) * this->imag());
}

template <ComplexNumber T>
double Complex<T>::magnitude() const {
    return std::sqrt(std::pow(this->real(), 2) + std::pow(this->imag(), 2));
}

template <ComplexNumber T>
T Complex<T>::real() const {
    return m_real;
}

template <ComplexNumber T>
T Complex<T>::imag() const {
    return m_imag;
}

template <ComplexNumber T>
bool Complex<T>::set_real(T x) {
    m_real = x;
    return true;
}

template <ComplexNumber T>
bool Complex<T>::set_imag(T y) {
    m_imag = y;
    return true;
}
