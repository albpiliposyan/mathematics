#ifndef FIELD_H
#define FIELD_H

#include "ring.h"
#include <stdexcept>

template <class T>
concept Divisible =std::regular<T> && requires(T x, T y) {
    { x + y } -> std::convertible_to<T>;
    { x - y } -> std::convertible_to<T>;
    { x * y } -> std::convertible_to<T>;
    { x / y } -> std::convertible_to<T>;
};

template <Divisible T = int>
class Field : public Ring<T>
{
public:
    Field(T x) : Ring<T>(x) {}
    Field(const Field<T>& other) : Ring<T>(other.get_value()) {}
    Field(Field<T>&& other) : Ring<T>(std::move(other.get_value()))
                { other.set_value(T()); }
    ~Field() = default;
public:
    virtual Field& operator=(const Field& rhs) {
        this->set_value(rhs.get_value());
        return *this;
    }
    virtual Field& operator+=(const Field& rhs) {
        this->set_value(this->get_value() + rhs.get_value());
        return *this;
    }
    virtual Field& operator-=(const Field& rhs) {
        this->set_value(this->get_value() - rhs.get_value());
        return *this;
    }
    virtual Field& operator*=(const Field& rhs) {
        this->set_value(this->get_value() * rhs.get_value());
        return *this;
    }
    virtual Field& operator/=(const Field& rhs) {
        if (rhs.get_value() == rhs.identity()) {
            throw std::runtime_error("Cannot divide by zero.");
        }
        this->set_value(this->get_value() / rhs.get_value());
        return *this;
    }
    virtual bool operator==(const Field& rhs) {
        return this->get_value() == rhs.get_value();
    }
    friend Field operator+(Field lhs, const Field& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend Field operator-(Field lhs, const Field& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend Field operator-(Field lhs) {
        lhs.set_value(-1 * lhs.get_value());
        return lhs;
    }
    friend Field operator*(Field lhs, const Field& rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend Field operator/(Field lhs, const Field& rhs) {
        lhs /= rhs;
        return lhs;
    }
    friend std::istream& operator>>(std::istream& in, Field& g) {
        T x;
        in >> x;
        g.set_value(x);
        return in;
    }
    friend std::ostream& operator<<(std::ostream& out, const Field& g) {
        out << g.get_value();
        return out;
    }
};

#endif
