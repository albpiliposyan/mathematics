#ifndef RING_H
#define RING_H

#include "group.h"
template <class T>
concept Multiplicable =std::regular<T> && requires(T x, T y) {
    { x + y } -> std::convertible_to<T>;
    { x - y } -> std::convertible_to<T>;
    { x * y } -> std::convertible_to<T>;
};

template <Multiplicable T = int>
class Ring : public Group<T>
{
public:
    Ring(T x) : Group<T>(x) {}
    Ring(const Ring<T>& other) : Group<T>(other.get_value()) {}
    Ring(Ring<T>&& other) : Group<T>(std::move(other.get_value()))
                { other.set_value(T()); }
    ~Ring() = default;
public:
    virtual Ring& operator=(const Ring& rhs) {
        this->set_value(rhs.get_value());
        return *this;
    }
    virtual Ring& operator+=(const Ring& rhs) {
        this->set_value(this->get_value() + rhs.get_value());
        return *this;
    }
    virtual Ring& operator-=(const Ring& rhs) {
        this->set_value(this->get_value() - rhs.get_value());
        return *this;
    }
    virtual Ring& operator*=(const Ring& rhs) {
        this->set_value(this->get_value() * rhs.get_value());
        return *this;
    }
    virtual bool operator==(const Ring& rhs) {
        return this->get_value() == rhs.get_value();
    }
    friend Ring operator+(Ring lhs, const Ring& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend Ring operator-(Ring lhs, const Ring& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend Ring operator-(Ring lhs) {
        lhs.set_value(-1 * lhs.get_value());
        return lhs;
    }
    friend Ring operator*(Ring lhs, const Ring& rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend std::istream& operator>>(std::istream& in, Ring& g) {
        T x;
        in >> x;
        g.set_value(x);
        return in;
    }
    friend std::ostream& operator<<(std::ostream& out, const Ring& g) {
        out << g.get_value();
        return out;
    }
};

#endif
