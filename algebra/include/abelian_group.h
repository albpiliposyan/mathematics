#include "group.h"

template <class T>
concept Abelian =std::regular<T> && requires(T x, T y) {
    { x + y } -> std::convertible_to<T>;
    { x - y } -> std::convertible_to<T>;
    { y + x } -> std::convertible_to<T>;
    { y - x } -> std::convertible_to<T>;
};

template <Abelian T = int>
class AbelianGroup : public Group<T>
{
public:
    AbelianGroup(T x) : Group<T>(x) {}
    AbelianGroup(const AbelianGroup<T>& other) : Group<T>(other) {}
    AbelianGroup(AbelianGroup<T>&& other) : Group<T>(std::move(other.get_value()))
                { other.set_value(T()); }
    ~AbelianGroup() = default;
public:
    virtual AbelianGroup& operator+=(const AbelianGroup& rhs) {
        this->set_value(this->get_value() + rhs.get_value());
        return *this;
    }
    virtual AbelianGroup& operator-=(const AbelianGroup& rhs) {
        this->set_value(this->get_value() - rhs.get_value());
        return *this;
    }
    virtual AbelianGroup& operator=(const AbelianGroup& rhs) {
        this->set_value(rhs.get_value());
        return *this;
    }
    virtual bool operator==(const AbelianGroup& rhs) {
        return this->get_value() == rhs.get_value();
    }
    friend AbelianGroup operator+(AbelianGroup lhs, const AbelianGroup& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend AbelianGroup operator-(AbelianGroup lhs, const AbelianGroup& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend AbelianGroup operator-(AbelianGroup lhs) {
        lhs.set_value(-1 * lhs.get_value());
        return lhs;
    }
    friend std::istream& operator>>(std::istream& in, AbelianGroup& g) {
        T x;
        in >> x;
        g.set_value(x);
        return in;
    }
    friend std::ostream& operator<<(std::ostream& out, const AbelianGroup& g) {
        out << g.get_value();
        return out;
    }
};
