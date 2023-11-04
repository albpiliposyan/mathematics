#ifndef GROUP_H
#define GROUP_H

template <class T>
concept Addible =std::regular<T> && requires(T x, T y) {
    { x +  y } -> std::convertible_to<T>;
    { x - y } -> std::convertible_to<T>;
};

template <Addible T = int>
class Group
{
public:
    Group(T x) : m_value(x) {}
    Group(const Group<T>& other) : m_value(other.get_value()) {}
    Group(Group<T>&& other) : m_value(std::move(other.get_value()))
                { other.set_value(T()); }
    virtual ~Group() = default;
public:
    virtual Group& operator+=(const Group& rhs) {
        m_value += rhs.get_value();
        return *this;
    }
    virtual Group& operator-=(const Group& rhs) {
        m_value -= rhs.get_value();
        return *this;
    }
    virtual Group& operator=(const Group& rhs) {
        m_value = rhs.get_value();
        return *this;
    }
    virtual bool operator==(const Group& rhs) {
        return m_value == rhs.get_value();
    }

public:
    virtual Group& initialize_group() { return *this; }
    virtual T identity() const  { return 0; };
    virtual T inverse(const T& element) const { return 0; };
public:
    friend Group operator+(Group lhs, const Group& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend Group operator-(Group lhs, const Group& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend Group operator-(Group lhs) {
        lhs.set_value(-1 * lhs.get_value());
        return lhs;
    }
    friend std::istream& operator>>(std::istream& in, Group& g) {
        T x;
        in >> x;
        g.set_value(x);
        return in;
    }
    friend std::ostream& operator<<(std::ostream& out, const Group& g) {
        out << g.get_value();
        return out;
    }
public:
    virtual T get_value() const {
        return m_value;
    }
    virtual bool set_value(T x) {
        m_value = x;
        return true;
    }
private:
    T m_value;
};

#endif
