#pragma once

#include "headers.hpp"

template <class U>
struct Edge {
    int v1;
    int v2;
    U value;
};


template <class U>
using adjacency_list_container = std::vector<std::list<U>>;

template <class T>
using matrix = std::vector<std::vector<T>>;



template <class T = char, class U = int>
class AbstractGraph
{
public:
    AbstractGraph();
    explicit AbstractGraph(int num);
    AbstractGraph(const AbstractGraph<T, U>& other);
    AbstractGraph(AbstractGraph<T, U>&& other) noexcept;
    virtual ~AbstractGraph();
public:
    virtual bool add_edge(const int&, const int&) = 0;
    virtual bool remove_edge(const int, const int) { return false; }
    virtual bool add_vertex(const T&);
    virtual bool add_vertex();
    // remove vertex
    bool print();
public:
    AbstractGraph<T, U>& operator=(const AbstractGraph<T, U>&) noexcept;
    AbstractGraph<T, U>& operator=(AbstractGraph<T, U>&&) noexcept;
public:
    [[nodiscard]] int get_number_of_vertices() const noexcept;
    [[nodiscard]] int get_number_of_edges() const noexcept;
    virtual bool set_number_of_vertices(int) noexcept;
    std::vector<T>& get_vertices() noexcept;
    std::vector<T> get_vertices() const noexcept;
    adjacency_list_container<U>& get_adjacency_list() noexcept;
    const adjacency_list_container<U> get_adjacency_list() const noexcept;
private:
    int m_number_of_vertices;
    adjacency_list_container<U> m_adjacency_list;
    std::vector<T> m_vertices;
};



template <class T, class U>
AbstractGraph<T, U>::AbstractGraph() : m_number_of_vertices(0) {}


template <class T, class U>
AbstractGraph<T, U>::AbstractGraph(int num) : m_number_of_vertices(num) {
    this->get_vertices().resize(num);
    this->get_adjacency_list().resize(num);
}


template <class T, class U>
AbstractGraph<T, U>::AbstractGraph(const AbstractGraph<T, U>& other)
    : m_number_of_vertices(other.m_number_of_vertices),
    m_adjacency_list(other.m_adjacency_list),
    m_vertices(other.m_vertices) {}


template <class T, class U>
AbstractGraph<T, U>::AbstractGraph(AbstractGraph<T, U>&& other) noexcept
    : m_number_of_vertices(other.m_number_of_vertices),
    m_adjacency_list(std::move(other.m_adjacency_list)),
    m_vertices(std::move(other.m_vertices)) {
        other.set_number_of_vertices(0);
}


template <class T, class U>
AbstractGraph<T, U>::~AbstractGraph() {
    m_adjacency_list.clear();
    m_vertices.clear();
}


template <class T, class U>
bool AbstractGraph<T, U>::add_vertex(const T& item) {
    set_number_of_vertices(get_number_of_vertices() + 1);
    get_vertices().emplace_back(item);
    return true;
}


template <class T, class U>
bool AbstractGraph<T, U>::add_vertex() {
    set_number_of_vertices(get_number_of_vertices() + 1);
    get_vertices().resize(get_number_of_vertices());
    return true;
}


template <class T, class U>
bool AbstractGraph<T, U>::print() {
    adjacency_list_container<U> tmp = get_adjacency_list();
    int i = 0;
    for (const auto& vec : tmp) {
        std::cout << i << " -> ";
        for(const auto& el : vec) {
            std::cout << el << " ";
        }
        ++i;
        std::cout << std::endl;
    }
    return true;
}


template <class T, class U>
AbstractGraph<T, U>& AbstractGraph<T, U>::operator=(const AbstractGraph<T, U>& other) noexcept {
    if (this != &other) {
        set_number_of_vertices(other.get_number_of_vertices());
        get_vertices().clear();
        get_vertices() = other.get_vertices();
        get_adjacency_list().clear();
        get_adjacency_list() = other.get_adjacency_list();
    }
    return *this;
}


template <class T, class U>
AbstractGraph<T, U>& AbstractGraph<T, U>::operator=(AbstractGraph<T, U>&& other) noexcept {
    set_number_of_vertices(other.get_number_of_vertices());
    get_vertices().clear();
    get_vertices() = std::move(other.get_vertices());
    get_adjacency_list().clear();
    get_adjacency_list() = std::move(other.get_adjacency_list());
    other.get_adjacency_list().clear();
    return *this;
}



template <class T, class U>
int AbstractGraph<T, U>::get_number_of_vertices() const noexcept {
    return m_number_of_vertices;
}


template <class T, class U>
int AbstractGraph<T, U>::get_number_of_edges() const noexcept {
    int counter = 0;
    const adjacency_list_container<U>& tmp = get_adjacency_list();
    for (const auto& vec : tmp) {
        counter += vec.size();
    }
    return counter;
}


template <class T, class U>
bool AbstractGraph<T, U>::set_number_of_vertices(const int num) noexcept {
    m_number_of_vertices = num;
    m_adjacency_list.resize(num);
    m_vertices.resize(num);
    return true;
}


template <class T, class U>
adjacency_list_container<U>& AbstractGraph<T, U>::get_adjacency_list() noexcept {
    return m_adjacency_list;
}


template <class T, class U>
const adjacency_list_container<U> AbstractGraph<T, U>::get_adjacency_list() const noexcept {
    return m_adjacency_list;
}


template <class T, class U>
std::vector<T> AbstractGraph<T, U>::get_vertices() const noexcept {
    return m_vertices;
}


template <class T, class U>
std::vector<T>& AbstractGraph<T, U>::get_vertices() noexcept {
    return m_vertices;
}


