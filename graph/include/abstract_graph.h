#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <stdexcept>
#include <algorithm>

using adjacency_list_container = std::vector<std::list<unsigned int>>;


template <class T>
class AbstractGraph
{
public:
    AbstractGraph();
    AbstractGraph(int num);
    AbstractGraph(const AbstractGraph<T>& other);
    AbstractGraph(const AbstractGraph<T>&& other);
    virtual ~AbstractGraph();
public:
    constexpr virtual bool add_edge(const unsigned int& ver1, const unsigned int& ver2) = 0;
    // remove edge
    constexpr virtual bool add_vertice(const T& item);
    constexpr virtual bool add_vertice();
    // remove vertice
    constexpr virtual bool print() noexcept;
public:
    constexpr AbstractGraph<T>& operator=(const AbstractGraph<T>& other) noexcept;
    constexpr AbstractGraph<T>& operator=(const AbstractGraph<T>&& other) noexcept;
public:
    constexpr unsigned int get_number_of_vertices() const noexcept;
    constexpr bool set_number_of_vertices(unsigned int num) noexcept;
    constexpr std::vector<T>& get_vertices() noexcept;
    constexpr std::vector<T> get_vertices() const noexcept;
    constexpr adjacency_list_container& get_adjacency_list() noexcept;
    constexpr adjacency_list_container get_adjacency_list() const noexcept;
private:
    unsigned int m_number_of_vertices;
    std::vector<T> m_vertices;
    adjacency_list_container m_adjacency_list;
};



template <class T>
AbstractGraph<T>::AbstractGraph() : m_number_of_vertices(0) {}

template <class T>
AbstractGraph<T>::AbstractGraph(int num) : m_number_of_vertices(num) {
    this->get_vertices().resize(num);
    this->get_adjacency_list().resize(num);
}

template <class T>
AbstractGraph<T>::AbstractGraph(const AbstractGraph<T>& other)
    : m_number_of_vertices(other.m_number_of_vertices),
            m_adjacency_list(other.m_adjacency_list),
            m_vertices(other.m_vertices) {}

template <class T>
AbstractGraph<T>::AbstractGraph(const AbstractGraph<T>&& other) 
    : m_number_of_vertices(other.m_number_of_vertices),
            m_adjacency_list(std::move(other.m_adjacency_list)),
            m_vertices(std::move(other.m_vertices)) {
    other.set_number_of_vertices(0);
}

template <class T>
AbstractGraph<T>::~AbstractGraph() {
    m_adjacency_list.clear();
    m_vertices.clear();
}

template <class T>
constexpr bool AbstractGraph<T>::add_vertice(const T& item) {
    set_number_of_vertices(get_number_of_vertices() + 1);
    get_vertices().emplace_back(item);
    return true;
}

template <class T>
constexpr bool AbstractGraph<T>::add_vertice() {
    set_number_of_vertices(get_number_of_vertices() + 1);
    get_vertices().resize(get_number_of_vertices());
    return true;
}


template <class T>
constexpr bool AbstractGraph<T>::print() noexcept {
    adjacency_list_container tmp = get_adjacency_list();
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

template <class T>
constexpr AbstractGraph<T>& AbstractGraph<T>::operator=(const AbstractGraph<T>& other) noexcept {
    if (this != &other) {
        set_number_of_vertices(other.get_number_of_vertices());
        get_vertices().clear();
        get_vertices() = other.get_vertices();
        get_adjacency_list().clear();
        get_adjacency_list() = other.get_adjacency_list();
    }
    return *this;
}

template <class T>
constexpr AbstractGraph<T>& AbstractGraph<T>::operator=(const AbstractGraph<T>&& other) noexcept {
    set_number_of_vertices(other.get_number_of_vertices());
    get_vertices().clear();
    get_vertices() = std::move(other.get_vertices());
    get_adjacency_list().clear();
    get_adjacency_list() = std::move(other.get_adjacency_list());
    other.get_adjacency_list().clear();
    return *this;
}



template <class T>
constexpr unsigned int AbstractGraph<T>::get_number_of_vertices() const noexcept {
    return m_number_of_vertices;
}

template <class T>
constexpr bool AbstractGraph<T>::set_number_of_vertices(unsigned int num) noexcept {
    m_number_of_vertices = num;
    m_adjacency_list.resize(num);
    m_vertices.resize(num);
    return true;
}

template <class T>
constexpr adjacency_list_container AbstractGraph<T>::get_adjacency_list() const noexcept {
    return m_adjacency_list;
}

template <class T>
constexpr adjacency_list_container& AbstractGraph<T>::get_adjacency_list() noexcept {
    return m_adjacency_list;
}

template <class T>
constexpr std::vector<T> AbstractGraph<T>::get_vertices() const noexcept {
    return m_vertices;
}

template <class T>
constexpr std::vector<T>& AbstractGraph<T>::get_vertices() noexcept {
    return m_vertices;
}





