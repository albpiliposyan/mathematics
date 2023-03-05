#pragma once

#include "abstract_graph.h"

template <class T = int>
class DirectedGraph : public AbstractGraph<T>
{
public:
    DirectedGraph();
    DirectedGraph(int num);
    DirectedGraph(const DirectedGraph<T>& other);
    DirectedGraph(const DirectedGraph<T>&& other);
    ~DirectedGraph() = default;
public:
    constexpr DirectedGraph<T>& operator=(const DirectedGraph<T>& other) noexcept;
    constexpr DirectedGraph<T>& operator=(const DirectedGraph<T>&& other) noexcept;
public:
    constexpr bool add_edge(const unsigned int& ver1, const unsigned int& ver2) override;
};



template <class T>
DirectedGraph<T>::DirectedGraph() : AbstractGraph<T>() {}

template <class T>
DirectedGraph<T>::DirectedGraph(int num) : AbstractGraph<T>(num) {}

template <class T>
DirectedGraph<T>::DirectedGraph(const DirectedGraph<T>& other)
    : AbstractGraph<T>(other) {}

template <class T>
DirectedGraph<T>::DirectedGraph(const DirectedGraph<T>&& other)
    : AbstractGraph<T>(other) {}

template <class T>
constexpr DirectedGraph<T>& DirectedGraph<T>::operator=(const DirectedGraph<T>& other) noexcept {
    if (this != &other) {
        AbstractGraph<T>::operator=(other);
    }
    return *this;
}

template <class T>
constexpr DirectedGraph<T>& DirectedGraph<T>::operator=(const DirectedGraph<T>&& other) noexcept {
    if (this != &other) {
        AbstractGraph<T>::operator=(std::move(other));
    }
    return *this;
}

template <class T>
constexpr bool DirectedGraph<T>::add_edge(const unsigned int& ver1, const unsigned int& ver2) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    auto& list_to_insert_1 = this->get_adjacency_list()[ver1];
    auto it = std::lower_bound(list_to_insert_1.begin(), list_to_insert_1.end(), ver2);
    list_to_insert_1.insert(it, ver2);
    return true;
}
