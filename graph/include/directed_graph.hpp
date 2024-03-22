#pragma once

#include "abstract_graph.hpp"

template <class T = char>
class DirectedGraph : public AbstractGraph<T, int>
{
public:
    DirectedGraph();
    explicit DirectedGraph(int num);
    DirectedGraph(const DirectedGraph<T>& other);
    DirectedGraph(DirectedGraph<T>&& other) noexcept;
    ~DirectedGraph() = default;
public:
    DirectedGraph<T>& operator=(const DirectedGraph<T>& other) noexcept;
    DirectedGraph<T>& operator=(DirectedGraph<T>&& other) noexcept;
public:
    bool add_edge(const int& ver1, const int& ver2) override;
    bool remove_edge(int ver1, int ver2) override;
};



template <class T>
DirectedGraph<T>::DirectedGraph() : AbstractGraph<T, int>() {}


template <class T>
DirectedGraph<T>::DirectedGraph(int num) : AbstractGraph<T, int>(num) {}


template <class T>
DirectedGraph<T>::DirectedGraph(const DirectedGraph<T>& other)
    : AbstractGraph<T, int>(other) {}


template <class T>
DirectedGraph<T>::DirectedGraph(DirectedGraph<T>&& other) noexcept
    : AbstractGraph<T, int>(other) {}


template <class T>
DirectedGraph<T>& DirectedGraph<T>::operator=(const DirectedGraph<T>& other) noexcept {
    if (this != &other) {
        AbstractGraph<T, int>::operator=(other);
    }
    return *this;
}


template <class T>
DirectedGraph<T>& DirectedGraph<T>::operator=(DirectedGraph<T>&& other) noexcept {
    if (this != &other) {
        AbstractGraph<T, int>::operator=(std::move(other));
    }
    return *this;
}


template <class T>
bool DirectedGraph<T>::add_edge(const int& ver1, const int& ver2) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    auto& list_to_insert_1 = this->get_adjacency_list()[ver1];
    auto it = std::lower_bound(list_to_insert_1.begin(), list_to_insert_1.end(), ver2);
    list_to_insert_1.insert(it, ver2);
    return true;
}


template <class T>
bool DirectedGraph<T>::remove_edge(int ver1, int ver2) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    this->get_adjacency_list()[ver1].remove(ver2);
    return true;
}


