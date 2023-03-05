#pragma once

#include "abstract_graph.h"

template <class T>
class Graph : public AbstractGraph<T>
{
public:
    Graph();
    Graph(int num);
    Graph(const Graph<T>& other);
    Graph(const Graph<T>&& other);
    ~Graph() = default;
public:
    constexpr Graph<T>& operator=(const Graph<T>& other) noexcept;
    constexpr Graph<T>& operator=(const Graph<T>&& other) noexcept;
public:
    constexpr bool add_edge(const unsigned int& ver1, const unsigned int& ver2) override;
};



template <class T>
Graph<T>::Graph() : AbstractGraph<T>() {}

template <class T>
Graph<T>::Graph(int num) : AbstractGraph<T>(num) {}

template <class T>
Graph<T>::Graph(const Graph<T>& other) : AbstractGraph<T>(other) {}

template <class T>
Graph<T>::Graph(const Graph<T>&& other) : AbstractGraph<T>(other) {}

template <class T>
constexpr Graph<T>& Graph<T>::operator=(const Graph<T>& other) noexcept {
    if (this != &other) {
        AbstractGraph<T>::operator=(other);
    }
    return *this;
}

template <class T>
constexpr Graph<T>& Graph<T>::operator=(const Graph<T>&& other) noexcept {
    if (this != &other) {
        AbstractGraph<T>::operator=(std::move(other));
    }
    return *this;
}

template <class T>
constexpr bool Graph<T>::add_edge(const unsigned int& ver1, const unsigned int& ver2) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    auto& list_to_insert_1 = this->get_adjacency_list()[ver1];
    auto it = std::lower_bound(list_to_insert_1.begin(), list_to_insert_1.end(), ver2);
    if (it != list_to_insert_1.end() && *it == ver2) {
        throw std::invalid_argument("Edge already exists");
    }
    list_to_insert_1.insert(it, ver2);


    auto& list_to_insert_2 = this->get_adjacency_list()[ver2];
    it = std::lower_bound(list_to_insert_2.begin(), list_to_insert_2.end(), ver1);
    if (it != list_to_insert_2.end() && *it == ver1) {
        throw std::invalid_argument("Edge already exists");
    }
    list_to_insert_2.insert(it, ver1);
    return true;
}
