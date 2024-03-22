#pragma once

#include "abstract_graph.hpp"

template <class T>
class Graph : public AbstractGraph<T, int>
{
public:
    Graph();
    explicit Graph(int num);
    Graph(const Graph<T>& other);
    Graph(Graph<T>&& other) noexcept;
    ~Graph() = default;
public:
    Graph<T>& operator=(const Graph<T>& other) noexcept;
    Graph<T>& operator=(Graph<T>&& other) noexcept;
public:
    bool add_edge(const int& ver1, const int& ver2) override;
    bool remove_edge(const int&, const int&) override;
    std::vector<int> cayley_representation();
};


template <class T>
Graph<T>::Graph() : AbstractGraph<T, int>() {}


template <class T>
Graph<T>::Graph(int num) : AbstractGraph<T, int>(num) {}


template <class T>
Graph<T>::Graph(const Graph<T>& other) : AbstractGraph<T, int>(other) {}


template <class T>
Graph<T>::Graph(Graph<T>&& other) noexcept :
                    AbstractGraph<T, int>(other) {}


template <class T>
Graph<T>& Graph<T>::operator=(const Graph<T>& other) noexcept {
    if (this != &other) {
        AbstractGraph<T, int>::operator=(other);
    }
    return *this;
}


template <class T>
Graph<T>& Graph<T>::operator=(Graph<T>&& other) noexcept {
    if (this != &other) {
        AbstractGraph<T, int>::operator=(std::move(other));
    }
    return *this;
}


template <class T>
bool Graph<T>::add_edge(const int& ver1, const int& ver2) {
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


template <class T>
bool Graph<T>::remove_edge(const int ver1, const int ver2) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    this->get_adjacency_list()[ver1].remove(ver2);
    this->get_adjacency_list()[ver2].remove(ver1);
    return true;
}


