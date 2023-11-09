#pragma once

#include "abstract_graph.hpp"

template <class T>
class Graph : public AbstractGraph<T, unsigned int>
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
    constexpr bool remove_edge(const unsigned int ver1, const unsigned int ver2) override;
    constexpr bool is_tree() const noexcept;
    std::vector<unsigned int> cayley_representation() const;
};


template <class T>
std::vector<unsigned int> Graph<T>::cayley_representation() const {
    std::vector<unsigned int> result;
    if (this->get_number_of_vertices() <= 2) {
        return result; 
    }
    if (!this->is_tree()) {
        // throw an error
        return result;
    }
    Graph tmp(*this);
    unsigned int n = this->get_number_of_vertices();
    result.reserve(n - 2);
    std::vector<bool> check_used(n, false);
    for(unsigned int i = 0; i < n - 2; ++i) {
        unsigned int j = 0;
        for (j = 0; j < n; ++j) {
            if (tmp.get_adjacency_list()[j].size() == 1) {
                break;
            }
        }
        result.push_back(j);
        unsigned int index = tmp.get_adjacency_list()[j].front();
        tmp.remove_edge(j, index);
    }
    return result;
}



template <class T>
Graph<T>::Graph() : AbstractGraph<T, unsigned int>() {}

template <class T>
Graph<T>::Graph(int num) : AbstractGraph<T, unsigned int>(num) {}

template <class T>
Graph<T>::Graph(const Graph<T>& other) : AbstractGraph<T, unsigned int>(other) {}

template <class T>
Graph<T>::Graph(const Graph<T>&& other) : AbstractGraph<T, unsigned int>(other) {}

template <class T>
constexpr Graph<T>& Graph<T>::operator=(const Graph<T>& other) noexcept {
    if (this != &other) {
        AbstractGraph<T, unsigned int>::operator=(other);
    }
    return *this;
}

template <class T>
constexpr Graph<T>& Graph<T>::operator=(const Graph<T>&& other) noexcept {
    if (this != &other) {
        AbstractGraph<T, unsigned int>::operator=(std::move(other));
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

template <class T>
constexpr bool Graph<T>::remove_edge(const unsigned int ver1, const unsigned int ver2) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    this->get_adjacency_list()[ver1].remove(ver2);
    this->get_adjacency_list()[ver2].remove(ver1);
    return true;
}

template <class T>
constexpr bool Graph<T>::is_tree() const noexcept {
    // G(m, n)
    // check n = m - 1
    // there is no 
    // TODO
    return true;
    for (auto& el : this->get_adjacency_list()) {
        if (el.size() > 2) {
            return false;
        }
    }
    return true;
}
