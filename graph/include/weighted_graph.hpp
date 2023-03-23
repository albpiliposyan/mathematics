#pragma once

#include "abstract_graph.hpp"


template <class T, class U>
class WeightedGraph : public AbstractGraph<T, std::pair<unsigned int, U> >
{
public:
    WeightedGraph();
    WeightedGraph(int num);
    WeightedGraph(const WeightedGraph<T, U>& other);
    WeightedGraph(const WeightedGraph<T, U>&& other);
    ~WeightedGraph() = default;
public:
    constexpr WeightedGraph<T, U>& operator=(const WeightedGraph<T, U>& other) noexcept;
    constexpr WeightedGraph<T, U>& operator=(const WeightedGraph<T, U>&& other) noexcept;
public:
    constexpr bool print() noexcept;
    constexpr bool add_edge(const unsigned int& ver1, const unsigned int& ver2)  override;
    constexpr bool add_edge(const unsigned int& ver1, const unsigned int& ver2, const U& weight);
};

template <class T, class U>
constexpr bool WeightedGraph<T, U>::add_edge(const unsigned int& ver1, const unsigned int& ver2) {
    throw std::invalid_argument("The weighted class requires an edge weight as third argument.");
    return false;
}

template <class T, class U>
constexpr bool WeightedGraph<T, U>::add_edge(const unsigned int& ver1, const unsigned int& ver2, const U& weight) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    auto& list_to_insert_1 = this->get_adjacency_list()[ver1];
    auto it = std::lower_bound(list_to_insert_1.begin(), list_to_insert_1.end(), std::make_pair(ver2, weight));
    if (it != list_to_insert_1.end() && it->first == ver2) {
        throw std::invalid_argument("Edge already exists");
    }
    list_to_insert_1.insert(it, std::make_pair(ver2, weight));

    auto& list_to_insert_2 = this->get_adjacency_list()[ver2];
    it = std::lower_bound(list_to_insert_2.begin(), list_to_insert_2.end(), std::make_pair(ver1, weight));
    if (it != list_to_insert_2.end() && it->first == ver1) {
        throw std::invalid_argument("Edge already exists");
    }
    list_to_insert_2.insert(it, std::make_pair(ver1, weight));
    return true;
}



template <class T, class U>
WeightedGraph<T, U>::WeightedGraph() : AbstractGraph<T, std::pair<unsigned int, U> >() {}

template <class T, class U>
WeightedGraph<T, U>::WeightedGraph(int num) : AbstractGraph<T, std::pair<unsigned int, U> >(num) {}

template <class T, class U>
WeightedGraph<T, U>::WeightedGraph(const WeightedGraph<T, U>& other) 
            : AbstractGraph<T, std::pair<unsigned int, U> >(other) {}


template <class T, class U>
WeightedGraph<T, U>::WeightedGraph(const WeightedGraph<T, U>&& other) 
            : AbstractGraph<T, std::pair<unsigned int, U> >(other) {}


template <class T, class U>
constexpr WeightedGraph<T, U>& WeightedGraph<T, U>::operator=(const WeightedGraph<T, U>& other) noexcept {
    if (this != &other) {
        this->set_number_of_vertices(other.get_number_of_vertices());
        this->get_vertices().clear();
        this->get_vertices() = other.get_vertices();
        this->get_adjacency_list().clear();
        this->get_adjacency_list() = other.get_adjacency_list();
        return *this;
    }
}

template <class T, class U>
constexpr WeightedGraph<T, U>& WeightedGraph<T, U>::operator=(const WeightedGraph<T, U>&& other) noexcept {
    if (this != &other) {
        this->set_number_of_vertices(other.get_number_of_vertices());
        this->get_vertices.clear();
        this->get_vertices() = std::move(other.get_vertices());
        this->get_adjacency_list().clear();
        this->get_adjacency_list() = std::move(other.get_adjacency_list());
        return *this;
    }
}

template <class T, class U>
constexpr bool WeightedGraph<T, U>::print() noexcept {
    adjacency_list_container<std::pair<unsigned int, U>> tmp = this->get_adjacency_list();
    int i = 0;
    for (const auto& vec : tmp) {
        std::cout << i << " -> ";
        for(const auto& el : vec) {
            std::cout << el.first << "[" << el.second << "]" << " ";
        }
        ++i;
        std::cout << std::endl;
    }
    return true;
}




