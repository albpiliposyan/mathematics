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
    // TODO
    constexpr bool add_edge(const unsigned int& ver1, const unsigned int& ver2) override;
    constexpr bool set_number_of_vertices(const unsigned int num) noexcept override;
};

// template <class T, class U>
// constexpr bool WeightedGraph<T, U>::add_edge(const unsigned int& ver1, const unsigned int& ver2) {return true;}
// 
// template <class T, class U>
// constexpr bool WeightedGraph<T, U>::set_number_of_vertices(const unsigned int num) noexcept {return true;}

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




