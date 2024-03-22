#ifndef DIRECTED_WEIGHTED_GRAPH_H
#define DIRECTED_WEIGHTED_GRAPH_H

#include "weighted_graph.hpp"

template <class T, class U>
class DirectedWeightedGraph : public WeightedGraph<T, U>
{
public:
    DirectedWeightedGraph();
    explicit DirectedWeightedGraph(int num);
    DirectedWeightedGraph(const DirectedWeightedGraph<T, U>& other);
    DirectedWeightedGraph(DirectedWeightedGraph<T, U>&& other) noexcept ;
    ~DirectedWeightedGraph() = default;
public:
    bool add_edge(const int& ver1, const int& ver2)  override;
    bool add_edge(const int& ver1, const int& ver2, const U& weight) override;
    bool remove_edge(int ver1, int ver2) override;
public:
    bool has_cycle();
private:
    bool has_cycle_helper(int v, std::vector<bool> visited, int parent);
};


template <class T, class U>
bool DirectedWeightedGraph<T, U>::has_cycle_helper(int v, std::vector<bool> visited, int parent) {
    visited[v] = true;
    for (const auto& neighbor : this->get_adjacency_list()[v]) {
        int next_vertex = neighbor.first;

        if (!visited[next_vertex]) {
            if (has_cycle_helper(next_vertex, visited, v)) {
                return true;
            }
            else if (next_vertex != parent) {
                return true;
            }
        }
    }
    return false;
    // if (visited[v]) {
    //     return true;
    // }
    // visited[v] = true;
    // bool FLAG = false;
    // for (auto& el : this->get_adjacency_list()[v]) {
    //     FLAG = has_cycle_helper(el.first, visited);
    //     if (FLAG) {
    //         return true;
    //     }

    // }
    // return false;
}

template <class T, class U>
bool DirectedWeightedGraph<T, U>::has_cycle() {
    std::vector<bool> visited(this->get_number_of_vertices(), false);

    for (int i = 0; i < this->get_number_of_vertices(); ++i) {
        if (!visited[i]) {
            if (has_cycle_helper(i, visited, -1)) {
                return true;
            }
        }
    }
    return false;
}


template <class T, class U>
DirectedWeightedGraph<T, U>::DirectedWeightedGraph() : WeightedGraph<T, U>() {}

template <class T, class U>
DirectedWeightedGraph<T, U>::DirectedWeightedGraph(int num)
        : WeightedGraph<T, U>(num) {}

template <class T, class U>
DirectedWeightedGraph<T, U>::DirectedWeightedGraph(const DirectedWeightedGraph<T, U>& other)
        : WeightedGraph<T, U>(other) {}


template <class T, class U>
DirectedWeightedGraph<T, U>::DirectedWeightedGraph(DirectedWeightedGraph<T, U>&& other) noexcept
        : WeightedGraph<T, U>(other) {}


template <class T, class U>
bool DirectedWeightedGraph<T, U>::add_edge(const int& ver1, const int& ver2) {
    throw std::invalid_argument("The weighted class requires an edge weight as third argument.");
}


template <class T, class U>
bool DirectedWeightedGraph<T, U>::add_edge(const int& ver1, const int& ver2, const U& weight) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    auto& list_to_insert_1 = this->get_adjacency_list()[ver1];
    auto it = std::lower_bound(list_to_insert_1.begin(), list_to_insert_1.end(), std::make_pair(ver2, weight));
    if (it != list_to_insert_1.end() && it->first == ver2) {
        throw std::invalid_argument("Edge already exists");
    }
    list_to_insert_1.insert(it, std::make_pair(ver2, weight));

    return true;
}


template <class T, class U>
bool DirectedWeightedGraph<T, U>::remove_edge(int ver1, int ver2) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    for (auto& it : this->get_adjacency_list()[ver1]) {
        if (it.first == ver2) {
            this->get_adjacency_list()[ver1].remove(it);
            break;
        }
    }
    return true;
}

#endif //DIRECTED_WEIGHTED_GRAPH_H
