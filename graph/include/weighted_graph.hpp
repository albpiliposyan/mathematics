#pragma once

#include "abstract_graph.hpp"


template <class T = char, class U = int>
class WeightedGraph : public AbstractGraph<T, std::pair<int, U> >
{
public:
    WeightedGraph();
    explicit WeightedGraph(int num);
    WeightedGraph(const WeightedGraph<T, U>& other);
    WeightedGraph(WeightedGraph<T, U>&& other) noexcept;
    ~WeightedGraph() = default;
public:
    WeightedGraph<T, U>& operator=(const WeightedGraph<T, U>& other) noexcept;
    WeightedGraph<T, U>& operator=(WeightedGraph<T, U>&& other) noexcept;
public:
    bool add_edge(const int& ver1, const int& ver2)  override;
    virtual bool add_edge(const int& ver1, const int& ver2, const U& weight);
    bool remove_edge(int ver1, int ver2) override;
    bool edge_exists(int v1, int v2);
    std::vector<Edge<U>> get_edges();
    matrix<U> get_adjacency_matrix();
    bool print();
public:
    bool are_vertices_connected(int v, int u);
    WeightedGraph<T, U> minimal_spanning_tree();
    U distance(int startpoint, int endpoint);
    std::vector<int> shortest_path(int startpoint, int endpoint);
    matrix<U> all_pairs_shortest_distances();
    matrix<int> all_pairs_shortest_paths();
    std::vector<WeightedGraph<T, U>> connected_components();
    int number_of_connected_components();
private:
    bool union_sets(std::vector<int>& parent, int x, int y);
    int find(std::vector<int>& parent, int v);
    std::pair<U, std::vector<int>> dijkstra_algorithm(int startpoint, int endpoint);
    std::pair<matrix<U>, matrix<int>> floyd_warshall_algorithm();
    WeightedGraph<T, U> kruskal_algorithm();
    void bfs(int startpoint);
    matrix<int> connected_vertices();
    std::unordered_map<int, int> connected_components_mapping(const matrix<int>&);
};

template <class T, class U>
std::vector<WeightedGraph<T, U>> WeightedGraph<T, U>::connected_components() {
    std::vector<WeightedGraph<T, U>> result;
    const int size = this->get_number_of_vertices();
    if (size == 0) {
        return result;
    }

    matrix<int> components = this->connected_vertices();

    std::unordered_map<int, int> mapping = connected_components_mapping(components);

    for (const std::vector<int>& comp : components) {
        WeightedGraph<T, U> tmp_graph(comp.size());
        for (int c : comp) {
            auto lst = this->get_adjacency_list()[c];
            for (const std::pair<int, U>& el : lst) {
                int mapped_v1 = mapping[c];
                int mapped_v2 = mapping[el.first];
                if (!tmp_graph.edge_exists(mapped_v1, mapped_v2)) {
                    tmp_graph.add_edge(mapped_v1, mapped_v2, el.second);
                }
            }
        }
        result.push_back(tmp_graph);
    }
    return result;
}

template <class T, class U>
std::unordered_map<int, int> WeightedGraph<T, U>::connected_components_mapping(const
                                matrix<int>& components) {
    std::unordered_map<int, int> mapping;
    int counter = 0;
    for (const auto& comp : components) {
        counter =  0;
        for (const int& c : comp) {
            mapping[c] = counter;
            ++counter;
        }
    }
    return mapping;
}


template <class T, class U>
matrix<int> WeightedGraph<T, U>::connected_vertices() {
    matrix<int> components;
    const int size = this->get_number_of_vertices();
    std::vector<bool> visited(size, false);
    int counter = 0;

    while (counter != size) {
        std::vector<int> comp;
        int v = counter;
        visited[v] = true;
        std::queue<int> Q;
        Q.push(counter);
        while (!Q.empty()) {
            v = Q.front();
            Q.pop();
            comp.push_back(v);
            for (std::pair<int, U> el: this->get_adjacency_list()[v]) {
                if (!visited[el.first]) {
                    visited[el.first] = true;
                    Q.push(el.first);
                }
            }
        }
        while (visited[counter]) {
            ++counter;
        }
        components.push_back(comp);
    }
    return components;
}

template <class T, class U>
int WeightedGraph<T, U>::number_of_connected_components() {
    return connected_vertices().size();
}


template <class T, class U>
void WeightedGraph<T, U>::bfs(int startpoint) {
    const int size = this->get_number_of_vertices();
    std::queue<int> Q;
    std::vector<bool> visited(size, false);
    Q.push(startpoint);
    int v = startpoint;
    visited[v] = true;
    while (!Q.empty()) {
        v = Q.front();
        Q.pop();
        std::cout << v << " ";
        for (std::pair<int, U> el: this->get_adjacency_list()[v]) {
            if (!visited[el.first]) {
                visited[el.first] = true;
                Q.push(el.first);
            }
        }
    }
}

template <class T, class U>
bool WeightedGraph<T, U>::are_vertices_connected(int v, int u) {
    int startpoint = v;
    const int size = this->get_number_of_vertices();
    std::queue<int> q;
    std::vector<bool> visited(size, false);
    q.push(startpoint);
    visited[v] = true;
    while (!q.empty()) {
        v = q.front();
        if (v == u) {
            return true;
        }
        q.pop();
        for (std::pair<int, U> el: this->get_adjacency_list()[v]) {
            if (!visited[el.first]) {
                visited[el.first] = true;
                q.push(el.first);
            }
        }
    }
    return false;
}

template <class T, class U>             // TODO print the shortest path
std::pair<U, std::vector<int>> WeightedGraph<T, U>::dijkstra_algorithm(int startpoint, int endpoint) {
    const int n = this->get_number_of_vertices();
    std::vector<bool> visited(n, false);
    std::vector<U> distance(n, std::numeric_limits<U>::max());
    std::vector<int> previous(n, -1);
    distance[startpoint] = 0;
    previous[startpoint] = startpoint;

    for (int i = 0; i < n; ++i) {
        int p = -1;

        for (int j = 0; j < n; ++j) {
            if (!visited[j] && (p == -1 || distance[j] < distance[p])) {
                p = j;
            }
        }
        visited[p] = true;

        for (const std::pair<int, U> el : this->get_adjacency_list()[p]) {
            int v = el.first;
            if (!visited[v]) {
                if (distance[p] + el.second <= distance[v]) {
                    previous[v] = p;
                }
                distance[v] = std::min(distance[p] + el.second, distance[v]);
            }
        }
    }

    std::vector<int> path;
    int v = endpoint;
    while (v != startpoint) {
        path.emplace_back(v);
        v = previous[v];
    }
    path.emplace_back(v);
    std::reverse(path.begin(), path.end());
    return std::make_pair(distance[endpoint], path);
}


template <class T, class U>
std::pair<matrix<U>, matrix<int>> WeightedGraph<T, U>::floyd_warshall_algorithm() {
    const int n = this->get_number_of_vertices();
    matrix<U> adj_matrix = get_adjacency_matrix();
    matrix<U> distance = adj_matrix;
    matrix<U> distance_prev;
    matrix<int> predecessors(n, std::vector<int>(n, std::numeric_limits<int>::max()));
    matrix<int> predecessors_prev;

    U numeric_max = std::numeric_limits<U>::max();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j || distance[i][j] == numeric_max) {
            }
            else {
                predecessors[i][j] = i;
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        distance_prev = distance;
        predecessors_prev = predecessors;
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                if (distance_prev[j][k] <= distance_prev[j][i] + distance_prev[i][k]) {
                    distance[j][k] = distance_prev[j][k];
                    predecessors[j][k] = predecessors_prev[j][k];
                }
                else {
                    distance[j][k] = distance_prev[j][i] + distance_prev[i][k];
                    predecessors[j][k] = predecessors_prev[i][k];
                }
            }
        }
    }

    return {distance, predecessors};
}



template <class T, class U>
U WeightedGraph<T, U>::distance(int startpoint, int endpoint) {
    bool connected = are_vertices_connected(startpoint, endpoint);
    if (!connected) {
        std::string message = "Cannot find the shortest path from " +  std::to_string(startpoint) + " to "
                + std::to_string(endpoint) + ". Vertices are not connected.";
        throw std::invalid_argument(message);
    }
    return this->dijkstra_algorithm(startpoint, endpoint).first;
}


template <class T, class U>
std::vector<int> WeightedGraph<T, U>::shortest_path(int startpoint, int endpoint) {
    return this->dijkstra_algorithm(startpoint, endpoint).second;
}


template <class T, class U>
matrix<U> WeightedGraph<T, U>::all_pairs_shortest_distances() {
    std::pair<matrix<U>, matrix<int>> p = this->floyd_warshall_algorithm();
    return p.first;
}


template <class T, class U>
matrix<int> WeightedGraph<T, U>::all_pairs_shortest_paths() {
    std::pair<matrix<U>, matrix<int>> p = this->floyd_warshall_algorithm();
    return p.second;
}


template <class T, class U>
WeightedGraph<T, U> WeightedGraph<T, U>::kruskal_algorithm() {
    std::vector<Edge<U>> edges = this->get_edges();

    std::sort(edges.begin(), edges.end(), [](const Edge<U>& a, const Edge<U>& b) {
        if (a.value == b.value) {
            return a.v1 < b.v1;
        }
        return a.value < b.value;
    });

    WeightedGraph<T, U> spanning_tree(this->get_number_of_vertices());

    std::vector<int> parent(this->get_number_of_vertices(), -1);
    for (const auto& edge : edges) {
        int u = edge.v1;
        int v = edge.v2;

        int set_u = find(parent, u);
        int set_v = find(parent, v);

        if (set_u != set_v) {
            spanning_tree.add_edge(u, v, edge.value);
            union_sets(parent, set_u, set_v);
        }
    }
    return spanning_tree;
}


template <class T, class U>
WeightedGraph<T, U> WeightedGraph<T, U>::minimal_spanning_tree() {
    if (number_of_connected_components() != 1) {
        std::cerr << "Graph is not connected and does not have a minimal spanning tree.";
        exit(1);
    }
    return kruskal_algorithm();
}


template <class T, class U>
int WeightedGraph<T, U>::find(std::vector<int>& parent, int v) {
    if (parent[v] == -1) {
        return v;
    }
    return find(parent, parent[v]);
}


template<class T, class U>
bool WeightedGraph<T, U>::union_sets(std::vector<int> &parent, int x, int y) {
    int root_x = find(parent, x);
    int root_y = find(parent, y);
    parent[root_x] = root_y;
    return true;
}


template<class T, class U>
std::vector<Edge<U>> WeightedGraph<T, U>::get_edges() {
    std::vector<Edge<U>> edges;
    edges.reserve(this->get_number_of_edges());
    int i = 0;
    for (auto& vec : this->get_adjacency_list()) {
        for (auto& el : vec) {
            Edge<U> tmp_edge;
            tmp_edge.v1 = i;
            tmp_edge.v2 = el.first;
            tmp_edge.value = el.second;
            edges.emplace_back(tmp_edge);
        }
        ++i;
    }
    return edges;
}


template <class T, class U>
matrix<U> WeightedGraph<T, U>::get_adjacency_matrix() {
    const int n = this->get_number_of_vertices();
    matrix<U> adj_matrix(n, std::vector<U>(n, std::numeric_limits<U>::max()));
    for (int i = 0; i < n; ++i) {
        adj_matrix[i][i] = U();
    }
    int i = 0;
    for (const std::list<std::pair<int, U>> v : this->get_adjacency_list()) {
        for (const std::pair<int, U> u : v) {
            adj_matrix[i][u.first] = u.second;
        }
        ++i;
    }

    return adj_matrix;
}


template <class T, class U>
bool WeightedGraph<T, U>::add_edge(const int& ver1, const int& ver2) {
    throw std::invalid_argument("The weighted class requires an edge weight as third argument.");
}


template <class T, class U>
bool WeightedGraph<T, U>::add_edge(const int& ver1, const int& ver2, const U& weight) {
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
bool WeightedGraph<T, U>::remove_edge(int ver1, int ver2) {
    if (ver1 >= this->get_number_of_vertices() || ver2 >= this->get_number_of_vertices()) {
        throw std::invalid_argument("Invalid vertices provided");
    }
    for (auto& it : this->get_adjacency_list()[ver1]) {
        if (it.first == ver2) {
            this->get_adjacency_list()[ver1].remove(it);
            break;
        }
    }
    for (auto& it : this->get_adjacency_list()[ver2]) {
        if (it.first == ver1) {
            this->get_adjacency_list()[ver2].remove(it);
            break;
        }
    }
    return true;
}


template <class T, class U>
bool WeightedGraph<T, U>::edge_exists(int v1, int v2) {
    auto& adj_list = this->get_adjacency_list()[v1];
    for (auto it = adj_list.begin(); it != adj_list.end(); it++) {
        if (it->first == v2) {
            return true;
        }
    }
    return false;
}


template <class T, class U>
bool WeightedGraph<T, U>::print() {
    adjacency_list_container<std::pair<int, U>> tmp = this->get_adjacency_list();
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


template <class T, class U>
WeightedGraph<T, U>::WeightedGraph()
                    : AbstractGraph<T, std::pair<int, U> >() {}


template <class T, class U>
WeightedGraph<T, U>::WeightedGraph(int num)
                    : AbstractGraph<T, std::pair<int, U> >(num) {}


template <class T, class U>
WeightedGraph<T, U>::WeightedGraph(const WeightedGraph<T, U>& other)
            : AbstractGraph<T, std::pair<int, U> >(other) {}


template <class T, class U>
WeightedGraph<T, U>::WeightedGraph(WeightedGraph<T, U>&& other) noexcept
            : AbstractGraph<T, std::pair<int, U> >(other) {}


template <class T, class U>
WeightedGraph<T, U>& WeightedGraph<T, U>::operator=(const WeightedGraph<T, U>& other) noexcept {
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
WeightedGraph<T, U>& WeightedGraph<T, U>::operator=(WeightedGraph<T, U>&& other) noexcept {
    if (this != &other) {
        this->set_number_of_vertices(other.get_number_of_vertices());
        this->get_vertices.clear();
        this->get_vertices() = std::move(other.get_vertices());
        this->get_adjacency_list().clear();
        this->get_adjacency_list() = std::move(other.get_adjacency_list());
        return *this;
    }
}

