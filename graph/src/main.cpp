
#include <iostream>
#include "weighted_graph.hpp"


int main() {
    std::cout << "Not Connected Weighted Graph:\n";
    int n = 6;
    WeightedGraph<int, double> not_connected_wg(n);                   // Not Connected Graph Example
    not_connected_wg.add_edge(0, 1, 1);
    not_connected_wg.add_edge(0, 2, 4);
    not_connected_wg.add_edge(1, 5, 1);
    not_connected_wg.add_edge(2, 5, 1);
    not_connected_wg.add_edge(4, 3, 3);
    not_connected_wg.print();
    std::cout << "\n";

    std::cout << "Connected components:\n";
    int ix = 1;
    std::vector<WeightedGraph<int, double>> components =  not_connected_wg.connected_components();     // Connected Components
    for (WeightedGraph<int, double>& comp : components) {
        std::cout << "Component " << ix++ << ":\n";
        comp.print();
    }
    std::cout << "\n\n";



    std::cout << "Weighted Graph:\n";
    n = 12;
    WeightedGraph<int, double> wg(n);					// Connected Graph Example
    wg.add_edge(0, 1, 1);
    wg.add_edge(0, 4, 4);
    wg.add_edge(1, 2, 3);
    wg.add_edge(1, 4, 1);
    wg.add_edge(2, 5, 2);
    wg.add_edge(3, 6, 1);
    wg.add_edge(3, 7, 2);
    wg.add_edge(4, 5, 6);
    wg.add_edge(4, 8, 1);
    wg.add_edge(5, 6, 3);
    wg.add_edge(5, 9, 1);
    wg.add_edge(5, 10, 2);
    wg.add_edge(6, 7, 1);
    wg.add_edge(6, 11, 3);
    wg.add_edge(7, 11, 4);
    wg.add_edge(8, 9, 3);
    wg.add_edge(8, 10, 2);
    wg.add_edge(9, 10, 2);
    wg.add_edge(10, 11, 3);
    wg.print();
    std::cout << "\n";

    std::cout << "Minimal Spanning tree:\n";
    WeightedGraph<int, double> spanning_tree = wg.minimal_spanning_tree();              // Minimal Spanning Tree
    spanning_tree.print();

    std::cout << "\n";                                                                  // Distance Dijkstra
    const int startpoint = 0;
    const int endpoint = 3;
    std::cout << "Distance from " << startpoint << " to " << endpoint << ": ";
    std::cout << wg.distance(startpoint, endpoint) << "\n";

    std::cout << "Shortest path from " << startpoint << " to " << endpoint << ": ";     // Shortest Path
    const std::vector<int> shortest_path = wg.shortest_path(startpoint, endpoint);
    for (int v : shortest_path) {
        std::cout << v << " ";
    }
    std::cout << "\n\n";

    matrix<double> distance = wg.all_pairs_shortest_distances();                       	// all pairs shortest path Floyd_Warshall algorithm
    matrix<int> predecessors = wg.all_pairs_shortest_paths();
    double numeric_max = std::numeric_limits<double>::max();

    std::cout << "All Pairs Shortest Distance:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (distance[i][j] == numeric_max) {
                std::cout << "i ";
            }
            else {
                std::cout << distance[i][j] << " ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "All Pairs Shortest Predecessors:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (predecessors[i][j] == std::numeric_limits<int>::max()) {
                std::cout << "inf ";
            }
            else {
                std::cout << predecessors[i][j] << " ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";


    return 0;
}
