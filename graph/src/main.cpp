#include <iostream>
#include <string>
#include "graph.hpp"
#include "weighted_graph.hpp"
#include "directed_graph.hpp"


int main() {
    std::cout << "Graph" << std::endl;
    Graph<char> g(9);
    g.add_edge(0, 1);
    g.add_edge(1, 2);
    g.add_edge(2, 3);
    g.add_edge(2, 4);
    g.add_edge(2, 5);
    g.add_edge(4, 6);
    g.add_edge(5, 7);
    g.add_edge(7, 8);
    g.print();
    std::cout << "Cayley representation: ";
    std::vector<unsigned int> cayley_result = g.cayley_representation();
    for (unsigned int& el : cayley_result) {
        std::cout << el;
    }
    std::cout << std::endl;
    std::cout << std::endl;


    std::cout << "Directed Graph" << std::endl;
    DirectedGraph<char> dg(5);
    dg.add_edge(0, 1);
    dg.add_edge(0, 2);
    dg.add_edge(0, 4);
    dg.add_edge(1, 4);
    dg.add_edge(2, 3);
    dg.add_edge(4, 2);
    dg.print();
    std::cout << std::endl;

    std::cout << "Weighted Graph" << std::endl;
    WeightedGraph<int, std::string> wg(5);
    wg.add_edge(0, 1, "test1");
    wg.add_edge(0, 2, "test2");
    wg.add_edge(1, 2, "test3");
    wg.add_edge(3, 4, "test4");

    wg.print();

    return 0;
}
