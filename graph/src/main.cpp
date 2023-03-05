#include <iostream>
#include "graph.h"
#include "directed_graph.h"


int main() {
    Graph<char> graph(5);
    graph.add_edge(0, 2);
    graph.add_edge(0, 4);
    graph.add_edge(1, 3);
    graph.add_edge(2, 3);
    graph.add_vertice('a');
    graph.print();
    std::cout << std::endl;

    DirectedGraph<char> directed_graph(5);
    directed_graph.add_edge(0, 1);
    directed_graph.add_edge(0, 2);
    directed_graph.add_edge(0, 4);
    directed_graph.add_edge(1, 4);
    directed_graph.add_edge(2, 3);
    directed_graph.add_edge(4, 2);
    directed_graph.print();


    return 0;
}