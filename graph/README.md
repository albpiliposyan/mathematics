# Graph Project

This is a C++ project for working with graphs.

Note: This project is still a work in progress and more functionality is being added.

## Requirements

The Matrix class requires a C++ compiler that supports C++20 or later.

## Project Structure

The project contains the following header files:

- `abstract_graph.h`: Defines an abstract graph class that can be inherited by concrete graph implementations.
- `graph.h`: Implements a concrete undirected graph class that uses an adjacency list representation.
- `directed_graph.h`: Implements a concrete directed graph class that uses an adjacency list representation.

## Usage

To use the project, simply include the necessary header files in your C++ code and instantiate the desired graph class. The following is an example of creating a `Graph` object with `int` vertex data type:

```c++
#include "graph.h"

int main() {
    Graph<int> graph;
    graph.add_vertex(1);
    graph.add_vertex(2);
    graph.add_edge(1, 2);

    return 0;
}
```
