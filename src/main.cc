#include "bdgecc.hh"
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <experimental/string_view>

int main(int argc, char *argv[])
{
    /* simple use */
    int indexFile = 1;
    if (argc == 1){
        std::cerr << "This programm needs a file\n";
	return -1;
    }
    
    FILE *file;
    bool print_ecc = false;
    if (std::experimental::string_view(argv[1]) == "--print"){
        print_ecc = true;
	indexFile++;
    }

    file = fopen (argv[indexFile] , "r");

    igraph_t graph;

    igraph_read_graph_edgelist(&graph, file,
                               0, false);

    fclose(file);
    std::cout << "Graph created \n";
    std::vector<int> eccs = bdgecc(graph);
    if (print_ecc){
	std::cout << "List of eccentricities:\n";
        for (auto& i: eccs)
            std::cout << i << ' ';
    }
    return 0;
 }
