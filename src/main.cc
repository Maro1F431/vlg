#include "bdgecc.hh"
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <algorithm>

int main()
{
    igraph_t g;
    igraph_vector_t v1;

    /* simple use */
    igraph_vector_init(&v1, 10);
    VECTOR(v1)[0] = 0;
    VECTOR(v1)[1] = 1;
    VECTOR(v1)[2] = 1;
    VECTOR(v1)[3] = 2;
    VECTOR(v1)[4] = 2;
    VECTOR(v1)[5] = 3;
    VECTOR(v1)[6] = 2;
    VECTOR(v1)[7] = 2;
    VECTOR(v1)[8] = 2;

    FILE *file;

    file = fopen ("../email-Enron.txt" , "r");

    igraph_t graph;

    igraph_read_graph_edgelist(&graph, file,
                               0, false);

    std::cout << "Graph created";
    fclose(file);
    igraph_create(&g, &v1, 0, 0);
    std::vector<int> eccs = bdgecc(graph);
    for (auto& i: eccs)
        std::cout << i << ' ';
    igraph_vector_destroy(&v1);
    std::cout << '\n' << "Calculated :" << *max_element(eccs.begin(), eccs.end());
    igraph_real_t res;

    igraph_diameter(&graph, &res,
                    NULL, NULL,
                    NULL,
                    false, true);
    std::cout << '\n' << "Ground truth: " << res;

    igraph_destroy(&g);
    igraph_destroy(&graph);
 }
