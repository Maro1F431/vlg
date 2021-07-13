#include "bdgecc.hh"
#include <iostream>
#include <vector>
int main()
{
    igraph_t g;
    igraph_vector_t v1, v2;

    /* simple use */
    igraph_vector_init(&v1, 16);
    VECTOR(v1)[0] = 0;
    VECTOR(v1)[1] = 1;
    VECTOR(v1)[2] = 1;
    VECTOR(v1)[3] = 2;
    VECTOR(v1)[4] = 2;
    VECTOR(v1)[5] = 3;
    VECTOR(v1)[6] = 2;
    VECTOR(v1)[7] = 2;
    VECTOR(v1)[8] = 2;

    igraph_create(&g, &v1, 0, 0);
    std::vector<int> eccs = bdgecc(g);
    for (auto& i: eccs)
        std::cout << i << ' ';
    igraph_destroy(&g);
}
