#include <igraph/igraph.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <algorithm>



int bdgecc(igraph_t V)
{
    std::vector<int> eccs(igraph_vcount(&V), 0);
    std::vector<int> eccs_lower(igraph_vcount(&V), std::numeric_limits<int>::min());
    std::vector<int> eccs_upper(igraph_vcount(&V), std::numeric_limits<int>::max());
    int delta_L = std::numeric_limits<int>::min(); 
    int delta_U = std::numeric_limits<int>::max();
    std::vector<int> candidate(igraph_vcount(&V), 1);
    while((std::find(candidate.begin(), candidate.end(), 1)) != candidate.end())
    {
        igraph_vector_t ecc_vec;
	igraph_vector_init(&ecc_vec, 1);
        int v_id = std::find(candidate.begin(), candidate.end(), 1) - candidate.begin(); // FIXME need real strategy
       	igraph_eccentricity(&V, &ecc_vec, igraph_vss_1(v_id), IGRAPH_OUT);
        int v_ecc = igraph_vector_e(&ecc_vec, 0);
        eccs[v_id] = v_ecc;

	delta_L = std::max(delta_L, v_ecc);
	delta_U = std::max(delta_U, 2 * v_ecc);

        igraph_vs_t all = igraph_vss_all();
        igraph_vit_t iterator;
        igraph_vit_create(&V, all, &iterator);
        while(!IGRAPH_VIT_END(iterator))
        {
            int w_id = IGRAPH_VIT_GET(iterator);
	    IGRAPH_VIT_NEXT(iterator);
            igraph_vs_t v = igraph_vss_1(v_id);
            igraph_vs_t w = igraph_vss_1(w_id);
            igraph_matrix_t mat_res;
            igraph_matrix_init(&mat_res, 1, 1);
            igraph_shortest_paths(&V, &mat_res , v, w, IGRAPH_OUT);
            int d_vw = MATRIX(mat_res, 0, 0);
            eccs_lower[w_id] = std::max(eccs_lower[w_id], std::max(v_ecc - d_vw, d_vw));
            eccs_upper[w_id] = std::min(eccs_upper[w_id], v_ecc + d_vw);
            if (eccs_lower[w_id] == eccs_upper[w_id] || (eccs_upper[w_id] <= delta_L && eccs_lower[w_id] <= delta_U / 2))
            {
                eccs[w_id] = eccs_lower[w_id];
		candidate[w_id] = 0;
                /*
                FIXME: deleting an element in the graph
                might change the ids of the vertices and so our
                eccs lists would stop reference the right vertice.
                A fix would be to only append in eccs (tuple id, ecc)
                and then remove w from both eccs_upper and eccs_lower
                (hoping that removing a vertex from a graph only does a
                shift in the indexes, like in a vector).
                */
            }
        }

    }
    return delta_L;
}
