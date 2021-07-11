#include <igraph/igraph.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>


//FIXME Have a better selection strategy and use igraph_bfs to compute eccentricity and distance at the same time 

std::vector<int> bdgecc(igraph_t V)
{
    std::vector<int> eccs(igraph_vcount(&V), 0); // FIXME may be useless
    std::vector<int> eccs_lower(igraph_vcount(&V), std::numeric_limits<int>::min());
    std::vector<int> eccs_upper(igraph_vcount(&V), std::numeric_limits<int>::max());
    int number_of_bfs = 0;
    std::vector<int> candidate(igraph_vcount(&V), 1);
    int candidates = candidate.size();
    while((std::find(candidate.begin(), candidate.end(), 1)) != candidate.end())
    {
        igraph_vector_t ecc_vec;
	igraph_vector_init(&ecc_vec, 1);
        int v_id = std::find(candidate.begin(), candidate.end(), 1) - candidate.begin(); // FIXME need better strategy
       	igraph_eccentricity(&V, &ecc_vec, igraph_vss_1(v_id), IGRAPH_OUT); // FIXME code ourself eccentricity

	number_of_bfs++;
        int v_ecc = igraph_vector_e(&ecc_vec, 0);
        eccs[v_id] = v_ecc;

        for (int w_id = 0; w_id < candidate.size(); w_id += 1)
        {
	    if (!candidate[w_id])
	        continue;
            igraph_vs_t v = igraph_vss_1(v_id);
            igraph_vs_t w = igraph_vss_1(w_id);
            igraph_matrix_t mat_res;
            igraph_matrix_init(&mat_res, 1, 1);
            igraph_shortest_paths(&V, &mat_res , v, w, IGRAPH_OUT);
            int d_vw = MATRIX(mat_res, 0, 0);
            eccs_lower[w_id] = std::max(eccs_lower[w_id], std::max(v_ecc - d_vw, d_vw));
            eccs_upper[w_id] = std::min(eccs_upper[w_id], v_ecc + d_vw);
            if (eccs_lower[w_id] == eccs_upper[w_id])
            {
                eccs[w_id] = eccs_lower[w_id];
		candidate[w_id] = 0;
            }
        }

    }
    //std::cout << number_of_bfs;
    return eccs;
}
