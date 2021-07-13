#include <igraph/igraph.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>
#include <stdbool.h>


//FIXME Have a better selection strategy and use igraph_bfs to compute eccentricity and distance at the same time 

void free_complist(igraph_vector_ptr_t *complist, long int index) {
    long int i;
    for (i = 0; i < igraph_vector_ptr_size(complist); i++) {
	if (i == index)
	    continue;
        igraph_destroy((igraph_t*)VECTOR(*complist)[i]);
        igraph_free((igraph_t*)VECTOR(*complist)[i]);
    }
}

std::vector<int> bdgecc(igraph_t V)
{
    // Picking the right component
    igraph_vector_ptr_t complist;
    igraph_vector_ptr_init(&complist, 0);
    igraph_decompose(&V, &complist, IGRAPH_WEAK, -1, 2);
    long int max_index = 0;
    igraph_integer_t max_vertice = 0;

    long int n = igraph_vector_ptr_size(&complist);
    for (long int i = 0; i < n; i++) {
        igraph_t part = *((igraph_t*)(VECTOR(complist)[i]));
	auto size = igraph_vcount(&part);

	if (size > max_vertice){
	    max_vertice = size;
	    max_index = i;
	}
    }

    //std::cout << max_index;

    free_complist(&complist, max_index);

    V = *((igraph_t*)(VECTOR(complist)[max_index]));

    std::vector<int> eccs(igraph_vcount(&V), 0); // FIXME may be useless
    std::vector<int> eccs_lower(igraph_vcount(&V), std::numeric_limits<int>::min());
    std::vector<int> eccs_upper(igraph_vcount(&V), std::numeric_limits<int>::max());
    int number_of_bfs = 0;
    std::vector<int> candidate(igraph_vcount(&V), 1);
    int candidates = candidate.size();
    int v_id = -1;
    bool high = true;
    int maxuppernode = -1;
    int minlowernode = -1;

    while(candidates > 0)
    {
	high = !high;
        if(v_id == -1) { // only in the first round, select node with highest degree
            v_id = 0;
	    auto v_neigh = -1;
            for (int i = 0; i < candidate.size(); i++){
	        igraph_vector_t i_id_neigh;
	        igraph_vector_init(&i_id_neigh, 1);
	        igraph_neighborhood_size(&V, &i_id_neigh,
                             igraph_vss_1(v_id), 1,
                             IGRAPH_ALL,
                             1);
		auto i_neigh = VECTOR(i_id_neigh)[0];
		if (v_neigh == i_neigh)
	            v_neigh = i_neigh;

                if(i_neigh > v_neigh) {
                    v_id = i;
		    v_neigh = i_neigh;
                }
            }
        }
        else if(high) // select node with highest upper bound
            v_id = maxuppernode;
        else // select node with lowest lower bound
            v_id = minlowernode;

	maxuppernode = -1;
	minlowernode = -1;
       	//igraph_eccentricity(&V, &ecc_vec, igraph_vss_1(v_id), IGRAPH_OUT); // FIXME code ourself eccentricity
	
	igraph_vector_t dist;


    	/* Initialize the vectors where the result will be stored. Any of these
     	* can be omitted and replaced with a null pointer when calling
     	* igraph_bfs() */
    	igraph_vector_init(&dist, 0);

    	/* Now call the BFS function */
    	igraph_bfs(&V, /*root=*/v_id, /*roots=*/ NULL, /*neimode=*/ IGRAPH_ALL,
               /*unreachable=*/ 0, /*restricted=*/ 0,
               NULL, NULL, NULL, NULL, NULL, &dist,
               /*callback=*/ 0, /*extra=*/ 0);

	int v_ecc = (int) igraph_vector_max(&dist);
	number_of_bfs++;
        eccs[v_id] = v_ecc;

        for (int w_id = 0; w_id < candidate.size(); w_id += 1)
        {
	    if (!candidate[w_id])
	        continue;
            int d_vw = VECTOR(dist)[w_id];
            eccs_lower[w_id] = std::max(eccs_lower[w_id], std::max(v_ecc - d_vw, d_vw));
            eccs_upper[w_id] = std::min(eccs_upper[w_id], v_ecc + d_vw);
            if (eccs_lower[w_id] == eccs_upper[w_id]){
                eccs[w_id] = eccs_lower[w_id];
		candidate[w_id] = 0;
		candidates--;
            }

	    // Compute max upper node and min lower node for next selection 
	   
	    igraph_vector_t w_id_neigh;
	    igraph_vector_init(&w_id_neigh, 1);
	    igraph_neighborhood_size(&V, &w_id_neigh,
                             igraph_vss_1(w_id), 1,
                             IGRAPH_ALL,
                             1);


            if(minlowernode == -1)
                minlowernode = w_id;
            else if(eccs_lower[w_id] == eccs_lower[minlowernode]){
                igraph_vector_t minlowernode_neigh;
	    	igraph_vector_init(&minlowernode_neigh, 1);
	    	igraph_neighborhood_size(&V, &minlowernode_neigh,
                             igraph_vss_1(minlowernode), 1,
                             IGRAPH_ALL,
                             1);
		if (VECTOR(w_id_neigh)[0] > VECTOR(minlowernode_neigh)[0])
			minlowernode = w_id;
	    }
            else if(eccs_lower[w_id] < eccs_lower[minlowernode])
                minlowernode = w_id;
            if(maxuppernode == -1)
                maxuppernode = w_id;
            else if(eccs_upper[w_id] == eccs_upper[maxuppernode]){
                igraph_vector_t maxuppernode_neigh;
	    	igraph_vector_init(&maxuppernode_neigh, 1);
	    	igraph_neighborhood_size(&V, &maxuppernode_neigh,
                             igraph_vss_1(maxuppernode), 1,
                             IGRAPH_ALL,
                             1);
		if (VECTOR(w_id_neigh)[0] > VECTOR(maxuppernode_neigh)[0])
			minlowernode = w_id;
	    }
            else if(eccs_upper[w_id] > eccs_upper[maxuppernode])
                maxuppernode = w_id;
        }

    }
    //std::cout << number_of_bfs;
    return eccs;
}
