#include <igraph/igraph.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>
#include <stdbool.h>
#include <numeric>

/*Function used to free igraph_vector_ptr in order to have 0 memory leaks*/
void free_complist(igraph_vector_ptr_t *complist, long int index) {
    long int i;
    for (i = 0; i < igraph_vector_ptr_size(complist); i++) {
	if (i == index)
	    continue;
        igraph_destroy((igraph_t*)VECTOR(*complist)[i]);
        igraph_free((igraph_t*)VECTOR(*complist)[i]);
    }
}

void computeCommunities(igraph_t* graph, igraph_vector_t* membership, igraph_integer_t *nb_cluster){
    igraph_rng_seed(igraph_rng_default(), 0);

    igraph_community_leiden(graph, NULL, NULL, 0.05, 0.01, 0, membership, nb_cluster, NULL);
}

int selectFromCommunitiesCycleOnce(std::vector<std::pair<long int, long int>> delegate, int* commuIndex, int* nb_cluster){
    if (*commuIndex > *nb_cluster)
        return -1;
    else{
    	auto res = delegate[*commuIndex].first;
	*commuIndex += 1;
	return res;
    }
        
}


int selectFrom(igraph_t V, int v_id, int maxuppernode, int minlowernode, std::vector<int> candidate, bool high, bool biggest, igraph_vector_t* membership, igraph_integer_t* nb_cluster, int* commuIndex, std::vector<std::pair<long int, long int>> delegate){ 
    int res_id = 0;
    /* only in the first round, select node with highest degree*/
    if(v_id == -1) {
    	auto res_neigh = -1;
        for (int i = 0; i < candidate.size(); i++){
	    igraph_vector_t i_id_neigh;
	    igraph_vector_init(&i_id_neigh, 1);
	    igraph_neighborhood_size(&V, &i_id_neigh,
                        igraph_vss_1(res_id), 1,
                        IGRAPH_ALL,
                        1);
            auto i_neigh = VECTOR(i_id_neigh)[0];
	    if (res_neigh == i_neigh)
	        res_neigh = i_neigh;

            if(i_neigh > res_neigh) {
                res_id = i;
		res_neigh = i_neigh;
            }
	    igraph_vector_destroy(&i_id_neigh);
         }
	 return res_id;
    }
    else{
	/*high is used to switch between lower bound and higher bound for better result*/
        res_id = selectFromCommunitiesCycleOnce(delegate, commuIndex, nb_cluster);	
	if (res_id == -1){
            if(high) /*select node with highest upper bound*/
                return maxuppernode;
            else /*select node with lowest lower bound*/
                return minlowernode;
	}
	return res_id;
    }
}

/* Main function */
std::vector<int> bdgecc(igraph_t V)
{
    igraph_vector_ptr_t complist;
    igraph_vector_ptr_init(&complist, 0);

    /* Calculation of LWCC*/
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

    free_complist(&complist, max_index);

    V = *((igraph_t*)(VECTOR(complist)[max_index]));

    /*Initialisation of the variables for the algorithm*/

    std::vector<int> eccs(igraph_vcount(&V), 0);
    std::vector<int> eccs_lower(igraph_vcount(&V), std::numeric_limits<int>::min());
    std::vector<int> eccs_upper(igraph_vcount(&V), std::numeric_limits<int>::max());
    int number_of_bfs = 0;
    std::vector<int> candidate(igraph_vcount(&V), 1);
    int candidates = candidate.size();
    int v_id = -1;
    int maxuppernode = -1;
    int minlowernode = -1;
    int commuIndex = 0;

    bool high = true;
    bool biggest = true;

    igraph_vector_t membership;
    igraph_vector_init(&membership, igraph_vcount(&V));
    igraph_integer_t nb_cluster;
    computeCommunities(&V, &membership, &nb_cluster);
    
    
    std::vector<std::pair<long int, long int>> delegate(nb_cluster, std::make_pair(-1, -1));

    for(int i = 0; i < igraph_vcount(&V); i++){
        igraph_vector_t i_neigh;
	igraph_vector_init(&i_neigh, 1);
	igraph_neighborhood_size(&V, &i_neigh,
                             igraph_vss_1(i), 1,
                             IGRAPH_ALL,
                             1);
	auto currentCommu = VECTOR(membership)[i];
	auto p = delegate[currentCommu];
	if (VECTOR(i_neigh)[0] > p.second){
	    delegate[currentCommu] = std::make_pair(i, VECTOR(i_neigh)[0]);
	}
    }

    /*Takes-Kost algorithm */
    while(candidates > 0)
    {
	high = !high;
	biggest = !biggest;
        v_id = selectFrom(V, v_id, maxuppernode, minlowernode, candidate, high, biggest, &membership, &nb_cluster, &commuIndex, delegate);
	maxuppernode = -1;
	minlowernode = -1;
	
	igraph_vector_t dist;


    	/* Computation of eccentricity using BFS*/
    	igraph_vector_init(&dist, 0);

    	igraph_bfs(&V, /*root=*/v_id, /*roots=*/ NULL, /*neimode=*/ IGRAPH_ALL,
               /*unreachable=*/ 0, /*restricted=*/ 0,
               NULL, NULL, NULL, NULL, NULL, &dist,
               /*callback=*/ 0, /*extra=*/ 0);

	number_of_bfs++;

	int v_ecc = (int) igraph_vector_max(&dist);
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

	    /* Compute max upper node and min lower node for next selection */
	   
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

		igraph_vector_destroy(&minlowernode_neigh);
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
		igraph_vector_destroy(&maxuppernode_neigh);
	    }
            else if(eccs_upper[w_id] > eccs_upper[maxuppernode])
                maxuppernode = w_id;
	    igraph_vector_destroy(&w_id_neigh);    
        }
	igraph_vector_destroy(&dist);

    }

    /* Destroying in order to avoid memory leaks */
    igraph_destroy((igraph_t*)VECTOR(*(&complist))[max_index]);
    igraph_free((igraph_t*)VECTOR(*(&complist))[max_index]);
    igraph_vector_ptr_destroy(&complist);
    igraph_vector_destroy(&membership);

    std::cout << '\n' << "Number of bfs: " << number_of_bfs << "\n\n";
    return eccs;
}
