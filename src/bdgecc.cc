#include <igraph/igraph.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>
#include <stdbool.h>
#include <numeric>



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
    //igraph_vector_t  degree;
    //igraph_real_t quality;

    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);

    /* Perform Leiden algorithm using CPM */
    //igraph_vector_init(&membership, igraph_vcount(&graph)); TODO Outside of this fonction
    igraph_community_leiden(graph, NULL, NULL, 0.05, 0.01, 0, membership, nb_cluster, NULL);

    //printf("Membership: ");
    //igraph_vector_print(membership);
    //printf("\n");

    /* Start from existing membership to improve it further */
    //igraph_community_leiden(&graph, NULL, NULL, 0.05, 0.01, 1, &membership, &nb_clusters, &quality);

    //printf("Iterated Leiden, using CPM (resolution parameter 0.05), quality is %.4f.\n", quality);
    printf("Membership: ");
    igraph_vector_print(membership);
    std::cout << " \n Mathieu  : " << igraph_vector_max(membership) << "nb_cluster: " << *nb_cluster;
    printf("\n");

    /* Initialize degree vector to use for optimizing modularity */
    //igraph_vector_init(&degree, igraph_vcount(&graph));
    //igraph_degree(&graph, &degree, igraph_vss_all(), IGRAPH_ALL, 1);

    /* Perform Leiden algorithm using modularity */
    //igraph_community_leiden(&graph, NULL, &degree, 1.0 / (2 * igraph_ecount(&graph)), 0.01, 0, &membership, &nb_clusters, &quality);

    //printf("Leiden found %" IGRAPH_PRId " clusters using modularity, quality is %.4f.\n", nb_clusters, quality);
    //printf("Membership: ");
    //igraph_vector_print(&membership);
    //printf("\n");

    //igraph_vector_destroy(&degree);
    //igraph_vector_destroy(&membership);
    //igraph_destroy(&graph);
}


int selectFromCommunities(igraph_t* V, igraph_vector_t* membership, igraph_integer_t* nb_cluster, bool biggest,std::vector<int> candidate){
    std::vector<long int> commuSizes(*nb_cluster, 0);
    std::vector<long int> commuDelegate(*nb_cluster, -1);

    for (long int i = 0; i < igraph_vcount(V); i++){
	commuSizes[VECTOR(*membership)[i]] += 1;
	if (commuDelegate[VECTOR(*membership)[i]] == -1){
	    if (candidate[i] == 1)
	        commuDelegate[VECTOR(*membership)[i]] = i;
	}
    }
    if (biggest){
	for (long int i = 0; i < *nb_cluster; i++){
	   if (commuDelegate[i] == -1){
               commuSizes[i] = std::numeric_limits<int>::min();
	   }
	}	
	long int maxCommuIndex = max_element(commuSizes.begin(), commuSizes.end()) - commuSizes.begin();
	return commuDelegate[maxCommuIndex];
    }
    else{
	for (long int i = 0; i < *nb_cluster; i++){
	    if (commuDelegate[i] == -1)
		 commuSizes[i] = std::numeric_limits<int>::max();
	}
	long int minCommuIndex = min_element(commuSizes.begin(), commuSizes.end()) - commuSizes.begin();
	return commuDelegate[minCommuIndex];
    }    
}

int selectFromCommunities_cycle(igraph_t* V, igraph_vector_t* membership, igraph_integer_t* nb_cluster, bool biggest,std::vector<int> candidate, int* commuIndex){
    std::vector<long int> commuSizes(*nb_cluster, 0);
    std::vector<long int> commuDelegate(*nb_cluster, -1);

    for (long int i = 0; i < igraph_vcount(V); i++){
	if (candidate[i] == 1){
	    commuSizes[VECTOR(*membership)[i]] += 1;
	    if (commuDelegate[VECTOR(*membership)[i]] == -1)
	        commuDelegate[VECTOR(*membership)[i]] = i;
	}
    }
    auto remainingClusters = *nb_cluster - std::count(commuDelegate.begin(), commuDelegate.end(), -1); 
    auto remainingClustersSize = std::accumulate(commuSizes.begin(), commuSizes.end(), 0);
    if (((float)remainingClustersSize /(remainingClusters * igraph_vcount(V))) > 0.1)
	    return -1;
    while (commuDelegate[*commuIndex] == -1){
        *commuIndex += 1;
	if (*commuIndex >= *nb_cluster)
	    *commuIndex = 0;
    }
    return commuDelegate[*commuIndex];
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
    if(v_id == -1) { // only in the first round, select node with highest degree
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
        res_id = selectFromCommunitiesCycleOnce(delegate, commuIndex, nb_cluster);	
//	*commuIndex += 1;
//	if (*commuIndex >= *nb_cluster)
//	    *commuIndex = 0;
	if (res_id == -1){
            if(high) // select node with highest upper bound
                return maxuppernode;
            else // select node with lowest lower bound
                return minlowernode;
	}
	return res_id;
    }
    /*else if(high) // select node with highest upper bound
        return maxuppernode;
    else // select node with lowest lower bound
        return minlowernode;*/
}

std::vector<int> bdgecc(igraph_t V)
{
    // Picking the right component
    igraph_vector_ptr_t complist;
    igraph_vector_ptr_init(&complist, 0);
    igraph_decompose(&V, &complist, IGRAPH_WEAK, -1, 2);
    long int max_index = 0;
    igraph_integer_t max_vertice = 0;

    //std::cout << max_index;
    //

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

    std::vector<int> eccs(igraph_vcount(&V), 0); // FIXME may be useless
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

    igraph_vector_t membership; // TODO Free
    igraph_vector_init(&membership, igraph_vcount(&V));
    igraph_integer_t nb_cluster; //TODO Free 
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

	//TODO free

	
    }
    //TODO Prune here 
    while(candidates > 0)
    {
	high = !high;
	biggest = !biggest;
        v_id = selectFrom(V, v_id, maxuppernode, minlowernode, candidate, high, biggest, &membership, &nb_cluster, &commuIndex, delegate);
	maxuppernode = -1;
	minlowernode = -1;
	
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
    igraph_destroy((igraph_t*)VECTOR(*(&complist))[max_index]);
    igraph_free((igraph_t*)VECTOR(*(&complist))[max_index]);
    //igraph_destroy(&V);
    //igraph_free(&V);
    igraph_vector_ptr_destroy(&complist);
    igraph_vector_destroy(&membership);

    //igraph_vector_ptr_destroy(&complist);
    std::cout << '\n' << "Num bfs: " << number_of_bfs << '\n';
    return eccs;
}
