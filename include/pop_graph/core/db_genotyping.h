#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>

void get_node_multiplicities(dBNode** branch1, int len_branch1, dBNode** branch2, int len_branch2,
                             int** br1_multiplicity_in_br1, int** br2_multiplicity_in_br2,
                             int** br1_multiplicity_in_br2, int** br2_multiplicity_in_br1,
                             Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*) );


