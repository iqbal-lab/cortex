#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>


typedef enum
 {
    full_flank_format   = 0,
    glf  = 1,
 } Variant_File_Format ;



void get_node_multiplicities(dBNode** branch1, int len_branch1, dBNode** branch2, int len_branch2,
                             int** br1_multiplicity_in_br1, int** br2_multiplicity_in_br2,
                             int** br1_multiplicity_in_br2, int** br2_multiplicity_in_br1,
			     boolean only_count_nodes_with_edge_in_specified_colour_func,
                             Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*) );


