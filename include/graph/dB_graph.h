
/*
  hash_table.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <open_hash/hash_table.h>
#include <stdio.h>

typedef HashTable dBGraph;





//print the supernode where the element is placed
void db_graph_print_supernode(FILE * file, dBNode * node, dBGraph * db_graph);

int db_graph_clip_tip(dBNode * node, int limit,dBGraph * db_graph);

char * get_seq_from_elem_to_end_of_supernode(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle, char * seq, int max_length);
        
dBNode* db_graph_get_first_node_in_supernode_containing_given_node(dBNode* node,  dBGraph* db_graph);

dBNode* db_graph_get_next_node_in_supernode(dBNode* node, Orientation orientation, Orientation* next_orientation,  dBGraph* db_graph);

int db_graph_get_perfect_path(dBNode * node, Orientation orientation, dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, boolean * is_cycle, int limit, NodeStatus status,dBGraph * db_graph);
boolean db_graph_detect_perfect_bubble(dBNode * node,
				       Orientation * orientation,
				       Nucleotide * base1, Nucleotide * base2, 
				       Nucleotide * labels,
				       dBNode ** end_node, Orientation * end_orientation,
				       dBGraph * db_graph);

void db_graph_get_supernode_length_marking_it_as_visited(dBGraph* db_graph, Element* node, int** array_of_supernode_lengths, int length_of_array);


void db_graph_traverse_with_array(void (*f)(HashTable*, Element *, int**, int),HashTable * hash_table, int** array, int length_of_array);

int db_graph_get_N50_of_supernodes(dBGraph* db_graph);
int int_cmp(const void *a, const void *b);


#endif /* DB_GRAPH_H_ */
