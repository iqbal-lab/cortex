/*
  dB_graph.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <hash_table.h>
#include <stdio.h>

typedef HashTable dBGraph;

//functions applies to Node
void db_graph_print_supernode(FILE * file, dBNode * node, dBGraph * db_graph);

dBNode * db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, 
				Orientation * next_orientation,
				Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph);

//char * get_seq_from_elem_to_end_of_supernode(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle);

//Functions applies to whole graph
//void db_graph_clip_tip(dBNode * node, int limit,dBGraph * db_graph);
void db_graph_set_all_visited_nodes_to_status_none(dBGraph* hash_table);


#endif /* DB_GRAPH_H_ */
