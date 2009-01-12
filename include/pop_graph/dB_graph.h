/*
  dB_graph.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <hash_table.h>
#include <stdio.h>

typedef HashTable dBGraph;

//print the supernode where the element is placed
void db_graph_print_supernode(FILE * file, dBNode * node, dBGraph * db_graph);
void db_graph_print_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index);

//void db_graph_clip_tip(dBNode * node, int limit,dBGraph * db_graph);

char * get_seq_from_elem_to_end_of_supernode(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle);
char * get_seq_from_elem_to_end_of_supernode_for_specific_person_or_pop(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle, EdgeArrayType type, int index);
#endif /* DB_GRAPH_H_ */
