/*
  hash_table.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <hash_table.h>
#include <stdio.h>

typedef HashTable dBGraph;

//print the supernode where the element is placed
void db_graph_print_supernode(FILE * file, dBNode * node, dBGraph * db_graph);

void db_graph_clip_tip(dBNode * node, int limit,dBGraph * db_graph);

#endif /* DB_GRAPH_H_ */
