/*
  dB_graph.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <hash_table.h>
#include <stdio.h>

typedef HashTable dBGraph;

//pays no attention to whether there is an edge joining current_node to the node you would get by adding this nucleotide.
//just checksto see if such a node is in the graph
dBNode * db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, 
				Orientation * next_orientation,
				Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph);




//Functions applying to whole graph

int db_graph_clip_tip(dBNode * node, int limit,dBGraph * db_graph);

void db_graph_set_all_visited_nodes_to_status_none(dBGraph* hash_table);


#endif /* DB_GRAPH_H_ */
