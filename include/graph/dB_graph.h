
/*
  db_graph.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <open_hash/hash_table.h>
#include <stdio.h>

typedef HashTable dBGraph;



int db_graph_clip_tip(dBGraph * db_graph);


int db_graph_get_perfect_path(dBNode * node, Orientation orientation, int limit, void (*node_action)(dBNode * node),
			      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
			      boolean * is_cycle, dBGraph * db_graph);

boolean db_graph_detect_perfect_bubble(dBNode * node,
				       Orientation * orientation, 
				       boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node), 
				       Nucleotide * base1, Nucleotide * base2, 
				       Nucleotide * labels,dBNode * * end_node,Orientation * end_orientation,
				       dBGraph * db_graph);

boolean db_graph_db_node_has_precisely_n_edges_with_status(dBNode * node,Orientation orientation,NodeStatus status,int n,
							   dBNode * * next_node, Orientation * next_orientation, Nucleotide * next_base,
							   dBGraph * db_graph);


int db_graph_supernode(dBNode * node,int limit,
		       boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node), 
		       char * string, dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
		       dBGraph * db_graph);


int db_graph_db_node_clip_tip(dBNode * node, int limit,
			     boolean (*condition)(dBNode * node),  void (*node_action)(dBNode * node),
			      dBGraph * db_graph);


boolean db_graph_is_condition_true_for_all_nodes_in_supernode(dBNode * node,int limit, boolean (*condition_on_initial_node_only)(dBNode * node), boolean (*condition_for_all_nodes)(dBNode * node),
                                                              void (*node_action)(dBNode * node),
                                                              char * string,dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, int* path_length,
                                                              dBGraph * db_graph);

boolean db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(dBNode * node,int limit, boolean (*condition_on_initial_node_only)(dBNode * node), boolean (*condition_for_all_nodes)(dBNode * node),  
										    void (*node_action)(dBNode * node), int min_start, int min_end,
										    char * string,dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, int* path_length,
										    dBGraph * db_graph);


void db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required, 
                                                                                  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index);



void db_graph_print_supernodes(dBGraph * db_graph);

void db_graph_detect_snps(dBGraph * db_graph);

void db_graph_print_coverage(dBGraph * db_graph);

void db_graph_clip_tips(dBGraph * db_graph);


#endif /* DB_GRAPH_H_ */
