
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
			      char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
			      boolean * is_cycle, dBGraph * db_graph);


boolean db_graph_db_node_has_precisely_n_edges_with_status(dBNode * node,Orientation orientation,NodeStatus status,int n,
							   dBNode * * next_node, Orientation * next_orientation, Nucleotide * next_base,
							   dBGraph * db_graph);



boolean db_graph_detect_bubble(dBNode * node,
			       Orientation orientation,
			       int limit, int delta,
			       void (*node_action)(dBNode * node), 
			       int * length1, Nucleotide * base1, dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,char * seq1,
			       int * length2, Nucleotide * base2, dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,char * seq2,
			       dBGraph * db_graph);

boolean db_graph_db_node_smooth_bubble(dBNode * node, Orientation orientation, int limit,int delta,int coverage_limit,double ratio_threshold, 
				       void (*node_action)(dBNode * node),
				       dBGraph * db_graph);

boolean db_graph_db_node_prune_low_coverage(dBNode * node, int coverage,
					    void (*node_action)(dBNode * node),
					    dBGraph * db_graph);


// limit is the max length
// min_coverage, max_coverage and avg_coveragte refer to the internal nodes
int db_graph_supernode(dBNode * node,int limit,
		       void (*node_action)(dBNode * node), 
		       dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
		       char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,boolean * is_cycle,
		       dBGraph * db_graph);

int db_graph_db_node_clip_tip(dBNode * node, int limit,
			      void (*node_action)(dBNode * node),
			      dBGraph * db_graph);



boolean db_graph_is_condition_true_for_all_nodes_in_supernode(dBNode * node,int limit, 						
							      boolean (*condition_for_all_nodes)(dBNode * node),  
							      void (*node_action)(dBNode * node),
							      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, int* path_length,
							      char * string, double * avg_coverage,int * min,int * max, boolean * is_cycle, 
							      dBGraph * db_graph);

boolean db_graph_is_condition_true_for_at_least_one_node_in_supernode(dBNode * node,int limit, 
								      boolean (*condition_for_all_nodes)(dBNode * node),  
								      void (*node_action)(dBNode * node),
								      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,  int* path_length,
								      char * string, double * avg_coverage,int * min,int * max, boolean * is_cycle, 
								      dBGraph * db_graph);


boolean db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(dBNode * node,int limit,
										    boolean (*condition_for_all_nodes)(dBNode * node),  
										    void (*node_action)(dBNode * node), 
										    int min_start, int min_end, int min_diff,
										    dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,int * path_length,
										    char * string, double * avg_coverage,int * min,int * max, boolean * is_cycle,
										    dBGraph * db_graph);


void db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required, FILE* fout,
										  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index);

void db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required, FILE* fout,
											  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index);

void db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required,
                                                                                                       int min_start, int min_end, int min_diff, FILE* fout,
                                                                                                       boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index);


void db_graph_print_supernodes(char * filename,int max_length, dBGraph * db_graph);

void db_graph_print_coverage(dBGraph * db_graph);

void db_graph_clip_tips(dBGraph * db_graph);

void db_graph_remove_low_coverage_nodes(int coverage, dBGraph * db_graph);
void db_graph_dump_binary(char * filename, boolean (*condition)(dBNode * node), dBGraph * db_graph);
void db_graph_smooth_bubbles(int coverage,int limit, int delta, dBGraph * db_graph);


void db_graph_detect_vars(int delta, int max_length, dBGraph * db_graph);

void db_graph_traverse_with_array(void (*f)(HashTable*, Element *, int**, int),HashTable * hash_table, int** array, int length_of_array);
int db_graph_get_N50_of_supernodes(dBGraph* db_graph);


int int_cmp(const void *a, const void *b);


dBNode* db_graph_get_first_node_in_supernode_containing_given_node(dBNode* node,  dBGraph* db_graph);

dBNode* db_graph_get_next_node_in_supernode(dBNode* node, Orientation orientation, Orientation* next_orientation,  dBGraph* db_graph);

void db_graph_get_supernode_length_marking_it_as_visited(dBGraph* db_graph, Element* node, int** array_of_supernode_lengths, int length_of_array);


#endif /* DB_GRAPH_H_ */
