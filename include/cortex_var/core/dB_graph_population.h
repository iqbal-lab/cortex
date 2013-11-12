/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  dB_graph_population.h

  wrapper for hash_table for people who want to think of the hash table as 
  many de bruijn graphs drawn on top of each other.

*/

#ifndef DB_GRAPH_POPULATION_H_
#define DB_GRAPH_POPULATION_H_

#include "global.h"
#include "element.h"
#include "open_hash/hash_table.h"
#include "seq.h"
#include "dB_graph.h"
#include "dB_graph_supernode.h"
#include "db_variants.h"
#include "graph_info.h"
#include "model_selection.h"
#include "db_complex_genotyping.h"


typedef struct {
  double mean;
}BasicStats;

dBNode * db_graph_get_next_node_for_specific_person_or_pop(dBNode * current_node, Orientation current_orientation,
                                                           Orientation * next_orientation,
                                                           Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, int index);

//This function does not  check that it there is such an edge in the specified person/colour - but it does check if the target node is in the specific person.
//if you want to be sure the dge exists in that colour, then check it before calling this function
//The last argument allows you to apply the operation to some subgraph - eg you might take the unuiion of colours 2 and 3, or of all colours.
//If for example you wanted to get the next node in the graph irrespective of colour, get_colour would return the union (bitwise AND) of all edges in a node.
dBNode * db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(dBNode * current_node, Orientation current_orientation, 
								       Orientation * next_orientation,
								       Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, 
								       Edges (*get_colour)(const dBNode*)
								       );


boolean db_graph_db_node_has_precisely_n_edges_with_status_for_specific_person_or_pop(
	dBNode * node,Orientation orientation,NodeStatus status,int n,
  dBNode * * next_node, Orientation * next_orientation, Nucleotide * next_base,
  dBGraph * db_graph, int index);


// clip a tip in the graph (the tip starts in node)
// limit is max length for tip
// node_action is applied to all the elements in the tip
// returns the length of the tip (0 means no length)
int db_graph_db_node_clip_tip_for_specific_person_or_pop(dBNode * node, int limit,
							 void (*node_action)(dBNode * node),
							 dBGraph * db_graph, int index);



int db_graph_db_node_clip_tip_in_subgraph_defined_by_func_of_colours(dBNode * node, int limit,
								     void (*node_action)(dBNode * node),
								     dBGraph * db_graph, 
								     Edges (*get_colour)(const dBNode*),
								     void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
								     void (*apply_reset_to_colour)(dBNode*)
								     );



void apply_reset_to_specific_edge_in_union_of_all_colours(dBNode* node, Orientation or, Nucleotide nuc);
void apply_reset_to_all_edges_in_union_of_all_colours(dBNode* node );


void db_graph_clip_tips_for_specific_person_or_pop(dBGraph * db_graph, int index);

void db_graph_clip_tips_in_subgraph_defined_by_func_of_colours(dBGraph * db_graph,
							       Edges (*get_colour)(const dBNode*),
							       void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
							       void (*apply_reset_to_colour)(dBNode*));

void db_graph_clip_tips_in_union_of_all_colours(dBGraph* db_graph);





// 1. the argument sum_of_covgs_in_desired_colours allows you to choose which colour (or union of colours) you want to apply this to
// eg you might want  to "remove" (as defined by the action) any node that has coverage <= your threshold in the UNION of all colours, or in colour red or whatever
// SO - this func returns the sum of coverages in the colours you care about
// 2. the argument get_edge_of_interest is a function that gets the "edge" you are interested in - may be a single edge/colour from the graph, or might be a union of some edges
// 3 Pass apply_reset_to_specified_edges which applies reset_one_edge to whichever set of edges you care about, 
// 4 Pass apply_reset_to_specified_edges_2 which applies db_node_reset_edges to whichever set of edges you care about, 
boolean db_graph_db_node_prune_low_coverage(dBNode * node, Covg coverage,
					    void (*node_action)(dBNode * node),
					    dBGraph * db_graph, 
					    Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
					    Edges (*get_edge_of_interest)(const Element*),
					    void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide),
					    void (*apply_reset_to_specified_edges_2)(dBNode*) );


boolean db_graph_db_node_prune_low_coverage_ignoring_colours(dBNode * node, Covg coverage,
							     void (*node_action)(dBNode * node),
							     dBGraph * db_graph);

boolean db_graph_db_node_prune_without_condition(dBNode * node, 
						 void (*node_action)(dBNode * node),
						 dBGraph * db_graph, 
						 Edges (*get_edge_of_interest)(const Element*),
						 void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide),
						 void (*apply_reset_to_specified_edges_2)(dBNode*) );



//if the node has covg <= coverage (arg2) and its supernode has length <=kmer+1 AND all the interiro nodes of the supernode have this low covg, then
//prune the whole of the interior of the supernode
//note the argument supernode_len is set to -1 if the passed in node has status != none
// returns true if the supernode is pruned, otherwise false
// if returns false, cannot assume path_nodes etc contain anything useful
boolean db_graph_remove_supernode_containing_this_node_if_looks_like_induced_by_error(
	dBNode* node, Covg coverage, dBGraph * db_graph, int max_expected_sup_len,
	Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
	Edges (*get_edge_of_interest)(const Element*), 
	void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
	void (*apply_reset_to_specified_edges_2)(dBNode*),
	dBNode** path_nodes, Orientation* path_orientations, Nucleotide* path_labels, 
	char* supernode_str, int* supernode_len);

// traverse graph. At each node, if covg <= arg1, get its supernode. If ALL interior nodes have covg <= arg1 
// then prune the node, and the interior nodes of the supernode.
// returns the number of pruned supernodes
long long db_graph_remove_errors_considering_covg_and_topology(
	Covg coverage, dBGraph * db_graph,
	Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
	Edges (*get_edge_of_interest)(const Element*), 
	void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
	void (*apply_reset_to_specified_edges_2)(dBNode*),
	int max_expected_sup);


boolean db_graph_remove_supernode_containing_this_node_if_more_likely_error_than_sampling(
  dBNode* node, int num_haploid_chroms, 
  double total_depth_of_covg, int read_len, 
  double error_rate_per_base,
  int max_length_allowed_to_prune,
  dBGraph* db_graph, int max_expected_sup,
  Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
  Edges (*get_edge_of_interest)(const Element*), 
  void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
  void (*apply_reset_to_specified_edges_2)(dBNode*),
  dBNode** path_nodes, Orientation* path_orientations, Nucleotide* path_labels,
  char* supernode_string, int* supernode_len, Covg covg_thresh);


long long db_graph_remove_supernodes_more_likely_errors_than_sampling(
  dBGraph * db_graph, GraphInfo* ginfo, GraphAndModelInfo* model_info,
  int max_length_allowed_to_prune, 
  Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
  Edges (*get_edge_of_interest)(const Element*), 
  void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
  void (*apply_reset_to_specified_edges_2)(dBNode*),
  int max_expected_sup, int manual_threshold);


int db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(
  dBNode *node, Orientation orientation, int limit, 
  Nucleotide fst_nucleotide,
  void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *seq, double *avg_coverage,Covg *min_coverage, Covg *max_coverage,
  boolean *is_cycle, dBGraph *db_graph, int index);

int db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(
  dBNode *node, Orientation orientation, int limit, 
  Nucleotide fst_nucleotide, void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *seq, double *avg_coverage, Covg *min_coverage, Covg *max_coverage,
  boolean *is_cycle, dBGraph *db_graph, 
  Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*));

int db_graph_get_perfect_path_for_specific_person_or_pop(
  dBNode *node, Orientation orientation, int limit, 
  void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *seq, double *avg_coverage, Covg *min_coverage, Covg *max_coverage,
  boolean *is_cycle, dBGraph *db_graph, int index);


int db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(
  dBNode *node, Orientation orientation, int limit, 
  void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *seq, double *avg_coverage, Covg *min_coverage, Covg *max_coverage,
  boolean *is_cycle, dBGraph *db_graph, 
  Edges (*get_colour)(const dBNode*),
  Covg (*get_covg)(const dBNode*)) ;


boolean db_graph_detect_bubble_for_specific_person_or_population(dBNode *node,
  Orientation orientation, int limit,
  void (*node_action)(dBNode *node), 
  int *length1, dBNode **path_nodes1, Orientation *path_orientations1, Nucleotide *path_labels1,
  char *seq1, double *avg_coverage1, Covg *min_coverage1, Covg *max_coverage1,
  int *length2,dBNode **path_nodes2, Orientation *path_orientations2, Nucleotide *path_labels2,
  char *seq2, double *avg_coverage2, Covg *min_coverage2, Covg *max_coverage2,
  dBGraph *db_graph, int index);


boolean db_graph_detect_bubble_in_subgraph_defined_by_func_of_colours(boolean (*start_condition)(dBNode*),								      
								      dBNode *node, Orientation orientation, int limit,
								      void (*node_action)(dBNode *node), 
								      int *length1,dBNode **path_nodes1, Orientation *path_orientations1, Nucleotide *path_labels1,
								      char *seq1, double *avg_coverage1, Covg *min_coverage1, Covg *max_coverage1,
								      int *length2,dBNode **path_nodes2, Orientation *path_orientations2, Nucleotide *path_labels2,
								      char *seq2, double *avg_coverage2, Covg *min_coverage2, Covg *max_coverage2,
								      dBGraph *db_graph, Edges (*get_colour)(const dBNode*)  , Covg (*get_covg)(const dBNode*),
								      boolean apply_special_action_to_branches, void (*special_action)(dBNode *node) );


boolean db_graph_db_node_smooth_bubble_for_specific_person_or_pop(
  dBNode *node, Orientation orientation, int limit, Covg coverage_limit,
  void (*node_action)(dBNode *node), dBGraph *db_graph, int index);


int db_graph_supernode_for_specific_person_or_pop(
  dBNode *node,int limit,void (*node_action)(dBNode *node), 
  dBNode **path_nodes, Orientation *path_orientations,
  Nucleotide *path_labels, char *supernode_str, 
  double *avg_coverage, Covg *min, Covg *max, boolean *is_cycle, 
  dBGraph *db_graph, int index);


int db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(
  dBNode *node,int limit,void (*node_action)(dBNode *node), 
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *supernode_str, double *avg_coverage, Covg *min_covg, Covg *max_covg,
  boolean *is_cycle, int *query_node_posn, dBGraph *db_graph, int index);
  
int db_graph_supernode_returning_query_node_posn_in_subgraph_defined_by_func_of_colours(
  dBNode *node,int limit,void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  char *supernode_str, double *avg_coverage, Covg *min_covg, Covg *max_covg,
  boolean *is_cycle, int *query_node_posn, dBGraph *db_graph,
  Edges (*get_colour)(const dBNode*),
  Covg (*get_covg)(const dBNode*));



int db_graph_supernode_in_subgraph_defined_by_func_of_colours(
	dBNode *node,int limit,void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations,
  Nucleotide *path_labels, char *supernode_str,
  double *avg_coverage, Covg *min_covg, Covg *max_covg,
  boolean *is_cycle, dBGraph *db_graph, 
  Edges (*get_colour)(const dBNode*),
  Covg (*get_covg)(const dBNode*));


boolean db_graph_is_condition_true_for_all_nodes_in_supernode(
  dBNode *node,int limit,
  boolean (*condition_for_all_nodes)(dBNode *node),
  void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels, int*path_length,
  char *string, double *avg_coverage, Covg *min_covg, Covg *max_covg,
  boolean *is_cycle, dBGraph *db_graph, int index);

boolean db_graph_is_condition_true_for_at_least_one_node_in_supernode(
  dBNode *node, int limit,
  boolean (*condition_for_all_nodes)(dBNode *node),
  void (*node_action)(dBNode *node),
  dBNode **path_nodes, Orientation *path_orientations,
  Nucleotide *path_labels,  int*path_length,
  char *string, double *avg_coverage, Covg *min_covg, Covg *max_covg,
  boolean *is_cycle, dBGraph *db_graph, int index);

boolean db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(
  dBNode *node,int limit, boolean (*condition_for_all_nodes)(dBNode *node),  
  void (*node_action)(dBNode *node), int min_start, int min_end, int min_diff,
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide *path_labels,
  int *path_length, char *string,
  double *avg_coverage, Covg *min_covg, Covg *max_covg,
  boolean *is_cycle, dBGraph *db_graph, int index);

void db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(
  dBGraph *db_graph, boolean (*condition)(dBNode * node), Covg min_covg_required, FILE* fout,  
  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index, int index);


void db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(
  dBGraph *db_graph, boolean (*condition)(dBNode *node), Covg min_covg_required,
  FILE *fout, boolean is_for_testing, char** for_test_array_of_supernodes,
  int *for_test_index, int index);

// can use this toprint potential indels.
// min_start and min_end are the number of nodes of overlap with the reference that  you want at the start and  end of the supernode
void db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(
  dBGraph *db_graph, boolean (*condition)(dBNode *node), Covg min_covg_required,
  int min_start, int min_end, int min_diff, FILE *fout,
  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index,
  int index);

//routine to DETECT/DISCOVER variants directly from the graph - reference-free (unless you have put the reference in the graph!)
// last argument is a condition which you apply to the flanks and branches to decide whether to call.
// e.g. this might be some constraint on the coverage of the branches, or one might have a condition that one branch
//      was in one colour  and the other in a different colour, or maybe that both branches are in the same colour
void db_graph_detect_vars(FILE* fout, /* FILE* fout_gls, */ int max_length, dBGraph * db_graph, 
			  boolean (*condition)(VariantBranchesAndFlanks*), 
			  void (*action_branches)(dBNode*),
			  void (*action_flanks)(dBNode*),
			  Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*),
			  void (*print_extra_info)(VariantBranchesAndFlanks*, FILE*),
			  boolean apply_model_selection, 
			  boolean (*model_selection_condition)(AnnotatedPutativeVariant*),
			  GraphAndModelInfo* model_info,
			  boolean (*start_condition)(dBNode*) );

void db_graph_detect_vars_after_marking_vars_in_reference_to_be_ignored(
									FILE* fout, /* FILE* fout_gls, */ 
									int max_length, dBGraph *db_graph, 
									boolean (*condition)(VariantBranchesAndFlanks*),
									Edges (*get_colour_ref)(const dBNode*),
									Covg (*get_covg_ref)(const dBNode*),
									Edges (*get_colour_indiv)(const dBNode*),
									Covg (*get_covg_indiv)(const dBNode*),
									void (*print_extra_info)(VariantBranchesAndFlanks*, FILE*));

boolean detect_vars_condition_always_true(VariantBranchesAndFlanks*);
boolean detect_vars_condition_branches_not_marked_to_be_ignored(VariantBranchesAndFlanks* var);
boolean detect_vars_condition_always_false(VariantBranchesAndFlanks*);
boolean detect_vars_condition_flanks_at_least_3(VariantBranchesAndFlanks*);

boolean detect_vars_condition_is_hom_nonref_given_colour_funcs_for_ref_and_indiv(
	VariantBranchesAndFlanks* var, 
  Covg (*get_covg_colourfunc1)(const dBNode*),
  Covg (*get_covg_colourfunc2)(const dBNode*));

boolean detect_vars_condition_is_hom_nonref_with_min_covg_given_colour_funcs_for_ref_and_indiv(
	VariantBranchesAndFlanks* var, 
  Covg (*get_covg_ref)(const dBNode*),
  Covg (*get_covg_indiv)(const dBNode*),
  Covg min_covg);




//given two lists of colours, we want to call variants where one branch is entirely in the 
// union of the colours of the first list, but not in the second, and vice-versa for the second branch
//UNLESS both lists are identical, in which case we just call bubbles in the union
//A consequence is that if you call bubbles on 1,2/1,3 you'll get nothing.
//if exclude_ref_bubbles_first==false, just pass in NULL, NULL for the last two arguments
void db_graph_detect_vars_given_lists_of_colours(FILE* fout, /* FILE* fout_gls, */ int max_length, dBGraph * db_graph, 
						 int* first_list, int len_first_list,
						 int* second_list,  int len_second_list,
						 boolean (*extra_condition)(VariantBranchesAndFlanks* var),
						 void (*print_extra_info)(VariantBranchesAndFlanks*, FILE*),
						 boolean exclude_ref_bubbles_first, 
						 Edges (*get_colour_ref)(const dBNode*),
						 Covg (*get_covg_ref)(const dBNode*),
						 boolean apply_model_selection, 
						 boolean (*model_selection_condition)(AnnotatedPutativeVariant*, GraphAndModelInfo*),
						 GraphAndModelInfo* model_info);





void db_graph_print_supernodes_for_specific_person_or_pop(
	char *filename_sups, char *filename_sings, int max_length, dBGraph *db_graph, int index,
  void (*print_extra_info)(dBNode**, Orientation*, int, FILE*));


long long  db_graph_count_error_supernodes(int max_length, dBGraph *db_graph, int index, 
  dBNode **path_nodes, Orientation *path_orientations, Nucleotide* path_labels, char *seq,
  boolean (*condition_node_not_visited)(dBNode *n), 
  boolean (*condition_is_error_supernode)(dBNode **path, int len, 
                                          int *count_error_nodes,
                                          Covg *count_error_covg),
  void (*action_set_visited)(dBNode* e));


void db_graph_print_supernodes_defined_by_func_of_colours(char * filename_sups, char* filename_sings, int max_length, 
							  dBGraph * db_graph, Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*),
							  void (*print_extra_info)(dBNode**, Orientation*, int, FILE*));

void db_graph_print_supernodes_defined_by_func_of_colours_given_condition(
  char * filename_sups, char* filename_sings, int max_length, 
  dBGraph * db_graph, Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*),
  void (*print_extra_info)(dBNode**, Orientation*, int, FILE*),
  boolean (*condition)(dBNode** path, Orientation* ors, int len));

void db_graph_print_novel_supernodes(char* outfile, int max_length, dBGraph * db_graph, 
				     int* first_list, int len_first_list,
				     int* second_list,  int len_second_list,
				     int min_contig_length_bp, int min_percentage_novel,
				     void (*print_extra_info)(dBNode**, Orientation*, int, FILE*));

void db_graph_print_coverage_for_specific_person_or_pop(dBGraph * db_graph, int index);

void print_covg_stats_for_timestamps_for_supernodes(char* outfile, dBGraph * db_graph, int max_expected_sup_len);



// 1. sum_of_covgs_in_desired_colours returns the sum of the coverages for the colours you are interested in
// 1. the argument get_edge_of_interest is a function that gets the "edge" you are interested in - may be a single edge/colour from the graph, or might be a union of some edges 
// 2 Pass apply_reset_to_specified_edges which applies reset_one_edge to whichever set of edges you care about,
// 3 Pass apply_reset_to_specified_edges_2 which applies db_node_reset_edges to whichever set of edges you care about,
void db_graph_remove_low_coverage_nodes(Covg coverage, dBGraph * db_graph,
					Covg (*sum_of_covgs_in_desired_colours)(const Element *), 
					Edges (*get_edge_of_interest)(const Element*), 
					void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
					void (*apply_reset_to_specified_edges_2)(dBNode*) );

void db_graph_remove_low_coverage_nodes_ignoring_colours(Covg coverage, dBGraph *db_graph);

int db_graph_dump_binary(char *filename, boolean (*condition)(dBNode *node),
                         dBGraph *db_graph, GraphInfo *db_graph_info, int version);

void db_graph_dump_single_colour_binary_of_colour0(char *filename,
                                                   boolean (*condition)(dBNode *node),
                                                   dBGraph *db_graph, 
						                                       GraphInfo *db_graph_info,
                                                   int version);

void db_graph_dump_single_colour_binary_of_specified_colour(char * filename, boolean (*condition)(dBNode * node), 
							    dBGraph * db_graph, GraphInfo* db_graph_info, int colour,
							    int version);

boolean db_node_is_supernode_end(dBNode * element,Orientation orientation, int edge_index, dBGraph* db_graph);




//***
// functions that are direct extensions of those in hash_table.h
// ****

//corresponds to hash_table_find
dBNode *  db_graph_find_node_restricted_to_specific_person_or_population(Key key, dBGraph * hash_table, int index);


// - not implemented: dBNode *  db_graph_find_or_insert_node_restricted_to_specific_person_or_population(Key key, boolean * found, dBGraph * hash_table, int index);



void db_graph_traverse_with_array(void (*f)(HashTable*, Element *, int**, int, int),
                                  HashTable * hash_table, int** array, int length_of_array, int index);

void db_graph_traverse_with_array_of_uint64(void (*f)(HashTable*, Element *, uint64_t*, uint32_t, int),
					    HashTable * hash_table, uint64_t* array, uint32_t length_of_array, int index);
 
void db_graph_get_covg_distribution(char* filename, dBGraph* db_graph, 
				    int index,  boolean (*condition)(dBNode* elem) );

long long  db_graph_count_covg1_kmers_in_func_of_colours(dBGraph* db_graph, Covg (*get_covg)(const dBNode*) );
long long  db_graph_count_covg2_kmers_in_func_of_colours(dBGraph* db_graph, Covg (*get_covg)(const dBNode*) );

void db_graph_get_supernode_length_marking_it_as_visited(
  dBGraph* db_graph, Element* node, int** array_of_supernode_lengths,
  int length_of_array, int index);

int db_graph_get_N50_of_supernodes(dBGraph* db_graph, int index);


void db_graph_traverse_to_gather_statistics_about_people(
  void (*f)(HashTable*, Element *, int**, int),
  HashTable *hash_table, int** array, int num_people);



//TODO - get rid of next two
dBNode* db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(
  dBNode* node, int index, dBGraph* db_graph);

dBNode* db_graph_get_next_node_in_supernode_for_specific_person_or_pop(
  dBNode* node, Orientation orientation, Orientation* next_orientation,
  int index, dBGraph* db_graph);


// min_covg_for_pop_supernode means number of people!  Not coverage of nodes
void  db_graph_find_population_consensus_supernode_based_on_given_node(
  Sequence* pop_consensus_supernode, int max_length_of_supernode, dBNode* node, 
  int min_covg_for_pop_supernode, int min_length_for_pop_supernode, dBGraph* db_graph);

int db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(
  char* subsection, dBNode* node, int start, int end, int index, dBGraph* db_graph);

void  db_graph_get_best_sub_supernode_given_min_covg_and_length_for_specific_person_or_pop(
  dBNode* first_node_in_supernode, int* index_for_start_of_sub_supernode,
  int* length_of_best_sub_supernode, int min_people_coverage,
  int min_length_of_sub_supernode, int index, dBGraph* db_graph);

void print_node_to_file_according_to_how_many_people_share_it(
  HashTable* db_graph, dBNode * node, FILE** list_of_file_ptrs);

void find_out_how_many_individuals_share_this_node_and_add_to_statistics(
  HashTable* db_graph, dBNode * node, int** array_of_counts, int number_of_people);

int db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta
                                                  (FILE* chrom_fptr, 
						   int number_of_nodes_to_load, int number_of_nodes_loaded_last_time,
						   int length_of_arrays, 
						   dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
						   Sequence* seq, KmerSlidingWindow* kmer_window, 
						   boolean expecting_new_fasta_entry, boolean last_time_was_not_start_of_entry, 
						   dBGraph* db_graph ) ;

int db_node_addr_cmp(const void *a, const void *b);

void get_coverage_from_array_of_nodes(
  dBNode** array, int length, Covg* min_coverage, Covg* max_coverage,
  double* avg_coverage, Covg* mode_coverage,
  double* percent_nodes_having_modal_value, int index);

void get_percent_novel_from_array_of_nodes(dBNode** array, int length,
                                           double* percent_novel,
                                           int index_of_reference_in_array_of_edges);

void get_covg_of_nodes_in_one_but_not_other_of_two_arrays(
  dBNode** array1, dBNode** array2, int length1, int length2,
  int* num_nodes_in_array_1not2, int * num_nodes_in_array_2not1,
  Covg** covgs_in_1not2, Covg** covgs_in_2not1,
  dBNode** reused_working_array1, dBNode** reused_working_array2, int index);


boolean make_reference_path_based_sv_calls_condition_always_true(VariantBranchesAndFlanks* var, int colour_ref, int colour_indiv);
boolean make_reference_path_based_sv_calls_condition_always_true_for_func_of_colours(VariantBranchesAndFlanks* var,  int colour_of_ref,
										     Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*) );


void print_no_extra_info(VariantBranchesAndFlanks* var, FILE* fout);
void print_no_extra_supernode_info(dBNode** node_array, Orientation* or_array, int len, FILE* fout);


int db_graph_make_reference_path_based_sv_calls(FILE* chrom_fasta_fptr,
  int index_for_indiv_in_edge_array, int index_for_ref_in_edge_array,
  int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, int max_anchor_span,
  Covg min_covg, Covg max_covg, int max_expected_size_of_supernode,
  int length_of_arrays, dBGraph* db_graph, FILE* output_file,
  int max_desired_returns,
  char** return_flank5p_array, char** return_trusted_branch_array, char** return_variant_branch_array, 
  char** return_flank3p_array, int** return_variant_start_coord,
  boolean (*condition)(VariantBranchesAndFlanks* var,  int colour_of_ref,  int colour_of_indiv),
  void (*action_for_branches_of_called_variants)(VariantBranchesAndFlanks* var),
  void (*print_extra_info)(VariantBranchesAndFlanks* var, FILE* fout), GraphAndModelInfo* model_info,
  AssumptionsOnGraphCleaning assump);

boolean make_reference_path_based_sv_calls_condition_always_true_in_subgraph_defined_by_func_of_colours(
  VariantBranchesAndFlanks* var, int colour_of_ref,
  Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*));

int db_graph_make_reference_path_based_sv_calls_in_subgraph_defined_by_func_of_colours(
  FILE* chrom_fasta_fptr,
  Edges (*get_colour)(const dBNode*),
  Covg (*get_covg)(const dBNode*),
  int ref_colour, //int index_for_ref_in_edge_array,
  int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, 
  int max_anchor_span, Covg min_covg, Covg max_covg, 
  int max_expected_size_of_supernode, int length_of_arrays, dBGraph* db_graph,
  FILE* output_file, int max_desired_returns,
  char** return_flank5p_array, char** return_trusted_branch_array, char** return_variant_branch_array, 
  char** return_flank3p_array, int** return_variant_start_coord,
  boolean (*condition)(VariantBranchesAndFlanks* var,  int colour_of_ref,  
  Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*)),
  void (*action_for_branches_of_called_variants)(VariantBranchesAndFlanks* var),
  void (*print_extra_info)(VariantBranchesAndFlanks* var, FILE* fout),
  GraphAndModelInfo* model_info, AssumptionsOnGraphCleaning assump,int start_variant_numbering_with_this
  );


int db_graph_make_reference_path_based_sv_calls_given_list_of_colours_for_indiv(
  int* list, int len_list,
  FILE* chrom_fasta_fptr, int ref_colour,
  int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, 
  int max_anchor_span, Covg min_covg, Covg max_covg, 
  int max_expected_size_of_supernode, int length_of_arrays, dBGraph* db_graph,
  FILE* output_file, int max_desired_returns,
  char** return_flank5p_array, char** return_trusted_branch_array,
  char** return_variant_branch_array, 
  char** return_flank3p_array, int** return_variant_start_coord,
  boolean (*condition)(VariantBranchesAndFlanks* var,  int colour_of_ref,  
  Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*)),
  void (*action_for_branches_of_called_variants)(VariantBranchesAndFlanks* var),
  void (*print_extra_info)(VariantBranchesAndFlanks* var, FILE* fout),
  GraphAndModelInfo* model_info, 
  AssumptionsOnGraphCleaning assump, int start_numbering_vars_from_this_number
  );

boolean condition_always_true(dBNode** flank_5p, int len5p, dBNode** ref_branch,
                              int len_ref, dBNode** var_branch, int len_var,
			                        dBNode** flank_3p, int len3p, int colour_of_ref,
                              int colour_of_indiv);


void apply_to_all_nodes_in_path_defined_by_fasta(
  void (*func)(dBNode*), FILE* fasta_fptr, int chunk_size, dBGraph* db_graph);

// use text_describing_comparison_with_other_path to allow printing coverages of
// nodes in this path but not in some specific other path
void print_fasta_from_path_for_specific_person_or_pop(FILE *fout,
  char * name, int length, double avg_coverage,
  Covg min_coverage, Covg max_coverage, Covg modal_coverage,
  double percent_nodes_with_modal_coverage, double percent_novel,
  dBNode * fst_node, Orientation fst_orientation,
  dBNode * lst_node, Orientation lst_orientation,
  char* text_describing_comparison_with_other_path, // may be NULL
  Covg* coverages_nodes_in_this_path_but_not_some_other, // may be NULL
  int length_of_coverage_array,
  char * string, // labels of paths
  int kmer_size, boolean include_first_kmer, int index);


// use text_describing_comparison_with_other_path to allow printing coverages
// of nodes in this path but not in some specific other path
void print_fasta_with_all_coverages_from_path_for_specific_person_or_pop(FILE *fout,
  char * name, int length, double avg_coverage,
  Covg min_coverage, Covg max_coverage, Covg modal_coverage,
  double percent_nodes_with_modal_coverage, double percent_novel,
  dBNode * fst_node, Orientation fst_orientation,
  dBNode * lst_node, Orientation lst_orientation,
  char* text_describing_comparison_with_other_path, // may be NULL
  Covg* coverages_nodes_in_this_path_but_not_some_other, // may be NULL
  int length_of_coverage_array,//refers to prev argument
  Covg* coverages_nodes_in_this_path, 
  Covg* coverages_in_ref_nodes_in_this_path,
  int number_nodes_in_this_path, char * string, //labels of paths
  int kmer_size, boolean include_first_kmer, int index);



boolean does_this_path_exist_in_this_colour(dBNode** array_nodes,
                                            Orientation* array_orientations,
                                            int len, Edges (*get_colour)(const dBNode*),
                                            dBGraph* db_graph);





void print_standard_extra_supernode_info(dBNode** node_array,
                                         Orientation* or_array, int len, FILE* fout);

void print_median_extra_supernode_info(dBNode** node_array,
				       Orientation* or_array,
				       int len, CovgArray* working_ca, FILE* fout);

void print_standard_extra_info(VariantBranchesAndFlanks* var, FILE* fout);

void print_median_covg_extra_info(VariantBranchesAndFlanks* var, CovgArray* working_ca,FILE* fout);

long long db_graph_health_check(boolean fix, dBGraph * db_graph);
long long db_graph_clean_orphan_edges(dBGraph * db_graph);

void db_graph_wipe_colour(int colour, dBGraph* db_graph);
void db_graph_wipe_two_colours_in_one_traversal(int colour1, int colour2, dBGraph* db_graph);

void db_graph_print_colour_overlap_matrix(int* first_col_list, int num1,
                                          int* second_col_list, int num2,
					  dBGraph* db_graph);

void print_call_given_var_and_modelinfo(VariantBranchesAndFlanks* var,
                                        FILE* fout,
                                        GraphAndModelInfo* model_info,
					DiscoveryMethod which_caller,
                                        dBGraph* db_graph,
                                        void (*print_extra_info)(VariantBranchesAndFlanks*, FILE*),
					AssumptionsOnGraphCleaning assump,
                                        GenotypingWorkingPackage* gwp,
                                        LittleHashTable* little_dbg, CovgArray* working_ca);


void db_graph_get_stats_of_supernodes_that_split_two_colour(int max_length, int colour1, int colour2,
							    dBGraph * db_graph, Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*),
							    boolean (*condition)(dBNode**, int, int*), int* bins_1_not_2, int* bins_2_not_1);

void db_graph_get_covgs_in_all_colours_of_col0union1_sups(
  int max_length, dBGraph * db_graph,
  Covg* pgf, Covg* cox, 
  Covg* data1, Covg* data2, Covg* data3, Covg* data4, Covg* data5);

void db_graph_get_proportion_of_cvg_on_each_sup(int max_length, dBGraph * db_graph,
						int tot_pgf, int tot_cox, 
						int tot_data1, int tot_data2, int tot_data3, int tot_data4, int tot_data5,
						boolean (*condition)(dBNode**, int, int*), FILE* fout);

void call_bubbles_distinguishing_cox_pgf(dBGraph* db_graph, int max_length,
                                         FILE* fout, GraphAndModelInfo* model_info);

#endif /* DB_GRAPH_POPULATION_H_ */
