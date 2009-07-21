/*
  dB_graph_population.h

  wrapper for hash_table for people who want to think of the hash table as 
  many de bruijn graphs drawn on top of each other.

*/

#ifndef DB_GRAPH__POPULATION_H_
#define DB_GRAPH_POPULATION_H_

#include <element.h>
#include <open_hash/hash_table.h>
#include <seq.h>
#include <dB_graph.h>
#include <dB_graph_supernode.h>

//***
// functions that are direct extensions of those in hash_table.h
// ****

//corresponds to hash_table_find
dBNode *  db_graph_find_node_restricted_to_specific_person_or_population(Key key, dBGraph * hash_table, EdgeArrayType type, int index);

//corresponds to hash_table_find_or_insert
dBNode *  db_graph_find_or_insert_node_restricted_to_specific_person_or_population(Key key, boolean * found, dBGraph * hash_table, EdgeArrayType type, int index);


//**************
// functions pulled over from Mario's graph code
// ************

//int db_graph_get_perfect_path_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit, void (*node_action)(dBNode * node),
//							 dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
//							 boolean * is_cycle, dBGraph * db_graph, EdgeArrayType type, int index);


int db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit,
                                                                         Nucleotide fst_nucleotide,
                                                                         void (*node_action)(dBNode * node),
                                                                         dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
                                                                         char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
                                                                         boolean * is_cycle, dBGraph * db_graph, EdgeArrayType type, int index);



int db_graph_get_perfect_path_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit,
                                                         void (*node_action)(dBNode * node),
                                                         dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
                                                         char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
                                                         boolean * is_cycle, dBGraph * db_graph, EdgeArrayType type, int index);




boolean db_graph_detect_bubble_for_specific_person_or_population(dBNode * node,
                                                                 Orientation orientation,
                                                                 int limit, int delta,
                                                                 void (*node_action)(dBNode * node),
                                                                 int * length1, Nucleotide * base1, dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,char * seq1,
                                                                 int * length2, Nucleotide * base2, dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,char * seq2,
                                                                 dBGraph * db_graph, EdgeArrayType type, int index);




//int db_graph_supernode_for_specific_person_or_pop(dBNode * node,int limit, boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node),
//						  char * string,dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
//						  dBGraph * db_graph, EdgeArrayType type, int index);

//int db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(dBNode * node,int limit, boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node),
//                                                                           char * string,dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, int* query_node_posn,
//                                                                          dBGraph * db_graph, EdgeArrayType type, int index);


int db_graph_supernode_for_specific_person_or_pop(dBNode * node,int limit,void (*node_action)(dBNode * node), 
						  dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
						  char * supernode_str, double * avg_coverage,int * min,int * max, boolean * is_cycle,
						  dBGraph * db_graph, EdgeArrayType type, int index);



int db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(dBNode * node,int limit,void (*node_action)(dBNode * node),
                                                                            dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
                                                                            char * supernode_str, double * avg_coverage,int * min,int * max, boolean * is_cycle,
                                                                            int* query_node_posn,
                                                                            dBGraph * db_graph, EdgeArrayType type, int index);



//***********************
//functions applied to node, with respect to an individual or population
//***********************



//char * get_seq_from_elem_to_end_of_supernode_for_specific_person_or_pop(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle, char * seq, int max_length, 
//									EdgeArrayType type, int index);

//void print_supernode_for_specific_person_or_pop(HashTable* db_graph, dBNode * node,EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test);
void db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop(HashTable* db_graph, dBNode * node, long* supernode_count, EdgeArrayType type, int index,
										    boolean is_for_testing, char** for_test, int* index_for_test);

void db_graph_choose_output_filename_and_print_potential_transloc_for_specific_person_or_pop(HashTable* db_graph, dBNode * node, long* supernode_count, EdgeArrayType type, int index,
                                                                                             int min_required_covg, int max_required_covg,
                                                                                             boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2);

void db_graph_choose_output_filename_and_print_potential_inversion_for_specific_person_or_pop(HashTable* db_graph, dBNode * node, long* supernode_count, EdgeArrayType type, int index,
                                                                                              int min_required_covg, int max_required_covg,
                                                                                              boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2);


void db_graph_print_chrom_intersections_for_supernode(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test);


// functions applied to a person/population's graph
void db_graph_traverse_specific_person_or_pop_for_supernode_printing(void (*f)(HashTable*, Element *, long*, EdgeArrayType, int, boolean, char**, int*),HashTable * hash_table, long* supernode_count, 
								     EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test);

void db_graph_traverse_specific_person_or_pop_for_supernode_and_chromosome_overlap_printing(void (*f)(HashTable*, Element *, long* , EdgeArrayType, int, int, int, boolean, char**, char**, int*, int*),
                                                                                            HashTable * hash_table, long* supernode_count, EdgeArrayType type, int index, int min_covg, int max_covg,
											    boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2);

void db_graph_traverse_to_gather_statistics_about_people(void (*f)(HashTable*, Element *, int**, int),HashTable *, int**, int);

void db_graph_print_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test );

void db_graph_print_chrom_intersections_for_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index,
										 boolean is_for_testing, char** for_test, int* index_for_test);

void db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population(dBGraph* hash_table, EdgeArrayType type, int index);


dBNode * db_graph_get_next_node_for_specific_person_or_pop(dBNode * current_node, Orientation current_orientation,
                                                           Orientation * next_orientation,
                                                           Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, EdgeArrayType type, int index);

//dBNode * db_node_get_next_node_for_specific_person_or_pop(dBNode * current_node, Orientation current_orientation,
//							   Orientation * next_orientation,
//							   Nucleotide edge, Nucleotide * reverse_edge, HashTable * db_graph, EdgeArrayType type, int index);

//returns true if the node side defined by the orientation is a conflict 
//or doesn't have any outgoing edge
boolean db_node_is_supernode_end(dBNode* node, Orientation orientation, EdgeArrayType edge_type, int edge_index, dBGraph* db_graph);


// given a node, find the start of the supernode that contains it. Then return the sequence for the subsection of that supernode starting at index start, and ending at end
// in the pre-malloc-ed char* which was passed in as the frst argument
int db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(char* subsection, dBNode* node, int start, int end, EdgeArrayType type, int index, dBGraph* db_graph);
void db_graph_get_subsection_of_supernode_containing_given_node_as_supernode_object(dBSupernode* subsection, dBNode* node, int start, int end, EdgeArrayType type, int index, dBGraph* db_graph);


//TODO - implement these!!
void apply_function_to_every_node_in_supernode_object(dBSupernode*, void (*f)(dBNode*));
void apply_function_to_every_node_in_supernode_containing_given_node(void (*f)(dBNode*, EdgeArrayType, int, dBGraph*), dBNode* node);

void db_graph_set_status_of_all_nodes_in_supernode(dBNode* node, NodeStatus status, EdgeArrayType type, int index, dBGraph* db_graph);

void db_graph_get_min_and_max_covg_of_nodes_in_supernode_for_specific_person_or_pop(dBNode* node, /*NodeStatus status,*/ EdgeArrayType type, int index,  dBGraph* dbgraph, int* min_covg, int* max_covg);

dBNode* db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(dBNode* node, EdgeArrayType type, int index, dBGraph* db_graph);
dBNode* db_graph_get_next_node_in_supernode_for_specific_person_or_pop(dBNode* node, Orientation o, Orientation* next_orientation, EdgeArrayType type, int index, dBGraph* db_graph);


void db_graph_print_supernode_if_is_potential_transloc_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index,
                                                                                  int min_required_covg, int max_required_covg,
                                                                                  boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2 );

void db_graph_print_supernode_if_is_potential_inversion_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index,
                                                                                   int min_required_covg, int max_required_covg,
                                                                                   boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2 );

boolean db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome(dBNode* node, EdgeArrayType type, int index, dBGraph* dbgraph, int* total_number_of_different_chromosomes_intersected);



//returns for this function are in the two int*'s.
void  db_graph_get_best_sub_supernode_given_min_covg_and_length_for_specific_person_or_pop(dBNode* first_node_in_supernode,  int* index_for_start_of_sub_supernode,
											   int* length_of_best_sub_supernode, int min_people_coverage, 
											   int min_length_of_sub_supernode, EdgeArrayType type, int index, dBGraph* db_graph);

void  db_graph_find_population_consensus_supernode_based_on_given_node(Sequence* pop_consensus_supernode, int max_length_of_supernode, dBNode* node, 
								       int min_covg_for_pop_supernode, int min_length_for_pop_supernode, dBGraph* db_graph);



void print_node_to_file_according_to_how_many_people_share_it(HashTable* db_graph, dBNode * node, FILE** list_of_file_ptrs);
void find_out_how_many_individuals_share_this_node_and_add_to_statistics(HashTable* db_graph, dBNode * node, int** array_of_counts, int number_of_people);



//the intention is that this function would be called repeatedly, each time returning more nodes corresponding to the path through
//the graph produced by the fasta.
int db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(FILE* chrom_fptr,
                                                                                                     int number_of_nodes_to_load, int number_of_nodes_loaded_last_time,
                                                                                                     int length_of_arrays,
                                                                                                     dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
                                                                                                     Sequence* seq, KmerSlidingWindow* kmer_window,
                                                                                                     boolean expecting_new_fasta_entry, boolean expecting_seq_to_have_extra_kmer_at_start,
                                                                                                     dBGraph* db_graph );

int db_graph_make_reference_path_based_sv_calls(FILE* chrom_fasta_fptr, EdgeArrayType which_array_holds_indiv, int index_for_indiv_in_edge_array,
                                                int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, int max_anchor_span, int min_covg, int max_covg,
                                                int max_expected_size_of_supernode, int length_of_arrays, dBGraph* db_graph, FILE* output_file);


#endif
