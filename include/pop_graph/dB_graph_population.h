/*
  dB_graph_population.h

  wrapper for hash_table for people who want to think of the hash table as 
  many de bruijn graphs drawn on top of each other.

*/

#ifndef DB_GRAPH__POPULATION_H_
#define DB_GRAPH_POPULATION_H_

#include <element.h>
#include <hash_table.h>
#include <seq.h>
#include <dB_graph_supernode.h>

//***********************
//functions applied to node, with respect to an individual or population
//***********************
boolean db_graph_is_this_node_in_this_person_or_populations_graph(dBNode* node, EdgeArrayType type, int index);

dBNode * db_graph_get_next_node_for_specific_person_or_pop(dBNode * current_node, Orientation current_orientation,
							   Orientation * next_orientation,
							   Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, EdgeArrayType type, int index);

dBNode *  db_graph_find_node_restricted_to_specific_person_or_population(Key key, dBGraph * hash_table, EdgeArrayType type, int index);

char * get_seq_from_elem_to_end_of_supernode_for_specific_person_or_pop(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle, EdgeArrayType type, int index);

void print_supernode_for_specific_person_or_pop(HashTable* db_graph, dBNode * node,EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test);


// functions applied to a person/population's graph
void db_graph_traverse_specific_person_or_pop(void (*f)(HashTable*, Element *, EdgeArrayType, int, boolean, char**, int*),HashTable * hash_table, EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test);
void db_graph_print_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test );

void db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population(dBGraph* hash_table, EdgeArrayType type, int index);


// given a node, find the start of the supernode that contains it. Then return the sequence for the subsection of that supernode starting at index start, and ending at end
// in the pre-malloc-ed char* which was passed in as the frst argument
void db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(char* subsection, dBNode* node, int start, int end, EdgeArrayType type, int index, dBGraph* db_graph);
void db_graph_get_subsection_of_supernode_containing_given_node_as_supernode_object(dBSupernode* subsection, dBNode* node, int start, int end, EdgeArrayType type, int index, dBGraph* db_graph);

void apply_function_to_every_node_in_supernode_object(dBSupernode*, void (*f)(dBNode*));
void apply_function_to_every_node_in_supernode_containing_given_node(void (*f)(dBNode*, EdgeArrayType, int, dBGraph*), dBNode* node);


void return_node_to_status_prior_to_tmp_touched_X(dBNode*);
dBNode* db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(dBNode* node, EdgeArrayType type, int index, dBGraph* db_graph);
dBNode* db_graph_get_next_node_in_supernode_for_specific_person_or_pop(dBNode* node, Orientation o, Orientation* next_orientation, EdgeArrayType type, int index, dBGraph* db_graph);


//returns for this function are in the two int*'s.
void  db_graph_get_best_sub_supernode_given_min_covg_and_length_for_specific_person_or_pop(const dBNode* first_node_in_supernode,  int* index_for_start_of_sub_supernode,
											   int* length_of_best_sub_supernode, int min_people_coverage, 
											   int min_length_of_sub_supernode, EdgeArrayType type, int index, dBGraph* db_graph);

void  db_graph_find_population_consensus_supernode_based_on_given_node(Sequence* pop_consensus_supernode, dBNode* node, int min_covg_for_pop_supernode, int min_length_for_pop_supernode, dBGraph* db_graph);



#endif
