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
  used when genotyping. one colour per allele, plus one for reference plus two spare for working
*/

#ifndef GENOTYPING_ELEMENT_H_
#define GENOTYPING_ELEMENT_H_


#include <binary_kmer.h>
#include <global.h>
#include <stdio.h>
#include <stdint.h>


//Bubble and PD callers only call biallelic for now
#define MAX_ALLELES_SUPPORTED_FOR_STANDARD_GENOTYPING 2




//type definitions


typedef struct{
  BinaryKmer kmer;
  //one for each allele + all your samples/colours + 2 working colours - assume one of NUMBER_OF_COLOURS is the ref
  //preserve all the same colour numbers in the main dBGraph. So colour 3 in a normal element
  // will be colour 3 in the GenotypingElement. AFTER them will come the 2 allele colours and the
  // 2 working colours
  int        coverage[MAX_ALLELES_SUPPORTED_FOR_STANDARD_GENOTYPING+NUMBER_OF_COLOURS+2];
  Edges      individual_edges[MAX_ALLELES_SUPPORTED_FOR_STANDARD_GENOTYPING+NUMBER_OF_COLOURS+2];
  char       status; //will case a NodeStatus to char
} GenotypingElement;


typedef GenotypingElement dBGenotypingNode;


GenotypingElement* new_genotyping_element();
void free_genotyping_element(GenotypingElement** genotyping_element);
void genotyping_element_assign(GenotypingElement* e1, GenotypingElement* e2);
void genotyping_element_initialise_from_normal_element(GenotypingElement* e1, 
							     Element* e2,
							     boolean set_status_to_none);

Edges* db_genotyping_node_get_edge(GenotypingElement, EdgeArrayType, int); 
Edges  db_genotyping_node_get_edge_copy(const GenotypingElement e, EdgeArrayType type,int index); 
Edges db_genotyping_node_get_union_of_edges(GenotypingElement e);
void  db_genotyping_node_add_edges(GenotypingElement*, EdgeArrayType, int, Edges);
void  db_genotyping_node_set_edges(GenotypingElement*, EdgeArrayType, int, Edges);
void  db_genotyping_node_reset_one_edge(GenotypingElement* e, Orientation orientation, Nucleotide nucleotide, EdgeArrayType type, int index);



Edges genotyping_element_get_colour_union_of_all_colours(const GenotypingElement*);
Edges genotyping_element_get_colour0(const GenotypingGenotyping_Element* e);
Edges genotyping_element_get_colour1(const GenotypingElement* e);
Edges genotyping_element_get_last_colour(const GenotypingElement* e);

int genotyping_element_get_covg_union_of_all_covgs(const dBGenotypingNode*);
int genotyping_element_get_covg_colour0(const dBGenotypingNode* e);
int genotyping_element_get_covg_colour1(const dBGenotypingNode* e);
int genotyping_element_get_covg_last_colour(const dBGenotypingNode* e);


int genotyping_element_get_number_of_people_or_pops_containing_this_element(GenotypingElement* e, EdgeArrayType type, int index);


boolean genotyping_element_smaller(GenotypingElement,GenotypingElement);
BinaryKmer* genotyping_element_get_kmer(GenotypingElement *);
boolean genotyping_element_is_key(Key,GenotypingElement, short kmer_size);
Key genotyping_element_get_key(BinaryKmer*,short kmer_size, Key preallocated_key);
void genotyping_element_initialise(GenotypingElement *,Key, short kmer_size);
void genotyping_element_initialise_kmer_covgs_edges_and_status_to_zero(GenotypingElement * e);

void genotyping_element_set_kmer(GenotypingElement * e, Key kmer, short kmer_size);


Orientation db_genotyping_node_get_orientation(BinaryKmer*, dBGenotypingNode *, short kmer_size);

//add an edge between nodes -- NB: it adds both edges: forward and reverse
boolean db_genotyping_node_add_edge(dBGenotypingNode *, dBGenotypingNode *, Orientation, Orientation, short kmer_size, EdgeArrayType edge_type, int edge_index); 


//returns yes if the label defined by the nucleotide coresponds to an 
//outgoing edge in the side defined by the orientation.   
boolean db_genotyping_node_edge_exist(dBGenotypingNode *, Nucleotide, Orientation, EdgeArrayType edge_type, int edge_index);

//final argument f is a function that returns an Edge that is a function of the different colured edges in a node.                                                                                       // eg we might be interested in the union of all the coloured edges in the graph, or just the colour/edge for the first person,                                                                          //    or the union of all edges except that corresponding to the reference.
boolean db_genotyping_node_edge_exist_within_specified_function_of_coloured_edges(dBGenotypingNode * genotyping_element,Nucleotide base,Orientation orientation, Edges (*f)(const GenotypingElement* ));


//returns the label of the first outgoing edge -- leaving from the side 
//defined by orientation. 
boolean db_genotyping_node_has_precisely_one_edge(dBGenotypingNode *, Orientation, Nucleotide *, EdgeArrayType edge_type, int edge_index);


boolean db_genotyping_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(dBGenotypingNode * node, Orientation orientation, Nucleotide * nucleotide, 
									      Edges (*get_colour)(const dBGenotypingNode*) );

boolean db_genotyping_node_has_precisely_one_edge_in_union_graph_over_all_people(dBGenotypingNode * node, Orientation orientation, Nucleotide * nucleotide);

boolean db_genotyping_node_has_precisely_two_edges(dBGenotypingNode * node, Orientation orientation, Nucleotide * nucleotide1, Nucleotide * nucleotide2, EdgeArrayType type, int index);

void db_genotyping_node_reset_all_edges_for_all_people_and_pops_to_zero(GenotypingElement* e);

//forgets about the edges
void db_genotyping_node_reset_edges(dBGenotypingNode *, EdgeArrayType, int  );

void db_genotyping_node_reset_edge(dBGenotypingNode *, Orientation, Nucleotide, EdgeArrayType, int );



//TODO - maybe do not need to export this:
void db_genotyping_node_reset_specified_edges(dBGenotypingNode * node, Orientation orientation, Nucleotide nucleotide, void (*f)(dBGenotypingNode*, Orientation, Nucleotide)  );


//check that the edges are 0's
boolean db_genotyping_node_edges_reset(dBGenotypingNode * node, EdgeArrayType edge_type, int edge_index);

boolean db_genotyping_node_check_status(dBGenotypingNode * node, NodeStatus status);
boolean db_genotyping_node_check_status_not_pruned(dBGenotypingNode * node);
boolean db_genotyping_node_check_status_not_pruned_or_visited(dBGenotypingNode * node);
boolean db_genotyping_node_check_status_to_be_dumped(dBGenotypingNode * node);
boolean db_genotyping_node_check_status_is_not_visited(dBGenotypingNode* node);


void db_genotyping_node_set_status(dBGenotypingNode * node,NodeStatus status);
void db_genotyping_node_trio_aware_set_pruned_status(dBGenotypingNode * node, int index);
void db_genotyping_node_set_status_to_none(dBGenotypingNode * node);



//actions and conditions 

void db_genotyping_node_action_set_status_none(dBGenotypingNode * node);
void db_genotyping_node_action_set_status_of_unpruned_to_none(dBGenotypingNode * node);

void db_genotyping_node_action_set_status_pruned(dBGenotypingNode * node);
void db_genotyping_node_action_set_status_visited(dBGenotypingNode * node);
void db_genotyping_node_action_set_status_ignore_this_node(dBGenotypingNode * node);

void db_genotyping_node_action_set_status_visited_or_visited_and_exists_in_reference(dBGenotypingNode * node);

void db_genotyping_node_action_unset_status_visited_or_visited_and_exists_in_reference(dBGenotypingNode * node);

void db_genotyping_node_action_unset_status_visited_or_visited_and_exists_in_reference_or_ignore_this_node(dBGenotypingNode * node);


void db_genotyping_node_action_do_nothing(dBGenotypingNode * node);

boolean db_genotyping_node_check_status_none(dBGenotypingNode * node);
boolean db_genotyping_node_check_for_flag_ALL_OFF(dBGenotypingNode * node);


boolean db_genotyping_node_check_status_visited(dBGenotypingNode * node);

boolean db_genotyping_node_check_status_exists_in_reference(dBGenotypingNode * node);

boolean db_genotyping_node_check_status_visited_and_exists_in_reference(dBGenotypingNode * node);

boolean db_genotyping_node_check_status_is_not_exists_in_reference(dBGenotypingNode * node);

boolean db_genotyping_node_check_status_is_not_visited_or_visited_and_exists_in_reference(dBGenotypingNode * node);

boolean db_genotyping_node_condition_always_true(dBGenotypingNode* node);





void db_genotyping_node_increment_coverage(dBGenotypingNode* e, EdgeArrayType type, int index);
void db_genotyping_node_decrement_coverage(dBNode* e, EdgeArrayType type, int index);
void db_genotyping_node_update_coverage(dBGenotypingNode* e, EdgeArrayType type, int index, int update);
int db_genotyping_node_get_coverage(const dBGenotypingNode* const e, EdgeArrayType type, int index);
void db_genotyping_node_set_coverage(dBGenotypingNode* e, EdgeArrayType type, int colour, int covg);
int db_genotyping_node_get_coverage_in_subgraph_defined_by_func_of_colours(const dBGenotypingNode* const e, int (*get_covg)(const dBGenotypingNode*) );



//check if node doesn't have any edges in a given orientation
boolean db_genotyping_node_is_blunt_end(dBGenotypingNode * node, Orientation orientation, EdgeArrayType edge_type, int edge_index);
boolean db_genotyping_node_is_blunt_end_in_subgraph_given_by_func_of_colours(dBGenotypingNode * node, Orientation orientation,  Edges (*get_colour)(const dBGenotypingNode*) );


boolean db_genotyping_node_is_this_node_in_this_person_or_populations_graph(dBGenotypingNode* node, EdgeArrayType type, int index);


boolean db_genotyping_node_is_this_node_in_subgraph_defined_by_func_of_colours(dBGenotypingNode* node, Edges (*get_colour)(const dBGenotypingNode*) );


//functions for binary format
void db_genotyping_node_print_multicolour_binary(FILE * fp, dBGenotypingNode * node);

void db_genotyping_node_print_single_colour_binary_of_colour0(FILE * fp, dBGenotypingNode * node);
void db_genotyping_node_print_single_colour_binary_of_specified_colour(FILE * fp, dBGenotypingNode * node, int colour);

//reading multicolour binaries
boolean db_genotyping_node_read_multicolour_binary(FILE * fp, short kmer_size, dBGenotypingNode * node, int num_colours_in_binary, int binversion_in_binheader);
//boolean db_genotyping_node_read_multicolour_binary_with_less_colours(FILE * fp, short kmer_size, dBGenotypingNode * node, int num_colours_in_binary);



//read a binary for an individual person, as dumped by the target "graph"
// the edge array type and index tell you which person you should load this data into
boolean db_genotyping_node_read_single_colour_binary(FILE * fp, short kmer_size, dBGenotypingNode * node, EdgeArrayType type, int index, int binversion_in_binheader);


boolean db_genotyping_node_check_read_start(dBGenotypingNode* node, Orientation ori);

void db_genotyping_node_set_read_start_status(dBGenotypingNode* node, Orientation ori);

boolean db_genotyping_node_check_duplicates(dBGenotypingNode* node1, Orientation o1, dBGenotypingNode* node2, Orientation o2);

//we have a read that starts at node in direction o1, and we want to know if a previous read started at that node in that direction
boolean db_genotyping_node_check_single_ended_duplicates(dBGenotypingNode* node1, Orientation o1);

#endif /* GENOTYPING_ELEMENT_H_ */
