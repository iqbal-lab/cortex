/*
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
  element.h defines the interface for the de Bruijn graph node. The implementation is complemented by 
  a hash table that stores every node indexed by kmers (BinaryKmers). 

  The element routines, ie the one required by hash_table/priority queue, are prefixed with element_ 
  The de Bruijn based routines are prefixed with db_node
*/


#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <binary_kmer.h>
#include <global.h>
#include <stdio.h>

//type definitions

typedef char Edges;




typedef unsigned short Flags;
#define  ALL_OFF  		       0

#define  USED			 (1 << 0) 
#define  VISITED		 (1 << 1) 
#define  PRUNED		         (1 << 2) 
#define  STOP_PATH		 (1 << 3) 
#define  VISITED_FF		 (1 << 4) 
#define  VISITED_RR		 (1 << 5) 
#define  VISITED_CLEANING	 (1 << 6) 




typedef struct{
	BinaryKmer kmer;
	int coverage;
	Flags flags; 
	Edges edges; // less significant nibble forward
	
} Element;


typedef Element dBNode;

typedef BinaryKmer* Key;


typedef Element GraphNode;


void element_assign(Element* e1, Element* e2);

//reverse orientation
Orientation opposite_orientation(Orientation);

boolean element_is_key(Key,Element, short kmer_size);

Key element_get_key(BinaryKmer*,short kmer_size, Key preallocated_key);

boolean element_smaller(Element,Element);

void element_initialise(Element *,BinaryKmer* kmer, short kmer_size);

BinaryKmer* element_get_kmer(Element *);

int element_get_coverage(Element *);
int element_update_coverage(Element *, int);

Orientation db_node_get_orientation(BinaryKmer*, dBNode *, short kmer_size);

//add an edge between nodes -- NB: it adds both edges: forward and reverse
boolean db_node_add_edge(dBNode *, dBNode *, Orientation, Orientation, short kmer_size); 

//returns true if the node side defined by the orientation is a conflict 
//or doesn't have any outgoing edge
boolean db_node_is_supernode_end(dBNode *, Orientation);

//returns yes if the label defined by the nucleotide coresponds to an 
//outgoing edge in the side defined by the orientation.   
boolean db_node_edge_exist(dBNode *, Nucleotide, Orientation);

//returns the label of the first outgoing edge -- leaving from the side 
//defined by orientation. 
boolean db_node_has_precisely_one_edge(dBNode *, Orientation, Nucleotide *);

//returns the label of one edge (the first that exists in order A,C,G,T), returns false if that edge doesn't exitst
boolean db_node_get_one_edge(dBNode *, Orientation, Nucleotide *);

//returns the label of the "two edges"
//defined by orientation. 
boolean db_node_has_precisely_two_edges(dBNode *, Orientation, Nucleotide *, Nucleotide *);

void db_node_unset_flag(dBNode * node, Flags f) ;

Edges db_node_get_edges(dBNode * node);
//char db_node_get_edges_coverage(dBNode * node);
//forgets about the edges
void db_node_reset_edges(dBNode * );

void db_node_reset_edge(dBNode *, Nucleotide, Orientation);

//check that the edges are 0's
boolean db_node_edges_reset(dBNode * );

//set every edge in 'edges' 
void db_node_set_edges(dBNode * node, Edges edges);
boolean db_node_edge_exist(dBNode * element,Nucleotide base,Orientation orientation);

//boolean db_node_edge_has_single_coverage(dBNode * element,Nucleotide base,Orientation orientation);



//check if node doesn't have any edges in a given orientation
boolean db_node_is_blunt_end(dBNode * node, Orientation orientation);


void db_node_print_binary(FILE * fp, dBNode * node, int); 

boolean db_node_read_binary(FILE * fp, short kmer_size, dBNode * node); 


//actions and conditions 

void db_node_action_set_flag_none(dBNode * node);

void db_node_action_set_flag_pruned(dBNode * node);

void db_node_action_set_flag_visited(dBNode * node);

void db_node_action_set_flag_visited_cleaning(dBNode * node);

void db_node_action_do_nothing(dBNode * node);

boolean db_node_check_flag_not_pruned(dBNode * node);

boolean db_node_check_flag_none(dBNode * node);

boolean db_node_check_nothing(dBNode * node);

boolean db_node_check_flag_visited(dBNode * node);

boolean db_node_update_edge_count(dBNode * node,Orientation o, Nucleotide base);

int db_node_edges_count(dBNode * node, Orientation orientation);

boolean db_node_condition_always_true(dBNode* node);


void db_node_set_flag(dBNode * node, Flags f);
Flags db_node_get_flags(dBNode * node, Flags f);
boolean db_node_check_for_flag(dBNode * node, Flags flag);
boolean db_node_check_for_flag_ALL_OFF(dBNode * node);

#endif /* ELEMENT_H_ */
