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


//type definitions

typedef char Edges;

typedef enum{
  none    = 0,
  visited = 1,
  pruned  = 2,
} NodeStatus;

typedef struct{
	BinaryKmer kmer;
	Edges edges; // less significant nibble forward
	NodeStatus status;
} Element;

typedef Element dBNode;

typedef BinaryKmer Key;

typedef enum{
  forward = 0,
  reverse = 1
} Orientation;

typedef Element GraphNode;

//reverse orientation
Orientation opposite_orientation(Orientation);

boolean element_is_key(Key,Element, short kmer_size);

Key element_get_key(BinaryKmer,short kmer_size);

boolean element_smaller(Element,Element);

void element_initialise(Element *,Key, short kmer_size);

BinaryKmer element_get_kmer(Element *);

Orientation db_node_get_orientation(BinaryKmer, dBNode *, short kmer_size);

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

//forgets about the edges
void db_node_reset_edges(dBNode * );

void db_node_reset_edge(dBNode *, Orientation, Nucleotide );

//check that the edges are 0's
boolean db_node_edges_reset(dBNode * );


boolean db_node_check_status(dBNode * node, NodeStatus status);

void db_node_set_status(dBNode * node,NodeStatus status);

//check if node doesn't have any edges in a given orientation
boolean db_node_is_blunt_end(dBNode * node, Orientation orientation);


#endif /* ELEMENT_H_ */
