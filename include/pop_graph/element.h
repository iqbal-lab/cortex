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

extern int NUMBER_OF_INDIVIDUALS_PER_POPULATION;
extern int NUMBER_OF_POPULATIONS;



//type definitions

typedef char Edges;

typedef enum
  {
    none    = 0,
    visited = 1,
    pruned_from_NA12878  = 2,
    pruned_from_NA12891  = 3,
    pruned_from_NA12892  = 4,
    pruned_from_NA12878_and_NA12891  = 5,
    pruned_from_NA12878_and_NA12892  = 6,
    pruned_from_NA12891_and_NA12892  = 7,
    pruned_from_NA12878_and_NA12891_and_NA12892  = 8,

  } NodeStatus;


//we will be dealing with two arrays of edges. One in an array of edges, one for each individual in  population
// and the other is an array of edges, one per population. However as far as our algorithms go, the edges
// can be treated identically, so we set up this enum, to allow you to specify which type of array you are talking about 

typedef enum{
  individual_edge_array = 0,
  population_edge_array = 1,  
} EdgeArrayType;

typedef struct{
  BinaryKmer kmer;
  char population_proportions;
  Edges* individual_edges; // array of edges of individuals in population that is currently being loaded/analysed
                           //These will be allocated for all individuals in a population, then reused for the next population
                           //In some sense these are *temporary*. You load all individuals for a population, process the graph to choose specific
                           //nodes/supernodes that you like, and then use them to specify the population_edge char for that population,
                           //You then move on to the next population, and reuse this array for the next population's individuals.
  Edges* population_edges; //array of edges, one per population, each one is union of edges in that population that satisfy some constraint.
  NodeStatus status;
} Element;


typedef Element dBNode;
typedef BinaryKmer Key;

typedef enum{
  forward = 0,
  reverse = 1
} Orientation;

typedef Element GraphNode;

Element* new_element();
void free_element(Element** element);


//utility function for getting the desired edge char, by specifying if talking about a population or an individual
// and giving the appropriate index in the relevant array

Edges* get_edge(Element, EdgeArrayType, int); //gets pointer to actual edge, so you can modify it
Edges get_edge_copy(const Element e, EdgeArrayType type,int index); //gets copy of edge
Edges get_union_of_edges(Element e);
void add_edges(Element*, EdgeArrayType, int, Edges);
void set_edges(Element*, EdgeArrayType, int, Edges);
void reset_one_edge(Element* e, Orientation orientation, Nucleotide nucleotide, EdgeArrayType type, int index);

int element_get_number_of_people_or_pops_containing_this_element(Element* e, EdgeArrayType type, int index);


boolean element_smaller(Element,Element);
BinaryKmer element_get_kmer(Element *);
boolean element_is_key(Key,Element, short kmer_size);
Key element_get_key(BinaryKmer,short kmer_size);
void element_initialise(Element *,Key, short kmer_size);

//reverse orientation
Orientation opposite_orientation(Orientation);
Orientation db_node_get_orientation(BinaryKmer, dBNode *, short kmer_size);

//add an edge between nodes -- NB: it adds both edges: forward and reverse
boolean db_node_add_edge(dBNode *, dBNode *, Orientation, Orientation, short kmer_size, EdgeArrayType edge_type, int edge_index); 


//returns yes if the label defined by the nucleotide coresponds to an 
//outgoing edge in the side defined by the orientation.   
boolean db_node_edge_exist(dBNode *, Nucleotide, Orientation, EdgeArrayType edge_type, int edge_index);

//returns the label of the first outgoing edge -- leaving from the side 
//defined by orientation. 
boolean db_node_has_precisely_one_edge(dBNode *, Orientation, Nucleotide *, EdgeArrayType edge_type, int edge_index);

boolean db_node_has_precisely_one_edge_in_union_graph_over_all_people(dBNode * node, Orientation orientation, Nucleotide * nucleotide);


void db_node_reset_all_edges_for_all_people_and_pops_to_zero(Element* e);

//forgets about the edges
void db_node_reset_edges(dBNode *, EdgeArrayType, int  );

void db_node_reset_edge(dBNode *, Orientation, Nucleotide, EdgeArrayType, int );

//check that the edges are 0's
boolean db_node_edges_reset(dBNode * node, EdgeArrayType edge_type, int edge_index);

boolean db_node_check_status(dBNode * node, NodeStatus status);
boolean db_node_check_status_not_pruned(dBNode * node);


void db_node_set_status(dBNode * node,NodeStatus status);
void db_node_trio_aware_set_pruned_status(dBNode * node, int index);
void db_node_set_status_to_none(dBNode * node);


//check if node doesn't have any edges in a given orientation
boolean db_node_is_blunt_end(dBNode * node, Orientation orientation, EdgeArrayType edge_type, int edge_index);


boolean db_node_is_this_node_in_this_person_or_populations_graph(dBNode* node, EdgeArrayType type, int index);

#endif /* ELEMENT_H_ */
