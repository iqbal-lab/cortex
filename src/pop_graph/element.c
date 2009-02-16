/*
 element.c -- implements the nodes of the dBruijn graph
 */

#include <element.h>
#include <stdlib.h>
#include <global.h>
#include <stdio.h>


int NUMBER_OF_INDIVIDUALS_PER_POPULATION = 5;
int NUMBER_OF_POPULATIONS = 2;


//currently noone calls this in normal use
// In normal use, the priority queue allocates space to put the eloement directly within,
// and calls element_initialise which allocs the arrays which the elemtn has pointers to
Element* new_element()
{
  Element* e = malloc(sizeof(Element));

  if (e==NULL)
    {
      printf("Unable to allocate a new element");
      exit(1);
    }
  
  e->kmer=0;
  e->population_proportions=0;


  e->individual_edges = malloc(sizeof(Edges)*NUMBER_OF_INDIVIDUALS_PER_POPULATION);
  if (e->individual_edges == NULL)
    {
      printf("cannot alloc the edges for individuals for this element");
      exit(1);
    }

  int i;
  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
    }


  e->population_edges = malloc(sizeof(Edges)*NUMBER_OF_POPULATIONS);
  if (e->population_edges == NULL)
    {
      printf("cannot alloc the edges for populations for this element");
      exit(1);
    }
  for (i=0; i< NUMBER_OF_POPULATIONS; i++)
    {
      e->population_edges[i]=0;
    }


  e->status=none;

  return e;
}


void free_element(Element** element)
{
  free((*element)->individual_edges);
  (*element)->individual_edges = NULL;
  free((*element)->population_edges);
  (*element)->population_edges = NULL;
  free(*element);
  *element=NULL;
}

//gets you a pointer to the edge you are referring to
Edges* get_edge(Element e, EdgeArrayType type,int index)
{

  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}
      return &e.individual_edges[index];
    }
  else if (type == population_edge_array)
    {
      if (index >= NUMBER_OF_POPULATIONS)
	{
	  exit(1);
	}
      return &e.population_edges[index];
    }

  exit(1);
}


//return a copy of the edge you are referring to
Edges get_edge_copy(const Element e, EdgeArrayType type,int index)
{

  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}
      return e.individual_edges[index];
    }
  else if (type == population_edge_array)
    {
      if (index >= NUMBER_OF_POPULATIONS)
	{
	  exit(1);
	}
      return e.population_edges[index];
    }

  exit(1);
}



Edges get_union_of_edges(Element e)
{

  int i;
  Edges edges=0;

  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      edges |= e.individual_edges[i];
    }

  for (i=0; i< NUMBER_OF_POPULATIONS; i++)
    {
      edges |=  e.population_edges[i];
    }

  return edges;
}








//adds edges from edge_char to the appropriate person/population edgeset, without removing existing edges
void add_edges(Element* e, EdgeArrayType type, int index, Edges edge_char)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}
      e->individual_edges[index] |= edge_char;
    }
  else if (type == population_edge_array)
    {
      if (index >= NUMBER_OF_POPULATIONS)
	{
	  exit(1);
	}
      e->population_edges[index] |= edge_char;
    }
  
}


void set_edges(Element* e, EdgeArrayType type, int index, Edges edge_char)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}
      e->individual_edges[index] = edge_char;
    }
  else if (type == population_edge_array)
    {
      if (index >= NUMBER_OF_POPULATIONS)
	{
	  exit(1);
	}
      e->population_edges[index] = edge_char;
    }
  
}


void db_node_reset_all_edges_for_all_people_and_pops_to_zero(Element* e)
{
  int i;

    for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
    }


  for (i=0; i<NUMBER_OF_POPULATIONS; i++)
    {
      e->population_edges[i]=0;
    }

}
void reset_one_edge(Element* e, Orientation orientation, Nucleotide nucleotide, EdgeArrayType type, int index)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}

      char edge = 1 << nucleotide;      
      if (orientation == reverse){
	edge <<= 4;
      }
      //toggle 1->0 0->1
      edge ^= (unsigned char) 0xFF; //xor with all 1's, ie 00010000 -> 11101111
      
      e->individual_edges[index] &= edge; //reset one edge
      

    }
  else if (type == population_edge_array)
    {
      if (index >= NUMBER_OF_POPULATIONS)
	{
	  exit(1);
	}
      char edge = 1 << nucleotide;

      if (orientation == reverse){
	edge <<= 4;
      }

      //toggle 1->0 0->1

      edge ^= (unsigned char) 0xFF; //xor with all 1's, ie 00010000 -> 11101111

      e->population_edges[index] &= edge; //reset one edge

    }
  
}


int element_get_number_of_people_or_pops_containing_this_element(Element* e, EdgeArrayType type, int index)
{
  int i;
  int count=0;
  if (type == individual_edge_array)
    {
      for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
	{
	  if ( (e->individual_edges)[i] != 0)
	    {
	      count++;
	    }
	}
    }
  else if (type == population_edge_array)
    {
      for (i=0; i< NUMBER_OF_POPULATIONS; i++)
        {
          if ( (e->population_edges)[i] != 0)
            {
	      count++;
            }
        }

    }

  return count;
}

boolean element_smaller(Element  e1, Element e2){
 
  return get_union_of_edges(e1)  <  get_union_of_edges(e2);
  
}



BinaryKmer element_get_kmer(Element * e){
  return e->kmer;
}

boolean element_is_key(Key key, Element e, short kmer_size){
  return key == e.kmer;
}

Key element_get_key(BinaryKmer kmer, short kmer_size){
  
  BinaryKmer rev_kmer = binary_kmer_reverse_complement(kmer,kmer_size);
  
  if (rev_kmer < kmer){
    kmer = rev_kmer;
  }

  return kmer;

}

//this allocs the population and person arrays
void element_initialise(Element * e, Key kmer, short kmer_size){

  e->kmer = element_get_key(kmer, kmer_size);

  e->individual_edges = malloc(sizeof(Edges)*NUMBER_OF_INDIVIDUALS_PER_POPULATION);
  if (e->individual_edges == NULL)
    {
      printf("cannot alloc the edges for individuals for this element");
      exit(1);
    }
  int i;
  for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
    }

  e->population_edges = malloc(sizeof(Edges)*NUMBER_OF_POPULATIONS);
  if (e->population_edges == NULL)
    {
      printf("cannot alloc the edges for populations for this element");
      exit(1);
    }

  for (i=0; i<NUMBER_OF_POPULATIONS; i++)
    {
      e->population_edges[i]=0;
    }

  e->population_proportions=0;
  e->status = none;
}


Orientation opposite_orientation(Orientation o){
  return o ^ 1;
  
}

Orientation db_node_get_orientation(BinaryKmer k, dBNode * e, short kmer_size){

  if (e->kmer == k){
    return forward;
  }
  
  if (e->kmer == binary_kmer_reverse_complement(k,kmer_size)){
    return reverse;
  }

  exit(1);
  
}





//After specifying which individual or population you are talking about, this function
//adds one edge ("arrow") to the appropriate edge in the appropriate array in the element -- basically sets a bit in the correct edges char
void db_node_add_labeled_edge(dBNode * e, Orientation o, Nucleotide base, EdgeArrayType edge_type, int edge_index){

  //set edge 
  char edge = 1 << base; // A (0) -> 0001, C (1) -> 0010, G (2) -> 0100, T (3) -> 1000
  
  if (o == reverse){
    edge <<= 4; //move to next nibble 
  }

  //update node
  add_edges(e, edge_type, edge_index, edge);
  
}


//adding an edge between two nodes implies adding two labeled edges (one in each direction)
//be aware that in the case of self-loops in palindromes the two labeled edges collapse in one

boolean db_node_add_edge(dBNode * src_e, dBNode * tgt_e, Orientation src_o, Orientation tgt_o, short kmer_size, EdgeArrayType edge_type, int edge_index){

  BinaryKmer src_k, tgt_k; 

  src_k = src_e->kmer;
  tgt_k = tgt_e->kmer;
 
  if (src_o == reverse){
    src_k = binary_kmer_reverse_complement(src_k,kmer_size);
  }
    
  if (tgt_o == reverse){
    tgt_k = binary_kmer_reverse_complement(tgt_k,kmer_size);
  }
    
  
  if (DEBUG){
    char dummy1[kmer_size];
    char dummy2[kmer_size];

    printf("add edge %s -%c-> %s to edge type %d, and edge index %d\n",binary_kmer_to_seq(src_k,kmer_size, dummy1),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(tgt_k)),
	   binary_kmer_to_seq(tgt_k,kmer_size, dummy2), edge_type, edge_index);
  }

  db_node_add_labeled_edge(src_e,src_o,binary_kmer_get_last_nucleotide(tgt_k), edge_type, edge_index);

  if (DEBUG){
    char dummy3[kmer_size];
    char dummy4[kmer_size];
    printf("add edge %s -%c-> %s to edge type %d, and edge index %d\n",binary_kmer_to_seq(tgt_k,kmer_size,dummy3),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(src_k,kmer_size))),binary_kmer_to_seq(src_k,kmer_size, dummy4),  edge_type, edge_index);
  }

  db_node_add_labeled_edge(tgt_e,opposite_orientation(tgt_o),binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(src_k,kmer_size)), edge_type, edge_index );

  return true;
}



boolean db_node_edge_exist(dBNode * element,Nucleotide base,Orientation orientation, EdgeArrayType edge_type, int edge_index){

  //get the edge char for this specific person or pop:
  char edge = get_edge_copy(*element, edge_type, edge_index);


  edge >>= base;
  if (orientation == reverse){
    edge >>= 4;
  }
  
  edge &= 1;
  
  if (edge == 1){
    return true;
  }
  else{
    return false;
  }
}




void db_node_reset_edges(dBNode * node,EdgeArrayType edge_type, int edge_index){
  set_edges(node, edge_type, edge_index, 0);
}

void db_node_reset_edge(dBNode * node, Orientation orientation, Nucleotide nucleotide, EdgeArrayType edge_type, int edge_index){
  reset_one_edge(node, orientation, nucleotide, edge_type, edge_index);
}




boolean db_node_edges_reset(dBNode * node, EdgeArrayType edge_type, int edge_index){
  return get_edge_copy(*node,edge_type,edge_index) == 0;
}


boolean db_node_has_precisely_one_edge(dBNode * node, Orientation orientation, Nucleotide * nucleotide, EdgeArrayType edge_type, int edge_index){
  
  Nucleotide n;
  Edges edges = get_edge_copy(*node,edge_type,edge_index);
  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
 
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      *nucleotide = n;
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 1);
  
}


boolean db_node_has_precisely_one_edge_in_union_graph_over_all_people(dBNode * node, Orientation orientation, Nucleotide * nucleotide){
  
  Nucleotide n;
  Edges edges = get_union_of_edges(*node);
  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
 
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      *nucleotide = n;
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 1);
  
}


boolean db_node_is_blunt_end(dBNode * node, Orientation orientation, EdgeArrayType edge_type, int edge_index){
  
  Edges edges = get_edge_copy(*node, edge_type, edge_index);


  if (orientation == reverse){
    edges >>= 4;
  }
  
  edges &= 15; // AND with 00001111 so that we only look at the 4 least significant bits
  
  return edges == 0;
}

boolean db_node_check_status(dBNode * node, NodeStatus status){
  return node->status == status;
}
boolean db_node_check_status_not_pruned(dBNode * node){
  if ( db_node_check_status(node, none) || db_node_check_status(node,visited))
    {
      return true;
    }
  return false;
}


void db_node_set_status(dBNode * node,NodeStatus status){
  node->status = status;
}
void db_node_set_status_to_none(dBNode * node){
  node->status = none;
}

//assumes index 0 = NA12878, index 1 = NA12891, index 2 = NA12892
//semantics - call this when pruning node from person defined by index
void db_node_trio_aware_set_pruned_status(dBNode * node, int index)
{
  if (db_node_check_status(node,none) || db_node_check_status(node,visited))
    {

  	  if (index==0)//NA12878
	    {
	      db_node_set_status(node, pruned_from_NA12878);
	    }
	  else if (index==1) //NA12891
	    {
	      db_node_set_status(node, pruned_from_NA12891);
	    }
	  else if (index==2)//NA12892
	    {
	      db_node_set_status(node, pruned_from_NA12892);
	    }
    }

  if (db_node_check_status(node, pruned_from_NA12878))
    {
      if (index==0)//NA12878
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==1) //NA12891
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891);
	}
      else if (index==2)//NA12892
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12892);
	}
    }
  
  else if (db_node_check_status(node, pruned_from_NA12891))
    {
      if (index==0)//NA12878
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891);
	  
	}
      else if (index==1) //NA12891
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==2)//NA12892
	{
	  db_node_set_status(node, pruned_from_NA12891_and_NA12892);
	}
    }
  else if (db_node_check_status(node, pruned_from_NA12892))
    {
      if (index==0)//NA12878
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12892);
	  
	}
      else if (index==1) //NA12891
	{
	  db_node_set_status(node, pruned_from_NA12891_and_NA12892);
	}
      else if (index==2)//NA12892
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
    }
  else if (db_node_check_status(node, pruned_from_NA12878_and_NA12891))
    {
      if (index==0)//NA12878
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==1) //NA12891
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==2)//NA12892
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891_and_NA12892);
	}
    }
  
  else if (db_node_check_status(node, pruned_from_NA12878_and_NA12892))
    {
      if (index==0)//NA12878
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==1) //NA12891
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891_and_NA12892);
	}
      else if (index==2)//NA12892
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
    }
  else if (db_node_check_status(node, pruned_from_NA12891_and_NA12892))
    {
      if (index==0)//NA12878
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891_and_NA12892);
	}
      else if (index==1) //NA12891
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==2)//NA12892
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
    }
}



boolean db_node_is_this_node_in_this_person_or_populations_graph(dBNode* node, EdgeArrayType type, int index)
{

  if (node ==NULL)
    {
      return false;
    }

  Edges edge_for_this_person_or_pop = get_edge_copy(*node, type, index);

  if (edge_for_this_person_or_pop == 0)
    {

      return false;
    }
  else
    {
      return true;
    }
 
}
