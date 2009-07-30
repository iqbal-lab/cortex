/*
 Elemenx.c -- implements the nodes of the dBruijn graph
 */

#include <element.h>
#include <stdlib.h>
#include <global.h>
#include <stdio.h>
#include <string.h>


//const int NUMBER_OF_INDIVIDUALS_PER_POPULATION = 5;
//const int NUMBER_OF_POPULATIONS = 2;

const unsigned char mask1 = 255-3; //11111100
const unsigned char mask2 = 255-12;//11110011
const unsigned char mask3 = 255-48;//11001111
const unsigned char mask4 = 255-192;//00111111


//currently noone calls this in normal use
// In normal use, the priority queue allocates space to put the eloement directly within,
// and calls element_initialise
Element* new_element()
{


  Element* e = malloc(sizeof(Element));

  if (e==NULL)
    {
      printf("Unable to allocate a new element");
      exit(1);
    }
  
  e->kmer=0;

  int i;
  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
      e->coverage[i]=0;
    }


  e->status=none;
  e->kmer=0;

  return e;
}


void free_element(Element** element)
{
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
 else 
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array, but we are getting type %d", type);
      exit(1);
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
  else 
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array, but we are getting type %d", type);
      exit(1);
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

  return edges;
}








//adds edges from edge_char to the appropriate person/population edgeset, without removing existing edges
void add_edges(Element* e, EdgeArrayType type, int index, Edges edge_char)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  printf("in element's add_edges function. index is %d, and should be at most %d", index, NUMBER_OF_INDIVIDUALS_PER_POPULATION);
	  exit(1);
	}
      e->individual_edges[index] |= edge_char;
    }

  else
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array, but we are getting type %d", type);
      exit(1);
    }
  
}


void set_edges(Element* e, EdgeArrayType type, int index, Edges edge_char)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  printf("in element's set_edges function. index is %d,and should be at most %d", index, NUMBER_OF_INDIVIDUALS_PER_POPULATION);
	  exit(1);
	}
      e->individual_edges[index] = edge_char;
    }

  else
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array, but we are getting type %d", type);
      exit(1);
    }
  
}


void db_node_reset_all_edges_for_all_people_and_pops_to_zero(Element* e)
{
  int i;

    for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
    }

}

void reset_one_edge(Element* e, Orientation orientation, Nucleotide nucleotide, EdgeArrayType type, int index)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  printf("in element's reset_one_edge function. index is %d,and should be at most %d", index, NUMBER_OF_INDIVIDUALS_PER_POPULATION);
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
  else
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array, but we are getting type %d", type);
      exit(1);
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
  else
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array, but we are getting type %d", type);
      exit(1);
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


void element_initialise(Element * e, Key kmer, short kmer_size){

  e->kmer = element_get_key(kmer, kmer_size);

  //has table has calloc-ed all elements, so already initialised to zero.
  int i;
  for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
      e->coverage[i]=0;
    }

  db_node_set_status(e, none);

}


void db_node_increment_coverage(dBNode* e, EdgeArrayType type, int index)
{
  e->coverage[index]=e->coverage[index]+1;
}

void db_node_update_coverage(dBNode* e, EdgeArrayType type, int index, short update)
{

  e->coverage[index]=db_node_get_coverage_as_short(e,type, index) + update;

}

//coverage stored as short, but we want to deal with it as an int. So access it via this getter
int db_node_get_coverage(dBNode* e, EdgeArrayType type, int index)
{
  short c = e->coverage[index];
  return (int) c;
}

int db_node_get_coverage_as_short(dBNode* e, EdgeArrayType type, int index)
{
  return e->coverage[index];
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

  printf("Exiting from  db_node_get_orientation because the binary kmer given is neither fw nor rev complement of the kmer associated with this element");
  char tmp_seq[kmer_size];
  printf("Kmer is %s\n, and element kmer is %s\n", binary_kmer_to_seq(k, kmer_size, tmp_seq), binary_kmer_to_seq(e->kmer, kmer_size, tmp_seq) );
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
  char tmp_seq[kmer_size];
 
  if (src_o == reverse){
    src_k = binary_kmer_reverse_complement(src_k,kmer_size);
  }
    
  if (tgt_o == reverse){
    tgt_k = binary_kmer_reverse_complement(tgt_k,kmer_size);
  }
    
  
  if (DEBUG){
    printf("add edge %s -%c-> %s to edge type %d, and edge index %d\n",binary_kmer_to_seq(src_k,kmer_size, tmp_seq),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(tgt_k)),
	   binary_kmer_to_seq(tgt_k,kmer_size, tmp_seq), edge_type, edge_index);
  }

  db_node_add_labeled_edge(src_e,src_o,binary_kmer_get_last_nucleotide(tgt_k), edge_type, edge_index);

  if (DEBUG){

    printf("add edge %s -%c-> %s to edge type %d, and edge index %d\n",binary_kmer_to_seq(tgt_k,kmer_size,tmp_seq),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(src_k,kmer_size))),binary_kmer_to_seq(src_k,kmer_size, tmp_seq),  edge_type, edge_index);
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

//a conflict - bifurcation
boolean db_node_has_precisely_two_edges(dBNode * node, Orientation orientation, Nucleotide * nucleotide1, Nucleotide * nucleotide2, EdgeArrayType type, int index){
  
  Nucleotide n;

  Edges edges = get_edge_copy(*node,type,index);

  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      if (edges_count == 0){
	*nucleotide1 = n;
      }
      
      if (edges_count == 1){
	*nucleotide2 = n;
      }
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 2);
  
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


boolean db_node_check_status_not_pruned_or_visited(dBNode * node)
{
  if ( db_node_check_status(node,visited) || db_node_check_status(node,visited_and_exists_in_reference) || 
	 db_node_check_status(node, pruned_from_NA12878) ||  db_node_check_status(node, pruned_from_NA12891) || db_node_check_status(node, pruned_from_NA12892) ||
	 db_node_check_status(node, pruned_from_NA12878_and_NA12891) || db_node_check_status(node, pruned_from_NA12878_and_NA12892) || db_node_check_status(node, pruned_from_NA12891_and_NA12892)
	 ) 
    {
      return false;
    }
  else
    {
      return true;
    }
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
	  //printf("WARNING. Pruning a node that is already pruned");
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
	  //printf("WARNING. Pruning a node that is already pruned");
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
	  //printf("WARNING. Pruning a node that is already pruned");
	}
    }
  else if (db_node_check_status(node, pruned_from_NA12878_and_NA12891))
    {
      if (index==0)//NA12878
	{
	  //printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==1) //NA12891
	{
	  //printf("WARNING. Pruning a node that is already pruned");
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
	  //printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==1) //NA12891
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891_and_NA12892);
	}
      else if (index==2)//NA12892
	{
	  //printf("WARNING. Pruning a node that is already pruned");
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
	  //printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==2)//NA12892
	{
	  //printf("WARNING. Pruning a node that is already pruned");
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


void db_node_print_binary(FILE * fp, dBNode * node)
{

  BinaryKmer kmer = element_get_kmer(node);
  short covg[NUMBER_OF_INDIVIDUALS_PER_POPULATION];
  Edges individual_edges[NUMBER_OF_INDIVIDUALS_PER_POPULATION]; 

  int i;
  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)                                                     
    {                                                                                                         
      covg[i] = db_node_get_coverage_as_short(node, individual_edge_array, i);
      individual_edges[i]= get_edge_copy(*node, individual_edge_array, i);
    }      

  fwrite(&kmer,  sizeof(BinaryKmer), 1, fp);
  fwrite(covg, sizeof(short), NUMBER_OF_INDIVIDUALS_PER_POPULATION, fp); 
  fwrite(individual_edges, sizeof(Edges), NUMBER_OF_INDIVIDUALS_PER_POPULATION, fp);

  
}


boolean db_node_read_sv_trio_binary(FILE * fp, short kmer_size, dBNode * node){

  BinaryKmer kmer = element_get_kmer(node);
  short covg[NUMBER_OF_INDIVIDUALS_PER_POPULATION];
  Edges individual_edges[NUMBER_OF_INDIVIDUALS_PER_POPULATION]; 

  int read;
  
  read = fread(&kmer,sizeof(BinaryKmer),1,fp);

  if (read>0){


    read = fread(covg, sizeof(short), NUMBER_OF_INDIVIDUALS_PER_POPULATION, fp);    
    if (read==0){
      puts("error with input file\n");
      exit(1);
    }

    read = fread(individual_edges, sizeof(Edges), NUMBER_OF_INDIVIDUALS_PER_POPULATION, fp);
    if (read==0){
      puts("error with input file\n");
      exit(1);
    }



  }
  else{
    return false;
  }

  element_initialise(node,kmer,kmer_size);

  int i;
  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      node->coverage[i]         = covg[i];
      node->individual_edges[i] = individual_edges[i];
    }

  return true;
}

//read a binary for an individual person, as dumped by the target "graph"
// the edge array type and index tell you which person you should load this data into
boolean db_node_read_graph_binary(FILE * fp, short kmer_size, dBNode * node, EdgeArrayType type, int index)
{

  if ( (index<0) || (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION))
    {
      printf("Invalid index for which person to load binary into: %d. Exiting.", index);
      exit(1);
    }

  BinaryKmer kmer;
  Edges edges;
  short coverage;
  int read;
  
  read = fread(&kmer,sizeof(BinaryKmer),1,fp);

  if (read>0){
    read = fread(&coverage,sizeof(short),1,fp);    
    if (read==0){
      puts("error with input file\n");
      exit(1);
    }
    read = fread(&edges,sizeof(Edges),1,fp);
    if (read==0){
      puts("error with input file\n");
      exit(1);
    }   
  }
  else{
    return false;
  }

  element_initialise(node,kmer,kmer_size);
  node->individual_edges[index]    = edges;
  node->coverage[index] = coverage;
  return true;
}


void db_node_action_set_status_none(dBNode * node){
  db_node_set_status(node,none);
}

//void db_node_action_set_status_pruned(dBNode * node){
//  db_node_set_status(node,pruned);
//}


void db_node_action_set_status_visited(dBNode * node){
  db_node_set_status(node,visited);
}

void db_node_action_set_status_visited_or_visited_and_exists_in_reference(dBNode * node){

  if (db_node_check_status(node, exists_in_reference))
    {
      db_node_set_status(node,visited_and_exists_in_reference);      
    }
  //WARNING. Need special case for pruned?
  else if (!db_node_check_status(node, unassigned)) 
    {
      db_node_set_status(node,visited);
    }

}


void db_node_action_unset_status_visited_or_visited_and_exists_in_reference(dBNode * node){
  if (db_node_check_status_visited_and_exists_in_reference(node))
  {
    db_node_set_status(node,exists_in_reference);
  }
  else if (db_node_check_status(node, visited))
  {
    db_node_set_status(node, none);
  }
      
}



void db_node_action_do_nothing(dBNode * node){
  
}


boolean db_node_check_status_none(dBNode * node){
  return db_node_check_status(node,none);
}

boolean db_node_check_status_visited(dBNode * node){
  return db_node_check_status(node,visited);
}

boolean db_node_check_status_exists_in_reference(dBNode * node){
  return db_node_check_status(node,exists_in_reference);
}

boolean db_node_check_status_visited_and_exists_in_reference(dBNode * node){
  return db_node_check_status(node,visited_and_exists_in_reference);
}

boolean db_node_check_status_is_not_exists_in_reference(dBNode * node){

  if (db_node_check_status(node,exists_in_reference)==true )
    {
      return false;
    }
  else
    {
      return true;
    }
}

boolean db_node_check_status_is_not_visited(dBNode* node)
{
   if (db_node_check_status(node,visited)  )
    {
      return false;
    }
  else
    {
      return true;
    }
  
}

boolean db_node_check_status_is_not_visited_or_visited_and_exists_in_reference(dBNode * node){
  if (db_node_check_status(node,visited_and_exists_in_reference) || db_node_check_status(node,visited)  )
    {
      return false;
    }
  else
    {
      return true;
    }
}


boolean db_node_condition_always_true(dBNode* node)
{
  return true;
}

