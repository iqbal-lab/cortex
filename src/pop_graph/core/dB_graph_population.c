/*
  dB_graph_population.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>

#include <element.h>
#include <open_hash/hash_table.h>

#include <dB_graph.h>
#include <dB_graph_population.h>
#include <pqueue_pop.h>
#include <seq.h>
#include <string.h>

//it doesn't check that it is a valid arrow
dBNode * db_graph_get_next_node_for_specific_person_or_pop(dBNode * current_node, Orientation current_orientation, 
							   Orientation * next_orientation,
							   Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, EdgeArrayType type, int index){
  
  BinaryKmer kmer = element_get_kmer(current_node);
  dBNode * next_node;
  BinaryKmer rev_kmer = binary_kmer_reverse_complement(kmer,db_graph->kmer_size);
  
  if (current_orientation == reverse){   
    *reverse_edge = binary_kmer_get_last_nucleotide(kmer);
    kmer = rev_kmer;
  }
  else{
    *reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
  }

  
  kmer = binary_kmer_add_nucleotide_shift(kmer,edge, db_graph->kmer_size);

   //get node from table
  next_node = hash_table_find(element_get_key(kmer,db_graph->kmer_size),db_graph);

 

  if (next_node != NULL){
    *next_orientation = db_node_get_orientation(kmer,next_node,db_graph->kmer_size);
  }

  return next_node;
}


// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)


int db_graph_get_perfect_path_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit, void (*node_action)(dBNode * node),
			      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
			      boolean * is_cycle, dBGraph * db_graph, EdgeArrayType type, int index){

  Orientation  current_orientation,next_orientation;
  dBNode * current_node;
  dBNode * next_node;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  
  //sanity checks
  if (node == NULL){
    printf("db_graph_get_perfect_path_for_specific_person_or_pop: can't pass a null node\n");
    exit(1);
  }
  else if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      printf("db_graph_get_perfect_path_for_specific_person_or_pop: this node does not exist in this persons graph");
      exit(1);
    }



  current_node        = node;
  current_orientation = orientation;

  *is_cycle = false;
  
  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  

  if (DEBUG){
    printf("\nNode %i in path: %s\n", length, binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
   }

  

  if (db_node_has_precisely_one_edge(node, orientation,&nucleotide, type, index)){ //first node special case
    
    do{ 
      if (length>0){
	node_action(current_node);
      }

      next_node =  db_graph_get_next_node_for_specific_person_or_pop(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph, type, index); 
      
      //sanity check
      if(next_node == NULL){
	fprintf(stderr,"dB_graph_get_perfect_path_for_specific_person_or_pop: didnt find next node in hash table after this one: %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
	exit(1);
      }

      path_labels[length]        = nucleotide;
      length++;
      path_nodes[length]         = next_node;
      path_orientations[length]  = next_orientation;

      if (DEBUG){
	printf("\nNode %i in path: %s\n", length, binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq));
      }

      current_node        = next_node;
      current_orientation = next_orientation;
    
    } while (length<limit && 
	   !((next_node == node) && (next_orientation == orientation)) && //loop
	   db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2, type, index) && //multiple entries
	   db_node_has_precisely_one_edge(current_node, current_orientation,&nucleotide, type, index)); //has one next edge only


    if ((next_node == node) && (next_orientation == orientation)){
      *is_cycle = true;
    }
  }
  
  if (DEBUG){
    printf("\nLast node in path: %s %i length: %i\n", binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq),get_edge_copy(*next_node,type, index),length);
  }
  
  return length;
  
}

//a perfect bubble starts in a node with only two outgoing edges in the same orientation
// every branch of the bubble is free of conflicts until it joins with the main branch
// the length of every branch is the same -- smaller than limit!

boolean db_graph_detect_perfect_bubble_for_specific_person_or_pop(dBNode * node,
				       Orientation * orientation, 
				       boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node), 
				       Nucleotide * base1, Nucleotide * base2, 
				       Nucleotide * labels,dBNode * * end_node,Orientation * end_orientation,
				       dBGraph * db_graph, EdgeArrayType type, int index){

  

  dBNode * next_node1;
  dBNode * next_node2;
  dBNode * nodes1[db_graph->kmer_size];
  dBNode * nodes2[db_graph->kmer_size];
  Orientation orientations1[db_graph->kmer_size], orientations2[db_graph->kmer_size];
  Nucleotide labels1[db_graph->kmer_size],labels2[db_graph->kmer_size];

  Orientation next_orientation1, next_orientation2;
  Nucleotide rev_base1, rev_base2;
  boolean is_cycle1, is_cycle2;
  int length1, length2;
  boolean ret = false;
  *end_node = NULL;

  if (condition(node)){

    if (db_node_has_precisely_two_edges(node,forward,base1,base2, type, index)){
    
      *orientation = forward;
      next_node1 = db_graph_get_next_node_for_specific_person_or_pop(node,forward,&next_orientation1,*base1,&rev_base1,db_graph, type, index);
      next_node2 = db_graph_get_next_node_for_specific_person_or_pop(node,forward,&next_orientation2,*base2,&rev_base2,db_graph, type, index);
      
      if (next_node1 == NULL || next_node2 == NULL){
	puts("error!");
	exit(1);
      }

      length1  = db_graph_get_perfect_path_for_specific_person_or_pop(next_node1,next_orientation1,
					   db_graph->kmer_size,node_action,
					   nodes1,orientations1,labels1,&is_cycle1,db_graph, type, index);

      length2  = db_graph_get_perfect_path_for_specific_person_or_pop(next_node2,next_orientation2,
					   db_graph->kmer_size,node_action,
					   nodes2,orientations2,labels2,&is_cycle2,db_graph, type, index);

    }
    else{
      if (db_node_has_precisely_two_edges(node,reverse,base1,base2, type, index)){

	*orientation = reverse;
	next_node1 = db_graph_get_next_node_for_specific_person_or_pop(node,reverse,&next_orientation1,*base1,&rev_base1,db_graph, type, index);
	next_node2 = db_graph_get_next_node_for_specific_person_or_pop(node,reverse,&next_orientation2,*base2,&rev_base2,db_graph, type, index);
	
	if (next_node1 == NULL || next_node2 == NULL){
	  puts("error!");
	  exit(1);
	}
	
	length1 = db_graph_get_perfect_path_for_specific_person_or_pop(next_node1,next_orientation1,db_graph->kmer_size,node_action,nodes1,orientations1,labels1,&is_cycle1,db_graph, type, index);
	length2 = db_graph_get_perfect_path_for_specific_person_or_pop(next_node2,next_orientation2,db_graph->kmer_size,node_action,nodes2,orientations2,labels2,&is_cycle2,db_graph, type, index);

      }
    }
  
    node_action(node);

    ret = (length1 == db_graph->kmer_size && 
	   length2 == db_graph->kmer_size &&
	   nodes1[length1] == nodes2[length2]);
    
    if (ret == true){
      *end_node        = nodes1[length1];
      *end_orientation = orientations1[length1];
      int i;
      for(i=0;i<db_graph->kmer_size;i++){
	labels[i] = labels1[i];
      }
    }
  }

  //sanity check
  if (ret == true && *end_node == NULL){
    printf("db_graph_detect_perfect_bubble_for_specific person or pop: error inconsistent state\n");
    exit(1);
  }
  return ret;
}




// limit is the max number of nodes in the supernode. But remember that the first one corresponds to kmer_size bases, while each subsequent
// one corresponds to an extra base. Therefore:
//string has to support limit+db_graph->kmer_size+1 (+1 as you need a space for the \0 at the end)

int db_graph_supernode_for_specific_person_or_pop(dBNode * node,int limit, boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node),
		       char * string,dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
		       dBGraph * db_graph, EdgeArrayType type, int index){

  if (node==NULL)
    {
      printf("do not call db_graph_supernode_for_specific_person_or_pop with NULL node\n");
      exit(1);
    }

  dBNode * nodes_reverse[limit];
  Orientation orientation_reverse[limit];
  Nucleotide labels_reverse[limit];
  boolean is_cycle;
  int length_reverse;
  int length = 0;

  char tmp_seq[db_graph->kmer_size+1];


  if (condition(node)==true){   
        
    //compute the reverse path until the end of the supernode
    //return is_cycle_reverse == true if the path closes a loop    
 
    length_reverse = db_graph_get_perfect_path_for_specific_person_or_pop(node,reverse,limit,db_node_action_do_nothing,
					       nodes_reverse,orientation_reverse,labels_reverse,
					       &is_cycle,db_graph, type, index);

    //apply action to the first node of supernode
    node_action(nodes_reverse[length_reverse]);

    //we are at the end of a supernode
    length = db_graph_get_perfect_path_for_specific_person_or_pop(nodes_reverse[length_reverse],opposite_orientation(orientation_reverse[length_reverse]),limit,node_action,
				       path_nodes,path_orientations,path_labels,
				       &is_cycle,db_graph, type, index);
       

    //apply action to the last node
    node_action(path_nodes[length]);

    //get string
    char tmp_seq[db_graph->kmer_size+1];
    BinaryKmer kmer = element_get_kmer(path_nodes[0]);
    
    if (path_orientations[0]==reverse){
      kmer = binary_kmer_reverse_complement(kmer, db_graph->kmer_size);
    }     
    binary_kmer_to_seq(kmer,db_graph->kmer_size,tmp_seq);
    
    //build the string
    
    //the first node kmer
    int i;
    for(i=0;i<db_graph->kmer_size;i++){
      string[i] = tmp_seq[i];
    }
    
    //the path
    int j;
    for(j=0;j<length;j++){
      string[i] = binary_nucleotide_to_char(path_labels[j]);
      i++;
    }
    string[i] = '\0';
  }
  
  else
    {
      return 0;
    }
  
  return db_graph->kmer_size + length;
}




boolean db_node_is_supernode_end(dBNode * element,Orientation orientation, EdgeArrayType edge_type, int edge_index, dBGraph* db_graph)
{
  if (element==NULL)
    {
      //printf("Sending null pointer to db_node_is_supernode_end\n");
      return false;
    }
  char edges = get_edge_copy(*element, edge_type, edge_index);
  
  if (orientation == reverse)
    {
      //shift along so the 4 most significant bits become the 4 least - we've passed out argument by copy so not altering the original
      edges >>= 4;
    }

  edges &= 15; // AND with 00001111 so that we only look at the 4 least significant bits 
   
  //is supernode end EITHER if it has 0 edges out
  if (edges == 0)
    {
      return true;
    }
  else if ((edges != 1) && (edges != 2) && (edges != 4) && (edges != 8))  // or if it has too many edges, ie >1
    {
      return true;
    }
  else  //or next node has more than one arrow in
    {
      Orientation next_orientation;
      
      //we know this element has only one edge out. What is it? The function db_node_has_precisely_one_edge fills the answer into argument 3
      Nucleotide nuc, reverse_nuc;
      if (db_node_has_precisely_one_edge(element, orientation, &nuc, edge_type, edge_index))
	{

	  dBNode* next_node = db_graph_get_next_node_for_specific_person_or_pop(element, orientation, &next_orientation, nuc, &reverse_nuc,db_graph, edge_type, edge_index);

	  if ( (next_node==element) && (next_orientation==orientation) )
	    {
	      return true; //infinite self-loop
	    }
	  else if (db_node_has_precisely_one_edge(next_node, opposite_orientation(next_orientation), &nuc, edge_type, edge_index))
	    {
	      return false;//successor node also has only one edge coming in
	    }
	  else //successor has multiple edges
	    {
	      return true;      
	    }
	}
    }

  //to keep compiler happy
  printf("We have got the end of of the is_supernode_end function - should not reach here");
  exit(1);
  return true;
}

/*
dBNode * db_node_get_next_node_for_specific_person_or_pop(dBNode * current_node, Orientation current_orientation, 
			       Orientation * next_orientation,
			       Nucleotide edge, Nucleotide * reverse_edge,HashTable * db_graph, EdgeArrayType type, int index){
  
  BinaryKmer kmer = element_get_kmer(current_node);
  dBNode * next_node;
  BinaryKmer rev_kmer = binary_kmer_reverse_complement(kmer,db_graph->kmer_size);
  
  if (current_orientation == reverse){   
    *reverse_edge = binary_kmer_get_last_nucleotide(kmer);
    kmer = rev_kmer;
  }
  else{
    *reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
  }


  //Mario - we are going to get confused by things like this. here edge variable is a nucleotide,
  //but elsewhere edge is th char that describes all the edges in a node.
  kmer = binary_kmer_add_nucleotide_shift(kmer,edge, db_graph->kmer_size);

  //get node from table 
  next_node =  db_graph_find_node_restricted_to_specific_person_or_population(element_get_key(kmer,db_graph->kmer_size),db_graph,type,index);

  if (next_node != NULL){
    *next_orientation = db_node_get_orientation(kmer,next_node,db_graph->kmer_size);
  }

  return next_node;
}

*/

// wrapper for hash_table_find, which allows you to look in the hash table
// specifically for nodes related to a specific person or population
// person or population  specified by which edge array type
// which person or pop specified by index

dBNode *  db_graph_find_node_restricted_to_specific_person_or_population(Key key, dBGraph * hash_table, EdgeArrayType type, int index)

{

  dBNode *  e = hash_table_find(key, hash_table);

  if (e==NULL)
    {
      return NULL;
    }

  //ASSUMING read length is strictly greater than kmer-length
  //then you should never see a kmer (node) which is unconnected to another (even if the other is itself)
  //ie to check if this kmer is seen in a given person/pop, it is enough to check if it has an edge

  Edges edge_for_this_person_or_pop = get_edge_copy(*e, type, index);

  if (edge_for_this_person_or_pop == 0)
    {
      return NULL;
    }
  else
    {
      return e;
    }
  
}


//WARNING: this marks as visited any node that it walks over (and is in the supernode)
char * get_seq_from_elem_to_end_of_supernode_for_specific_person_or_pop(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle, char * seq, int max_length, 
									EdgeArrayType type, int index)
{
  Nucleotide nucleotide1, nucleotide2, rev_nucleotide;

  Orientation original_orientation, next_orientation;
  dBNode * original_node;
  dBNode * next_node;
  int seq_length = 0;
  char tmp_seq[db_graph->kmer_size];

  original_node = node;
  original_orientation = orientation; 
  

  //Mark the element we're starting at as visited
  db_node_set_status(node,visited);

  *is_cycle = false;

  
  while(db_node_has_precisely_one_edge(node,orientation,&nucleotide1, type, index)) {
 
    if ((seq_length+1)>max_length){
      printf("cannot allocate a sequence longer than the max length %d\n", max_length);
      exit(1);
    }

    next_node =  db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,nucleotide1,&rev_nucleotide,db_graph, type, index);
    
    if(next_node == NULL){
      printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size, tmp_seq));
      exit(1);
    }	         
    
    if (DEBUG){
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide1),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size, tmp_seq) :  binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size),db_graph->kmer_size, tmp_seq));
	     

    }
    
    //check for multiple entry edges 
    if (db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2, type, index)){
      seq[seq_length] =  binary_nucleotide_to_char(nucleotide1);
      seq_length++;
      if (DEBUG){
	printf("ADD %c\n",binary_nucleotide_to_char(nucleotide1));
      }
    }
    else{
      if (DEBUG){
	printf("Multiple entries\n");
      }
      break;
    }
    
    
    //loop
    if ((next_node == original_node) && (next_orientation == original_orientation)){      
      *is_cycle = true;
      //remove last addition
      seq_length--;
      break;
    }
      
    db_node_set_status(next_node,visited);
    
    node = next_node;
    orientation = next_orientation;      
  }

  seq[seq_length] = '\0';
  return seq;
}

//formerly void print_supernode_for_specific_person_or_pop(HashTable* db_graph, dBNode * node,EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test)
void db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop(HashTable* db_graph, dBNode * node, long* supernode_count, EdgeArrayType type, int index, 
									   boolean is_for_testing, char** for_test, int* index_for_test)
{

  FILE * fout;
  
  char filename [200];
  if (*supernode_count % 100000000 == 0)
    {
      int num = *supernode_count / 100000000;
      
      if (*supernode_count !=0)
	{
	  fclose(fout);
	}

      if (type == individual_edge_array)
	{
	  sprintf(filename,"out_nodes_kmer_%i_person_%i_subset_%i",db_graph->kmer_size,index,num);
	}
      else
	{
	  sprintf(filename,"out_nodes_kmer_%i_population_%i_subset_%i",db_graph->kmer_size,index,num);
	}
      //fprintf(stderr,"opening file %s\n",filename);
      fout = fopen(filename,"w");
    }
  *supernode_count = *supernode_count+1;
  db_graph_print_supernode_for_specific_person_or_pop(fout,node,db_graph, type,index, is_for_testing,  for_test,index_for_test);

}




void db_graph_choose_output_filename_and_print_potential_transloc_for_specific_person_or_pop(HashTable* db_graph, dBNode * node, long* supernode_count, EdgeArrayType type, int index,
											     int min_required_covg, int max_required_covg,
											     boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2)
{

  //**debug only zam
  //printf("call db_graph_choose_output_filename_and_print_potential_transloc_for_specific_person_or_pop\n");
  //char tmp_seq[db_graph->kmer_size];
  //printf("Check if this node %s has been visited\n",binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq));

  //if (db_node_check_status(node,visited))
  // {
  //    printf("Node is already visited\n");
  //  }
  //else
  // {
  //    printf("Node is not visited\n");
  //  }

      //end of debug only
  FILE * fout;
  
  char filename [200];
  if (*supernode_count % 100000000 == 0)
    {
      int num = *supernode_count / 100000000;
      
      if (*supernode_count !=0)
	{
	  fclose(fout);
	}

      if (type == individual_edge_array)
	{
	  sprintf(filename,"translocations_person_%i_mincovg_%i_maxcovg_%i_subset_%i",index, min_required_covg, max_required_covg, num);
	}
      else
	{
	  sprintf(filename,"translocations_population_%i_mincovg_%i_maxcovg_%i_subset_%i", index, min_required_covg, max_required_covg, num);
	}
      //fprintf(stderr,"opening file %s\n",filename);
      fout = fopen(filename,"w");
    }
  *supernode_count = *supernode_count+1;
  db_graph_print_supernode_if_is_potential_transloc_for_specific_person_or_pop(fout,node,db_graph, type,index, min_required_covg, max_required_covg,
									       is_for_testing,  for_test1, for_test2, index_for_test1, index_for_test2);

}


void db_graph_choose_output_filename_and_print_potential_inversion_for_specific_person_or_pop(HashTable* db_graph, dBNode * node, long* supernode_count, EdgeArrayType type, int index,
											      int min_required_covg, int max_required_covg,
											      boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2)
{

  FILE * fout;
  
  char filename [200];
  if (*supernode_count % 100000000 == 0)
    {
      int num = *supernode_count / 100000000;
      
      if (*supernode_count !=0)
	{
	  fclose(fout);
	}

      if (type == individual_edge_array)
	{
	  sprintf(filename,"inversions_person_%i_mincovg_%i_maxcovg_%i_subset_%i",index,min_required_covg, max_required_covg, num);
	}
      else
	{
	  sprintf(filename,"inversions_population_%i_mincovg_%i_maxcovg_%i_subset_%i",index, min_required_covg, max_required_covg, num);
	}
      //fprintf(stderr,"opening file %s\n",filename);
      fout = fopen(filename,"w");
    }
  *supernode_count = *supernode_count+1;
  db_graph_print_supernode_if_is_potential_inversion_for_specific_person_or_pop(fout,node,db_graph, type,index, min_required_covg, max_required_covg,
										is_for_testing,  for_test1, for_test2, index_for_test1, index_for_test2);

}





// *********************************************
// functions applied to a person/pop's graph
// *********************************************


void db_graph_traverse_specific_person_or_pop_for_supernode_printing(void (*f)(HashTable*, Element *, long* , EdgeArrayType, int, boolean, char**, int*),HashTable * hash_table, long* supernode_count, 
								     EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test){


  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(hash_table, &hash_table->table[i], supernode_count, type,index, is_for_testing, for_test, index_for_test);
    }
  }


}

void db_graph_traverse_specific_person_or_pop_for_supernode_and_chromosome_overlap_printing(void (*f)(HashTable*, Element *, long* , EdgeArrayType, int, int, int, boolean, char**, char**, int*, int*),
											    HashTable * hash_table, long* supernode_count, EdgeArrayType type, int index, int min_covg, int max_covg, 
											     boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2){

  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(hash_table, &hash_table->table[i], supernode_count, type, index, min_covg, max_covg, is_for_testing, for_test1, for_test2, index_for_test1, index_for_test2  );
    }
  }



}



void db_graph_traverse_to_gather_statistics_about_people(void (*f)(HashTable*, Element *, int**, int),HashTable * hash_table, int** array, int num_people )
{

  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(hash_table, &hash_table->table[i], array, num_people);
    }
  }




}


//will print only nodes for specific person/population

// ******NOTE****** if you want to go on to print for other people, you need to set status to none for all the nodes again after going thrrough the whole
// graph printing supernodes - otherwise some will be left marked visited.

//if you pass in true for the boolean is_for_testing,  then instead of printing to the file, it adds the supernode to the array of char*'s 
//this enables the caller to sort the supernodes lexicographically, and test that the right supernodes are being found
//this is not scalable to a ful genome, but sufficient for detailed unit tests

void db_graph_print_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test ){

  if (!db_node_is_this_node_in_this_person_or_populations_graph(node, type, index))
    {
      return;
    }
  
  const int              max_expected_supernode_length=12000;
  dBNode *         nodes_path[max_expected_supernode_length];
  Orientation      orientations_path[max_expected_supernode_length];
  Nucleotide       labels_path[max_expected_supernode_length];
  char             seq[max_expected_supernode_length];
  int              j;

  //initialise
  for (j=0; j<max_expected_supernode_length; j++)
    {
      seq[j]='0';
    }


  //get the supernode.
  int length = db_graph_supernode_for_specific_person_or_pop(node,max_expected_supernode_length,&db_node_check_status_not_pruned_or_visited,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							     seq,nodes_path,orientations_path, labels_path,db_graph, type, index);
  if (length==0) //means the condition db_node_check_status_not_pruned_or_visited was FALSE, so don't want to print anything
    {
      return;
    }

  if (!is_for_testing) 
    {
      fprintf(file,">NODE\n%s\n",seq);
    }
  else
    {
      
      for_test[*index_for_test] = (char*) calloc(length,sizeof(char));
      if (for_test[*index_for_test]==NULL)
	{
	  printf("Unable to calloc for supernode");
	  exit(1);
	}

      int j;

      for (j=0; j<length; j++)
	{
	  for_test[*index_for_test][j]=seq[j];
	}
      for_test[*index_for_test][length]='\0';
      *index_for_test=*index_for_test+1;
      
	//Now make sure you are using the smaller of the sequence and its rev complement

	  //for the moment, only do this for short supernodes :-(
	  
	//	if (length_of_supernode<32)
	// {
	//   BinaryKmer snode_kmer = seq_to_binary_kmer(for_test[*index_for_test],length_of_supernode);
	//   BinaryKmer rev_snode_kmer =  binary_kmer_reverse_complement(snode_kmer, length_of_supernode);
	    
	//    if (rev_snode_kmer<snode_kmer)
	//     {
	//	binary_kmer_to_seq(rev_snode_kmer,length_of_supernode, for_test[*index_for_test]);
	//     }
	// }
	//TODO - fix this - I was trying to find reverse complement by using binary_kmer_reverse complement, and this assumes k<31, so can use long long. But supernodes can be much longer than 31 bases.
	// this is only an issue because I want to print out the smaller of supernode and rev_comp(supernde), so not critical.
	

    }
  
  
}

//TODO - move thiscode out into /sv_trio
//the final argument returns the number of chromosomes intersected, but only gives the right answer if the function returns true - ie no node intersected >1 chromosome
boolean db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome(dBNode* node, EdgeArrayType type, int index, dBGraph* dbgraph, int* total_number_of_different_chromosomes_intersected)
{
 
  *total_number_of_different_chromosomes_intersected=0;
  char tmp_seq[dbgraph->kmer_size];

  int chrom_list[24]; //keep track of which chromosomes have been hit by this supernode
  int i;
  for (i=0; i<24; i++)
    {
      chrom_list[i]=-1;
    }
  int chrom_ctr=0; //points to entry in array available for next chromosome we find

  //first check the node exists in this person's graph
  if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      if (node==NULL)
        {
          printf("Bloody node is null so of course cant get first node. Should not be calling this function with null nodes");
	  exit(1);
        }
      else
        {
          printf("In db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome: This person %d does not have this node %s\n", index, binary_kmer_to_seq(node->kmer, dbgraph->kmer_size, tmp_seq));
	  //the only case when this should happen is when a read contains k-mer length of bases, and then N's. So the node does not have any edges in that person, so we think it's not in their graph.
        }
      //exit(1);
      return false;
    }

  dBNode* first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, type, index, dbgraph);
  //  db_node_set_status(first_node, visited);

  int which_chromosome=-99;
  if (! db_node_has_at_most_one_intersecting_chromosome(first_node, &which_chromosome) )
    {
      return false;
    }

 
  if (which_chromosome>0)
    {
      chrom_list[0]=which_chromosome;
      chrom_ctr=1;
    }
  else if (which_chromosome == -99)
    {
      printf("programming error counting chrom overlaps");
      exit(1);
    }
  else if (which_chromosome ==0)
    {
      //  ignore. no overlapping chromosome
      
    }
  else
    {
      printf("Programming error. which_chromosome is %d", which_chromosome);
    }


  dBNode* current_node;
  dBNode* next_node;
  current_node=first_node;
  Orientation start_orientation, current_orientation, next_orientation;

  //work out which direction to leave supernode in. Function is_supernode_end will also return true if is an infinite self-loop
  if (db_node_is_supernode_end(first_node,forward, type,index, dbgraph))
    {
      if (db_node_is_supernode_end(first_node,reverse, type,index, dbgraph))
	{
	  //singleton
	  return true;
	}
      else
	{
	  start_orientation=reverse;
	}
    }
  else
    {
      start_orientation=forward;
    }

  current_orientation=start_orientation;
  


  //unfortunately, this means applying is_supernode_end twice altogether to the start_node. TODO - improve this 
  while (!db_node_is_supernode_end(current_node,current_orientation, type,index, dbgraph))
    {
      next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(current_node, current_orientation, &next_orientation, type, index, dbgraph);
      //      db_node_set_status(next_node, visited);
      if (! db_node_has_at_most_one_intersecting_chromosome(next_node, &which_chromosome) )
	{
	  return false;
	}
      if ((next_node==first_node) && (next_orientation==start_orientation))//back to the start - will loop forever if not careful ;-0
	{
	  break;
	}

      boolean new_chromosome=false;

      if (which_chromosome>0) //will be zero if the node intersects no chromosomes
	{
	  new_chromosome=true;
	  
	  if (chrom_ctr>0)
	    {
	      for (i=0; i<chrom_ctr; i++)
		{
		  if (chrom_list[i]==0)
		    {
		      printf("Programming error. i is %d,  Should not have zero entries in this array until AFTER %d\n", i,chrom_ctr);
		      printf("Current values in chrom_list are \n");
		      int j;
		      for (j=0; j<=i ; j++)
			{
			  printf("%d ", chrom_list[j]);
			}
		      exit(1);
		}
		  else if (chrom_list[i]==which_chromosome)
		    {
		      //we have seen it already
		      new_chromosome=false;
		    }
		}
	    }
	}
	  
      if (new_chromosome)
	{
	  chrom_list[chrom_ctr]=which_chromosome;
	  chrom_ctr++;
	}
      current_node=next_node;
      current_orientation=next_orientation;
	
    }

  *total_number_of_different_chromosomes_intersected=chrom_ctr;
  //printf("\nTOTAL number of chromosomes intersected is %d\n", chrom_ctr);
  return true;


}

//unlike the previous function, this works even when some nodes intersect multiple chromosomes
//int db_graph_get_number_of_chromosomes_intersected_by_this_supernode(dBNode* node, EdgearrayType type, int index, dBGraph* dbgraph)
//{
  
//}


//assume we have checked and no node has >1 chrom intwersection
void db_graph_print_chrom_intersections_for_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, 
						      boolean is_for_testing, char** for_test, int* index_for_test)
{

  char tmp_seq[db_graph->kmer_size];

  //ignore if not in this persons graph
  if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      printf("In db_graph_print_chrom_intersections_for_supernode : This person %d does not have this node %s\n", index, binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq));
      return;
    }

  int len_chrom_string_for_test=0;


  dBNode* first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, type, index, db_graph);

  int which_chromosome=0;

  //if it overlaps >1 chromosomes, which_chromosome will be -1
  db_node_has_at_most_one_intersecting_chromosome(first_node, &which_chromosome);

  dBNode* current_node;
  dBNode* next_node;
  current_node=first_node;
  Orientation start_orientation, current_orientation, next_orientation;
  Overlap next_overlap;
  boolean is_singleton=false;
  
  //work out which direction to leave supernode in. Function is_supernode_end will also return true if is an infinite self-loop
  if (db_node_is_supernode_end(first_node,forward, type,index, db_graph))
    {
      if (db_node_is_supernode_end(first_node,reverse, type,index, db_graph))
	{
	  fprintf(file, "is singleton ");
	  is_singleton=true;
	  start_orientation=forward; //somewhat moot. The supernode is just this node.
	}
      else
	{
	  start_orientation=reverse;
	}
    }
  else
    {
      start_orientation=forward;
    }


  current_orientation=start_orientation;

  Overlap first_overlap=db_node_get_direction_through_node_in_which_chromosome_passes(first_node,which_chromosome);
  char tmp_char[2];
  overlap_to_char(first_overlap,tmp_char);

  char direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node[2];
  compare_chrom_overlap_and_supernode_direction(first_overlap, start_orientation, direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node);


  if (!is_for_testing)
    {
      fprintf(file, "%d%s ", which_chromosome, direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node);
    }
  else
    {

      for_test[*index_for_test] = (char*) calloc(1000,sizeof(char));//TODO - don't have hardcded 1000, use length of supernode
      if (for_test[*index_for_test]==NULL)
	{
	  printf("Unable to calloc for supernode");
	  exit(1);
	}
      for_test[*index_for_test][0]='\0';
      char chrom_as_string[3];
      sprintf(chrom_as_string, "%d%s ",which_chromosome, direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node);
      strcat(for_test[*index_for_test],chrom_as_string);
      len_chrom_string_for_test=len_chrom_string_for_test+strlen(chrom_as_string);
    }

  if (is_singleton)
    {
      
      if (!is_for_testing)
	{
	  fprintf(file, " singleton\n");
	}
      else
	{
	  for_test[*index_for_test][len_chrom_string_for_test]='\0';
	  *index_for_test=*index_for_test+1;
	  
	}
      return;
      
    }
  

  //OK - have done the first node - now do the rest in a loop.
  //unfortunately, this means applying is_supernode_end twice altogether to the start_node. TODO - improve this 
  while (!db_node_is_supernode_end(current_node,current_orientation, type,index, db_graph))
    {
      next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(current_node, current_orientation, &next_orientation, type, index, db_graph);
      db_node_has_at_most_one_intersecting_chromosome(next_node, &which_chromosome);

      next_overlap=db_node_get_direction_through_node_in_which_chromosome_passes(next_node,which_chromosome);
      if (which_chromosome==-1)
	{
	  printf("Called db_graph_print_chrom_intersections_for_supernode_for_specific_person_or_pop with a node that intersects >1 chromosome. Should not do that\n");
	  exit(1);
	}

      overlap_to_char(next_overlap,tmp_char);
      compare_chrom_overlap_and_supernode_direction(next_overlap, next_orientation, direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node);
      if ((next_node==first_node) && (next_orientation==start_orientation))//back to the start - will loop forever if not careful ;-0
	{
	  if (!is_for_testing)
	    {
	      fprintf(file, " Looped back to start node");
	    }

	  break;
	}

      if (!is_for_testing)
	{
	  if (which_chromosome==0)
	    {
	      //fprintf(file, "00 ");
	      if (strcmp("0", direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node))
		{
		  printf("prg error. Since which_chromosome=0 there should be no intersection.");
		  printf(" but direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node is not zero!");
		  exit(1);
		}
	    }
	  //else
	  //  {
	  fprintf(file, "%d%s ", which_chromosome, direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node);
	  //  }

	}
      else
	{
	  char chrom_as_string2[4];
	  sprintf(chrom_as_string2, "%d%s ",which_chromosome, direction_of_chromosome_passing_through_node_compared_with_direction_of_supernode_passing_through_node);
	  strcat(for_test[*index_for_test],chrom_as_string2);
	  len_chrom_string_for_test=len_chrom_string_for_test+strlen(chrom_as_string2);
	}
      current_node=next_node;
      current_orientation=next_orientation;
	
    }


  if (!is_for_testing)
    {
      //ok - so current_node is now the end of the supernode.
      fprintf(file, "\n");
    }
  else
    {
      for_test[*index_for_test][len_chrom_string_for_test]='\0';
      *index_for_test=*index_for_test+1;

    }
  return;
}




//TODO - implement db_graph_apply_function_to_all_nodes-in_supernode_stopping_as_soon_as_one_returns_false


void db_graph_print_supernode_if_is_potential_transloc_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, 
										  int min_required_covg, int max_required_covg, 
										  boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2 )
{
  
  // ignore if visited
  if (db_node_check_status(node, visited))
    {
      //char tmp_seq[db_graph->kmer_size];
      //printf("ignoring visited node %s\n",binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq));
      return;
    }

  //ignore if not in this persons graph
  if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      return;
    }

  //ignore if pruned from this person's graph - TODO

  int total_number_of_different_chromosomes_intersected=0;

  //ignore unless for all nodes in supernode, has <=1 chrom intersection
  if  ( db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome(node, type, index, db_graph, &total_number_of_different_chromosomes_intersected))
    {

      int min=0;
      int max=0;
      db_graph_get_min_and_max_covg_of_nodes_in_supernode_for_specific_person_or_pop(node, type, index, db_graph,&min, &max);
      if ( (min>=min_required_covg) && (max<= max_required_covg))
        {
	  
	  if (total_number_of_different_chromosomes_intersected==2)
	    {
	      //then this is a supernode which is a potential sv locus.
	      
	      //first print out the supernode itself
	      db_graph_print_supernode_for_specific_person_or_pop(file, node, db_graph, type, index, is_for_testing, for_test1, index_for_test1 );
	      //then print out the chromosome intersections
	      db_graph_print_chrom_intersections_for_supernode_for_specific_person_or_pop(file, node, db_graph, type, index, is_for_testing, for_test2,  index_for_test2);
	    }
	}
    }


  
  //now mark all nodes in supernode as visited
  db_graph_set_status_of_all_nodes_in_supernode(node, visited, type, index, db_graph);
}



void db_graph_print_supernode_if_is_potential_inversion_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, 
										   int min_required_covg, int max_required_covg,
										   boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2 )
{
  
  // ignore if visited 
  if (db_node_check_status(node, visited))
    {
      //char tmp_seq[db_graph->kmer_size];
      //printf("ignoring visited node %s\n",binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq));
      return;
    }

  //ignore if not in this persons graph
  if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      return;
    }

  //ignore if pruned from this person's graph - TODO

  int total_number_of_different_chromosomes_intersected=0;

  //ignore unless for all nodes in supernode, has <=1 chrom intersection
  if  ( db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome(node, type, index, db_graph, &total_number_of_different_chromosomes_intersected))
    {

      int min=0;
      int max=0;
      db_graph_get_min_and_max_covg_of_nodes_in_supernode_for_specific_person_or_pop(node, type, index, db_graph,&min, &max);
      if ( (min>=min_required_covg) && (max<= max_required_covg))
	{
	  
	  
	  if (total_number_of_different_chromosomes_intersected==1)
	    {
	      //then this is a supernode which is a potential inversion.
	      
	      //first print out the supernode itself
	      db_graph_print_supernode_for_specific_person_or_pop(file, node, db_graph, type, index, is_for_testing, for_test1, index_for_test1 );
	      //then print out the chromosome intersections
	      db_graph_print_chrom_intersections_for_supernode_for_specific_person_or_pop(file, node, db_graph, type, index, is_for_testing, for_test2,  index_for_test2);
	    }
	}
    

    }

  //now mark all nodes in supernode as visited
  db_graph_set_status_of_all_nodes_in_supernode(node, visited, type, index, db_graph);
  
}


void db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population(dBGraph* hash_table, EdgeArrayType type, int index)
{
  printf("not implemented yet");
  exit(1);
}

void db_graph_set_status_of_all_nodes_in_supernode(dBNode* node, NodeStatus status, EdgeArrayType type, int index,  dBGraph* dbgraph)
{

  dBNode* first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, type, index, dbgraph);
  dBNode* current_node=first_node;

  dBNode* next_node;
  Orientation current_orientation, next_orientation, start_orientation;
  

  db_node_set_status(current_node, status);


  //work out which direction to leave supernode in. Function is_supernode_end will also return true if is an infinite self-loop
  if (db_node_is_supernode_end(current_node,forward, type,index, dbgraph))
    {
      if (db_node_is_supernode_end(current_node,reverse, type,index, dbgraph))
	{
	  //singleton
	  return ;
	}
      else
	{
	  start_orientation=reverse;
	}
    }
  else
    {
      start_orientation=forward;
    }

  current_orientation=start_orientation;

  //unfortunately, this means applying is_supernode_end twice altogether to the start_node. TODO - improve this 
  while (!db_node_is_supernode_end(current_node,current_orientation, type,index, dbgraph))
    {
      next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(current_node, current_orientation, &next_orientation, type, index, dbgraph);

      if ((next_node==first_node) && (next_orientation==start_orientation))//back to the start - will loop forever if not careful ;-0
	{
	  break;
	}
      db_node_set_status(next_node, status);

      current_node=next_node;
      current_orientation=next_orientation;
	
    }

  return;


}


void db_graph_get_min_and_max_covg_of_nodes_in_supernode_for_specific_person_or_pop(dBNode* node, /* NodeStatus status,*/  EdgeArrayType type, int index,  dBGraph* dbgraph, int* min_covg, int* max_covg)
{

  if (node==NULL)
    {
      printf("Do not call get_min_and_max_covg on a NULL node");
      exit(1);
    }

  
  dBNode* first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, type, index, dbgraph);
  dBNode* current_node=first_node;

  int min = db_node_get_coverage(first_node, type, index);
  int max = db_node_get_coverage(first_node, type, index);

  dBNode* next_node;
  Orientation current_orientation, next_orientation, start_orientation;
  



  //work out which direction to leave supernode in. Function is_supernode_end will also return true if is an infinite self-loop
  if (db_node_is_supernode_end(current_node,forward, type,index, dbgraph))
    {
      if (db_node_is_supernode_end(current_node,reverse, type,index, dbgraph))
	{
	  //singleton
	  *min_covg=min;
	  *max_covg=max;
	  return ;
	}
      else
	{
	  start_orientation=reverse;
	}
    }
  else
    {
      start_orientation=forward;
    }

  current_orientation=start_orientation;

  //unfortunately, this means applying is_supernode_end twice altogether to the start_node. TODO - improve this 
  while (!db_node_is_supernode_end(current_node,current_orientation, type,index, dbgraph))
    {
      next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(current_node, current_orientation, &next_orientation, type, index, dbgraph);

      if ((next_node==first_node) && (next_orientation==start_orientation))//back to the start - will loop forever if not careful ;-0
	{
	  break;
	}

      if (db_node_get_coverage(next_node,type,index)>max)
	{
	  max = db_node_get_coverage(next_node,type,index);
	}
      if (db_node_get_coverage(next_node,type,index)<min)
	{
	  min = db_node_get_coverage(next_node,type,index);
	}
      current_node=next_node;
      current_orientation=next_orientation;
	
    }

  *min_covg=min;
  *max_covg=max;
  return;


}










dBNode* db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(dBNode* node, EdgeArrayType type, int index, dBGraph* db_graph)
{
  
  char tmp_seq[db_graph->kmer_size];

  if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      if (node==NULL)
	{
	  //printf("Bloody node is null so of course cant get first node");
	}
      else
	{
	  //printf("This person %d does not have this node %s\n", index, binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq));
	}
      return NULL;
    }

  if (!db_node_check_status_not_pruned(node)  )
    {
      //don't waste time with pruned nodes.
      return NULL;
    }
  

  boolean is_cycle;
  Nucleotide nucleotide1, nucleotide2, rev_nucleotide;
  Orientation original_orientation, next_orientation, orientation;
  dBNode * original_node=node;
  dBNode * next_node;


  //First node in supernode is, almost by definition, what you get if you go in the Reverse direction (with respect to the Node)
  // as far as you can go.
  original_orientation = reverse; 
  orientation = reverse;
  is_cycle = false;


  while(db_node_has_precisely_one_edge(node,orientation,&nucleotide1, type, index)) {
 

    next_node =  db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,nucleotide1,&rev_nucleotide,db_graph, type, index);
    
    if(next_node == NULL){
      printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size, tmp_seq));
      exit(1);
    }	         

    if (DEBUG)
      {
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide1),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size, tmp_seq) :  
	     binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size),db_graph->kmer_size, tmp_seq));
	     
      }
    

    //check for multiple entry edges 
    if (db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2, type, index))
      {
      }
    else
      {
	if (DEBUG)
	  {
	    printf("Multiple entries\n");
	  }
	//break;
	 
	//printf("returning this first nodem, with kmer %s\n", binary_kmer_to_seq(node->kmer, db_graph->kmer_size));
      return node; //we have gone as far as we can go - the next node has multiple entries. So we are now at the first node of the supernode
      }
    
    
    //loop
    if ((next_node == original_node) && (next_orientation == original_orientation))
      {      

	is_cycle = true;
	//break;
	
	//printf("We have a loop, so original node will do, with kmer %s\n", binary_kmer_to_seq(original_node->kmer, db_graph->kmer_size, tmp_seq));
	return original_node; //we have a loop that returns to where we start. Might as well consider ourselves as at the fiurst node of the supernode right at the beginning
      }
    
    //printf("zam8\n");
    node = next_node;
    orientation = next_orientation;      
  }
  //printf("We have found the first node, it is %s\n", binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq));
  return node;

}


dBNode* db_graph_get_next_node_in_supernode_for_specific_person_or_pop(dBNode* node, Orientation orientation, Orientation* next_orientation, EdgeArrayType type, int index, dBGraph* db_graph)
{
  char tmp_seq[db_graph->kmer_size];

  if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      //printf("\nThis node is not in the graph of this person\n");
      return NULL;
    }
  else if (!db_node_check_status_not_pruned(node))
    {
      printf("ignore pruned node");
      //don't waste time with pruned nodes.
      return NULL;
    }
  else if (db_node_is_supernode_end(node, orientation, type, index, db_graph))
    {
      if (DEBUG)
	{
	  printf("this node is at the end of the supernode, in this orientation, so cant return the next one\n");
	}
      return NULL;
    }

  Nucleotide nucleotide_for_only_edge, reverse_nucleotide_for_only_edge;

  
  db_node_has_precisely_one_edge(node,orientation,&nucleotide_for_only_edge, type, index);//gives us nucleotide
  
  dBNode* next_node =  db_graph_get_next_node_for_specific_person_or_pop(node,orientation ,next_orientation,nucleotide_for_only_edge,&reverse_nucleotide_for_only_edge, db_graph, type, index);
  
  if(next_node == NULL){
    printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size, tmp_seq));
    exit(1);
  }	         

  if (DEBUG)
    {
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide_for_only_edge),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size, tmp_seq) :  
	     binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size),db_graph->kmer_size, tmp_seq));
      
    }
    

  //check for multiple entry edges 
  Nucleotide nuc;
  if (db_node_has_precisely_one_edge(next_node,opposite_orientation(*next_orientation),&nuc, type, index))
    {
    }
  else
    {
      //double check
      if (node ==NULL)
	{
	  printf("programming error. returning null node when my model in my head says impossible\n");
	  exit(1);
	}
      return node; //we have gone as far as we can go - the next node has multiple entries. So we are now at the first node of the supernode
    }
    
    
  //loop
  if ((next_node == node) && (*next_orientation == orientation))
    {      
      //double check
      if (node ==NULL)
	{
          printf("programming error. returning null node when my model in my head says impossible\n");
	  exit(1);
        }

	return node; //we have a kmer that loops back on itself
    }
  
  return next_node;

}


// Given a node, will look at each person's supernode containing that node. 
// For each of these supernodes, independently,  find the longest contiguous subsection that has min people-coverage > min_covg_for_pop_supernode.
// (Note this may not contain the original node.)
// Then compare across people, and define your consensus as the sub_supernode which has the most people backing it.
// If two people have equally good sub_supernodes, choose that of the person with highest index.
// Pass in pre-malloced Sequence* for answer, of length max_length_for_supernode (arg 2)

void  db_graph_find_population_consensus_supernode_based_on_given_node(Sequence* pop_consensus_supernode, int max_length_of_supernode, dBNode* node, int min_covg_for_pop_supernode, int min_length_for_pop_supernode, 
								       dBGraph* db_graph)
{
  int length_of_best_sub_supernode_in_each_person[NUMBER_OF_INDIVIDUALS_PER_POPULATION];
  int index_of_start_of_best_sub_supernode_in_each_person[NUMBER_OF_INDIVIDUALS_PER_POPULATION];

  int i;
  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      length_of_best_sub_supernode_in_each_person[i]=0;
      index_of_start_of_best_sub_supernode_in_each_person[i]=0;
    }

  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
        //get the first node in the supernode for this edge, within this person't graph
      dBNode* first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, individual_edge_array, i, db_graph);
      
      //note that node may be null, and also may not be in this person's graph. In this case, ignore.
      if (first_node ==NULL)
	{
	  //printf("first node for %d is nul, so ignore \n",i);
	  index_of_start_of_best_sub_supernode_in_each_person[i]=0;
	  length_of_best_sub_supernode_in_each_person[i]=0;
	}
      else
	{
	  //this returns 0 in length_of_best_sub_supernode_in_each_person if no subsupernode matches constraints
	  //printf("\nFind best sub supernode in person %d based on first node in supernode, which is %s\n", i, binary_kmer_to_seq(first_node->kmer, db_graph->kmer_size) );
	  db_graph_get_best_sub_supernode_given_min_covg_and_length_for_specific_person_or_pop(first_node,  &index_of_start_of_best_sub_supernode_in_each_person[i], &length_of_best_sub_supernode_in_each_person[i], 
											       min_covg_for_pop_supernode,  min_length_for_pop_supernode, individual_edge_array, i, db_graph); 
	  //printf("OK - that sub supernode in person %d starts at %d ends at %d\n", i, index_of_start_of_best_sub_supernode_in_each_person[i], length_of_best_sub_supernode_in_each_person[i]);
	}
    }

  
  //check that at least one person has a subsupernode that matches criteria
  boolean noone_has_decent_sub_supernode=true;
  for (i=0; (i< NUMBER_OF_INDIVIDUALS_PER_POPULATION) && (noone_has_decent_sub_supernode==true); i++)
    {
      if ( length_of_best_sub_supernode_in_each_person[i] >0 )
	{
	  noone_has_decent_sub_supernode=false;
	}
    }
  
  if (noone_has_decent_sub_supernode)
    {
      pop_consensus_supernode->seq[0]='\0';
      //printf("noone has a decent supernode\n");
      return;
    }


  //Now find which person has the bext sub_supernode nucleated at this node 
  int person_with_best_sub_supernode;
  int max=0;
  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      if ( length_of_best_sub_supernode_in_each_person[i] > max )
	{
	  person_with_best_sub_supernode=i;
	  max = length_of_best_sub_supernode_in_each_person[i];
	}
    }

  //printf("We think the person with best supernodenode is %d and length is %d\n",  person_with_best_sub_supernode, length_of_best_sub_supernode_in_each_person[person_with_best_sub_supernode]); 
  if (max==0)
    {
      printf("This should be impossible. Max size of sub supernode over all people is 0. Since we start wth a kmer that is in the graph, at least one person ought to have it\n");
      exit(1);
    }
  
  else
    {
      if (max > max_length_of_supernode)
	{
	  printf("pop_consensus_supernode not big enough - only alloced %d and we need %d\n", max_length_of_supernode, max);
	  exit(1);
	}


      int start = index_of_start_of_best_sub_supernode_in_each_person[person_with_best_sub_supernode];
      int end =    start+ length_of_best_sub_supernode_in_each_person[person_with_best_sub_supernode]-1;
      if ((start<0) || (end<0))
	{
	  printf("This is wrong. start is %d and end is %d for person %d and their length is %d\n", start,end, person_with_best_sub_supernode, length_of_best_sub_supernode_in_each_person[person_with_best_sub_supernode]);
	  exit(1);
	}      
      //      if (start==end)
      //	{
      //	  //then the best we have anaged is a single kmer. No need to bother getting that subsection
      //	  pop_consensus_supernode->seq = ???
      //	}
  
      else if (db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(pop_consensus_supernode->seq, node, start, end, individual_edge_array, person_with_best_sub_supernode, db_graph)==1)
	{
	  printf("Aborting - something wrong with gettig the subsection");
	  exit(1);
	}
      //printf("OKOK found bvest subsection for is %s\n", pop_consensus_supernode->seq);
    }

  

}

//subsection (argument1) is allocated by the caller
//returns 0 if successfully, 1 otherwise
int db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(char* subsection, dBNode* node, int start, int end, EdgeArrayType type, int index, dBGraph* db_graph)
{
  char tmp_seq[db_graph->kmer_size];

  //printf("CALL GET SUBSECTION With start %d and end %d\n\n", start, end);
  if ( (start<0) || (end<0) || (end-start<0) )
    {
      if (DEBUG)
	{
	  printf("bad args for getting subsection start %d and end %d", start, end);
	}
      subsection="ERROR";
      return 1;
    }

  
  dBNode* first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, type, index, db_graph);
  dBNode* current_node=first_node;
  dBNode* next_node;
  Orientation correct_direction_to_go;

  if (db_node_is_supernode_end(first_node,forward, type, index, db_graph))
    {
      correct_direction_to_go=reverse;
    }
  else if (db_node_is_supernode_end(first_node,reverse, type, index, db_graph))
    {
      correct_direction_to_go=forward; 
    }
  else
    {
      printf("Something has gone wrong. This node has been given as first node of a supernode, but db_node_is_supernode_end thinks it is not in either direction\n");
      return 1;
    }

  

 int i;
 Orientation current_orientation, next_orientation;
 current_orientation = correct_direction_to_go;

  for (i=0; i<start; i++)
    {
      next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(current_node, current_orientation, &next_orientation, type, index, db_graph);
      if (next_node==NULL)
	{
	  if (DEBUG)
	    {
	      printf("You're asking for the section from node %d to node %d in supernode that is not even that long", start, end);
	    }
	  return 1;
	}
      current_node=next_node;
      current_orientation=next_orientation;
    }

  char* first_kmer_in_subsection; 
  //first kmer depends on which direction you start in with respect to the first node OF THE SUBSECTION, not the first node of the supernode
  if (current_orientation == forward) //i.e. orientation of first node of subsection
    {
      first_kmer_in_subsection = binary_kmer_to_seq(current_node->kmer, db_graph->kmer_size, tmp_seq);
    }
  else
    {
      first_kmer_in_subsection = binary_kmer_to_seq(   binary_kmer_reverse_complement(current_node->kmer, db_graph->kmer_size) , db_graph->kmer_size, tmp_seq);
    }
  

   for(i=0; i<db_graph->kmer_size; i++)
   {
      subsection[i] = first_kmer_in_subsection[i];
   }


   for (i=db_graph->kmer_size; i<=db_graph->kmer_size + end-start-1; i++)
   {
      next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(current_node, current_orientation, &next_orientation, type, index, db_graph);
      if (next_node==NULL)
        {
	  if (DEBUG)
	    {
	      printf("Youre asking for the section from node %d to node %d in supernode that is not even that long", start, end);
	    }
          return 1;
        }
      Nucleotide nuc_for_next_edge;
      db_node_has_precisely_one_edge(current_node,current_orientation,&nuc_for_next_edge, type, index);
      subsection[i]=binary_nucleotide_to_char(nuc_for_next_edge);

      current_node=next_node;
      current_orientation=next_orientation;

      
   }
   subsection[db_graph->kmer_size + end-start]='\0';

   return 0;
}


// walk through the supernode, from one end to the other. As soon as we hit a node that has people-covg >= min_covg_for_pop_supernode
// then note down the coordinate (how far through the supernode), and keep going until we meet a node with people-covg < min_covg_for_pop_supernode
// Work out the length of the chunk you found, which had sufficient covg. If > min_length_for_pop_supernode, then you have found your first
// candidate for this person's best sub-supernode. Log the start-point and length. Then continue walkign through the supernode until you
// find another node that has people_covg sufficiently high. repeat the process for this point, and if its candidate is longer than our
// current candidate, then replace it with this new one. etc.

//Note that if none of the nodes match criteria, return 0 in length of best sub_supernode
void  db_graph_get_best_sub_supernode_given_min_covg_and_length_for_specific_person_or_pop(dBNode* first_node_in_supernode,  int* index_for_start_of_sub_supernode,
                                                                                           int* length_of_best_sub_supernode, int min_people_coverage,
                                                                                           int min_length_of_sub_supernode, EdgeArrayType type, int index, dBGraph* db_graph)
{

  char tmp_seq[db_graph->kmer_size];

  if (first_node_in_supernode==NULL)
    {
      printf("do not give null pointer to get_best_sub_supernode function");
      exit(1);
    }

  // OK - which direction do we walk?
  // Either - this node is a singleton supernode (either totally unconnected or too many edges in both directions)
  // Or    - we can go forward/reverse only, because the other direction has 0 or >1 edges.
  // Any other possibility breaks the condition that we are at a supernode end.

  Orientation correct_direction_to_go;

  if (db_node_is_supernode_end(first_node_in_supernode,forward, type, index, db_graph))
    {
    correct_direction_to_go=reverse;
    }
  else if (db_node_is_supernode_end(first_node_in_supernode,reverse, type,  index, db_graph))
    {
      correct_direction_to_go=forward; 
    }
  else
    {
      //must therefore be the case that this node is not in this person's graph. Check this
      if ( ! (db_node_is_this_node_in_this_person_or_populations_graph(first_node_in_supernode, type, index)) )
	{
	  *index_for_start_of_sub_supernode=0;
	  *length_of_best_sub_supernode=0;
	}
      else
	{
	  printf("Programming error. This node is supposed to be at the start of a supernode for this person, but is not in either direction. However it IS in their graph\n");
	}
    }


  
  int start_of_best_section_so_far=0;
  int length_of_best_section_so_far=0;
  boolean in_middle_of_a_good_section=false;
  int current_start=0;
  int current=0;

  dBNode* current_node=first_node_in_supernode;
  // printf("starting at first node in supernode: %s\n", binary_kmer_to_seq(current_node->kmer,db_graph->kmer_size, tmp_seq));
  dBNode* next_node;
  Orientation current_orientation = correct_direction_to_go;
  Orientation next_orientation;

  boolean reached_end=false;

  //does the start_node have enough coverage?

  if (element_get_number_of_people_or_pops_containing_this_element(current_node, type, index) < min_people_coverage)
    {
    }
  else
    {
      in_middle_of_a_good_section=true;      
    }
  while(!reached_end)
    {
      next_node=db_graph_get_next_node_in_supernode_for_specific_person_or_pop(current_node, current_orientation, &next_orientation, type, index, db_graph);
      if (next_node==NULL)
	{
	  //printf("Reached end of supernode\n");
	  reached_end=true;
	  continue;
	}
      else if (next_node==first_node_in_supernode) 
	{
	  //printf("reached end of supernode - in this case back at the start so it stops just before becoming a loop\n");
	  reached_end=true;
	  continue;
	}


      int next_people_cov= element_get_number_of_people_or_pops_containing_this_element(next_node, type, index);

      if  ( next_people_cov < min_people_coverage)
	{

	  current++;
	  current_node=next_node;
	  current_orientation=next_orientation;
	  current_start=current;
	  in_middle_of_a_good_section=false;

	  //rest of this if clause should be DEBUG only
	  char* next_kmer;
	  if (next_orientation==forward)
	    {
	      next_kmer = binary_kmer_to_seq(next_node->kmer,db_graph->kmer_size, tmp_seq);
	    }
	  else
	    {
	      next_kmer=binary_kmer_to_seq( binary_kmer_reverse_complement(next_node->kmer,db_graph->kmer_size), db_graph->kmer_size, tmp_seq );
	    }
	  
	  
	  //printf("Looking for best subsection Next node is %s\n", next_kmer);
	  //printf("Too little peope coverage on this node - only %d\n", next_people_cov);
	  
	 
	}
      
      else //there is a next node, and it has enough people coverage
	{

	  current++;
          current_node=next_node;
          current_orientation=next_orientation;

	  if (in_middle_of_a_good_section)
	    {
	      
	    }
	  else
	    {
	      current_start=current;
	      in_middle_of_a_good_section=true;
	    }

	  //printf("Looking for best subsection Next node is %s\n", binary_kmer_to_seq(next_node->kmer,db_graph->kmer_size, tmp_seq));
	
	  if (current-current_start+1 > length_of_best_section_so_far)
	    {
	      start_of_best_section_so_far=current_start;
	      length_of_best_section_so_far = current-current_start+1; //remember is length in nodes not bases
	    }
	  //printf("People covg is sufficient, at %d\nCurrent is %d, current_start is %d, best length sofar s %d\n", next_people_cov, current, current_start, length_of_best_section_so_far);
	  
	}
    }

  //length of best section is counted in supernodes. min length of sub supernode is in bases
  if ( ( length_of_best_section_so_far + db_graph->kmer_size - 1) >= min_length_of_sub_supernode)
    {
      //printf("Good. store best section start as %d, and length %d", start_of_best_section_so_far, length_of_best_section_so_far);
      *index_for_start_of_sub_supernode=start_of_best_section_so_far;
      *length_of_best_sub_supernode=length_of_best_section_so_far;
      return;
    }
  else
    {
      //printf("This person has nothing. Min length is %d but we only manage %d\n", min_length_of_sub_supernode, length_of_best_section_so_far);
      *index_for_start_of_sub_supernode=0;
      *length_of_best_sub_supernode=0;
      return;
    }
}



void print_node_to_file_according_to_how_many_people_share_it(HashTable* db_graph, dBNode * node, FILE** list_of_file_ptrs)
{
  int i;
  char tmp_seq[db_graph->kmer_size];
  
  int number_of_individuals_with_this_node=0;

  for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      if (get_edge_copy(*node, individual_edge_array, i) ==0 )
	{
	}
      else
	{
	  number_of_individuals_with_this_node++;
	}
    }

  char* kmer_as_string = binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq);
  fprintf(list_of_file_ptrs[number_of_individuals_with_this_node], "%s\n", kmer_as_string);
  
}

//array_of_counts totals up the number of kmers that are shared by 0,1,2,... individuals. Obviously the first element should be zero (no element
//in the graph is shared by no people) - use this as a crude check.
//number of people loaded is a param
void find_out_how_many_individuals_share_this_node_and_add_to_statistics(HashTable* db_graph, dBNode * node, int** array_of_counts, int number_of_people)
{

  char tmp_seq[db_graph->kmer_size];
  int i;

  if (number_of_people>NUMBER_OF_INDIVIDUALS_PER_POPULATION)
    {
      printf("Cannot call find_out_how_many_individuals_share_this_node_and_add_to_statistics with number_of_people = %d, as it's bigger than the NUMBER per pop, %d", number_of_people,NUMBER_OF_INDIVIDUALS_PER_POPULATION);
      exit(1);
    }
  int number_of_individuals_with_this_node=0;

  for (i=0; i<number_of_people; i++)
    {
      
      if (get_edge_copy(*node, individual_edge_array, i) ==0 )
	{
        }
      else
        {
          number_of_individuals_with_this_node++;
        }
    }

  //char* kmer_as_string = binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq);
  //printf("There are %d people with node %s\n", number_of_individuals_with_this_node,kmer_as_string);
  

  if (number_of_individuals_with_this_node>0)
    {
      (*(array_of_counts[number_of_individuals_with_this_node]))++;
    }
}


