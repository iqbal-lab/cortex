/*
  dB_graph_population.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>

#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
// #include <pqueue_pop.h>
#include <seq.h>
#include <string.h>
#include <limits.h>
#include <file_reader.h>

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

  //  char tmp_seq[db_graph->kmer_size+1];


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




// identical code to db_graph_supernode_for_specific_person_or_pop but this returns the index of the query node (first argument) in th supernode array
// will make a significant performance gain compared with getting the path_nodes array back and then searching it
int db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(dBNode * node,int limit, boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node),
									    char * string,dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, int* query_node_posn, 
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

  //  char tmp_seq[db_graph->kmer_size+1];


  if (condition(node)==true){   
        
    //compute the reverse path until the end of the supernode
    //return is_cycle_reverse == true if the path closes a loop    
 
    length_reverse = db_graph_get_perfect_path_for_specific_person_or_pop(node,reverse,limit,db_node_action_do_nothing,
					       nodes_reverse,orientation_reverse,labels_reverse,
					       &is_cycle,db_graph, type, index);

    //returns the index of node within the supernode array that we are creating
    *query_node_posn=length_reverse;

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


/*
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

*/



// *********************************************
// functions applied to a person/pop's graph
// *********************************************


void db_graph_traverse_specific_person_or_pop_for_supernode_printing(void (*f)(HashTable*, Element *, long* , EdgeArrayType, int, boolean, char**, int*),HashTable * hash_table, 
								     long* supernode_count, 
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
      
      for_test[*index_for_test] = (char*) calloc(length+1,sizeof(char));
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

  printf("Stop using this");
  exit(1);

}



//assume we have checked and no node has >1 chrom intwersection
void db_graph_print_chrom_intersections_for_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, 
						      boolean is_for_testing, char** for_test, int* index_for_test)
{
  //deprecated
  printf("Stop using this");
  exit(1);
}




void db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population(dBGraph* hash_table, EdgeArrayType type, int index)
{
  printf("not implemented yet");
  exit(1);
}

void db_graph_set_status_of_all_nodes_in_supernode(dBNode* node, NodeStatus status, EdgeArrayType type, int index,  dBGraph* dbgraph)
{

  // printf("reimplement db_graph_set_status_of_all_nodes_in_supernode using db_graph_supernode_for_specific_person_or_pop\n");
  
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


  //TODO - fix this to only look at interior nodes in supernode
  // ALSO fix this to call db_graph_supernode

  printf("reimplement to use db_graph_supernode - currently NULL\n");
  

  /*

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
  */

}










dBNode* db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(dBNode* node, EdgeArrayType type, int index, dBGraph* db_graph)
{

  //  printf("Reimplement db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop using db_grapg_supernode");

    
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

void  db_graph_find_population_consensus_supernode_based_on_given_node(Sequence* pop_consensus_supernode, int max_length_of_supernode, dBNode* node, 
								       int min_covg_for_pop_supernode, int min_length_for_pop_supernode, dBGraph* db_graph)
{
  
  //  printf("Reimplement using db_graph_supernode\n");


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


  //TODO - reimplement this whole section using db_graph_supernode. 

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

  //  char tmp_seq[db_graph->kmer_size];
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



// This assumes the given fasta has already been loaded into the graph. All we do here is get the sequence of nodes
// that correspond to the sequence of bases in the file.

// Note that we have to handle the N's we expect to see in a reference fasta. The most reliable thing to do, for these purposes, is to put a NULL node
// in the array for every kmer in the fasta that contains an N

//Note we need to know whether we expect a new entry this time, AND whether there was a new fasta entry last time:
// When about to load a new set of nodes, decide whether to preload seq with the last kmer of the previous set. Do this iff you are not expecting the start of a new entry.
// Only tricky bit is determining the index of the last kmer within seq (from the last time around), when you want to move it to the start of seq just before loading a new set of nodes.
// If last time, we preloaded seq with a kmer at the start (indices 0 to k-1), and we loaded N nodes, then the last kmer in seq is from indices N to N+k-1
// If last time we did not preload seq, and we loaded N nodes, then the last kmer is from N-1 to N+k-2.
// This boils down to knowing whether the last time was a new fasta entry


//Note what happens as you repeatedly call this, with number_of_nodes_to_load = length_of_arrays/2. The first time, you load nodes into the back (right hand) end of the array, and the front end is empty.
//    The next time, these are pushed to the front, and new ones are loaded into the back.
//  the first time you call this, expecting_new_fasta_entry=true,  last_time_was_not_start_of_entry=false
//  the next time,                expecting_new_fasta_entry=false, last_time_was_not_start_of_entry=false,
//  after that                    expecting_new_fasta_entry=false, last_time_was_not_start_of_entry=true

// returns the number of nodes loaded into the array
int db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(FILE* chrom_fptr, 
												     int number_of_nodes_to_load, int number_of_nodes_loaded_last_time,
												     int length_of_arrays, 
												     dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
												     Sequence* seq, KmerSlidingWindow* kmer_window, 
												     boolean expecting_new_fasta_entry, boolean last_time_was_not_start_of_entry, 
												     dBGraph* db_graph ) 
{

 if (length_of_arrays%2 !=0)
    {
      printf("Must only call db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta  with even length_of_arrays\n");
      exit(1);
    }

  if (number_of_nodes_to_load>length_of_arrays)
    {
      printf("Insufficient space in arrays, length %d, to load %d nodes, in db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta\n", 
	     length_of_arrays, number_of_nodes_to_load);
      exit(1);
    }


  //push everything in the arrays left by number_of_nodes_to_load. 
  int i;
  for(i=0; i<length_of_arrays-number_of_nodes_to_load; i++)
    {
      path_nodes[i]        = path_nodes[i+number_of_nodes_to_load];
      path_orientations[i] = path_orientations[i+number_of_nodes_to_load];
      path_labels[i]       = path_labels[i+number_of_nodes_to_load];
      path_string[i]       = path_string[i+number_of_nodes_to_load];
    } 
  for(i=length_of_arrays-number_of_nodes_to_load; i<length_of_arrays; i++)
    {
      path_nodes[i]        = NULL;
      path_orientations[i] = forward;
      path_labels[i]     =  Undefined;
      path_string[i]     = 'N';
    }
  path_string[length_of_arrays-number_of_nodes_to_load]='\0';

  //move a kmer-worth of bases to front of seq, if appropriate
  if(!expecting_new_fasta_entry)
    {

      //sanity
      if (number_of_nodes_loaded_last_time==0)
	{
	  printf("Not expecting a new fasta entry, but did not load any nodes last time - programming error\n");
	  exit(1);
	}
	
      // When about to load a new set of nodes, decide whether to preload seq with the last kmer of the previous set. Do this iff you are not expecting the start of a new entry.
      // Only tricky bit is determining the index of the last kmer within seq (from the last time around), when you want to move it to the start of seq just before loading a new set of nodes. 
      // If last time, we preloaded seq with a kmer at the start (indices 0 to k-1), and we loaded N nodes, then the last kmer in seq is from indices N to N+k-1
      // If last time we did not preload seq, and we loaded N nodes, then the last kmer is from N-1 to N+k-2.
      // This boils down to knowing whether the last time was a new fasta entry
      
      
      if (last_time_was_not_start_of_entry==true) //which will be true most of the time
	{
	  for (i=0; i<db_graph->kmer_size; i++)
	    {
	      seq->seq[i]=seq->seq[number_of_nodes_loaded_last_time+i]; 
	    }
	}
      else
	{
	  for (i=0; i<db_graph->kmer_size; i++)
	    {
	      seq->seq[i]=seq->seq[number_of_nodes_loaded_last_time-1+i]; 
	    }
	}
      seq->seq[db_graph->kmer_size]='\0';
    }

  return load_seq_into_array(chrom_fptr, number_of_nodes_to_load, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  
}



// New improved reference-based SV calling algorithm
// Walks along chromosome path, comparing with supernodes, and for each supernode, sees where it attaches.
// max_anchor_span is the biggest gap we allow/look for between where the 3' and 5' anchors attach. If you want to be able to find, say 10kb deletions, then set this to 10000, etc
// min_fiveprime_flank_anchor is counted in number of nodes
// length_of_arrays MUST BE EVEN
// returns number of variants found
int db_graph_make_reference_path_based_sv_calls(FILE* chrom_fasta_fptr, EdgeArrayType which_array_holds_indiv, int index_for_indiv_in_edge_array,
						int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, int max_anchor_span, int min_covg, int max_covg, 
						int length_of_arrays, dBGraph* db_graph, FILE* output_file)
{
  
  int num_variants_found=0;

  //makes life much simpler to insist the array is even length.
  if (length_of_arrays%2 !=0)
    {
      printf("Must only call db_graph_make_reference_path_based_sv_calls  with even length_of_arrays\n");
      exit(1);
    }



  //prepare arrays. Four will hold 10,000 nodes/orientations/etc, corresponding to the consecutive 10000 bases of the chromosome
  //                Another set of four will be reused to hold supernodes


  int number_of_nodes_to_load=length_of_arrays/2;
  //int number_of_nodes_to_load=max_anchor_span/2;//each time you load, you pull in half an array's worth
  //int length_of_arrays=max_anchor_span;




  //*************************************
  // malloc and initialising
  //*************************************
  dBNode**     chrom_path_array        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); //everything these dBNode*'s point to will be in the hash table - ie owned by the hash table
  Orientation* chrom_orientation_array = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
  Nucleotide*  chrom_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  char*        chrom_string            = (char*) malloc(sizeof(char)*length_of_arrays+1); //+1 for \0

  dBNode**     current_supernode       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  Orientation* curr_sup_orientations   = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  Nucleotide*  curr_sup_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  char*        supernode_string        = (char*) malloc(sizeof(char)*length_of_arrays+1); //+1 for \0



  int n;
  for (n=0; n<length_of_arrays; n++)
    {
      chrom_path_array[n]=NULL;
      chrom_orientation_array[n]=forward;
      chrom_labels[n]=Undefined;
      chrom_string[n]='N';
      current_supernode[n]=NULL;
      curr_sup_orientations[n]=forward;
      curr_sup_labels[n]=Undefined;
      supernode_string[n]='N';
    }
  chrom_string[0]='\0';
  supernode_string[0]='\0';

  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,number_of_nodes_to_load+db_graph->kmer_size+1,LINE_MAX);
  seq->seq[0]='\0';


  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*1000);    //*(number_of_nodes_to_load + db_graph->kmer_size));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  
  kmer_window->nkmers=0;

  // ***********************************************************
  //********** end of malloc and initialise ********************


  //load
  int ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fasta_fptr, number_of_nodes_to_load, 0, 
													     length_of_arrays,
													     chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													     seq, kmer_window, 
													     true, false,
													     db_graph);
  if (ret !=number_of_nodes_to_load)
    {
      printf("db_graph_make_reference_path_based_sv_calls failed on loading first batch. Exit. Sequence loaded  is %s, length %d, and we expect it to be length %d\n", seq->seq, ret, number_of_nodes_to_load);
      exit(1);
    }

  char tmp_seqzam[db_graph->kmer_size+1];
  tmp_seqzam[db_graph->kmer_size]='\0';

  int i;
 
  //one more batch, then array is full, and ready for the main loop:
  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fasta_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, false,
													 db_graph);
  if (ret !=number_of_nodes_to_load)
    {
      printf("db_graph_make_reference_path_based_sv_calls failed on loading second batch. Exit. Sequence loaded  is %s, length %d, and we expect it to be length %d\n", seq->seq, ret, number_of_nodes_to_load);
      exit(1);
    }



  printf("WARNING _ WHAT ABOUT LOOPS?");   
  
  //each call   will push left the nodes/etc in the various arrays by max_anchor_span
  // and then put the ndew nodes etc in on the right of that
  //TODO - change to do--while?
  do
    {
      int start_node_index=0;//this is a position in chrom_path_array 

      while (start_node_index<=max_anchor_span/2)
	{

	  //if this chromosome node does not exist in the person's graph, move on to the next node
	  if (db_node_is_this_node_in_this_person_or_populations_graph(chrom_path_array[start_node_index], which_array_holds_indiv, index_for_indiv_in_edge_array) == false)
	    {
	      start_node_index++;
              break;
	    }

	  int index_of_query_node_in_supernode_array=-1;
	  int length_curr_supernode = db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(chrom_path_array[start_node_index], max_anchor_span, &db_node_condition_always_true, &db_node_action_do_nothing,
										    supernode_string, current_supernode, curr_sup_orientations, curr_sup_labels, &index_of_query_node_in_supernode_array, db_graph,
										    which_array_holds_indiv, index_for_indiv_in_edge_array);


	  if (index_of_query_node_in_supernode_array==-1)
	    {
	      printf("Warning - red alert!!!  - failed to get index of query\n");
	      start_node_index++;
	      break;
	    }

	  //Exclude/ignore hubs and singletons as well as supernodes without enough length to generate our anchors
	  // check if there is room between query node and either end of supernode to fit both anchors
	  if (length_curr_supernode - index_of_query_node_in_supernode_array < min_fiveprime_flank_anchor+min_threeprime_flank_anchor)
	    {
	      start_node_index++;
	      break;
	    }



	  // we now have an array for our supernode, and we know where our query node is within it,  but we do not know 
	  // whether we want to traverse the array backwards or forwards

	  boolean traverse_sup_left_to_right;

	  if (curr_sup_orientations[index_of_query_node_in_supernode_array]==chrom_orientation_array[start_node_index])
	    {
	      //we traverse the supernode array forwards - ie with increasing index
	      traverse_sup_left_to_right=true;
	    }
	  else
	    {
	      traverse_sup_left_to_right=false;
	    }


	  //now see how far along the chromosome you go before the supernode breaks off and differs
	  int first_index_in_chrom_where_supernode_differs_from_chromosome = start_node_index;
	  int index_in_supernode_where_supernode_differs_from_chromosome=index_of_query_node_in_supernode_array;

	  //find how far along can go on reference before supernode branches
	  while ( (chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome]==current_supernode[index_in_supernode_where_supernode_differs_from_chromosome])
		  && (index_in_supernode_where_supernode_differs_from_chromosome < length_curr_supernode ) && (index_in_supernode_where_supernode_differs_from_chromosome>=0)
		  && (first_index_in_chrom_where_supernode_differs_from_chromosome<max_anchor_span) 
		  )
	    {
	      if (traverse_sup_left_to_right==true)
		{
		  index_in_supernode_where_supernode_differs_from_chromosome++;
		}
	      else
		{
		  index_in_supernode_where_supernode_differs_from_chromosome--;
		}
	      first_index_in_chrom_where_supernode_differs_from_chromosome++;
	    }

	  if (  abs(index_in_supernode_where_supernode_differs_from_chromosome-index_of_query_node_in_supernode_array) <= min_fiveprime_flank_anchor-1)
	    {
	      //does not have sufficient 5' anchor
	      start_node_index=first_index_in_chrom_where_supernode_differs_from_chromosome;
	      break;
	    }

	  //note that we do one last increment/decrement at the end of the previous loop, taking us to the point beyond where they last agree
	  if (    (index_in_supernode_where_supernode_differs_from_chromosome==length_curr_supernode)
		  || (index_in_supernode_where_supernode_differs_from_chromosome==-1) 
	     )
	    {
	      //then actually the whole supernode matches the reference exactly.
	      start_node_index = index_in_supernode_where_supernode_differs_from_chromosome+1;
	      break;
	    }

	  //We now have  decent 5-prime anchor, and we know we don't match the ref exactly.
	  // Now see if the END of the supernode matches anywhere in our current chuk of reference chromosome

	  int start_of_3prime_anchor_in_sup;

	  if (traverse_sup_left_to_right==true)
	    {
	      start_of_3prime_anchor_in_sup=length_curr_supernode-min_threeprime_flank_anchor ;
	    }
	  else
	    {
	      start_of_3prime_anchor_in_sup = 0 + min_threeprime_flank_anchor;
	    }

	  int start_of_3prime_anchor_in_chrom = first_index_in_chrom_where_supernode_differs_from_chromosome+1;
	  boolean found_other_anchor=false;

	  //walk along entire chromosome array, and see if can attach the 3prime anchor
	  // it may attach multiple times, but we will only find the closest
	  while ( (found_other_anchor==false) && ( start_of_3prime_anchor_in_chrom < max_anchor_span-min_threeprime_flank_anchor))
	    {
	      boolean potential_anchor=true;
	      int j;
	      for (j=0 ; (j< min_threeprime_flank_anchor) && (potential_anchor==true); j++)
		{
		  int k=j;
		  if (traverse_sup_left_to_right==false)
		    {
		      k=-j;
		    }

		  if (chrom_path_array[start_of_3prime_anchor_in_chrom + j] != current_supernode[start_of_3prime_anchor_in_sup + k])
		    {
		      potential_anchor=false;
		    }
		  if (potential_anchor==true)
		    {
		      found_other_anchor=true;
		    }

		}
	      start_of_3prime_anchor_in_chrom++;
	    }
	  
	  
	  if (found_other_anchor==false)
	    {
	      start_node_index++;
	      break;
	    }
	  else
	    {
	      //we have found a potential SV locus.
	      num_variants_found++;

	      //print reference section
	      printf("VARIANT, 5' anchor %d, 3' anchor %d, ref branch length %d, individual branch length %d \n", first_index_in_chrom_where_supernode_differs_from_chromosome-start_node_index, min_threeprime_flank_anchor, 
		     start_of_3prime_anchor_in_chrom-first_index_in_chrom_where_supernode_differs_from_chromosome, abs(start_of_3prime_anchor_in_sup-index_in_supernode_where_supernode_differs_from_chromosome));
	      printf("Ref\n");
	      int k;
	      for (k=start_node_index; k< start_of_3prime_anchor_in_chrom+min_threeprime_flank_anchor; k++)
		{
		  printf("%c", chrom_string[k]);
		}
	      printf("Supernode\n");

	      if (traverse_sup_left_to_right==true)
		{
		  
		  for (k=index_of_query_node_in_supernode_array; k< start_of_3prime_anchor_in_sup+min_threeprime_flank_anchor; k++)
		    {
		      printf("%c", supernode_string[k]);
		    }
		}
	      else
		{
		  for (k=index_of_query_node_in_supernode_array; k>start_of_3prime_anchor_in_sup-min_threeprime_flank_anchor; k--)
		    {
		      printf("%c", supernode_string[k]);
		    }
		}
	      
	    }

	}

      

    }  while (db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fasta_fptr, number_of_nodes_to_load, number_of_nodes_to_load,
													       length_of_arrays, 
													       chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													       seq, kmer_window, 
													       false, true, db_graph)>0);





  return num_variants_found;


}
