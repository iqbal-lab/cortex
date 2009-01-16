/*
  dB_graph_population.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>

#include <element.h>
#include <hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <pqueue_pop.h>

#include <string.h>



boolean db_graph_is_this_node_in_this_person_or_populations_graph(dBNode* node, EdgeArrayType type, int index)
{
 
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





// wrapper for hash_table_find, which allows you to look in the hash table
// specifically for nodes related to a specific person or population
// person or population  specified by which edge array type
// which person or pop specified by index

dBNode *  db_graph_find_node_restricted_to_specific_person_or_population(Key key, dBGraph * hash_table, EdgeArrayType type, int index)

{

  dBNode *  e = hash_table_find(key, hash_table);

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


char * get_seq_from_elem_to_end_of_supernode_for_specific_person_or_pop(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle, EdgeArrayType type, int index){
  Nucleotide nucleotide1, nucleotide2, rev_nucleotide;
  char * seq = NULL;
  int max = 1500;
  Orientation original_orientation, next_orientation;
  dBNode * original_node;
  dBNode * next_node;
  int seq_length = 0;
  
  

  seq = malloc(1500*sizeof(char)); 

  if (seq == NULL){
    puts("dB_graph: cannot assign memory for sequence\n");
    exit(1);
  }

  original_node = node;
  original_orientation = orientation; 
  

  //Mark the element we're starting at as visited
  db_node_set_status(node,visited);

  *is_cycle = false;

  
  while(db_node_has_precisely_one_edge(node,orientation,&nucleotide1, type, index)) {
 
    if ((seq_length+1)>max){
      max+=1000;
      seq = (char*)realloc(seq, sizeof(char) * max);
    }

    if (seq == NULL){
      printf("dB_graph: cannot assign memory\n");
      exit(1);
    }

   
    next_node =  db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,nucleotide1,&rev_nucleotide,db_graph, type, index);
    
    if(next_node == NULL){
      printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
      exit(1);
    }	         
    
    if (DEBUG){
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide1),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size) :  binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size),db_graph->kmer_size));
	     

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

  //clever thing about this next line is that if it turns out that the very first node that you handed in as an argument
  //had an edge back to itself - so it was a cycle, then the first character in seq is \0 and when the person who calls this prints what you return, they print nothing
  seq[seq_length] = '\0';
  return seq;
}


void print_supernode_for_specific_person_or_pop(HashTable* db_graph, dBNode * node,EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test)
{

  FILE * fout;
  long count=0;
  
  char filename [200];
  if (count % 100000000 == 0)
    {
      int num = count / 100000000;
      
      if (count !=0)
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
  
  db_graph_print_supernode_for_specific_person_or_pop(fout,node,db_graph, type,index, is_for_testing,  for_test,index_for_test);
 


}


// *********************************************
// functions applied to a person/pop's graph
// *********************************************


void db_graph_traverse_specific_person_or_pop(void (*f)(HashTable*, Element *, EdgeArrayType, int, boolean, char**, int*),HashTable * hash_table, EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test){

  int i;
  for(i=0;i<hash_table->number_buckets;i++){
    pqueue_traverse_specific_person_or_pop(f,hash_table, &(hash_table->table[i]), type,index, is_for_testing, for_test, index_for_test);
  }

}

//will print only nodes for specific person/population

// ******NOTE****** if you want to go on to print for other people, you need to set status to none for all the nodes again after going thrrough the whole
// graph printing supernodes - otherwise some will be left marked visited.

//if you pass in true for the boolean is_for_testing,  then instead of printing to the file, it adds the supernode to the array of char*'s 
//this enables the caller to sort the supernodes lexicographically, and test that the right supernodes are being found
//this is not scalable to a ful genome, but sufficient for detailed unit tests

void db_graph_print_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test ){


  //don't do anything if this node does not occur in the graph for this person/population
  if (!db_graph_is_this_node_in_this_person_or_populations_graph(node, type, index))
    {
      return;
    }

  char * seq = NULL;
  char * seq_forward = NULL;
  char * seq_reverse = NULL; 
  char * seq_reverse_reversed = NULL;
  int i;
  int length_reverse = 0;
  boolean is_cycle_forward, is_cycle_reverse;

  //printf("PRINT SUPERNODE for person %d\n", index);
  if (! db_node_check_status(node,visited) && ! db_node_check_status(node,pruned))
    {
      
      //this function binary_kmer_to_seq does the malloc for what it returns.
      seq = binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size);
      
      if (DEBUG){
      printf("\nSTART Supernode %s\n",seq);    
      printf("go forward\n");
      }
      
      //compute the forward path until the end of the supernode
      //mark the nodes in the path as visited.
 
     //return is_cycle_forward == true if the path closes a loop
 
      //this function mallocs the sequence that it returns, and will always initialise the is_cycle_forward variable, that is currently uninitialised
      //note that if going forward there is a cycle, then seq_forward will contain garbage, so don't use it. 
     seq_forward = get_seq_from_elem_to_end_of_supernode_for_specific_person_or_pop(node,forward,db_graph,&is_cycle_forward, type, index);
      
     if (DEBUG){
       printf("Z NODE c %s\n",seq); 
       if (is_cycle_forward)
	 {
	   printf("NODE f - none since is cycle forward\n");
	 }
       else
	 {
	   printf("NODE f %s\n",seq_forward);
	 }
     }
      
      if (! is_cycle_forward)
	{
    	  if (DEBUG){
	  printf("go reverse\n");
	  }
	  
	  //compute the reverse path...
	  //this function mallocs the sequence that it returns, andwill always initialise the is_cycle_reverse variable, that is currently uninitialised
	  //printf("getting rev seq for person %d", index);
	  seq_reverse = get_seq_from_elem_to_end_of_supernode_for_specific_person_or_pop(node,reverse,db_graph,&is_cycle_reverse, type, index);
	  
	  if (is_cycle_reverse){
	    puts("cycle in the reverse orientation without cycle in the forward orientation - impossible!\n");
	    exit(1);
	  }
	
	  if (DEBUG)
	    {
	      printf("Z NODE r %s\n",seq_reverse);
	    }
	  
	  length_reverse = strlen(seq_reverse);
	}

      seq_reverse_reversed = malloc((length_reverse+1)*sizeof(char));


      if (seq_reverse_reversed==NULL)
	{
	  printf("Cannot mallo cseq rev rev");
	  exit(1);
	}
      
      
      //reverse the reverse sequence
      for(i=0;i<length_reverse;i++){
	seq_reverse_reversed[i] = reverse_char_nucleotide(seq_reverse[length_reverse-i-1]);
	//printf("ZAMZAM seq rev rev next char is %c\n", seq_reverse_reversed[i]);
      }
      seq_reverse_reversed[length_reverse]='\0';
      //so if we had a cycle in the forward/reverse directions, this would only contain the single character \0

      if (DEBUG){
	if (length_reverse>0)
	  {
	    printf("Z NODE rr %s\n", seq_reverse_reversed);
	  }
      }
      
      //printf(">NODE\n%s%s%s\n",seq_reverse_reversed,seq,seq_forward); 
      
      if (!is_for_testing) 
	{
	  fprintf(file,">NODE\n%s%s%s\n",seq_reverse_reversed,seq,seq_forward); 
	}
      else
	{
	  int length_of_supernode = strlen(seq_reverse_reversed)+strlen(seq) +strlen(seq_forward) ;
	  if (length_of_supernode==0)
	    {
	      printf("Null supernode");
	      exit(1);
	    }


	  //assume the called has malloced space for output

	  for_test[*index_for_test] = (char*) calloc(length_of_supernode+1,sizeof(char));
	  if (for_test[*index_for_test]==NULL)
	    {
	      printf("Unable to calloc for supernode");
	      exit(1);
	    }
	  for_test[*index_for_test][0]='\0';
	  strcat(for_test[*index_for_test],seq_reverse_reversed);
	  strcat(for_test[*index_for_test],seq);
	  strcat(for_test[*index_for_test],seq_forward);
	  for_test[*index_for_test][length_of_supernode]='\0';
	  
	  (*index_for_test)++;

	}


      free(seq);
      free(seq_forward);
      free(seq_reverse);
      free(seq_reverse_reversed);
	  
    }
  else{
    if (DEBUG){
      if ( db_node_check_status(node,visited)){
	printf("\n%s: visited\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
      }
      if ( db_node_check_status(node,pruned)){
	printf("\n%s: pruned\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
      }
    }
  }


}




void db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population(dBGraph* hash_table, EdgeArrayType type, int index)
{
  printf("not implemented yet");
  exit(1);
}

