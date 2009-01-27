
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
#include <seq.h>
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


//WARNING: this marks as visited any node that it walks over (and is in the supernode)
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
	  
	  //Now make sure you are using the smaller of the sequence and its rev complement
	  BinaryKmer snode_kmer = seq_to_binary_kmer(for_test[*index_for_test],length_of_supernode);
	  BinaryKmer rev_snode_kmer =  binary_kmer_reverse_complement(snode_kmer, length_of_supernode);
	  
	  if (rev_snode_kmer<snode_kmer)
	    {
	      free(for_test[*index_for_test]);
	      for_test[*index_for_test] = binary_kmer_to_seq(rev_snode_kmer,length_of_supernode);
	    }
	  
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


void return_node_to_status_prior_to_tmp_touched_X(dBNode* node)
{
  if (node->status == tmp_touched_previously_none)
    {
      node->status = none;
    }
  else if (node->status == tmp_touched_previously_visited)
    {
      node->status = visited;
    }
  else if (node->status == tmp_touched_previously_pruned)
    {
      node->status = pruned;
    }
}


//scratch_node_array is prealloced. caller should ignore its contents - is used purely internally.
dBNode* db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(dBNode* node, EdgeArrayType type, int index, dBGraph* db_graph, dBNodeArray* scratch_node_array)
{
  printf("start of get first node function\n");
  if (! (db_graph_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      return NULL;
    }

  if (db_node_check_status(node,pruned))
    {
      //don't waste time with pruned nodes.
      return NULL;
    }
  
  //initialise scratch node array length to 0;
  scratch_node_array->length=0;

  boolean is_cycle;
  Nucleotide nucleotide1, nucleotide2, rev_nucleotide;
  Orientation original_orientation, next_orientation, orientation;
  dBNode * original_node=node;
  dBNode * next_node;
  int j;

  //First node in supernode is, almost by definition, what you get if you go in the Reverse direction (with respect to the Node)
  // as far as you can go.
  original_orientation = reverse; 
  orientation = reverse;
  

  //Mark the element we're starting at as visited

  if (db_node_check_status(original_node,none))
    {
      db_node_set_status(original_node,tmp_touched_previously_none);
    }
  else if (db_node_check_status(original_node,visited))
    {
      db_node_set_status(original_node,tmp_touched_previously_visited);
    }
  else 
    {
      printf("This node %s has status that is not visited or none or pruned",binary_kmer_to_seq(original_node->kmer, db_graph->kmer_size)); //creates a mem leak but about to die anyway
      exit(1);
    }
  

  //add this node to scratch array, so that we can unset its status afterwards
  if ((scratch_node_array->length)+1 > scratch_node_array->max)
    {
      printf("Out of mem in scratch node array");
      exit(1);
    }
  (scratch_node_array->nodes)[scratch_node_array->length]=original_node;
  scratch_node_array->length= (scratch_node_array->length) +1;



  is_cycle = false;

  
  while(db_node_has_precisely_one_edge(node,orientation,&nucleotide1, type, index)) {
 
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
      
      for (j=0; j< scratch_node_array->length; j++)
	{
	  return_node_to_status_prior_to_tmp_touched_X((scratch_node_array->nodes)[j]);
	}
      printf("returning this first nodem, with kmer %s\n", binary_kmer_to_seq(node->kmer, db_graph->kmer_size));
      return node; //we have gone as far as we can go - the next node has multiple entries. So we are now at the first node of the supernode
      }
    
    
    //loop
    if ((next_node == original_node) && (next_orientation == original_orientation))
      {      
	is_cycle = true;
	//break;
	
	for (j=0; j< scratch_node_array->length; j++)
	  {
	    return_node_to_status_prior_to_tmp_touched_X((scratch_node_array->nodes)[j]);
	  }
      printf("returning this first nodem, with kmer %s\n", binary_kmer_to_seq(original_node->kmer, db_graph->kmer_size));
	return original_node; //we have a loop that returns to where we start. Might as well consider ourselves as at the fiurst node of the supernode right at the beginning
      }
    
    
    //   db_node_set_status(next_node,visited);
    if (db_node_check_status(next_node,none))
      {
	db_node_set_status(next_node,tmp_touched_previously_none);
      }
    else if (db_node_check_status(next_node,visited))
      {
	db_node_set_status(next_node,tmp_touched_previously_visited);
      }
    else
      {
	printf("This node %s has status that is not visited or none or pruned",binary_kmer_to_seq(next_node->kmer, db_graph->kmer_size)); //creates a mem leak but about to die anyway
	exit(1);
      }
    //add this node to scratch array, so that we can unset its status afterwards
    if (scratch_node_array->length+1 > scratch_node_array->max)
      {
	printf("Out of mem in scratch node array");
	exit(1);
      }
    (scratch_node_array->nodes)[scratch_node_array->length]=next_node;
    scratch_node_array->length= (scratch_node_array->length) +1;

    node = next_node;
    orientation = next_orientation;      
  }


  //if we get here then something has gone wrong
  printf("Reached end of function without finding first node in supernode");
  exit(1);
  return NULL; //to keep compiler happy

}



// Given a node, will look at each person's supernode containing that node. 
// For each of these supernodes, independently,  find the longest contiguous subsection that has min people-coverage > min_covg_for_pop_supernode.
// (Note this may not contain the original node.)
// Then compare across people, and define your consensus as the sub_supernode which has the most people backing it.
// If two people have equally good sub_supernodes, choose that of the person with highest index.
// Pass in pre-malloced Sequence* for answer
// pass in pre-malloced dBNodeArray for internal workings - do not look at contents afterwards.

void  db_graph_find_population_consensus_supernode_based_on_given_node(Sequence* pop_consensus_supernode, dBNode* node, int min_covg_for_pop_supernode, int min_length_for_pop_supernode, 
								       dBGraph* db_graph, dBNodeArray* scratch_node_array)
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
      //      Edges person_edges = get_edge_copy(*node, individual_edge_array, i);

      //get the first node in the supernode for this edge, within this person't graph
      dBNode* first_node = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, individual_edge_array, i, db_graph, scratch_node_array);
      //clean up those nodes that were set with status tmp_touched_previously_X to X
      
//db_graph_return_supernode_statuses_to_previous_state(node,individual_edge_array, i);

      // walk through the supernode, from one end to the other. As soon as we hit a node that has people-covg >= min_covg_for_pop_supernode
      // then note down the coordinate (how far through the supernode), and keep going until we meet a node with people-covg < min_covg_for_pop_supernode
      // Work out the length of the chunk you found, which had sufficient covg. If > min_length_for_pop_supernode, then you have found your first
      // candidate for this person's best sub-supernode. Log the start-point and length. Then continue walkign through the supernode until you
      // find another node that has people_covg sufficiently high. repeat the process for this point, and if its candidate is longer than our
      // current candidate, then replace it with this new one. etc.
      
 //  dB_graph_get_best_sub_supernode_given_min_covg_and_length_for_specific_person_or_pop(first_node,  &index_of_start_of_best_sub_supernode_in_each_person[i], &length_of_best_sub_supernode_in_each_person[i], 
      //											   min_covg_for_pop_supernode,  min_length_for_pop_supernode, individual_edge_array, i); 

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

  if (max==0)
    {
      printf("This should be impossible. Max size of sub supernode over all people is 0. Since we start wth a kmer that is in the graph, at least one person ought to have it\n");
    }
  
  else
    {
      if (max > pop_consensus_supernode->max)
	{
	  pop_consensus_supernode->seq = (char*) realloc(pop_consensus_supernode->seq, sizeof(char*) * max);
	  if (pop_consensus_supernode->seq == NULL)
	    {
	      printf("OOM getting pop cons supernode\n");
	      exit(1);
	    }

	  pop_consensus_supernode->max=max;
	}
      int j;
      int start = index_of_start_of_best_sub_supernode_in_each_person[person_with_best_sub_supernode];
      int end =    start+ length_of_best_sub_supernode_in_each_person[person_with_best_sub_supernode];
      db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(pop_consensus_supernode->seq, node, start, end, individual_edge_array, person_with_best_sub_supernode);
    }

}


void db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(char* subsection, dBNode* node, int start, int end, EdgeArrayType type, int index)
{
}
