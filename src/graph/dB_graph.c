/*
  dB_graph.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <string.h>

dBNode * db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, 
			       Orientation * next_orientation,
			       Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph){
  
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


void db_graph_remove_path(dBNode * node, Orientation orientation, int length, dBGraph * db_graph){

  Nucleotide nucleotide, rev_nucleotide;
  dBNode * next_node;
  Orientation next_orientation;
  int i;

  for(i = 0; i < length; i++){

    if (db_node_has_precisely_one_edge(node,orientation,&nucleotide)){

      if (DEBUG){
	printf("CLIPPING %c\n",binary_nucleotide_to_char(nucleotide));
      }
      
      next_node = db_graph_get_next_node(node,orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph);

      if (DEBUG){
	printf("RESET %c BACK\n",binary_nucleotide_to_char(rev_nucleotide));
      }
 
      

      if(next_node == NULL){
	printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
	exit(1);
      }	           
      
      
      db_node_set_status(node,pruned);
      db_node_reset_edges(node);

      db_node_reset_edge(next_node,opposite_orientation(next_orientation),rev_nucleotide);

    }
    else{
      printf("dB_graph: assume unique path kmer:%s\n",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
      exit(1);
    }
    
    node = next_node;
    orientation = next_orientation;
  }
}

int db_graph_detect_tip(dBNode * node, Orientation orientation, int limit,dBGraph * db_graph){ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end(node, opposite_orientation(orientation))){
    
    boolean join_main_trunk = false;

     while(db_node_has_precisely_one_edge(node,orientation,&nucleotide)) {
       Orientation next_orientation;
       dBNode * next_node;
       
       next_node = db_graph_get_next_node(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
       
       if(next_node == NULL){
	 printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
	 exit(1);
       }	           
       
       length ++;

       
       if (length>limit){
	 break;
       }

       //want to stop when we join another trunk
       if (!db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide)){
	 join_main_trunk = true;
	 break;
       }

       
       node = next_node;
       orientation = next_orientation;
      
    
     }
  
  
     if (! join_main_trunk){
       length = 0;
     }

  }

  return length;
}



char * get_seq_from_elem_to_end_of_supernode(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle){
  Nucleotide nucleotide1, nucleotide2, rev_nucleotide;
  char * seq = NULL;
  int max = 1000;
  Orientation original_orientation, next_orientation;
  dBNode * original_node;
  dBNode * next_node;
  int seq_length = 0;
  
  

  seq = malloc(1000*sizeof(char)); 

  if (seq == NULL){
    puts("dB_graph: cannot assign memory for sequence\n");
    exit(1);
  }

  original_node = node;
  original_orientation = orientation; 
  

  //Mark the element we're starting at as visited
  db_node_set_status(node,visited);

  *is_cycle = false;

  
  while(db_node_has_precisely_one_edge(node,orientation,&nucleotide1)) {
 
    if ((seq_length+1)>max){
      max+=1000;
      seq = (char*)realloc(seq, sizeof(char) * max);
    }

    if (seq == NULL){
      printf("dB_graph: cannot assign memory\n");
      exit(1);
    }

   
    next_node =  db_graph_get_next_node(node,orientation,&next_orientation,nucleotide1,&rev_nucleotide,db_graph);
    
    if(next_node == NULL){
      printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
      exit(1);
    }	         
    
    if (DEBUG){
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide1),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size) :  binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size),db_graph->kmer_size));
	     

    }
    
    //check for multiple entry edges 
    if (db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2)){
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

void db_graph_print_supernode(FILE * file, dBNode * node, dBGraph * db_graph){

  char * seq = NULL;
  char * seq_forward = NULL;
  char * seq_reverse = NULL; 
  char * seq_reverse_reversed = NULL;
  int i;
  int length_reverse = 0;
  boolean is_cycle_forward, is_cycle_reverse;

  if (! db_node_check_status(node,visited) && ! db_node_check_status(node,pruned)){
  
    seq = binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size);

    if (DEBUG){
      printf("\nSTART Supernode %s\n",seq);    
      printf("go forward\n");
    }
    
    //compute the forward path until the end of the supernode
    //mark the nodes in the path as visited.
    //return is_cycle_forward == true if the path closes a loop
    seq_forward = get_seq_from_elem_to_end_of_supernode(node,forward,db_graph,&is_cycle_forward);
    
    if (DEBUG){
      printf("NODE c %s\n",seq); 
      printf("NODE f %s\n",seq_forward);
    }
    
    if (! is_cycle_forward){
      
      if (DEBUG){
	printf("go reverse\n");
      }
      
      //compute the reverse path...
      seq_reverse = get_seq_from_elem_to_end_of_supernode(node,reverse,db_graph,&is_cycle_reverse);
      
      if (is_cycle_reverse){
	puts("cycle reverse orientation without cycle in the forward orientation\n");
	exit(1);
      }
      
      if (DEBUG){
	printf("NODE r %s\n",seq_reverse);
      }
      
      length_reverse = strlen(seq_reverse);
    }
    
    seq_reverse_reversed = malloc((length_reverse+1)*sizeof(char));
    
    //reverse the reverse sequence
    for(i=0;i<length_reverse;i++){
      seq_reverse_reversed[i] = reverse_char_nucleotide(seq_reverse[length_reverse-i-1]);
      
    }
    seq_reverse_reversed[length_reverse]='\0';
    
    if (DEBUG){
      printf("NODE rr %s\n",seq_reverse_reversed);
    }
    
    fprintf(file,">NODE\n%s%s%s\n",seq_reverse_reversed,seq,seq_forward); 
    
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

void db_graph_clip_tip(dBNode * node, int limit,dBGraph * db_graph){

  int length_tip;

  

  if (DEBUG){
    printf("Node %s\n",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
  }

  if (!db_node_check_status(node,pruned)){

    length_tip = db_graph_detect_tip(node,forward,limit,db_graph);
    if (DEBUG){
      printf("\tforward tip: %i\n",length_tip);
    }
    
    if (length_tip!=0){
      
      db_graph_remove_path(node,forward,length_tip,db_graph);
      
    }
    else{
      length_tip = db_graph_detect_tip(node,reverse,limit,db_graph);
      if (DEBUG){
	printf("\treverse tip: %i\n",length_tip);
      }
      
      if (length_tip!=0){
	db_graph_remove_path(node,reverse,length_tip,db_graph);
      
      }
    }     
  }
  else{
    if (DEBUG){
      if ( db_node_check_status(node,pruned)){
	printf("\n%s: pruned\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size));
      }
    }
  }
}

