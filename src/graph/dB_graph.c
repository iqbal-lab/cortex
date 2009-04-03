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

  char seq[db_graph->kmer_size];
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
	printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
	exit(1);
      }	           
      
      
      db_node_set_status(node,pruned);
      db_node_reset_edges(node);

      db_node_reset_edge(next_node,opposite_orientation(next_orientation),rev_nucleotide);

    }
    else{
      printf("dB_graph: assume unique path kmer:%s\n",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
      exit(1);
    }
    
    node = next_node;
    orientation = next_orientation;
  }
}

int db_graph_detect_tip(dBNode * node, Orientation orientation, int limit,dBGraph * db_graph){ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  char seq[db_graph->kmer_size];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end(node, opposite_orientation(orientation))){
    
    boolean join_main_trunk = false;

     while(db_node_has_precisely_one_edge(node,orientation,&nucleotide)) {
       Orientation next_orientation;
       dBNode * next_node;
       
       next_node = db_graph_get_next_node(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
       
       if(next_node == NULL){
	 printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
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


int db_graph_clip_tip_with_orientation(dBNode * node, Orientation orientation, int limit,dBGraph * db_graph){ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  int i;
  dBNode * nodes[limit];
  Orientation next_orientation;
  dBNode * next_node;
  char seq[db_graph->kmer_size];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end(node, opposite_orientation(orientation))){
    
    boolean join_main_trunk = false;

     while(db_node_has_precisely_one_edge(node,orientation,&nucleotide)) {
  
       nodes[length] = node;

       next_node = db_graph_get_next_node(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
       
       if(next_node == NULL){
	 printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
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

       //keep track of node
       

       node = next_node;
       orientation = next_orientation;      
    
     }
  
  
     if (! join_main_trunk){
       length = 0;
     }
     else{//clear edges and mark nodes as pruned
       for(i=0;i<length;i++){

	 if (DEBUG){
	   printf("CLIPPING node: %s\n",binary_kmer_to_seq(element_get_kmer(nodes[i]),db_graph->kmer_size,seq));
	 }

	 db_node_set_status(nodes[i],pruned);
	 db_node_reset_edges(nodes[i]);
       }

       if (DEBUG){
	printf("RESET %c BACK\n",binary_nucleotide_to_char(reverse_nucleotide));
      }
       db_node_reset_edge(next_node,opposite_orientation(next_orientation),reverse_nucleotide);
     }

  }

  return length;
}

char * get_seq_from_elem_to_end_of_supernode(dBNode * node, Orientation orientation, dBGraph * db_graph, boolean * is_cycle,char * seq, int max_length){
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

  
  while(db_node_has_precisely_one_edge(node,orientation,&nucleotide1)) {
 
    if ((seq_length+1)> max_length){
      fprintf (stderr,"cannot allocate a sequence longer than max length [%i]\n",max_length);
      exit(1);
    }

    next_node =  db_graph_get_next_node(node,orientation,&next_orientation,nucleotide1,&rev_nucleotide,db_graph);
    
    if(next_node == NULL){
      fprintf(stderr,"dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq));
      exit(1);
    }	         
    
    if (DEBUG){
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide1),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq) :  binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size),db_graph->kmer_size,tmp_seq));
	     

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

  char seq[db_graph->kmer_size];
  char seq_forward[5000];
  char seq_reverse[5000]; 
  char seq_reverse_reversed[5000];
  int i;
  int length_reverse = 0;
  boolean is_cycle_forward, is_cycle_reverse;
 

  if (! db_node_check_status(node,visited) && ! db_node_check_status(node,pruned)){
  
    binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq);

    if (DEBUG){
      printf("\nSTART Supernode %s\n",seq);    
      printf("go forward\n");
    }
    
    //compute the forward path until the end of the supernode
    //mark the nodes in the path as visited.
    //return is_cycle_forward == true if the path closes a loop
    get_seq_from_elem_to_end_of_supernode(node,forward,db_graph,&is_cycle_forward,seq_forward,5000);
    
    if (DEBUG){
      printf("NODE c %s\n",seq); 
      printf("NODE f %s\n",seq_forward);
    }
    
    if (! is_cycle_forward){
      
      if (DEBUG){
	printf("go reverse\n");
      }
      
      //compute the reverse path...
      get_seq_from_elem_to_end_of_supernode(node,reverse,db_graph,&is_cycle_reverse,seq_reverse,5000);
      
      if (is_cycle_reverse){
	puts("cycle reverse orientation without cycle in the forward orientation\n");
	exit(1);
      }
      
      if (DEBUG){
	printf("NODE r %s\n",seq_reverse);
      }
      
      length_reverse = strlen(seq_reverse);
    }
    
   
    
    //reverse the reverse sequence
    for(i=0;i<length_reverse;i++){
      seq_reverse_reversed[i] = reverse_char_nucleotide(seq_reverse[length_reverse-i-1]);
      
    }
    seq_reverse_reversed[length_reverse]='\0';
    
    if (DEBUG){
      printf("NODE rr %s\n",seq_reverse_reversed);
    }
    
    fprintf(file,">NODE\n%s%s%s\n",seq_reverse_reversed,seq,seq_forward); 
    
  }
  else{
    if (DEBUG){
      if ( db_node_check_status(node,visited)){
	printf("\n%s: visited\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
      }
      if ( db_node_check_status(node,pruned)){
	printf("\n%s: pruned\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
      }
    }
  }
}

int db_graph_clip_tip(dBNode * node, int limit,dBGraph * db_graph){

  int length_tip = 0;
  char seq[db_graph->kmer_size];

  if (DEBUG){
    printf("START Node %s\n",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
  }

  if (!db_node_check_status(node,pruned)){

    length_tip = db_graph_clip_tip_with_orientation(node,forward,limit,db_graph);
   
    
    if (length_tip!=0){
      if (DEBUG){
	 printf("\tforward tip: %i\n",length_tip);
       }
    }
    else{
      length_tip = db_graph_clip_tip_with_orientation(node,reverse,limit,db_graph);
      
      if (length_tip!=0){
	if (DEBUG){
	 printf("\treverse tip: %i\n",length_tip);
	}
      }     
    }
  }
  else{
    if (DEBUG){
      if ( db_node_check_status(node,pruned)){
	printf("\n%s: pruned\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
      }
    }
  }
  
  return length_tip;
}


// perfect path -- no conflict no cycle -- returns length
// path node_1 edge_1 node_2 edge_2 ... node_n edge_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n+1 array from 0..n with all the labels for the edges (see above)
int db_graph_get_perfect_path(dBNode * node, Orientation orientation, dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, boolean * is_cycle, int limit, NodeStatus status, dBGraph * db_graph){
  Orientation  current_orientation,next_orientation;
  dBNode * current_node;
  dBNode * next_node;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size];

  //sanity check
  if (node == NULL){
    printf("db_graph_get_perfect_path: can't pass a null node\n");
    exit(1);
  }

  current_node        = node;
  current_orientation = orientation;

  *is_cycle = false;

  path_nodes[0]         = node;
  path_orientations[0]  = orientation;
  db_node_set_status(node,status);

  if (DEBUG){
     printf("\nNode in path: %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
   }

  do{ 

    if (!db_node_has_precisely_one_edge(current_node, current_orientation,&nucleotide)){
      break;
    }
    else{
      next_node =  db_graph_get_next_node(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph);
      
      //sanity check
      if(next_node == NULL){
	fprintf(stderr,"dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
	exit(1);
      }

      db_node_set_status(next_node,status);

      if (DEBUG){
	printf("\nNode in path: %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
      }
      
      path_labels[length]        = nucleotide;
      length++;
      path_nodes[length]         = next_node;
      path_orientations[length]  = next_orientation;

      current_node        = next_node;
      current_orientation = next_orientation;
    }
  } while (length<limit && 
	   !((next_node == node) && (next_orientation == orientation)) && //loop
	   db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2)); //multiple entries
	 
  if (DEBUG){
    printf("\nLast node in path: %s %i\n", binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq),next_node->edges);
  }
  
  return length;

}

//a perfect bubble starts in a node with only two outgoing edges in the same orientation
// every branch of the bubble is free of conflicts until it joins with the main branch
// the length of every branch is the same -- smaller than limit!

boolean db_graph_detect_perfect_bubble(dBNode * node,
				       Orientation * orientation,
				       Nucleotide * base1, Nucleotide * base2, 
				       Nucleotide * labels,dBNode * * end_node,Orientation * end_orientation,
				       dBGraph * db_graph){

  

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

  if (! db_node_check_status(node,visited)){

    if (db_node_has_precisely_two_edges(node,forward,base1,base2)){
    
      *orientation = forward;
      next_node1 = db_graph_get_next_node(node,forward,&next_orientation1,*base1,&rev_base1,db_graph);
      next_node2 = db_graph_get_next_node(node,forward,&next_orientation2,*base2,&rev_base2,db_graph);
      
      if (next_node1 == NULL || next_node2 == NULL){
	puts("error!");
	exit(1);
      }

      length1  = db_graph_get_perfect_path(next_node1,next_orientation1,nodes1,orientations1,labels1,&is_cycle1,db_graph->kmer_size,visited,db_graph);
      length2  = db_graph_get_perfect_path(next_node2,next_orientation2,nodes2,orientations2,labels2,&is_cycle2,db_graph->kmer_size,visited,db_graph);

    }
    else{
      if (db_node_has_precisely_two_edges(node,reverse,base1,base2)){

	*orientation = reverse;
	next_node1 = db_graph_get_next_node(node,reverse,&next_orientation1,*base1,&rev_base1,db_graph);
	next_node2 = db_graph_get_next_node(node,reverse,&next_orientation2,*base2,&rev_base2,db_graph);
	
	if (next_node1 == NULL || next_node2 == NULL){
	  puts("error!");
	  exit(1);
	}
	
	length1 = db_graph_get_perfect_path(next_node1,next_orientation1,nodes1,orientations1,labels1,&is_cycle1,db_graph->kmer_size,visited,db_graph);
	length2 = db_graph_get_perfect_path(next_node2,next_orientation2,nodes2,orientations2,labels2,&is_cycle2,db_graph->kmer_size,visited,db_graph);

      }
    }
  

    db_node_set_status(node,visited);

   
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
    printf("db_graph_detect_perfect_bubble: error inconsistent state\n");
    exit(1);
  }
  return ret;
}
