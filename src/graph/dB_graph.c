

/*
  dB_graph.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <string.h>
#include <limits.h>

void print_fasta_from_path(FILE *fout,
			   char * name,
			   int length,
			   double avg_coverage,
			   int min_coverage,
			   int max_coverage,
			   dBNode * fst_node,
			   Orientation fst_orientation,
			   dBNode * lst_node,
			   Orientation lst_orientation,
			   char * string, //labels of paths
			   int kmer_size,
			   boolean include_first_kmer);

//it doesn't check that it is a valid arrow -- it just assumes the arrows is fine
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


boolean db_graph_db_node_has_precisely_n_edges_with_status(dBNode * node,Orientation orientation,NodeStatus status,int n,
						    dBNode * * next_node, Orientation * next_orientation, Nucleotide * next_base,
						    dBGraph * db_graph){

  int count = 0;
  boolean ret = false;

  void check_edge(Nucleotide base){

    dBNode * tmp_next_node;
    Orientation tmp_next_orientation;
    Nucleotide rev_base;

   

    if (db_node_edge_exist(node,base,orientation)){

      tmp_next_node = db_graph_get_next_node(node,orientation,&tmp_next_orientation,
					 base,&rev_base,db_graph);

      if (db_node_check_status(tmp_next_node,status)){
	next_base[count] = base;
	next_node[count] = tmp_next_node;
	
	next_orientation[count] = tmp_next_orientation;
	count++;
      }
      
    }
    
  }
  
  
  nucleotide_iterator(&check_edge);

  if (count == n){
    ret = true;
  }

  return ret;

}



int db_graph_db_node_clip_tip_with_orientation(dBNode * node, Orientation orientation, int limit,
					       void (*node_action)(dBNode * node),dBGraph * db_graph){ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  int i;
  dBNode * nodes[limit];
  Orientation next_orientation;
  dBNode * next_node;
  char seq[db_graph->kmer_size+1];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end(node, opposite_orientation(orientation))){
    
    boolean join_main_trunk = false;

     while(db_node_has_precisely_one_edge(node,orientation,&nucleotide)) {
  
       nodes[length] = node;

       next_node = db_graph_get_next_node(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
       
       if(next_node == NULL){
	 printf("dB_graph_db_node_clip_tip_with_orientation: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
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
	 
	 node_action(nodes[i]);
	 //perhaps we want to move this inside the node action?
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


// clip a tip in the graph (the tip starts in node)
// limit is max length for tip
// node_action is applied to all the elements in the tip
// returns the length of the tip (0 means no length)
int db_graph_db_node_clip_tip(dBNode * node, int limit,
			      void (*node_action)(dBNode * node),
			      dBGraph * db_graph){

  int length_tip = 0;

  
  length_tip = db_graph_db_node_clip_tip_with_orientation(node,forward,limit,node_action,db_graph);
  
  if (length_tip==0){
    length_tip = db_graph_db_node_clip_tip_with_orientation(node,reverse,limit,node_action,db_graph);
    
  }
  
  return length_tip;
}


//remove any node under x coverage (and the arrows)

boolean db_graph_db_node_prune_low_coverage(dBNode * node, int coverage,
			   void (*node_action)(dBNode * node),
			   dBGraph * db_graph){

  boolean ret = false;

  if (element_get_coverage(node)<=coverage){
    ret = true;
  
    void nucleotide_action(Nucleotide nucleotide){
      Orientation next_orientation;
      Nucleotide reverse_nucleotide;
      dBNode * next_node;
      
      if (db_node_edge_exist(node,nucleotide,forward)){
	next_node = db_graph_get_next_node(node,forward,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
	db_node_reset_edge(next_node,opposite_orientation(next_orientation),reverse_nucleotide);
	  
      }
	
      if (db_node_edge_exist(node,nucleotide,reverse)){
	next_node = db_graph_get_next_node(node,reverse,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
	db_node_reset_edge(next_node,opposite_orientation(next_orientation),reverse_nucleotide);
	
      }
    }

    nucleotide_iterator(&nucleotide_action);

    node_action(node);
    db_node_reset_edges(node);	
   

  }

  return ret;
}


// computes a perfect path starting from a node and an edge
// ie the starting node can have  multiple exits
// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only
// returns length of the path (as number of labels)

int db_graph_get_perfect_path_with_first_edge(dBNode * node, Orientation orientation, int limit, 
					      Nucleotide fst_nucleotide,
					      void (*node_action)(dBNode * node),
					      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
					      char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
					      boolean * is_cycle, dBGraph * db_graph){

  Orientation  current_orientation,next_orientation;
  dBNode * current_node = NULL;
  dBNode * next_node = NULL;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  int sum_coverage = 0;
  int coverage  = 0;

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
  *max_coverage         = 0;
  *min_coverage         = INT_MAX;
 
  if (DEBUG){
    printf("\nNode %i in path: %s\n", length, binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
   }
    
  
  //first edge defined
  nucleotide = fst_nucleotide;

  do{ 
    if (length>0){
      node_action(current_node);
      sum_coverage += coverage;
      *max_coverage = *max_coverage < coverage ? coverage : *max_coverage;
      *min_coverage = *min_coverage > coverage ? coverage : *min_coverage;
    }

    next_node =  db_graph_get_next_node(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph);

      
    //sanity check
    if(next_node == NULL){
      fprintf(stderr,"dB_graph: didnt find node in hash table: %s %c %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),binary_nucleotide_to_char(nucleotide), current_orientation == forward ? "forward" : "reverse");
      exit(1);
    }

    
    path_labels[length]        = nucleotide;
    seq[length]                = binary_nucleotide_to_char(nucleotide);
    coverage                   = element_get_coverage(next_node);
    
    length++;

    //printf("perfect_path_with_edge %i\n",length);
    path_nodes[length]         = next_node;
    path_orientations[length]  = next_orientation;
    
    if (DEBUG){
      printf("\nNode %i in path: %s\n", length, binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq));
    }
    
    current_node        = next_node;
    current_orientation = next_orientation;
    
  } while (length<limit && 
	   !((next_node == node) && (next_orientation == orientation)) && //loop
	   db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2) && //multiple entries
	   db_node_has_precisely_one_edge(current_node, current_orientation,&nucleotide)); //has one next edge only
  
  
  if ((next_node == node) && (next_orientation == orientation)){
    *is_cycle = true;
  }
  
  
  
  if (DEBUG){
    printf("\nLast node in path: %s %i length: %i\n", binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq),next_node->edges,length);
  }

  seq[length] = '\0';
  *avg_coverage = (length-1<=0) ? 0 : (double) sum_coverage/(double) (length-1);

  if (*min_coverage == INT_MAX)
    {
      *min_coverage = 0;
    };

  return length;
  
}


// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only

int db_graph_get_perfect_path(dBNode * node, Orientation orientation, int limit, 
			      void (*node_action)(dBNode * node),
			      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
			      char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
			      boolean * is_cycle, dBGraph * db_graph)
{

  int length =0;
  Nucleotide nucleotide;

  //sanity check
  if (node == NULL){
    printf("db_graph_get_perfect_path: can't pass a null node\n");
    exit(1);
  }

  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  
 

  if (db_node_has_precisely_one_edge(node,orientation,&nucleotide))
    {
    
      length= db_graph_get_perfect_path_with_first_edge(node,orientation, limit, nucleotide,
							node_action,
							path_nodes,path_orientations,path_labels,
							seq,avg_coverage,min_coverage,max_coverage,
							is_cycle,db_graph);
    }
  else{
    *max_coverage         = 0;
    *min_coverage         = 0;
    seq[0] = '\0';
  }

  return length;

}


// a bubble starts in a node with only two outgoing edges in the same orientation
// every branch of the bubble is free of conflicts until it joins with the main branch
// the algorithms works by trying to find supernodes for all possible edges in the given orientation.
// the action is applied to all the internal nodes in the supernodes (even if they don't form a bubble).
// at the moment returns only one bubble, notice that potentially there could be situations where a node forms two bubbles (something to implement). 


//TODO - could have triallelic
boolean db_graph_detect_bubble(dBNode * node,
			       Orientation orientation,
			       int limit,
			       void (*node_action)(dBNode * node), 
			       int * length1,dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,
			       char * seq1, double * avg_coverage1, int * min_coverage1, int * max_coverage1,
			       int * length2,dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,
			       char * seq2, double * avg_coverage2, int * min_coverage2, int * max_coverage2,
			       dBGraph * db_graph){

  
  
  boolean ret = false;
  
  dBNode * path_nodes[4][limit+1];
  Orientation path_orientations[4][limit+1];
  Nucleotide path_labels[4][limit];
  char seq[4][limit+1];
  double avg_coverage[4];
  int min_coverage[4];
  int max_coverage[4];
  int lengths[4];
  int i=0;
  
  void check_nucleotide(Nucleotide n){

    boolean is_cycle;      
    if (db_node_edge_exist(node,n,orientation)){
	lengths[i] = db_graph_get_perfect_path_with_first_edge(node,orientation,
							       limit,n,
							       node_action,
							       path_nodes[i],path_orientations[i],path_labels[i],
							       seq[i],&avg_coverage[i],&min_coverage[i],&max_coverage[i],
							       &is_cycle,db_graph);
	
	i++;
    }
  }
  
  
  nucleotide_iterator(check_nucleotide);
 
  int j = 0;
  int k = 0;

  while (!ret && j<i){
    k = j+1;
    while (!ret && k<i){
      
      if (path_nodes[j][lengths[j]] == path_nodes[k][lengths[k]]){
	*length1 = lengths[j];
	*length2 = lengths[k];

	*avg_coverage1 = avg_coverage[j];
	*avg_coverage2 = avg_coverage[k];
	
	*min_coverage1 = min_coverage[j];
	*min_coverage2 = min_coverage[k];
		
	*max_coverage1 = max_coverage[j];
	*max_coverage2 = max_coverage[k];
	  
	int l;
	for(l=0;l<lengths[j]; l++){
	  path_nodes1[l] = path_nodes[j][l];
	  path_orientations1[l] = path_orientations[j][l];
	  path_labels1[l] = path_labels[j][l];
	  seq1[l] = seq[j][l];
	}
	path_nodes1[lengths[j]] = path_nodes[j][lengths[j]];
	path_orientations1[lengths[j]] = path_orientations[j][lengths[j]];
	seq1[lengths[j]] = seq[j][lengths[j]];

	for(l=0;l<lengths[k]; l++){
	  path_nodes2[l] = path_nodes[k][l];
	  path_orientations2[l] = path_orientations[k][l];
	  path_labels2[l] = path_labels[k][l];
	  seq2[l] = seq[k][l];
	}
	path_nodes2[lengths[k]] = path_nodes[k][lengths[k]];
	path_orientations2[lengths[k]] = path_orientations[k][lengths[k]];
	seq2[lengths[k]] = seq[k][lengths[k]];
	
	ret = true;
	  
      }
      else{
	k++;
      }
      
    }
    j++;
  }
    
  return ret;
}



//clip the branch with smaller coverage in the first kmer --> this can be more sophisticated
//dont want to flatten repeats -- use coverage_limit for this
boolean db_graph_db_node_smooth_bubble(dBNode * node, Orientation orientation, 
				       int limit,int coverage_limit,
				       void (*node_action)(dBNode * node),
				       dBGraph * db_graph){

  boolean ret = false;
  int length1, length2;
  dBNode * path_nodes1[limit+1];
  dBNode * path_nodes2[limit+1];
  Orientation path_orientations1[limit+1];
  Orientation path_orientations2[limit+1];
  Nucleotide path_labels1[limit];
  Nucleotide path_labels2[limit];
  
  dBNode * * path_nodes_tmp;
  Orientation * path_orientations_tmp;
  Nucleotide * path_labels_tmp;

  char seq1[limit+1];
  double avg_coverage1;
  int min_coverage1,max_coverage1;

  char seq2[limit+1];
  double avg_coverage2;
  int min_coverage2,max_coverage2;


  int length_tmp;
  Nucleotide reverse_nucleotide;
  Orientation next_orientation;


  if (db_graph_detect_bubble(node,orientation,limit,&db_node_action_do_nothing,
			     &length1,path_nodes1,path_orientations1,path_labels1,
			     seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
			     &length2,path_nodes2,path_orientations2,path_labels2,
			     seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
			     db_graph)){
 

    if ((double)coverage_limit>=avg_coverage1 && 
        (double)coverage_limit>=avg_coverage2){

     
      //prune
      int i;
     
      ret = true;
      if (avg_coverage1<avg_coverage2){
	path_nodes_tmp = path_nodes1;
	path_orientations_tmp = path_orientations1;
	path_labels_tmp  = path_labels1;
	length_tmp     = length1;
	db_node_reset_edge(node,orientation,path_labels1[0]);
      }
      else
	{
	  path_nodes_tmp = path_nodes2;
	  length_tmp     = length2;
	  path_orientations_tmp = path_orientations2;
	  path_labels_tmp       = path_labels2;
	  db_node_reset_edge(node,orientation,path_labels2[0]);
	}
      
      for(i=1;i<length_tmp;i++){
	node_action(path_nodes_tmp[i]);
	db_node_reset_edges(path_nodes_tmp[i]);
      }
      
      db_graph_get_next_node(path_nodes_tmp[length_tmp-1],path_orientations_tmp[length_tmp-1],&next_orientation,path_labels_tmp[length_tmp-1],&reverse_nucleotide,db_graph);
      db_node_reset_edge(path_nodes_tmp[length_tmp],opposite_orientation(next_orientation),reverse_nucleotide);
    }
    
  }

  return ret;
}


// it returns the supernode containing 'node'  
// string has to support limit+1 (+1 as you need a space for the \0 at the end)
// node_action has to be idempotent as it can be applied to the same node twice!!
// supernode_str returns the string made of the labels of the path (doesn't include first kmer). 
// returns the length of supernode 


int db_graph_supernode(dBNode * node,int limit,void (*node_action)(dBNode * node), 
		       dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
		       char * supernode_str, double * avg_coverage,int * min,int * max, boolean * is_cycle, 
		       dBGraph * db_graph){

  
  //use the allocated space as a temporary space
  dBNode * * nodes_reverse = path_nodes;
  Orientation * orientations_reverse = path_orientations;
  Nucleotide * labels_reverse = path_labels;
     
  boolean is_cycler;
  int length_reverse;
  int length = 0;
  
  int minr,maxr;
  double avg_coverager;


  //compute the reverse path until the end of the supernode
  //return is_cycle_reverse == true if the path closes a loop    
  
  length_reverse = db_graph_get_perfect_path(node,reverse,limit,&db_node_action_do_nothing,
					     nodes_reverse,orientations_reverse,labels_reverse,
					     supernode_str,&avg_coverager,&minr,&maxr,
					     &is_cycler,db_graph);
  
  
  
  if (length_reverse>0){
    //let's re do the last step, we need to do that because the last node could have had multiple entries
    
    Nucleotide label;
    Orientation next_orientation;
    
    dBNode * lst_node = db_graph_get_next_node(nodes_reverse[length_reverse-1],orientations_reverse[length_reverse-1],
					       &next_orientation, labels_reverse[length_reverse-1],&label,db_graph);
    
    //sanity check
    if (lst_node != nodes_reverse[length_reverse]){
      puts("db_graph_supernode broken!\n");
      exit(1);
    }
    
    
    length = db_graph_get_perfect_path_with_first_edge(nodes_reverse[length_reverse],
						       opposite_orientation(orientations_reverse[length_reverse]),
						       limit,label,
						       node_action,
						       path_nodes,path_orientations,path_labels,
						       supernode_str,avg_coverage,min,max,
						       is_cycle,db_graph);
    
  }
  else{
    length = db_graph_get_perfect_path(node,forward,
				       limit,
				       node_action,
				       path_nodes,path_orientations,path_labels,
				       supernode_str,avg_coverage,min,max,
				       is_cycle,db_graph);
  }
  
  
  

  //apply action to the fst and last node
  node_action(path_nodes[0]);
  node_action(path_nodes[length]);
  
  return length;
}



// The idea is this. First of all you check the node you are given with condition_on_initial_node_only.
// If this is false, forget the whole supernode - ignore it. In that case *path_length==0.
// If condition_on_initial_node_only is true, then apply big OR / big AND of condition_for_all_nodes (depending on flag)
// apply the action to all nodes in the supernode, 
// whether true or false, will return the supernode in path_nodes, path_orientations, path_labels
// node_action MUST BE IDEMPOTENT

boolean db_graph_logical_operation_on_condition_for_all_nodes_in_supernode(dBNode * node,int limit, 
									   boolean (*condition_for_all_nodes)(dBNode * node),  
									   void (*node_action)(dBNode * node),
									   dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, int* path_length,
									   char * string, double * avg_coverage,int * min,int * max, boolean * is_cycle, 
									   boolean big_and, //otherwise is big or
									   dBGraph * db_graph)
{

  boolean return_val;
  if (big_and)
    {
      return_val = true;
    }
  else
    {
      return_val=false;
    }

  void action_check_condition_on_all_nodes(dBNode* n)
    {
      if ( big_and == true)
	{
	  return_val = return_val && condition_for_all_nodes(n);
	}
      else
	{
	  return_val = return_val || condition_for_all_nodes(n);
	}
    }

       
  *path_length = db_graph_supernode(node, limit, 
				    action_check_condition_on_all_nodes, 
				    path_nodes, path_orientations, path_labels,
				    string, avg_coverage, min, max, is_cycle,
				    db_graph);

  
  //now apply action to all nodes
  //notice that a node could be met twice with different orientations, that's why we don't apply the action in the function above (passed to db_graph_supernode) as it could interact with the condition we are checking for. 
  int i;
  for (i=0; i<= *path_length; i++)
    {
      node_action(path_nodes[i]);
    }
  
  return return_val;
  
}


// The idea is this. First of all you check the node you are given with condition_on_initial_node_only. If this is false, forget the whole supernode - ignore it. In that case *path_length==0.
// If condition_on_initial_node_only  is true, then you get the supernode and apply the action to all nodes in the supernode, and also you check
// if condition_for_all_nodes is true for all nodes, and return true, else false.
// whether true or false, will return the supernode in path_nodes, path_orientations, path_labels
// node_action MUST BE IDEMPOTENT
boolean db_graph_is_condition_true_for_all_nodes_in_supernode(dBNode * node,int limit, 						
							      boolean (*condition_for_all_nodes)(dBNode * node),  
							      void (*node_action)(dBNode * node),
							      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, int* path_length,
							      char * string, double * avg_coverage,int * min,int * max, boolean * is_cycle, 
							      dBGraph * db_graph)
{

  return db_graph_logical_operation_on_condition_for_all_nodes_in_supernode(node,limit,
									    condition_for_all_nodes,
									    node_action,
									    path_nodes, path_orientations, path_labels, path_length, 
									    string, avg_coverage, min, max, is_cycle,
									    true, db_graph);
									   
}


// The idea is this. First of all you check the node you are given with condition_on_initial_node_only.
// Initial node dependes on the order of traversal (we cannot determine that beforehand, depends on hash_table size/parameters)
// If this is false, forget the whole supernode - ignore it. In that case *path_length==0.
// If condition_on_initial_node_only  is true, then you get the supernode and apply the action to all nodes in the supernode, and also you check
// if condition_for_all_nodes is true for any of the nodes, and return true, else false.
// whether true or false, will return the supernode in path_nodes, path_orientations, path_labels
// node_action MUST BE IDEMPOTENT



boolean db_graph_is_condition_true_for_at_least_one_node_in_supernode(dBNode * node,int limit, 
								      boolean (*condition_for_all_nodes)(dBNode * node),  
								      void (*node_action)(dBNode * node),
								      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,  int* path_length,
								      char * string, double * avg_coverage,int * min,int * max, boolean * is_cycle, 
								      dBGraph * db_graph)
{
  return db_graph_logical_operation_on_condition_for_all_nodes_in_supernode(node,limit,
									    condition_for_all_nodes,
									    node_action,
									    path_nodes, path_orientations, path_labels, path_length,
									    string, avg_coverage, min, max, is_cycle,
									    false, db_graph);

}


// The idea is this. First of all you check the node you are given with condition_on_initial_node_only. 
// If this is false, forget the whole supernode - ignore it. In that case *path_length==0.
// If condition_on_initial_node_only  is true, then you get the supernode;
// if condition_for_all_nodes is true for >= min_start nodes at the start and >=min_end at the end of the supernode, but NOT ALL, then return true, else false.
// whether true or false, will return the supernode in path_nodes, path_orientations, path_labels
// apply the action to all nodes in the supernode, and also you check
// node_action MUST BE IDEMPOTENT

boolean db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(dBNode * node,int limit,
										    boolean (*condition_for_all_nodes)(dBNode * node),  
										    void (*node_action)(dBNode * node), 
										    int min_start, int min_end, int min_diff,
										    dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,int * path_length,
										    char * string, double * avg_coverage,int * min,int * max, boolean * is_cycle,
										    dBGraph * db_graph)
{

  boolean condition_is_true_for_all_nodes_in_supernode=true;

  void action_check_condition_on_all_nodes(dBNode* n)
    {
      //printf("db_graph_is_condition_true..:  applying action to node %s\n", binary_kmer_to_seq(element_get_kmer(n),db_graph->kmer_size,tmp_seq) );
      condition_is_true_for_all_nodes_in_supernode = condition_is_true_for_all_nodes_in_supernode  &&   condition_for_all_nodes(n);
    }

  
  *path_length = db_graph_supernode(node, limit, 
				    action_check_condition_on_all_nodes, 
				    path_nodes, path_orientations, path_labels,
				    string, avg_coverage, min, max, is_cycle,
				    db_graph);
  
  int num_nodes_at_start_where_condition_is_true=0;
  int num_nodes_at_end_where_condition_is_true=0;
  
  if (*path_length>0)
    {
	   //now see if have enoug nodes at start and end with condition true.
      if (!condition_is_true_for_all_nodes_in_supernode)
	{
	  int j;
	  boolean flag=true;
	  for (j=0; j<=*path_length; j++)
	    {
	      if ( (flag==true) && (condition_for_all_nodes(path_nodes[j])==true) )
		{
		  num_nodes_at_start_where_condition_is_true++;
		}
	      else
		{
		  flag=false;
		}
	    }
	  
	  
	  flag=true;
	  for (j=*path_length; j>=0; j--)
	    {
	      if ( (flag==true) && (condition_for_all_nodes(path_nodes[j])==true) )
		{
		  num_nodes_at_end_where_condition_is_true++;
		}
	      else
		{
		  flag=false;
		}
	    }
	  
	}
    }

  //now apply action to all nodes
  int i;
  for (i=0; i<= *path_length; i++)
    {
      node_action(path_nodes[i]);
    }
  
  boolean ret_val = false;
  if ( (num_nodes_at_start_where_condition_is_true>=min_start) && (num_nodes_at_end_where_condition_is_true>=min_end) 
       && ((*path_length+1 - num_nodes_at_start_where_condition_is_true - num_nodes_at_end_where_condition_is_true)>min_diff) )
    {
      ret_val = true;
    }

  return ret_val;
  
}




void db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required, FILE* fout,  
										  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index)
{

  int count_nodes=0;
  void print_supernode(dBNode * e)
    {

      dBNode * nodes_path[5000];
      Orientation orientations_path[5000];
      Nucleotide labels_path[5000];
      char seq[5000+1];
      double avg_coverage;
      int min_coverage;
      int  max_coverage;
      int length_path=0;
      boolean is_cycle;
      
      if (db_node_check_status_is_not_visited_or_visited_and_exists_in_reference(e)){

	
	if (db_graph_is_condition_true_for_all_nodes_in_supernode(e, 5000, condition,
								  &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
								  nodes_path,orientations_path, labels_path, &length_path, 
								  seq, &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
								  db_graph))
	  {
	
	   

	    if (min_coverage>=min_covg_required)
	      {
		if (!is_for_testing)
		  {
		    if (fout == NULL)
		      {
			printf(">node_%i length: %i min covg: %d max covg: %d\n",count_nodes,length_path+db_graph->kmer_size, min_coverage, max_coverage);
			count_nodes++;
			printf("%s\n",seq);
		      }
		    else
		      {
			fprintf(fout, ">node_%i length: %i min covg: %d max covg: %d\n",count_nodes,length_path+db_graph->kmer_size, min_coverage, max_coverage);
			count_nodes++;
			fprintf(fout, "%s\n",seq);
		      }
		  }
		else
		  {
		    //return the supernode in the preallocated array, so the test can check it.
		    for_test_array_of_supernodes[*for_test_index][0]='\0';
		    strcat(for_test_array_of_supernodes[*for_test_index], seq);
		    *for_test_index=*for_test_index+1;
		  }
	      }
	  }
      }
    }
  hash_table_traverse(&print_supernode,db_graph); 
}

/*
//not yet implemented
void db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode_AND_print_supernodes_it_connects_to(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required, FILE* fout,
														      boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index)
{
  exit(1);//remove when implemented


  int count_nodes=0;
  void print_supernode_and_neighbours(dBNode * e)
    {

      dBNode * nodes_path[5000];
      Orientation orientations_path[5000];
      Nucleotide labels_path[5000];
      char seq[5000+db_graph->kmer_size+1];
      int length_path=0;
      
      if (db_graph_is_condition_true_for_all_nodes_in_supernode(e, 5000, &db_node_check_status_is_not_visited_or_visited_and_exists_in_reference, condition,
								&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
								seq,nodes_path,orientations_path, labels_path, &length_path, db_graph))
	{
	  if (length_path>0) // db_graph_is_condition_true_for_all_nodes_in_supernode returns zero length if the condition for initial node only is false. Usually this condition is checking if visited
	    {
	      //work out min and max covg for this supernode
	      int min=element_get_coverage(nodes_path[0]);
	      int max=element_get_coverage(nodes_path[0]);
	      int j;
	      for (j=0; j<= length_path; j++)
		{
		  int cov = element_get_coverage(nodes_path[j]);
		  if (cov<min)
		    {
		      min=cov;
		    }
		  if (cov>max)
		    {
		      max=cov;
		    }
		}
	      
	      if (min>=min_covg_required)
		{
		  if (!is_for_testing)
		    {
		      if (fout == NULL)
			{
			  printf(">node_%i length: %i min covg: %d max covg: %d\n",count_nodes,length_path+db_graph->kmer_size, min, max);
			  count_nodes++;
			  printf("%s\n",seq);

			  //Now print out the supernodes BEFORE this supernode, which connect to nodes_path[0]

			  
			  //Now print out the supernodes AFTER this supernode, which nodes_path[length+1] connects to  ..is that right??

			}
		      else
			{
			  fprintf(fout, ">node_%i length: %i min covg: %d max covg: %d\n",count_nodes,length_path+db_graph->kmer_size, min, max);
                          count_nodes++;
                          fprintf(fout, "%s\n",seq);
			}
		    }
		  else
		    {
		      //return the supernode in the preallocated array, so the test can check it.
		      for_test_array_of_supernodes[*for_test_index][0]='\0';
		      strcat(for_test_array_of_supernodes[*for_test_index], seq);
		      *for_test_index=*for_test_index+1;
		    }
		}
	    }
	}
    }
  hash_table_traverse(&print_supernode,db_graph); 
  
}

*/



void db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required, FILE* fout, 
										  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index)
{
  
  int count_nodes=0;
  void print_supernode(dBNode * e)
  {
    
    dBNode * nodes_path[5000];
    Orientation orientations_path[5000];
    Nucleotide labels_path[5000];
    char seq[5000+1];
    int length_path=0;
    double avg_coverage;
    int min_coverage;
    int  max_coverage;
    boolean is_cycle;
    
    
    if (db_node_check_status_is_not_visited_or_visited_and_exists_in_reference(e))
      {
    
	if (db_graph_is_condition_true_for_at_least_one_node_in_supernode(e, 5000, condition,
									  &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
									  nodes_path,orientations_path, labels_path, &length_path, 
									  seq, &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
									  db_graph))
	{

	  
	  if (min_coverage>=min_covg_required)
	    {
	      if (!is_for_testing)
		{
		  if (fout==NULL)
		    {
		      printf(">node_%i length: %i min covg: %d max covg: %d\n",count_nodes,length_path+db_graph->kmer_size, min_coverage, max_coverage);
		      count_nodes++;
		      printf("%s\n",seq);
		    }
		  else
		    {
		      fprintf(fout, ">node_%i length: %i min covg: %d max covg: %d\n",count_nodes,length_path+db_graph->kmer_size, min_coverage, max_coverage);
		      count_nodes++;
		      fprintf(fout, "%s\n",seq);
		    }
	      }
	      else
		{
		  //return the supernode in the preallocated array, so the test can check it.
		  for_test_array_of_supernodes[*for_test_index][0]='\0';
		  strcat(for_test_array_of_supernodes[*for_test_index], seq);
		  *for_test_index=*for_test_index+1;
		}
	    }
	}
      }
  }
  hash_table_traverse(&print_supernode,db_graph); 
}




// can use this toprint potential indels.
// min_start and min_end are the number of nodes of overlap with the reference that  you want at the start and  end of the supernode
void db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required,
												       int min_start, int min_end, int min_diff, FILE* fout,
												       boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index)
{

  int count_nodes=0;
  void print_supernode(dBNode * e)
    {

      dBNode * nodes_path[5000];
      Orientation orientations_path[5000];
      Nucleotide labels_path[5000];
      char seq[5000+1];
      double avg_coverage;
      int min_coverage;
      int  max_coverage;
      int length_path=0;
      boolean is_cycle;

      if (db_node_check_status_is_not_visited_or_visited_and_exists_in_reference(e))
	{
	  if (db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(e, 5000, condition,
											  &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
											  min_start, min_end,min_diff,
											  nodes_path,orientations_path, labels_path, &length_path,
											  seq, &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
											  db_graph))
	{
	  
	  
	      if (min_coverage>=min_covg_required)
		{
		  
		  if (!is_for_testing)
		    {
		      
		      if (fout == NULL)
			{
			  printf(">potential_sv_node_%i length: %i min covg: %d max covg: %d\n",count_nodes,length_path+db_graph->kmer_size, min_coverage, max_coverage);
			  count_nodes++;
			  printf("%s\n",seq);
			}
		      else
			{
			  fprintf(fout, ">potential_sv_node_%i length: %i min covg: %d max covg: %d\n",count_nodes,length_path+db_graph->kmer_size, min_coverage, max_coverage);
                          count_nodes++;
                          fprintf(fout, "%s\n",seq);

			}
		    }
		  else
		    {
		      //return the supernode in the preallocated array, so the test can check it.
		      for_test_array_of_supernodes[*for_test_index][0]='\0';
		      strcat(for_test_array_of_supernodes[*for_test_index], seq);
		      *for_test_index=*for_test_index+1;
		    }
		}
	      else
		{
		  //printf("Too low covg of only %d and we require %d\n",min, min_covg_required); 
		}
	}
	}
    }
  hash_table_traverse(&print_supernode,db_graph); 
}


//routine to get short indels and SNPs (delta = 0 necessary condition, not sufficient as we could also find balancing duplications or similar)

void db_graph_detect_vars(int delta, int max_length, dBGraph * db_graph){
  
  int count_vars = 0; 
	
  int flanking_length = 1000; 

  
  void get_vars(dBNode * node){
   
    void get_vars_with_orientation(dBNode * node, Orientation orientation){
      
      int length1, length2;
      dBNode * path_nodes1[max_length+1];
      dBNode * path_nodes2[max_length+1];
      Orientation path_orientations1[max_length+1];
      Orientation path_orientations2[max_length+1];
      Nucleotide path_labels1[max_length];
      Nucleotide path_labels2[max_length];
      char seq1[max_length+1];
      char seq2[max_length+1];
    
      double avg_coverage1;
      int min_coverage1,max_coverage1;
      
      double avg_coverage2;
      int min_coverage2,max_coverage2;
      
      dBNode * current_node = node;
   
      do{
	
	if (db_graph_detect_bubble(current_node,orientation,max_length,&db_node_action_set_status_visited,
				   &length1,path_nodes1,path_orientations1,path_labels1,
				   seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
				   &length2,path_nodes2,path_orientations2,path_labels2,
				   seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
				   db_graph)){
	  int length_flank5p = 0;	
	  int length_flank3p = 0;
	  dBNode * nodes5p[flanking_length];
	  dBNode * nodes3p[flanking_length];
	  Orientation orientations5p[flanking_length];
	  Orientation orientations3p[flanking_length];
	  Nucleotide labels_flank5p[flanking_length];
	  Nucleotide labels_flank3p[flanking_length];
	  char seq5p[flanking_length+1];
	  char seq3p[flanking_length+1];
	  boolean is_cycle5p, is_cycle3p;
	  
	  double avg_coverage5p;
	  int min5p,max5p;
	  double avg_coverage3p;

	  int min3p,max3p;
	  
	  printf("\nVARIATION: %i - coverage: %d\n",count_vars,element_get_coverage(current_node));
	  count_vars++;
	  
	  //compute 5' flanking region       	
	  int length_flank5p_reverse = db_graph_get_perfect_path(current_node,opposite_orientation(orientation),
								 flanking_length,&db_node_action_set_status_visited,
								 nodes5p,orientations5p,labels_flank5p,
								 seq5p,&avg_coverage5p,&min5p,&max5p,
								 &is_cycle5p,db_graph);
	  
	  if (length_flank5p_reverse>0){
	    Nucleotide label;
	    Orientation next_orientation;
	    
	    dBNode * lst_node = db_graph_get_next_node(nodes5p[length_flank5p_reverse-1],orientations5p[length_flank5p_reverse-1],
						       &next_orientation, labels_flank5p[length_flank5p_reverse-1],&label,db_graph);
	    	    
	    length_flank5p = db_graph_get_perfect_path_with_first_edge(nodes5p[length_flank5p_reverse],
								       opposite_orientation(orientations5p[length_flank5p_reverse]),
								       flanking_length,label,
								       &db_node_action_set_status_visited,
								       nodes5p,orientations5p,labels_flank5p,
								       seq5p,&avg_coverage5p,&min5p,&max5p,
								       &is_cycle5p,db_graph);
	  }
	  else{
	    length_flank5p = 0;
	  }

	  printf("length 5p flank: %i avg_coverage:%5.2f \n",length_flank5p,avg_coverage5p);
	  
	  
	  //compute 3' flanking region
	  length_flank3p = db_graph_get_perfect_path(path_nodes2[length2],path_orientations2[length2],
						     flanking_length,&db_node_action_set_status_visited,
						     nodes3p,orientations3p,labels_flank3p,
						     seq3p,&avg_coverage3p,&min3p,&max3p,
						     &is_cycle3p,db_graph);
	  
	  printf("length 3p flank: %i avg_coverage:%5.2f\n",length_flank3p,avg_coverage3p);
	  
	  printf("length branch 1: %i - avg_coverage: %5.2f\n",length1,avg_coverage1);
	  printf("length branch 2: %i - avg_coverage: %5.2f\n",length2,avg_coverage2);
	  
	  
	  char name[100];
	  
	  //print flank5p - 
	  sprintf(name,"var_5p_flank_%i",count_vars);
	  print_fasta_from_path(stdout,name,length_flank5p,avg_coverage5p,min5p,max5p,
				nodes5p[0],orientations5p[0],				
				nodes5p[length_flank5p],orientations5p[length_flank5p],				
				seq5p,
				db_graph->kmer_size,false);	
	  
	  //print branches
	  sprintf(name,"branch_%i_1",count_vars);
	  print_fasta_from_path(stdout,name,length1,
				avg_coverage1,min_coverage1,max_coverage1,
				path_nodes1[0],path_orientations1[0],path_nodes1[length1],path_orientations1[length1],
				seq1,db_graph->kmer_size,false);
	  
	  sprintf(name,"branch_%i_2",count_vars);
	  print_fasta_from_path(stdout,name,length2,
				avg_coverage2,min_coverage2,max_coverage2,
				path_nodes2[0],path_orientations2[0],path_nodes2[length2],path_orientations2[length2],
				seq2,db_graph->kmer_size,false);
	  
	  //print flank3p
	  sprintf(name,"var_3p_flank_%i",count_vars);
	  print_fasta_from_path(stdout,name,length_flank3p,avg_coverage3p,min3p,max3p,
				nodes3p[0],orientations3p[0],nodes3p[length_flank3p],orientations3p[length_flank3p],
				seq3p,db_graph->kmer_size,false);
	  

	  db_node_action_set_status_visited(path_nodes2[length2]);
	  
	  current_node = path_nodes2[length2];
	  orientation = path_orientations2[length2];
	}
      } while (current_node != node //to avoid cycles
	       && db_node_check_status_none(current_node));
    }
    
  
    if (db_node_check_status_none(node)){     
      db_node_action_set_status_visited(node);
      get_vars_with_orientation(node,forward);
      get_vars_with_orientation(node,reverse);
    }
  }
 
  hash_table_traverse(&get_vars,db_graph); 
    
}



void db_graph_smooth_bubbles(int coverage,int limit,dBGraph * db_graph){
  
  void smooth_bubble(dBNode * node){
    if (db_node_check_status_none(node)){
      db_graph_db_node_smooth_bubble(node,forward,limit,coverage,
				     &db_node_action_set_status_pruned,db_graph);
      db_graph_db_node_smooth_bubble(node,reverse,limit,coverage,
				     &db_node_action_set_status_pruned,db_graph);
      
    }
  }
  
  hash_table_traverse(&smooth_bubble,db_graph); 

}




void db_graph_clip_tips(dBGraph * db_graph){
  
  void clip_tips(dBNode * node){
    
    if (db_node_check_status_none(node)){
      db_graph_db_node_clip_tip(node, db_graph->kmer_size*2,&db_node_action_set_status_pruned,db_graph);
    }
  }

  hash_table_traverse(&clip_tips,db_graph);
  
}


void db_graph_print_supernodes(char * filename, int max_length, dBGraph * db_graph){

  FILE * fout; //binary output
  fout= fopen(filename, "w"); 

  int count_nodes=0;
  
  dBNode * *    path_nodes;
  Orientation * path_orientations;
  Nucleotide *  path_labels;
  char * seq;
  boolean is_cycle;
  double avg_coverage;
  int min,max;
  
  
  path_nodes        = calloc(max_length,sizeof(dBNode*));
  path_orientations = calloc(max_length,sizeof(Orientation));
  path_labels       = calloc(max_length,sizeof(Nucleotide));
  seq               = calloc(max_length+1+db_graph->kmer_size,sizeof(char));
  
  

  long long count_kmers = 0;
  long long count_sing  = 0;

  void print_supernode(dBNode * node){
    
   

    count_kmers++;
    
    char name[100];

    if (db_node_check_status_none(node) == true){
      int length = db_graph_supernode(node,max_length,&db_node_action_set_status_visited,
				      path_nodes,path_orientations,path_labels,
				      seq,&avg_coverage,&min,&max,&is_cycle,
				      db_graph);
    
 

      if (length>0){	
	sprintf(name,"node_%i",count_nodes);

	print_fasta_from_path(fout,name,length,avg_coverage,min,max,path_nodes[0],path_orientations[0],path_nodes[length],path_orientations[length],seq,db_graph->kmer_size,true);
	if (length==max_length){
	  printf("contig length equals max length [%i] for node_%i\n",max_length,count_nodes);
	}
	count_nodes++;
      }
      else{
	sprintf(name,"node_%qd",count_sing);
	print_fasta_from_path(stdout,name,length,avg_coverage,min,max,path_nodes[0],path_orientations[0],path_nodes[length],path_orientations[length],seq,db_graph->kmer_size,true);
	
	count_sing++;
      }
    
    }
  }
  
  hash_table_traverse(&print_supernode,db_graph); 
  printf("%qd nodes visted [%qd singletons]\n",count_kmers,count_sing);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(seq);
  fclose(fout);
}





void db_graph_print_coverage(dBGraph * db_graph){
  long long count_kmers=0;

  void print_coverage(dBNode * node){
    char seq[db_graph->kmer_size];

    printf("%qd\t%s\t%i\n",count_kmers,binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq),element_get_coverage(node));
    count_kmers++;;
  }
  hash_table_traverse(&print_coverage,db_graph); 
}

void db_graph_remove_low_coverage_nodes(int coverage, dBGraph * db_graph){
  
  void prune_node(dBNode * node){
    db_graph_db_node_prune_low_coverage(node,coverage,
			   &db_node_action_set_status_pruned,
			   db_graph);
  }

  hash_table_traverse(&prune_node,db_graph); 
}




void db_graph_dump_binary(char * filename, boolean (*condition)(dBNode * node), dBGraph * db_graph){
  FILE * fout; //binary output
  fout= fopen(filename, "w"); 
  
  long long count=0;
  //routine to dump graph
  void print_node_binary(dBNode * node){   
    if (condition(node)){
      count++;
      db_node_print_binary(fout,node);
    }
  }

  hash_table_traverse(&print_node_binary,db_graph); 
  fclose(fout);

  printf("%qd kmers dumped\n",count);
}

void db_graph_traverse_with_array(void (*f)(HashTable*, Element *, int**, int),HashTable * hash_table, int** array, int length_of_array){

  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(hash_table, &hash_table->table[i], array, length_of_array);
    }
  }


}

dBNode* db_graph_get_first_node_in_supernode_containing_given_node(dBNode* node,  dBGraph* db_graph)
{
  char tmp_seq[db_graph->kmer_size];

  if (node==NULL)
    {
      printf("db_graph_get_first_node_in_supernode_containing_given_node: node is null so of course cant get first node");
      exit(1);
    }
    
  if (db_node_check_status(node, pruned) )
    {
      //don't waste time with pruned nodes.
      printf ("Warning. Getting first node in supernode of a pruned node");
      exit(1);
    }
  

  Nucleotide nucleotide1, nucleotide2, rev_nucleotide;
  Orientation original_orientation, next_orientation, orientation;
  dBNode * original_node=node;
  dBNode * next_node;

  //are we already at the end of a supernode
  if (db_node_is_supernode_end(original_node, forward))
    {
      return original_node;
    }
  else if (db_node_is_supernode_end(original_node, reverse))
    {
      return original_node;
    }

  //arbitrary choice - go in reverse direction until you reach the end of the supernode
  original_orientation = reverse; 
  orientation = reverse;

  while(db_node_has_precisely_one_edge(node,orientation,&nucleotide1))
  {
    next_node =  db_graph_get_next_node(node,orientation,&next_orientation,nucleotide1,&rev_nucleotide,db_graph);
    
    if(next_node == NULL){
      printf("dB_graph_get_first_node_in_supernode_containing_given_node: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size, tmp_seq));
      exit(1);
    }	         

    if (DEBUG)
      {
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide1),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size, tmp_seq) :  
	     binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size),db_graph->kmer_size, tmp_seq));
	     
      }
    

    //check for multiple entry edges 
    if (db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2))
      {
      }
    else
      {
	return node; //we have gone as far as we can go - the next node has multiple entries. So we are now at the first node of the supernode
      }
    
    
    //loop
    if ((next_node == original_node) && (next_orientation == original_orientation))
      {      
	return original_node; //we have a loop that returns to where we start. Might as well consider ourselves as at the fiurst node of the supernode right at the beginning
      }
    
    node = next_node;
    orientation = next_orientation;      
  }
  //printf("We have found the first node, it is %s\n", binary_kmer_to_seq(node->kmer, db_graph->kmer_size, tmp_seq));
  return node;

}

dBNode* db_graph_get_next_node_in_supernode(dBNode* node, Orientation orientation, Orientation* next_orientation,  dBGraph* db_graph)
{
  char tmp_seq[db_graph->kmer_size];

  if (node==NULL)
    {
      printf("do not try to get next node of a NULL node");
      exit(1);
    }
  else if (db_node_check_status(node, pruned))
    {
      printf("ignore pruned node");
      //don't waste time with pruned nodes.
      exit(1);
    }
  else if (db_node_is_supernode_end(node, orientation))
    {
      if (DEBUG)
	{
	  printf("this node is at the end of the supernode, in this orientation, so cant return the next one\n");
	}
      next_orientation=NULL;
      return NULL;
    }

  Nucleotide nucleotide_for_only_edge, reverse_nucleotide_for_only_edge;

  
  db_node_has_precisely_one_edge(node,orientation,&nucleotide_for_only_edge);//gives us nucleotide
  
  dBNode* next_node =  db_graph_get_next_node(node,orientation ,next_orientation,nucleotide_for_only_edge,&reverse_nucleotide_for_only_edge, db_graph );
 
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
  if (db_node_has_precisely_one_edge(next_node,opposite_orientation(*next_orientation),&nuc))
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


//the length of the supernode is passed to the caller in the array of supernode lengths that is passed in.
//the idea is that this function is passed into the traverse function, and the caller of that
void db_graph_get_supernode_length_marking_it_as_visited(dBGraph* db_graph, Element* node, int** array_of_supernode_lengths, int length_of_array)
{

  if (db_node_check_status(node, visited) || db_node_check_status(node,pruned))
    {
      return;
    }

  int length_of_supernode=1;
  dBNode* first_node=db_graph_get_first_node_in_supernode_containing_given_node(node, db_graph);
  dBNode* current_node;
  dBNode* next_node;
  current_node=first_node;
  Orientation start_orientation, current_orientation, next_orientation;
  db_node_set_status(current_node, visited);

  //work out which direction to leave supernode in. 
  if (db_node_is_supernode_end(first_node,forward))
    {
      if (db_node_is_supernode_end(first_node,reverse))
	{
	  //singleton
	  *(array_of_supernode_lengths[length_of_supernode])=*(array_of_supernode_lengths[length_of_supernode])+1;
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
  while (!db_node_is_supernode_end(current_node,current_orientation))
    {
      next_node = db_graph_get_next_node_in_supernode(current_node, current_orientation, &next_orientation, db_graph);

      if ((next_node==first_node) && (next_orientation==start_orientation))//back to the start - will loop forever if not careful ;-0
	{
	  break;
	}

      length_of_supernode++;
      current_node=next_node;
      current_orientation=next_orientation;
      db_node_set_status(current_node, visited);

    }


  if (length_of_supernode>length_of_array-1)
    {
      printf("We have a supernode of length %d, but have only allocated space to record lengths up to %d", length_of_supernode, length_of_array);
      exit(1);
    }
  *(array_of_supernode_lengths[length_of_supernode])=*(array_of_supernode_lengths[length_of_supernode])+1;
  return ;


  
}


int db_graph_get_N50_of_supernodes(dBGraph* db_graph)
{

  int numbers_of_supernodes_of_specific_lengths[10000]; //i-th entry is number of supernodes of length i.
  int* numbers_of_supernodes_of_specific_lengths_ptrs[10000];
  int lengths_of_supernodes[10000];
  
  int i;
  for (i=0; i<10000; i++)
    {
      numbers_of_supernodes_of_specific_lengths[i]=0;
      numbers_of_supernodes_of_specific_lengths_ptrs[i]=&(numbers_of_supernodes_of_specific_lengths[i]);
      lengths_of_supernodes[i]=0;

    }


  db_graph_traverse_with_array(&db_graph_get_supernode_length_marking_it_as_visited, db_graph,  numbers_of_supernodes_of_specific_lengths_ptrs, 10000);

  //we now have an array containing the number of supernodes of each length from 1 to 10,000 nodes.
  // for example it might be 0,100,1,0,0,0,23,4. note there are no supernodes of length 0, so first entry should always be 0

  //let's make another array, containing the lengths
  // this would , in the above example, be: 1,2,6,7 
  int ctr=0;

  for (i=0; i<10000; i++)
    {
      if (numbers_of_supernodes_of_specific_lengths[i]>0)
	{
	  lengths_of_supernodes[ctr]=i;
	  ctr++;
	}
    }

  //1 sort the list of lengths, to give in our example 0,0,....0,1,2,6,7
  size_t numbers_len=sizeof(numbers_of_supernodes_of_specific_lengths)/sizeof(int);
  qsort( lengths_of_supernodes, numbers_len, sizeof(int), int_cmp); 

  //2 add all of lengths*number of supernodes up - get total length
  int total=0;
  
  for(i=0; i<10000; i++)
    {
      total=total+i*numbers_of_supernodes_of_specific_lengths[i];
    }

  //3. Start at the max supernode length and work your way down until the length so far is half the total. Where you are is the n50
  if (ctr<0)
    {
      printf("negative ctr");
      exit(1);
    }


  int current_cumulative_len=0;
  //printf("Total length of supernodes is %d and ctr is %d\n", total, ctr);

  for (i=9999; i>=9999-ctr+1; i--)
    {
      current_cumulative_len += numbers_of_supernodes_of_specific_lengths[lengths_of_supernodes[i]] * lengths_of_supernodes[i] ;
      //printf("Looking at supernodes whose length is %d  - there are %d of them, i is %d, current cum length is %d\n", lengths_of_supernodes[i],  numbers_of_supernodes_of_specific_lengths[lengths_of_supernodes[i]], i, current_cumulative_len);

      if (current_cumulative_len >= total/2)
	{
	  //  printf("total is %d and n50 is %d\n", total, lengths_of_supernodes[i]);
	  // then we have found the N50.
	  return lengths_of_supernodes[i];
	}

    }
  
  //should not get here
  printf("Should  not get here. Exiting in N50 calc code");
  exit(1);
  return 1;
}

int int_cmp(const void *a, const void *b)
{
  const int *ia = (const int *)a; // casting pointer types
  const int *ib = (const int *)b;
  return *ia  - *ib; 
  /* integer comparison: returns negative if b > a 
     and positive if a > b */
}


void compute_label(dBNode * node, Orientation o, char * label){

  int i=0; 

  void check_nucleotide(Nucleotide n){

    if (db_node_edge_exist(node,n,o)){
	label[i] =  binary_nucleotide_to_char(n);
	i++;
      }
  }
  
  nucleotide_iterator(check_nucleotide);

  label[i]='\0';
}

void print_fasta_from_path(FILE *fout,
			   char * name,
			   int length,
			   double avg_coverage,
			   int min_coverage,
			   int max_coverage,
			   dBNode * fst_node,
			   Orientation fst_orientation,
			   dBNode * lst_node,
			   Orientation lst_orientation,
			   char * string, //labels of paths
			   int kmer_size,
			   boolean include_first_kmer){

    

  char fst_f[5], fst_r[5],lst_f[5],lst_r[5];

  compute_label(fst_node,forward,fst_f);
  compute_label(fst_node,reverse,fst_r);
  compute_label(lst_node,forward,lst_f);
  compute_label(lst_node,reverse,lst_r);

  char fst_seq[kmer_size+1], lst_seq[kmer_size+1];
 
  BinaryKmer fst_kmer = element_get_kmer(fst_node);
  if (fst_orientation==reverse){
    fst_kmer = binary_kmer_reverse_complement(fst_kmer,kmer_size);
  } 
  binary_kmer_to_seq(fst_kmer,kmer_size,fst_seq);

  BinaryKmer lst_kmer = element_get_kmer(lst_node);
  if (lst_orientation==reverse){
    lst_kmer = binary_kmer_reverse_complement(lst_kmer,kmer_size);
  } 
  binary_kmer_to_seq(lst_kmer,kmer_size,lst_seq);

  fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s\n", name,
	  (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage,
	  element_get_coverage(fst_node),
	  fst_seq,
	  (fst_orientation == forward ? fst_r : fst_f),
	  (fst_orientation == forward ? fst_f : fst_r),
	  element_get_coverage(lst_node),
	  lst_seq,
	  (lst_orientation == forward ? lst_r : lst_f),
	  (lst_orientation == forward ? lst_f : lst_r));


  if (include_first_kmer){    
    
    fprintf(fout,"%s",binary_kmer_to_seq(fst_kmer,kmer_size,fst_seq));
  }
  
  fprintf(fout,"%s\n",string);
  
}
