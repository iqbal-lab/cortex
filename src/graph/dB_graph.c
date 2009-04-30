
/*
  dB_graph.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <string.h>

//it doesn't check that it is a valid arrow
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


int db_graph_db_node_clip_tip(dBNode * node, int limit,
			     boolean (*condition)(dBNode * node),  void (*node_action)(dBNode * node),
			      dBGraph * db_graph){

  int length_tip = 0;
  char seq[db_graph->kmer_size+1];

  if (DEBUG){
    printf("START Node %s\n",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
  }

  if (condition(node)){

    length_tip = db_graph_db_node_clip_tip_with_orientation(node,forward,limit,node_action,db_graph);
   
    
    if (length_tip!=0){
      if (DEBUG){
	 printf("\tforward tip: %i\n",length_tip);
       }
    }
    else{
      length_tip = db_graph_db_node_clip_tip_with_orientation(node,reverse,limit,node_action,db_graph);
      
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


//remove any node under x coverage (and the arrows)

boolean db_graph_db_node_prune(dBNode * node, int coverage,
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

// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)


int db_graph_get_perfect_path(dBNode * node, Orientation orientation, int limit, void (*node_action)(dBNode * node),
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
  sum_coverage          += element_get_coverage(node);
  *min_coverage          = sum_coverage;
  *max_coverage          = sum_coverage;

  if (DEBUG){
    printf("\nNode %i in path: %s\n", length, binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
   }

  

  if (db_node_has_precisely_one_edge(node, orientation,&nucleotide)){ //first node special case
    

    do{ 
      if (length>0){
	node_action(current_node);
      }

      next_node =  db_graph_get_next_node(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph);

      
      //sanity check
      if(next_node == NULL){
	fprintf(stderr,"dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
	exit(1);
      }

           
      path_labels[length]        = nucleotide;
      seq[length]                = binary_nucleotide_to_char(nucleotide);
      int coverage               = element_get_coverage(next_node);
      sum_coverage               += coverage;
      *max_coverage = *max_coverage < coverage ? coverage : *max_coverage;
      *min_coverage = *min_coverage > coverage ? coverage : *min_coverage;

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
	   db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2) && //multiple entries
	   db_node_has_precisely_one_edge(current_node, current_orientation,&nucleotide)); //has one next edge only


    if ((next_node == node) && (next_orientation == orientation)){
      *is_cycle = true;
    }

  }
  
  if (DEBUG){
    printf("\nLast node in path: %s %i length: %i\n", binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq),next_node->edges,length);
  }

  seq[length] = '\0';
  *avg_coverage = (double) sum_coverage/(double) (length+1);
  return length;
  
}

//a bubble starts in a node with only two outgoing edges in the same orientation
// every branch of the bubble is free of conflicts until it joins with the main branch
// the differente between the length of the brances is <= delta

boolean db_graph_detect_bubble(dBNode * node,
			       Orientation orientation,
			       int limit, int delta,
			       void (*node_action)(dBNode * node), 
			       int * length1, Nucleotide * base1, dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,char * seq1,
			       int * length2, Nucleotide * base2, dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,char * seq2,
			       dBGraph * db_graph){

  

  
  dBNode * next_node1;
  dBNode * next_node2;
  Orientation next_orientation1, next_orientation2;
  Nucleotide rev_base1, rev_base2;
  boolean is_cycle1, is_cycle2;
  boolean ret = false;

  
  if (db_node_has_precisely_two_edges(node,orientation,base1,base2)){
  
    next_node1 = db_graph_get_next_node(node,orientation,&next_orientation1,*base1,&rev_base1,db_graph);
    next_node2 = db_graph_get_next_node(node,orientation,&next_orientation2,*base2,&rev_base2,db_graph);
    
  
    if (next_node1 == NULL || next_node2 == NULL){
      puts("error!");
      exit(1);
    }
    
    double avg_coverage;
    int min_coverage;
    int max_coverage;
    
    *length1  = db_graph_get_perfect_path(next_node1,next_orientation1,
					  limit,node_action,
					  path_nodes1,path_orientations1,path_labels1,
					  seq1,&avg_coverage,&min_coverage,&max_coverage,
					  &is_cycle1,db_graph);

   
    *length2  = db_graph_get_perfect_path(next_node2,next_orientation2,
					  limit,node_action,
					  path_nodes2,path_orientations2,path_labels2,
					  seq2,&avg_coverage,&min_coverage,&max_coverage,
					  &is_cycle2,db_graph);
  
    //action to the begining and end of branches
    node_action(next_node1);
    node_action(next_node2);
    node_action(path_nodes1[*length1]);
    node_action(path_nodes2[*length2]);

    ret = (abs(*length1-*length2)<=delta && 
	   path_nodes1[*length1] == path_nodes2[*length2]);
  }
  return ret;
}



boolean db_graph_detect_perfect_bubble(dBNode * node,
				       Orientation * orientation, 
				       void (*node_action)(dBNode * node), 
				       Nucleotide * base1, Nucleotide * base2, 
				       Nucleotide * labels,dBNode * * end_node,Orientation * end_orientation,
				       dBGraph * db_graph){
  
  boolean ret = false;
  int length1, length2;
  dBNode * path_nodes1[db_graph->kmer_size+1];
  dBNode * path_nodes2[db_graph->kmer_size+1];
  Orientation path_orientations1[db_graph->kmer_size+1];
  Orientation path_orientations2[db_graph->kmer_size+1];
  Nucleotide path_labels1[db_graph->kmer_size];
  Nucleotide path_labels2[db_graph->kmer_size];
  char seq1[db_graph->kmer_size+1];
  char seq2[db_graph->kmer_size+1];

  node_action(node);

 
  if (db_graph_detect_bubble(node,forward,db_graph->kmer_size,0,node_action,
			     &length1,base1,path_nodes1,path_orientations1,path_labels1,seq1,
			     &length2,base2,path_nodes2,path_orientations2,path_labels2,seq2,
			     db_graph)){
    *orientation = forward;
    ret = true;
  }
  else{
    if (db_graph_detect_bubble(node,reverse,db_graph->kmer_size,0,node_action,
			       &length1,base1,path_nodes1,path_orientations1,path_labels1,seq1,
			       &length2,base2,path_nodes2,path_orientations2,path_labels2,seq2,
			       db_graph)){

      *orientation = reverse;
      ret = true;
      
    }
      
  }

  if (ret){
    int i;
    for(i=0;i<length1;i++){
      labels[i] = path_labels1[i];
    }
    
    //sanity check
    if (path_nodes1[length1] != path_nodes2[length2]){
      puts("db_graph_detect_perfect_bubble problem\n");
      exit(1);
    }

    *end_node        = path_nodes1[length1];
    *end_orientation = path_orientations1[length1];
    

  }
  return ret;
}

//clip the branch with smaller coverage in the first kmer --> this can be more sophisticated
//dont want to flatten repeats -- use coverage_limit for this
boolean db_graph_db_node_smooth_bubble(dBNode * node, Orientation orientation, int limit,int delta,int coverage_limit,
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
  char seq2[limit+1];

  int length_tmp;
  Nucleotide reverse_nucleotide;
  Orientation next_orientation;
  Nucleotide base1,base2;

  if (db_graph_detect_bubble(node,orientation,limit,delta,&db_node_action_do_nothing,
			     &length1,&base1,path_nodes1,path_orientations1,path_labels1,seq1,
			     &length2,&base2,path_nodes2,path_orientations2,path_labels2,seq2,
			     db_graph)){
 
    if (coverage_limit>=element_get_coverage(path_nodes1[0]) && 
	coverage_limit>=element_get_coverage(path_nodes2[0])){
          
      ret = true;
      //prune
      int i;

      if (element_get_coverage(path_nodes1[0])<element_get_coverage(path_nodes2[0])){
	path_nodes_tmp = path_nodes1;
	path_orientations_tmp = path_orientations1;
	path_labels_tmp  = path_labels1;
	length_tmp     = length1;
	db_node_reset_edge(node,orientation,base1);
      }
      else{
	path_nodes_tmp = path_nodes2;
	length_tmp     = length2;
	path_orientations_tmp = path_orientations2;
	path_labels_tmp       = path_labels2;
	db_node_reset_edge(node,orientation,base2);
      }
	
      for(i=0;i<length_tmp;i++){
	 node_action(path_nodes_tmp[i]);
	 db_node_reset_edges(path_nodes_tmp[i]);
      }
	
      db_graph_get_next_node(path_nodes_tmp[length_tmp-1],path_orientations_tmp[length_tmp-1],&next_orientation,path_labels_tmp[length_tmp-1],&reverse_nucleotide,db_graph);
      db_node_reset_edge(path_nodes_tmp[length_tmp],opposite_orientation(next_orientation),reverse_nucleotide);

      
    }
    
   
  }  

  return ret;
}

//string has to support limit+db_graph->kmer_size+1 (+1 as you need a space for the \0 at the end)
//node_action has to be idempotent as it can be applied to the same node twice!!

int db_graph_supernode(dBNode * node,int limit, boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node),
		       char * string,dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
		       dBGraph * db_graph){

  dBNode * nodes_reverse[limit];
  Orientation orientation_reverse[limit];
  Nucleotide labels_reverse[limit];
  boolean is_cycle;
  int length_reverse;
  int length = 0;
  char seq[limit+1];
  int min,max;
  double avg_coverage;


  if (condition(node)){    
        
    //compute the reverse path until the end of the supernode
    //return is_cycle_reverse == true if the path closes a loop    
 
    length_reverse = db_graph_get_perfect_path(node,reverse,limit,&db_node_action_do_nothing,
					       nodes_reverse,orientation_reverse,labels_reverse,
					       seq,&avg_coverage,&min,&max,
					       &is_cycle,db_graph);

    //apply action
    node_action(nodes_reverse[length_reverse]);

    //we are at the end of a supernode
    length = db_graph_get_perfect_path(nodes_reverse[length_reverse],opposite_orientation(orientation_reverse[length_reverse]),limit,node_action,
				       path_nodes,path_orientations,path_labels,
				       seq,&avg_coverage,&min,&max,
				       &is_cycle,db_graph);
    
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
  
  return length;
}



//routine to get SNPS 

void db_graph_detect_snps(dBGraph * db_graph){
  
  /* int count_snps = 0; */
		 
/*   void get_snps(dBNode * node){ */
    
/*     Orientation orientation,end_orientation; */
/*     Nucleotide base1,base2; */
/*     Nucleotide labels[db_graph->kmer_size]; */
/*     dBNode * end_node; */

/*     if (db_node_check_status_none(node)){ */

/*       if (db_graph_detect_perfect_bubble(node,&orientation,&db_node_action_set_status_visited,&base1,&base2,labels,&end_node,&end_orientation,db_graph)){ */
      
/* 	int length_flank5p = 0; */
/* 	int length_flank3p = 0; */
/* 	dBNode * nodes5p[100]; */
/* 	dBNode * nodes3p[100]; */
/* 	Orientation orientations5p[100]; */
/* 	Orientation orientations3p[100]; */
/* 	Nucleotide labels_flank5p[100];  */
/* 	Nucleotide labels_flank3p[100]; */
/* 	boolean is_cycle5p, is_cycle3p; */
/* 	char tmp_seq[db_graph->kmer_size+1]; */
/* 	seq */
/* 	int i; */
	
/* 	printf("SNP: %i - coverage: %d\n",count_snps,element_get_coverage(node)); */
/* 	count_snps++; */
	
/* 	printf("five prime end\n"); */
/* 	length_flank5p = db_graph_get_perfect_path(node,opposite_orientation(orientation),100,db_node_action_set_status_visited, */
/* 						   nodes5p,orientations5p,labels_flank5p, */
/* 						   seq,&avg_coverage,&min,&max, */
/* 						   &is_cycle5p,db_graph);   */
/* 	printf("length 5p flank: %i\n",length_flank5p); */
	
/* 	printf("three prime end\n"); */
/* 	length_flank3p = db_graph_get_perfect_path(end_node,end_orientation,100,db_node_action_set_status_visited, */
/* 						   nodes3p,orientations3p,labels_flank3p, */
/* 						   seq,&avg_coverage,&min,&max, */
/* 						   &is_cycle3p,db_graph);     */
/* 	printf("length 3p flank: %i\n",length_flank3p); */
	
/* 	//print flank5p */
/* 	for(i=length_flank5p-1;i>=0;i--){ */
/* 	  printf("%c",reverse_char_nucleotide(binary_nucleotide_to_char(labels_flank5p[i]))); */
/* 	} */

/* 	//print the initial node */
/* 	printf(" "); */
/* 	if (orientation == forward){ */
/* 	  printf ("%s ",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq)); */
/* 	} */
/* 	else{ */
/* 	  printf ("%s ",binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(node),db_graph->kmer_size),db_graph->kmer_size,tmp_seq)); */
/* 	} */
	
/* 	printf (" [%c %c] ",binary_nucleotide_to_char(base1),binary_nucleotide_to_char(base2)); */
	
/* 	//print bubble */
/* 	for(i=0;i<db_graph->kmer_size;i++){ */
/* 	  printf("%c",binary_nucleotide_to_char(labels[i])); */
/* 	} */
	
/* 	printf(" "); */
/* 	//print flank3p */
/* 	for(i=0;i<length_flank3p;i++){ */
/* 	  printf("%c",binary_nucleotide_to_char(labels_flank3p[i])); */
/* 	} */
	
/* 	printf("\n"); */
/*       } */
/*     }     */
/*   } */
/*   hash_table_traverse(&get_snps,db_graph); */
}

//routine to get SNPS 

void db_graph_smooth_bubbles(int coverage,int limit, int delta, dBGraph * db_graph){
  
  void smooth_bubble(dBNode * node){
    if (db_node_check_status_none(node)){
      db_graph_db_node_smooth_bubble(node,forward,limit,delta,coverage,
				     &db_node_action_set_status_pruned,db_graph);
      db_graph_db_node_smooth_bubble(node,reverse,limit,delta,coverage,
				     &db_node_action_set_status_pruned,db_graph);
      
    }
  }
  
  hash_table_traverse(&smooth_bubble,db_graph); 

}




void db_graph_clip_tips(dBGraph * db_graph){
  
  void clip_tips(dBNode * node){
    
    db_graph_db_node_clip_tip(node, db_graph->kmer_size*2,&db_node_check_status_none,&db_node_action_set_status_pruned,db_graph);
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
  char  tmp_seq[db_graph->kmer_size+1];
  BinaryKmer bin_kmer;
  boolean is_cycle = false;
  double avg_coverage;
  int min,max;
  

  path_nodes        = calloc(max_length,sizeof(dBNode*));
  path_orientations = calloc(max_length,sizeof(Orientation));
  path_labels       = calloc(max_length,sizeof(Nucleotide));
  seq               = calloc(max_length+1,sizeof(char));
  

  void print_supernode(dBNode * node){
    
    int count;

    if (db_node_check_status(node,none)){

      int length = 0;
      
      bin_kmer = element_get_kmer(node);
      
      count = db_node_edges_count(node,forward);
      if (count!=1){
	length = db_graph_get_perfect_path(node,reverse,max_length,&db_node_action_set_status_visited,
					   path_nodes, path_orientations,path_labels,
					   seq,&avg_coverage,&min,&max,
					   &is_cycle,db_graph);
	
	bin_kmer = binary_kmer_reverse_complement(bin_kmer, db_graph->kmer_size);
      }
      else{

	count = db_node_edges_count(node,reverse);
	if (count!=1){
	  length = db_graph_get_perfect_path(node,forward,max_length,&db_node_action_set_status_visited,
					     path_nodes, path_orientations,path_labels,
					     seq,&avg_coverage,&min,&max,
					     &is_cycle,db_graph);
	  
	}

      }

      if (length>0){
	db_node_set_status(path_nodes[length],visited);


	fprintf(fout,">node_%i length:%i min_coverage:%i max_coverage:%i average_coverage:%5.2f fst:%i lstr:%i lstf:%i\n",count_nodes,length+db_graph->kmer_size,min,max,avg_coverage,count,db_node_edges_count(path_nodes[length],path_orientations[length]),db_node_edges_count(path_nodes[length],opposite_orientation(path_orientations[length])));
	count_nodes++;
	fprintf(fout,"%s%s\n",binary_kmer_to_seq(bin_kmer,db_graph->kmer_size,tmp_seq),seq);
      }
    }
  }
  
  hash_table_traverse(&print_supernode,db_graph); 
  
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

void db_graph_prune_low_coverage_nodes(int coverage, dBGraph * db_graph){
  
  void prune_node(dBNode * node){
    db_graph_db_node_prune(node,coverage,
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
