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


void db_graph_remove_path(dBNode * node, Orientation orientation, int length, dBGraph * db_graph){

  char seq[db_graph->kmer_size+1];
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


// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)


int db_graph_get_perfect_path(dBNode * node, Orientation orientation, int limit, void (*node_action)(dBNode * node),
			      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
			      boolean * is_cycle, dBGraph * db_graph){

  Orientation  current_orientation,next_orientation;
  dBNode * current_node;
  dBNode * next_node;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  
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
  
  return length;
  
}

//a perfect bubble starts in a node with only two outgoing edges in the same orientation
// every branch of the bubble is free of conflicts until it joins with the main branch
// the length of every branch is the same -- smaller than limit!

boolean db_graph_detect_perfect_bubble(dBNode * node,
				       Orientation * orientation, 
				       boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node), 
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

  if (condition(node)){

    if (db_node_has_precisely_two_edges(node,forward,base1,base2)){
    
      *orientation = forward;
      next_node1 = db_graph_get_next_node(node,forward,&next_orientation1,*base1,&rev_base1,db_graph);
      next_node2 = db_graph_get_next_node(node,forward,&next_orientation2,*base2,&rev_base2,db_graph);
      
      if (next_node1 == NULL || next_node2 == NULL){
	puts("error!");
	exit(1);
      }

      length1  = db_graph_get_perfect_path(next_node1,next_orientation1,
					   db_graph->kmer_size,node_action,
					   nodes1,orientations1,labels1,&is_cycle1,db_graph);

      length2  = db_graph_get_perfect_path(next_node2,next_orientation2,
					   db_graph->kmer_size,node_action,
					   nodes2,orientations2,labels2,&is_cycle2,db_graph);

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
	
	length1 = db_graph_get_perfect_path(next_node1,next_orientation1,db_graph->kmer_size,node_action,nodes1,orientations1,labels1,&is_cycle1,db_graph);
	length2 = db_graph_get_perfect_path(next_node2,next_orientation2,db_graph->kmer_size,node_action,nodes2,orientations2,labels2,&is_cycle2,db_graph);

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
    printf("db_graph_detect_perfect_bubble: error inconsistent state\n");
    exit(1);
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



  if (condition(node)){

    node_action(node);
        
    //compute the reverse path until the end of the supernode
    //return is_cycle_reverse == true if the path closes a loop    
 
    length_reverse = db_graph_get_perfect_path(node,reverse,limit,node_action,
					       nodes_reverse,orientation_reverse,labels_reverse,
					       &is_cycle,db_graph);

    //apply action
    node_action(nodes_reverse[length_reverse]);

    //we are at the end of a supernode
    length = db_graph_get_perfect_path(nodes_reverse[length_reverse],opposite_orientation(orientation_reverse[length_reverse]),limit,node_action,
				       path_nodes,path_orientations,path_labels,
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
  
  int count_snps = 0;
		 
  void get_snps(dBNode * node){
    
    Orientation orientation,end_orientation;
    Nucleotide base1,base2;
    Nucleotide labels[db_graph->kmer_size];
    dBNode * end_node;

    if (db_graph_detect_perfect_bubble(node,&orientation,&db_node_check_status_none,&db_node_action_set_status_visited,&base1,&base2,labels,&end_node,&end_orientation,db_graph)){
      
      int length_flank5p = 0;
      int length_flank3p = 0;
      dBNode * nodes5p[100];
      dBNode * nodes3p[100];
      Orientation orientations5p[100];
      Orientation orientations3p[100];
      Nucleotide labels_flank5p[100]; 
      Nucleotide labels_flank3p[100];
      boolean is_cycle5p, is_cycle3p;
      char tmp_seq[db_graph->kmer_size+1];
      int i;

      printf("SNP: %i - coverage: %d\n",count_snps,element_get_coverage(node));
      count_snps++;

      printf("five prime end\n");
      length_flank5p = db_graph_get_perfect_path(node,opposite_orientation(orientation),100,db_node_action_set_status_visited,
						 nodes5p,orientations5p,labels_flank5p,
						 &is_cycle5p,db_graph);  
      printf("length 5p flank: %i\n",length_flank5p);

      printf("three prime end\n");
      length_flank3p = db_graph_get_perfect_path(end_node,end_orientation,100,db_node_action_set_status_visited,
						 nodes3p,orientations3p,labels_flank3p,
						 &is_cycle3p,db_graph);    
      printf("length 3p flank: %i\n",length_flank3p);

      //print flank5p
      for(i=length_flank5p-1;i>=0;i--){
	printf("%c",reverse_char_nucleotide(binary_nucleotide_to_char(labels_flank5p[i])));
      }

      //print the initial node
      printf(" ");
      if (orientation == forward){
	printf ("%s ",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq));
      }
      else{
	printf ("%s ",binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(node),db_graph->kmer_size),db_graph->kmer_size,tmp_seq));
      }

      printf (" [%c %c] ",binary_nucleotide_to_char(base1),binary_nucleotide_to_char(base2));
      
      //print bubble
      for(i=0;i<db_graph->kmer_size;i++){
	printf("%c",binary_nucleotide_to_char(labels[i]));
      }

      printf(" ");
      //print flank3p
      for(i=0;i<length_flank3p;i++){
	printf("%c",binary_nucleotide_to_char(labels_flank3p[i]));
      }

      printf("\n");
    }
  }    
  hash_table_traverse(&get_snps,db_graph);
}

void db_graph_clip_tips(dBGraph * db_graph){
  
  void clip_tips(dBNode * node){
    
    db_graph_db_node_clip_tip(node, db_graph->kmer_size*2,&db_node_check_status_none,&db_node_action_set_status_pruned,db_graph);
  }

  hash_table_traverse(&clip_tips,db_graph);
  
}


void db_graph_print_supernodes(dBGraph * db_graph){
  int count_nodes=0;
  void print_supernode(dBNode * node){
    dBNode * nodes_path[5000];
    Orientation orientations_path[5000];
    Nucleotide labels_path[5000];
    char seq[5000+db_graph->kmer_size+1];


    int length = db_graph_supernode(node,5000,&db_node_check_status_none,&db_node_action_set_status_visited,
				    seq,nodes_path,orientations_path, labels_path,db_graph);
    if (length>0){
      printf(">node_%i %i\n",count_nodes,length+db_graph->kmer_size);
      count_nodes++;
      printf("%s\n",seq);
    }
  }
  hash_table_traverse(&print_supernode,db_graph); 
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
