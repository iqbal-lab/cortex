
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


void print_fasta_from_path_for_specific_person_or_pop(FILE *fout,
						      char * name,
						      int length,
						      double avg_coverage,
						      int min_coverage,
						      int max_coverage,
                                                      int modal_coverage,
                                                      double  percent_nodes_with_modal_coverage,
						      double percent_novel, 
						      dBNode * fst_node,
						      Orientation fst_orientation,
						      dBNode * lst_node,
						      Orientation lst_orientation,
						      char * string, //labels of paths
						      int kmer_size,
						      boolean include_first_kmer,
						      EdgeArrayType type,
						      int index
						      );


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

 

  if (next_node != NULL)
    {
      *next_orientation = db_node_get_orientation(kmer,next_node,db_graph->kmer_size);
    }
  
  //need to check the node is in this person's graph
  if (! (db_node_is_this_node_in_this_person_or_populations_graph(next_node, type, index)))
    {
      return NULL;
    }
 
  return next_node;
}


/*
  
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

  
   //first node special case
   if (db_node_has_precisely_one_edge(node, orientation,&nucleotide, type, index)){ 
   
   do{ 
   if (length>0)
   {
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

*/





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


int db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit, 
									 Nucleotide fst_nucleotide,
									 void (*node_action)(dBNode * node),
									 dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
									 char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
									 boolean * is_cycle, dBGraph * db_graph, EdgeArrayType type, int index){

  Orientation  current_orientation,next_orientation;
  dBNode * current_node = NULL;
  dBNode * next_node = NULL;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  int sum_coverage = 0;
  int coverage  = 0;

  //sanity checks
  if (node == NULL)
    {
      printf("db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop: can't pass a null node\n");
      exit(1);
    }
  else if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      //printf("\nThis node is not in the graph of this person - in db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop\n");
      return false;
    }

  current_node        = node;
  current_orientation = orientation;

  *is_cycle = false;
  
  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  
  *max_coverage         = 0;
  *min_coverage         = INT_MAX;
 
  if (DEBUG){
    printf("\n ZAM Node %i in path: %s\n", length, binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq));
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

    //this will return NULL if the next node is not in the person's graph. It does NOT check if the edge is in the graph...
    next_node =  db_graph_get_next_node_for_specific_person_or_pop(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph, type, index);

      


    //sanity check
    if(next_node == NULL)
      {
	fprintf(stderr,"dB_graph: didnt find node in hash table: %s %c %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),binary_nucleotide_to_char(nucleotide), current_orientation == forward ? "forward" : "reverse");
	exit(1);
      }


    path_labels[length]        = nucleotide;
    seq[length]                = binary_nucleotide_to_char(nucleotide);
    coverage                   = db_node_get_coverage(next_node, type, index);
    
    length++;

    path_nodes[length]         = next_node;
    path_orientations[length]  = next_orientation;
    
    if (DEBUG){
      printf("\n Length of path so far is  %i : %s\n", length, binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size,tmp_seq));
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

  //debug
  /*
  if ((next_node == node) && (next_orientation == orientation))
    {
      printf("stopped becaise of loop\n");
    }
  else if (!db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide2, type, index))
    {
      printf("stopped because next node has >1 edge in - multple entries\n");
    }
  else if (!db_node_has_precisely_one_edge(current_node, current_orientation,&nucleotide, type, index))
    {
      printf("stopped because current node has >1 edge\n");
    }
  else if (length>=limit)
    {
      printf("stopped becase limit exceeded: length %d and limit %d\n", length, limit);
    }
  else
    {
      printf("WARNING IMPOSSIBLE ZAM");
    }
  */
  
  
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

int db_graph_get_perfect_path_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit, 
							 void (*node_action)(dBNode * node),
							 dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
							 char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
							 boolean * is_cycle, dBGraph * db_graph, EdgeArrayType type, int index)
{

  int length =0;
  Nucleotide nucleotide;

  //sanity check
  if (node == NULL){
    printf("db_graph_get_perfect_path_for_specific_person_or_pop: can't pass a null node\n");
    exit(1);
  }

  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  
 

  if (db_node_has_precisely_one_edge(node,orientation,&nucleotide, type, index))
    {
    
      length= db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(node,orientation, limit, nucleotide,
										   node_action,
										   path_nodes,path_orientations,path_labels,
										   seq,avg_coverage,min_coverage,max_coverage,
										   is_cycle,db_graph, type, index);
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
// the differente between the length of the branches is <= delta

boolean db_graph_detect_bubble_for_specific_person_or_population(dBNode * node,
								 Orientation orientation,
								 int limit, int delta,
								 void (*node_action)(dBNode * node), 
								 int * length1, Nucleotide * base1, dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,char * seq1,
								 int * length2, Nucleotide * base2, dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,char * seq2,
								 dBGraph * db_graph, EdgeArrayType type, int index){


  if (node==NULL)
    {
      printf("Do not call db_graph_detect_bubble_for_specific_person_or_population with NULL node. Exiting\n");
      exit(1);
    }
  else if (! (db_node_is_this_node_in_this_person_or_populations_graph(node, type, index)))
    {
      //printf("\nThis node is not in the graph of this person\n");
      return false;
    }
  
  dBNode * next_node1;
  dBNode * next_node2;
  Orientation next_orientation1, next_orientation2;
  Nucleotide rev_base1, rev_base2;
  boolean is_cycle1, is_cycle2;
  boolean ret = false;


  if (db_node_has_precisely_two_edges(node,orientation,base1,base2,type,index)){
  
    next_node1 = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation1,*base1,&rev_base1,db_graph,type,index);
    next_node2 = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation2,*base2,&rev_base2,db_graph,type,index);
    
  
    if (next_node1 == NULL || next_node2 == NULL){
      puts("error in db_graph_detect_bubble_for_specific_person_or_population! next_node 1 or 2 is NULL ");
      exit(1);
    }
    
    double avg_coverage;
    int min_coverage;
    int max_coverage;
    
    *length1  = db_graph_get_perfect_path_for_specific_person_or_pop(next_node1,next_orientation1,
									    limit,node_action,
									    path_nodes1,path_orientations1,path_labels1,
									    seq1,&avg_coverage,&min_coverage,&max_coverage,
									    &is_cycle1,db_graph, type, index);

   
    *length2  = db_graph_get_perfect_path_for_specific_person_or_pop(next_node2,next_orientation2,
									    limit,node_action,
									    path_nodes2,path_orientations2,path_labels2,
									    seq2,&avg_coverage,&min_coverage,&max_coverage,
									    &is_cycle2,db_graph, type, index);
  

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


/*


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

*/



// it returns the supernode containing 'node'  
// string has to support limit+1 (+1 as you need a space for the \0 at the end)
// node_action has to be idempotent as it can be applied to the same node twice!!
// supernode_str returns the string made of the labels of the path (doesn't include first kmer). 
// returns the length of supernode 


int db_graph_supernode_for_specific_person_or_pop(dBNode * node,int limit,void (*node_action)(dBNode * node), 
						  dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char * supernode_str, 
						  double * avg_coverage,int * min,int * max, boolean * is_cycle, 
						  dBGraph * db_graph, EdgeArrayType type, int index){

  
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
  
  length_reverse = db_graph_get_perfect_path_for_specific_person_or_pop(node,reverse,limit,&db_node_action_do_nothing,
									nodes_reverse,orientations_reverse,labels_reverse,
									supernode_str,&avg_coverager,&minr,&maxr,
									&is_cycler,db_graph, type, index);
  
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
    
    
    length = db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(nodes_reverse[length_reverse],
										  opposite_orientation(orientations_reverse[length_reverse]),
										  limit,label,
										  node_action,
										  path_nodes,path_orientations,path_labels,
										  supernode_str,avg_coverage,min,max,
										  is_cycle,db_graph, type, index);
    
  }
  else{
    length = db_graph_get_perfect_path_for_specific_person_or_pop(node,forward,
								   limit,
								   node_action,
								   path_nodes,path_orientations,path_labels,
								   supernode_str,avg_coverage,min,max,
								   is_cycle,db_graph, type, index);
  }
  
  
  

  //apply action to the fst and last node
  node_action(path_nodes[0]);
  node_action(path_nodes[length]);
  
  return length;
}






// identical code to db_graph_supernode_for_specific_person_or_pop but this returns the index of the query node (first argument) in th supernode array
// will make a significant performance gain compared with getting the path_nodes array back and then searching it

int db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(dBNode * node,int limit,void (*node_action)(dBNode * node), 
									    dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
									    char * supernode_str, double * avg_coverage,int * min,int * max, boolean * is_cycle, 
									    int* query_node_posn,
									    dBGraph * db_graph, EdgeArrayType type, int index){
  
  
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
  
  length_reverse = db_graph_get_perfect_path_for_specific_person_or_pop(node,reverse,limit,&db_node_action_do_nothing,
									nodes_reverse,orientations_reverse,labels_reverse,
									supernode_str,&avg_coverager,&minr,&maxr,
									&is_cycler,db_graph, type, index);

  *query_node_posn = length_reverse;
  
  
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
    
    
    length = db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(nodes_reverse[length_reverse],
										  opposite_orientation(orientations_reverse[length_reverse]),
										  limit,label,
										  node_action,
										  path_nodes,path_orientations,path_labels,
										  supernode_str,avg_coverage,min,max,
										  is_cycle,db_graph, type, index);
    
  }
  else{
    length = db_graph_get_perfect_path_for_specific_person_or_pop(node,forward,
								   limit,
								   node_action,
								   path_nodes,path_orientations,path_labels,
								   supernode_str,avg_coverage,min,max,
								   is_cycle,db_graph, type, index);
  }
  
  
  

  //apply action to the fst and last node
  node_action(path_nodes[0]);
  node_action(path_nodes[length]);
  
  return length;
}




/*
int db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop(dBNode * node,int limit, boolean (*condition)(dBNode * node), void (*node_action)(dBNode * node),
									    char * string,dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, int* query_node_posn, 
									    dBGraph * db_graph, EdgeArrayType type, int index){

  if (node==NULL)
    {
      printf("do not call db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop with NULL node\n");
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

*/





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



/*


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

*/







//TODO - move thiscode out into /sv_trio
//the final argument returns the number of chromosomes intersected, but only gives the right answer if the function returns true - ie no node intersected >1 chromosome
boolean db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome(dBNode* node, EdgeArrayType type, int index, dBGraph* dbgraph, int* total_number_of_different_chromosomes_intersected)
{

  printf("Stop using db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome");
  exit(1);

}



//assume we have checked and no node has >1 chrom intwersection
void db_graph_print_chrom_intersections_for_supernode_for_specific_person_or_pop(FILE * file, dBNode * node, dBGraph * db_graph, EdgeArrayType type, int index, 
						      boolean is_for_testing, char** for_test, int* index_for_test)
{
  //deprecated
  printf("Stop using db_graph_print_chrom_intersections_for_supernode_for_specific_person_or_pop");
  exit(1);
}




void db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population(dBGraph* hash_table, EdgeArrayType type, int index)
{
  printf("not implemented db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population yet");
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
  int person_with_best_sub_supernode=-1;
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

  Orientation correct_direction_to_go=forward;

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
int db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta
                                                  (FILE* chrom_fptr, 
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
      path_string[i]     = 'X';//so we spot a problem if we start seeing X's.
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



int int_cmp(const void *a, const void *b)
{
  const int *ia = (const int *)a; // casting pointer types
  const int *ib = (const int *)b;
  return *ia  - *ib;
  /* integer comparison: returns negative if b > a
     and positive if a > b */
}


void get_coverage_from_array_of_nodes(dBNode** array, int length, int* min_coverage, int* max_coverage, double* avg_coverage, int* mode_coverage, double*  percent_nodes_having_modal_value, EdgeArrayType type, int index)
{

  
  int sum_coverage=0;

  *max_coverage         = 0;
  *min_coverage         = INT_MAX;

  int i;
  int coverages[length];

  for (i=0; i< length; i++)
    {
      int this_covg = 0;
      
      if (array[i]!= NULL)
	{
	  this_covg = db_node_get_coverage(array[i], type, index);
	}
      
      coverages[i]=this_covg; //will use this later, for the mode

      sum_coverage += this_covg;
      *max_coverage = *max_coverage < this_covg ? this_covg : *max_coverage;
      *min_coverage = *min_coverage > this_covg ? this_covg : *min_coverage;
      //printf("i is %d, this node has coverage %d, min is %d, max is %d\n", i, this_covg, *min_coverage, *max_coverage);

    }  

  if (*min_coverage==INT_MAX)
    {
      *min_coverage=0;
    }
  *avg_coverage = sum_coverage/length;


  qsort( coverages, length, sizeof(int), int_cmp);
  int covg_seen_most_often=coverages[0];
  int number_of_nodes_with_covg_seen_most_often=1;
  int current_run_of_identical_adjacent_covgs=1;

  for (i=1; i< length; i++)
    {
      if (coverages[i]==coverages[i-1])
	{
	  current_run_of_identical_adjacent_covgs++;
	}
      else
	{
	  current_run_of_identical_adjacent_covgs=1;
	}

      if (current_run_of_identical_adjacent_covgs > number_of_nodes_with_covg_seen_most_often)
	{
	  number_of_nodes_with_covg_seen_most_often = current_run_of_identical_adjacent_covgs;
	  covg_seen_most_often=coverages[i];
	}
    }

  *mode_coverage = covg_seen_most_often;
  *percent_nodes_having_modal_value = 100* number_of_nodes_with_covg_seen_most_often/length;

}


void get_percent_novel_from_array_of_nodes(dBNode** array, int length, double* percent_novel, EdgeArrayType type_for_reference, int index_of_reference_in_array_of_edges)
{
  int sum_novel=0;

  int i;
  for (i=0; i< length; i++)
    {
      if (array[i] != NULL)
	{
	  //boolean this_node_is_novel = !(db_node_check_status_exists_in_reference(array[i]) || db_node_check_status_visited_and_exists_in_reference(array[i]) );
	  boolean this_node_is_novel = !db_node_is_this_node_in_this_person_or_populations_graph(array[i], type_for_reference, index_of_reference_in_array_of_edges);
	  if (this_node_is_novel==true)
	    {
	      sum_novel++;
	    }
	}
    }

  *percent_novel =  100*sum_novel/length;
}



// *******************************************************************************
// New SV calling algorithm based on comparing a supernode with a "trusted path". This trusted path might be the reference, or some contig that
// we have obtained by bootstrapping.
// ********************************************************************************
// Walks along chromosome path, comparing with supernodes, and for each supernode that touches the trusted path, sees where it attaches.
// max_anchor_span is the biggest gap we allow/look for between start of 5' and end of 3' anchors
//       If you want to be able to find, say 10kb deletions, then set this to 10000 + sum of your desired minimum anchors/flanks
// Length of arrays should be double the max_anchor_span, and max_expected_size)of_supernode must be < max_anchor_span

// min_three/fiveprime_flank_anchor is counted in number of nodes
// length_of_arrays MUST BE EVEN
// returns number of variants found. 
// if max_desired_returns>0, then the first max_desired_returns results are returned in the preallocated arrays branch1_array and branch2_array
// In normal use, this should be zero - we'll find far too many variants. But can be used for testing.
// the edgearraytype and index for the reference are purely used for checking if nodes exist in the reference at all, or are novel. The trusted path is, in general, not necessarily the reference.
// The trusted path comes entirely from chrom_fptr, and doe not need to be the same as the reference, as specified in arguments 4,5
int db_graph_make_reference_path_based_sv_calls(FILE* chrom_fasta_fptr, EdgeArrayType which_array_holds_indiv, int index_for_indiv_in_edge_array,
						EdgeArrayType which_array_holds_ref, int index_for_ref_in_edge_array,
						int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, int max_anchor_span, int min_covg, int max_covg, 
						int max_expected_size_of_supernode, int length_of_arrays, dBGraph* db_graph, FILE* output_file,
						int max_desired_returns,
						char** return_flank5p_array, char** return_trusted_branch_array, char** return_variant_branch_array, 
						char** return_flank3p_array, int** return_variant_start_coord)
{

  int num_variants_found=0;
  
  //makes life much simpler to insist the array is even length.
  if (length_of_arrays%2 !=0)
    {
      printf("Must only call db_graph_make_reference_path_based_sv_calls  with even length_of_arrays\n");
      exit(1);
    }

  //insist max_anchor_span=length_of_arrays/2
  if (max_anchor_span!=length_of_arrays/2)
    {
      printf("You have max_anchor_span = %d, and length_of_arrays = %d. If calling db_graph_make_reference_path_based_sv_calls, must have max_anchor span as half of length_of_arrays.\n",
	     max_anchor_span, length_of_arrays);
      exit(1);
    }

  //insist max length of supernodes <= max_anchor_span
  if (max_expected_size_of_supernode>max_anchor_span)
    {
      printf("You have called db_graph_make_reference_path_based_sv_calls  with arguments: max_expected_size_of_supernode=%d, and max_anchor_span=%d", max_expected_size_of_supernode, max_anchor_span);
      exit(1);
    }

  int number_of_nodes_to_load=length_of_arrays/2;





  //*************************************
  // malloc and initialising
  //*************************************
  dBNode**     chrom_path_array        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); //everything these dBNode*'s point to will be in the hash table - ie owned by the hash table
  Orientation* chrom_orientation_array = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
  Nucleotide*  chrom_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  char*        chrom_string            = (char*) malloc(sizeof(char)*length_of_arrays+1); //+1 for \0

  dBNode**     current_supernode       = (dBNode**) malloc(sizeof(dBNode*)*(max_expected_size_of_supernode+ db_graph->kmer_size));
  Orientation* curr_sup_orientations   = (Orientation*) malloc(sizeof(Orientation)*(max_expected_size_of_supernode+ db_graph->kmer_size));
  Nucleotide*  curr_sup_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_expected_size_of_supernode+ db_graph->kmer_size));
  char*        supernode_string        = (char*) malloc(sizeof(char)*((max_expected_size_of_supernode+ db_graph->kmer_size)+1)); //+1 for \0



  int n;
  for (n=0; n<length_of_arrays; n++)
    {
      chrom_path_array[n]=NULL;
      chrom_orientation_array[n]=forward;
      chrom_labels[n]=Undefined;
      chrom_string[n]='N';
    }
  for (n=0; n< (max_expected_size_of_supernode+ db_graph->kmer_size); n++)
    {
      current_supernode[n]=NULL;
      curr_sup_orientations[n]=forward;
      curr_sup_labels[n]=Undefined;
      supernode_string[n]='N';
    }
  chrom_string[0]='\0';
  chrom_string[length_of_arrays]='\0';
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
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*length_of_arrays);    //*(number_of_nodes_to_load + db_graph->kmer_size));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  
  kmer_window->nkmers=0;



  //some strings in which to hold the long sequences we will print.
  char* trusted_branch = malloc(sizeof(char)*(length_of_arrays+1));
  char* variant_branch = malloc(sizeof(char)*(max_expected_size_of_supernode+1));
  
  if ( (trusted_branch==NULL) || (variant_branch==NULL) )
    {
      printf("OOM. Unable to malloc trusted and variant branches. Exit. \n");
      exit(1);
    }

  int k;
  for (k=0; k<length_of_arrays; k++)
    {
      trusted_branch[k]=0;
    }
  for (k=0; k<max_expected_size_of_supernode; k++)
    {
      variant_branch[k]=0;
    }
  trusted_branch[0]='\0';
  trusted_branch[length_of_arrays]=0;
  variant_branch[0]='\0';
  variant_branch[max_expected_size_of_supernode]=0;





  // ***********************************************************
  //********** end of malloc and initialise ********************


  //in order to be able to report the coordinates of the SV we find, with respect to the trusted path fasta (usually the fasta of a reference chromosome)
  // we keep track of what coordinate in that fasta is the first base in the array.
  int coord_of_start_of_array_in_trusted_fasta=0;


  //load a set of nodes into thr back end (right hand/greatest indices) of array
  // each call   will push left the nodes/etc in the various arrays by length_of_arrays/2=max_anchor_span
  // and then put the new nodes etc in on the right of that

  int ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fasta_fptr, number_of_nodes_to_load, 0, 
													     length_of_arrays,
													     chrom_path_array, chrom_orientation_array, 
													     chrom_labels, chrom_string,
													     seq, kmer_window, 
													     true, false,
													     db_graph);



  if (ret ==number_of_nodes_to_load)
    {
      //one more batch, then array is full, and ready for the main loop:
      ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fasta_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													     length_of_arrays,
													     chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													     seq, kmer_window, 
													     false, false,
													     db_graph);
    }



  coord_of_start_of_array_in_trusted_fasta=0;


  char tmp_zam[db_graph->kmer_size];




  // ********************************************************
  /* the core loop of this function: basic structure is
     **********************************************************

  do
    {
    while (start_node_index<max_anchor_span)
       {
  
           look at the supernode in individual which contains chrom_path_array[start_node_index], and see if we can find 5' and 3' anchors
           that match thechromsome/reference.
	   If yes, print it. If no, increment start_node_index.
       }
    }
    while (load another max_anchor_span bases into the rear end of our array)



    The idea is that we will find any variant where there is a change/variation mid-supernode, PROVIDED the supernode
    is no longer than max_anchor_span

                                                                         
  */
                                                                   




  do
    {
      int start_node_index=0;//this is a position in chrom_path_array 

      while (start_node_index<max_anchor_span)
	{

	  //NULL node means you have an N, so move on, or you're at the end of the array, which is full of NULLs
	  if (chrom_path_array[start_node_index]==NULL)
	    {
	      start_node_index++;
	      continue;
	    }

	  //if this chromosome node does not exist in the person's graph, move on to the next node
	  if (db_node_is_this_node_in_this_person_or_populations_graph(chrom_path_array[start_node_index], which_array_holds_indiv, index_for_indiv_in_edge_array) == false)
	    {
	      start_node_index++;
	      continue;
	    }

	  if (db_node_check_status_visited(chrom_path_array[start_node_index])==true)
	  {
	    start_node_index++;
	    continue;
	  }

	  int index_of_query_node_in_supernode_array=-1;

	  
	  //coverage variables:
	  double avg_coverage=0;
	  int min_coverage=0;
	  int max_coverage=0;
	  boolean is_cycle=false;
	  
	  
	  int length_curr_supernode = db_graph_supernode_returning_query_node_posn_for_specific_person_or_pop
	                                      (chrom_path_array[start_node_index], max_expected_size_of_supernode, &db_node_action_set_status_visited,
					       current_supernode, curr_sup_orientations, curr_sup_labels, supernode_string, 
					       &avg_coverage, &min_coverage, &max_coverage, &is_cycle, 
					       &index_of_query_node_in_supernode_array, 
					       db_graph, which_array_holds_indiv, index_for_indiv_in_edge_array);

	  if (length_curr_supernode>max_anchor_span)
	    {
	      printf("Warning - bad choice of params. Current supernode is length %d, but max anchro size is %d\n", length_curr_supernode, max_anchor_span);
	      exit(1);
	    }

	  //char tmp_seqzam[db_graph->kmer_size];
	  //printf("Start looking at the supernode in indiv nucleated at query node %s\n Query node position is %d, Supernode is %s\n", 
	  //	 binary_kmer_to_seq(chrom_path_array[start_node_index]->kmer, db_graph->kmer_size, tmp_seqzam), index_of_query_node_in_supernode_array, supernode_string);


	  if (index_of_query_node_in_supernode_array==-1)
	    {
	      printf("Warning - red alert!!!  - failed to get index of query\n");
	      start_node_index++;
	      exit(1);
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



	  // Exclude/ignore hubs and singletons as well as supernodes without enough length to generate our anchors
	  // check if there is room between query node and either end of supernode to fit both anchors

	  if (traverse_sup_left_to_right)
	    {
	      if (length_curr_supernode - index_of_query_node_in_supernode_array < min_fiveprime_flank_anchor+min_threeprime_flank_anchor)
		{
		  //printf("Insuff room on supernode for anchors at start_node_index %d, corresponding to kmer %s . length of current supernode is %d, index of query node in supernode is %d,Move to next position\n", 
		  // start_node_index, binary_kmer_to_seq(chrom_path_array[start_node_index]->kmer, db_graph->kmer_size, tmp_zam),length_curr_supernode, index_of_query_node_in_supernode_array);
		  start_node_index++;
		  continue;
		}
	    }
	  else
	    {
	      if (index_of_query_node_in_supernode_array < min_fiveprime_flank_anchor+min_threeprime_flank_anchor)
		{
		  //printf("Insuff room on supernode for anchors at start_node_index %d, corresponding to kmer %s . length of current supernode is %d, index of query node in supernode is %d, traversing supernode from right to leftMove to next position\n", 
		  //	 start_node_index, binary_kmer_to_seq(chrom_path_array[start_node_index]->kmer, db_graph->kmer_size, tmp_zam),length_curr_supernode, index_of_query_node_in_supernode_array);
		  start_node_index++;
		  continue;
		}

	    }


	  //now see how far along the chromosome you go before the supernode breaks off and differs
	  int first_index_in_chrom_where_supernode_differs_from_chromosome = start_node_index;
	  int index_in_supernode_where_supernode_differs_from_chromosome=index_of_query_node_in_supernode_array;


	  // complete paranoia - might delete these checks
	  if (chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome]!=current_supernode[index_in_supernode_where_supernode_differs_from_chromosome])
	    {
	      printf("WARNING 1. two arrays do not meet at start.");
	      exit(1);
	    }
	  else if (index_in_supernode_where_supernode_differs_from_chromosome > length_curr_supernode )  
	    {
	      printf("WARNING 2.index is %d and length of sup is %d", index_in_supernode_where_supernode_differs_from_chromosome, length_curr_supernode);
	      exit(1);
	    }
	  else if (index_in_supernode_where_supernode_differs_from_chromosome<0)
	    {
	      printf("WARNING 3. Index is %d\n", index_in_supernode_where_supernode_differs_from_chromosome);
	      exit(1);
	    }
	  else if (first_index_in_chrom_where_supernode_differs_from_chromosome-start_node_index >= max_anchor_span)
	    {
	      printf("WARNING 4\n");
	      exit(1);
	    }


	  //find how far along can go on reference before supernode branches off and starts to differ from the ref/chromosome
	  if (traverse_sup_left_to_right==true)
	    {
	      while ( (index_in_supernode_where_supernode_differs_from_chromosome < length_curr_supernode )  
		  && (first_index_in_chrom_where_supernode_differs_from_chromosome-start_node_index<max_anchor_span) 
		  && (chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome]
		      ==current_supernode[index_in_supernode_where_supernode_differs_from_chromosome])		  
		  )
		{
		  index_in_supernode_where_supernode_differs_from_chromosome++;
		  first_index_in_chrom_where_supernode_differs_from_chromosome++;
		}
	      
	    }
	  else
	    {
	      while ( 
		     (index_in_supernode_where_supernode_differs_from_chromosome>=0)
		     && (first_index_in_chrom_where_supernode_differs_from_chromosome-start_node_index<max_anchor_span) 
		     && (chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome]
			 ==current_supernode[index_in_supernode_where_supernode_differs_from_chromosome])		  
		     
		     )
		{
		  index_in_supernode_where_supernode_differs_from_chromosome--;
		  first_index_in_chrom_where_supernode_differs_from_chromosome++;
		}

	    }


	  if (index_in_supernode_where_supernode_differs_from_chromosome==index_of_query_node_in_supernode_array)
	    {
	      printf("WARNING Did not even go into loop working way along ref - so supernode only touches ref at initial node - I don't think this should be possible, should go into that loop once\n");
	      start_node_index++;
	      exit(1);
	      continue;
	    }

	  if ( ( (traverse_sup_left_to_right) && (index_in_supernode_where_supernode_differs_from_chromosome>=length_curr_supernode) )
	    ||
	       ( (!traverse_sup_left_to_right) && (index_in_supernode_where_supernode_differs_from_chromosome<=0) )
	       )
	    {
	      //then actually the whole supernode matches the reference exactly, and first_index_in_chrom_where_supernode_differs_from_chromosome is actually just one base beyond the length of the supernode
	      start_node_index = first_index_in_chrom_where_supernode_differs_from_chromosome; 
	      continue;
	    }

	  if ( chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome]==NULL)
	    {
	      //then the only reason they differ is an N in the trusted path - will create a false variant. Ignore and move on.
	      start_node_index = first_index_in_chrom_where_supernode_differs_from_chromosome+1;
	      continue;
	    }

	  if (  abs(index_in_supernode_where_supernode_differs_from_chromosome-index_of_query_node_in_supernode_array) <= min_fiveprime_flank_anchor-1)
	    {
	      //does not have sufficient 5' anchor
	      start_node_index=first_index_in_chrom_where_supernode_differs_from_chromosome;
	      // printf("supernode differs from ref/chromosome too soon. Insufficient space for 5' anchor. index_in_supernode_where_supernode_differs_from_chromosome is %d, and index of query node in supernode array is %d - move on, and set start_node_index to %d\n", 
	      //   index_in_supernode_where_supernode_differs_from_chromosome, index_of_query_node_in_supernode_array, start_node_index);
	      continue;
	    }
	      
	  
	  //We now have  decent 5-prime anchor, and we know we don't match the ref exactly. Note that the point at which they differ may be an N in the trusted path, so may be a NULL ptr
	  // Now see if the END of the supernode matches anywhere in our current array of reference chromosome
	  
	  int start_of_3prime_anchor_in_sup;
	  
	  if (traverse_sup_left_to_right==true)
	    {
	      start_of_3prime_anchor_in_sup=length_curr_supernode-min_threeprime_flank_anchor; 
	      //printf("trav sup left to right, length supernodes is %d, so start of 3' anchor is %d\n", length_curr_supernode, start_of_3prime_anchor_in_sup);

	      //this needs to be *after* the position where supernode first differs from trusted path
	      if (start_of_3prime_anchor_in_sup<= index_in_supernode_where_supernode_differs_from_chromosome)
		{
		  start_node_index=first_index_in_chrom_where_supernode_differs_from_chromosome;
		  continue;
		}

	    }
	  else
	    {
	      start_of_3prime_anchor_in_sup = 0 + min_threeprime_flank_anchor;
	      //printf("trav sup right to left, length supernodes is %d, and start of 3' anchor is %d\n", length_curr_supernode, start_of_3prime_anchor_in_sup);

	      //this needs to be *after* the position where supernode first differs from trusted path, in the direction in which we are walking the supernode (right to left)
	      if (start_of_3prime_anchor_in_sup>= index_in_supernode_where_supernode_differs_from_chromosome)
		{
		  start_node_index=first_index_in_chrom_where_supernode_differs_from_chromosome;
		  continue;
		}

	      

	    }



	  //walk along entire chromosome array, and see if can attach the 3prime anchor
	  // it may attach multiple times, but we will only find the closest
	  
	  int start_of_3prime_anchor_in_chrom = first_index_in_chrom_where_supernode_differs_from_chromosome;
	  boolean found_other_anchor=false;
	  
	  while ( (found_other_anchor==false) && ( start_of_3prime_anchor_in_chrom < length_of_arrays-min_threeprime_flank_anchor))
	    {
	      start_of_3prime_anchor_in_chrom++;
	      //printf("Start of 3' anchor in chrom is %d\n", start_of_3prime_anchor_in_chrom);

	      boolean potential_anchor=true;
	      int j;
	      for (j=0 ; ((j< min_threeprime_flank_anchor) && (potential_anchor==true)) ; j++)
		{
		  int k=j;
		  if (traverse_sup_left_to_right==false)
		    {
		      k=-j;
		    }
		  
		  
		  //cannot have a null ptr in supernode, so this removes cases of N in chrom array
		  if (chrom_path_array[start_of_3prime_anchor_in_chrom + j] != current_supernode[start_of_3prime_anchor_in_sup + k])
		    {
		      //printf("chrom node %d does not work as start of anchor\n", start_of_3prime_anchor_in_chrom );
		      potential_anchor=false;
		    }
		  else if (chrom_path_array[start_of_3prime_anchor_in_chrom + j]==NULL)//paranoia - impossible
		    {
		      //printf("chrom node %d and sup node %d are both NULL. So this cann't be an anchor\n", start_of_3prime_anchor_in_chrom + j, start_of_3prime_anchor_in_sup + k);
		      potential_anchor=false;
		    }


		  //debug
		  //if ( (chrom_path_array[start_of_3prime_anchor_in_chrom + j] != NULL)&& (current_supernode[start_of_3prime_anchor_in_sup + k]!=NULL) )
		  // {
		  //   char tmp_dbg1[db_graph->kmer_size];
		  //   char tmp_dbg2[db_graph->kmer_size];
		      //printf("Compare chrom node %d, %s,  and supernode node %d, %s,\n", start_of_3prime_anchor_in_chrom + j, 
		      //	     binary_kmer_to_seq(chrom_path_array[start_of_3prime_anchor_in_chrom + j]->kmer, db_graph->kmer_size, tmp_dbg1),
		      //	     start_of_3prime_anchor_in_sup + k,
		      //	     binary_kmer_to_seq(current_supernode[start_of_3prime_anchor_in_sup + k]->kmer, db_graph->kmer_size, tmp_dbg2));
		  // }


		}
	      
	      if (potential_anchor==true)
		{
		  found_other_anchor=true;
		}

	    }
	      
	      
	  if (found_other_anchor==false)
	    {
	      //printf("Did not find 3' anchor\n");
	      start_node_index++;
	      continue;
	    }
	  else
	    {

	      //we have found a potential SV locus.
	      num_variants_found++;


	      // Note if max_desired_returns>0, then we are going to return the first max_desired_returns  results in the arguments passed in for this purpose, return_branch1_array etc
	      // We can check this easily - if num_variants_found<max_desired_returns, then we add this result to the arrays. This is essentically used only for testing.


	      
	      // this is length of 5p flank in number of edges in the 5p flank.
	      // (in addition there will be one kmer of bases of course. The print function we call will handle this for us)
	      int length_5p_flank = first_index_in_chrom_where_supernode_differs_from_chromosome-start_node_index-1;


	      //work out start coordinate and print it out
	      int start_coord_of_variant_in_trusted_path_fasta = coord_of_start_of_array_in_trusted_fasta + first_index_in_chrom_where_supernode_differs_from_chromosome + db_graph->kmer_size;
	      fprintf(output_file, "VARIATION: %d, start coordinate %d\n", num_variants_found, start_coord_of_variant_in_trusted_path_fasta);
	      
	      if (num_variants_found<=max_desired_returns) //if we are asking for some of the things we find to be returned within the arguments passed in. See comments at start of this whole function
		{
		  *(return_variant_start_coord[num_variants_found-1]) = start_coord_of_variant_in_trusted_path_fasta;
		}



	      int length_3p_flank = min_threeprime_flank_anchor;//we can do better than this. Probably, the two branches will be identical for some long stretch. Also quite possibly the 3prime anchor
	                                                        //will extend further in the 3prime direction. We avoided doing this in the main search loop for efficiency's sake 
	      
	      
	      
	      //we are going to see if we can extend the 3prime anchor in the 3prime and 5prime directions. Which variables are we going to update, and what do we
	      // have to pay attention to if we do that? 
	      // will update primarily start_of_3prime_anchor_in_chrom, and start_of_3prime_anchor_in_sup.
	      // consequently length_3p_flank will change. I think that's it! No problem huh?

	      int how_many_steps_in_3prime_dir_can_we_extend=0;
	      int how_many_steps_in_5prime_dir_can_we_extend=0;

	      //check in 5prime dir first. 
	      int j=1;
	      boolean can_extend_further_in_5prime_dir=true;
	      
	      if (
		  ( (traverse_sup_left_to_right)&&(start_of_3prime_anchor_in_sup-1==index_in_supernode_where_supernode_differs_from_chromosome) )
		  ||
		  ( (!traverse_sup_left_to_right)&&(start_of_3prime_anchor_in_sup+1==index_in_supernode_where_supernode_differs_from_chromosome) )
		  )
		{
		  can_extend_further_in_5prime_dir=false;
		}


	      
	      //Here's an example of why you have to be careful
	      // Supernode:      KingZamIqbalIsTheBossOfEngland
	      // Ref:            KingZamFrenchFrenchFrenchFrenchZamIqbalIsTheBossOfEngland

	      //Our algorithm would say: OK, KingZam is the 5prime anchor, let's see if we can find a 3prime anchor that is 7 characters long (say). Aha! England will do. OK - can we extend backwards? Yes, we can extend backwards 
	      // all the way back to Zam. But WAIT!! Zam is before the split!! That's fine, it means there is a repeat, just don't try and extend beyond there. Hence the else if condition marked below with ***

	      //(By the way this is not a pathologival case - happens in Human chromosome 1)

	      while (can_extend_further_in_5prime_dir==true)  
		{
		  int k=j;

		  if (traverse_sup_left_to_right==false)
		    {
		      k=-j;
		    }
		  
		  
		  //debug
		  //char tmp_dbg1[db_graph->kmer_size];
		  //char tmp_dbg2[db_graph->kmer_size];
		  //if (chrom_path_array[start_of_3prime_anchor_in_chrom - j]==NULL)
		  //  {
		      //printf("Cannot extend further - have hit an N in the trsuted path\n");
		  //  }
		  //else
		  // {
		      //printf("Extending 3prime anchor in 5prime dir. Compare chrom node %d, %s,  and supernode node %d, %s,\n", 
		      //	     start_of_3prime_anchor_in_chrom - j,
		      //	     binary_kmer_to_seq(chrom_path_array[start_of_3prime_anchor_in_chrom - j]->kmer, db_graph->kmer_size, tmp_dbg1),
		      //	     start_of_3prime_anchor_in_sup - k,
		      //	     binary_kmer_to_seq(current_supernode[start_of_3prime_anchor_in_sup - k]->kmer, db_graph->kmer_size, tmp_dbg2));
		  //  }
		  //end debug

		  if (chrom_path_array[start_of_3prime_anchor_in_chrom - j] != current_supernode[start_of_3prime_anchor_in_sup - k])
		    {
		      //this condition will cover cases when there is an N/Null in the trusted path/chrom_path_array
		      can_extend_further_in_5prime_dir=false;
		    }
		  else if (  (traverse_sup_left_to_right) && (start_of_3prime_anchor_in_sup - j <= index_in_supernode_where_supernode_differs_from_chromosome+1) ) // *** see comment above this while loop
		    {
		      can_extend_further_in_5prime_dir=false;
		    }
		  else if ( (!traverse_sup_left_to_right) && (start_of_3prime_anchor_in_sup + j>= index_in_supernode_where_supernode_differs_from_chromosome-1) )
		    {
		      can_extend_further_in_5prime_dir=false;
		    }
		  else if (start_of_3prime_anchor_in_chrom-how_many_steps_in_5prime_dir_can_we_extend<=first_index_in_chrom_where_supernode_differs_from_chromosome+1)
		    {
		      can_extend_further_in_5prime_dir=false;
		    }
		  else
		    {
		      how_many_steps_in_5prime_dir_can_we_extend++;
		      j++;
		    }

		  //make sure we stop before the end of the supernode. traverse_sup_left_to_right==true => 5prime direction is in direction of decreasing index
		  //this is surely unnecessary - leaving in out of paranoia
		  if (    ((traverse_sup_left_to_right) &&  (start_of_3prime_anchor_in_sup - k -1==0) )
			  ||
			  ((!traverse_sup_left_to_right) && (start_of_3prime_anchor_in_sup - k +1==length_curr_supernode-1) ) //zam changed this from -2 to -1
			  )
		  {
		    can_extend_further_in_5prime_dir=false;
		  }
		  

		  //sanity
		  if (start_of_3prime_anchor_in_chrom-how_many_steps_in_5prime_dir_can_we_extend<=first_index_in_chrom_where_supernode_differs_from_chromosome)
		    {
		      //something has gone wrong
		      printf("We are extending the 3prime anchor in 5prime direction and have gone back beyond the point at which the variant starts in the TRUSTED array!!\n");
		      printf("start_of_3prime_anchor_in_chrom = %d\nhow_many_steps_in_5prime_dir_can_we_extend=%d\n, first_index_in_chrom_where_supernode_differs_from_chromosome=%d\n",
		  	     start_of_3prime_anchor_in_chrom, how_many_steps_in_5prime_dir_can_we_extend,  first_index_in_chrom_where_supernode_differs_from_chromosome);
		      exit(1);
		    }
		  if (traverse_sup_left_to_right)
		    {
		      if (start_of_3prime_anchor_in_sup - how_many_steps_in_5prime_dir_can_we_extend <= index_in_supernode_where_supernode_differs_from_chromosome)
			{
			  printf("We are extending the 3prime anchor in 5prime direction and have gone back beyond the point in the supernode where it separated from the 5prime flank. Should not reach here.\n");
			  printf("start_of_3prime_anchor_in_sup = %d\nhow_many_steps_in_5prime_dir_can_we_extend=%d\n, index_in_supernode_where_supernode_differs_from_chromosome=%d\n",
				 start_of_3prime_anchor_in_sup, how_many_steps_in_5prime_dir_can_we_extend, index_in_supernode_where_supernode_differs_from_chromosome);
			  printf("We are traversing supernode from left to right. Exit.\n");
			  exit(1);
			}
		    }
		  else
		    {
		      if (start_of_3prime_anchor_in_sup + how_many_steps_in_5prime_dir_can_we_extend >= index_in_supernode_where_supernode_differs_from_chromosome)
			{
			  printf("We are extending the 3prime anchor in 5prime direction and have gone back beyond the point in the supernode where it separated from the 5prime flank. Should not reach here.\n");
			  printf("start_of_3prime_anchor_in_sup = %d\nhow_many_steps_in_5prime_dir_can_we_extend=%d\n, index_in_supernode_where_supernode_differs_from_chromosome=%d\n",
				 start_of_3prime_anchor_in_sup, how_many_steps_in_5prime_dir_can_we_extend, index_in_supernode_where_supernode_differs_from_chromosome);
			  printf("We are traversing supernode from right to left. Exit.\n");
			  exit(1);
			}
		      
		    }

		}		  


	      //check in 3prime dir 
	      j=length_3p_flank+1;
	      boolean can_extend_further_in_3prime_dir=false;

	      if (  ( ((traverse_sup_left_to_right) &&  (start_of_3prime_anchor_in_sup+length_3p_flank<length_curr_supernode) )
		    ||
		     ((!traverse_sup_left_to_right) &&  (start_of_3prime_anchor_in_sup-length_3p_flank>0) ))
		   &&
		   (start_of_3prime_anchor_in_chrom+length_3p_flank<length_of_arrays)
		   )
		{
		  can_extend_further_in_3prime_dir=true;
		}


	      while (can_extend_further_in_3prime_dir==true) 
		{
		  int k=j;

		  if (traverse_sup_left_to_right==false)
		    {
		      k=-j;
		    }
		  

		  
		  //char tmp_dbg1[db_graph->kmer_size];
		  //char tmp_dbg2[db_graph->kmer_size];
		  //printf("Extending 3prime anchor in 3prime dir. Compare chrom node %d, %s,  and supernode node %d, %s,\n", 
		  //	 start_of_3prime_anchor_in_chrom + j,
		  //	 binary_kmer_to_seq(chrom_path_array[start_of_3prime_anchor_in_chrom + j]->kmer, db_graph->kmer_size, tmp_dbg1),
		  //	 start_of_3prime_anchor_in_sup + k,
		  //	 binary_kmer_to_seq(current_supernode[start_of_3prime_anchor_in_sup + k]->kmer, db_graph->kmer_size, tmp_dbg2));
		  

		    //we don't need to worry about N's in the chromosome array - there will never be any in the supernode array
		  if (chrom_path_array[start_of_3prime_anchor_in_chrom + j] != current_supernode[start_of_3prime_anchor_in_sup + k])
		    {
		      can_extend_further_in_3prime_dir=false;
		    }
		  else
		    {
		      how_many_steps_in_3prime_dir_can_we_extend++;
		      j++;
		    }

		  //make sure we stop if we're at end of supernode
		  if ( (start_of_3prime_anchor_in_sup + k==1) || (start_of_3prime_anchor_in_sup + k==length_curr_supernode-2) )
		  {
		    can_extend_further_in_3prime_dir=false;
		  }

		  //sanity
		  if (start_of_3prime_anchor_in_chrom + length_3p_flank+ how_many_steps_in_3prime_dir_can_we_extend >= length_of_arrays-1)
		    {
		      //something has gone wrong
		      printf("We are extending the 3prime anchor in 3prime direction and have reached the end of our array. I don't believe it!!\n");
		      printf("start_of_3prime_anchor_in_chrom = %d\nhow_many_steps_in_3prime_dir_can_we_extend=%d\n, length_of_arrays=%d\n",
			     start_of_3prime_anchor_in_chrom, how_many_steps_in_3prime_dir_can_we_extend,  length_of_arrays);
		      exit(1);
		    }

		}		  



	      // collect info about how much bigger our 3p flank can be made:
	      length_3p_flank+= how_many_steps_in_5prime_dir_can_we_extend + how_many_steps_in_3prime_dir_can_we_extend;
	      start_of_3prime_anchor_in_chrom -= how_many_steps_in_5prime_dir_can_we_extend;

	      if (traverse_sup_left_to_right)
		{
		  start_of_3prime_anchor_in_sup -= how_many_steps_in_5prime_dir_can_we_extend;
		}
	      else
		{
		  start_of_3prime_anchor_in_sup += how_many_steps_in_5prime_dir_can_we_extend;
		}


	      printf("Collected following info:\n");
	      printf("start_of_3prime_anchor_in_sup = %d\n", start_of_3prime_anchor_in_sup);
	      printf("start_of_3prime_anchor_in_chrom = %d\n", start_of_3prime_anchor_in_chrom);
	      printf("Length 3p flank is  %d\n", length_3p_flank);
	      printf("Length 5p flank is  %d\n", length_5p_flank);
	      printf("Length of arrays is %d\n", length_of_arrays);
	      printf("Start node index is %d\n", start_node_index);
	      printf("first_index_in_chrom_where_supernode_differs_from_chromosome is %d\n", first_index_in_chrom_where_supernode_differs_from_chromosome);
	      printf("index_of_query_node_in_supernode_array is %d\n", index_of_query_node_in_supernode_array);
	      
	      
	      //sanity
	      if (start_of_3prime_anchor_in_chrom<= first_index_in_chrom_where_supernode_differs_from_chromosome)
		{
		  printf("Exit. start_of_3prime_anchor_in_chrom is %d and is less than first_index_in_chrom_where_supernode_differs_from_chromosome is %d. This is terrible. ", 
			 start_of_3prime_anchor_in_chrom, first_index_in_chrom_where_supernode_differs_from_chromosome);
		  printf("Start node index is %d, and coordinate in fasta file is %d\n", start_node_index, start_coord_of_variant_in_trusted_path_fasta);
		  exit(1);
		}
	      


	      char flank5p[length_5p_flank+1];
	      flank5p[0] = '\0';
	      strncpy(flank5p, chrom_string+start_node_index, length_5p_flank);
	      flank5p[length_5p_flank]='\0';
	      
	      int flank5p_max_covg=0;
	      int flank5p_min_covg=0;
	      double flank5p_avg_covg=0;
	      int flank5p_mode_coverage=0;
	      double flank5p_percent_nodes_with_modal_covg=0;
	      double flank5p_percent_novel=0; //if trusted path=reference, this is zero by definition, but we want to leave room for trusted path not the reference
                                              // a node is novel iff it is not in the reference genome

	      //get coverage of the 5p flank, and %novel sequence
	      if (traverse_sup_left_to_right==true)
		{
		  get_coverage_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array, 
						   length_5p_flank, &flank5p_min_covg, &flank5p_max_covg, &flank5p_avg_covg, &flank5p_mode_coverage, &flank5p_percent_nodes_with_modal_covg,
						   which_array_holds_indiv, index_for_indiv_in_edge_array);

		  //actually looks at the nodes and sees if they are in the graph of the person specified in the last 2 arguments - in this case, the reference => hence we find if novel
		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array, 
							length_5p_flank, &flank5p_percent_novel,
							which_array_holds_ref, index_for_ref_in_edge_array);
				



		}
	      else
		{
		  get_coverage_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank, 
						   length_5p_flank, &flank5p_min_covg, &flank5p_max_covg, &flank5p_avg_covg, &flank5p_mode_coverage, &flank5p_percent_nodes_with_modal_covg,
						   which_array_holds_indiv, index_for_indiv_in_edge_array);


		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank, 
							length_5p_flank, &flank5p_percent_novel,
							which_array_holds_ref, index_for_ref_in_edge_array);



		}

	      
	      

	      //print 5pflank
	      char name[300];
	      sprintf(name,"var_%i_5p_flank",num_variants_found);
	      print_fasta_from_path_for_specific_person_or_pop(output_file, name, length_5p_flank, 
							       flank5p_avg_covg,flank5p_min_covg,flank5p_max_covg, flank5p_mode_coverage, flank5p_percent_nodes_with_modal_covg, flank5p_percent_novel,
							       chrom_path_array[start_node_index],
							       chrom_orientation_array[start_node_index],
							       chrom_path_array[start_node_index+length_5p_flank-1],
							       chrom_orientation_array[start_node_index+length_5p_flank-1],
							       flank5p,
							       db_graph->kmer_size,
							       true,
							       which_array_holds_indiv, index_for_indiv_in_edge_array);

	      if (num_variants_found<=max_desired_returns)
		{
		  //add first kmer onto flank5p - this si done automatically by print_fasta_from_path
		  if (chrom_orientation_array[start_node_index]==forward)
		    {
		      strcat(return_flank5p_array[num_variants_found-1],
			     binary_kmer_to_seq(chrom_path_array[start_node_index]->kmer, db_graph->kmer_size, tmp_zam));
		    }
		  else
		    {
		      char tmp_zam2[db_graph->kmer_size+1];
		      tmp_zam2[db_graph->kmer_size]='\0';

		      strcat(return_flank5p_array[num_variants_found-1],
			     seq_reverse_complement(binary_kmer_to_seq(chrom_path_array[start_node_index]->kmer, db_graph->kmer_size, tmp_zam), 
						    db_graph->kmer_size,tmp_zam2)
			     );
		    }
		      
		  strcat(return_flank5p_array[num_variants_found-1],flank5p);
		}

		    

	      //print first branch, which is is from the trusted path
	      sprintf(name,"var_%i_branch1_trusted", num_variants_found);

	      trusted_branch[length_of_arrays]='\0';
	      trusted_branch[0]='\0';

	      int len_trusted_branch = start_of_3prime_anchor_in_chrom  - (first_index_in_chrom_where_supernode_differs_from_chromosome-1);
	      
	      //sanity
	      if ( (len_trusted_branch<=0) || (len_trusted_branch>length_of_arrays-1) )
		{
		  printf("len_trusted_branch is %d, and length of arrays is %d - this should never happen. Exit.\n", len_trusted_branch, length_of_arrays);
		  exit(1);;
		}

	      strncpy(trusted_branch, chrom_string+first_index_in_chrom_where_supernode_differs_from_chromosome-1, len_trusted_branch);
	      
	      trusted_branch[len_trusted_branch]='\0';
	      int trusted_branch_min_covg=0;
	      int trusted_branch_max_covg=0;
	      double trusted_branch_avg_covg=0;
	      int trusted_branch_mode_coverage=0;
	      double trusted_branch_percent_nodes_with_modal_covg=0;
	      double trusted_branch_percent_novel=0;

	      //get coverage on the trusted branch. 
	      get_coverage_from_array_of_nodes(chrom_path_array+start_node_index+length_5p_flank,
					       len_trusted_branch, &trusted_branch_min_covg, &trusted_branch_max_covg, &trusted_branch_avg_covg, &trusted_branch_mode_coverage, &trusted_branch_percent_nodes_with_modal_covg,
					       which_array_holds_indiv, index_for_indiv_in_edge_array);

	      get_percent_novel_from_array_of_nodes(chrom_path_array+start_node_index+length_5p_flank, 
						    len_trusted_branch, &trusted_branch_percent_novel,
						    which_array_holds_ref, index_for_ref_in_edge_array);
	      

	      print_fasta_from_path_for_specific_person_or_pop(output_file, name, len_trusted_branch, 
							       trusted_branch_avg_covg, trusted_branch_min_covg, trusted_branch_max_covg, 
							       trusted_branch_mode_coverage, trusted_branch_percent_nodes_with_modal_covg, trusted_branch_percent_novel,
							       chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
							       chrom_orientation_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
							       chrom_path_array[start_of_3prime_anchor_in_chrom],
							       chrom_orientation_array[start_of_3prime_anchor_in_chrom],
							       trusted_branch,
							       db_graph->kmer_size,
							       false,
							       which_array_holds_indiv, index_for_indiv_in_edge_array);
	      

	      if (num_variants_found<=max_desired_returns)
		{
		  strcat(return_trusted_branch_array[num_variants_found-1],trusted_branch);
		}

	      

	      //print second branch - this is the variant - which is the part of the supernode between the anchors/flanking regions
	      sprintf(name,"var_%i_branch2_variant", num_variants_found);

	      int branch2_min_covg=0;
	      int branch2_max_covg=0;
	      double branch2_avg_covg=0;
	      int branch2_mode_coverage=0;
	      double branch2_percent_nodes_with_modal_covg=0;
	      double branch2_percent_novel=0;

	      int len_branch2;

	      if (traverse_sup_left_to_right==true)
		{
		  len_branch2 = start_of_3prime_anchor_in_sup +1 - index_in_supernode_where_supernode_differs_from_chromosome;

		  //sanity. Note it is legitimate for len_branch2 to be zero
		  if (len_branch2<0)
		    {
		      printf("Traversing supernode left to right, but len_branch2 is %d, which is <=0 - this should be impossible. Exit\n", len_branch2);
		      exit(1);
		    }
		  else if (len_branch2> length_curr_supernode)
		    {
		      printf("Traversing supernode left to right, but len_branch2 is %d, which is > length of supernode %d - this should be impossible. Exit\n", len_branch2, length_curr_supernode);
		      exit(1);
		    }

		  //sanity
		  if (index_of_query_node_in_supernode_array+length_5p_flank+len_branch2>max_expected_size_of_supernode-length_3p_flank)
		    {
		      printf("Programming error. index_of_query_node_in_supernode_array+length_5p_flank+len_branch2 = %d is > max_expected_size_of_supernode-length_3p_flank %d, which should never happen.\n",
			     index_of_query_node_in_supernode_array+length_5p_flank+len_branch2,
			     max_expected_size_of_supernode-length_3p_flank );
		      printf("index_of_query_node_in_supernode_array is %d, length_5p_flank is %d, len_branch2 is %d, \n", index_of_query_node_in_supernode_array, length_5p_flank, len_branch2);
		      printf(" max_expected_size_of_supernode is %d, length_3p_flank is %d",  max_expected_size_of_supernode, length_3p_flank);
		      exit(1);
		    }

		  strncpy(variant_branch, supernode_string + index_of_query_node_in_supernode_array+length_5p_flank, len_branch2);
		  variant_branch[len_branch2]='\0';

		  
		  //get coverages
		  get_coverage_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank,
                                                   len_branch2, &branch2_min_covg, &branch2_max_covg, &branch2_avg_covg, &branch2_mode_coverage, &branch2_percent_nodes_with_modal_covg,
                                                   which_array_holds_indiv, index_for_indiv_in_edge_array);

		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank, 
							len_branch2, &branch2_percent_novel,
							which_array_holds_ref, index_for_ref_in_edge_array);




		  print_fasta_from_path_for_specific_person_or_pop(output_file, name, len_branch2, 
								   branch2_avg_covg, branch2_min_covg, branch2_max_covg, branch2_mode_coverage, branch2_percent_nodes_with_modal_covg, branch2_percent_novel, 
								   current_supernode[index_of_query_node_in_supernode_array+length_5p_flank], 
								   curr_sup_orientations[index_of_query_node_in_supernode_array+length_5p_flank],
								   current_supernode[start_of_3prime_anchor_in_sup], 
								   curr_sup_orientations[start_of_3prime_anchor_in_sup],
								   variant_branch, db_graph->kmer_size, false, 
								   which_array_holds_indiv, index_for_indiv_in_edge_array
								   );
		  if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_variant_branch_array[num_variants_found-1],variant_branch);
		    }
		  


		}
	      else //traversing supernode right to left
		{
		  len_branch2 = index_in_supernode_where_supernode_differs_from_chromosome +1 - start_of_3prime_anchor_in_sup ;


		  //sanity. Note it is legitimate for len_branch2 to be zero
		  if (len_branch2<0)
		    {
		      printf("Traversing supernode right to left, but len_branch2 is %d, which is <=0 - this should be impossible.\n", len_branch2);
		      printf("index_in_supernode_where_supernode_differs_from_chromosome is %d, and start_of_3prime_anchor_in_sup is %d. Exit.\n", 
			     index_in_supernode_where_supernode_differs_from_chromosome, start_of_3prime_anchor_in_sup);
		      exit(1);
		    }
		  else if (len_branch2> length_curr_supernode)
		    {
		      printf("Traversing supernode right to left, but len_branch2 is %d, which is > length of supernode %d - this should be impossible. Exit\n", len_branch2, length_curr_supernode);
		      exit(1);
		    }


		  // ok - need to take reverse complement of the full supernode sequence, inlcuding first kmer
		  char temp[max_expected_size_of_supernode+db_graph->kmer_size+1];
		  temp[0]='\0';
		  temp[max_expected_size_of_supernode+db_graph->kmer_size]='\0';

		  if (curr_sup_orientations[0]==forward)
		    {
		      strcat(temp, binary_kmer_to_seq(current_supernode[0]->kmer, db_graph->kmer_size, tmp_zam));
		    }
		  else
		    {
		      strcat(temp, binary_kmer_to_seq(binary_kmer_reverse_complement(current_supernode[0]->kmer, db_graph->kmer_size), 
						      db_graph->kmer_size, tmp_zam) );
		    }

		  strcat(temp, supernode_string);//temp is now the full sequence of the supernode.
		                                 //remember the query node, where the supernode touches the trusted path at start_node_index, 
                                                 // is not the 0-th element of the supernode. (definitely not, as we are traversing right to left!)



		  //sanity
		  if ( (index_of_query_node_in_supernode_array-length_5p_flank-len_branch2<0) || (index_of_query_node_in_supernode_array-length_5p_flank-len_branch2> max_expected_size_of_supernode) )
		    {
		      printf("Programming error. Either index_of_query_node_in_supernode_array-length_5p_flank-len_branch2<0\n");
		      printf(" or index_of_query_node_in_supernode_array-length_5p_flank-len_branch2> max_expected_size_of_supernode.\n");
		      printf("Neither should be possible. \n");
		      printf("index_of_query_node_in_supernode_array is %d, length_5p_flank is %d, len_branch2 is %d, max_expected_size_of_supernode is %d", 
			     index_of_query_node_in_supernode_array, length_5p_flank, len_branch2, max_expected_size_of_supernode );
		      exit(1);
		    }
		  strncpy(variant_branch, temp+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2, len_branch2);

		  variant_branch[len_branch2]='\0';
		  char rev_variant_branch[len_branch2+1];
		  seq_reverse_complement(variant_branch, len_branch2, rev_variant_branch);
		  rev_variant_branch[len_branch2]='\0';


		  //get coverages
		  get_coverage_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2,
                                                   len_branch2, &branch2_min_covg, &branch2_max_covg, &branch2_avg_covg, &branch2_mode_coverage, &branch2_percent_nodes_with_modal_covg,
                                                   which_array_holds_indiv, index_for_indiv_in_edge_array);



		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2, 
							len_branch2, &branch2_percent_novel,
							which_array_holds_ref, index_for_ref_in_edge_array);




		  print_fasta_from_path_for_specific_person_or_pop(output_file, name, len_branch2, 
								   branch2_avg_covg, branch2_min_covg, branch2_max_covg, branch2_mode_coverage, branch2_percent_nodes_with_modal_covg, branch2_percent_novel,
								   current_supernode[index_of_query_node_in_supernode_array-length_5p_flank], 
								   curr_sup_orientations[index_of_query_node_in_supernode_array-length_5p_flank],
								   current_supernode[start_of_3prime_anchor_in_sup], 
								   curr_sup_orientations[start_of_3prime_anchor_in_sup],
								   rev_variant_branch, db_graph->kmer_size, false, 
								   which_array_holds_indiv, index_for_indiv_in_edge_array
								   );


		  if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_variant_branch_array[num_variants_found-1],rev_variant_branch);
		    }
	
		}


	      //print 3p flank
	      char flank3p[length_3p_flank+1];
	      flank3p[0] = '\0';

	      //sanity
	      if (start_of_3prime_anchor_in_chrom<0)
		{
		  printf("Start of 3rime anchor in chrom <0 - impossible!\n");
		  exit(1);
		}
	      else if (start_of_3prime_anchor_in_chrom+ length_3p_flank>length_of_arrays)
		{
		  printf("Problem. start_of_3prime_anchor_in_chrom is %d and length_3p_flank is %d, and they add up to more thn length of arrays: %d\n", start_of_3prime_anchor_in_chrom, length_3p_flank, length_of_arrays);
		  exit(1);
		}
	      strncpy(flank3p, chrom_string+start_of_3prime_anchor_in_chrom, length_3p_flank);
	      flank3p[length_3p_flank]='\0';

	      int flank3p_max_covg=0;
	      int flank3p_min_covg=0;
	      double flank3p_avg_covg=0;
	      int flank3p_mode_coverage=0;
	      double flank3p_percent_nodes_with_modal_covg=0;
	      double flank3p_percent_novel=0;


	      if (traverse_sup_left_to_right==true)
		{
		  get_coverage_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank+len_branch2, 
						   length_3p_flank, &flank3p_min_covg, &flank3p_max_covg, &flank3p_avg_covg, &flank3p_mode_coverage, &flank3p_percent_nodes_with_modal_covg,
						   which_array_holds_indiv, index_for_indiv_in_edge_array);

		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank+len_branch2, 
							length_3p_flank, &flank3p_percent_novel,
							which_array_holds_ref, index_for_ref_in_edge_array);

				

		}
	      else
		{
		  get_coverage_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2-length_3p_flank, 
						   length_3p_flank, &flank3p_min_covg, &flank3p_max_covg, &flank3p_avg_covg,  &flank3p_mode_coverage, &flank3p_percent_nodes_with_modal_covg,
						   which_array_holds_indiv, index_for_indiv_in_edge_array);

		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2-length_3p_flank, 
							length_3p_flank, &flank3p_percent_novel,
							which_array_holds_ref, index_for_ref_in_edge_array);

			


		}


							       
	      sprintf(name,"var_%i_3p_flank",num_variants_found);
	      print_fasta_from_path_for_specific_person_or_pop(output_file, name, length_3p_flank, 
							       flank3p_avg_covg, flank3p_min_covg, flank3p_max_covg, flank3p_mode_coverage, flank3p_percent_nodes_with_modal_covg, flank3p_percent_novel, 
							       chrom_path_array[start_of_3prime_anchor_in_chrom],
							       chrom_orientation_array[start_of_3prime_anchor_in_chrom],
							       chrom_path_array[start_of_3prime_anchor_in_chrom+length_3p_flank],
							       chrom_orientation_array[start_of_3prime_anchor_in_chrom+length_3p_flank],
							       flank3p,
							       db_graph->kmer_size,
							       false,
							       which_array_holds_indiv, index_for_indiv_in_edge_array);

	      if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_flank3p_array[num_variants_found-1],flank3p);
		    }
	

	      start_node_index = start_of_3prime_anchor_in_chrom+length_3p_flank+1;
 
	    }
	 
	}

      
      
      coord_of_start_of_array_in_trusted_fasta+=number_of_nodes_to_load;
      
    }  while (db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta
	                    (chrom_fasta_fptr, number_of_nodes_to_load, number_of_nodes_to_load,
			     length_of_arrays, 
			     chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
			     seq, kmer_window, 
			     false, true, db_graph)>0);


  //free malloc-ed variables

  free(chrom_path_array);
  free(chrom_orientation_array);
  free(chrom_labels);
  free(chrom_string);
  free(current_supernode);
  free(curr_sup_orientations);
  free(curr_sup_labels);
  free(supernode_string);
  free_sequence(&seq);
  free(kmer_window->kmer);
  free(kmer_window);
  free(trusted_branch);
  free(variant_branch);
													   
  return num_variants_found;


}


void compute_label(dBNode * node, Orientation o, char * label, EdgeArrayType type, int index){

  int i=0; 

  void check_nucleotide(Nucleotide n){

    if (db_node_edge_exist(node,n,o, type, index)){
	label[i] =  binary_nucleotide_to_char(n);
	i++;
      }
  }
  
  nucleotide_iterator(check_nucleotide);

  label[i]='\0';
}


void print_fasta_from_path_for_specific_person_or_pop(FILE *fout,
						      char * name,
						      int length,
						      double avg_coverage,
						      int min_coverage,
						      int max_coverage,
						      int modal_coverage,
						      double percent_nodes_with_modal_coverage, 
						      double percent_novel,
						      dBNode * fst_node,
						      Orientation fst_orientation,
						      dBNode * lst_node,
						      Orientation lst_orientation,
						      char * string, //labels of paths
						      int kmer_size,
						      boolean include_first_kmer,
						      EdgeArrayType type,
						      int index
						      ){


  if (fout==NULL)
    {
      printf("Exiting - have passed a null file pointer to print_fasta_from_path_for_specific_person_or_pop\n");
      exit(1);
    }
  
  if ( (fst_node==NULL) || (lst_node==NULL) )
    {
      fprintf(fout, "WARNING - print_fasta command as have been given a NULL node as first or last node - will not print fst/lst kmers.\n");
      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i mode_coverage: %i percent_nodes_with_modal_covg: %5.2f percent_novel: %5.2f\n", 
              name, (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage, modal_coverage, percent_nodes_with_modal_coverage, percent_novel);
      fprintf(fout,"%s\n",string);
    }
  else
    {
      char fst_f[5], fst_r[5],lst_f[5],lst_r[5];
      
      compute_label(fst_node,forward,fst_f, type, index);
      compute_label(fst_node,reverse,fst_r, type, index);
      compute_label(lst_node,forward,lst_f, type, index);
      compute_label(lst_node,reverse,lst_r, type, index);
      
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
      
      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i mode_coverage: %i percent_nodes_with_modal_covg: %5.2f percent_novel: %5.2f fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s\n", name,
	      (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage, modal_coverage, percent_nodes_with_modal_coverage, percent_novel,
	      db_node_get_coverage(fst_node, type, index),
	      fst_seq,
	      (fst_orientation == forward ? fst_r : fst_f),
	      (fst_orientation == forward ? fst_f : fst_r),
	      db_node_get_coverage(lst_node, type, index),
	      lst_seq,
	      (lst_orientation == forward ? lst_r : lst_f),
	      (lst_orientation == forward ? lst_f : lst_r));


      if (include_first_kmer){    
	
	fprintf(fout,"%s",binary_kmer_to_seq(fst_kmer,kmer_size,fst_seq));
      }
      
      fprintf(fout,"%s\n",string);
      
    }


}


/*
void get_covg_of_nodes_in_one_but_not_other_of_two_arrays(const dBNode** const array1, const dBNode** const array2, int length1, int length2, 
							  int* num_nodes_in_array_1not2, int * num_nodes_in_array_2not1, int** covgs_in_1not2, int** covgs_in_2not1, //results
							  const dBNode** const array_for_working1, const dBNode** const array_for_working2)
{

  //get copies of these arrays into the arrays which we shall sort
  int i;
  for (i=0; i<length1; i++)
    {
      array_for_working1[i] = array1[i];
    }
  for (i=0; i<length2; i++)
    {
      array_for_working2[i] = array2[i];
    }  
  qsort (array_for_working1, length1, sizeof(dBNode*), int_cmp); 
  qsort (array_for_working2, length2, sizeof(dBNode*), int_cmp); 

  int j=0;
  i=0;

  for (i=0; i<length1; i++)
    {
      for (j=0; j<length2; j++)
	{
	  if (int_cmp(i,j)==0)
	    {
	    }
	  else if (int_cmp(i,j)>0)
	    {
iqbal
	      
	    }
	}
    }

}


*/
