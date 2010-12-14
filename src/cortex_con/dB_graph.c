/*
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */

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
#include <file_reader.h>

void next_node_sanity_check(dBNode * current_node, Nucleotide nucleotide, Orientation current_orientation,dBNode * next_node, dBGraph * db_graph){
  char tmp_seq[db_graph->kmer_size+1]; 
  if(next_node == NULL){
    fprintf(stderr,"dB_graph: didnt find node in hash table: %s %c %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),binary_nucleotide_to_char(nucleotide), current_orientation == forward ? "forward" : "reverse");
    exit(1);
  }
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
			   boolean include_first_kmer);




//it doesn't check that it is a valid arrow -- it just assumes the arrows is fine
dBNode * db_graph_get_next_node(dBNode * current_node, Orientation current_orientation, 
			       Orientation * next_orientation,
			       Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph){

  BinaryKmer local_copy_of_kmer;
  binary_kmer_assignment_operator(local_copy_of_kmer, current_node->kmer);

  BinaryKmer tmp_kmer;
  dBNode * next_node=NULL;

  // after the following line tmp_kmer and rev_kmer are pointing to the same B Kmer
  BinaryKmer* rev_kmer = binary_kmer_reverse_complement(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer);

  
  if (current_orientation == reverse){   
    *reverse_edge = binary_kmer_get_last_nucleotide(&local_copy_of_kmer);
    binary_kmer_assignment_operator(local_copy_of_kmer,*rev_kmer);
  }
  else{
    *reverse_edge = binary_kmer_get_last_nucleotide(rev_kmer);
  }


  binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&local_copy_of_kmer, edge, db_graph->kmer_size);

   //get node from table
  next_node = hash_table_find(element_get_key(&local_copy_of_kmer,db_graph->kmer_size, &tmp_kmer),db_graph);

  if (next_node != NULL){
    *next_orientation = db_node_get_orientation(&local_copy_of_kmer,next_node,db_graph->kmer_size);
  }

  return next_node;
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
       next_node_sanity_check(node,nucleotide,orientation,next_node,db_graph);

       //if(next_node == NULL){
       // printf("dB_graph_db_node_clip_tip_with_orientation: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
       // exit(1);
       //}	           
       
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
       db_node_reset_edge(next_node,reverse_nucleotide,opposite_orientation(next_orientation));
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
	next_node_sanity_check(node,nucleotide,forward,next_node,db_graph);
	db_node_reset_edge(next_node,reverse_nucleotide,opposite_orientation(next_orientation));
	  
      }
	
      if (db_node_edge_exist(node,nucleotide,reverse)){
	next_node = db_graph_get_next_node(node,reverse,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
	next_node_sanity_check(node,nucleotide,reverse,next_node,db_graph);
	db_node_reset_edge(next_node,reverse_nucleotide,opposite_orientation(next_orientation));
	
      }
    }

    nucleotide_iterator(&nucleotide_action);

    node_action(node);
    db_node_reset_edges(node);	
   

  }

  return ret;
}




// computes a random path until a stop signal is found
// ie the starting node can have  multiple exits
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only
// returns length of the path (as number of labels)

int db_graph_get_random_path(dBNode * node, Orientation orientation, int limit, 
			     void (*node_action)(dBNode * node),
			     dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
			     char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
			     int * count_fresh_nodes,
			     boolean * is_cycle, dBGraph * db_graph){

  Orientation  current_orientation,next_orientation;
  dBNode * current_node = NULL;
  dBNode * next_node = NULL;
  Nucleotide nucleotide,rev_nucleotide;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  tmp_seq[db_graph->kmer_size]='\0';
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
    

  if (db_node_get_one_edge(node,orientation,&nucleotide)==true){

    do{ 
      if (length>0){
	if (db_node_check_flag_visited(current_node) == false){
	  (*count_fresh_nodes)++;
	}

	db_node_set_flag(current_node,VISITED);

	if (current_orientation==forward){
	  db_node_set_flag(current_node,VISITED_FF);
	}
	else{
	  db_node_set_flag(current_node,VISITED_RR);
	}

	sum_coverage += coverage;
	*max_coverage = (*max_coverage < coverage) ? coverage : *max_coverage;
	*min_coverage = (*min_coverage > coverage) ? coverage : *min_coverage;
      }
      
      //printf("Path %i!\n",nucleotide);
      next_node =  db_graph_get_next_node(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph);
      next_node_sanity_check(current_node,nucleotide,current_orientation,next_node,db_graph);
      
      //sanity check
      //if(next_node == NULL){
      //	fprintf(stderr,"dB_graph: didnt find node in hash table: %s %c %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),binary_nucleotide_to_char(nucleotide), current_orientation == forward ? "forward" : "reverse");
      //exit(1);
      //}
      
      
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
	     db_node_check_for_flag(next_node, STOP_PATH)==false &&
	     !(db_node_check_for_flag(next_node,VISITED_FF)==true && next_orientation==forward) &&
	     !(db_node_check_for_flag(next_node,VISITED_RR)==true && next_orientation==reverse) &&
	     db_node_get_one_edge(next_node,next_orientation,&nucleotide)==true);
  
    if ((next_node == node) && (next_orientation == orientation)){
      *is_cycle = true;
    }

    int i;
    for(i=1;i<length;i++){
      db_node_unset_flag(path_nodes[i], VISITED_FF | VISITED_RR);
      node_action(path_nodes[i]);
     
    }
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
  tmp_seq[db_graph->kmer_size]='\0';
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
      *max_coverage = (*max_coverage < coverage) ? coverage : *max_coverage;
      *min_coverage = (*min_coverage > coverage) ? coverage : *min_coverage;
    }

    next_node =  db_graph_get_next_node(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph);
    next_node_sanity_check(current_node,nucleotide,current_orientation,next_node,db_graph);
      
    //sanity check
    //if(next_node == NULL){
    // fprintf(stderr,"dB_graph: didnt find node in hash table: %s %c %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),binary_nucleotide_to_char(nucleotide), current_orientation == forward ? "forward" : "reverse");
    // exit(1);
    //}

    
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
      
      if (path_nodes[j][lengths[j]] == path_nodes[k][lengths[k]] && 
	  path_orientations[j][lengths[j]] == path_orientations[k][lengths[k]]){
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






// a dirty bubble starts in a node with at least two outgoing edges in the same orientation
// every branch of the bubble is not necesssary a perfect path 
// we just look for a couple of path that meet. 

boolean db_graph_detect_dirty_bubble(dBNode * node,
				     Orientation orientation,
				     int depth, int breadth, //breadth is the number of paths allowed 
				     void (*node_action)(dBNode * node), 
				     int * length1,dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,
				     char * seq1, double * avg_coverage1, int * min_coverage1, int * max_coverage1,
				     int * length2,dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,
				     char * seq2, double * avg_coverage2, int * min_coverage2, int * max_coverage2,
				     dBGraph * db_graph){



  int lengths[breadth]; //lengths of paths in the search space
  boolean full_path[breadth]; //paths that have been completed (ie when reaching a fork a path is completed
  dBNode * path_nodes[breadth][depth];
  Orientation path_orientations[breadth][depth];
  Nucleotide path_labels[breadth][depth];
  int new_paths;
  int number_paths=0;
  int cdepth = 1;

  
  void extend_path(int i, Nucleotide n){
    //at this point we know path can be extended with n
    //printf("\textend path %i with %c\n",i,binary_nucleotide_to_char(n));
    Orientation next_orientation;
    Nucleotide rev_nucleotide;
    dBNode * next_node;
    //printf("\t\tget next node \n");
    next_node =  db_graph_get_next_node(path_nodes[i][lengths[i]],path_orientations[i][lengths[i]],&next_orientation,n,&rev_nucleotide,db_graph);
    //printf("\t\tdone with get next node\n");
    lengths[i]++;
    path_nodes[i][lengths[i]] = next_node;
    path_orientations[i][lengths[i]] = next_orientation;
    path_labels[i][lengths[i]-1] = n;
    //printf("\tdone with extesion\n");
    
  }

  void do_one_path(int i){
    int edge_count_f = db_node_edges_count(path_nodes[i][lengths[i]],path_orientations[i][lengths[i]]);
    int edge_count_r = db_node_edges_count(path_nodes[i][lengths[i]],opposite_orientation(path_orientations[i][lengths[i]]));
    int pos = lengths[i];
    //printf("\tdo one path %i %i - path %i pos %i\n",edge_count_f,edge_count_r, i , pos);
    if (edge_count_f==0){
      full_path[i]=true;
    }
    else{
      if (edge_count_f==1 && edge_count_r ==1){ //extend
	void do_one_step(Nucleotide n){
	  if (db_node_edge_exist(path_nodes[i][pos],n,path_orientations[i][pos])){
	    extend_path(i,n);
	  }
	}	
	nucleotide_iterator(do_one_step);
      }
      else{ //create copies
	void create_path_and_extend(Nucleotide n){
	  if (db_node_edge_exist(path_nodes[i][pos],n,path_orientations[i][pos])){
	    //copy path to number_path+new_paths index
	    //printf("CHECK %i %i %i\n",i,number_paths,pos);
	    int j;
	    for(j=0;j<=pos;j++){	 
	      path_nodes[number_paths+new_paths][j] =  path_nodes[i][j];
	      path_orientations[number_paths+new_paths][j] =  path_orientations[i][j];
	      if (j<pos){path_labels[number_paths+new_paths][j] =  path_labels[i][j];}
	    }
	    lengths[number_paths+new_paths]=lengths[i];
	    extend_path(number_paths+new_paths,n);
	    new_paths++;
	  }
	}
	
	nucleotide_iterator(create_path_and_extend);
	full_path[i]=true;
      }
    }

  }
  

  void create_path(int i, int *length, dBNode ** nodes, Orientation * orientations, Nucleotide * labels, char * seq, double * avg_coverage, int * min_coverage, int * max_coverage){
    *max_coverage = 0;
    *min_coverage = INT_MAX;
    *length = lengths[i];
    int sum_coverage;
    int j;
    for(j=0;j<=lengths[i];j++){
      nodes[j]=path_nodes[i][j];
      node_action(path_nodes[i][j]);
      orientations[j] = path_orientations[i][j];
      if (j<lengths[i]){labels[j] = path_labels[i][j];}
      //printf("path %i pos %i\n",i,j);
      if (j<lengths[i]){seq[j] = binary_nucleotide_to_char(path_labels[i][j]);}

      int coverage = element_get_coverage(path_nodes[i][j]);
      sum_coverage += coverage;
      *max_coverage = (*max_coverage < coverage) ? coverage : *max_coverage;
      *min_coverage = (*min_coverage > coverage) ? coverage : *min_coverage;
    }

    nodes[j]=path_nodes[i][j];
    orientations[j]=path_orientations[i][j];
    seq[j-1]='\0';

    *avg_coverage = (length-1<=0) ? 0 : (double) sum_coverage/(double) (*length-1);
    if (*min_coverage == INT_MAX)
      {
	*min_coverage = 0;
      };
  }
  



  boolean bubble_found = false;
  if (db_node_edges_count(node,orientation)>1){   

    void create_first_step(Nucleotide n){
      if (db_node_edge_exist(node,n,orientation)){

	path_nodes[number_paths][0]=node;
	lengths[number_paths]=0;
	path_orientations[number_paths][0]=orientation;
	extend_path(number_paths,n);
	full_path[number_paths]=false;
	number_paths++;
      }
    }
	
    nucleotide_iterator(create_first_step);


    int j,k;
    
    while( (!bubble_found) && (number_paths<=breadth) && (cdepth<=depth) ){
    
      new_paths=0;
      int i;
      for(i=0;i<number_paths;i++){
	//extend/copy paths
	if (!full_path[i]){
	  do_one_path(i);
	}
      }
      //printf("\tfinish one step -  %i\n",new_paths);
      number_paths+=new_paths;
      //printf("\tLOOK1 %i\n",cdepth);
      //check for bubble
      //printf("\tcheck for bubbles\n");
      for(j=0;j<number_paths;j++){
	
	for(k=j+1;k<number_paths;k++){
	  if ((path_nodes[j][lengths[j]] == path_nodes[k][lengths[k]]) &&
	      (path_orientations[j][lengths[j]] == path_orientations[k][lengths[k]])){
	    //printf("\tbubble found!!\n");
	    //found bubble
	    //printf("\tcreate path1\n");
	    create_path(j,length1,path_nodes1,path_orientations1,path_labels1,seq1,avg_coverage1,min_coverage1,max_coverage1);
	    //printf("\tcreate path2\n");
	    create_path(k,length2,path_nodes2,path_orientations2,path_labels2,seq2,avg_coverage2,min_coverage2,max_coverage2);
	    
	    bubble_found = true;
	  }
	}
      }
      cdepth++;
      //printf("\tLOOK2 %i\n",cdepth);
    }
  }
   
 return bubble_found;
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
    //let's redo the last step, we need to do that because the last node could have had multiple entries
    
    Nucleotide label;
    Orientation next_orientation;
    
    dBNode * lst_node = db_graph_get_next_node(nodes_reverse[length_reverse-1],orientations_reverse[length_reverse-1],
					       &next_orientation, labels_reverse[length_reverse-1],&label,db_graph);

    //sanity checks
    next_node_sanity_check(nodes_reverse[length_reverse-1],labels_reverse[length_reverse-1],orientations_reverse[length_reverse-1],lst_node,db_graph);
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


void db_graph_detect_vars(int max_length, dBGraph * db_graph, int choice);

void db_graph_detect_vars_clean_bubbles(int max_length, dBGraph * db_graph){
  db_graph_detect_vars(max_length,db_graph,0);
}

void db_graph_detect_vars_dirty_bubbles(int max_length, dBGraph * db_graph){
  db_graph_detect_vars(max_length,db_graph,1);
}


//routine to get short indels and SNPs (delta = 0 necessary condition, not sufficient as we could also find balancing duplications or similar)

void db_graph_detect_vars(int max_length, dBGraph * db_graph, int clean_dirty){
  
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
	
	boolean go_ahead = false;

	if (clean_dirty == 0){
	 
	  go_ahead = db_graph_detect_bubble(current_node,orientation,max_length,&db_node_action_set_flag_visited, 
					    &length1,path_nodes1,path_orientations1,path_labels1,
					    seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
					    &length2,path_nodes2,path_orientations2,path_labels2,
					    seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
					    db_graph);
	}
	else{	  

	  go_ahead = db_graph_detect_dirty_bubble(current_node,
						  orientation,
						  max_length,10,
						  &db_node_action_set_flag_visited,
						  &length1,path_nodes1,path_orientations1,path_labels1,
						  seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
						  &length2,path_nodes2,path_orientations2,path_labels2,
						  seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
						  db_graph);
	}

	if (go_ahead == false){
	}
	else{
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
	  count_vars++;
	  printf("\nVARIATION: %i - coverage: %d\n",count_vars,element_get_coverage(current_node));
	  
	  //compute 5' flanking region       	
	  int length_flank5p_reverse = db_graph_get_perfect_path(current_node,opposite_orientation(orientation),
								 flanking_length,&db_node_action_set_flag_visited,
								 nodes5p,orientations5p,labels_flank5p,
								 seq5p,&avg_coverage5p,&min5p,&max5p,
								 &is_cycle5p,db_graph);
	  
	  if (length_flank5p_reverse>0){
	    Nucleotide label;
	    Orientation next_orientation;
	    dBNode * next_node;

	    next_node = db_graph_get_next_node(nodes5p[length_flank5p_reverse-1],orientations5p[length_flank5p_reverse-1],
					       &next_orientation, labels_flank5p[length_flank5p_reverse-1],&label,db_graph);
	    
	    //sanity check
	    next_node_sanity_check(nodes5p[length_flank5p_reverse-1],labels_flank5p[length_flank5p_reverse-1],orientations5p[length_flank5p_reverse-1],next_node,db_graph);
	    
	    length_flank5p = db_graph_get_perfect_path_with_first_edge(nodes5p[length_flank5p_reverse],
								       opposite_orientation(orientations5p[length_flank5p_reverse]),
								       flanking_length,label,
								       &db_node_action_set_flag_visited,
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
						     flanking_length,&db_node_action_set_flag_visited,
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
	  

	  db_node_action_set_flag_visited(path_nodes2[length2]);
	  
	  current_node = path_nodes2[length2];
	  orientation = path_orientations2[length2];
	}
      } while (current_node != node //to avoid cycles
	       && db_node_check_flag_visited(current_node)==false);
    }
    
  
    if (db_node_check_flag_visited(node)==false){     
      db_node_action_set_flag_visited(node);
      get_vars_with_orientation(node,forward);
      get_vars_with_orientation(node,reverse);
    }
  }
 
  hash_table_traverse(&get_vars,db_graph); 
    
}




int db_graph_remove_bubbles(int limit, dBGraph * db_graph){


  int removed_nodes=0;
  int count_bubbles=0;

  void remove_bubble_with_orientation(dBNode * node, Orientation orientation){
    int length1, length2;
    dBNode * path_nodes1[limit+1];
    dBNode * path_nodes2[limit+1];
    
    Orientation path_orientations1[limit+1];
    Orientation path_orientations2[limit+1];

    Nucleotide path_labels1[limit];
    Nucleotide path_labels2[limit];

    char seq1[limit+1];
    char seq2[limit+1];

    double avg_coverage1;
    int min_coverage1,max_coverage1;	  
    double avg_coverage2;
    int min_coverage2,max_coverage2;
	  
    if (db_node_edges_count(node,orientation)==2 &&
	db_graph_detect_bubble(node,orientation,limit,&db_node_action_do_nothing,
			       &length1,path_nodes1,path_orientations1,path_labels1,
			       seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
			       &length2,path_nodes2,path_orientations2,path_labels2,
			       seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
			       db_graph) == true){

      count_bubbles++;
      //prune branch -- select branch with less coverage
      dBNode * * path_to_prune;
      Orientation * orientations_to_prune;
      Nucleotide * labels_to_prune;
      int length_to_prune;

      if (element_get_coverage(path_nodes1[1])>element_get_coverage(path_nodes2[1])){
	length_to_prune = length2;
	path_to_prune = path_nodes2;
	labels_to_prune = path_labels2;
	orientations_to_prune = path_orientations2;
      }
      else{
	length_to_prune = length1;
	path_to_prune = path_nodes1;
	labels_to_prune = path_labels1;
	orientations_to_prune = path_orientations1;
      }


      int i;
      dBNode * next_node;
      Orientation next_orientation;
      Nucleotide reverse_nucleotide;

      //remove connection at start of loop
      db_node_reset_edge(path_to_prune[0],labels_to_prune[0],orientations_to_prune[0]);

      //remove connection at end of loop
      next_node = db_graph_get_next_node(path_to_prune[length_to_prune-1],
					 orientations_to_prune[length_to_prune-1],
					 &next_orientation,
					 labels_to_prune[length_to_prune-1],
					 &reverse_nucleotide,db_graph);

      //sanity checks
      next_node_sanity_check(path_to_prune[length_to_prune-1],labels_to_prune[length_to_prune-1],orientations_to_prune[length_to_prune-1],next_node,db_graph);

      if (next_node != path_to_prune[length_to_prune]){
	printf("dB_graph_remove_bubbles: didnt find node\n");
      }
      if (db_node_edge_exist(path_to_prune[length_to_prune],reverse_nucleotide,opposite_orientation(orientations_to_prune[length_to_prune]))==false){
	printf("dB_graph_remove_bubbles: didnt find edge in hash table - length1: %i length2: %i\n",length1,length2);
      }
      
      db_node_reset_edge(path_to_prune[length_to_prune],reverse_nucleotide,opposite_orientation(orientations_to_prune[length_to_prune]));
      

      for(i=1;i<length_to_prune;i++){
	removed_nodes++;
	db_node_reset_edges(path_to_prune[i]);
	db_node_action_set_flag_pruned(path_to_prune[i]);
      }        
    }

  }

  void remove_bubble(dBNode * node){
    remove_bubble_with_orientation(node,forward);
    remove_bubble_with_orientation(node,reverse);
  }


  hash_table_traverse(&remove_bubble,db_graph);
  printf("%i bubbles removed\n",count_bubbles);
  return removed_nodes;
}

//threshold is the length of the tips
int db_graph_clip_tips(int threshold, dBGraph * db_graph){
  int nodes_clipped = 0;

  void clip_tips(dBNode * node){
    
    if (db_node_check_for_flag(node,PRUNED)==false){
      nodes_clipped += db_graph_db_node_clip_tip(node, threshold,&db_node_action_set_flag_pruned,db_graph);
    }
  }

  hash_table_traverse(&clip_tips,db_graph);
  
  return nodes_clipped; 
}

//threshold is the coverage 
int db_graph_clip_low_coverage_supernodes(int threshold,int threshold_length, int limit,dBGraph * db_graph){
  
  int supernodes_removed=0;
  int nodes_removed=0;
  
  int length;
  dBNode * path_nodes[limit+1];
  Orientation path_orientations[limit+1];
  Nucleotide path_labels[limit];
  char seq[limit+1];
  boolean is_cycle;
  double avg_coverage;
  int min;
  int max;
 
  void clip_supernode(dBNode * node){
 
    if (db_node_check_for_flag(node,VISITED_CLEANING)==false && element_get_coverage(node)<=threshold){

      length = db_graph_supernode(node,limit,&db_node_action_set_flag_visited_cleaning,
				  path_nodes,path_orientations,path_labels,
				  seq,&avg_coverage,&min,&max,&is_cycle,
				  db_graph);

      if ((is_cycle == false) && (length >1) && (length<threshold_length)){
	//check that all internal nodes in the supernode have coverage smaller than threshold
	int i;
	boolean remove = true;
	for(i=1;i<length;i++){
	  if (element_get_coverage(path_nodes[i])>threshold){
	    remove = false;
	    break;
	  }
	  
	}
	
	if (remove == true){
	   dBNode * next_node;
	   Orientation next_orientation;
	   Nucleotide reverse_nucleotide;
	  
	  //remove supernode 
	  supernodes_removed++;
	  //remove connection from first node
	  db_node_reset_edge(path_nodes[0],path_labels[0],path_orientations[0]);
	  
	  
	  //remove connection at end of path
	  next_node = db_graph_get_next_node(path_nodes[length-1],
					     path_orientations[length-1],
					     &next_orientation,
					     path_labels[length-1],
					     &reverse_nucleotide,db_graph);
	  

	  //sanity checks
	  next_node_sanity_check(path_nodes[length-1], path_labels[length-1], path_orientations[length-1],next_node,db_graph);

	  if (next_node != path_nodes[length]){
	    printf("db_graph_clip_low_coverage_supernodes: didnt find node\n");
	  }
	  else{
	    if (db_node_edge_exist(path_nodes[length],reverse_nucleotide,opposite_orientation(path_orientations[length]))==false){
	      printf("db_graph_clip_low_coverage_supernodes: didnt find edge in hash table\n");
	    }
	    else{
	      if (next_orientation != path_orientations[length]){
		printf("db_graph_clip_low_coverage_supernodes: inconsistent orientation %s %s\n",next_orientation == forward ? "forward" : "reverse",opposite_orientation(path_orientations[length]) == forward ? "forward" : "reverse");
	      }
	      else{
		db_node_reset_edge(path_nodes[length],reverse_nucleotide,opposite_orientation(path_orientations[length]));
	      }
	    }
	  }

	  for(i=1;i<length;i++){
	    nodes_removed++;
	    db_node_reset_edges(path_nodes[i]);
	    db_node_action_set_flag_pruned(path_nodes[i]);	  
	  }
	}
      }
    }
  } 
	  

  hash_table_traverse(&clip_supernode,db_graph);
  printf("%i supernodes removed %i kmers_removed\n",supernodes_removed,nodes_removed);
  return nodes_removed;
}


void db_graph_print_paths(char * filename, int max_length,boolean with_coverages, int coverage_threshold, dBGraph * db_graph){

  FILE * fout; 
  FILE * fout_cov; 

  fout= fopen(filename, "w"); 

  if (with_coverages){
    char filename_cov[strlen(filename)+10];
    sprintf(filename_cov,"%s_cov",filename);
    fout_cov = fopen(filename_cov,"w");
  }

  int count_nodes=0;
  
  dBNode * *    path_nodes_r;
  dBNode * *    path_nodes_f;
  Orientation * path_orientations;
  Nucleotide *  path_labels;
  char * seq_aux;
  char * seq;
  boolean is_cycle;
  double avg_coverage_r,avg_coverage_f;
  int min_r,max_r,min_f,max_f;
  
  
  path_nodes_r        = calloc(max_length,sizeof(dBNode*));
  path_nodes_f        = calloc(max_length,sizeof(dBNode*));
  path_orientations = calloc(max_length,sizeof(Orientation));
  path_labels       = calloc(max_length,sizeof(Nucleotide));
  seq_aux           = calloc(max_length+1,sizeof(char));
  seq               = calloc(2*max_length+1,sizeof(char));
  
  

  long long count_kmers = 0;


  void print_path(dBNode * node){
   

    count_kmers++;
    
    char name[100];

    if ((db_node_check_flag_visited(node) == false) && (element_get_coverage(node)<=coverage_threshold)){  
      
      int count_fresh_nodes=0;
      db_node_action_set_flag_visited(node);

      sprintf(name,"node_%i",count_nodes);
      

      //printf("get reverse\n");
      int length_r = db_graph_get_random_path(node,reverse,max_length,&db_node_action_set_flag_visited,
					      path_nodes_r,path_orientations,path_labels,
					      seq_aux,&avg_coverage_r,&min_r,&max_r,&count_fresh_nodes,&is_cycle,
					      db_graph);
      
      dBNode * fst_node = path_nodes_r[length_r];
      Orientation fst_orientation =  opposite_orientation(path_orientations[length_r]);
      db_node_action_set_flag_visited(fst_node);
      
  
      //printf("compute reverse complement %i %s\n",length_r,name);
      seq_reverse_complement(seq_aux, strlen(seq_aux), seq);
      //fprintf(fout,"%s",seq_f);     
      seq_aux[0]='\0';


      //get first kmer
      BinaryKmer fst_kmer;
      char fst_seq[db_graph->kmer_size+1];
      binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(node)));
      binary_kmer_to_seq(&fst_kmer,db_graph->kmer_size,fst_seq);
      //printf("R: %s\n",seq);
      //printf("K: %s\n",fst_seq);
      strcat(seq,fst_seq);
      
      //printf("get forward\n");
      int length_f = db_graph_get_random_path(node,forward,max_length,&db_node_action_set_flag_visited,
					      path_nodes_f,path_orientations,path_labels,
					      seq_aux,&avg_coverage_f,&min_f,&max_f,&count_fresh_nodes,&is_cycle,
					      db_graph);

      db_node_action_set_flag_visited(path_nodes_f[length_f]);
     
      //printf("F: %s\n",seq_aux);
      strcat(seq,seq_aux);
      if (count_fresh_nodes>db_graph->kmer_size+1){

	if (with_coverages){
	  fprintf(fout_cov,">%s length:%i length_string:%i\n",name,length_r+length_f+1,strlen(seq));
	  int i;
	  for(i=length_r;i>=0;i--){
	    fprintf(fout_cov,"%i ",element_get_coverage(path_nodes_r[i]));
	  }
	  
	  for(i=1;i<length_f;i++){
	    fprintf(fout_cov,"%i ",element_get_coverage(path_nodes_f[i]));
	  }
	  fprintf(fout_cov,"\n");
	}

	print_fasta_from_path(fout,name,length_r+length_f+db_graph->kmer_size,(avg_coverage_f+avg_coverage_r)/2,
			      ((min_r<min_f && min_r!=0))?min_r:min_f,
			      (max_r>max_f)?max_r:max_f,
			      fst_node,fst_orientation,
			      path_nodes_f[length_f],path_orientations[length_f],
			      seq,db_graph->kmer_size,
			      false);
	
	//printf("node %s fresh nodes %i length %i\n", name, count_fresh_nodes, length_r+length_f);
	count_nodes++;
      }
      seq_aux[0]='\0';
      seq[0]='\0';

    }
  }
    
  hash_table_traverse(&print_path,db_graph); 

  free(path_nodes_r);
  free(path_nodes_f);
  free(path_orientations);
  free(path_labels);
  free(seq_aux);
  free(seq);
  fclose(fout);
  if (with_coverages){
    fclose(fout_cov);
  }
}



//prints perfect paths

void db_graph_print_supernodes(char * filename, int max_length, boolean with_coverages, dBGraph * db_graph){

  FILE * fout; 
  FILE * fout_cov; 
  fout= fopen(filename, "w"); 

  if (with_coverages){
    char filename_cov[strlen(filename)+10];
    sprintf(filename_cov,"%s_cov",filename);
    fout_cov = fopen(filename_cov,"w");
  }

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
  char name[100];
  void print_supernode(dBNode * node){
    
    Nucleotide n;

    count_kmers++;
    
    

    if (db_node_check_flag_visited(node) == false){
      if (db_node_has_precisely_one_edge(node,forward,&n)==false || db_node_has_precisely_one_edge(node,forward,&n)==false){
	db_node_set_flag(node,STOP_PATH);
	db_node_set_flag(node,VISITED);
      }
      else{
	int length = db_graph_supernode(node,max_length,&db_node_action_set_flag_visited,
					path_nodes,path_orientations,path_labels,
					seq,&avg_coverage,&min,&max,&is_cycle,
					db_graph);
	
 
	
	if (length>0){	
	  sprintf(name,"node_%i",count_nodes);
	  if (with_coverages){
	    fprintf(fout_cov,">%s\n",name);
	    int i;
	    for(i=0;i<db_graph->kmer_size-1;i++){
	      fprintf(fout_cov,"%i ",element_get_coverage(path_nodes[0]));
	    }
	    
	    for(i=db_graph->kmer_size;i<length;i++){
	      fprintf(fout_cov,"%i ",element_get_coverage(path_nodes[i-db_graph->kmer_size+1]));
	    }
	    fprintf(fout_cov,"\n");
	  }
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
  }


  void print_complex_connections(dBNode * node){
    Orientation next_o;
    Nucleotide rev_n;
     
      void check_nucleotide_orientation(Nucleotide n,Orientation o){
      if (db_node_edge_exist(node,n,o)){
	dBNode * next_node = db_graph_get_next_node(node,o,&next_o,n,&rev_n,db_graph);
	//sanity checks
	next_node_sanity_check(node,n,o,next_node,db_graph);

	if (db_node_check_for_flag(next_node,STOP_PATH)==true){
	  //print path
	  sprintf(name,"node_%i",count_nodes);
	  char kmer_str[2];
	  kmer_str[0]=binary_nucleotide_to_char(n);
	  kmer_str[1]='\0';

	  print_fasta_from_path(fout,name,1,0,0,0,node,o,next_node,next_o,kmer_str,db_graph->kmer_size,true);
	  count_nodes++;

	    
	}
	
      }
    }

    if (db_node_check_for_flag(node,STOP_PATH)){
      nucleotide_iterator_orientation(&check_nucleotide_orientation);
      
    }
  }
 
  hash_table_traverse(&print_supernode,db_graph); 
  printf("%qd nodes visted [%qd singletons]\n",count_kmers,count_sing);
  
  hash_table_traverse(&print_complex_connections,db_graph); 

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(seq);
  fclose(fout);
  if (with_coverages){
    fclose(fout_cov);
  }
}


long long db_graph_remove_weak_edges(int threshold,dBGraph * db_graph){

  long long edges_removed=0;

  void find_weak_edge_with_orientation(dBNode * node,Orientation orientation){
    dBNode * next_node1;
    dBNode * next_node2;
    Nucleotide nucleotide1, nucleotide2, reverse_nucleotide1,reverse_nucleotide2;
    Orientation next_orientation1, next_orientation2;

    if (db_node_has_precisely_two_edges(node,orientation,&nucleotide1,&nucleotide2) && db_node_edges_count(node,opposite_orientation(orientation))==1){
       next_node1 = db_graph_get_next_node(node,orientation,&next_orientation1,nucleotide1,&reverse_nucleotide1,db_graph);
       next_node_sanity_check(node,nucleotide1,orientation,next_node1,db_graph);
       
       next_node2 = db_graph_get_next_node(node,orientation,&next_orientation2,nucleotide2,&reverse_nucleotide2,db_graph);
       next_node_sanity_check(node,nucleotide2,orientation,next_node2,db_graph);
       
       if (element_get_coverage(node)>threshold){
	 if (element_get_coverage(node)>element_get_coverage(next_node1) && 
	     element_get_coverage(node)-element_get_coverage(next_node1)==1 &&
	     element_get_coverage(next_node2)>threshold
	     ){
	   //remove arrow node->next_node2
	   edges_removed++;
	   printf("coverage %i coverage1 %i coverage2 %i\n",element_get_coverage(node),element_get_coverage(next_node1),element_get_coverage(next_node2));
	   db_node_reset_edge(node,nucleotide2,orientation);
	   db_node_reset_edge(next_node2,reverse_nucleotide2,opposite_orientation(next_orientation2));
	   
	 }
	 else{
	   if (element_get_coverage(node)>element_get_coverage(next_node2) && 
	       element_get_coverage(node)-element_get_coverage(next_node2)==1 &&
	       element_get_coverage(next_node1)>threshold
	       ){
	    //remove arrow node->next_node1
	     edges_removed++;
	     printf("coverage %i coverage1 %i coverage2 %i\n",element_get_coverage(node),element_get_coverage(next_node1),element_get_coverage(next_node2));
	     db_node_reset_edge(node,nucleotide1,orientation);
	     db_node_reset_edge(next_node1,reverse_nucleotide1,opposite_orientation(next_orientation1));
	   }
	 }
	 
       }
    }
  }

  void remove_weak_edge(dBNode * node){
    
    find_weak_edge_with_orientation(node,forward);
    find_weak_edge_with_orientation(node,reverse);
  }
  hash_table_traverse(&remove_weak_edge,db_graph); 
  return edges_removed;
}


long long db_graph_find_double_Ys(int max_length, dBGraph * db_graph){
  dBNode * *    path_nodes;
  Orientation * path_orientations;
  Nucleotide *  path_labels;
  char * seq;
  boolean is_cycle;
  double avg_coverage;
  int min,max; 
  long long count_double_Ys=0;
  
  path_nodes        = calloc(max_length,sizeof(dBNode*));
  path_orientations = calloc(max_length,sizeof(Orientation));
  path_labels       = calloc(max_length,sizeof(Nucleotide));
  seq               = calloc(max_length+1+db_graph->kmer_size,sizeof(char));
  
  
  boolean db_node_is_double_Y_with_orientation(dBNode * node, Orientation orientation){
    boolean ret = false;
    int length;

    if (db_node_edges_count(node,orientation)>1){
      length = db_graph_get_perfect_path(node,opposite_orientation(orientation),
					 max_length,&db_node_action_do_nothing,
					 path_nodes,path_orientations,path_labels,
					 seq,&avg_coverage,&min,&max,
					 &is_cycle,db_graph);
      
      if (length>0){
	if (is_cycle == false && 
	    db_node_edges_count(path_nodes[length],opposite_orientation(path_orientations[length]))==1 &&
	    db_node_edges_count(path_nodes[length],path_orientations[length])>1){
	  
	  db_node_set_flag(node,STOP_PATH);
	  ret = true;
	}
      }
      else{
	if (db_node_edges_count(node,opposite_orientation(orientation))>1){
	  db_node_set_flag(node,STOP_PATH);
	  ret = true;
	}
      }
    }

    return ret;
  }

  void db_node_is_double_Y(dBNode * node){

    if (db_node_is_double_Y_with_orientation(node,forward)==true){
      count_double_Ys++;
    }
    else{
      if (db_node_is_double_Y_with_orientation(node,reverse)==true){
	count_double_Ys++;
      }
    }

   
  }

  hash_table_traverse(&db_node_is_double_Y,db_graph); 

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(seq);
  return count_double_Ys;

}





//check all edges in graph
long long db_graph_health_check(boolean fix, dBGraph * db_graph){
  dBNode * next_node;
  Nucleotide reverse_nucleotide;
  Orientation next_orientation;
  long long count_nodes;
  char tmp_seq[db_graph->kmer_size+1]; 

  void check_node_with_orientation(dBNode * node, Orientation orientation){
  

    void check_base(Nucleotide n){     
      if (db_node_edge_exist(node,n,orientation)==true){
	next_node = db_graph_get_next_node(node,orientation,&next_orientation,n,&reverse_nucleotide,db_graph);
	
	if(next_node == NULL){
	  printf("Health check problem -  didnt find node in hash table: %s %c %s\n",
		 binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq),
		 binary_nucleotide_to_char(n), 
		 orientation == forward ? "forward" : "reverse");
	  if (fix){
	      db_node_reset_edge(node,n,orientation);
	    }
	}
	else{
	  if (db_node_edge_exist(next_node,reverse_nucleotide,opposite_orientation(next_orientation))==false){
	    printf("Health check problem - inconsitency return arrow missing: %s %c %s\n",
		   binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq),
		   binary_nucleotide_to_char(n), 
		   orientation == forward ? "forward" : "reverse");
	    if (fix){
	      db_node_reset_edge(node,n,orientation);
	    }
	  }
	}
      
      }
    }
    nucleotide_iterator(&check_base);
  }
  
  void check_node(dBNode * node){
    check_node_with_orientation(node,forward);
    check_node_with_orientation(node,reverse);
    count_nodes++;
  }
    
  hash_table_traverse(&check_node,db_graph); 
  printf("%qd nodes checked\n",count_nodes);
  return count_nodes;
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

int db_graph_remove_low_coverage_nodes(int coverage, dBGraph * db_graph){
  int pruned_nodes;

  void prune_node(dBNode * node){
    if (db_graph_db_node_prune_low_coverage(node,coverage,
			   &db_node_action_set_flag_pruned,
					    db_graph)){
      pruned_nodes++;
    }
  }

  hash_table_traverse(&prune_node,db_graph); 
  return pruned_nodes;
}




void db_graph_dump_binary(char * filename, boolean (*condition)(dBNode * node), dBGraph * db_graph){
  FILE * fout; //binary output
  fout= fopen(filename, "w"); 
  
  print_binary_signature(fout,db_graph->kmer_size,1);

  

  long long count=0;
  //routine to dump graph
  void print_node_binary(dBNode * node){   
    if (condition(node)){
      count++;
      db_node_print_binary(fout,node,db_graph->kmer_size);
    }
  }

  hash_table_traverse(&print_node_binary,db_graph); 
  fclose(fout);

  printf("%qd kmers dumped\n",count);
}

void db_graph_dump_hash_table(char * filename, dBGraph * db_graph){
  FILE * fout; //binary output
  fout= fopen(filename, "w"); 
  
  print_hash_table_signature(fout,db_graph);

  hash_table_dump_to_file(fout,db_graph);
  
  fclose(fout);

  printf("%qd kmers dumped (full table)\n",hash_table_get_capacity(db_graph));
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
  fst_seq[kmer_size]='\0';
  lst_seq[kmer_size]='\0';

  BinaryKmer fst_kmer;
  BinaryKmer tmp_kmer;

  binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)));

  if (fst_orientation==reverse){
    binary_kmer_reverse_complement(&fst_kmer,kmer_size, &tmp_kmer);
    binary_kmer_assignment_operator(fst_kmer, tmp_kmer);
  } 
  binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq);

  BinaryKmer lst_kmer;
  binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)));

  if (lst_orientation==reverse){
    binary_kmer_reverse_complement(&lst_kmer,kmer_size, &tmp_kmer);
    binary_kmer_assignment_operator(lst_kmer, tmp_kmer);
  } 
  binary_kmer_to_seq(&lst_kmer,kmer_size,lst_seq);

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


  //print string
  int start = 0;
  if (include_first_kmer){   
    binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq);
    if (kmer_size<80){
      fprintf(fout,"%s",fst_seq);
    }
    else{
      while(start+80<kmer_size){
	fprintf(fout,"%.*s\n",80,&fst_seq[start]);
	start+=80;	      
      }
      if (kmer_size % 80>0){
	fprintf(fout,"%.*s",kmer_size % 80,&fst_seq[start]);
      }
    }

    //complete tis line with sequence
    fprintf(fout,"%.*s\n",80-(kmer_size % 80),&string[0]);
    start = 80-(kmer_size % 80);
  }
  
  int str_length = strlen(string);

  while(start+80<str_length){
    fprintf(fout,"%.*s\n",80,&string[start]);
    start+=80;
  }

  if (str_length-start>0){
    fprintf(fout,"%.*s\n",str_length-start,&string[start]);
  }

  

}


