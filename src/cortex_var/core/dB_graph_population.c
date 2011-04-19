/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo  
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
  dB_graph_population.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>

#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
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
						      char* text_describing_comparison_with_other_path, //may be NULL - use to allow printing coverages of nodes in this path but not in some specific other path
						      int* coverages_nodes_in_this_path_but_not_some_other, //may be NULL
						      int length_of_coverage_array,
						      char * string, //labels of paths
						      int kmer_size,
						      boolean include_first_kmer,
						      EdgeArrayType type,
						      int index
						      );

void print_fasta_from_path_in_subgraph_defined_by_func_of_colours(FILE *fout,
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
								  char* text_describing_comparison_with_other_path, 
								     //text_describing_comparison_with_other_path may be NULL 
								     // - use to allow printing coverages of nodes in this path but not in some specific other path
								  int* coverages_nodes_in_this_path_but_not_some_other, //may be NULL
								  int length_of_coverage_array,
								  char * string, //labels of paths
								  int kmer_size,
								  boolean include_first_kmer,
								  Edges (*get_colour)(const dBNode*),
								  int (*get_covg)(const dBNode*)
								  );

void print_minimal_fasta_from_path_for_specific_person_or_pop(FILE *fout,
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
							      boolean include_first_kmer,
							      EdgeArrayType type,
							      int index
							      );

void print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(FILE *fout,
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
									  boolean include_first_kmer,
									  Edges (*get_colour)(const dBNode*),
									  int (*get_covg)(const dBNode*)
									  );

                                                                          
//This function does not  check that it there is such an edge in the specified person/colour - but it does check if the target node is in the specific person.
//if you want to be sure the dge exists in that colour, then check it before calling this function
dBNode * db_graph_get_next_node_for_specific_person_or_pop(dBNode * current_node, Orientation current_orientation, 
							   Orientation * next_orientation,
							   Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, EdgeArrayType type, int index)
{

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

  if (next_node != NULL)
    {
      *next_orientation = db_node_get_orientation(&local_copy_of_kmer,next_node,db_graph->kmer_size);
    }
  else
    {
      //no else
    }


  //need to check the node is in this person's graph
  if (! (db_node_is_this_node_in_this_person_or_populations_graph(next_node, type, index)))
    {
      return NULL;
    }
 
  return next_node;
  

}



//This function does not  check that it there is such an edge in the specified person/colour - but it does check if the target node is in the specific person.
//if you want to be sure the dge exists in that colour, then check it before calling this function
//The last argument allows you to apply the operation to some subgraph - eg you might take the unuiion of colours 2 and 3, or of all colours.
//If for example you wanted to get the next node in the graph irrespective of colour, get_colour would return the union (bitwise AND) of all edges in a node.
dBNode * db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(dBNode * current_node, Orientation current_orientation, 
								       Orientation * next_orientation,
								       Nucleotide edge, Nucleotide * reverse_edge,dBGraph * db_graph, 
								       Edges (*get_colour)(const dBNode*)
								       )
{

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

  if (next_node != NULL)
    {
      *next_orientation = db_node_get_orientation(&local_copy_of_kmer,next_node,db_graph->kmer_size);
    }
  else
    {
      //no else
    }


  //need to check the node is in the specified subgraph graph
  if (! (db_node_is_this_node_in_subgraph_defined_by_func_of_colours(next_node, get_colour)) )
    {
      return NULL;
    }
 
  return next_node;
  

}


//MARIO Mario - noode seems to call or use this?
boolean db_graph_db_node_has_precisely_n_edges_with_status_for_specific_person_or_pop(dBNode * node,Orientation orientation,NodeStatus status,int n,
										      dBNode * * next_node, Orientation * next_orientation, Nucleotide * next_base,
										      dBGraph * db_graph, EdgeArrayType type, int index)
{

  if ( (n>4)||(n<0))
    {
      printf("WARNING - calling db_graph_db_node_has_precisely_n_edges_with_status with n=%d", n);
      return false;
    }

  int count = 0;
  boolean ret = false;

  void check_edge(Nucleotide base){

    dBNode * tmp_next_node;
    Orientation tmp_next_orientation;
    Nucleotide rev_base;

   

    if (db_node_edge_exist(node,base,orientation, type, index)){

      tmp_next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&tmp_next_orientation,
									base,&rev_base,db_graph, type, index);

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


int db_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit,
									  void (*node_action)(dBNode * node),dBGraph * db_graph, EdgeArrayType type, int index)
{ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  int i;
  dBNode * nodes[limit];
  Orientation next_orientation;
  dBNode * next_node;
  char seq[db_graph->kmer_size+1];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end(node, opposite_orientation(orientation), type, index)){
    
    boolean join_main_trunk = false;

    while(db_node_has_precisely_one_edge(node,orientation,&nucleotide, type, index)) {
  
       nodes[length] = node;

       next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph, type, index);
       
       if(next_node == NULL){
	 printf("dB_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
	 exit(1);
       }	           
       
       length ++;

       
       if (length>limit){
	 break;
       }

       //want to stop when we join another trunk
       if (!db_node_has_precisely_one_edge(next_node,opposite_orientation(next_orientation),&nucleotide, type, index)){
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
	 db_node_reset_edges(nodes[i], type, index);
       }

       if (DEBUG){
	printf("RESET %c BACK\n",binary_nucleotide_to_char(reverse_nucleotide));
      }
       db_node_reset_edge(next_node,opposite_orientation(next_orientation),reverse_nucleotide, type, index);
     }

  }

  return length;
}




// clip a tip in the graph (the tip starts in node)
// limit is max length for tip
// node_action is applied to all the elements in the tip
// returns the length of the tip (0 means no length)
int db_graph_db_node_clip_tip_for_specific_person_or_pop(dBNode * node, int limit,
							 void (*node_action)(dBNode * node),
							 dBGraph * db_graph, EdgeArrayType type, int index)
{

  int length_tip = 0;

  
  length_tip = db_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop(node,forward,limit,node_action,db_graph, type, index);
  
  if (length_tip==0){
    length_tip = db_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop(node,reverse,limit,node_action,db_graph, type, index);
    
  }
  
  return length_tip;
}




int db_graph_db_node_clip_tip_with_orientation_in_subgraph_defined_by_func_of_colours(dBNode * node, Orientation orientation, int limit,
										       void (*node_action)(dBNode * node),dBGraph * db_graph, 
										       Edges (*get_colour)(const dBNode*),
										      void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
										       void (*apply_reset_to_colour)(dBNode*)
										       )
{ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  int i;
  dBNode * nodes[limit];
  Orientation next_orientation;
  dBNode * next_node;
  char seq[db_graph->kmer_size+1];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end_in_subgraph_given_by_func_of_colours(node, opposite_orientation(orientation), get_colour))
    {
      boolean join_main_trunk = false;

      while(db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(node,orientation,&nucleotide, get_colour)) 
	{
	  nodes[length] = node;
	  
	  next_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph, get_colour);
	  
	  if(next_node == NULL){
	    printf("dB_graph_db_node_clip_tip_with_orientation_for_specific_person_or_pop: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
	    exit(1);
	  }	           
	  
	  length ++;
	  
	  
	  if (length>limit){
	    break;
	  }
	  
	  //want to stop when we join another trunk
	  if (!db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(next_node,opposite_orientation(next_orientation),&nucleotide, get_colour)){
	    join_main_trunk = true;
	    break;
	  }

	  //keep track of node
	  
	  node = next_node;
	  orientation = next_orientation;        
	}
      
      if (! join_main_trunk)
	{
	  length = 0;
	}
      else
	{//clear edges and mark nodes as pruned
	  for(i=0;i<length;i++)
	    {
	      if (DEBUG){
		printf("CLIPPING node: %s\n",binary_kmer_to_seq(element_get_kmer(nodes[i]),db_graph->kmer_size,seq));
	      }
	      
	      node_action(nodes[i]);
	      //perhaps we want to move this inside the node action?
	      apply_reset_to_colour(nodes[i]);
	    }
	  
	  if (DEBUG){
	    printf("RESET %c BACK\n",binary_nucleotide_to_char(reverse_nucleotide));
	  }
	  apply_reset_to_specific_edge_in_colour(next_node,opposite_orientation(next_orientation),reverse_nucleotide);
     }
      
    }
  
  return length;
}




// clip a tip in the graph (the tip starts in node)
// limit is max length for tip
// node_action is applied to all the elements in the tip
// returns the length of the tip (0 means no length)
int db_graph_db_node_clip_tip_in_subgraph_defined_by_func_of_colours(dBNode * node, int limit,
								     void (*node_action)(dBNode * node),
								     dBGraph * db_graph, 
								     Edges (*get_colour)(const dBNode*),
								     void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
								     void (*apply_reset_to_colour)(dBNode*)
								     )
{

  int length_tip = 0;

  
  length_tip = db_graph_db_node_clip_tip_with_orientation_in_subgraph_defined_by_func_of_colours(node,forward,limit,node_action,db_graph, get_colour, apply_reset_to_specific_edge_in_colour, apply_reset_to_colour );
  
  if (length_tip==0){
    length_tip = db_graph_db_node_clip_tip_with_orientation_in_subgraph_defined_by_func_of_colours(node,reverse,limit,node_action,db_graph, get_colour, apply_reset_to_specific_edge_in_colour, apply_reset_to_colour);
    
  }
  
  return length_tip;
}






// 1. the argument sum_of_covgs_in_desired_colours allows you to choose which colour (or union of colours) you want to apply this to
// eg you might want  to "remove" (as defined by the action) any node that has coverage <= your threshold in the UNION of all colours, or in colour red or whatever
// SO - this func returns the sum of coverages in the colours you care about
// 2. the argument get_edge_of_interest is a function that gets the "edge" you are interested in - may be a single edge/colour from the graph, or might be a union of some edges
// 3 Pass apply_reset_to_specified_edges which applies reset_one_edge to whichever set of edges you care about, 
// 4 Pass apply_reset_to_specified_edges_2 which applies db_node_reset_edges to whichever set of edges you care about, 
boolean db_graph_db_node_prune_low_coverage(dBNode * node, int coverage,
					    void (*node_action)(dBNode * node),
					    dBGraph * db_graph, 
					    int (*sum_of_covgs_in_desired_colours)(const Element *), 
					    Edges (*get_edge_of_interest)(const Element*),
					    void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide),
					    void (*apply_reset_to_specified_edges_2)(dBNode*) )
{

  boolean ret = false;

  if (sum_of_covgs_in_desired_colours(node)<=coverage){
    ret = true;
  
    void nucleotide_action(Nucleotide nucleotide){
      Orientation next_orientation;
      Nucleotide reverse_nucleotide;
      dBNode * next_node;


      //remove evidence in adjacent nodes where they point to this one
      if (db_node_edge_exist_within_specified_function_of_coloured_edges(node,nucleotide,forward, get_edge_of_interest)){
	next_node = db_graph_get_next_node(node,forward,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
	db_node_reset_specified_edges(next_node,opposite_orientation(next_orientation),reverse_nucleotide, apply_reset_to_specified_edges);
	  
      }
	
      if (db_node_edge_exist_within_specified_function_of_coloured_edges(node,nucleotide,reverse, get_edge_of_interest)){
	next_node = db_graph_get_next_node(node,reverse,&next_orientation,nucleotide,&reverse_nucleotide,db_graph);
	db_node_reset_specified_edges(next_node,opposite_orientation(next_orientation),reverse_nucleotide, apply_reset_to_specified_edges);
	
      }
    }

    nucleotide_iterator(&nucleotide_action);

    node_action(node);
    //remove all edges from this node, in colours we care about
    apply_reset_to_specified_edges_2(node);   

  }

  return ret;
}


boolean db_graph_db_node_prune_low_coverage_ignoring_colours(dBNode * node, int coverage,
							     void (*node_action)(dBNode * node),
							     dBGraph * db_graph)
{

  Edges get_edge_of_interest(const Element* node)
  {
    return get_union_of_edges(*node);
  }
  
  Edges set_to_zero(Edges edge)
  {
    return 0;
  }

  void apply_reset_to_specified_edges(dBNode* node , Orientation or , Nucleotide nuc)
  {
      int i;
      for (i=0; i< NUMBER_OF_COLOURS; i++)
	{
	  reset_one_edge(node, or, nuc, individual_edge_array, i);
	}
  }
  
  void apply_reset_to_specified_edges_2(dBNode* node)
  {
      int i;
      for (i=0; i< NUMBER_OF_COLOURS; i++)
	{
	  db_node_reset_edges(node, individual_edge_array, i);
	}
  }

  int sum_of_covgs_in_desired_colours(const dBNode * node)
  {
    int total=0;
    int i;
    for (i=0; i< NUMBER_OF_COLOURS; i++)
      {
	total += db_node_get_coverage(node, individual_edge_array,i);
      }
    return total;
  }

  return db_graph_db_node_prune_low_coverage(node, coverage, node_action, db_graph,
					     &sum_of_covgs_in_desired_colours,
					     &get_edge_of_interest,
					     &apply_reset_to_specified_edges,
					     &apply_reset_to_specified_edges_2);


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
  tmp_seq[db_graph->kmer_size]='\0';
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
      *max_coverage = (*max_coverage < coverage) ? coverage : *max_coverage;
      *min_coverage = (*min_coverage > coverage) ? coverage : *min_coverage;
    }

    //this will return NULL if the next node is not in the person's graph. It does NOT check if the edge is in the graph...
    next_node =  db_graph_get_next_node_for_specific_person_or_pop(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph, type, index);

      


    //sanity check
    if(next_node == NULL)
      {
	fprintf(stderr,"db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop: didnt find node in hash table: %s %c %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),binary_nucleotide_to_char(nucleotide), current_orientation == forward ? "forward" : "reverse");
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
      printf("WARNING IMPOSSIBLE");
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

//get_colour returns an edge which is a function of the edges in a node. eg union of all edges, or of colours 1 and 2
// get covg returns coverge desired - most ikely is sum of covgs in each of colours which are summed/unioned in get_colour
int db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(dBNode * node, Orientation orientation, int limit, 
										     Nucleotide fst_nucleotide,
										     void (*node_action)(dBNode * node),
										     dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
										     char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
										     boolean * is_cycle, dBGraph * db_graph, 
										     Edges (*get_colour)(const dBNode*),
										     int (*get_covg)(const dBNode*) 
										     )
{

  Orientation  current_orientation,next_orientation;
  dBNode * current_node = NULL;
  dBNode * next_node = NULL;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  tmp_seq[db_graph->kmer_size]='\0';
  int sum_coverage = 0;
  int coverage  = 0;

  //sanity checks
  if (node == NULL)
    {
      printf("db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours: can't pass a null node\n");
      exit(1);
    }
  else if (! (db_node_is_this_node_in_subgraph_defined_by_func_of_colours(node, get_colour)))
    {
      printf("\nThis node is not in the graph of this person - in db_graph_get_perfect_path_\n");
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
      *max_coverage = (*max_coverage < coverage) ? coverage : *max_coverage;
      *min_coverage = (*min_coverage > coverage) ? coverage : *min_coverage;
    }

    //this will return NULL if the next node is not in the subgraph specified by get_colour. 
    next_node =  db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(current_node,current_orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph, get_colour);

      


    //sanity check
    if(next_node == NULL)
      {
	fprintf(stderr,"db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours: didnt find node in hash table: %s %c %s\n", binary_kmer_to_seq(element_get_kmer(current_node),db_graph->kmer_size,tmp_seq),binary_nucleotide_to_char(nucleotide), current_orientation == forward ? "forward" : "reverse");
	exit(1);
      }


    path_labels[length]        = nucleotide;
    seq[length]                = binary_nucleotide_to_char(nucleotide);
    coverage                   = db_node_get_coverage_in_subgraph_defined_by_func_of_colours(next_node, get_covg);
    
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
	   db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(next_node,opposite_orientation(next_orientation),&nucleotide2, get_colour) && //multiple entries
	   db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(current_node, current_orientation,&nucleotide, get_colour)); //has one next edge only
  
  
  if ((next_node == node) && (next_orientation == orientation)){
    *is_cycle = true;
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


// perfect path -- no conflict no cycle -- returns length
// path: node_0 edge_0 node_1 edge_1 ... node_n-1 edge_n-1 node_n
// path_nodes is a n+1 array from 0..n with all the nodes in the path
// path_orientations is n+1 array from 0..n with the orientations of the node in the path
// path labels is n array from 0..n-1 with all the labels for the edges (see above)
// node_action only applied to internal nodes (not node_0 and node_n)
// seq is a string with all the labels concatenated (NB: it doesn't contain the kmer in the first node)
// avg_coverage, min_coverge, max_coverge -> refers to the internal nodes only

int db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(dBNode * node, Orientation orientation, int limit, 
								     void (*node_action)(dBNode * node),
								     dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
								     char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
								     boolean * is_cycle, dBGraph * db_graph, 
								     Edges (*get_colour)(const dBNode*),
								     int (*get_covg)(const dBNode*) )
{

  int length =0;
  Nucleotide nucleotide;

  //sanity check
  if (node == NULL){
    printf("db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours: can't pass a null node\n");
    exit(1);
  }

  path_nodes[0]         = node;
  path_orientations[0]  = orientation;  
 

  if (db_node_has_precisely_one_edge_in_subgraph_defined_by_func_of_colours(node,orientation,&nucleotide, get_colour))
    {
    
      length= db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(node,orientation, limit, nucleotide,
											       node_action,
											       path_nodes,path_orientations,path_labels,
											       seq,avg_coverage,min_coverage,max_coverage,
											       is_cycle,db_graph, 
											       get_colour, get_covg);
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

// the algorithms works by trying to find supernodes for all possible edges in the given orientation.
// the action is applied to all the internal nodes in the supernodes (even if they don't form a bubble).
// at the moment returns only one bubble, notice that potentially there could be situations where a node forms two bubbles (something to implement). 


//TODO - could have triallelic
boolean db_graph_detect_bubble_for_specific_person_or_population(dBNode * node,
								  Orientation orientation,
								  int limit,
								  void (*node_action)(dBNode * node), 
								  int * length1,dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,
								  char * seq1, double * avg_coverage1, int * min_coverage1, int * max_coverage1,
								  int * length2,dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,
								  char * seq2, double * avg_coverage2, int * min_coverage2, int * max_coverage2,
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
    if (db_node_edge_exist(node,n,orientation, type, index)){
	lengths[i] = db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(node,orientation,
											 limit,n,
											 node_action,
											 path_nodes[i],path_orientations[i],path_labels[i],
											 seq[i],&avg_coverage[i],&min_coverage[i],&max_coverage[i],
											 &is_cycle,db_graph, type, index);
	
	i++;
    }
  }
  
  
  nucleotide_iterator(check_nucleotide);
 
  int j = 0;
  int k = 0;

  while (!ret && j<i){
    k = j+1;
    while (!ret && k<i){
      
      if ( (path_nodes[j][lengths[j]] == path_nodes[k][lengths[k]])
	   &&
	   (path_orientations[j][lengths[j]] == path_orientations[k][lengths[k]])
	   )
	{
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



// a bubble starts in a node with only two outgoing edges in the same orientation
// every branch of the bubble is free of conflicts until it joins with the main branch
// the differente between the length of the branches is <= delta

// the algorithms works by trying to find supernodes for all possible edges in the given orientation.
// the action is applied to all the internal nodes in the supernodes (even if they don't form a bubble).
// at the moment returns only one bubble, notice that potentially there could be situations where a node forms two bubbles (something to implement). 

// sometimes, you want to mark branches of a bubble specially, in which case pass true for apply_special_action_to_branches, and a unction pointer to special_action.
// otherwise, pass false, NULL for the last two actions.

//TODO - could have triallelic
boolean db_graph_detect_bubble_in_subgraph_defined_by_func_of_colours(dBNode * node,
								      Orientation orientation,
								      int limit,
								      void (*node_action)(dBNode * node), 
								      int * length1,dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,
								      char * seq1, double * avg_coverage1, int * min_coverage1, int * max_coverage1,
								      int * length2,dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,
								      char * seq2, double * avg_coverage2, int * min_coverage2, int * max_coverage2,
								      dBGraph * db_graph, Edges (*get_colour)(const dBNode*),  int (*get_covg)(const dBNode*),
								      boolean apply_special_action_to_branches, void (*special_action)(dBNode * node))
{

  if (node==NULL)
    {
      printf("Do not call db_graph_detect_bubble_in_subgraph_defined_by_func_of_colours with NULL node. Exiting\n");
      exit(1);
    }
  else if (! (db_node_is_this_node_in_subgraph_defined_by_func_of_colours(node, get_colour)))
    // this will include pruned nodes
    {
      //printf("\nThis node is not in the graph of this person\n");
      return false;
    }

  
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
    //check in all four directions, and in the process apply action (usually setting to visited)
    if (db_node_edge_exist_within_specified_function_of_coloured_edges(node,n,orientation, get_colour)){
	lengths[i] = db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(node,orientation,
												      limit,n,
												      node_action,
												      path_nodes[i],path_orientations[i],path_labels[i],
												      seq[i],&avg_coverage[i],&min_coverage[i],&max_coverage[i],
												      &is_cycle,db_graph, get_colour, get_covg);
	
	i++;
    }
  }

  // i is now a count of how many edges there are  
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
	//OK - we are done here. We have found a bubble, and put its data into the relevant arrays, so happy to return
	//First we apply the special action to the branches of the bubble we have found


	if (apply_special_action_to_branches==true)
	  {
	    for(l=0;l<lengths[j]; l++)
	      {
		//printf("Applying special action %d\n",l); 
		special_action(path_nodes1[l]);
	      }
	    for(l=0;l<lengths[k]; l++)
	      {
		//printf("Applying special action %d\n",l); 
		special_action(path_nodes2[l]);
	      }
	  }
	  
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


boolean db_graph_db_node_smooth_bubble_for_specific_person_or_pop(dBNode * node, Orientation orientation, 
								  int limit,int coverage_limit,
								  void (*node_action)(dBNode * node),
								  dBGraph * db_graph, EdgeArrayType type, int index){

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


  if (db_graph_detect_bubble_for_specific_person_or_population(node,orientation,limit,&db_node_action_do_nothing,
							       &length1,path_nodes1,path_orientations1,path_labels1,
							       seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
							       &length2,path_nodes2,path_orientations2,path_labels2,
							       seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
							       db_graph, type, index)){
    

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
	db_node_reset_edge(node,orientation,path_labels1[0], type, index);
      }
      else
	{
	  path_nodes_tmp = path_nodes2;
	  length_tmp     = length2;
	  path_orientations_tmp = path_orientations2;
	  path_labels_tmp       = path_labels2;
	  db_node_reset_edge(node,orientation,path_labels2[0], type, index);
	}
      
      for(i=1;i<length_tmp;i++){
	node_action(path_nodes_tmp[i]);
	db_node_reset_edges(path_nodes_tmp[i], type, index);
      }
      
      db_graph_get_next_node_for_specific_person_or_pop(path_nodes_tmp[length_tmp-1],path_orientations_tmp[length_tmp-1],&next_orientation,path_labels_tmp[length_tmp-1],&reverse_nucleotide,
							db_graph, type, index);
      db_node_reset_edge(path_nodes_tmp[length_tmp],opposite_orientation(next_orientation),reverse_nucleotide, type, index);
    }
    
  }

  return ret;
}




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


// it returns the supernode containing 'node'  
// string has to support limit+1 (+1 as you need a space for the \0 at the end)
// node_action has to be idempotent as it can be applied to the same node twice!!
// supernode_str returns the string made of the labels of the path (doesn't include first kmer). 
// returns the length of supernode 
// last two arguments mean you find supernode in a subgraph of the main graph defined by the edges as returned by get_colour.
int db_graph_supernode_in_subgraph_defined_by_func_of_colours(dBNode * node,int limit,void (*node_action)(dBNode * node), 
							      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char * supernode_str, 
							      double * avg_coverage,int * min,int * max, boolean * is_cycle, 
							      dBGraph * db_graph, 
							      Edges (*get_colour)(const dBNode*),
							      int (*get_covg)(const dBNode*)  ){

  
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
  
  length_reverse = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(node,reverse,limit,&db_node_action_do_nothing,
										    nodes_reverse,orientations_reverse,labels_reverse,
										    supernode_str,&avg_coverager,&minr,&maxr,
										    &is_cycler,db_graph, get_colour, get_covg);
  
  if (length_reverse>0){
    //let's re do the last step, we need to do that because the last node could have had multiple entries
    
    Nucleotide label;
    Orientation next_orientation;
    
    dBNode * lst_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(nodes_reverse[length_reverse-1],orientations_reverse[length_reverse-1],
										      &next_orientation, labels_reverse[length_reverse-1],&label,db_graph, get_colour);
    
    //sanity check
    if (lst_node != nodes_reverse[length_reverse]){
      puts("db_graph_supernode broken!\n");
      exit(1);
    }
    
    
    length = db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(nodes_reverse[length_reverse],
											      opposite_orientation(orientations_reverse[length_reverse]),
											      limit,label,
											      node_action,
											      path_nodes,path_orientations,path_labels,
											      supernode_str,avg_coverage,min,max,
											      is_cycle,db_graph, get_colour, get_covg);
    
  }
  else{
    length = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(node,forward,
									      limit,
									      node_action,
									      path_nodes,path_orientations,path_labels,
									      supernode_str,avg_coverage,min,max,
									      is_cycle,db_graph, get_colour, get_covg);
  }
  
  
  

  //apply action to the fst and last node
  node_action(path_nodes[0]);
  node_action(path_nodes[length]);
  
  return length;
}





boolean condition_always_true(dBNode** flank_5p, int len5p, dBNode** ref_branch, int len_ref, dBNode** var_branch, int len_var,
			      dBNode** flank_3p, int len3p, int colour_of_ref, int colour_of_indiv)
{
  return true;
}


//this is the default condition used by the trsuted-path/reference-based SV caller
// it requires of a potential variant site that (after normalising by factoring out coverage due to the rest of the genome)
// there is coverage of variant nodes that are not in the ref-allele/trusted-path branch
// this assumes that we are doing a DIPLOID assembly, and requires you to pass in the expected coverage per chromosome/haplotype
// i.e if you (total reads*read_length)/length_of_genome = 30x coverage, then pass in 15, as you expect 15x per copy of a chromosome.
boolean condition_default(dBNode** flank_5p, int len5p, dBNode** ref_branch, int len_ref, dBNode** var_branch, int len_var,
			      dBNode** flank_3p, int len3p, int colour_of_ref, int colour_of_indiv)
{
  
  //basic idea is simple. Take the nodes which are in the variant and not the trusted path
  
  return true;

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


int db_graph_supernode_returning_query_node_posn_in_subgraph_defined_by_func_of_colours(dBNode * node,int limit,void (*node_action)(dBNode * node), 
											dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
											char * supernode_str, double * avg_coverage,int * min,int * max, boolean * is_cycle, 
											int* query_node_posn,
											dBGraph * db_graph, 
											Edges (*get_colour)(const dBNode*),
											int (*get_covg)(const dBNode*) )
{
  
  
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


  length_reverse = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(node,reverse,limit,&db_node_action_do_nothing,
										    nodes_reverse,orientations_reverse,labels_reverse,
										    supernode_str,&avg_coverager,&minr,&maxr,
										    &is_cycler,db_graph, get_colour, get_covg);

  *query_node_posn = length_reverse;
  
  
  if (length_reverse>0){
    //let's re do the last step, we need to do that because the last node could have had multiple entries
    
    Nucleotide label;
    Orientation next_orientation;

    dBNode * lst_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(nodes_reverse[length_reverse-1],orientations_reverse[length_reverse-1],
										      &next_orientation, labels_reverse[length_reverse-1],&label,db_graph, get_colour);
    
    //sanity check
    if (lst_node != nodes_reverse[length_reverse]){
      puts("db_graph_supernode broken in db_graph_supernode_returning_query_node_posn_in_subgraph_defined_by_func_of_colours\n");
      exit(1);
    }
    

    length = db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(nodes_reverse[length_reverse],
											      opposite_orientation(orientations_reverse[length_reverse]),
											      limit,label,
											      node_action,
											      path_nodes,path_orientations,path_labels,
											      supernode_str,avg_coverage,min,max,
											      is_cycle,db_graph, 
											      get_colour, get_covg);
    
  }
  else{
    length = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(node,forward,
									      limit,
									      node_action,
									      path_nodes,path_orientations,path_labels,
									      supernode_str,avg_coverage,min,max,
									      is_cycle,db_graph, get_colour, get_covg);
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
									   dBGraph * db_graph, EdgeArrayType type, int index)
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

       
  *path_length = db_graph_supernode_for_specific_person_or_pop(node, limit, 
							       action_check_condition_on_all_nodes, 
							       path_nodes, path_orientations, path_labels,
							       string, avg_coverage, min, max, is_cycle,
							       db_graph, type, index);

  
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
							      dBGraph * db_graph, EdgeArrayType type, int index)
{

  return db_graph_logical_operation_on_condition_for_all_nodes_in_supernode(node,limit,
									    condition_for_all_nodes,
									    node_action,
									    path_nodes, path_orientations, path_labels, path_length, 
									    string, avg_coverage, min, max, is_cycle,
									    true, db_graph, type, index);
									   
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
								      dBGraph * db_graph, EdgeArrayType type, int index)
{
  return db_graph_logical_operation_on_condition_for_all_nodes_in_supernode(node,limit,
									    condition_for_all_nodes,
									    node_action,
									    path_nodes, path_orientations, path_labels, path_length,
									    string, avg_coverage, min, max, is_cycle,
									    false, db_graph, type, index);

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
										    dBGraph * db_graph, EdgeArrayType type, int index)
{

  boolean condition_is_true_for_all_nodes_in_supernode=true;

  void action_check_condition_on_all_nodes(dBNode* n)
    {
      //printf("db_graph_is_condition_true..:  applying action to node %s\n", binary_kmer_to_seq(element_get_kmer(n),db_graph->kmer_size,tmp_seq) );
      condition_is_true_for_all_nodes_in_supernode = condition_is_true_for_all_nodes_in_supernode  &&   condition_for_all_nodes(n);
    }

  
  *path_length = db_graph_supernode_for_specific_person_or_pop(node, limit, 
							       action_check_condition_on_all_nodes, 
							       path_nodes, path_orientations, path_labels,
							       string, avg_coverage, min, max, is_cycle,
							       db_graph, type, index);

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
										  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index, EdgeArrayType type, int index)
{

  int count_nodes=0;
  void print_supernode(dBNode * e)
    {
      // TODO  - change this magic 5000 to global max-expected supernode size
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
								  db_graph, type, index))
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



void db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(dBGraph * db_graph, boolean (*condition)(dBNode * node), int min_covg_required, FILE* fout, 
											  boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index, 
											  EdgeArrayType type, int index)
{
  
  int count_nodes=0;
  void print_supernode(dBNode * e)
  {
    
    dBNode * nodes_path[5000];
    Orientation orientations_path[5000];
    Nucleotide labels_path[5000];
    char seq[5000+1];
    seq[0]='\0';
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
									  db_graph, type, index))
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
												       boolean is_for_testing, char** for_test_array_of_supernodes, int* for_test_index,
												       EdgeArrayType type, int index)
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
											  db_graph, type, index))
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



boolean detect_vars_condition_always_true(VariantBranchesAndFlanks* var)
{
  return true;
}

boolean detect_vars_condition_branches_not_marked_to_be_ignored(VariantBranchesAndFlanks* var)
{
  int i;
  for (i=0; i<var->len_one_allele; i++)
    {
      if (db_node_check_status((var->one_allele)[i], ignore_this_node)==true)
	{
	  return false;
	}
    }

  for (i=0; i<var->len_other_allele; i++)
    {
      if (db_node_check_status((var->other_allele)[i], ignore_this_node)==true)
	{
	  return false;
	}
    }

  return true;
}





boolean detect_vars_condition_always_false(VariantBranchesAndFlanks* var)
{
  return false;
}


boolean detect_vars_condition_flanks_at_least_3(VariantBranchesAndFlanks* var)
{
  if ((var->len_flank5p >=3) && (var->len_flank3p >=3))
    {
      return true;
    }
  else
    {
      return false;
    }
}


//typical usecase is, colourfunc1 = colour of ref, colourfunc2 = colour of individual ---> calling hom nonref
boolean detect_vars_condition_is_hom_nonref_given_colour_funcs_for_ref_and_indiv(VariantBranchesAndFlanks* var, 
										 int (*get_covg_colourfunc1)(const dBNode*), int (*get_covg_colourfunc2)(const dBNode*) )
{
  int covg_threshold = 1;
  int i;
  int count_how_many_nodes_in_one_allele_have_covg_by_colourfunc2=0;
  int count_how_many_nodes_in_other_allele_have_covg_by_colourfunc2=0;
  int count_how_many_nodes_in_one_allele_have_covg_by_colourfunc1=0;
  int count_how_many_nodes_in_other_allele_have_covg_by_colourfunc1=0;
  
  for (i=0; i< var->len_one_allele; i++)
    {
      if ( get_covg_colourfunc2((var->one_allele)[i]) >= covg_threshold )
	{
	  count_how_many_nodes_in_one_allele_have_covg_by_colourfunc2++;
	}
      if ( get_covg_colourfunc1((var->one_allele)[i]) > 0 )
	{
	  count_how_many_nodes_in_one_allele_have_covg_by_colourfunc1++;
	}
    }
  for (i=0; i< var->len_other_allele; i++)
    {
      if ( get_covg_colourfunc2((var->other_allele)[i]) >= covg_threshold )
	{
	  count_how_many_nodes_in_other_allele_have_covg_by_colourfunc2++;
	}
      if ( get_covg_colourfunc1((var->other_allele)[i]) > 0 )
	{
	  count_how_many_nodes_in_other_allele_have_covg_by_colourfunc1++;
	}
    }


  if (//colourfunc2/individual has branch1 but not branch2 
      (count_how_many_nodes_in_one_allele_have_covg_by_colourfunc2==var->len_one_allele)
      &&
      (count_how_many_nodes_in_other_allele_have_covg_by_colourfunc2<=1)//last node of the two branches is same
      &&
      //colourfunc1/reference has branch2 only
      (count_how_many_nodes_in_one_allele_have_covg_by_colourfunc1<=1)//Mario - do you agree?
      &&
      (count_how_many_nodes_in_other_allele_have_covg_by_colourfunc1==var->len_other_allele)
      )
    {
      return true;
    }
  else if (//individual has branch2 but not branch1
	   (count_how_many_nodes_in_one_allele_have_covg_by_colourfunc2<=1)
	   &&
	   (count_how_many_nodes_in_other_allele_have_covg_by_colourfunc2==var->len_other_allele)
	   &&
	   //reference has branch1 only
	   (count_how_many_nodes_in_one_allele_have_covg_by_colourfunc1==var->len_one_allele)
	   &&
	   (count_how_many_nodes_in_other_allele_have_covg_by_colourfunc1<=1)
	   )
    {
      //printf("Return true\n");
      return true;
    }
  else
    {
      //printf("Return false\n");
      return false;
    }
  
}






boolean detect_vars_condition_is_hom_nonref_with_min_covg_given_colour_funcs_for_ref_and_indiv(VariantBranchesAndFlanks* var, int (*get_covg_ref)(const dBNode*), int (*get_covg_indiv)(const dBNode*), int min_covg )
{
  //Assumes the reference is colour 0 and the individual is colour 1
  int covg_threshold = min_covg;
  int i;
  int count_how_many_nodes_in_one_allele_have_covg_by_indiv=0;
  int count_how_many_nodes_in_other_allele_have_covg_by_indiv=0;
  int count_how_many_nodes_in_one_allele_have_covg_by_ref=0;
  int count_how_many_nodes_in_other_allele_have_covg_by_ref=0;
  
  for (i=0; i< var->len_one_allele; i++)
    {
      if ( get_covg_indiv((var->one_allele)[i]) >= covg_threshold )
	{
	  count_how_many_nodes_in_one_allele_have_covg_by_indiv++;
	}
      if ( get_covg_ref((var->one_allele)[i]) > 0 )
	{
	  count_how_many_nodes_in_one_allele_have_covg_by_ref++;
	}
    }
  for (i=0; i< var->len_other_allele; i++)
    {
      if ( get_covg_indiv((var->other_allele)[i]) >= covg_threshold )
	{
	  count_how_many_nodes_in_other_allele_have_covg_by_indiv++;
	}
      if ( get_covg_ref((var->other_allele)[i]) > 0 )
	{
	  count_how_many_nodes_in_other_allele_have_covg_by_ref++;
	}
    }


  if (//individual has branch1 but not branch2 
      (count_how_many_nodes_in_one_allele_have_covg_by_indiv==var->len_one_allele)
      &&
      (count_how_many_nodes_in_other_allele_have_covg_by_indiv<=1)//last node of the two branches is same
      &&
      //reference has branch2 only
      (count_how_many_nodes_in_one_allele_have_covg_by_ref<=1)//Mario - do you agree?
      &&
      (count_how_many_nodes_in_other_allele_have_covg_by_ref==var->len_other_allele)
      )
    {
      return true;
    }
  else if (//individual has branch2 but not branch1
	   (count_how_many_nodes_in_one_allele_have_covg_by_indiv<=1)
	   &&
	   (count_how_many_nodes_in_other_allele_have_covg_by_indiv==var->len_other_allele)
	   &&
	   //reference has branch1 only
	   (count_how_many_nodes_in_one_allele_have_covg_by_ref==var->len_one_allele)
	   &&
	   (count_how_many_nodes_in_other_allele_have_covg_by_ref<=1)
	   )
    {
      //printf("Return true\n");
      return true;
    }
  else
    {
      //printf("Return false\n");
      return false;
    }
  
}





//routine to DETECT/DISCOVER variants directly from the graph - reference-free (unless you have put the reference in the graph!)
// "condition" argument is a condition which you apply to the flanks and branches to decide whether to call.
// e.g. this might be some constraint on the coverage of the branches, or one might have a condition that one branch
//      was in one colour  and the other in a different colour, or maybe that both branches are in the same colour
// last argument get_colour specifies some combination of colours, defining the graph within which we look for bubbles.
// most obvious choices are: colour/edge with index (say)2, or union of all edges, or union of ll except that which is the reference genome

void db_graph_detect_vars(FILE* fout, int max_length, dBGraph * db_graph, 
			  boolean (*condition)(VariantBranchesAndFlanks*),
			  void (*action_branches)(dBNode*),
			  void (*action_flanks)(dBNode*),
			  Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*),
			  void (*print_extra_info)(VariantBranchesAndFlanks*, FILE*))
{

  int count_vars = 0; 
  int flanking_length = 1000; 

  //allocate arrays
  dBNode** path_nodes1 = (dBNode**) malloc(sizeof(dBNode*)*(max_length+1));
  dBNode** path_nodes2 = (dBNode**) malloc(sizeof(dBNode*)*(max_length+1));
  Orientation* path_orientations1 = (Orientation*) malloc(sizeof(Orientation)*(max_length+1));
  Orientation* path_orientations2 = (Orientation*) malloc(sizeof(Orientation)*(max_length+1));
  Nucleotide* path_labels1 = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_length+1) );
  Nucleotide* path_labels2 = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_length+1) );
  char* seq1 = (char*) malloc(sizeof(char)*(max_length+1));
  char* seq2 = (char*) malloc(sizeof(char)*(max_length+1));

  if ( (path_nodes1==NULL) || (path_nodes2==NULL) || (path_orientations1==NULL) || (path_orientations2==NULL) 
       || (path_labels1==NULL) || (path_labels2==NULL) || (seq1==NULL) || (seq2==NULL) )
    {
      printf("Could not allocate arrays in db_graph_detect_vars. Out of memory, or you asked for unreasonably big max branch size: %d\n", max_length);
      exit(1);
    }
  
  void get_vars(dBNode * node){
   
    void get_vars_with_orientation(dBNode * node, Orientation orientation){
      
      int length1, length2;
      //dBNode * path_nodes1[max_length+1];
      //dBNode * path_nodes2[max_length+1];
      //Orientation path_orientations1[max_length+1];
      //Orientation path_orientations2[max_length+1];
      //Nucleotide path_labels1[max_length];
      //Nucleotide path_labels2[max_length];
      //char seq1[max_length+1];
      //char seq2[max_length+1];

      double avg_coverage1;
      int min_coverage1,max_coverage1;
      
      double avg_coverage2;
      int min_coverage2,max_coverage2;
      
      dBNode * current_node = node;
      VariantBranchesAndFlanks var;//will reuse this as we traverse the graph
   
      do{
	

	//The idea is that db_graph_detect_bubble_in_subgraph_defined_by_func_of_colours will mark anything is sees with action_flanks
	//However if it does find a bubble, and if you want it to (penultimate argument=true)  then the branches are  marked with action_branches
	if (db_graph_detect_bubble_in_subgraph_defined_by_func_of_colours(current_node,orientation,max_length,action_flanks,
									  &length1,path_nodes1,path_orientations1,path_labels1,
									  seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
									  &length2,path_nodes2,path_orientations2,path_labels2,
									  seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
									  db_graph, get_colour, get_covg,
									  true, action_branches))
	  {

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
	    
	    
	    //compute 5' flanking region       	
	    int length_flank5p_reverse = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(current_node,opposite_orientation(orientation),
													  flanking_length,action_flanks,
													  nodes5p,orientations5p,labels_flank5p,
													  seq5p,&avg_coverage5p,&min5p,&max5p,
													  &is_cycle5p,db_graph, get_colour, get_covg);
	    
	    if (length_flank5p_reverse>0){
	      Nucleotide label;
	      Orientation next_orientation;
	     
	      dBNode * lst_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(nodes5p[length_flank5p_reverse-1],orientations5p[length_flank5p_reverse-1],
												&next_orientation, labels_flank5p[length_flank5p_reverse-1],&label,
												db_graph, get_colour);
	      
	      length_flank5p = db_graph_get_perfect_path_with_first_edge_in_subgraph_defined_by_func_of_colours(nodes5p[length_flank5p_reverse],
														opposite_orientation(orientations5p[length_flank5p_reverse]),
														flanking_length,label,
														action_flanks,
														nodes5p,orientations5p,labels_flank5p,
														seq5p,&avg_coverage5p,&min5p,&max5p,
														&is_cycle5p,db_graph, get_colour, get_covg);
	    }
	    else{
	      length_flank5p = 0;
	    }
	    
	    
	    
	    //compute 3' flanking region
	    length_flank3p = db_graph_get_perfect_path_in_subgraph_defined_by_func_of_colours(path_nodes2[length2],path_orientations2[length2],
											      flanking_length, action_flanks,
											      nodes3p,orientations3p,labels_flank3p,
											      seq3p,&avg_coverage3p,&min3p,&max3p,
											      &is_cycle3p,db_graph, get_colour, get_covg);


	    
	    char name[100];
	    
	    //warning - array of 5prime nodes, oprientations is in reverse order to what you would expect - it is never used in what follows
	    set_variant_branches_and_flanks(&var, nodes5p, orientations5p, length_flank5p, path_nodes1, path_orientations1, length1, 
					    path_nodes2, path_orientations2, length2, nodes3p, orientations3p, length_flank3p, unknown);
	    

	    if (condition(&var)==true) 
	      {
		//printf("\nPassed condition - found VARIATION: %i\n",count_vars);
		count_vars++;
		
		//printf("length 5p flank: %i avg_coverage:%5.2f \n",length_flank5p,avg_coverage5p);	    
		//printf("length branch 1: %i - avg_coverage: %5.2f\n",length1,avg_coverage1);
		//printf("length branch 2: %i - avg_coverage: %5.2f\n",length2,avg_coverage2);
		//printf("length 3p flank: %i avg_coverage:%5.2f\n",length_flank3p,avg_coverage3p);


		//print flank5p - 
		sprintf(name,"var_%i_5p_flank",count_vars);

		print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout,name,length_flank5p,avg_coverage5p,min5p,max5p,
										     nodes5p[0],orientations5p[0],			
										     nodes5p[length_flank5p],orientations5p[length_flank5p],				
										     seq5p,
										     db_graph->kmer_size,true, get_colour, get_covg);	
		
		//print branches
		sprintf(name,"branch_%i_1",count_vars);
		print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout,name,length1,
										     avg_coverage1,min_coverage1,max_coverage1,
										     path_nodes1[0],path_orientations1[0],path_nodes1[length1],path_orientations1[length1],
										     seq1,
										     db_graph->kmer_size,false, get_colour, get_covg);
		
		sprintf(name,"branch_%i_2",count_vars);
		print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout,name,length2,
										     avg_coverage2,min_coverage2,max_coverage2,
										     path_nodes2[0],path_orientations2[0],path_nodes2[length2],path_orientations2[length2],
										     seq2,
										     db_graph->kmer_size,false, get_colour, get_covg);
		
		//print flank3p
		sprintf(name,"var_%i_3p_flank",count_vars);
		print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout,name,length_flank3p,avg_coverage3p,min3p,max3p,
										     nodes3p[0],orientations3p[0],nodes3p[length_flank3p],orientations3p[length_flank3p],
										     seq3p,
										     db_graph->kmer_size,false, get_colour, get_covg);
		//ZAMZAM TEMPORARY FOR SIMS sprintf(name,"variant_%i", count_vars);
		//ZAMZAM TEMPORARY FOR SIM PARSING OF OUTPUT fprintf(fout,"\n%s - extra information\n", name);
		print_extra_info(&var, fout);

		
	      }

	    //db_node_action_set_status_visited(path_nodes2[length2]);
	    action_branches(path_nodes2[length2]);

	    current_node = path_nodes2[length2];
	    orientation = path_orientations2[length2];
	  }
      } while (current_node != node //to avoid cycles
	       && db_node_check_status_none(current_node));
    }
    
    
    if (db_node_check_status_none(node)){     
      //db_node_action_set_status_visited(node);
      action_flanks(node);
      get_vars_with_orientation(node,forward);
      get_vars_with_orientation(node,reverse);
    }
  }
 
  hash_table_traverse(&get_vars,db_graph); 

  //cleanup
  free(path_nodes1);
  free(path_nodes2);
  free(path_orientations1);
  free(path_orientations2);
  free(path_labels1);
  free(path_labels2);
  free(seq1);
  free(seq2);

    
}


//looks for bubbles in the individual's colour, after having marked all bubbles found purely in the reference colour (ie ignore those)
//condition is only applied to the second set of bubble,s called in the individual
//Leaves the bubbles in the ref marked as to be ignored at the end
void db_graph_detect_vars_after_marking_vars_in_reference_to_be_ignored(FILE* fout, int max_length, dBGraph * db_graph, 
									boolean (*condition)(VariantBranchesAndFlanks*),
									Edges (*get_colour_ref)(const dBNode*), int (*get_covg_ref)(const dBNode*) ,
									Edges (*get_colour_indiv)(const dBNode*), int (*get_covg_indiv)(const dBNode*),
									void (*print_extra_info)(VariantBranchesAndFlanks*, FILE*) )
{

  //first detect bubbles in the ref colour, but do not print them out, so they get marked as visited
  db_graph_detect_vars(NULL, max_length, db_graph, 
		       &detect_vars_condition_always_false,
		       &db_node_action_set_status_ignore_this_node,//mark branches to be ignored
		       &db_node_action_set_status_visited, //mark everything else as visited
		       get_colour_ref, get_covg_ref, print_extra_info);

  //unset the nodes marked as visited, but not those marked as to be ignored
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	

  //then start again, detecting variants in the individual.
  printf("Now see if can detect anything afer marking stuff to be ignored\n");
  db_graph_detect_vars(fout, max_length, db_graph, 
		       &detect_vars_condition_branches_not_marked_to_be_ignored,//ignore anything you find that is marked to be ignored
		       &db_node_action_set_status_visited,
		       &db_node_action_set_status_visited,
		       get_colour_indiv, get_covg_indiv, print_extra_info);


  //unset the nodes marked as visited, but not those marked as to be ignored
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	

}




//utility function, do not export
boolean are_two_lists_identical(const int const * list1, int len_list1, const int const * list2, int len_list2)
{
  if (len_list1 != len_list2)
    {
      return false;
    }
  //sort
  
  int copy_list1[len_list1];
  int copy_list2[len_list2];
  int i;
  for (i=0; i< len_list1; i++)
    {
      copy_list1[i]=list1[i];
    }
  for (i=0; i< len_list2; i++)
    {
      copy_list2[i]=list2[i];
    }
  
  qsort(copy_list1,len_list1, sizeof(int), int_cmp); 
  qsort(copy_list2,len_list2, sizeof(int), int_cmp); 
  
  boolean ret=true;

  for (i=0; i<len_list1; i++)
    {
      if (copy_list1[i] != copy_list2[i])
	{
	  ret = false;
	}
    }
  
  return ret;
}

//given two lists of colours, we want to call variants where one branch is entirely in the 
// union of the colours of the first list, but not in the second, and vice-versa for the second branch
//UNLESS both lists are identical, in which case we just call bubbles in the union
//A consequence is that if you call bubbles on 1,2/1,3 you'll get nothing.
//if exclude_ref_bubbles_first==false, just pass in NULL, NULL for the last two arguments
void db_graph_detect_vars_given_lists_of_colours(FILE* fout, int max_length, dBGraph * db_graph, 
						 int* first_list, int len_first_list,
						 int* second_list,  int len_second_list,
						 boolean (*extra_condition)(VariantBranchesAndFlanks* var),
						 void (*print_extra_info)(VariantBranchesAndFlanks*, FILE*),
						 boolean exclude_ref_bubbles_first, 
						 Edges (*get_colour_ref)(const dBNode*), int (*get_covg_ref)(const dBNode*))
{


  Edges get_union_first_list_colours(const dBNode* e)
  {
    int i;
    Edges edges=0;
    
    for (i=0; i< len_first_list; i++)
      {
	edges |= e->individual_edges[first_list[i]];
      }
    return edges;    
  }

  Edges get_union_second_list_colours(const dBNode* e)
  {
    int i;
    Edges edges=0;
    
    for (i=0; i< len_second_list; i++)
      {
	edges |= e->individual_edges[second_list[i]];
      }
    return edges;    
  }

  Edges get_union_first_and_second_list_colours(const dBNode* e)
  {
    int i;
    Edges edges=0;
    //in thise case, do not have to worry about intersection of lists, as is just OR-ing
    for (i=0; i< len_first_list; i++)
      {
	edges |= e->individual_edges[first_list[i]];
      }    
    for (i=0; i< len_second_list; i++)
      {
	edges |= e->individual_edges[second_list[i]];
      }
    return edges;    
  }

  

  int get_covg_of_union_first_list_colours(const dBNode* e)
  {
    int i;
    int covg=0;
  
    for (i=0; i< len_first_list; i++)
      {
	covg += e->coverage[first_list[i]];
      }
    return covg;
  }
  int get_covg_of_union_second_list_colours(const dBNode* e)
  {
    int i;
    int covg=0;
  
    for (i=0; i< len_second_list; i++)
      {
	covg += e->coverage[second_list[i]];
      }
    return covg;
  }

  int get_covg_of_union_first_and_second_list_colours(const dBNode* e)
  {
    int i;
    int covg=0;

    /*
    printf("Got list 1:\n");
    for (i=0; i< len_first_list; i++)
      {
	printf("%d ", first_list[i]);
      }
    printf("\n");
    printf("Got list 2:\n");
    for (i=0; i< len_second_list; i++)
      {
	printf("%d ", second_list[i]);
      }
    printf("\n");
    */

    //concatenate the two lists
    int full_list[len_first_list+len_second_list];
    for (i=0; i< len_first_list; i++)
      {
	full_list[i]=first_list[i];
      }
    for (i=0; i< len_second_list; i++)
      {
	full_list[len_first_list+i]=second_list[i];
      }
    //sort 
    qsort(full_list,len_first_list+len_second_list, sizeof(int), int_cmp); 

    for (i=0; i< len_first_list+len_second_list; i++)
      {
	if ( (i>0) && (full_list[i] != full_list[i-1]) )
	  {
	    //printf("Sorted, merged lists 1 and 2 contains : %d\n", i);
	    covg += e->coverage[full_list[i]];
	  }
	else if (i==0)
	  {
	    //printf("Sorted, merged lists 1 and 2 contains : %d\n", i);
	    covg += e->coverage[full_list[i]];
	  }
      }
    return covg;
  }

  boolean condition_two_branches_lie_in_opposite_lists(VariantBranchesAndFlanks* var)
  {
    return detect_vars_condition_is_hom_nonref_given_colour_funcs_for_ref_and_indiv(var, &get_covg_of_union_first_list_colours, &get_covg_of_union_second_list_colours);
  }

  boolean both_conditions(VariantBranchesAndFlanks* var)
  {
    if ( (condition_two_branches_lie_in_opposite_lists(var)==true) && (extra_condition(var)==true) )
      {
	return true;
      }
    else
      {
	return false;
      }
  }

  boolean both_conditions_and_avoid_branches_marked_ignore(VariantBranchesAndFlanks* var)
  {
    if ( (condition_two_branches_lie_in_opposite_lists(var)==true) && (extra_condition(var)==true) 
	 && (detect_vars_condition_branches_not_marked_to_be_ignored(var)==true) )
      {
	return true;
      }
    else
      {
	return false;
      }
  }

  boolean extra_condition_and_avoid_branches_marked_ignore(VariantBranchesAndFlanks* var)
  {
    if ( (extra_condition(var)==true) && (detect_vars_condition_branches_not_marked_to_be_ignored(var)==true) )
      {
        return true;
      }
    else
      {
        return false;
      }
  }



  if (exclude_ref_bubbles_first==true)
    {

      //first detect bubbles in the ref colour, but do not print them out, so they get marked as visited
      db_graph_detect_vars(NULL, max_length, db_graph, 
			   &detect_vars_condition_always_false,
			   &db_node_action_set_status_ignore_this_node,//mark branches to be ignored
			   &db_node_action_set_status_visited, //mark everything else as visited
			   get_colour_ref, get_covg_ref, print_extra_info);
      
      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
    }


  if (are_two_lists_identical(first_list, len_first_list, second_list, len_second_list)==false)
    {
      
      //then we do demand that a bubble must have one branch in first list and NOT second list, and vice versa.
      if (exclude_ref_bubbles_first==false)
	{
	  db_graph_detect_vars(fout, max_length, db_graph, 
			       &both_conditions,
			       &db_node_action_set_status_visited,
			       &db_node_action_set_status_visited,
			       &get_union_first_and_second_list_colours, &get_covg_of_union_first_and_second_list_colours,
			       print_extra_info);
	}
      else
	{
	  db_graph_detect_vars(fout, max_length, db_graph, 
			       &both_conditions_and_avoid_branches_marked_ignore,
			       &db_node_action_set_status_visited,
			       &db_node_action_set_status_visited,
			       &get_union_first_and_second_list_colours, &get_covg_of_union_first_and_second_list_colours,
			       print_extra_info);
	}
    }
  else//the two lists of colours are identical, so we don't demand the branches take different colours!
    {
      if (exclude_ref_bubbles_first==false)
	{
	  db_graph_detect_vars(fout, max_length, db_graph, 
			       extra_condition,
			       &db_node_action_set_status_visited,
			       &db_node_action_set_status_visited,
			       &get_union_first_and_second_list_colours, &get_covg_of_union_first_and_second_list_colours,
			       print_extra_info);
	}
      else
	{
	  db_graph_detect_vars(fout, max_length, db_graph, 
			       &extra_condition_and_avoid_branches_marked_ignore,
			       &db_node_action_set_status_visited,
			       &db_node_action_set_status_visited,
			       &get_union_first_and_second_list_colours, &get_covg_of_union_first_and_second_list_colours,
			       print_extra_info);
	}
    }



  if (exclude_ref_bubbles_first==true)
    {
      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
    }
}


/*
//ie exclude bubbles found in the reference first, THEN find bubbles where the two branches lie in opposites lists
void db_graph_detect_vars_given_lists_of_colours_excluding_reference_bubbles(FILE* fout, int max_length, dBGraph * db_graph, 
									     int* first_list, int len_first_list,
									     int* second_list,  int len_second_list,
									     Edges (*get_colour_ref)(const dBNode*), int (*get_covg_ref)(const dBNode*) ,
									     boolean (*extra_condition)(VariantBranchesAndFlanks* var),
									     void (*print_extra_info)(VariantBranchesAndFlanks*, FILE*))
{







  //then start again, detecting variants.
  if (are_two_lists_identical(first_list, len_first_list, second_list, len_second_list)==false)
    {

      //then we do demand that a bubble must have one branch in first list and NOT second list, and vice versa. 
      db_graph_detect_vars_given_lists_of_colours(fout, max_length, db_graph, 
						  first_list, len_first_list,
						  second_list, len_second_list,
						  &detect_vars_condition_branches_not_marked_to_be_ignored,//ignore anything you find that is marked to be ignored
						  print_extra_info);
    }

  
}

*/








//do not export - discuss with Mario
void db_graph_smooth_bubbles_for_specific_person_or_pop(int coverage,int limit,dBGraph * db_graph, EdgeArrayType type, int index){
  
  void smooth_bubble(dBNode * node){
    if (db_node_check_status_none(node)){
      db_graph_db_node_smooth_bubble_for_specific_person_or_pop(node,forward,limit,coverage,
								&db_node_action_set_status_pruned,db_graph, type, index);
      db_graph_db_node_smooth_bubble_for_specific_person_or_pop(node,reverse,limit,coverage,
								&db_node_action_set_status_pruned,db_graph, type, index);
      
    }
  }
  
  hash_table_traverse(&smooth_bubble,db_graph); 

}


void db_graph_clip_tips_for_specific_person_or_pop(dBGraph * db_graph, EdgeArrayType type, int index)
{
  
  void clip_tips(dBNode * node){
    
    //use max length k+1, which is what you would get with a single base error - a bubble of that length
    if (db_node_check_status_none(node)){
      db_graph_db_node_clip_tip_for_specific_person_or_pop(node, 1+db_graph->kmer_size,&db_node_action_set_status_pruned,db_graph, type, index);
    }
  }

  hash_table_traverse(&clip_tips,db_graph);
  
}





void db_graph_clip_tips_in_subgraph_defined_by_func_of_colours(dBGraph * db_graph,
							       Edges (*get_colour)(const dBNode*),
							       void (*apply_reset_to_specific_edge_in_colour)(dBNode*, Orientation, Nucleotide),
							       void (*apply_reset_to_colour)(dBNode*))
{
  
  void clip_tips(dBNode * node){
    
    //use max length k+1, which is what you would get with a single base error - a bubble of that length
    if (db_node_check_status_none(node)){
      db_graph_db_node_clip_tip_in_subgraph_defined_by_func_of_colours(node, 1+db_graph->kmer_size,&db_node_action_set_status_pruned,db_graph, 
								       get_colour, apply_reset_to_specific_edge_in_colour, apply_reset_to_colour);
    }
  }

  hash_table_traverse(&clip_tips,db_graph);
  
}


void apply_reset_to_specific_edge_in_union_of_all_colours(dBNode* node, Orientation or, Nucleotide nuc)
{
  int j;
  for (j=0; j<NUMBER_OF_COLOURS; j++)
    {
      reset_one_edge(node, or, nuc, individual_edge_array, j);
    }
}

void apply_reset_to_all_edges_in_union_of_all_colours(dBNode* node )
{
  int j;
  for (j=0; j<NUMBER_OF_COLOURS; j++)
    {
      db_node_reset_edges(node, individual_edge_array, j);
    }
}


void db_graph_clip_tips_in_union_of_all_colours(dBGraph* db_graph)
{
  db_graph_clip_tips_in_subgraph_defined_by_func_of_colours(db_graph,
							    &element_get_colour_union_of_all_colours,
							    &apply_reset_to_specific_edge_in_union_of_all_colours,
							    &apply_reset_to_all_edges_in_union_of_all_colours);
}



void db_graph_print_supernodes_for_specific_person_or_pop(char * filename_sups, char* filename_sings, int max_length, dBGraph * db_graph, EdgeArrayType type, int index, 
							  void (*print_extra_info)(dBNode**, Orientation*, int, FILE*)){

  FILE * fout1; //file to which we will write all supernodes which are longer than 1 node in fasta format
  fout1= fopen(filename_sups, "w"); 
  FILE * fout2; //file to which we will write all "singleton" supernodes, that are just  1 node, in fasta format
  fout2= fopen(filename_sings, "w"); 
  
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
      int length = db_graph_supernode_for_specific_person_or_pop(node,max_length,&db_node_action_set_status_visited,
								 path_nodes,path_orientations,path_labels,
								 seq,&avg_coverage,&min,&max,&is_cycle,
								 db_graph, type, index);
      
 

      if (length>0){	
	sprintf(name,"node_%i",count_nodes);

	print_minimal_fasta_from_path_for_specific_person_or_pop(fout1,name,length,avg_coverage,min,max,
								 path_nodes[0],path_orientations[0],path_nodes[length],path_orientations[length],seq,
								 db_graph->kmer_size,true, type, index);
	if (length==max_length){
	  printf("contig length equals max length [%i] for node_%i\n",max_length,count_nodes);
	}
	//fprintf(fout1, "extra information:\n");
	print_extra_info(path_nodes, path_orientations, length, fout1);
	count_nodes++;
      }
      else{
	sprintf(name,"node_%qd",count_sing);
	print_minimal_fasta_from_path_for_specific_person_or_pop(fout2,name,length,avg_coverage,min,max,
								 path_nodes[0],path_orientations[0],path_nodes[length],path_orientations[length],seq,
								 db_graph->kmer_size,true, type, index);
	//fprintf(fout2, "extra information:\n");
	print_extra_info(path_nodes, path_orientations, length, fout2);
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
  fclose(fout1);
  fclose(fout2);
}


void db_graph_print_supernodes_defined_by_func_of_colours(char * filename_sups, char* filename_sings, int max_length, 
							  dBGraph * db_graph, Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*),
							  void (*print_extra_info)(dBNode**, Orientation*, int, FILE*)){

  boolean do_we_print_singletons=true;//singletons are supernodes consisting of ONE node.

  FILE * fout1; //file to which we will write all supernodes which are longer than 1 node in fasta format
  fout1= fopen(filename_sups, "w"); 
  if (fout1==NULL)
    {
      printf("Cannot open file %s in db_graph_print_supernodes_defined_by_func_of_colours\n", filename_sups);
      exit(1);
    }

  FILE * fout2=NULL; //file to which we will write all "singleton" supernodes, that are just  1 node, in fasta format
  if ( strcmp(filename_sings, "")==0 )
    {
      printf("Only printing supernodes consisting of >1 node (ie contigs longer than %d bases)\n", db_graph->kmer_size);
      do_we_print_singletons=false;
    }
  else
    {
      fout2= fopen(filename_sings, "w"); 
      if (fout2==NULL)
	{
	  printf("Cannot open file %s in db_graph_print_supernodes_defined_by_func_of_colours\n", filename_sings);
	  exit(1);
	}
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

  void print_supernode(dBNode * node){
    
    count_kmers++;
    char name[100];

    if (db_node_check_status_none(node) == true){
      int length = db_graph_supernode_in_subgraph_defined_by_func_of_colours(node,max_length,&db_node_action_set_status_visited,
									     path_nodes,path_orientations,path_labels,
									     seq,&avg_coverage,&min,&max,&is_cycle,
									     db_graph, get_colour, get_covg);
      
      if (length>0){	
	sprintf(name,"node_%i",count_nodes);

	print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout1,name,length,avg_coverage,min,max,
									     path_nodes[0],path_orientations[0],path_nodes[length],path_orientations[length],seq,
									     db_graph->kmer_size,true, get_colour, get_covg);
	if (length==max_length){
	  printf("contig length equals max length [%i] for node_%i\n",max_length,count_nodes);
	}
	//fprintf(fout1, "extra information:\n");
	print_extra_info(path_nodes, path_orientations, length, fout1);
	count_nodes++;
      }
      else{
	count_sing++;
	if (do_we_print_singletons==true)
	  {
	    sprintf(name,"node_%qd",count_sing);
	    print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(fout2,name,length,avg_coverage,min,max,
										 path_nodes[0],path_orientations[0],path_nodes[length],path_orientations[length],seq,
										 db_graph->kmer_size,true, get_colour, get_covg);
	    //fprintf(fout2, "extra information:\n");
	    print_extra_info(path_nodes, path_orientations, length, fout2);
	  }
	

      }
    
    }
  }
  
  hash_table_traverse(&print_supernode,db_graph); 
  printf("%qd nodes visted [%qd singletons]\n",count_kmers,count_sing);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(seq);
  fclose(fout1);
  if (fout2 != NULL)
    {
      fclose(fout2);
    }
}



void db_graph_print_coverage_for_specific_person_or_pop(dBGraph * db_graph, EdgeArrayType type, int index){
  long long count_kmers=0;

  void print_coverage(dBNode * node){
    char seq[db_graph->kmer_size];

    printf("%qd\t%s\t%i\n",count_kmers,binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq),db_node_get_coverage(node, type, index));
    count_kmers++;;
  }
  hash_table_traverse(&print_coverage,db_graph); 
}


//if the node has covg <= coverage (arg2) and its supernode has length <=kmer+1 AND all the interiro nodes of the supernode have this low covg, then
//prune the whole of the interior of the supernode
//note the argument supernode_len is set to -1 if the passed in node has status != none
// returns true if the supernode is pruned, otherwise false
// if returns false, cannot assume path_nodes etc contain anything useful
boolean db_graph_remove_supernode_containing_this_node_if_looks_like_induced_by_singlebase_error(dBNode* node, int coverage, dBGraph * db_graph, int max_expected_sup_len,
											       int (*sum_of_covgs_in_desired_colours)(const Element *), 
											       Edges (*get_edge_of_interest)(const Element*), 
											       void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
											       void (*apply_reset_to_specified_edges_2)(dBNode*),
											       dBNode** path_nodes, Orientation* path_orientations, Nucleotide* path_labels, 
											       char* supernode_str, int* supernode_len)

{

  boolean is_supernode_pruned=true;

  if (node==NULL)
    {
      printf("Called db_graph_remove_supernode_containing_this_node_if_looks_like_induced_by_singlebase_error with a NULL node. Programming error\n");
      exit(1);
    }

  if (db_node_check_status(node, none)==false)//don't touch stuff that is visited or pruned, or whatever
    {
      //printf("Ignore as status us not none\n\n");
      *supernode_len=-1;//caller can check this
      is_supernode_pruned=false;
      return is_supernode_pruned;
    }


  if ( (sum_of_covgs_in_desired_colours(node)>0) && (sum_of_covgs_in_desired_colours(node)<=coverage) )
      {
	//get the supernode, setting nodes to visited
	double avg_cov;
	int min_cov;
	int max_cov;
	boolean is_cycle;
	
	//length_sup is the number of edges in the supernode
	int length_sup =  db_graph_supernode_in_subgraph_defined_by_func_of_colours(node,max_expected_sup_len,
										    &db_node_action_set_status_visited,
										    path_nodes, path_orientations, path_labels, supernode_str,
										    &avg_cov,&min_cov, &max_cov, &is_cycle,
										    db_graph, 
										    get_edge_of_interest,
										    sum_of_covgs_in_desired_colours);

	*supernode_len=length_sup;

        if (length_sup>db_graph->kmer_size+1)
	  {
	    is_supernode_pruned=false;
	    return is_supernode_pruned;//do nothing. This is too long to look like a supernode induced by a single-base error
	  }
	else if (length_sup <=1)
	  {
	    is_supernode_pruned=false;
	    return is_supernode_pruned;//do nothing. This supernode has no interior, is just 1 or 2 nodes, so cannot prune it
	  }
	else
	  {
	    int i;
	    //to look like an error, must all have actual coverage, caused by an actual errored read, BUT must have low covg, <=threshold
	    boolean interior_nodes_look_like_error=true;
	    for (i=1; (i<=length_sup-1) && (interior_nodes_look_like_error==true); i++)
	      {
		if (sum_of_covgs_in_desired_colours(path_nodes[i])>coverage)
		  {
		    interior_nodes_look_like_error=false;
		  }
	      }

	    if (interior_nodes_look_like_error==true)
	      {
		for (i=1; (i<=length_sup-1); i++)
		  {

		    db_graph_db_node_prune_low_coverage(path_nodes[i],coverage,
							&db_node_action_set_status_pruned,
							db_graph,
							sum_of_covgs_in_desired_colours, get_edge_of_interest, apply_reset_to_specified_edges, apply_reset_to_specified_edges_2
							);
		  }
	      }
	    else//interior nodes do not look like error
	      {
		//don't prune - some interior ode has high coverage
		is_supernode_pruned=false;
	      }
	  }
      }
  else//debug only
    {
      //printf("OK - query node has too much covg to consider removing\n");
      is_supernode_pruned=false;
    }
    return is_supernode_pruned;
}


// traverse graph. At each node, if covg <= arg1, get its supernode. If that supernode length is <= kmer-length, and ALL interior nodes have covg <= arg1 
// then prune the node, and the interior nodes of the supernode.
// returns the number of pruned supernodes
long long db_graph_remove_errors_considering_covg_and_topology(int coverage, dBGraph * db_graph,
							       int (*sum_of_covgs_in_desired_colours)(const Element *), 
							       Edges (*get_edge_of_interest)(const Element*), 
							       void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
							       void (*apply_reset_to_specified_edges_2)(dBNode*),
							       int max_expected_sup)
{

  //int max_expected_sup = 10000;
  dBNode**     path_nodes        = (dBNode**) malloc(sizeof(dBNode*)*max_expected_sup); 
  Orientation* path_orientations = (Orientation*) malloc(sizeof(Orientation)*max_expected_sup); 
  Nucleotide*  path_labels       = (Nucleotide*) malloc(sizeof(Nucleotide)*max_expected_sup);
  char*        supernode_string  = (char*) malloc(sizeof(char)*max_expected_sup+1); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (supernode_string==NULL) )
    {
      printf("Cannot malloc arrays for db_graph_remove_errors_considering_covg_and_topology");
      exit(1);
    }


  long long prune_supernode_if_it_looks_like_is_induced_by_singlebase_errors(dBNode* node)
  {

    int len;
    boolean is_sup_pruned;
    
    is_sup_pruned = db_graph_remove_supernode_containing_this_node_if_looks_like_induced_by_singlebase_error(node, coverage, db_graph, max_expected_sup,
													     sum_of_covgs_in_desired_colours,
													     get_edge_of_interest,
													     apply_reset_to_specified_edges, 
													     apply_reset_to_specified_edges_2,
													     path_nodes, path_orientations, path_labels,supernode_string,&len);
    if (is_sup_pruned==true)
      {
	return 1;
      }
    else
      {
	return 0;
      }

  }
  

  long long number_of_pruned_supernodes  = hash_table_traverse_returning_sum(&prune_supernode_if_it_looks_like_is_induced_by_singlebase_errors, db_graph);
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);
  return number_of_pruned_supernodes;
}

							  

// 1. sum_of_covgs_in_desired_colours returns the sum of the coverages for the colours you are interested in
// 1. the argument get_edge_of_interest is a function that gets the "edge" you are interested in - may be a single edge/colour from the graph, or might be a union of some edges 
// 2 Pass apply_reset_to_specified_edges which applies reset_one_edge to whichever set of edges you care about,
// 3 Pass apply_reset_to_specified_edges_2 which applies db_node_reset_edges to whichever set of edges you care about,
void db_graph_remove_low_coverage_nodes(int coverage, dBGraph * db_graph,
					int (*sum_of_covgs_in_desired_colours)(const Element *), 
					Edges (*get_edge_of_interest)(const Element*), 
					void (*apply_reset_to_specified_edges)(dBNode*, Orientation, Nucleotide), 
					void (*apply_reset_to_specified_edges_2)(dBNode*) )
{
  
  void prune_node(dBNode * node){
    db_graph_db_node_prune_low_coverage(node,coverage,
					&db_node_action_set_status_pruned,
					db_graph,
					sum_of_covgs_in_desired_colours, get_edge_of_interest, apply_reset_to_specified_edges, apply_reset_to_specified_edges_2
					);
  }

  hash_table_traverse(&prune_node,db_graph); 
} 


void db_graph_remove_low_coverage_nodes_ignoring_colours(int coverage, dBGraph * db_graph)
{
  
  void prune_node(dBNode * node){
    db_graph_db_node_prune_low_coverage_ignoring_colours(node,coverage,
							 &db_node_action_set_status_pruned,
							 db_graph
							 );
  }

  hash_table_traverse(&prune_node,db_graph); 
}



//if you don't want to/care about graph_info, pass in NULL
void db_graph_dump_binary(char * filename, boolean (*condition)(dBNode * node), dBGraph * db_graph, GraphInfo* db_graph_info){
  FILE * fout; //binary output
  fout= fopen(filename, "w"); 
  
  if (db_graph_info==NULL)
    {
      int i;
      int means[NUMBER_OF_COLOURS];
      long long tots[NUMBER_OF_COLOURS];
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  means[i]=0;
	  tots[i]=0;
	}
      print_binary_signature(fout, db_graph->kmer_size, NUMBER_OF_COLOURS, means, tots);
    }
  else
    {
      print_binary_signature(fout, db_graph->kmer_size, NUMBER_OF_COLOURS, db_graph_info->mean_read_length, db_graph_info->total_sequence);
    }

  long long count=0;
  //routine to dump graph
  void print_node_multicolour_binary(dBNode * node){   
    if (condition(node)){
      count++;
      db_node_print_multicolour_binary(fout,node);
    }
  }

  hash_table_traverse(&print_node_multicolour_binary,db_graph); 
  fclose(fout);

  printf("%qd kmers dumped\n",count);
}



void db_graph_dump_single_colour_binary_of_colour0(char * filename, boolean (*condition)(dBNode * node), dBGraph * db_graph, GraphInfo* db_graph_info){
  FILE * fout; //binary output
  fout= fopen(filename, "w"); 
  if (fout==NULL)
    {
      printf("Unable to open output file %s\n", filename);
      exit(1);
    }

  if (db_graph_info==NULL)
    {
      int means=0;
      long long tots=0;
      print_binary_signature(fout, db_graph->kmer_size,1, &means, &tots);
    }
  else
    {
      print_binary_signature(fout, db_graph->kmer_size, 1, &(db_graph_info->mean_read_length[0]), &(db_graph_info->total_sequence[0]) );
    }


  long long count=0;
  //routine to dump graph
  void print_node_single_colour_binary_of_colour0(dBNode * node){   
    if (condition(node) ){
      count++;
      db_node_print_single_colour_binary_of_colour0(fout,node);
    }
  }

  hash_table_traverse(&print_node_single_colour_binary_of_colour0,db_graph); 
  fclose(fout);

  //printf("%qd kmers dumped\n",count);
}





void db_graph_dump_single_colour_binary_of_specified_colour(char * filename, boolean (*condition)(dBNode * node), dBGraph * db_graph, GraphInfo* db_graph_info, int colour){

  if ( (colour<0) || (colour>=NUMBER_OF_COLOURS) )
    {
      printf("Cannot call db_graph_dump_single_colour_binary_of_specified_colour and try to dump a colour not in 0....%d when %d is the compile time max number of colours supported\n",
	     NUMBER_OF_COLOURS-1, NUMBER_OF_COLOURS);
      exit(1);
    }

  FILE * fout; //binary output
  fout= fopen(filename, "w"); 
  if (fout==NULL)
    {
      printf("Unable to open output file %s\n", filename);
      exit(1);
    }

  if (db_graph_info==NULL)
    {
      int means=0;
      long long tots=0;

      print_binary_signature(fout, db_graph->kmer_size,1, &means, &tots);
    }
  else
    {
      print_binary_signature(fout, db_graph->kmer_size, 1, &(db_graph_info->mean_read_length[colour]), &(db_graph_info->total_sequence[colour]) );
    }


  long long count=0;
  //routine to dump graph
  void print_node_single_colour_binary_of_specified_colour(dBNode * node){   
    if (condition(node)){
      count++;
      db_node_print_single_colour_binary_of_specified_colour(fout,node, colour);
    }
  }

  hash_table_traverse(&print_node_single_colour_binary_of_specified_colour,db_graph); 
  fclose(fout);

  //printf("%qd kmers dumped\n",count);
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
      //shift along so the 4 most significant bits become the 4 least - we've passed our argument by copy so not altering the original
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

  //ASSUMING read length is strictly greater than kmer-length
  //then you should never see a kmer (node) which is unconnected to another (even if the other is itself)
  //ie to check if this kmer is seen in a given person/pop, it is enough to check if it has an edge

  if (db_node_is_this_node_in_this_person_or_populations_graph(e, type, index)==false)
    {
      return NULL;
    }
  else
    {
      return e;
    }
  
}



void db_graph_traverse_with_array(void (*f)(HashTable*, Element *, int**, int, EdgeArrayType, int),HashTable * hash_table, int** array, int length_of_array, EdgeArrayType type, int index){
  
  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(hash_table, &hash_table->table[i], array, length_of_array, type, index);
    }
  }
  
  
}


void db_graph_traverse_with_array_of_longlongs(void (*f)(HashTable*, Element *, long long**, int, EdgeArrayType, int),
					       HashTable * hash_table, long long** array, int length_of_array, EdgeArrayType type, int index)
{
  
  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(hash_table, &hash_table->table[i], array, length_of_array, type, index);
    }
  }
}

void db_graph_get_covg_distribution(char* filename, dBGraph* db_graph, EdgeArrayType type, int index, boolean (*condition)(dBNode* elem) )
{
  int i;

  FILE* fout=fopen(filename, "w");
  if (fout==NULL)
    {
      printf("Cannot open %s\n", filename);
      exit(1);
    }


  long long* covgs = (long long*) malloc(sizeof(long long) * 10001);
  if (covgs==NULL)
    {
      printf("Could not alloc array to hold covg distrib\n");
      exit(1);
    }
  long long** covgs_ptrs = (long long**) malloc(sizeof(long long*) * 10001);
  if (covgs_ptrs==NULL)
    {
      printf("Could not alloc array to hold covg distrib\n");
      exit(1);
    }

  for (i=0; i<=10000; i++)
    {
      covgs[i]=0;
      covgs_ptrs[i]=&(covgs[i]);
    }

  void bin_covg_and_add_to_array(HashTable* htable, Element * e , long long** arr, int len, EdgeArrayType type, int colour)
  {
    if (condition(e)==true)
      {
	int bin = e->coverage[colour];
	if (bin>10000)
	  {
	    bin = 10000;
	  }
	*(arr[bin]) = *(arr[bin]) +1; 
      }
  }
  
  db_graph_traverse_with_array_of_longlongs(&bin_covg_and_add_to_array, db_graph, covgs_ptrs, 10001, type, index);

  for (i=0; i<=10000; i++)
    {
      fprintf(fout, "multiplicity:%d\tNumber:%qd\n",i,covgs[i]);
    }
  fclose(fout);
  free(covgs_ptrs);
  free(covgs);
}







void db_graph_get_covg_distribution_on_bubbles(dBGraph* db_graph, EdgeArrayType type, int index)
{
  int i;

  long long* covgs = (long long*) malloc(sizeof(long long) * 10001);
  if (covgs==NULL)
    {
      printf("Could not alloc array to hold covg distrib\n");
      exit(1);
    }
  long long** covgs_ptrs = (long long**) malloc(sizeof(long long*) * 10001);
  if (covgs_ptrs==NULL)
    {
      printf("Could not alloc array to hold covg distrib\n");
      exit(1);
    }

  for (i=0; i<=10000; i++)
    {
      covgs[i]=0;
      covgs_ptrs[i]=&(covgs[i]);
    }

  void bin_covg_and_add_to_array(HashTable* htable, Element * e , long long** arr, int len, EdgeArrayType type, int colour)
  {
    int bin = e->coverage[colour];
    if (bin>10000)
      {
	bin = 10000;
      }
    *(arr[bin]) = *(arr[bin]) +1; 
  }
  
  db_graph_traverse_with_array_of_longlongs(&bin_covg_and_add_to_array, db_graph, covgs_ptrs, 10001, type, index);

  for (i=0; i<=10000; i++)
    {
      printf("multiplicity:%d\tNumber:%qd\n",i,covgs[i]);
    }

  free(covgs_ptrs);
  free(covgs);
}




//TODO - update this to use db_graph_supernode
/*
//the length of the supernode is passed to the caller in the array of supernode lengths that is passed in.
//the idea is that this function is passed into the traverse function, and the caller of that
void db_graph_get_supernode_length_marking_it_as_visited(dBGraph* db_graph, Element* node, int** array_of_supernode_lengths, int length_of_array,
							 EdgeArrayType type, int index)
{

  if (db_node_check_status(node, visited) || db_node_check_status(node,pruned))
    {
      return;
    }

  int length_of_supernode=1;
  dBNode* first_node=db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(node, individual_edge_array, 0 , db_graph);
  dBNode* current_node;
  dBNode* next_node;
  current_node=first_node;
  Orientation start_orientation, current_orientation, next_orientation;
  db_node_set_status(current_node, visited);

  //work out which direction to leave supernode in. 
  if (db_node_is_supernode_end(first_node,forward, type, index, db_graph))
    {
      if (db_node_is_supernode_end(first_node,reverse, type, index,db_graph))
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
  while (!db_node_is_supernode_end(current_node,current_orientation, type, index, db_graph))
    {
      next_node = db_graph_get_next_node_in_supernode_for_specific_person_or_pop(current_node, current_orientation, &next_orientation, type, index, db_graph);

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





int db_graph_get_N50_of_supernodes(dBGraph* db_graph, EdgeArrayType type, int index)
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


  db_graph_traverse_with_array(&db_graph_get_supernode_length_marking_it_as_visited, db_graph,  numbers_of_supernodes_of_specific_lengths_ptrs, 10000, type, index);

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
*/




  /* get rid of this

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
  */


/* delete?
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
*/




void db_graph_traverse_to_gather_statistics_about_people(void (*f)(HashTable*, Element *, int**, int),HashTable * hash_table, int** array, int num_people )
{

  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(hash_table, &hash_table->table[i], array, num_people);
    }
  }




}






void db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population(dBGraph* hash_table, EdgeArrayType type, int index)
{
  printf("not implemented db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population yet");
  exit(1);
}




//TODO - get rid of this
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
      printf("db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size, tmp_seq));
      exit(1);
    }	         

    if (DEBUG)
      {
	BinaryKmer tmp_kmer;
	printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide1),
	       next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size, tmp_seq) :  
	       binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size, &tmp_kmer),db_graph->kmer_size, tmp_seq));
	     
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



//TODO - get rid of this

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
      BinaryKmer tmp_kmer;
      printf("TRY TO ADD %c - next node %s\n",binary_nucleotide_to_char(nucleotide_for_only_edge),
	     next_orientation == forward ? binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size, tmp_seq) :  
	     binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(next_node),db_graph->kmer_size, &tmp_kmer),db_graph->kmer_size, tmp_seq));
      
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


  int length_of_best_sub_supernode_in_each_person[NUMBER_OF_COLOURS];
  int index_of_start_of_best_sub_supernode_in_each_person[NUMBER_OF_COLOURS];

  int i;
  for (i=0; i< NUMBER_OF_COLOURS; i++)
    {
      length_of_best_sub_supernode_in_each_person[i]=0;
      index_of_start_of_best_sub_supernode_in_each_person[i]=0;
    }

  for (i=0; i< NUMBER_OF_COLOURS; i++)
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
  for (i=0; (i< NUMBER_OF_COLOURS) && (noone_has_decent_sub_supernode==true); i++)
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
  for (i=0; i< NUMBER_OF_COLOURS; i++)
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
      first_kmer_in_subsection = binary_kmer_to_seq(&(current_node->kmer), db_graph->kmer_size, tmp_seq);
    }
  else
    {
      BinaryKmer tmp_kmer;
      first_kmer_in_subsection = binary_kmer_to_seq(   binary_kmer_reverse_complement(&(current_node->kmer), db_graph->kmer_size, &tmp_kmer) , db_graph->kmer_size, tmp_seq);
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
	      next_kmer = binary_kmer_to_seq(element_get_kmer(next_node),db_graph->kmer_size, tmp_seq);
	    }
	  else
	    {
	      BinaryKmer tmp_kmer;
	      next_kmer=binary_kmer_to_seq( binary_kmer_reverse_complement(&(next_node->kmer),db_graph->kmer_size, &tmp_kmer), db_graph->kmer_size, tmp_seq );
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

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (get_edge_copy(*node, individual_edge_array, i) ==0 )
	{
	}
      else
	{
	  number_of_individuals_with_this_node++;
	}
    }

  char* kmer_as_string = binary_kmer_to_seq(element_get_kmer(node), db_graph->kmer_size, tmp_seq);
  fprintf(list_of_file_ptrs[number_of_individuals_with_this_node], "%s\n", kmer_as_string);
  
}

//array_of_counts totals up the number of kmers that are shared by 0,1,2,... individuals. Obviously the first element should be zero (no element
//in the graph is shared by no people) - use this as a crude check.
//number of people loaded is a param
void find_out_how_many_individuals_share_this_node_and_add_to_statistics(HashTable* db_graph, dBNode * node, int** array_of_counts, int number_of_people)
{

  //  char tmp_seq[db_graph->kmer_size];
  int i;

  if (number_of_people>NUMBER_OF_COLOURS)
    {
      printf("Cannot call find_out_how_many_individuals_share_this_node_and_add_to_statistics with number_of_people = %d, as it's bigger than the NUMBER per pop, %d", number_of_people,NUMBER_OF_COLOURS);
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


// Note what happens as you repeatedly call this, with number_of_nodes_to_load = length_of_arrays/2. The first time, you load nodes into the back (right hand) 
// end of the array, and the front end is empty.
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


int db_node_addr_cmp(const void *a, const void *b)
{
  const dBNode** ia = (const dBNode **)a; // casting pointer types
  const dBNode** ib = (const dBNode **)b;
  return *ia  - *ib;
  //address comparison: returns negative if b > a
  // and positive if a > b 
  
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


void get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(dBNode** array, int length, 
									     int* min_coverage, int* max_coverage, double* avg_coverage, 
									     int* mode_coverage, double*  percent_nodes_having_modal_value,
									     Edges (*get_colour)(const dBNode*),
									     int (*get_covg)(const dBNode*)
									     )
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
	  this_covg = db_node_get_coverage_in_subgraph_defined_by_func_of_colours(array[i], get_covg);
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


//this is horribly inefficient - order n^2 to compare. But for now I prefer correctness/definitely right-ness over efficiency
// also remember that 95% of everything you find will be very short.
// I won't count null pointers that are in one but not the other
void get_covg_of_nodes_in_one_but_not_other_of_two_arrays(dBNode** array1, dBNode** array2, int length1, int length2, 
							  int* num_nodes_in_array_1not2, int * num_nodes_in_array_2not1, int** covgs_in_1not2, int** covgs_in_2not1,
							  dBNode** reused_working_array1, dBNode** reused_working_array2,
							  EdgeArrayType type, int index)
     
{
  int j=0;
  int i=0;
  *num_nodes_in_array_1not2=0;
  *num_nodes_in_array_2not1=0;

  //copy nodes into working arrays, which we then sort.
  for (i=0; i< length1; i++)
    {
      reused_working_array1[i]=array1[i];
    }
  for (i=0; i < length2; i++)
    {
      reused_working_array2[i]=array2[i];
    }
  
  qsort(reused_working_array1, length1, sizeof(dBNode*), db_node_addr_cmp);
  qsort(reused_working_array2, length2, sizeof(dBNode*), db_node_addr_cmp);


  for (i=0; i<length1; i++)
    {
      boolean elem_in_both_arrays=false;
      if (reused_working_array1[i]==NULL)
	{
	  continue;
	}
      if (i>0)
	{
	  if (reused_working_array1[i]==reused_working_array1[i-1])
	    {
	      continue;
	    }
	}
      for (j=0; (j<length2) && (elem_in_both_arrays==false); j++)
	{
	  if (reused_working_array1[i]==reused_working_array2[j])
	    {
	      //ignore this element - it is in both arrays
	      elem_in_both_arrays=true;
	    }
	}
      
      if (elem_in_both_arrays==false)
	{
	  //this element reused_working_array1[i] is in reused_working_array1 but not reused_working_array2, and is not NULL
	  *num_nodes_in_array_1not2= *num_nodes_in_array_1not2+1;
	  *(covgs_in_1not2[*num_nodes_in_array_1not2-1])=db_node_get_coverage(reused_working_array1[i], type, index);	  
	}

    }

  //now the other way around

  for (i=0; i<length2; i++)
    {
      boolean elem_in_both_arrays=false;
      if (reused_working_array2[i]==NULL)
	{
	  continue;
	}
      if (i>0)
	{
	  if (reused_working_array2[i]==reused_working_array2[i-1])
	    {
	      continue;
	    }
	}

      for (j=0; (j<length1) && (elem_in_both_arrays==false); j++)
	{
	  if (reused_working_array2[i]==reused_working_array1[j])
	    {
	      //ignore this element - it is in both arrays
	      elem_in_both_arrays=true;
	    }
	}
      
      if (elem_in_both_arrays==false)
	{
	  //this element reused_working_array2[i] is in reused_working_array2 but not reused_working_array1, and is not NULL
	  *num_nodes_in_array_2not1= *num_nodes_in_array_2not1+1;
	  *(covgs_in_2not1[*num_nodes_in_array_2not1-1])=db_node_get_coverage(reused_working_array2[i], type, index);	  
	}

    }




}


//this is horribly inefficient - order n^2 to compare. But for now I prefer correctness/definitely right-ness over efficiency
// also remember that 95% of everything you find will be very short.
// I won't count null pointers that are in one but not the other
void get_covg_of_nodes_in_one_but_not_other_of_two_arrays_in_subgraph_defined_by_func_of_colours(dBNode** array1, dBNode** array2, int length1, int length2, 
												 int* num_nodes_in_array_1not2, int * num_nodes_in_array_2not1, 
												 int** covgs_in_1not2, int** covgs_in_2not1,
												 dBNode** reused_working_array1, dBNode** reused_working_array2,
												 Edges (*get_colour)(const dBNode*),int (*get_covg)(const dBNode*))
  
{
  int j=0;
  int i=0;
  *num_nodes_in_array_1not2=0;
  *num_nodes_in_array_2not1=0;

  //copy nodes into working arrays, which we then sort.
  for (i=0; i< length1; i++)
    {
      reused_working_array1[i]=array1[i];
    }
  for (i=0; i < length2; i++)
    {
      reused_working_array2[i]=array2[i];
    }
  
  qsort(reused_working_array1, length1, sizeof(dBNode*), db_node_addr_cmp);
  qsort(reused_working_array2, length2, sizeof(dBNode*), db_node_addr_cmp);


  for (i=0; i<length1; i++)
    {
      boolean elem_in_both_arrays=false;
      if (reused_working_array1[i]==NULL)
	{
	  continue;
	}
      if (i>0)
	{
	  if (reused_working_array1[i]==reused_working_array1[i-1])
	    {
	      continue;
	    }
	}
      for (j=0; (j<length2) && (elem_in_both_arrays==false); j++)
	{
	  if (reused_working_array1[i]==reused_working_array2[j])
	    {
	      //ignore this element - it is in both arrays
	      elem_in_both_arrays=true;
	    }
	}
      
      if (elem_in_both_arrays==false)
	{
	  //this element reused_working_array1[i] is in reused_working_array1 but not reused_working_array2, and is not NULL
	  *num_nodes_in_array_1not2= *num_nodes_in_array_1not2+1;
	  *(covgs_in_1not2[*num_nodes_in_array_1not2-1])=db_node_get_coverage_in_subgraph_defined_by_func_of_colours(reused_working_array1[i], get_covg);	  
	}

    }

  //now the other way around

  for (i=0; i<length2; i++)
    {
      boolean elem_in_both_arrays=false;
      if (reused_working_array2[i]==NULL)
	{
	  continue;
	}
      if (i>0)
	{
	  if (reused_working_array2[i]==reused_working_array2[i-1])
	    {
	      continue;
	    }
	}

      for (j=0; (j<length1) && (elem_in_both_arrays==false); j++)
	{
	  if (reused_working_array2[i]==reused_working_array1[j])
	    {
	      //ignore this element - it is in both arrays
	      elem_in_both_arrays=true;
	    }
	}
      
      if (elem_in_both_arrays==false)
	{
	  //this element reused_working_array2[i] is in reused_working_array2 but not reused_working_array1, and is not NULL
	  *num_nodes_in_array_2not1= *num_nodes_in_array_2not1+1;
	  *(covgs_in_2not1[*num_nodes_in_array_2not1-1])=db_node_get_coverage_in_subgraph_defined_by_func_of_colours(reused_working_array2[i], get_covg);	  
	}

    }




}



boolean make_reference_path_based_sv_calls_condition_always_true(VariantBranchesAndFlanks* var, int colour_ref, int colour_indiv)
{
  return true;
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
						char** return_flank3p_array, int** return_variant_start_coord,
						boolean (*condition)(VariantBranchesAndFlanks* var,  int colour_of_ref,  int colour_of_indiv),
						void (*action_for_branches_of_called_variants)(VariantBranchesAndFlanks* var),
						void (*print_extra_info)(VariantBranchesAndFlanks* var, FILE* fout)
						)
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
  //kmer_window->kmer = (BinaryKmer*) malloc(sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER*length_of_arrays);    //*(number_of_nodes_to_load + db_graph->kmer_size));
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



  //finally, we need a couple of arrays that will hold dBNode*'s, which we wil reuse when comparing two branches
  dBNode** working_array1 = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 
  dBNode** working_array2 = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 

  int* covgs_in_trusted_not_variant = (int*) malloc(sizeof(int)*length_of_arrays);
  int* covgs_in_variant_not_trusted = (int*) malloc(sizeof(int)*max_expected_size_of_supernode);
  //int* covgs_of_indiv_on_trusted_path = (int*) malloc(sizeof(int)*length_of_arrays);
  //int* covgs_of_ref_on_trusted_path = (int*) malloc(sizeof(int)*length_of_arrays);
  //int* covgs_of_indiv_on_variant_path = (int*) malloc(sizeof(int)*length_of_arrays);
  //int* covgs_of_ref_on_variant_path = (int*) malloc(sizeof(int)*length_of_arrays);

  if ( (working_array1==NULL) || (working_array2==NULL)  ) 
    {
      printf("OOM at start - cannot alloc working arrays. Exit");
      exit(1);
    }

  if ( (covgs_in_trusted_not_variant==NULL) || (covgs_in_variant_not_trusted==NULL) )
    // || (covgs_of_indiv_on_trusted_path==NULL) || (covgs_of_ref_on_trusted_path==NULL)  || (covgs_of_indiv_on_variant_path==NULL) || (covgs_of_ref_on_variant_path==NULL) )
    {
      printf("OOM at start  - cannot malloc covg arrays. Exit\n");
      exit(1);
    }

  int** ptrs_to_covgs_in_trusted_not_variant = (int**) malloc(sizeof(int*)*length_of_arrays);
  int** ptrs_to_covgs_in_variant_not_trusted = (int**) malloc(sizeof(int*)*max_expected_size_of_supernode);

  if ( (ptrs_to_covgs_in_trusted_not_variant==NULL) || (ptrs_to_covgs_in_variant_not_trusted==NULL) )
    {
      printf("Cannot alloc ptr arrays for covgs in make_ref assisted calls fn");
      exit(1);
    }

  for (k=0; k< length_of_arrays; k++)
    {
      ptrs_to_covgs_in_trusted_not_variant[k]=&covgs_in_trusted_not_variant[k];
    }
  for (k=0; k< max_expected_size_of_supernode; k++)
    {
      ptrs_to_covgs_in_variant_not_trusted[k]=&covgs_in_variant_not_trusted[k];
    }


  VariantBranchesAndFlanks var; //we will reuse this


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
  BinaryKmer tmp_kmer;
  


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
	      //printf("ZAM DEBUG node number %d not in the graph of person %d\n", start_node_index, index_for_indiv_in_edge_array);

	      start_node_index++;
	      continue;
	    }

	  if (db_node_check_status_visited(chrom_path_array[start_node_index])==true)
	  {
	    //printf("ZAM DEBUG - node %d is visited\n", start_node_index);

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
	      //printf("ZAM DEBUG - found potential SV locus\n");
	      num_variants_found++;


	      // Note if max_desired_returns>0, then we are going to return the first max_desired_returns  results in the arguments passed in for this purpose, return_branch1_array etc
	      // We can check this easily - if num_variants_found<max_desired_returns, then we add this result to the arrays. This is essentically used only for testing.


	      
	      // this is length of 5p flank in number of edges in the 5p flank.
	      // (in addition there will be one kmer of bases of course. The print function we call will handle this for us)
	      int length_5p_flank = first_index_in_chrom_where_supernode_differs_from_chromosome-start_node_index-1;


	      //work out start coordinate and print it out
	      int start_coord_of_variant_in_trusted_path_fasta = coord_of_start_of_array_in_trusted_fasta + first_index_in_chrom_where_supernode_differs_from_chromosome + db_graph->kmer_size;
	      //fprintf(output_file, "VARIATION: %d, start coordinate %d\n", num_variants_found, start_coord_of_variant_in_trusted_path_fasta);
	      
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

	      //Our algorithm would say: OK, KingZam is the 5prime anchor, let's see if we can find a 3prime anchor that is 7 characters long (say). 
	      // Aha! England will do. OK - can we extend backwards? Yes, we can extend backwards 
	      // all the way back to Zam. But WAIT!! Zam is before the split!! That's fine, it means there is a repeat, just don't try and extend beyond there. 
	      //Hence the else if condition marked below with ***

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
		  else if (  (traverse_sup_left_to_right) && (start_of_3prime_anchor_in_sup - j <= index_in_supernode_where_supernode_differs_from_chromosome+1) ) 
		    // *** see comment above this while loop
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

	      
	      //printf("Collected following info:\n");
	      //printf("start_of_3prime_anchor_in_sup = %d\n", start_of_3prime_anchor_in_sup);
	      //printf("start_of_3prime_anchor_in_chrom = %d\n", start_of_3prime_anchor_in_chrom);
	      //printf("Length 3p flank is  %d\n", length_3p_flank);
	      //printf("Length 5p flank is  %d\n", length_5p_flank);
	      //printf("Length of arrays is %d\n", length_of_arrays);
	      //printf("Start node index is %d\n", start_node_index);
	      //printf("first_index_in_chrom_where_supernode_differs_from_chromosome is %d\n", first_index_in_chrom_where_supernode_differs_from_chromosome);
	      //printf("index_of_query_node_in_supernode_array is %d\n", index_of_query_node_in_supernode_array);
	      
	      
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

	      
	      

	      // 5pflank info now ready for printing. Will be done below, after have collected info on branches and 3p flank


	      //this for testing only (ie used by unit tests)
	      if (num_variants_found<=max_desired_returns)
		{
		  //add first kmer onto flank5p - this si done automatically by print_fasta_from_path
		  if (chrom_orientation_array[start_node_index]==forward)
		    {
		      strcat(return_flank5p_array[num_variants_found-1],
			     binary_kmer_to_seq(element_get_kmer(chrom_path_array[start_node_index]), db_graph->kmer_size, tmp_zam));
		    }
		  else
		    {
		      char tmp_zam2[db_graph->kmer_size+1];
		      tmp_zam2[db_graph->kmer_size]='\0';

		      strcat(return_flank5p_array[num_variants_found-1],
			     seq_reverse_complement(binary_kmer_to_seq(element_get_kmer(chrom_path_array[start_node_index]), db_graph->kmer_size, tmp_zam), 
						    db_graph->kmer_size,tmp_zam2)
			     );
		    }
		      
		  strcat(return_flank5p_array[num_variants_found-1],flank5p);
		}

		    


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
	      int num_nodes_on_trusted_but_not_variant_branch=0;
	      int num_nodes_on_variant_but_not_trusted_branch=0;
	      
	      //get coverage on the trusted branch. 
	      get_coverage_from_array_of_nodes(chrom_path_array+start_node_index+length_5p_flank,
					       len_trusted_branch, &trusted_branch_min_covg, &trusted_branch_max_covg, &trusted_branch_avg_covg, &trusted_branch_mode_coverage, &trusted_branch_percent_nodes_with_modal_covg,
					       which_array_holds_indiv, index_for_indiv_in_edge_array);

	      get_percent_novel_from_array_of_nodes(chrom_path_array+start_node_index+length_5p_flank, 
						    len_trusted_branch, &trusted_branch_percent_novel,
						    which_array_holds_ref, index_for_ref_in_edge_array);
	      



	      //almost ready to print trusted branch, but need to get coverages of nodes that are in trusted but not variant branch, and vice-versa. 
	      // will first set up variant branch, and then print the two together

	      

	      int branch2_min_covg=0;
	      int branch2_max_covg=0;
	      double branch2_avg_covg=0;
	      int branch2_mode_coverage=0;
	      double branch2_percent_nodes_with_modal_covg=0;
	      double branch2_percent_novel=0;

	      int len_branch2;

								 
              // ******************* variant branch - left to right
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

		  
		  //get coverages for variant branch
		  get_coverage_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank,
                                                   len_branch2, &branch2_min_covg, &branch2_max_covg, &branch2_avg_covg, &branch2_mode_coverage, &branch2_percent_nodes_with_modal_covg,
                                                   which_array_holds_indiv, index_for_indiv_in_edge_array);

		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank, 
							len_branch2, &branch2_percent_novel,
							which_array_holds_ref, index_for_ref_in_edge_array);





		  //now do comparison with trusted branch - look for nodes that are in one but not the other
		  // data from this will be used in the next two calls to print_fasta, for the two branches - trusted, and variant
		  get_covg_of_nodes_in_one_but_not_other_of_two_arrays(chrom_path_array+start_node_index+length_5p_flank, current_supernode+index_of_query_node_in_supernode_array+length_5p_flank, 
								       len_trusted_branch, len_branch2,
		    						       &num_nodes_on_trusted_but_not_variant_branch, &num_nodes_on_variant_but_not_trusted_branch, 
								       ptrs_to_covgs_in_trusted_not_variant, ptrs_to_covgs_in_variant_not_trusted,
		    						       working_array1, working_array2,
		    						       which_array_holds_indiv, index_for_indiv_in_edge_array);
		  //we are ready to print the trusted and variant paths, when variant path traversed left to right. Will do it after have collected info on 3p flank
		  


		  //for unit tests only
		  if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_trusted_branch_array[num_variants_found-1],trusted_branch);
		      strcat(return_variant_branch_array[num_variants_found-1],variant_branch);
		    }
		  


		}
	      else //************ variant branch - traversing supernode right to left
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
		      strcat(temp, binary_kmer_to_seq(element_get_kmer(current_supernode[0]), db_graph->kmer_size, tmp_zam));
		    }
		  else
		    {
		      strcat(temp, binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(current_supernode[0]), db_graph->kmer_size, &tmp_kmer), 
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


		  //get coverages
		  get_coverage_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2,
                                                   len_branch2, &branch2_min_covg, &branch2_max_covg, &branch2_avg_covg, &branch2_mode_coverage, &branch2_percent_nodes_with_modal_covg,
                                                   which_array_holds_indiv, index_for_indiv_in_edge_array);



		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2, 
							len_branch2, &branch2_percent_novel,
							which_array_holds_ref, index_for_ref_in_edge_array);


		  //now do comparison with trusted branch - look for nodes that are in one but not the other
		  get_covg_of_nodes_in_one_but_not_other_of_two_arrays(chrom_path_array+start_node_index+length_5p_flank, 
								       current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2, 
								       len_trusted_branch, len_branch2,
		    						       &num_nodes_on_trusted_but_not_variant_branch, &num_nodes_on_variant_but_not_trusted_branch, 
								       ptrs_to_covgs_in_trusted_not_variant, ptrs_to_covgs_in_variant_not_trusted,
		    						       working_array1, working_array2,
		    						       which_array_holds_indiv, index_for_indiv_in_edge_array);


		  //ready to print trusted and var branches - will do it below, when have got 5p flank


		  //for unit tests
		  if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_trusted_branch_array[num_variants_found-1],trusted_branch);
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

	      //oK - we have all we need to print the 3prime flank now - will do so below
							       
	      
	      //for unit tests only
	      if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_flank3p_array[num_variants_found-1],flank3p);
		    }
	      


	      int start_of_var_branch_in_supernode_array;//coord of first node, may be at lh or rh end of section of array corresponding to this branch
	      int left_most_coord_of_var_branch_in_supernode_array; //if traversing right to left, this is the start of the 3prime anchor
	      if (traverse_sup_left_to_right==true)
		{
		  start_of_var_branch_in_supernode_array           = index_of_query_node_in_supernode_array+length_5p_flank;
		  left_most_coord_of_var_branch_in_supernode_array = start_of_var_branch_in_supernode_array;
		}
	      else
		{
		  start_of_var_branch_in_supernode_array           = index_of_query_node_in_supernode_array-length_5p_flank; 
		  left_most_coord_of_var_branch_in_supernode_array = start_of_3prime_anchor_in_sup;
		}
	      

	      set_variant_branches_and_flanks(&var, 
					      &chrom_path_array[start_node_index], 
					      &chrom_orientation_array[start_node_index], 
					      length_5p_flank,
					      &chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome], 
					      &chrom_orientation_array[first_index_in_chrom_where_supernode_differs_from_chromosome], 
					      len_trusted_branch,
					      &current_supernode[left_most_coord_of_var_branch_in_supernode_array],
					      &curr_sup_orientations[left_most_coord_of_var_branch_in_supernode_array],
					      len_branch2,
					      &chrom_path_array[start_of_3prime_anchor_in_chrom],
					      &chrom_orientation_array[start_of_3prime_anchor_in_chrom],
					      length_3p_flank, 
					      first);//first tells it that the trusted branch is the reference


	      if (condition(&var, index_for_ref_in_edge_array, index_for_indiv_in_edge_array)==true)
		{
		  //in some situations you want to mark all branches of variants for future use after this function returns
		  action_for_branches_of_called_variants(&var);

		  if (output_file != NULL)
		    {
		      
		      char name[300]; 
		      sprintf(name,"var_%i_5p_flank",num_variants_found);
		      print_fasta_from_path_for_specific_person_or_pop(output_file, name, length_5p_flank, 
								       flank5p_avg_covg,flank5p_min_covg,flank5p_max_covg, flank5p_mode_coverage, 
								       flank5p_percent_nodes_with_modal_covg, flank5p_percent_novel,
								       chrom_path_array[start_node_index],
								       chrom_orientation_array[start_node_index],
								       chrom_path_array[start_node_index+length_5p_flank-1],
								       chrom_orientation_array[start_node_index+length_5p_flank-1],
								       NULL, NULL, 0, //these 3 args are for when you want to find covg of nodes in this path and not another - not needed here
								       flank5p,
								       db_graph->kmer_size,
								       true,
								       which_array_holds_indiv, index_for_indiv_in_edge_array);
		      
		      
		      if (traverse_sup_left_to_right==true)
			{
			  sprintf(name,"var_%i_trusted_branch", num_variants_found);
			  
			  
			  print_fasta_from_path_for_specific_person_or_pop(output_file, name, len_trusted_branch, 
									   trusted_branch_avg_covg, trusted_branch_min_covg, trusted_branch_max_covg, 
									   trusted_branch_mode_coverage, trusted_branch_percent_nodes_with_modal_covg, trusted_branch_percent_novel,
									   chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
									   chrom_orientation_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
									   chrom_path_array[start_of_3prime_anchor_in_chrom],
									   chrom_orientation_array[start_of_3prime_anchor_in_chrom],
									   "covgs of trusted not variant nodes: ", covgs_in_trusted_not_variant, num_nodes_on_trusted_but_not_variant_branch, 
									   trusted_branch,
									   db_graph->kmer_size,
									   false,
									   which_array_holds_indiv, index_for_indiv_in_edge_array);
			  
			  
			  
			  sprintf(name,"var_%i_variant_branch", num_variants_found);
			  
			  print_fasta_from_path_for_specific_person_or_pop(output_file, name, len_branch2, 
									   branch2_avg_covg, branch2_min_covg, branch2_max_covg, branch2_mode_coverage, 
									   branch2_percent_nodes_with_modal_covg, branch2_percent_novel, 
									   current_supernode[start_of_var_branch_in_supernode_array], 
									   curr_sup_orientations[start_of_var_branch_in_supernode_array],
									   current_supernode[start_of_3prime_anchor_in_sup], 
									   curr_sup_orientations[start_of_3prime_anchor_in_sup],
									   "covgs of variant not trusted nodes: ", covgs_in_variant_not_trusted, num_nodes_on_variant_but_not_trusted_branch,
									   variant_branch, db_graph->kmer_size, false, 
									   which_array_holds_indiv, index_for_indiv_in_edge_array
									   );
			  
			  
			}
		      else //traversing right to left
			{
			  
			  sprintf(name,"var_%i_trusted_branch", num_variants_found);
			  
			  print_fasta_from_path_for_specific_person_or_pop(output_file, name, len_trusted_branch, 
									   trusted_branch_avg_covg, trusted_branch_min_covg, trusted_branch_max_covg, 
									   trusted_branch_mode_coverage, trusted_branch_percent_nodes_with_modal_covg, trusted_branch_percent_novel,
									   chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
									   chrom_orientation_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
									   chrom_path_array[start_of_3prime_anchor_in_chrom],
									   chrom_orientation_array[start_of_3prime_anchor_in_chrom],
									   "covgs of trusted not variant nodes: ",covgs_in_trusted_not_variant, num_nodes_on_trusted_but_not_variant_branch, 
									   trusted_branch,
									   db_graph->kmer_size,
									   false,
									   which_array_holds_indiv, index_for_indiv_in_edge_array);
			  
			  
			  char rev_variant_branch[len_branch2+1];
			  seq_reverse_complement(variant_branch, len_branch2, rev_variant_branch);
			  rev_variant_branch[len_branch2]='\0';
			  sprintf(name,"var_%i_variant_branch", num_variants_found);
			  
			  print_fasta_from_path_for_specific_person_or_pop(output_file, name, len_branch2, 
									   branch2_avg_covg, branch2_min_covg, branch2_max_covg, branch2_mode_coverage, 
									   branch2_percent_nodes_with_modal_covg, branch2_percent_novel,
									   current_supernode[start_of_var_branch_in_supernode_array], 
									   curr_sup_orientations[start_of_var_branch_in_supernode_array],
									   current_supernode[start_of_3prime_anchor_in_sup], 
									   curr_sup_orientations[start_of_3prime_anchor_in_sup],
									   "covgs of variant not trusted nodes: ", covgs_in_variant_not_trusted, num_nodes_on_variant_but_not_trusted_branch,
									   rev_variant_branch, db_graph->kmer_size, false, 
									   which_array_holds_indiv, index_for_indiv_in_edge_array
									   );
			  
			  //for unit tests only
			  if (num_variants_found<=max_desired_returns)
			    {
			      strcat(return_variant_branch_array[num_variants_found-1],rev_variant_branch);
			    }
			  
			}
		      sprintf(name,"var_%i_3p_flank",num_variants_found);
		      print_fasta_from_path_for_specific_person_or_pop(output_file, name, length_3p_flank, 
								       flank3p_avg_covg, flank3p_min_covg, flank3p_max_covg, flank3p_mode_coverage, 
								       flank3p_percent_nodes_with_modal_covg, flank3p_percent_novel, 
								       chrom_path_array[start_of_3prime_anchor_in_chrom],
								       chrom_orientation_array[start_of_3prime_anchor_in_chrom],
								       chrom_path_array[start_of_3prime_anchor_in_chrom+length_3p_flank],
								       chrom_orientation_array[start_of_3prime_anchor_in_chrom+length_3p_flank],
								       NULL, NULL, 0,
								       flank3p,
								       db_graph->kmer_size,
								       false,
								       which_array_holds_indiv, index_for_indiv_in_edge_array);
		      
		      
		      
		    }
		 
		  //fprintf(output_file, "var_%i - extra information\n", num_variants_found);
		  print_extra_info(&var, output_file);
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
  free(ptrs_to_covgs_in_trusted_not_variant);
  free(ptrs_to_covgs_in_variant_not_trusted);
  free(covgs_in_variant_not_trusted);
  free(covgs_in_trusted_not_variant);
  //  free(covgs_of_indiv_on_trusted_path);
  //free(covgs_of_ref_on_trusted_path);
  //free(covgs_of_indiv_on_variant_path);
  //free(covgs_of_ref_on_variant_path);
  free(working_array1);
  free(working_array2);
  return num_variants_found;


}



boolean make_reference_path_based_sv_calls_condition_always_true_in_subgraph_defined_by_func_of_colours(VariantBranchesAndFlanks* var, 
													int colour_of_ref,
													Edges (*get_colour)(const dBNode*),
													int (*get_covg)(const dBNode*))
{
  return true;
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
int db_graph_make_reference_path_based_sv_calls_in_subgraph_defined_by_func_of_colours(FILE* chrom_fasta_fptr,
										       Edges (*get_colour)(const dBNode*),
										       int (*get_covg)(const dBNode*),
										       int ref_colour, //EdgeArrayType which_array_holds_ref, int index_for_ref_in_edge_array,
										       int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, 
										       int max_anchor_span, int min_covg, int max_covg, 
										       int max_expected_size_of_supernode, int length_of_arrays, dBGraph* db_graph, FILE* output_file,
										       int max_desired_returns,
										       char** return_flank5p_array, char** return_trusted_branch_array, char** return_variant_branch_array, 
										       char** return_flank3p_array, int** return_variant_start_coord,
										       boolean (*condition)(VariantBranchesAndFlanks* var,  int colour_of_ref,  
													    Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*)),
										       void (*action_for_branches_of_called_variants)(VariantBranchesAndFlanks* var),
										       void (*print_extra_info)(VariantBranchesAndFlanks* var, FILE* fout)
										       )
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
  //kmer_window->kmer = (BinaryKmer*) malloc(sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER*length_of_arrays);    //*(number_of_nodes_to_load + db_graph->kmer_size));
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



  //finally, we need a couple of arrays that will hold dBNode*'s, which we wil reuse when comparing two branches
  dBNode** working_array1 = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 
  dBNode** working_array2 = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 

  int* covgs_in_trusted_not_variant = (int*) malloc(sizeof(int)*length_of_arrays);
  int* covgs_in_variant_not_trusted = (int*) malloc(sizeof(int)*max_expected_size_of_supernode);
  //int* covgs_of_indiv_on_trusted_path = (int*) malloc(sizeof(int)*length_of_arrays);
  //int* covgs_of_ref_on_trusted_path = (int*) malloc(sizeof(int)*length_of_arrays);
  //int* covgs_of_indiv_on_variant_path = (int*) malloc(sizeof(int)*length_of_arrays);
  //int* covgs_of_ref_on_variant_path = (int*) malloc(sizeof(int)*length_of_arrays);

  if ( (working_array1==NULL) || (working_array2==NULL)  ) 
    {
      printf("OOM at start - cannot alloc working arrays. Exit");
      exit(1);
    }

  if ( (covgs_in_trusted_not_variant==NULL) || (covgs_in_variant_not_trusted==NULL) )
    // || (covgs_of_indiv_on_trusted_path==NULL) || (covgs_of_ref_on_trusted_path==NULL)  || (covgs_of_indiv_on_variant_path==NULL) || (covgs_of_ref_on_variant_path==NULL) )
    {
      printf("OOM at start  - cannot malloc covg arrays. Exit\n");
      exit(1);
    }

  int** ptrs_to_covgs_in_trusted_not_variant = (int**) malloc(sizeof(int*)*length_of_arrays);
  int** ptrs_to_covgs_in_variant_not_trusted = (int**) malloc(sizeof(int*)*max_expected_size_of_supernode);

  if ( (ptrs_to_covgs_in_trusted_not_variant==NULL) || (ptrs_to_covgs_in_variant_not_trusted==NULL) )
    {
      printf("Cannot alloc ptr arrays for covgs in make_ref assisted calls fn");
      exit(1);
    }

  for (k=0; k< length_of_arrays; k++)
    {
      ptrs_to_covgs_in_trusted_not_variant[k]=&covgs_in_trusted_not_variant[k];
    }
  for (k=0; k< max_expected_size_of_supernode; k++)
    {
      ptrs_to_covgs_in_variant_not_trusted[k]=&covgs_in_variant_not_trusted[k];
    }


  VariantBranchesAndFlanks var; //we will reuse this


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
  BinaryKmer tmp_kmer;
  


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
	  if (db_node_is_this_node_in_subgraph_defined_by_func_of_colours(chrom_path_array[start_node_index], get_colour) == false)
	    {
	      //printf("ZAM DEBUG node number %d not in the graph of person %d\n", start_node_index, index_for_indiv_in_edge_array);
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
	  
	  int length_curr_supernode = 
	    db_graph_supernode_returning_query_node_posn_in_subgraph_defined_by_func_of_colours(chrom_path_array[start_node_index], 
												max_expected_size_of_supernode, &db_node_action_set_status_visited,
												current_supernode, curr_sup_orientations, curr_sup_labels, supernode_string, 
												&avg_coverage, &min_coverage, &max_coverage, &is_cycle, 
												&index_of_query_node_in_supernode_array, 
												db_graph, get_colour, get_covg);
	  

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
	      printf("WARNING  - Did not even go into loop working way along ref - so supernode only touches ref at initial node - I don't think this should be possible, should go into that loop once\nContact Zam - this bug has not been seen for over a year");
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
	      //fprintf(output_file, "VARIATION: %d, start coordinate %d\n", num_variants_found, start_coord_of_variant_in_trusted_path_fasta);
	      
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

	      //Our algorithm would say: OK, KingZam is the 5prime anchor, let's see if we can find a 3prime anchor that is 7 characters long (say). 
	      // Aha! England will do. OK - can we extend backwards? Yes, we can extend backwards 
	      // all the way back to Zam. But WAIT!! Zam is before the split!! That's fine, it means there is a repeat, just don't try and extend beyond there. 
	      //Hence the else if condition marked below with ***

	      //(By the way this is not a pathologival case - happens in Human chromosome 1)

	      while (can_extend_further_in_5prime_dir==true)  
		{
		  int k=j;

		  if (traverse_sup_left_to_right==false)
		    {
		      k=-j;
		    }
		  

		  if (chrom_path_array[start_of_3prime_anchor_in_chrom - j] != current_supernode[start_of_3prime_anchor_in_sup - k])
		    {
		      //this condition will cover cases when there is an N/Null in the trusted path/chrom_path_array
		      can_extend_further_in_5prime_dir=false;
		    }
		  else if (  (traverse_sup_left_to_right) && (start_of_3prime_anchor_in_sup - j <= index_in_supernode_where_supernode_differs_from_chromosome+1) ) 
		    // *** see comment above this while loop
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

	      
	      //printf("Collected following info:\n");
	      //printf("start_of_3prime_anchor_in_sup = %d\n", start_of_3prime_anchor_in_sup);
	      //printf("start_of_3prime_anchor_in_chrom = %d\n", start_of_3prime_anchor_in_chrom);
	      //printf("Length 3p flank is  %d\n", length_3p_flank);
	      //printf("Length 5p flank is  %d\n", length_5p_flank);
	      //printf("Length of arrays is %d\n", length_of_arrays);
	      //printf("Start node index is %d\n", start_node_index);
	      //printf("first_index_in_chrom_where_supernode_differs_from_chromosome is %d\n", first_index_in_chrom_where_supernode_differs_from_chromosome);
	      //printf("index_of_query_node_in_supernode_array is %d\n", index_of_query_node_in_supernode_array);
	      
	      
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
		  get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(current_supernode+index_of_query_node_in_supernode_array, 
											  length_5p_flank, &flank5p_min_covg, &flank5p_max_covg, &flank5p_avg_covg, 
											  &flank5p_mode_coverage, &flank5p_percent_nodes_with_modal_covg,
											  get_colour, get_covg);

		  //actually looks at the nodes and sees if they are in the graph of the person specified in the last 2 arguments - in this case, the reference => hence we find if novel
		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array, 
							length_5p_flank, &flank5p_percent_novel,
							individual_edge_array, ref_colour);
			


		}
	      else
		{
		  get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank, 
											  length_5p_flank, &flank5p_min_covg, &flank5p_max_covg, &flank5p_avg_covg, 
											  &flank5p_mode_coverage, &flank5p_percent_nodes_with_modal_covg,
											  get_colour, get_covg);


		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank, 
							length_5p_flank, &flank5p_percent_novel,
							individual_edge_array, ref_colour);

		}

	      
	      

	      // 5pflank info now ready for printing. Will be done below, after have collected info on branches and 3p flank


	      //this for testing only (ie used by unit tests)
	      if (num_variants_found<=max_desired_returns)
		{
		  //add first kmer onto flank5p - this si done automatically by print_fasta_from_path
		  if (chrom_orientation_array[start_node_index]==forward)
		    {
		      strcat(return_flank5p_array[num_variants_found-1],
			     binary_kmer_to_seq(element_get_kmer(chrom_path_array[start_node_index]), db_graph->kmer_size, tmp_zam));
		    }
		  else
		    {
		      char tmp_zam2[db_graph->kmer_size+1];
		      tmp_zam2[db_graph->kmer_size]='\0';

		      strcat(return_flank5p_array[num_variants_found-1],
			     seq_reverse_complement(binary_kmer_to_seq(element_get_kmer(chrom_path_array[start_node_index]), db_graph->kmer_size, tmp_zam), 
						    db_graph->kmer_size,tmp_zam2)
			     );
		    }
		      
		  strcat(return_flank5p_array[num_variants_found-1],flank5p);
		}

		    


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
	      int num_nodes_on_trusted_but_not_variant_branch=0;
	      int num_nodes_on_variant_but_not_trusted_branch=0;
	      
	      //get coverage on the trusted branch. 
	      get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(chrom_path_array+start_node_index+length_5p_flank,
										      len_trusted_branch, &trusted_branch_min_covg, &trusted_branch_max_covg, &trusted_branch_avg_covg, 
										      &trusted_branch_mode_coverage, &trusted_branch_percent_nodes_with_modal_covg,
										      get_colour, get_covg);

	      get_percent_novel_from_array_of_nodes(chrom_path_array+start_node_index+length_5p_flank, 
						    len_trusted_branch, &trusted_branch_percent_novel,
						    individual_edge_array, ref_colour);
	      



	      //almost ready to print trusted branch, but need to get coverages of nodes that are in trusted but not variant branch, and vice-versa. 
	      // will first set up variant branch, and then print the two together

	      

	      int branch2_min_covg=0;
	      int branch2_max_covg=0;
	      double branch2_avg_covg=0;
	      int branch2_mode_coverage=0;
	      double branch2_percent_nodes_with_modal_covg=0;
	      double branch2_percent_novel=0;

	      int len_branch2;

								 
              // ******************* variant branch - left to right
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

		  
		  //get coverages for variant branch
		  get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank,
											  len_branch2, &branch2_min_covg, &branch2_max_covg, &branch2_avg_covg, 
											  &branch2_mode_coverage, &branch2_percent_nodes_with_modal_covg,
											  get_colour, get_covg);

		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank, 
							len_branch2, &branch2_percent_novel,
							individual_edge_array, ref_colour);





		  //now do comparison with trusted branch - look for nodes that are in one but not the other
		  // data from this will be used in the next two calls to print_fasta, for the two branches - trusted, and variant
		  get_covg_of_nodes_in_one_but_not_other_of_two_arrays_in_subgraph_defined_by_func_of_colours(chrom_path_array+start_node_index+length_5p_flank, 
													      current_supernode+index_of_query_node_in_supernode_array+length_5p_flank, 
													      len_trusted_branch, len_branch2,
													      &num_nodes_on_trusted_but_not_variant_branch, 
													      &num_nodes_on_variant_but_not_trusted_branch, 
													      ptrs_to_covgs_in_trusted_not_variant, 
													      ptrs_to_covgs_in_variant_not_trusted,
													      working_array1, working_array2,
													      get_colour, get_covg);
		  //we are ready to print the trusted and variant paths, when variant path traversed left to right. Will do it after have collected info on 3p flank
		  


		  //for unit tests only
		  if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_trusted_branch_array[num_variants_found-1],trusted_branch);
		      strcat(return_variant_branch_array[num_variants_found-1],variant_branch);
		    }
		  


		}
	      else //************ variant branch - traversing supernode right to left
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
		      strcat(temp, binary_kmer_to_seq(element_get_kmer(current_supernode[0]), db_graph->kmer_size, tmp_zam));
		    }
		  else
		    {
		      strcat(temp, binary_kmer_to_seq(binary_kmer_reverse_complement(element_get_kmer(current_supernode[0]), db_graph->kmer_size, &tmp_kmer), 
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


		  //get coverages
		  get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2,
											  len_branch2, &branch2_min_covg, &branch2_max_covg, &branch2_avg_covg, 
											  &branch2_mode_coverage, &branch2_percent_nodes_with_modal_covg,
											  get_colour, get_covg);

		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2, 
							len_branch2, &branch2_percent_novel,
							individual_edge_array, ref_colour);


		  //now do comparison with trusted branch - look for nodes that are in one but not the other
		  get_covg_of_nodes_in_one_but_not_other_of_two_arrays_in_subgraph_defined_by_func_of_colours(chrom_path_array+start_node_index+length_5p_flank, 
													      current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2, 
													      len_trusted_branch, len_branch2,
													      &num_nodes_on_trusted_but_not_variant_branch, 
													      &num_nodes_on_variant_but_not_trusted_branch, 
													      ptrs_to_covgs_in_trusted_not_variant, ptrs_to_covgs_in_variant_not_trusted,
													      working_array1, working_array2,
													      get_colour, get_covg);


		  //ready to print trusted and var branches - will do it below, when have got 5p flank


		  //for unit tests
		  if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_trusted_branch_array[num_variants_found-1],trusted_branch);
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
		  get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank+len_branch2, 
											  length_3p_flank, &flank3p_min_covg, &flank3p_max_covg, &flank3p_avg_covg, 
											  &flank3p_mode_coverage, &flank3p_percent_nodes_with_modal_covg,
											  get_colour, get_covg);

		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array+length_5p_flank+len_branch2, 
							length_3p_flank, &flank3p_percent_novel,
							individual_edge_array, ref_colour);

				

		}
	      else
		{
		  get_coverage_from_array_of_nodes_in_subgraph_defined_by_func_of_colours(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2-length_3p_flank, 
											  length_3p_flank, &flank3p_min_covg, &flank3p_max_covg, &flank3p_avg_covg,  
											  &flank3p_mode_coverage, &flank3p_percent_nodes_with_modal_covg,
											  get_colour, get_covg);
		  
		  get_percent_novel_from_array_of_nodes(current_supernode+index_of_query_node_in_supernode_array-length_5p_flank-len_branch2-length_3p_flank, 
							length_3p_flank, &flank3p_percent_novel,
							individual_edge_array, ref_colour);
		}

	      //oK - we have all we need to print the 3prime flank now - will do so below
							       
	      
	      //for unit tests only
	      if (num_variants_found<=max_desired_returns)
		    {
		      strcat(return_flank3p_array[num_variants_found-1],flank3p);
		    }
	      


	      int start_of_var_branch_in_supernode_array;//coord of first node, may be at lh or rh end of section of array corresponding to this branch
	      int left_most_coord_of_var_branch_in_supernode_array; //if traversing right to left, this is the start of the 3prime anchor
	      if (traverse_sup_left_to_right==true)
		{
		  start_of_var_branch_in_supernode_array           = index_of_query_node_in_supernode_array+length_5p_flank;
		  left_most_coord_of_var_branch_in_supernode_array = start_of_var_branch_in_supernode_array;
		}
	      else
		{
		  start_of_var_branch_in_supernode_array           = index_of_query_node_in_supernode_array-length_5p_flank; 
		  left_most_coord_of_var_branch_in_supernode_array = start_of_3prime_anchor_in_sup;
		}
	      

	      set_variant_branches_and_flanks(&var, 
					      &chrom_path_array[start_node_index], 
					      &chrom_orientation_array[start_node_index], 
					      length_5p_flank,
					      &chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome], 
					      &chrom_orientation_array[first_index_in_chrom_where_supernode_differs_from_chromosome], 
					      len_trusted_branch,
					      &current_supernode[left_most_coord_of_var_branch_in_supernode_array],
					      &curr_sup_orientations[left_most_coord_of_var_branch_in_supernode_array],
					      len_branch2,
					      &chrom_path_array[start_of_3prime_anchor_in_chrom],
					      &chrom_orientation_array[start_of_3prime_anchor_in_chrom],
					      length_3p_flank, 
					      first);//first tells it that the trusted branch is the reference


	      if (condition(&var, ref_colour, get_colour, get_covg)==true)
		{
		  //in some situations you want to mark all branches of variants for future use after this function returns
		  action_for_branches_of_called_variants(&var);

		  if (output_file != NULL)
		    {
		      

		      char name[300]; 
		      sprintf(name,"var_%i_5p_flank",num_variants_found);
		      print_fasta_from_path_in_subgraph_defined_by_func_of_colours(output_file, name, length_5p_flank, 
										   flank5p_avg_covg,flank5p_min_covg,flank5p_max_covg, flank5p_mode_coverage, 
										   flank5p_percent_nodes_with_modal_covg, flank5p_percent_novel,
										   chrom_path_array[start_node_index],
										   chrom_orientation_array[start_node_index],
										   chrom_path_array[start_node_index+length_5p_flank-1],
										   chrom_orientation_array[start_node_index+length_5p_flank-1],
										   NULL, NULL, 0, 
										      //the above 3 args are for when you want to find covg of nodes in this path and not another - not needed here
										   flank5p,
										   db_graph->kmer_size,
										   true,
										   get_colour, get_covg);
		      
		      
		      if (traverse_sup_left_to_right==true)
			{
			  sprintf(name,"var_%i_trusted_branch", num_variants_found);
			  
			  
			  print_fasta_from_path_in_subgraph_defined_by_func_of_colours(output_file, name, len_trusted_branch, 
									   trusted_branch_avg_covg, trusted_branch_min_covg, trusted_branch_max_covg, 
									   trusted_branch_mode_coverage, trusted_branch_percent_nodes_with_modal_covg, trusted_branch_percent_novel,
									   chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
									   chrom_orientation_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
									   chrom_path_array[start_of_3prime_anchor_in_chrom],
									   chrom_orientation_array[start_of_3prime_anchor_in_chrom],
									   "covgs of trusted not variant nodes: ", covgs_in_trusted_not_variant, num_nodes_on_trusted_but_not_variant_branch, 
									   trusted_branch,
									   db_graph->kmer_size,
									   false,
									   get_colour, get_covg);
			  
			  
			  
			  sprintf(name,"var_%i_variant_branch", num_variants_found);
			  
			  print_fasta_from_path_in_subgraph_defined_by_func_of_colours(output_file, name, len_branch2, 
									   branch2_avg_covg, branch2_min_covg, branch2_max_covg, branch2_mode_coverage, 
									   branch2_percent_nodes_with_modal_covg, branch2_percent_novel, 
									   current_supernode[start_of_var_branch_in_supernode_array], 
									   curr_sup_orientations[start_of_var_branch_in_supernode_array],
									   current_supernode[start_of_3prime_anchor_in_sup], 
									   curr_sup_orientations[start_of_3prime_anchor_in_sup],
									   "covgs of variant not trusted nodes: ", covgs_in_variant_not_trusted, num_nodes_on_variant_but_not_trusted_branch,
									   variant_branch, db_graph->kmer_size, false, 
									   get_colour, get_covg
									   );
			  
			  
			}
		      else //traversing right to left
			{
			  
			  sprintf(name,"var_%i_trusted_branch", num_variants_found);
			  
			  print_fasta_from_path_in_subgraph_defined_by_func_of_colours(output_file, name, len_trusted_branch, 
									   trusted_branch_avg_covg, trusted_branch_min_covg, trusted_branch_max_covg, 
									   trusted_branch_mode_coverage, trusted_branch_percent_nodes_with_modal_covg, trusted_branch_percent_novel,
									   chrom_path_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
									   chrom_orientation_array[first_index_in_chrom_where_supernode_differs_from_chromosome],
									   chrom_path_array[start_of_3prime_anchor_in_chrom],
									   chrom_orientation_array[start_of_3prime_anchor_in_chrom],
									   "covgs of trusted not variant nodes: ",covgs_in_trusted_not_variant, num_nodes_on_trusted_but_not_variant_branch, 
									   trusted_branch,
									   db_graph->kmer_size,
									   false,
									   get_colour, get_covg);
			  
			  
			  char rev_variant_branch[len_branch2+1];
			  seq_reverse_complement(variant_branch, len_branch2, rev_variant_branch);
			  rev_variant_branch[len_branch2]='\0';
			  sprintf(name,"var_%i_variant_branch", num_variants_found);
			  
			  print_fasta_from_path_in_subgraph_defined_by_func_of_colours(output_file, name, len_branch2, 
									   branch2_avg_covg, branch2_min_covg, branch2_max_covg, branch2_mode_coverage, 
									   branch2_percent_nodes_with_modal_covg, branch2_percent_novel,
									   current_supernode[start_of_var_branch_in_supernode_array], 
									   curr_sup_orientations[start_of_var_branch_in_supernode_array],
									   current_supernode[start_of_3prime_anchor_in_sup], 
									   curr_sup_orientations[start_of_3prime_anchor_in_sup],
									   "covgs of variant not trusted nodes: ", covgs_in_variant_not_trusted, num_nodes_on_variant_but_not_trusted_branch,
									   rev_variant_branch, db_graph->kmer_size, false, 
									   get_colour, get_covg
									   );
			  
			  //for unit tests only
			  if (num_variants_found<=max_desired_returns)
			    {
			      strcat(return_variant_branch_array[num_variants_found-1],rev_variant_branch);
			    }
			  
			}
		      sprintf(name,"var_%i_3p_flank",num_variants_found);
		      print_fasta_from_path_in_subgraph_defined_by_func_of_colours(output_file, name, length_3p_flank, 
								       flank3p_avg_covg, flank3p_min_covg, flank3p_max_covg, flank3p_mode_coverage, 
								       flank3p_percent_nodes_with_modal_covg, flank3p_percent_novel, 
								       chrom_path_array[start_of_3prime_anchor_in_chrom],
								       chrom_orientation_array[start_of_3prime_anchor_in_chrom],
								       chrom_path_array[start_of_3prime_anchor_in_chrom+length_3p_flank],
								       chrom_orientation_array[start_of_3prime_anchor_in_chrom+length_3p_flank],
								       NULL, NULL, 0,
								       flank3p,
								       db_graph->kmer_size,
								       false,
								       get_colour, get_covg);
		      
		      
		      
		    }
		 
		  //fprintf(output_file, "var_%i - extra information\n", num_variants_found);
		  print_extra_info(&var, output_file);
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
  free(ptrs_to_covgs_in_trusted_not_variant);
  free(ptrs_to_covgs_in_variant_not_trusted);
  free(covgs_in_variant_not_trusted);
  free(covgs_in_trusted_not_variant);
  //  free(covgs_of_indiv_on_trusted_path);
  //free(covgs_of_ref_on_trusted_path);
  //free(covgs_of_indiv_on_variant_path);
  //free(covgs_of_ref_on_variant_path);
  free(working_array1);
  free(working_array2);
  return num_variants_found;


}




void db_graph_make_reference_path_based_sv_calls_given_list_of_colours_for_indiv(int* list, int len_list,
										 FILE* chrom_fasta_fptr, int ref_colour,
										 int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, 
										 int max_anchor_span, int min_covg, int max_covg, 
										 int max_expected_size_of_supernode, int length_of_arrays, dBGraph* db_graph, FILE* output_file,
										 int max_desired_returns,
										 char** return_flank5p_array, char** return_trusted_branch_array, char** return_variant_branch_array, 
										 char** return_flank3p_array, int** return_variant_start_coord,
										 boolean (*condition)(VariantBranchesAndFlanks* var,  int colour_of_ref,  
												      Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*)),
										 void (*action_for_branches_of_called_variants)(VariantBranchesAndFlanks* var),
										 void (*print_extra_info)(VariantBranchesAndFlanks* var, FILE* fout)
										 )
{


  Edges get_union_list_colours(const dBNode* e)
  {
    int i;
    Edges edges=0;
    
    for (i=0; i< len_list; i++)
      {
	edges |= e->individual_edges[list[i]];
      }
    return edges;    
  }

  int get_covg_of_union_first_list_colours(const dBNode* e)
  {
    int i;
    int covg=0;
  
    for (i=0; i< len_list; i++)
      {
	covg += e->coverage[list[i]];
      }
    return covg;
  }

  db_graph_make_reference_path_based_sv_calls_in_subgraph_defined_by_func_of_colours(chrom_fasta_fptr,
										     &get_union_list_colours, &get_covg_of_union_first_list_colours,
										     ref_colour,
										     min_fiveprime_flank_anchor, min_threeprime_flank_anchor,
										     max_anchor_span, min_covg, max_covg,
										     max_expected_size_of_supernode, length_of_arrays, db_graph, output_file,
										     max_desired_returns,
										     return_flank5p_array, return_trusted_branch_array, return_variant_branch_array,
										     return_flank3p_array, return_variant_start_coord,
										     condition, action_for_branches_of_called_variants, print_extra_info);
										     



}







/*
int db_graph_make_reference_path_based_sv_calls_after_marking_vars_in_ref_to_be_ignored(char* chrom_fasta, EdgeArrayType which_array_holds_indiv, int index_for_indiv_in_edge_array,
											EdgeArrayType which_array_holds_ref, int  index_for_ref_in_edge_array,
											int min_fiveprime_flank_anchor, int min_threeprime_flank_anchor, 
											int max_anchor_span, int min_covg, int max_covg,
											int max_expected_size_of_supernode, int length_of_arrays, dBGraph* db_graph, FILE* output_file,
											int max_desired_returns,
											char** return_flank5p_array, char** return_trusted_branch_array, char** return_variant_branch_array,
											char** return_flank3p_array, int** return_variant_start_coord,
											boolean (*condition)(VariantBranchesAndFlanks* var,  int colour_of_ref,  int colour_of_indiv)
											)
{

  FILE* fptr = fopen(chrom_fasta, "r");
  if (fptr==NULL)
    {
      printf("Cannot open %s the first tme in  db_graph_make_reference_path_based_sv_calls_after_marking_vars_in_ref_to_be_ignored", chrom_fasta);
      exit(1);
    }
  //call variants in the reference colour
  int num_false_vars_masked = db_graph_make_reference_path_based_sv_calls(fptr, which_array_holds_ref, index_for_ref_in_edge_array, which_array_holds_ref, index_for_ref_in_edge_array, 
									  min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
									  max_expected_size_of_supernode, length_of_arrays, db_graph, NULL, 0,
									  NULL, NULL, NULL, NULL, NULL,
									  &make_reference_path_based_sv_calls_condition_always_true, &action_set_flanks_and_branches_to_be_ignored);
  fclose(fptr);

  printf("Masked %d false variants by calling in ref colour\n", num_false_vars_masked);
  //unset visited, but leaving nodes marked as to_be_ignored if they were in variants in the ref colour
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);


  fptr = fopen(chrom_fasta, "r");
  if (fptr==NULL)
    {
      printf("Cannot open %s the second tme in  db_graph_make_reference_path_based_sv_calls_after_marking_vars_in_ref_to_be_ignored", chrom_fasta);
      exit(1);
    }
  int num_vars_found = db_graph_make_reference_path_based_sv_calls(fptr, which_array_holds_indiv, index_for_indiv_in_edge_array, which_array_holds_ref, index_for_ref_in_edge_array, 
								   min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
								   max_expected_size_of_supernode, length_of_arrays, db_graph, output_file, max_desired_returns,
								   return_flank5p_array, return_trusted_branch_array, return_variant_branch_array, return_flank3p_array, return_variant_start_coord, 
								   condition, &db_variant_action_do_nothing);
  fclose(fptr);
  return num_vars_found;
  
}
*/


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

void compute_label_in_subgraph_defined_by_func_of_colours(dBNode * node, Orientation o, char * label, Edges (*get_colour)(const dBNode*)){

  int i=0; 

  void check_nucleotide(Nucleotide n){

    if (db_node_edge_exist_within_specified_function_of_coloured_edges(node,n,o, get_colour)){
	label[i] =  binary_nucleotide_to_char(n);
	i++;
      }
  }
  
  nucleotide_iterator(check_nucleotide);

  label[i]='\0';
}


//sometimes we want to apply a function to all nodes in a path, but the path is so long that we don't
//want to have it all in memory. So we read through the fasta defining the path.
// caller must have an idea of size of fasta, and suggest a chunk_size (ideally <0.5 of the size of the file) in bases.
//  So chunk size wil be the number of bases/nodes loaded each time we go through the internal file-reading loop in this function
//  We assume chunk size is less than half the entire file.
void apply_to_all_nodes_in_path_defined_by_fasta(void (*func)(dBNode*), FILE* fasta_fptr, int chunk_size, dBGraph* db_graph)
{

  int length_of_arrays=2*chunk_size;

  //malloc arrays
  
  dBNode**     path_array        = (dBNode**)    malloc(sizeof(dBNode*)*length_of_arrays); //everything these dBNode*'s point to will be in the hash table - ie owned by the hash table
  Orientation* orientation_array = (Orientation*)malloc(sizeof(Orientation)*length_of_arrays); 
  Nucleotide*  label_array       = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  char*        path_string       = (char*)       malloc(sizeof(char)*(length_of_arrays+1)); //+1 for \0


  int n;
  for (n=0; n<length_of_arrays; n++)
    {
      path_array[n]=NULL;
      orientation_array[n]=forward;
      label_array[n]=Undefined;
      path_string[n]='N';
    }

  path_string[0]='\0';
  path_string[length_of_arrays]='\0';


  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence in apply_to_all_nodes_in_path_defined_by_fasta\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,chunk_size+db_graph->kmer_size+1,LINE_MAX);
  seq->seq[0]='\0';


  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in apply_to_all_nodes_in_path_defined_by_fasta. Exit.\n");
      exit(1);
    }

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*length_of_arrays);   
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in apply_to_all_nodes_in_path_defined_by_fastadb_graph. Exit.\n");
      exit(1);
    }  
  kmer_window->nkmers=0;




  //load a set of nodes into thr back end (right hand/greatest indices) of array
  // each call   will push left the nodes/etc in the various arrays by length_of_arrays
  // and then put the new nodes etc in on the right of that

  int total_so_far=0;
  int ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(fasta_fptr, chunk_size, 0, 
													     length_of_arrays,
													     path_array, orientation_array, 
													     label_array, path_string,
													     seq, kmer_window, 
													     true, false,
													     db_graph);
  total_so_far+=ret;


  if (ret == chunk_size)
    {
      //one more batch, then array is full, and ready for the main loop:
      ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(fasta_fptr, chunk_size, chunk_size, 
													     length_of_arrays,
													     path_array, orientation_array, label_array, path_string,
													     seq, kmer_window, 
													     false, false,
													     db_graph);
      total_so_far+=ret;
    }



  //apply function to first half of array
  for (n=0; n<chunk_size; n++)
    {
      func(path_array[n]);
    }


  //now in this loop keep applying function to the back half of the array
  do
    {
      for (n=chunk_size; n<2*chunk_size; n++)
	{
	  func(path_array[n]);
	}

    }
  while (db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta
	 (fasta_fptr, chunk_size, chunk_size, length_of_arrays, path_array, orientation_array, label_array, path_string,
	  seq, kmer_window, false, true, db_graph)>0);




  free(path_array);
  free(orientation_array);
  free(label_array);
  free(path_string);
  free_sequence(&seq);
  free(kmer_window->kmer);
  free(kmer_window);





}



void print_no_extra_info(VariantBranchesAndFlanks* var, FILE* fout)
{
}
void print_no_extra_supernode_info(dBNode** node_array, Orientation* or_array, int len, FILE* fout)
{
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
                                                      char* text_describing_comparison_with_other_path,//may be NULL -use to allow printing coveragesof nodes in this path but not in some specific other path
                                                      int* coverages_nodes_in_this_path_but_not_some_other, //may be NULL
                                                      int length_of_coverage_array,
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
 
      BinaryKmer fst_kmer;
      BinaryKmer tmp_kmer;

      if (fst_orientation==reverse)
	{
	  binary_kmer_assignment_operator(fst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(fst_node),kmer_size, &tmp_kmer) ) );
	} 
      else
	{
	  binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)) );
	}
      
      binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq);
      
      BinaryKmer lst_kmer;

      if (lst_orientation==reverse){
	binary_kmer_assignment_operator(lst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(lst_node),kmer_size, &tmp_kmer) ) );
      } 
      else
	{
	  binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)) );
	}
      binary_kmer_to_seq(&lst_kmer,kmer_size,lst_seq);
      

      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i mode_coverage: %i percent_nodes_with_modal_covg: %5.2f percent_novel: %5.2f fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s ", name,
	      (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage, modal_coverage, percent_nodes_with_modal_coverage, percent_novel,
	      db_node_get_coverage(fst_node, type, index),
	      fst_seq,
	      (fst_orientation == forward ? fst_r : fst_f),
	      (fst_orientation == forward ? fst_f : fst_r),
	      db_node_get_coverage(lst_node, type, index),
	      lst_seq,
	      (lst_orientation == forward ? lst_r : lst_f),
	      (lst_orientation == forward ? lst_f : lst_r));

      // often want to print out coverages of nodes that are in this path, but not in some other path
      if ((text_describing_comparison_with_other_path!=NULL)&& (coverages_nodes_in_this_path_but_not_some_other!=NULL) && (length_of_coverage_array>0))
	{
	  fprintf(fout, "%s ", text_describing_comparison_with_other_path);
	  int i;
	  for (i=0; i<length_of_coverage_array; i++)
	    {
	      fprintf(fout, "%d ", coverages_nodes_in_this_path_but_not_some_other[i]);
	    } 
	  fprintf(fout, "number_of_such_nodes: %d", length_of_coverage_array);
	  
	}

      // now print the nodes of all of the nodes in the array, in order

      fprintf(fout, "\n");

      if (include_first_kmer){    
	
	fprintf(fout,"%s",binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq));
      }
      
      fprintf(fout,"%s\n",string);
      
    }


}

void print_fasta_from_path_in_subgraph_defined_by_func_of_colours(FILE *fout,
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
								  char* text_describing_comparison_with_other_path, 
								     //text_describing_comparison_with_other_path may be NULL 
								     // - use to allow printing coverages of nodes in this path but not in some specific other path
								  int* coverages_nodes_in_this_path_but_not_some_other, //may be NULL
								  int length_of_coverage_array,
								  char * string, //labels of paths
								  int kmer_size,
								  boolean include_first_kmer,
								  Edges (*get_colour)(const dBNode*),
								  int (*get_covg)(const dBNode*)
								  )
{
  if (fout==NULL)
    {
      printf("Exiting - have passed a null file pointer to print_fasta_from_path_for_specific_person_or_pop\n");
      exit(1);
    }
  
  if ( (fst_node==NULL) || (lst_node==NULL) )
    {
      fprintf(fout, "WARNING - print_fasta command has been given a NULL node as first or last node - will not print fst/lst kmers.\n");
      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i mode_coverage: %i percent_nodes_with_modal_covg: %5.2f percent_novel: %5.2f\n", 
              name, (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage, modal_coverage, percent_nodes_with_modal_coverage, percent_novel);
      fprintf(fout,"%s\n",string);
    }
  else
    {
      char fst_f[5], fst_r[5],lst_f[5],lst_r[5];
      
      compute_label_in_subgraph_defined_by_func_of_colours(fst_node,forward,fst_f, get_colour);
      compute_label_in_subgraph_defined_by_func_of_colours(fst_node,reverse,fst_r, get_colour);
      compute_label_in_subgraph_defined_by_func_of_colours(lst_node,forward,lst_f, get_colour);
      compute_label_in_subgraph_defined_by_func_of_colours(lst_node,reverse,lst_r, get_colour);

      char fst_seq[kmer_size+1], lst_seq[kmer_size+1];
 
      BinaryKmer fst_kmer;
      BinaryKmer tmp_kmer;

      if (fst_orientation==reverse)
	{
	  binary_kmer_assignment_operator(fst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(fst_node),kmer_size, &tmp_kmer) ) );
	} 
      else
	{
	  binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)) );
	}
      
      binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq);
      
      BinaryKmer lst_kmer;

      if (lst_orientation==reverse){
	binary_kmer_assignment_operator(lst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(lst_node),kmer_size, &tmp_kmer) ) );
      } 
      else
	{
	  binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)) );
	}
      binary_kmer_to_seq(&lst_kmer,kmer_size,lst_seq);
      

      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i mode_coverage: %i percent_nodes_with_modal_covg: %5.2f percent_novel: %5.2f fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s ", name,
	      (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage, modal_coverage, percent_nodes_with_modal_coverage, percent_novel,
	      db_node_get_coverage_in_subgraph_defined_by_func_of_colours(fst_node, get_covg),
	      fst_seq,
	      (fst_orientation == forward ? fst_r : fst_f),
	      (fst_orientation == forward ? fst_f : fst_r),
	      db_node_get_coverage_in_subgraph_defined_by_func_of_colours(lst_node, get_covg),
	      lst_seq,
	      (lst_orientation == forward ? lst_r : lst_f),
	      (lst_orientation == forward ? lst_f : lst_r));
      
      // often want to print out coverages of nodes that are in this path, but not in some other path
      if ((text_describing_comparison_with_other_path!=NULL)&& (coverages_nodes_in_this_path_but_not_some_other!=NULL) && (length_of_coverage_array>0))
	{
	  fprintf(fout, "%s ", text_describing_comparison_with_other_path);
	  int i;
	  for (i=0; i<length_of_coverage_array; i++)
	    {
	      fprintf(fout, "%d ", coverages_nodes_in_this_path_but_not_some_other[i]);
	    } 
	  fprintf(fout, "number_of_such_nodes: %d", length_of_coverage_array);
	  
	}

      // now print the nodes of all of the nodes in the array, in order

      fprintf(fout, "\n");

      if (include_first_kmer){    
	
	fprintf(fout,"%s",binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq));
      }
      
      fprintf(fout,"%s\n",string);
      
    }


}





void print_fasta_with_all_coverages_from_path_for_specific_person_or_pop(FILE *fout,
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
									 char* text_describing_comparison_with_other_path,//may be NULL -use to allow printing coveragesof nodes in this path but not in some specific other path
									 int* coverages_nodes_in_this_path_but_not_some_other, //may be NULL
									 int length_of_coverage_array,//refers to prev argument
									 int* coverages_nodes_in_this_path, 
									 int* coverages_in_ref_nodes_in_this_path,
									 int number_nodes_in_this_path,
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


      BinaryKmer fst_kmer;
      BinaryKmer tmp_kmer;

      if (fst_orientation==reverse)
	{
	  binary_kmer_assignment_operator(fst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(fst_node),kmer_size, &tmp_kmer) ) );
	} 
      else
	{
	  binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)) );
	}
      
      binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq);
      
      BinaryKmer lst_kmer;

      if (lst_orientation==reverse){
	binary_kmer_assignment_operator(lst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(lst_node),kmer_size, &tmp_kmer) ) );
      } 
      else
	{
	  binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)) );
	}
      binary_kmer_to_seq(&lst_kmer,kmer_size,lst_seq);

      
      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i mode_coverage: %i percent_nodes_with_modal_covg: %5.2f percent_novel: %5.2f fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s ", name,
	      (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage, modal_coverage, percent_nodes_with_modal_coverage, percent_novel,
	      db_node_get_coverage(fst_node, type, index),
	      fst_seq,
	      (fst_orientation == forward ? fst_r : fst_f),
	      (fst_orientation == forward ? fst_f : fst_r),
	      db_node_get_coverage(lst_node, type, index),
	      lst_seq,
	      (lst_orientation == forward ? lst_r : lst_f),
	      (lst_orientation == forward ? lst_f : lst_r));

      // often want to print out coverages of nodes that are in this path, but not in some other path
      if ((text_describing_comparison_with_other_path!=NULL)&& (coverages_nodes_in_this_path_but_not_some_other!=NULL) && (length_of_coverage_array>0))
	{
	  fprintf(fout, "%s ", text_describing_comparison_with_other_path);
	  int i;
	  for (i=0; i<length_of_coverage_array; i++)
	    {
	      fprintf(fout, "%d ", coverages_nodes_in_this_path_but_not_some_other[i]);
	    } 
	  fprintf(fout, "number_of_such_nodes: %d", length_of_coverage_array);
	  
	}

      // now print the nodes of all of the nodes in the array, in order
      fprintf(fout, "Coverages in indiv of all nodes in path in order:  ");
      int i;
      for (i=0; i<number_nodes_in_this_path; i++)
	{
	  fprintf(fout, "%d ", coverages_nodes_in_this_path[i]);
	} 
      fprintf(fout, "Coverages in ref of all nodes in path in order:  ");
      for (i=0; i<number_nodes_in_this_path; i++)
	{
          fprintf(fout, "%d ", coverages_in_ref_nodes_in_this_path[i]);
	}

      
      fprintf(fout, "\n");

      if (include_first_kmer){    
	
	fprintf(fout,"%s",binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq));
      }
      
      fprintf(fout,"%s\n",string);
      
    }


}


void print_minimal_fasta_from_path_for_specific_person_or_pop(FILE *fout,
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
							      boolean include_first_kmer,
							      EdgeArrayType type,
							      int index
							      )
{

  if (fout==NULL)
    {
      printf("Exiting - have passed a null file pointer to print_fasta_from_path_for_specific_person_or_pop\n");
      exit(1);
    }
  
  if ( (fst_node==NULL) || (lst_node==NULL) )
    {
      fprintf(fout, "WARNING - print_fasta command as have been given a NULL node as first or last node - will not print fst/lst kmers.\n");
      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i\n", 
              name, (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage);
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
 
      BinaryKmer fst_kmer;
      BinaryKmer tmp_kmer;

      if (fst_orientation==reverse)
	{
	  binary_kmer_assignment_operator(fst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(fst_node),kmer_size, &tmp_kmer) ) );
	} 
      else
	{
	  binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)) );
	}
      
      binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq);
      
      BinaryKmer lst_kmer;

      if (lst_orientation==reverse){
	binary_kmer_assignment_operator(lst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(lst_node),kmer_size, &tmp_kmer) ) );
      } 
      else
	{
	  binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)) );
	}
      binary_kmer_to_seq(&lst_kmer,kmer_size,lst_seq);
      

      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s ", name,
	      (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage, 
	      db_node_get_coverage(fst_node, type, index),
	      fst_seq,
	      (fst_orientation == forward ? fst_r : fst_f),
	      (fst_orientation == forward ? fst_f : fst_r),
	      db_node_get_coverage(lst_node, type, index),
	      lst_seq,
	      (lst_orientation == forward ? lst_r : lst_f),
	      (lst_orientation == forward ? lst_f : lst_r));



      fprintf(fout, "\n");

      if (include_first_kmer){    
	
	fprintf(fout,"%s",binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq));
      }
      
      fprintf(fout,"%s\n",string);
      
    }



}



void print_minimal_fasta_from_path_in_subgraph_defined_by_func_of_colours(FILE *fout,
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
									  boolean include_first_kmer,
									  Edges (*get_colour)(const dBNode*),
									  int (*get_covg)(const dBNode*)
									  )
{

  if (fout==NULL)
    {
      printf("Exiting - have passed a null file pointer to print_fasta_from_path_for_specific_person_or_pop\n");
      exit(1);
    }
  
  if ( (fst_node==NULL) || (lst_node==NULL) )
    {
      fprintf(fout, "WARNING - print_fasta command as have been given a NULL node as first or last node - will not print fst/lst kmers.\n");
      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i\n", 
              name, (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage);
      fprintf(fout,"%s\n",string);
    }
  else
    {
      char fst_f[5], fst_r[5],lst_f[5],lst_r[5];
      
      compute_label_in_subgraph_defined_by_func_of_colours(fst_node,forward,fst_f, get_colour);
      compute_label_in_subgraph_defined_by_func_of_colours(fst_node,reverse,fst_r, get_colour);
      compute_label_in_subgraph_defined_by_func_of_colours(lst_node,forward,lst_f, get_colour);
      compute_label_in_subgraph_defined_by_func_of_colours(lst_node,reverse,lst_r, get_colour);
      
      char fst_seq[kmer_size+1], lst_seq[kmer_size+1];
 
      BinaryKmer fst_kmer;
      BinaryKmer tmp_kmer;

      if (fst_orientation==reverse)
	{
	  binary_kmer_assignment_operator(fst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(fst_node),kmer_size, &tmp_kmer) ) );
	} 
      else
	{
	  binary_kmer_assignment_operator(fst_kmer, *(element_get_kmer(fst_node)) );
	}
      
      binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq);
      
      BinaryKmer lst_kmer;

      if (lst_orientation==reverse){
	binary_kmer_assignment_operator(lst_kmer, *(  binary_kmer_reverse_complement(element_get_kmer(lst_node),kmer_size, &tmp_kmer) ) );
      } 
      else
	{
	  binary_kmer_assignment_operator(lst_kmer, *(element_get_kmer(lst_node)) );
	}
      binary_kmer_to_seq(&lst_kmer,kmer_size,lst_seq);
      

      fprintf(fout,">%s length:%i average_coverage:%5.2f min_coverage:%i max_coverage:%i fst_coverage:%i fst_kmer:%s fst_r:%s fst_f:%s lst_coverage:%i lst_kmer:%s lst_r:%s lst_f:%s ", name,
	      (include_first_kmer ? length+kmer_size:length),avg_coverage,min_coverage,max_coverage, 
	      db_node_get_coverage_in_subgraph_defined_by_func_of_colours(fst_node, get_covg),
	      fst_seq,
	      (fst_orientation == forward ? fst_r : fst_f),
	      (fst_orientation == forward ? fst_f : fst_r),
	      db_node_get_coverage_in_subgraph_defined_by_func_of_colours(lst_node, get_covg),
	      lst_seq,
	      (lst_orientation == forward ? lst_r : lst_f),
	      (lst_orientation == forward ? lst_f : lst_r));



      fprintf(fout, "\n");

      if (include_first_kmer){    
	
	fprintf(fout,"%s",binary_kmer_to_seq(&fst_kmer,kmer_size,fst_seq));
      }
      
      fprintf(fout,"%s\n",string);
      
    }



}



boolean does_this_path_exist_in_this_colour(dBNode** array_nodes, Orientation* array_orientations,  int len, Edges (*get_colour)(const dBNode*), dBGraph* db_graph )
{
  boolean ret = true;
  int i;

  if (array_nodes[0]==NULL)
    {
      return false;
    }

  for (i=1; i<len; i++)
    {
      //get last base in kmer
      Nucleotide n;

      if (array_nodes[i]==NULL)
	{
	  return false;
	}

      if (array_orientations[i]==forward)
	{
	  n = binary_kmer_get_last_nucleotide(&(array_nodes[i]->kmer));
	}
      else
	{
	  BinaryKmer tmp_kmer;
	  n = binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(&(array_nodes[i]->kmer), db_graph->kmer_size, &tmp_kmer));
	}
      
      if (!(db_node_edge_exist_within_specified_function_of_coloured_edges(array_nodes[i-1], n, array_orientations[i-1], get_colour)) )
	{
	  ret=false;
	  break;
	}
      
    }

  return ret;
  
  
}


void print_standard_extra_supernode_info(dBNode** node_array, Orientation* or_array, int len, FILE* fout)
{
  
  int col;
  for (col=0; col<NUMBER_OF_COLOURS; col++)
    {
      
      fprintf(fout, "Covg in Colour %d:\n", col);
      int i;
      for (i=0; i<len; i++)
	{
	  if (node_array[i]!=NULL)
	    {
	      fprintf(fout, "%d ", node_array[i]->coverage[col]);
	    }
	  else
	    {
	      fprintf(fout, "0 ");
	    }
	}
      fprintf(fout, "\n");
      
    }
  
}

void print_standard_extra_info(VariantBranchesAndFlanks* var, FILE* fout)
{
  fprintf(fout, "\n");
  //print coverages:
  fprintf(fout, "branch1 coverages\n");
  print_standard_extra_supernode_info(var->one_allele, var->one_allele_or, var->len_one_allele, fout);
  fprintf(fout, "branch2 coverages\n");
  print_standard_extra_supernode_info(var->other_allele, var->other_allele_or, var->len_other_allele, fout);
  fprintf(fout, "\n\n");
}



//check all edges in graph
long long db_graph_health_check(boolean fix, dBGraph * db_graph){
  dBNode * next_node;
  Nucleotide reverse_nucleotide;
  Orientation next_orientation;
  long long count_nodes=0;
  char tmp_seq[db_graph->kmer_size+1]; 

  void check_node_with_orientation(dBNode * node, Orientation orientation){

    int j;
    for (j=0; j< NUMBER_OF_COLOURS; j++)
      {

	void check_base(Nucleotide n)
	{     
	  if (db_node_edge_exist(node,n,orientation, individual_edge_array, j)==true){

	    next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,n,&reverse_nucleotide,db_graph, individual_edge_array, j);
	    
	    if(next_node == NULL)
	      {
		printf("Health check problem -  didnt find node in hash table: %s %c %s\n",
		       binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq),
		       binary_nucleotide_to_char(n), 
		       orientation == forward ? "forward" : "reverse");
		if (fix)
		  {
		    db_node_reset_edge(node,n,orientation, individual_edge_array, j);
		  }
	      }
	    else
	      {
		if (db_node_edge_exist(next_node,reverse_nucleotide,opposite_orientation(next_orientation), individual_edge_array, j)==false)
		  {
		    printf("Health check problem - inconsitency return arrow missing: %s %c %s\n",
			   binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,tmp_seq),
			   binary_nucleotide_to_char(n), 
			   orientation == forward ? "forward" : "reverse");
		    if (fix)
		      {
			db_node_reset_edge(node,orientation, n, individual_edge_array, j);
		      }
		  }
	      }
	    
	  }
	}

	nucleotide_iterator(&check_base);
      }
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



//wipes a colour clean - for all nodes, sets covg=0, edge=0.
void db_graph_wipe_colour(int colour, dBGraph* db_graph)
{
  void wipe_node(dBNode* node)
  {
    node->individual_edges[colour]=0;
    node->coverage[colour]=0;
  }
  hash_table_traverse(&wipe_node, db_graph);
}

void db_graph_print_colour_overlap_matrix(int* first_col_list, int num1,
					  int* second_col_list, int num2,
					  dBGraph* db_graph)
{
  int i;
  int j;
  
  printf("\t");
  for (i=0; i<num1; i++)
    {
      printf("%d", first_col_list[i]);
      if (i<num1-1)
	{
	  printf("\t");
	}
      else
	{
	  printf("\n");
	}      
    }
  for (j=0; j<num2; j++)
    {
      printf("%d\t", second_col_list[j]);
      for (i=0; i<num1; i++)
	{
	  //local function
	  long long overlap_cols_i_and_j(Element* node)
	  {
	    char str[db_graph->kmer_size];

	    if ( 
		(db_node_is_this_node_in_this_person_or_populations_graph(node, individual_edge_array, first_col_list[i])==true)
		&&
		(db_node_is_this_node_in_this_person_or_populations_graph(node, individual_edge_array, second_col_list[j])==true)
		 )
	      {

		return 1;
	      }
	    else
	      {
		return 0;
	      }
	  }

	  long long  overlap = hash_table_traverse_returning_sum(&overlap_cols_i_and_j, db_graph);
	  printf("%qd", overlap);
	  if (i<num1-1)
	    {
	      printf("\t");
	    }
	  else
	    {
	      printf("\n");
	    }
	}
    }
}
