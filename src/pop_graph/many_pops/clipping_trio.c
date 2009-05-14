/*
  clipping_trio.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>

#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <pqueue_pop.h>
#include <clipping_trio.h>
#include <seq.h>
#include <string.h>


//internal function, do not export
void db_graph_remove_path_for_specific_person_or_population(dBNode * node, Orientation orientation, int length, dBGraph * db_graph, EdgeArrayType type, int index)
{

  char seq[db_graph->kmer_size];
  Nucleotide nucleotide, rev_nucleotide;
  dBNode * next_node;
  Orientation next_orientation;
  int i;

  for(i = 0; i < length; i++){

    if (db_node_has_precisely_one_edge(node,orientation,&nucleotide, type, index)){

      if (DEBUG){
	printf("CLIPPING %c\n",binary_nucleotide_to_char(nucleotide));
      }
      
      //the only reason we know this next node is connected by en edge, is because this function db_node_has_precisely_one_edge_in_union_graph_over_all_people returned it in argument 3
      next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,nucleotide,&rev_nucleotide,db_graph, type, index); 

      if (DEBUG){
	printf("RESET %c BACK\n",binary_nucleotide_to_char(rev_nucleotide));
      }
 
      
      if(next_node == NULL){
	printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
	exit(1);
      }	           
      
      
      //set status of node to pruned
      db_node_set_status(node, pruned);
      fprintf(stderr, "WARNING. Talk to Mario. Setting node to pruned when only pruned from one persons graph");

      //independent of all that, whichever people the node is now pruned for, the only thing we have just done is make it pruned for one more person - identified by type, index
      db_node_reset_edges(node, type, index);

      db_node_reset_edge(next_node,opposite_orientation(next_orientation),rev_nucleotide, type, index);

    }
    else{
      printf("dB_graph: assume unique path kmer:%s\n",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
      exit(1);
    }
    
    node = next_node;
    orientation = next_orientation;
  }
}

//internal function do not export
int db_graph_detect_tip_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit,dBGraph * db_graph, EdgeArrayType type, int index)
{ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  char seq[db_graph->kmer_size];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end(node, opposite_orientation(orientation), type, index)){
    
    boolean join_main_trunk = false;

     while(db_node_has_precisely_one_edge(node,orientation,&nucleotide,type, index)) {
       Orientation next_orientation;
       dBNode * next_node;
       
       next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph, type, index);
       
       if(next_node == NULL){
	 printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
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

       
       node = next_node;
       orientation = next_orientation;
      
    
     }
  
  
     if (! join_main_trunk){
       length = 0;
     }

  }

  return length;
}

//internal function, do not export
int db_graph_clip_tip_with_orientation_for_specific_person_or_pop(dBNode * node, Orientation orientation, int limit,dBGraph * db_graph, EdgeArrayType type, int index){ 

  Nucleotide nucleotide, reverse_nucleotide;
  int length = 0;
  int i;
  dBNode * nodes[limit];
  Orientation next_orientation;
  dBNode * next_node;
  char seq[db_graph->kmer_size];
  
  //starting in a blunt end also prevents full loops 
  if (db_node_is_blunt_end(node, opposite_orientation(orientation), type, index)){
    
    boolean join_main_trunk = false;

     while(db_node_has_precisely_one_edge(node,orientation,&nucleotide, type, index))
       {
  
	 nodes[length] = node;
	 
	 next_node = db_graph_get_next_node_for_specific_person_or_pop(node,orientation,&next_orientation,nucleotide,&reverse_nucleotide,db_graph, type, index);
	 
	 if(next_node == NULL){
	   printf("dB_graph: didnt find node in hash table: %s\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
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

	 db_node_set_status(nodes[i],pruned);
	 fprintf(stderr, "WARNING. Talk to Mario. Setting node to pruned when only pruned from one persons graph");
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




int db_graph_clip_tip_for_specific_person_or_pop(dBNode * node, int limit,dBGraph * db_graph, EdgeArrayType type, int index){

  int length_tip = 0;
  char seq[db_graph->kmer_size];

  if (DEBUG){
    printf("START Node %s\n",binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
  }

  if (db_node_check_status(node,none) || db_node_check_status(node,visited) )
    {

      length_tip = db_graph_clip_tip_with_orientation_for_specific_person_or_pop(node,forward,limit,db_graph, type, index);
      
      if (length_tip!=0){
	if (DEBUG){
	  printf("\tforward tip: %i\n",length_tip);
	}
      }
      else{
	length_tip = db_graph_clip_tip_with_orientation_for_specific_person_or_pop(node,reverse,limit,db_graph, type, index);
	
	if (length_tip!=0){
	  if (DEBUG){
	    printf("\treverse tip: %i\n",length_tip);
	  }
	}     
      }
    }
  else{
    if (DEBUG){
      if ( !(   db_node_check_status(node,none) || db_node_check_status(node,visited)  ) )
      {
	printf("\n%s: pruned\n", binary_kmer_to_seq(element_get_kmer(node),db_graph->kmer_size,seq));
      }
    }
  }
  
  return length_tip;
}
