/*
  dB_graph.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>
#include <binary_kmer.h>
#include <element.h>
#include <dB_graph.h>
#include <dB_graph_population.h>

#include <string.h>


//This gets the next node in the graph, and does not care about whether it
//is there for any specific person or population
//it just looks to see if it is there at all, and if so, gets it
// Also does not care if there is an edge joining this node to the one it returns.
// You simply tell it to add a certain nucleotide, see if there is such a node, and if so, return it.

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



void db_graph_set_all_visited_nodes_to_status_none(dBGraph* hash_table)
{
  hash_table_traverse(&db_node_set_status_to_none, hash_table);
}














