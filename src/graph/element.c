/*
 element.c -- implements the nodes of the dBruijn graph
 */

#include <element.h>
#include <stdlib.h>
#include <global.h>
#include <stdio.h>
#include <binary_kmer.h>

boolean element_is_key(Key key, Element e, short kmer_size)
{
  return key == e.kmer;
}

boolean element_smaller(Element  e1, Element e2){
	return e1.edges < e2.edges;
}


BinaryKmer element_get_kmer(Element * e){
  return e->kmer;
}

Key element_get_key(BinaryKmer kmer, short kmer_size){
  
  BinaryKmer rev_kmer = binary_kmer_reverse_complement(kmer,kmer_size);
  
  if (rev_kmer < kmer){
    kmer = rev_kmer;
  }

  return kmer;

}

void element_initialise(Element * e, Key kmer, short kmer_size){

  e->kmer = element_get_key(kmer, kmer_size);
  e->edges = 0;
  e->status = none;
}


Orientation db_node_get_orientation(BinaryKmer k, dBNode * e, short kmer_size){

  if (e->kmer == k){
    return forward;
  }
  
  if (e->kmer == binary_kmer_reverse_complement(k,kmer_size)){
    return reverse;
  }

  exit(1);
  
}


//add one edge to element -- basically sets a bit in the edges char
void db_node_add_labeled_edge(dBNode * e, Orientation o, Nucleotide base){

  //set edge 
  char edge = 1 << base; // A (0) -> 0001, C (1) -> 0010, G (2) -> 0100, T (3) -> 1000
  
  if (o == reverse){
    edge <<= 4; //move to next nibble 
  }

  //update node
  e->edges |= edge;
  
}


//adding an edge between two nodes implies adding two labeled edges (one in each direction)
//be aware that in the case of self-loops in palindromes the two labeled edges collapse in one

boolean db_node_add_edge(dBNode * src_e, dBNode * tgt_e, Orientation src_o, Orientation tgt_o, short kmer_size){

  BinaryKmer src_k, tgt_k; 

  src_k = src_e->kmer;
  tgt_k = tgt_e->kmer;
 
  if (src_o == reverse){
    src_k = binary_kmer_reverse_complement(src_k,kmer_size);
  }
    
  if (tgt_o == reverse){
    tgt_k = binary_kmer_reverse_complement(tgt_k,kmer_size);
  }
    
  
  if (DEBUG){
    printf("add edge %s -%c-> %s\n",binary_kmer_to_seq(src_k,kmer_size),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(tgt_k)),binary_kmer_to_seq(tgt_k,kmer_size));
  }

  db_node_add_labeled_edge(src_e,src_o,binary_kmer_get_last_nucleotide(tgt_k));

  if (DEBUG){
    printf("add edge %s -%c-> %s\n",binary_kmer_to_seq(tgt_k,kmer_size),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(src_k,kmer_size))),binary_kmer_to_seq(src_k,kmer_size));
  }

  db_node_add_labeled_edge(tgt_e,opposite_orientation(tgt_o),binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(src_k,kmer_size)));

  return true;
}



boolean db_node_edge_exist(dBNode * element,Nucleotide base,Orientation orientation){
  char edge = element->edges;
  edge >>= base;
  if (orientation == reverse){
    edge >>= 4;
  }
  
  edge &= 1;
  
  if (edge == 1){
    return true;
  }
  else{
    return false;
  }
}


boolean db_node_is_supernode_end(dBNode * element,Orientation orientation)
{
  char edges = element->edges;
  
  if (orientation == reverse)
    {
      //shift along so the 4 most significant bits become the 4 least - we've passed out argument by copy so not altering the original
      edges >>= 4;
    }

  edges &= 15; // AND with 00001111 so that we only look at the 4 least significant bits 
   
  //is supernode end EITHER if it has 0 edges out, or >1.
  return ((edges == 0) || ((edges != 1) && (edges != 2) && (edges != 4) && (edges != 8)));

}


Orientation opposite_orientation(Orientation o){
  return o ^ 1;
  
}


void db_node_reset_edges(dBNode * node){
  node->edges = 0;
}

void db_node_reset_edge(dBNode * node, Orientation orientation, Nucleotide nucleotide){

  char edge = 1 << nucleotide;

  if (orientation == reverse){
    edge <<= 4;
  }
  
  //toggle 1->0 0->1

  edge ^= (unsigned char) 0xFF; //xor with all 1's, ie 00010000 -> 11101111

  node->edges &= edge; //reset one edge
  
}




boolean db_node_edges_reset(dBNode * node){
  return node->edges == 0;
}


boolean db_node_has_precisely_one_edge(dBNode * node, Orientation orientation, Nucleotide * nucleotide){
  
  Nucleotide n;
  Edges edges = node->edges;
  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
 
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      *nucleotide = n;
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 1);
  
}

boolean db_node_is_blunt_end(dBNode * node, Orientation orientation){
  
  Edges edges = node->edges;


  if (orientation == reverse){
    edges >>= 4;
  }
  
  edges &= 15; // AND with 00001111 so that we only look at the 4 least significant bits
  
  return edges == 0;
}

boolean db_node_check_status(dBNode * node, NodeStatus status){
  return node->status == status;
}

void db_node_set_status(dBNode * node,NodeStatus status){
  node->status = status;
}

