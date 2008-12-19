#include <CUnit.h>
#include <Basic.h>
#include <dB_graph.h>
#include <dB_graph_node.h>
#include <binary_kmer.h>
#include <graph.h>
#include "test_dB_graph.h"

#include <stdio.h>
#include <assert.h>

void test_supernode_walking()
{

  //get output file 
  FILE* output = fopen("test_dB_graph_dir/test_log.txt","w");

  int kmer_size;
  int number_of_buckets;

  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a graph/hash table.
  //         The test then picks an element in the graph, and calls get_seq_from_elem_to_end_of_supernode
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two self loops, and a single edge joining them.
  // ****

  //first set up the hash/graph

  kmer_size = 3;
  number_of_buckets=5;
  HashTable* hash_table = new_hash_table(kmer_size,5);

  FILE* fp_test1 = fopen("test_dB_graph_dir/generates_graph_with_two_self_loops.fasta","r");
  load_fasta_data_into_graph(fp_test1, hash_table);

  // now start at GTA and get all the sequence from there to the end of the supernode, and see
  // if that is right.
  
  Element* test_element1;
  char* test_kmer1="GTA";
  test_element1->kmer = seq_to_bin(test_kmer1, kmer_size); //take responsibility for freeing this
  test_element1->edges=0;
  test_element1->visited=0;

  char* test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table);
  
  CU_ASSERT_STRING_EQUAL(test1_result, "GTAC");

  free(test_element1->kmer); 

  //now start at ACG and get all the sequence from there to the end of the supernode and check if that is right
  test_kmer1="ACG";
  test_element1->kmer = seq_to_bin(test_kmer1, kmer_size); //take responsibility for freeing this

  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table);
  CU_ASSERT_STRING_EQUAL(test1_result,"ACGTAC");
  free(test_element1->kmer);

  //now test the opposite orientation
  test_kmer1="ACG";
  test_element1->kmer = seq_to_bin(test_kmer1, kmer_size); //take responsibility for freeing this

  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, reverse, hash_table);
  CU_ASSERT_STRING_EQUAL(test1_result,"GCATGC");
  free(test_element1->kmer);

  



//1. Test realloc when sequence > max. Does realloc preserve the data you already had in there?


  fclose(output);
	
}		    
