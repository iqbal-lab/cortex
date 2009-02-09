#include <CUnit.h>
#include <Basic.h>
#include <dB_graph.h>
#include <element.h>
#include <binary_kmer.h>
//#include <graph.h>
#include <file_reader.h>
#include "test_dB_graph.h"
#include <global.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

void test_supernode_walking()
{

  ///first with the non efficient loading functions

  //get output file 
  FILE* output = fopen("../test/test_dB_graph_dir/test_log.txt","w");

  int kmer_size;
  int number_of_buckets;

  //first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  HashTable* hash_table = hash_table_new(number_of_buckets,kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a graph/hash table.
  //         The test then picks an element in the graph, and calls get_seq_from_elem_to_end_of_supernode
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****


  load_fasta_data_from_filename_into_graph("../test/test_dB_graph_dir/generates_graph_with_two_self_loops.fasta", hash_table);

 
  // now start at GTA and get all the sequence from there to the end of the supernode, and see
  // if that is right. Remember that once encoded GTA, might be in a kmer where to read GTA you must go
  // in the reverse orientation. 
  
  Element* test_element1 = hash_table_find(seq_to_binary_kmer("GTA", kmer_size) ,hash_table);

  CU_ASSERT(!(test_element1==NULL));


  boolean is_cycle=false;  

  char* test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table, &is_cycle);
  CU_ASSERT(is_cycle);

  //printf("Sequence after GTA and up to end of supernode is  %s\n",test1_result);
  //printf("the is_cycle variable is %d\n",is_cycle); 

  //depending on what orientation GTA is with respect to it's kmer, the answer you get will be either CGT or GT.
  boolean pass = ((!strcmp(test1_result,"CGT")) || (!strcmp(test1_result,"GT")) );
  //printf("test 1.1 result returned by get_seq_from_elem_to_end_of_supernode is %s",test1_result);
  CU_ASSERT(pass);


  hash_table_free(&hash_table);
  free(test1_result);

 
  // ****
  //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two outward/exiting edges 
  // ****



  //first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);
  load_fasta_data_from_filename_into_graph("../test/test_dB_graph_dir/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta", hash_table);
  
  test_element1 = hash_table_find(seq_to_binary_kmer("ACA", kmer_size) ,hash_table);

  is_cycle=false;

  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table, &is_cycle);
  CU_ASSERT(!is_cycle);

  //printf("Sequence after ACA and up to end of supernode is  %s\n",test1_result);
  //printf("the is_cycle variable is %d\n",is_cycle);

  //depending on what orientation ACA is with respect to it's kmer, the answer you get will be either TT or nothing
  pass = ((!strcmp(test1_result,"TT")) || (!strcmp(test1_result,"")) );
  //  printf("this time result is %s",test1_result);
  CU_ASSERT(pass);

  hash_table_free(&hash_table);
  free(test1_result);




  // ****
  //1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two INWARD edges in the opposite direction
  // ****

//first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);
  load_fasta_data_from_filename_into_graph("../test/test_dB_graph_dir/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta", hash_table);
  
  test_element1 = hash_table_find(seq_to_binary_kmer("ACA", kmer_size) ,hash_table);

  is_cycle=false;

  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table, &is_cycle);
  CU_ASSERT(!is_cycle);

  //printf("Sequence after ACA and up to end of supernode is  %s\n",test1_result);
  //printf("the is_cycle variable is %d\n",is_cycle);

  //depending on what orientation ACA is with respect to it's kmer, the answer you get will be either TT or nothing
  pass = ((!strcmp(test1_result,"TT")) || (!strcmp(test1_result,"")) );
  //  printf("this time result is %s",test1_result);
  CU_ASSERT(pass);

  hash_table_free(&hash_table);
  free(test1_result);



  // ****
  //1.4 Fasta file that generate a graph with an infinite loop at a single kmer
  //   
  // ****

  
//first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);
  load_fasta_data_from_filename_into_graph("../test/test_dB_graph_dir/generates_graph_with_infinite_loop.fasta", hash_table);
  
  test_element1 = hash_table_find(seq_to_binary_kmer("AAA", kmer_size) ,hash_table);

  is_cycle=false;

  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table, &is_cycle);
  CU_ASSERT(is_cycle);//this time it should find a cycle

  //printf("Sequence after AAA and up to end of supernode is  %s\n",test1_result);
  //printf("the is_cycle variable is %d\n",is_cycle);

  pass = (!strcmp(test1_result,""));
  //printf("this time result is %s and pass variable is %d",test1_result, pass);
  CU_ASSERT(pass);

  hash_table_free(&hash_table);
  free(test1_result);





  /// then with the efficient loading functions



  // Same Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a graph/hash table.
  //         The test then picks an element in the graph, and calls get_seq_from_elem_to_end_of_supernode
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****

  //first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);
  int num_reallocs=0;
  load_fasta_data_from_filename_into_graph_efficient("../test/test_dB_graph_dir/generates_graph_with_two_self_loops.fasta", hash_table, &num_reallocs);
  
  CU_ASSERT(num_reallocs==0);


  // now start at GTA and get all the sequence from there to the end of the supernode, and see
  // if that is right. Remember that once encoded GTA, might be in a kmer where to read GTA you must go
  // in the reverse orientation. 
  
  test_element1 = hash_table_find(seq_to_binary_kmer("GTA", kmer_size) ,hash_table);

  CU_ASSERT(!(test_element1==NULL));

  is_cycle=false;  
  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table, &is_cycle);
  CU_ASSERT(is_cycle);

  //printf("Sequence after GTA and up to end of supernode is  %s\n",test1_result);
  //printf("the is_cycle variable is %d\n",is_cycle); 

  //depending on what orientation GTA is with respect to it's kmer, the answer you get will be either CGT or GT.
  pass = ((!strcmp(test1_result,"CGT")) || (!strcmp(test1_result,"GT")) );
  //printf("test 1.1 result returned by get_seq_from_elem_to_end_of_supernode is %s",test1_result);
  CU_ASSERT(pass);


  hash_table_free(&hash_table);
  free(test1_result);

  // ****
  //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two outward/exiting edges 
  // ****



  //first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);
  num_reallocs=0;
  load_fasta_data_from_filename_into_graph_efficient("../test/test_dB_graph_dir/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta", hash_table, &num_reallocs);

  CU_ASSERT(num_reallocs==0);



  test_element1 = hash_table_find(seq_to_binary_kmer("ACA", kmer_size) ,hash_table);

  is_cycle=false;

  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table, &is_cycle);
  CU_ASSERT(!is_cycle);

  //printf("Sequence after ACA and up to end of supernode is  %s\n",test1_result);
  //printf("the is_cycle variable is %d\n",is_cycle);

  //depending on what orientation ACA is with respect to it's kmer, the answer you get will be either TT or nothing
  pass = ((!strcmp(test1_result,"TT")) || (!strcmp(test1_result,"")) );
  //  printf("this time result is %s",test1_result);
  CU_ASSERT(pass);

  hash_table_free(&hash_table);
  free(test1_result);




  // ****
  //1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two INWARD edges in the opposite direction
  // ****

//first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);
  num_reallocs=0;
  load_fasta_data_from_filename_into_graph_efficient("../test/test_dB_graph_dir/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta", hash_table, &num_reallocs);
  CU_ASSERT(num_reallocs==0);

  
  test_element1 = hash_table_find(seq_to_binary_kmer("ACA", kmer_size) ,hash_table);

  is_cycle=false;

  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table, &is_cycle);
  CU_ASSERT(!is_cycle);

  //printf("Sequence after ACA and up to end of supernode is  %s\n",test1_result);
  //printf("the is_cycle variable is %d\n",is_cycle);

  //depending on what orientation ACA is with respect to it's kmer, the answer you get will be either TT or nothing
  pass = ((!strcmp(test1_result,"TT")) || (!strcmp(test1_result,"")) );
  //  printf("this time result is %s",test1_result);
  CU_ASSERT(pass);

  hash_table_free(&hash_table);
  free(test1_result);



  // ****
  //1.4 Fasta file that generate a graph with an infinite loop at a single kmer
  //   
  // ****

  
//first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);
  num_reallocs=0;
  load_fasta_data_from_filename_into_graph_efficient("../test/test_dB_graph_dir/generates_graph_with_infinite_loop.fasta", hash_table, &num_reallocs);
  CU_ASSERT(num_reallocs==0);

  
  test_element1 = hash_table_find(seq_to_binary_kmer("AAA", kmer_size) ,hash_table);

  is_cycle=false;

  test1_result = get_seq_from_elem_to_end_of_supernode(test_element1, forward, hash_table, &is_cycle);
  CU_ASSERT(is_cycle);//this time it should find a cycle

  //printf("Sequence after AAA and up to end of supernode is  %s\n",test1_result);
  //printf("the is_cycle variable is %d\n",is_cycle);

  pass = (!strcmp(test1_result,""));
  //printf("this time result is %s and pass variable is %d",test1_result, pass);
  CU_ASSERT(pass);

  hash_table_free(&hash_table);
  free(test1_result);

  fclose(output);

}		    


