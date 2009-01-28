#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <dB_graph_population.h>
#include <dB_graph_supernode.h>
#include <element.h>
#include <hash_table.h>
#include <stdlib.h>
#include "supernode_cmp.h"


void test_find_first_node_in_supernode()
{


  //*********************************************************************************
  // 1. Test two simple fasta, one each for a different person in the same pop
  //
  //       >person1_read1
  //        AAAAAAAAAAAAAAA
  //        >person1_read2
  //        ACGTT
  //
  //
  //       >person2_read1
  //       GGGGGGGGGGGGGGGGGG
  //       >person2_read2
  //       TTGACG


  int kmer_size = 3;
  int number_of_buckets=5;
  HashTable* hash_table = hash_table_new(number_of_buckets,kmer_size);
  
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  int seq_loaded = load_population_as_fasta("../test/data/test_pop_load_and_print/two_individuals_simple.txt", hash_table);
  printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 44);

  // Now just see if it can correctly find the forst node in a supernode
  dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GGG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  if (query_node==NULL)
    {
      printf("query_node is NULL> exit");
      exit(1);
    }


  CU_ASSERT((query_node==NULL));

  dBNode* testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));
  char* answer = binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "GGG") || !strcmp(answer, "CCC") );


  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));

  answer = binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size);

  CU_ASSERT( !strcmp(answer, "ACG") || !strcmp(answer, "AAC") );



  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));

  answer = binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size);

  CU_ASSERT( !strcmp(answer, "AAA") || !strcmp(answer, "TTT") );




    
  


}
