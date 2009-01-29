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
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 44);

  // Now just see if it can correctly find the first node in a supernode


  // *******
  //first check all the kmers in person 2. all of them should only be found in person2, except ACG
  // ****


  // **** GGG first

  dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GGG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  dBNode* testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  
  //this kmer doesn't exist in person 1.
  CU_ASSERT(testnode==NULL);

  //try again for person 2:
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  char* answer = binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "GGG") || !strcmp(answer, "CCC") );


  // *** then TTG
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  
  //this kmer doesn't exist in person 1.
  CU_ASSERT(testnode==NULL);

  //try again for person 2:
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  answer = binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "TTG") || !strcmp(answer, "CAA") || !strcmp(answer, "CGT") || !strcmp(answer, "ACG") );
  

  // *****then TGA
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  
  //this kmer doesn't exist in person 1.
  CU_ASSERT(testnode==NULL);

  //try again for person 2:
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  answer = binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "TTG") || !strcmp(answer, "CAA") || !strcmp(answer, "CGT") || !strcmp(answer, "ACG") );


  //*** then GAC
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GAC",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT((query_node!=NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  
  //this kmer doesn't exist in person 1.
  CU_ASSERT(testnode==NULL);

  //try again for person 2:
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode != NULL);  
  answer = binary_kmer_to_seq(testnode ->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "TTG") || !strcmp(answer, "CAA") || !strcmp(answer, "CGT") || !strcmp(answer, "ACG") );

  
  // then ACG - this is the one that I expect to see in both person1 and person2

  //check person1
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));
  answer = binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "ACG") || !strcmp(answer, "AAC")  || !strcmp(answer, "CGT") || !strcmp(answer, "GTT") );

  //then person2
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(!(testnode==NULL));
  answer = binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "ACG") || !strcmp(answer, "CGT")  || !strcmp(answer, "TTG") || !strcmp(answer, "CAA") );


  //*****
  // then do the kmers you expect to find only in person1 (not person 2) (horrible having index from 0 to 1 and people labelled 1 and 2 - stupid of me)


  // *** AAA
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));
  answer = binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "AAA") || !strcmp(answer, "TTT") );

  //confirm not seen in person1
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode==NULL);



  // *** GTT

  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GTT",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 0, hash_table);
  CU_ASSERT(!(testnode==NULL));
  answer = binary_kmer_to_seq(testnode->kmer, hash_table->kmer_size);
  CU_ASSERT( !strcmp(answer, "ACG") || !strcmp(answer, "CGT")  ||  !strcmp(answer, "GTT") || !strcmp(answer, "AAC") );

  //confirm not seen in person1
  testnode = db_graph_get_first_node_in_supernode_containing_given_node_for_specific_person_or_pop(query_node, individual_edge_array, 1, hash_table);
  CU_ASSERT(testnode==NULL);
    
  


}



void test_correctly_find_subsection_of_supernode()
{
}

void test_find_best_subsection_of_supernode()
{
}

void test_get_population_consensus_supernode()
{
}
