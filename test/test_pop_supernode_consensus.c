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
    
  
  hash_table_free(&hash_table);

}



void test_correctly_find_subsection_of_supernode()
{

  int kmer_size = 5;
  int number_of_buckets=8;
  HashTable* hash_table = hash_table_new(number_of_buckets,kmer_size);
  
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  int seq_loaded = load_population_as_fasta("../test/data/pop/two_people_test_consensus.txt", hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 23);

  
  //have just loaded the following
  // >cons_person1_read
  //  AAGACTGGAGAT
  // >cons_person2_read
  // ATCCTGGAGAA




  //notice the people share CTGGA
  //notice that each person basically has one big supernode, so the "breaking of supernodes" only occurs when you compare people.


  char subsection[24];

  dBNode* node = hash_table_find(element_get_key(seq_to_binary_kmer("GGAGA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(node != NULL);

  //test the subsections of person1's supernode, startng from the beginning

  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,1,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACT");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,2,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,3,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,4,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGGA");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,5,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGGAG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,6,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGGAGA");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,7,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "AAGACTGGAGAT");

  //what happens if you go too far?
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,8,individual_edge_array, 0, hash_table)==1);
  //what happens with negative numbers?
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,-1,individual_edge_array, 0, hash_table)==1);
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, -1,8,individual_edge_array, 0, hash_table)==1);
  //what happens if start=end
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,0,individual_edge_array, 0, hash_table)==1);
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 1,1,individual_edge_array, 0, hash_table)==1);
  //what happens if start>end
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 5,4,individual_edge_array, 0, hash_table)==1);

  //now check a few others, not starting from 0
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 3,6,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ACTGGAGA");

  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 2,3,individual_edge_array, 0, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "GACTGG");


  //now do the same kind of thing for the other person

  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,1,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,2,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,3,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGGA");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,4,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGGAG");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,5,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGGAGA");
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,6,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "ATCCTGGAGAA");

  //what happens if you go too far?
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,18,individual_edge_array, 1, hash_table)==1);
  //what happens with negative numbers?
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,-1,individual_edge_array, 1, hash_table)==1);
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, -1,8,individual_edge_array, 1, hash_table)==1);
  //what happens if start=end
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 0,0,individual_edge_array, 1, hash_table)==1);
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 1,1,individual_edge_array, 1, hash_table)==1);
  //what happens if start>end
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 5,4,individual_edge_array, 1, hash_table)==1);
  
  //now check a few others, not starting from 0
  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 1,3,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "TCCTGGA");

  CU_ASSERT(db_graph_get_subsection_of_supernode_containing_given_node_as_sequence(subsection, node, 3,5,individual_edge_array, 1, hash_table)==0);
  CU_ASSERT_STRING_EQUAL(subsection, "CTGGAGA");
  

  hash_table_free(&hash_table);
}

void test_find_best_subsection_of_supernode()
{
}

void test_get_population_consensus_supernode()
{
}
