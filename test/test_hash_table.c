#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <element.h>
#include <hash_table.h>
#include <stdlib.h>
#include "test_hash_table.h"



void test_hash_table_find()
{
  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_buckets=5;
  HashTable* hash_table = hash_table_new(number_of_buckets,kmer_size);

  //Load the following fasta:
  //    >read1
  //    AAAAAAAA
  //    >read2
  //    GGCT
  //    >read3
  //    TAGG

  load_fasta_data_from_filename_into_graph("../test/data/basic/test_hash_table.fasta", hash_table);

  //all the kmers and their reverse complements from the reads
  Element* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size),kmer_size) ,hash_table);
  Element* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size),kmer_size) ,hash_table);
  Element* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size),kmer_size) ,hash_table);
  Element* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size),kmer_size) ,hash_table);
  Element* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size),kmer_size) ,hash_table);
  Element* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size),kmer_size) ,hash_table);
  Element* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size),kmer_size) ,hash_table);
  Element* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size),kmer_size) ,hash_table);
  Element* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size),kmer_size) ,hash_table);
  Element* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size),kmer_size) ,hash_table);

  //kmers that should not be in the graph

  Element* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size),kmer_size) ,hash_table);
  Element* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size),kmer_size) ,hash_table);
  Element* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size),kmer_size) ,hash_table);
  Element* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size),kmer_size) ,hash_table);
  Element* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size),kmer_size) ,hash_table);
  Element* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size),kmer_size) ,hash_table);
  Element* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size),kmer_size) ,hash_table);
  Element* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size),kmer_size) ,hash_table);
  Element* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size),kmer_size) ,hash_table);
  Element* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size),kmer_size) ,hash_table);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element3 != NULL);
  CU_ASSERT(test_element4 != NULL);
  CU_ASSERT(test_element5 != NULL);
  CU_ASSERT(test_element6 != NULL);
  CU_ASSERT(test_element7 != NULL);
  CU_ASSERT(test_element8 != NULL);
  CU_ASSERT(test_element9 != NULL);
  CU_ASSERT(test_element10 != NULL);

  CU_ASSERT(test_element11 == NULL);
  CU_ASSERT(test_element12 == NULL);
  CU_ASSERT(test_element13 == NULL);
  CU_ASSERT(test_element14 == NULL);
  CU_ASSERT(test_element15 == NULL);
  CU_ASSERT(test_element16 == NULL);
  CU_ASSERT(test_element17 == NULL);
  CU_ASSERT(test_element18 == NULL);
  CU_ASSERT(test_element19 == NULL);
  CU_ASSERT(test_element20 == NULL);

  
}
