#include <CUnit.h>
#include <Basic.h>
#include <dB_graph.h>
#include <element.h>
#include <binary_kmer.h>
#include <file_reader.h>
#include "test_dB_graph_population.h"
#include "dB_graph_population.h"
#include <global.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>


void test_getting_stats_of_how_many_indivduals_share_a_node()
{ 


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

  long long bad_reads=0;
  long long total_kmers=0;
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.txt", &total_kmers, &bad_reads, hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 44);
  CU_ASSERT(bad_reads==0);

  int stats[6]; //array of ints. stats[0] - number of kmers shared by noone; stats[1], number shared by one person, etc
  stats[0]=0; stats[1]=0; stats[2]=0; stats[3]=0; stats[4]=0; stats[5]=0; 

  int* array[3];
  array[0]=&stats[0];
  array[1]=&stats[1];
  array[2]=&stats[2];

  int num_people=2;

  db_graph_traverse_to_gather_statistics_about_people(&find_out_how_many_individuals_share_this_node_and_add_to_statistics , hash_table, array, num_people);

  CU_ASSERT(stats[0]==0);
  CU_ASSERT(stats[1]==6);
  CU_ASSERT(stats[2]==1);

  hash_table_free(&hash_table);
 
}

