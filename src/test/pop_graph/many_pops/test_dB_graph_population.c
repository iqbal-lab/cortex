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

void test_is_supernode_end()
{

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

  long long count_kmers=0;
  long long bad_reads=0;

  int seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_two_self_loops", &count_kmers, &bad_reads, hash_table);
  CU_ASSERT(seq_loaded==6);
  CU_ASSERT(bad_reads==0);

  //GTA is not the end of a supernode in either direction
  dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("GTA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));

  //CGT is not the end of a supernode in either direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("CGT",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));

  //ACG is not  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));

  hash_table_free(&hash_table);


  // ****
  //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two outward/exiting edges 
  // ****



  //first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;

  count_kmers=0;
  bad_reads=0;

  hash_table = hash_table_new(number_of_buckets,kmer_size);
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_one_long_supernode_with_conflict_at_end", &count_kmers, &bad_reads, hash_table);

  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(bad_reads==0);

  //ACA IS  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));
  

  //ATG is not a supernode end
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));


  //ATT IS  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATT",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));
  
  //TTC is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTC",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));

  //TTG is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));


  hash_table_free(&hash_table);


  // ****
  //1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two INWARD edges in the opposite direction
  // ****

//first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);

  count_kmers=0;
  bad_reads=0;

  seq_loaded=load_population_as_fasta("../data/test/pop_graph/supernode/one_person_one_long_supernode_with_inward_conflict_at_end",&count_kmers, &bad_reads, hash_table);

  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(bad_reads==0);

  //ACA IS  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));
  
  //ATG is not a supernode end
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));


  //ATT IS  the end of a supernode in the reverse direction
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ATT",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));
  
  //TTC is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTC",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));

  //AAT is a supernode end in reverse directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAT",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));

  //TAA is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TAA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));


  hash_table_free(&hash_table);



  // ****
  //1.4 Fasta file that generate a graph with an infinite loop at a single kmer
  //   
  // ****

  
//first set up the hash/graph
  kmer_size = 3;
  number_of_buckets=5;
  hash_table = hash_table_new(number_of_buckets,kmer_size);

  count_kmers=0;
  bad_reads=0;
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_infiniteloop", &count_kmers, &bad_reads, hash_table);
  CU_ASSERT(bad_reads==0);
  CU_ASSERT(seq_loaded==25);

  //AAA is a supernode end in both directions
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_is_supernode_end(query_node, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node, reverse, individual_edge_array, 0, hash_table));


  hash_table_free(&hash_table);


}


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

