#include <CUnit.h>
#include <Basic.h>
#include <dB_graph.h>
#include <element.h>
#include <binary_kmer.h>
#include <file_reader.h>
#include <test_dB_graph.h>

void test_hash_table_find()
{
  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_buckets = 4; // has to be power of 2 -- need to fix this
  long long count_kmers = 0;
  long long bad_reads = 0; 
  


  dBGraph * db_graph = hash_table_new(number_of_buckets,kmer_size);
  
  //Load the following fasta:
  //    >read1
  //    AAAAAAAA
  //    >read2
  //    GGCT
  //    >read3
  //    TAGG

  int seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/test_dB_graph.fasta", &count_kmers, &bad_reads, 20, db_graph);
  
  //length of total sequence
  CU_ASSERT(seq_length == 16);
  
  //number of kmers

  CU_ASSERT_EQUAL(count_kmers, 5);

  //bad reads
  CU_ASSERT_EQUAL(bad_reads, 0);

  //all the kmers and their reverse complements from the reads
  Element* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size),kmer_size) ,db_graph);
  Element* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size),kmer_size) ,db_graph);
  Element* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size),kmer_size) ,db_graph);
  Element* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size),kmer_size) ,db_graph);
  Element* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size),kmer_size) ,db_graph);
  Element* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size),kmer_size) ,db_graph);
  Element* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size),kmer_size) ,db_graph);
  Element* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size),kmer_size) ,db_graph);
  Element* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size),kmer_size) ,db_graph);
  Element* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size),kmer_size) ,db_graph);

  //kmers that should not be in the graph

  Element* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size),kmer_size) ,db_graph);
  Element* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size),kmer_size) ,db_graph);
  Element* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size),kmer_size) ,db_graph);
  Element* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size),kmer_size) ,db_graph);
  Element* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size),kmer_size) ,db_graph);
  Element* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size),kmer_size) ,db_graph);
  Element* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size),kmer_size) ,db_graph);
  Element* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size),kmer_size) ,db_graph);
  Element* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size),kmer_size) ,db_graph);
  Element* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size),kmer_size) ,db_graph);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);

  CU_ASSERT(test_element3 != NULL);
  CU_ASSERT(test_element4 != NULL);
  CU_ASSERT(test_element3 == test_element4);

  CU_ASSERT(test_element5 != NULL);
  CU_ASSERT(test_element6 != NULL);  
  CU_ASSERT(test_element5 == test_element6);

  CU_ASSERT(test_element7 != NULL);
  CU_ASSERT(test_element8 != NULL);
  CU_ASSERT(test_element7 == test_element8);

  CU_ASSERT(test_element9 != NULL);
  CU_ASSERT(test_element10 != NULL);
  CU_ASSERT(test_element9 == test_element10);  

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

  
  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
  
}




void test_supernode_walking()
{

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_buckets=1;
  int seq_length;
  long long count_kmers = 0;
  long long bad_reads = 0;

   dBGraph * db_graph = hash_table_new(number_of_buckets,kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a dB_graph.
  //         The test then picks an element in the graph, and calls get_seq_from_elem_to_end_of_supernode
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****


   seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta",&count_kmers, &bad_reads, 20,  db_graph);

   CU_ASSERT_EQUAL(seq_length,6);
   CU_ASSERT_EQUAL(count_kmers,2);
   CU_ASSERT_EQUAL(bad_reads,0);
   

  // now start at GTA and get all the sequence from there to the end of the supernode, and see
  // if that is right. Remember that once encoded GTA, might be in a kmer where to read GTA you must go
  // in the reverse orientation.
  
   Element* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size), kmer_size),db_graph);
   CU_ASSERT(test_element1!=NULL);
   
   
   boolean is_cycle=false;
   
   char test1_result[5000];
   
   get_seq_from_elem_to_end_of_supernode(test_element1, forward, db_graph, &is_cycle, test1_result,5000);
   CU_ASSERT(is_cycle);
   
   //depending on what orientation GTA is with respect to it's kmer, the answer you get will be either CGT or GT
   //Zam we know GTA<TAC!
   
   CU_ASSERT_STRING_EQUAL(test1_result,"CGT");
   
   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);
   
 
   /*   // **** */
   /*   //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end */
   /*   //   caused by two outward/exiting edges */
   /*   // **** */

   //first set up the hash/graph
   kmer_size = 3;
   number_of_buckets=4;
   count_kmers = 0;
   bad_reads = 0;

   db_graph = hash_table_new(number_of_buckets,kmer_size);
   seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta",&count_kmers,&bad_reads,20,db_graph);

   CU_ASSERT_EQUAL(seq_length,13);
   CU_ASSERT_EQUAL(count_kmers,5);
   CU_ASSERT_EQUAL(bad_reads,0);


   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size), kmer_size),db_graph);
   CU_ASSERT(test_element1!=NULL);

   //ACA < TGT so forward gives TT
   is_cycle=false;
   get_seq_from_elem_to_end_of_supernode(test_element1, forward, db_graph, &is_cycle, test1_result, 5000);
   CU_ASSERT(!is_cycle);   
   CU_ASSERT_STRING_EQUAL(test1_result,"TT");

   
   //ACA < TGT so backward gives ""
   is_cycle=false;
   get_seq_from_elem_to_end_of_supernode(test_element1, reverse, db_graph, &is_cycle, test1_result, 5000);
   CU_ASSERT(!is_cycle);   
   CU_ASSERT_STRING_EQUAL(test1_result,"");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);




   // ****
   //1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
   //   caused by two INWARD edges in the opposite direction
   // ****
   
   //first set up the hash/graph
   kmer_size = 3;
   number_of_buckets=8;
   count_kmers = 0;
   bad_reads = 0;

   db_graph = hash_table_new(number_of_buckets,kmer_size);
   
   seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta",&count_kmers, &bad_reads, 20,  db_graph);

   CU_ASSERT_EQUAL(seq_length,13);
   CU_ASSERT_EQUAL(count_kmers,5);
   CU_ASSERT_EQUAL(bad_reads,0);
   

   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size), kmer_size) ,db_graph);
   CU_ASSERT(test_element1!=NULL);

   //ACA < TGT so forward gives TT
   is_cycle=false;
   get_seq_from_elem_to_end_of_supernode(test_element1, forward, db_graph, &is_cycle, test1_result, 5000);
   CU_ASSERT(!is_cycle);   
   CU_ASSERT_STRING_EQUAL(test1_result,"TT");
   
   //ACA < TGT so backwar gives ""
   is_cycle=false;
   get_seq_from_elem_to_end_of_supernode(test_element1, reverse, db_graph, &is_cycle, test1_result, 5000);
   CU_ASSERT(!is_cycle);   
   CU_ASSERT_STRING_EQUAL(test1_result,"");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


   // ****
   //1.4 Fasta file that generate a graph with an infinite loop at a single kmer
   //
   // ****

  
   //first set up the hash/graph
   kmer_size = 3;
   number_of_buckets=32;
   count_kmers = 0;
   bad_reads = 0;

   db_graph = hash_table_new(number_of_buckets,kmer_size);
   seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_infinite_loop.fasta",&count_kmers,&bad_reads,30,db_graph);

   CU_ASSERT_EQUAL(seq_length,25);
   CU_ASSERT_EQUAL(count_kmers,1);
   CU_ASSERT_EQUAL(bad_reads,0);

   test_element1 = hash_table_find(seq_to_binary_kmer("AAA", kmer_size) ,db_graph);
   CU_ASSERT(test_element1!=NULL);

   //forward
   is_cycle=false;
   get_seq_from_elem_to_end_of_supernode(test_element1, forward, db_graph, &is_cycle, test1_result, 5000);
   CU_ASSERT(is_cycle);//this time it should find a cycle
   CU_ASSERT_STRING_EQUAL(test1_result,"");
  

   //backward
   is_cycle=false;
   get_seq_from_elem_to_end_of_supernode(test_element1, reverse, db_graph, &is_cycle, test1_result, 5000);
   CU_ASSERT(is_cycle);//this time it should find a cycle
   CU_ASSERT_STRING_EQUAL(test1_result,"");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


}


void test_writing_reading_graph(){

  dBNode node1, node2;
  int kmer_size = 3;
  boolean succ;

  char seq[kmer_size];

  seq_to_binary_kmer("AAA",3);

  
  element_initialise(&node1,seq_to_binary_kmer("AAA",kmer_size),kmer_size);
  node1.edges = 'a'; // ie 64 = (0110 0100)2


  FILE* fp1 = fopen("../data/test/graph/dump_element.bin", "w");
  db_node_print_binary(fp1,&node1);
  fclose(fp1);

  FILE* fp2 = fopen("../data/test/graph/dump_element.bin", "r");
  
  succ = db_node_read_binary(fp2,kmer_size,&node2);

  CU_ASSERT_EQUAL(succ, true);

  CU_ASSERT_EQUAL(node1.kmer, node2.kmer);
  CU_ASSERT_EQUAL(node1.edges, node2.edges);

  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(node1.kmer,kmer_size,seq));
  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(node2.kmer,kmer_size,seq));
  
  CU_ASSERT_EQUAL(node1.edges, node2.edges);
  CU_ASSERT_EQUAL(node1.edges, 'a');
  CU_ASSERT_EQUAL(node2.edges, 'a');
  
  succ = db_node_read_binary(fp2,kmer_size,&node2);
  CU_ASSERT_EQUAL(succ, false);

  fclose(fp2);
}
   
