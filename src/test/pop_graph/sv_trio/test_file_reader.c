#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <stdlib.h>
#include <test_file_reader.h>

void test_dump_load_sv_trio_binary(){

  int kmer_size = 3;
  int number_of_bits_pre = 4; 
  int number_of_bits_post = 8;
  int bucket_size = 5;
  long long bad_reads = 0; 
  FILE * fout = fopen("../data/test/pop_graph/dump_sv_trio_graph.bin", "w");
  int seq_length_pre,seq_length_post;
  dBGraph * db_graph_pre;
  dBGraph * db_graph_post;
  
  

  void print_node_binary(dBNode * node){
    db_node_print_binary(fout,node);
  }

  db_graph_pre = hash_table_new(number_of_bits_pre,bucket_size,10,kmer_size);
  seq_length_pre = load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/test_dB_graph.fasta", &bad_reads, 20, db_graph_pre, individual_edge_array,0);

  if (fout == NULL){
    fprintf(stderr,"cannot open ../data/test/pop_graph/dump_sv_trio_graph.bin");
    exit(1);
  }

  hash_table_traverse(&print_node_binary,db_graph_pre);
  
  hash_table_free(&db_graph_pre);
  CU_ASSERT(db_graph_pre==NULL);

  fclose(fout);

  db_graph_post = hash_table_new(number_of_bits_post,bucket_size,10,kmer_size);
  seq_length_post = load_sv_trio_binary_data_from_filename_into_graph("../data/test/pop_graph/dump_sv_trio_graph.bin", db_graph_post);

  printf("***************** I don't understand this\n");
  CU_ASSERT_EQUAL(seq_length_post,15);
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),5);

  //all the kmers and their reverse complements from the reads
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size),kmer_size) ,db_graph_post);


  //kmers that should not be in the graph
  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size),kmer_size) ,db_graph_post);
  dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size),kmer_size) ,db_graph_post);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);
  CU_ASSERT(db_node_get_coverage(test_element1,individual_edge_array,0)==1);

  CU_ASSERT(test_element3 != NULL);
  CU_ASSERT(test_element4 != NULL);
  CU_ASSERT(test_element3 == test_element4);
  CU_ASSERT(db_node_get_coverage(test_element3,individual_edge_array,0)==1);

  CU_ASSERT(test_element5 != NULL);
  CU_ASSERT(test_element6 != NULL);  
  CU_ASSERT(test_element5 == test_element6);
  CU_ASSERT(db_node_get_coverage(test_element5,individual_edge_array,0)==1);

  CU_ASSERT(test_element7 != NULL);
  CU_ASSERT(test_element8 != NULL);
  CU_ASSERT(test_element7 == test_element8);
  CU_ASSERT(db_node_get_coverage(test_element7,individual_edge_array,0)==1);

  CU_ASSERT(test_element9 != NULL);
  CU_ASSERT(test_element10 != NULL);
  CU_ASSERT(test_element9 == test_element10);  
  CU_ASSERT(db_node_get_coverage(test_element9,individual_edge_array,0)==1);

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


  //check arrows
  
  Nucleotide base;

  // AAA -A-> AAA
  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, individual_edge_array, 0)==true);
  CU_ASSERT_EQUAL(base,Adenine);

  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, individual_edge_array, 0)==true);
  CU_ASSERT_EQUAL(base,Thymine);

  //GGC -T-> GCT
  CU_ASSERT(db_node_has_precisely_one_edge(test_element3, reverse,&base, individual_edge_array, 0)==true);
  CU_ASSERT_EQUAL(base,Thymine);


  //TAG -G-> AGG
  CU_ASSERT(db_node_has_precisely_one_edge(test_element7, reverse,&base, individual_edge_array, 0)==true);
  CU_ASSERT_EQUAL(base,Guanine);

  //AGC -C-> GCC
  CU_ASSERT(db_node_has_precisely_one_edge(test_element6, forward,&base, individual_edge_array, 0)==true);
  CU_ASSERT_EQUAL(base,Cytosine);

  //CCT -A-> CTA
  CU_ASSERT(db_node_has_precisely_one_edge(test_element10, reverse,&base, individual_edge_array, 0)==true);
  CU_ASSERT_EQUAL(base,Adenine);

  //add one extra arrows by had -- this breaks the graph - it is only to check that the arrows get set correctly

  add_edges(test_element1,individual_edge_array,0,0x02);
  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, individual_edge_array, 0)==false);


  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, individual_edge_array, 0)==true);
  CU_ASSERT_EQUAL(base,Thymine);

  CU_ASSERT(db_node_edge_exist(test_element1, Adenine, forward, individual_edge_array, 0)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, forward, individual_edge_array, 0)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Guanine, forward, individual_edge_array, 0)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Thymine, forward, individual_edge_array, 0)==false);
  


  add_edges(test_element1, individual_edge_array,0,0x20);
  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, individual_edge_array, 0)==false);

  CU_ASSERT(db_node_edge_exist(test_element1, Adenine, forward, individual_edge_array, 0)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, forward, individual_edge_array, 0)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Guanine, forward, individual_edge_array, 0)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Thymine, forward, individual_edge_array, 0)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Adenine, reverse, individual_edge_array, 0)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, reverse, individual_edge_array, 0)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Guanine, reverse, individual_edge_array, 0)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Thymine, reverse, individual_edge_array, 0)==true);
  

  hash_table_free(&db_graph_post);
  CU_ASSERT(db_graph_post == NULL);



  ////************************************************************************************************
  // Now try a more complicated example. Load fasta for each member of the trio
  // Load chromosome overlap data also.
  // Dump binary, reload, and compare results.


  


  
}



void test_load_graph_binary()
{

  int kmer_size = 5;
  int number_of_bits = 7;
  int bucket_size = 5;
  long long bad_reads = 0;
  int seq_len=0;
  int max_retries=10;
  
  dBGraph* db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  //pre-prepared binary file built from this fasta
  // 
  seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_to_load_and_dump_as_graph_bin_kmer5.bin", db_graph, individual_edge_array, 0);

  CU_ASSERT(seq_len==6);//kmers loaded
  printf("Seq length is %d", seq_len);
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 6);
  printf("Has table says %lld", hash_table_get_unique_kmers(db_graph));


  //all the nodes and their rev complements from the grap
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AACGT",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACGTT",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGTTC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GAACG",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTCC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("GGAAC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTCA",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("TGAAC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("TTTTT",  kmer_size),kmer_size) ,db_graph);

  // nodes that should not be in the graph
  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("ATATA",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("TGGGG",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("AATAG",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("CTCTC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("GGCGG",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGGA",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TACTA",  kmer_size),kmer_size) ,db_graph);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);
  CU_ASSERT(db_node_get_coverage(test_element1,individual_edge_array,0)==4);
  printf("Warning - graph has a different defn of covg to sv_trio. It increments each time you see a kmer. sv_trio only incrmeents once per read\n");

  CU_ASSERT(test_element3 != NULL);
  CU_ASSERT(test_element4 != NULL);
  CU_ASSERT(test_element3 == test_element4);
  CU_ASSERT(db_node_get_coverage(test_element3,individual_edge_array,0)==2);

  CU_ASSERT(test_element5 != NULL);
  CU_ASSERT(test_element6 != NULL);
  CU_ASSERT(test_element5 == test_element6);
  CU_ASSERT(db_node_get_coverage(test_element5,individual_edge_array,0)==2);

  CU_ASSERT(test_element7 != NULL);
  CU_ASSERT(test_element8 != NULL);
  CU_ASSERT(test_element7 == test_element8);
  CU_ASSERT(db_node_get_coverage(test_element7,individual_edge_array,0)==1);

  CU_ASSERT(test_element9 != NULL);
  CU_ASSERT(test_element10 != NULL);
  CU_ASSERT(test_element9 == test_element10);
  CU_ASSERT(db_node_get_coverage(test_element10,individual_edge_array,0)==254);

  CU_ASSERT(test_element11 == NULL);
  CU_ASSERT(test_element12 == NULL);
  CU_ASSERT(test_element13 == NULL);
  CU_ASSERT(test_element14 == NULL);
  CU_ASSERT(test_element15 == NULL);
  CU_ASSERT(test_element16 == NULL);
  CU_ASSERT(test_element17 == NULL);



}
