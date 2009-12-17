#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <dB_graph_population.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <stdlib.h>
#include <test_file_reader.h>
#include <seq.h>

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

  //we need the following arguments for the API but we will not use them - for duplicate removal and homopolymer breaking
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  long long dup_reads=0;
  int homopolymer_cutoff=0;

  seq_length_pre = load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/test_dB_graph.fasta", &bad_reads, &dup_reads, 20, 
										      remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
										      db_graph_pre, individual_edge_array,0);

  CU_ASSERT(seq_length_pre==16);

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


  //load_sv_trio_binary_data_from_filename_into_graph returns total number of unique kmers loaded, times kmer_length
  CU_ASSERT_EQUAL(seq_length_post,15);
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),5);

  BinaryKmer tmp_kmer1, tmp_kmer2;

  //all the kmers and their reverse complements from the reads
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


  //kmers that should not be in the graph
  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);
  CU_ASSERT(db_node_get_coverage(test_element1,individual_edge_array,0)==6); //checking that we count KMER coverage. ie this kmer occurs many times in the same read

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





  
}



void test_load_graph_binary()
{

  int kmer_size = 5;
  int number_of_bits = 7;
  int bucket_size = 5;
  // long long bad_reads = 0;
  int seq_len=0;
  int max_retries=10;
  BinaryKmer tmp_kmer1, tmp_kmer2;
  
  dBGraph* db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  //pre-prepared binary file built from this fasta
  // 
  // >read
  // AACGTTCC
  // >read
  // AACGTTCC
  // >zam
  // GTTCA
  // >read1 - will repeat 300 times
  // AAAAAAA
  // >read1
  // AAAAAAA

  // note this contains 5 unique kmers.read 1 seems to contain 4 unique kmers but only contains 3, as AACGTT is just one kmer looped back on itself


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==1)
    {
      //seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_to_load_and_dump_as_graph_bin_kmer5.bin", db_graph, individual_edge_array, 0);
      seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_kmer5_NUM_BITFIELDS_1.ctx", db_graph, individual_edge_array, 0);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==2)
    {
      seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_kmer5_NUM_BITFIELDS_2.ctx", db_graph, individual_edge_array, 0);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==3)
    {
      seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_kmer5_NUM_BITFIELDS_3.ctx", db_graph, individual_edge_array, 0);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==4)
    {
      seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_kmer5_NUM_BITFIELDS_4.ctx", db_graph, individual_edge_array, 0);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==5)
    {
      seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_kmer5_NUM_BITFIELDS_5.ctx", db_graph, individual_edge_array, 0);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==6)
    {
      seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_kmer5_NUM_BITFIELDS_6.ctx", db_graph, individual_edge_array, 0);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==7)
    {
      seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_kmer5_NUM_BITFIELDS_7.ctx", db_graph, individual_edge_array, 0);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==8)
    {
      seq_len = load_individual_binary_data_from_filename_into_graph("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_kmer5_NUM_BITFIELDS_8.ctx", db_graph, individual_edge_array, 0);
    }

  CU_ASSERT(seq_len==25);//kmers loaded * length of kmer
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 5);


  //all the nodes and their rev complements from the graph
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AACGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACGTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGTTC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GAACG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("GGAAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("TGAAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("TTTTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);

  // nodes that should not be in the graph
  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("ATATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("TGGGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("AATAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("CTCTC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("GGCGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TACTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);
  CU_ASSERT(db_node_get_coverage(test_element1,individual_edge_array,0)==4);

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

  CU_ASSERT(test_element11 == NULL);
  CU_ASSERT(test_element12 == NULL);
  CU_ASSERT(test_element13 == NULL);
  CU_ASSERT(test_element14 == NULL);
  CU_ASSERT(test_element15 == NULL);
  CU_ASSERT(test_element16 == NULL);
  CU_ASSERT(test_element17 == NULL);

  hash_table_free(&db_graph);

}



void test_load_individual_binaries_into_sv_trio()
{

  //prepare a hash table

  int kmer_size = 31;
  int number_of_bits = 0;
  int bucket_size = 41;
  int seq_len=0;
  int max_retries=82;
  BinaryKmer tmp_kmer1, tmp_kmer2;
  
  dBGraph* db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);



  //prior to this test, using graph/ target, created 3 graph/ binaries which we will load into sv_trio.
  
  //Person 0: created binary from this fasta
  //  >read1 taken from an Alu
  //  GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC

  //Person 1:
  //  > read1 different line from same Alu as person1, so will be supernode of its own
  //  GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG

  //Person2:
  //  > read1 matches first kmer of person 0, followed by a final A not G
  //  GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGA



  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==1)
    {
      load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_graph_bins_to_sv_trio_if_num_bitfields_1", db_graph);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==2)
    {
     load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_graph_bins_to_sv_trio_if_num_bitfields_2", db_graph);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==3)
    {
     load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_graph_bins_to_sv_trio_if_num_bitfields_3", db_graph);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==4)
    {
     load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_graph_bins_to_sv_trio_if_num_bitfields_4", db_graph);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==5)
    {
     load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_graph_bins_to_sv_trio_if_num_bitfields_5", db_graph);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==6)
    {
     load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_graph_bins_to_sv_trio_if_num_bitfields_6", db_graph);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==7)
    {
     load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_graph_bins_to_sv_trio_if_num_bitfields_7", db_graph);
    }
  else if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==8)
    {
     load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_graph_bins_to_sv_trio_if_num_bitfields_8", db_graph);
    }

  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 41);


  //start with a kmer that should be in person0 and person 2 only

  dBNode* test_element1_person0 = db_graph_find_node_restricted_to_specific_person_or_population(
						element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
						db_graph, individual_edge_array, 0);

  //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
  dBNode* test_element2_person0 = db_graph_find_node_restricted_to_specific_person_or_population(
						element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
						db_graph, individual_edge_array, 0);

  dBNode* test_element1_person1 = db_graph_find_node_restricted_to_specific_person_or_population(
						element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
						db_graph, individual_edge_array, 1);

  //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
  dBNode* test_element2_person1 = db_graph_find_node_restricted_to_specific_person_or_population(
						element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
						db_graph, individual_edge_array, 1);


  dBNode* test_element1_person2 = db_graph_find_node_restricted_to_specific_person_or_population(
						element_get_key(seq_to_binary_kmer("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
						db_graph, individual_edge_array, 2);

  //look for reverse complement of GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG
  dBNode* test_element2_person2 = db_graph_find_node_restricted_to_specific_person_or_population(
						element_get_key(seq_to_binary_kmer("CAGAACTCCTGGACTCAAGCGATCCTCCCAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
						db_graph, individual_edge_array, 2);



  CU_ASSERT(test_element1_person0 !=NULL);
  CU_ASSERT(test_element2_person0 !=NULL);
  CU_ASSERT(test_element1_person0==test_element2_person0);
  CU_ASSERT(test_element1_person1 ==NULL);
  CU_ASSERT(test_element2_person1 ==NULL);
  CU_ASSERT(test_element1_person2 !=NULL);
  CU_ASSERT(test_element2_person2 !=NULL);
  CU_ASSERT(test_element1_person2==test_element2_person0);
  CU_ASSERT(test_element1_person0==test_element1_person2);


  //Now some kmers that exist in person 2 only

  dBNode* test_element3_person0 = db_graph_find_node_restricted_to_specific_person_or_population(element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
												 db_graph, individual_edge_array, 0);

  dBNode* test_element3_person1 = db_graph_find_node_restricted_to_specific_person_or_population(element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
												 db_graph, individual_edge_array, 1);

  dBNode* test_element3_person2 = db_graph_find_node_restricted_to_specific_person_or_population(element_get_key(seq_to_binary_kmer("GATCGGGTGTCCGCACTAAGTTCGGCATCAA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),
												 db_graph, individual_edge_array, 2);


  CU_ASSERT(test_element3_person0==NULL);
  CU_ASSERT(test_element3_person1!=NULL);
  CU_ASSERT(test_element3_person2==NULL);

  if (test_element3_person1==NULL)
    {
      printf("test_element3_person1 is NULL");
    }
  if (test_element1_person2==NULL)
    {
      printf("test_element1_person2 is NULL");
    }

  //now check the edges

  //in person1 we know this kmer GATCGGGTGTCCGCACTAAGTTCGGCATCAA has an Thymine edge. This is in the forward direction with respect to the node, as the reverse complement
  //  of GATCGGGTGTCCGCACTAAGTTCGGCATCAA starts with an T, and so is bigger
  Nucleotide base;
  CU_ASSERT(db_node_has_precisely_one_edge(test_element3_person1, forward,&base, individual_edge_array, 1)==true);
  CU_ASSERT(db_node_edge_exist(test_element3_person1, Thymine, forward, individual_edge_array, 1)==true);


  //and the coup de grace. Does this kmer GTGGGAGGATCGCTTGAGTCCAGGAGTTCTG node have two different edges, for person 0 and person 2. ie does person 0 have precisely a G,
  // and person 2 an A?

  CU_ASSERT(db_node_has_precisely_one_edge(test_element1_person0, reverse,&base, individual_edge_array, 0)==true);
  CU_ASSERT(db_node_edge_exist(test_element1_person0, Guanine, reverse, individual_edge_array, 0)==true);
  CU_ASSERT(db_node_has_precisely_one_edge(test_element1_person2, reverse,&base, individual_edge_array, 2)==true);
  CU_ASSERT(db_node_edge_exist(test_element1_person2, Adenine, reverse, individual_edge_array, 2)==true);


  //from now on we are pretty happy. It's loading the right nodes, and putting the edges in the right place.
  //the following is not ideal from a modular code point of view - I'm going to use code in db_graph_population to get supernodes, and see if they are what I expect.

  
  int max_expected_supernode_length=30;
  dBNode * nodes_path[max_expected_supernode_length];
  Orientation orientations_path[max_expected_supernode_length];
  Nucleotide labels_path[max_expected_supernode_length];
  char seq[max_expected_supernode_length+db_graph->kmer_size+1];

  double avg_coverage=0;
  int max_coverage=0;
  int min_coverage=0;
  boolean is_cycle=false;


  db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(test_element1_person0, reverse, max_expected_supernode_length, Guanine, &db_node_action_do_nothing,
								       nodes_path,orientations_path, labels_path, seq,
								       &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
								       db_graph, individual_edge_array, 0);

  CU_ASSERT_STRING_EQUAL("GGCTGTAGTGCGCTATGCC", seq);


  db_graph_get_perfect_path_with_first_edge_for_specific_person_or_pop(test_element1_person0, reverse, max_expected_supernode_length, Adenine, &db_node_action_do_nothing,
								       nodes_path,orientations_path, labels_path, seq,
								       &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
								       db_graph, individual_edge_array, 2);

  CU_ASSERT_STRING_EQUAL("A", seq);



                             


  hash_table_free(&db_graph);
  

}







void test_load_seq_into_array()
{

  int length_of_arrays=120;
  boolean expecting_new_fasta_entry=true;
  int max_kmer_size_used_in_this_test=31;


  // **************************
  // Example 1. Load fasta containing one read, GGGG into the graph. Then use same fasta and try to pull out the array of nodes corresponding to that path through the graph

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  long long seq_loaded=0;
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson1really", &bad_reads, db_graph);
  
  CU_ASSERT(seq_loaded==4);

  
  dBNode**     path_nodes       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); 
  Orientation* path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  Nucleotide*  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  char*        path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
    {
      printf("Failed to allocate arrays for test\n");
      exit(1);
    }

  int i=0;

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';

  char tmp_seq[5000];
  
                           
  int num_of_nodes_to_read=2;
                                                                                         
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*1000);
  kmer_window->nkmers=0;

  FILE* fptr = fopen("../data/test/pop_graph/simple1.fasta", "r");

  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple1.fasta\n");
      exit(1);
    }

  int retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string,
					 seq, kmer_window, expecting_new_fasta_entry, db_graph);


  int offset=length_of_arrays-num_of_nodes_to_read;


  CU_ASSERT(retvalue==2);
  CU_ASSERT_STRING_EQUAL("CCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==reverse);
  CU_ASSERT_STRING_EQUAL("CCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==reverse);
  CU_ASSERT_STRING_EQUAL("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZG", path_string);
  CU_ASSERT(path_labels[offset]==Guanine); 
  CU_ASSERT(path_labels[offset+1]==Undefined); 
  

  //we should now hit the end of the file, but it should not affect anything in our arrays
  CU_ASSERT(load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph)==0);
  CU_ASSERT_STRING_EQUAL("CCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==reverse);
  CU_ASSERT_STRING_EQUAL("CCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==reverse);
  CU_ASSERT_STRING_EQUAL("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZG", path_string);
  CU_ASSERT(path_labels[offset]==Guanine); 
  CU_ASSERT(path_labels[offset+1]==Undefined); 


  fclose(fptr);


  // **********************************
  // Example 2. 
  // Fasta file containing GGNG. We expect to get back two nodes, both containing the dummy kmer TTT, signifying that they contain an N


  //deliberately do not clean up db_graph. Should make no difference what other paths/nodes/edges there are in th graph, so long as the fasta whose path we are following
  // has itself been loaded - we are only ever following one  path.
  
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson2really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==4);
  fptr = fopen("../data/test/pop_graph/simple2.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple2.fasta\n");
      exit(1);
    }


  num_of_nodes_to_read=2;
  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;// 120 - 2

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';



  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);


  CU_ASSERT(retvalue==2);
  CU_ASSERT(path_nodes[offset]==NULL);
  CU_ASSERT(path_orientations[offset]==forward);
  CU_ASSERT(path_labels[offset]==Undefined);
  CU_ASSERT(path_string[offset]=='G');
  CU_ASSERT_STRING_EQUAL("G", path_string+offset);


  CU_ASSERT(path_nodes[offset+1]==NULL);
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset+1]==Undefined);


  //we should now hit the end of the file, but it should not affect anything in our arrays
  CU_ASSERT(load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph)
	    ==0);
  CU_ASSERT(path_nodes[offset]==NULL);
  CU_ASSERT(path_orientations[offset]==forward);
  CU_ASSERT(path_labels[offset]==Undefined);
  CU_ASSERT(path_string[offset]=='G');

  CU_ASSERT(path_nodes[offset+1]==NULL);
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset+1]==Undefined);




  fclose(fptr);



  // ***********************************************************************************
  // Example 3: read contains only GNGNGNGNGNG


  
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson3really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==11);
  fptr = fopen("../data/test/pop_graph/simple3.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple3.fasta\n");
      exit(1);
    }


  num_of_nodes_to_read=9;

  //need more space in seq 
  free_sequence(&seq);
  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);


  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read; //-db_graph->kmer_size+1);


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);
  
  for (i=0; i< num_of_nodes_to_read; i++)
    {
      CU_ASSERT(path_nodes[offset+i]==NULL);
      CU_ASSERT(path_orientations[offset+i]==forward);
      CU_ASSERT(path_labels[offset+i]==Undefined);
    }
  

  CU_ASSERT(path_string[offset]   == 'N');
  CU_ASSERT(path_string[offset+1] == 'G');
  CU_ASSERT(path_string[offset+2] == 'N');
  CU_ASSERT(path_string[offset+3] == 'G');
  CU_ASSERT(path_string[offset+4] == 'N');
  CU_ASSERT(path_string[offset+5] == 'G');
  CU_ASSERT(path_string[offset+6] == 'N');
  CU_ASSERT(path_string[offset+7] == 'G');
  CU_ASSERT(path_string[offset+8]=='\0');

  CU_ASSERT(load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, 
				db_graph)==0);

  for (i=0; i< num_of_nodes_to_read; i++)
    {

      CU_ASSERT(path_nodes[offset+i]==NULL);
      CU_ASSERT(path_orientations[offset+i]==forward);
      CU_ASSERT(path_labels[offset+i]==Undefined);
      }

  CU_ASSERT(path_string[offset]   == 'N');
  CU_ASSERT(path_string[offset+1] == 'G');
  CU_ASSERT(path_string[offset+2] == 'N');
  CU_ASSERT(path_string[offset+3] == 'G');
  CU_ASSERT(path_string[offset+4] == 'N');
  CU_ASSERT(path_string[offset+5] == 'G');
  CU_ASSERT(path_string[offset+6] == 'N');
  CU_ASSERT(path_string[offset+7] == 'G');
  CU_ASSERT(path_string[offset+8]=='\0');
  
  fclose(fptr);





  // Example 4: fasta file contains: ACGCGCGTTTACG

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson4really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==13);
  fptr = fopen("../data/test/pop_graph/simple4.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple4.fasta\n");
      exit(1);
    }


  num_of_nodes_to_read=11;

  //need more space in seq
  free_sequence(&seq);
  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  path_nodes       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
    {
      printf("Failed to allocate arrays for test\n");
      exit(1);
    }

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);


  //  ******** here is where we test the path_string ***************
  CU_ASSERT_STRING_EQUAL("CGCGTTTACG", path_string+offset); 



  CU_ASSERT(path_nodes[offset]!=NULL);
  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==forward);


  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==reverse);
  CU_ASSERT(path_labels[offset+1]==Guanine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==reverse);
  CU_ASSERT(path_labels[offset+3]==Guanine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+5]==reverse);
  CU_ASSERT(path_labels[offset+4]==Thymine);

  CU_ASSERT_STRING_EQUAL("AAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+6]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+6]==reverse);
  CU_ASSERT(path_labels[offset+5]==Thymine);

  //note that path_nodes[offset+7] is NOT NULL. even though it encoded a TTT. This is because we check for 1's in ALL NUMBER_OF_BITFIELDS_IN_BINARY_KMER*64 bits of the BinaryKmer,
  //whereas when you read a TTT you only get 1's in the first 6 bits.
  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+7]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+7]==reverse);
  CU_ASSERT(path_labels[offset+6]==Thymine);

  CU_ASSERT_STRING_EQUAL("TAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+8]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+8]==reverse);
  CU_ASSERT(path_labels[offset+7]==Adenine );

  CU_ASSERT_STRING_EQUAL("GTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+9]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+9]==reverse);
  CU_ASSERT(path_labels[offset+8]==Cytosine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+10]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+10]==forward);
  CU_ASSERT(path_labels[offset+9]==Guanine);

  
  //now reach the end of file
  CU_ASSERT(load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry,
				db_graph)==0);
  //arrays should be unaffected


  CU_ASSERT_STRING_EQUAL("CGCGTTTACG", path_string+offset); 

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==forward);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==reverse);
  CU_ASSERT(path_labels[offset+1]==Guanine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CGC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==reverse);
  CU_ASSERT(path_labels[offset+3]==Guanine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+5]==reverse);
  CU_ASSERT(path_labels[offset+4]==Thymine);

  CU_ASSERT_STRING_EQUAL("AAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+6]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+6]==reverse);
  CU_ASSERT(path_labels[offset+5]==Thymine);

  //note that path_nodes[offset+7] is NOT NULL. even though it encoded a TTT. This is because we check for 1's in ALL 64 bits of the kmer,
  //whereas when you read a TTT you only get 1's in the first 6 bits.
  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+7]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+7]==reverse);
  CU_ASSERT(path_labels[offset+6]==Thymine);

  CU_ASSERT_STRING_EQUAL("TAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+8]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+8]==reverse);
  CU_ASSERT(path_labels[offset+7]==Adenine );

  CU_ASSERT_STRING_EQUAL("GTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+9]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+9]==reverse);
  CU_ASSERT(path_labels[offset+8]==Cytosine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+10]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+10]==forward);
  CU_ASSERT(path_labels[offset+9]==Guanine);

  fclose(fptr);
  free_sequence(&seq);





  
  // ****************************************************************************************
  // Example 5: Fasta contains a chunk of chromosome 1
  
  //  >from chromosome 1
  // AACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTC
  // no N's in this


  //cleanup, as we now want a new graph
  hash_table_free(&db_graph);


  kmer_size = 31;
  number_of_bits = 7;
  bucket_size    = 5;
  bad_reads = 0;
  max_retries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson5really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==360);
  fptr = fopen("../data/test/pop_graph/simple5.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple5.fasta\n");
      exit(1);
    }


  num_of_nodes_to_read=330;
  length_of_arrays=400;


  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  path_nodes       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
    {
      printf("Failed to allocate arrays for test\n");
      exit(1);
    }

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);  

  //path_string should be the sequence loaded, but not the first kmer - ie matches the edges
  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTC", path_string+offset);


  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==forward);
  CU_ASSERT(path_labels[offset+1]==Adenine);


  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  
  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Cytosine);

  CU_ASSERT_STRING_EQUAL("GGGTTAGGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+5]==reverse);
  CU_ASSERT(path_labels[offset+4]==Cytosine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCCTAACCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+6]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+6]==forward);
  CU_ASSERT(path_labels[offset+5]==Thymine);

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+7]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+7]==forward);
  CU_ASSERT(path_labels[offset+6]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+8]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+8]==forward);
  CU_ASSERT(path_labels[offset+7]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+9]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+9]==forward);
  CU_ASSERT(path_labels[offset+8]==Cytosine);
  

  //and one node near the end
  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCCTAACCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+308]),db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+308]==forward);
  CU_ASSERT(path_labels[offset+307]==Cytosine);


  fclose(fptr);



  // ****************************************************************************************
  // Example 6. Fasta contains a small chunk of chromosome 1 with an N inserted
  
  //  >r
  //  AACCCTAACCCTAACCCTAACCCTAACCCCTAACCNTAACCCTAACCCTAACCCTAACCTAACCCTAA


  //Deliberately do not cleanup, as we do not want a new graph - we want to be sure we can follow our path without other stuff in graqph distracting


  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson6really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==68);
  fptr = fopen("../data/test/pop_graph/simple6.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple6.fasta\n");
      exit(1);
    }


  num_of_nodes_to_read=38;//num of bases is 68
  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';

  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  CU_ASSERT_STRING_EQUAL("AACCNTAACCCTAACCCTAACCCTAACCTAACCCTAA", path_string+offset);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==forward);
  CU_ASSERT(path_labels[offset+1]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Cytosine);


  for (i=0; i<31; i++)
    {
      CU_ASSERT(path_nodes[offset+5+i]==NULL);
    }

  CU_ASSERT_STRING_EQUAL("TAACCCTAACCCTAACCCTAACCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+36]),  db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+36]==forward);
  CU_ASSERT(path_labels[offset+35]==Undefined);  // ******** <<<<<<<<<<<< this is important!!!

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+37]),  db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+37]==forward);
  CU_ASSERT(path_labels[offset+36]==Adenine);






  // ****************************************************************************************
  // Example 7. Fasta contains a small chunk of chromosome 1 with several consecutive N's inserted
  
  // >zammo
  // AACCCTAACCCTAACCCTAACCCTAACCCCTAACCNNNNNTAACCCTAACCCTAACCCTAACCTAACCCTAA


  //Deliberately do not cleanup, as we do not want a new graph - we want to be sure we can follow our path without other stuff in graqph distracting


  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson7really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==72);
  fptr = fopen("../data/test/pop_graph/simple7.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple7.fasta\n");
      exit(1);
    }

  num_of_nodes_to_read=42;//number of bases is 72
  expecting_new_fasta_entry=true;
  offset=length_of_arrays-num_of_nodes_to_read;

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);


  CU_ASSERT(retvalue==num_of_nodes_to_read);

  CU_ASSERT_STRING_EQUAL("AACCNNNNNTAACCCTAACCCTAACCCTAACCTAACCCTAA", path_string+offset);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==forward);
  CU_ASSERT(path_labels[offset+1]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Cytosine);


  for (i=0; i<35; i++)
    {
      CU_ASSERT(path_nodes[offset+5+i]==NULL);
    }

  CU_ASSERT_STRING_EQUAL("TAACCCTAACCCTAACCCTAACCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+40]),  db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+40]==forward);
  CU_ASSERT(path_labels[offset+39]==Undefined);  // ******** <<<<<<<<<<<< this is important!!!

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+41]),  db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+41]==forward);
  CU_ASSERT(path_labels[offset+40]==Adenine);


  fclose(fptr);





  // ****************************************************************************************
  // Example 8: Fasta entry starts with Ns

  // >z
  // NNNNAAACGT



  //cleanup, as we now want a new graph - changing kmer size
  hash_table_free(&db_graph);


  kmer_size = 5;
  number_of_bits = 7;
  bucket_size    = 5;
  bad_reads = 0;
  max_retries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }


  num_of_nodes_to_read=6;
  length_of_arrays=40;
  offset=length_of_arrays-num_of_nodes_to_read;
  expecting_new_fasta_entry=true;

  free_sequence(&seq);
  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);

  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  path_nodes       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
    {
      printf("Failed to allocate arrays for test\n");
      exit(1);
    }


  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson8really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==10);
  fptr = fopen("../data/test/pop_graph/simple8.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple8.fasta\n");
      exit(1);
    }



  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);


  CU_ASSERT(retvalue==num_of_nodes_to_read);
  
  CU_ASSERT_STRING_EQUAL("AACGT", path_string+offset);

  CU_ASSERT(path_nodes[offset]==NULL);
  CU_ASSERT(path_nodes[offset+1]==NULL);
  CU_ASSERT(path_nodes[offset+2]==NULL);
  CU_ASSERT(path_nodes[offset+3]==NULL);


  CU_ASSERT_STRING_EQUAL("AAACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Undefined);

  CU_ASSERT_STRING_EQUAL("AACGT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+5]==forward);
  CU_ASSERT(path_labels[offset+4]==Thymine);

  fclose(fptr);

  


  // ****************************************************************************************
  // Example 9:

  // >zam
  // AAAANAANAAANAAAANCGTCT



  num_of_nodes_to_read=18;//22 bases
  length_of_arrays=40;
  offset=length_of_arrays-num_of_nodes_to_read;
  expecting_new_fasta_entry=true;

  free_sequence(&seq);
  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq, num_of_nodes_to_read+db_graph->kmer_size, 10);


  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson9really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==22);
  fptr = fopen("../data/test/pop_graph/simple9.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple9.fasta\n");
      exit(1);
    }

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';


  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);


  CU_ASSERT(retvalue==num_of_nodes_to_read);
  CU_ASSERT_STRING_EQUAL("AANAAANAAAANCGTCT", path_string+offset);
  
  CU_ASSERT(path_nodes[offset]==NULL);
  CU_ASSERT(path_nodes[offset+1]==NULL);
  CU_ASSERT(path_nodes[offset+2]==NULL);
  CU_ASSERT(path_nodes[offset+3]==NULL);
  CU_ASSERT(path_nodes[offset+4]==NULL);
  CU_ASSERT(path_nodes[offset+5]==NULL);
  CU_ASSERT(path_nodes[offset+6]==NULL);
  CU_ASSERT(path_nodes[offset+7]==NULL);
  CU_ASSERT(path_nodes[offset+8]==NULL);
  CU_ASSERT(path_nodes[offset+9]==NULL);
  CU_ASSERT(path_nodes[offset+10]==NULL);
  CU_ASSERT(path_nodes[offset+11]==NULL);
  CU_ASSERT(path_nodes[offset+12]==NULL);
  CU_ASSERT(path_nodes[offset+13]==NULL);
  CU_ASSERT(path_nodes[offset+14]==NULL);
  CU_ASSERT(path_nodes[offset+15]==NULL);
  CU_ASSERT(path_nodes[offset+16]==NULL);


  CU_ASSERT_STRING_EQUAL("AGACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+17]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+4]==forward);
  CU_ASSERT(path_labels[offset+3]==Undefined);

  free_sequence(&seq);
  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  hash_table_free(&db_graph);
  
  fclose(fptr);

  

  // So far we have tested that when we call this function once, we get the right nodes/Nucleotides etc. What remains is to check that we can call it in a loop
  // and it will successfully keep reading the next chunk of the file, and putting it in the right place in the array. Remember we are doing this:

  // array (........, X,X,X,X,...X)        The X's are the result of a previous call. Say there are 100 of them
  // then push the array 50 places to the left and call load_seq_into_array, loading another 50 into the right-hand-end
  // Does everything go into the right index of the array? Do we get the right nodes at the start of the second batch? At that point we will read an extra base, but need the previous
  // k-1 to construct the binary kmer. What happens if we had N's at the end of the previous batch?



  // ****************************************************************************************
  // Example 10 - simplest example of loading batches. We will have an array of length 8, and load 4 bases at a time.
  // 

  // >zam
  // ACGTACGTACGTACGT


  kmer_size = 3;
  number_of_bits = 7;
  bucket_size    = 5;
  bad_reads = 0;
  max_retries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }


  num_of_nodes_to_read=2;//let's get the first 4 bases first
  length_of_arrays=8;

  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq, length_of_arrays+db_graph->kmer_size, 10);


  path_nodes       = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
    {
      printf("Failed to allocate arrays for test\n");
      exit(1);
    }

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';



  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson10really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==16);
  fptr = fopen("../data/test/pop_graph/simple10.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple10.fasta\n");
      exit(1);
    }

  //initialise the arrays
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  
  //when we load the first chunk, we expect the start of an entry..
  expecting_new_fasta_entry=true;
  //remember offset is not passed into the function call, it just tells us where we expect the answers to be
  offset=length_of_arrays-num_of_nodes_to_read;

  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  for (i=0; i<offset; i++)
    {
      CU_ASSERT(path_nodes[i]==NULL);
    }

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset]==forward);
  CU_ASSERT_STRING_EQUAL("T", path_string+offset);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==reverse);
  CU_ASSERT(path_labels[offset]==Thymine);


  //now we want to read another 4 nodes
  num_of_nodes_to_read=4;

  //Now move the data we have currently in our arrays along by Number_of_nodes_to_read
  for (i=0; i<length_of_arrays-num_of_nodes_to_read; i++)
    {
      path_nodes[i]       = path_nodes[i+num_of_nodes_to_read];
      path_orientations[i]= path_orientations[i+num_of_nodes_to_read];
      path_labels[i]      = path_labels[i+num_of_nodes_to_read];
      path_string[i]      = path_string[i+num_of_nodes_to_read];
    }

  //check we have moved it correctly! Arrays are length 8, and we have loaded 4 bases=2 kmers so far
  for (i=0; i<2; i++)
    {
      CU_ASSERT(path_nodes[i]==NULL);
      CU_ASSERT(path_orientations[i]==forward);
      CU_ASSERT(path_labels[i]==Undefined);
    }
   CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[2]), db_graph->kmer_size, tmp_seq));
   CU_ASSERT_STRING_EQUAL("T", path_string+offset-num_of_nodes_to_read);//offset

   CU_ASSERT(path_orientations[2]==forward);
   CU_ASSERT(path_labels[1]==Undefined);
   CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[3]), db_graph->kmer_size, tmp_seq));
   CU_ASSERT(path_orientations[3]==reverse);
   CU_ASSERT(path_labels[2]==Thymine);
  
  //move final kmer to the start of seq. 
  for (i=0; i<db_graph->kmer_size; i++)
  {
  seq->seq[i]=seq->seq[num_of_nodes_to_read-db_graph->kmer_size+i];
   }
  
  //and of course offset alters:
  offset=offset-num_of_nodes_to_read;

  //OK - now ready to load next batch, this time of 4 nodes, remembering that now we no longer expect a new fasta entry.
  expecting_new_fasta_entry=false;
  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);
  CU_ASSERT(retvalue==num_of_nodes_to_read); 


  CU_ASSERT_STRING_EQUAL("TACGT", path_string+offset);


  //first four entries in array unchanged
  //check first 4 entries
  for (i=0; i<2; i++)
  {
  CU_ASSERT(path_nodes[i]==NULL);
  CU_ASSERT(path_orientations[i]==forward);
  CU_ASSERT(path_labels[i]==Undefined);
  }
  
    
  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[2]==forward);
  CU_ASSERT(path_labels[1]==Undefined);
  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[3]==reverse);
  CU_ASSERT(path_labels[2]==Thymine);

  //check last 4
  CU_ASSERT_STRING_EQUAL("GTA", binary_kmer_to_seq(element_get_kmer(path_nodes[4]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[4]==forward);
  CU_ASSERT(path_labels[3]==Adenine);

  CU_ASSERT_STRING_EQUAL("GTA", binary_kmer_to_seq(element_get_kmer(path_nodes[5]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[5]==reverse);
  CU_ASSERT(path_labels[4]==Cytosine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[6]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[6]==forward);
  CU_ASSERT(path_labels[5]==Guanine);

  CU_ASSERT_STRING_EQUAL("ACG", binary_kmer_to_seq(element_get_kmer(path_nodes[7]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[7]==reverse);
  CU_ASSERT(path_labels[6]==Thymine);



  free_sequence(&seq);
  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  hash_table_free(&db_graph);


  



  // Example 11: Harder case of loading batches, this time where batch size much bigger than kmer_size


    
//  First 20 lines of chromosome 1, contains 1140 bases, but I want to not have to worry about that and just keep loading chunks until I run out. 
 // Kmer_size 31. Take an array of length 800, and load 400 nodes at a time.

//  >1 dna:chromosome chromosome:NCBI36:1:1:247249719:1
//  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC
//  TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCT
//  AACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCT
//  AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCC
//  TAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAAC
//  CCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAAC
//  CCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTA
//  ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTC
//  AGCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCA
//  CCGAAATCTGTGCAGAGGACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCT
//  GAGGAGAACGCAACTCCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGC
//  AGACACATGCTAGCGCGTCGGGGTGGAGGCGTGGCGCAGGCGCAGAGAGGCGCGCCGCGC
//  CGGCGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCA
//  GAGAGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAG
//  GCGTGGCGCAGGCGCAGAGACGCAAGCCTACGGGCGGGGGTTGGGGGGGCGTGTGTTGCA
//  GGAGCAAAGTCGCACGGCGCCGGGCTGGGGCGGGGGGAGGGTGGCGCCGTGCACGCGCAG
//  AAACTCACGTCACGGTGGCGCGGCGCAGAGACGGGTAGAACCTCAGTAATCCGAAAAGCC
//  GGGATCGACCGCCCCTTGCTTGCAGCCGGGCACTACAGGACCCGCTTGCTCACGGTGCTG
//  TGCCAGGGCGCCCCCTGCTGGCGACTAGGGCAACTGCAGGGCTCTCTTGCTTAGAGTGGT
//  
//
  

  kmer_size = 31;
  number_of_bits = 12;
  bucket_size    = 20;
  bad_reads = 0;
  max_retries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (db_graph==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }


  num_of_nodes_to_read=400;
  length_of_arrays=800;;

  seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq, length_of_arrays+db_graph->kmer_size, 40);//40 is max length of read name


  path_nodes             = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays);
  path_orientations      = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays);
  path_labels            = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  path_string            = (char*) malloc(sizeof(char)*(length_of_arrays+ max_kmer_size_used_in_this_test +1) ); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL))
    {
      printf("Failed to allocate arrays for test\n");
      exit(1);
    }

  //initialise
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
    }
  for (i=0; i< length_of_arrays+max_kmer_size_used_in_this_test; i++)
    {
      path_string[i]='Z';
    }
  path_string[length_of_arrays+max_kmer_size_used_in_this_test]='\0';





  seq_loaded = load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson11really", &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==1140);
  fptr = fopen("../data/test/pop_graph/simple11.fasta", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/simple10.fasta\n");
      exit(1);
    }

  //initialise the arrays
  for (i=0; i<length_of_arrays; i++)
    {
      path_nodes[i]=NULL;
      path_orientations[i]=forward;
      path_labels[i]=Undefined;
      path_string[i]='N';
    }
  
  //when we load the first chunk, we expect the start of an entry..
  expecting_new_fasta_entry=true;

  //remember offset is not passed into the function call, it just tells us where we expect the answers to be
  offset=length_of_arrays-num_of_nodes_to_read;

  // *********** first batch load *************
  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);

  CU_ASSERT(retvalue==num_of_nodes_to_read);

  for (i=0; i<offset; i++)
    {
      CU_ASSERT(path_nodes[i]       ==NULL);
      CU_ASSERT(path_orientations[i]==forward);
      CU_ASSERT(path_labels[i]      ==Undefined);
      CU_ASSERT(path_string[i]      =='N');
    }
  

  //take first 430 bases of file = first 400 nodes, MINUS the first 31-mer
  CU_ASSERT_STRING_EQUAL( "AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCC" ,path_string+offset);


  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq)); //reverse complement of first kmer in fasta - TAACCCTAACCCTAACCCTAACCCTAACCCT
  CU_ASSERT(path_orientations[offset]==reverse);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+1]==forward);
  CU_ASSERT(path_labels[offset]==Adenine);

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+2]==forward);
  CU_ASSERT(path_labels[offset+1]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+3]==forward);
  CU_ASSERT(path_labels[offset+2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+60]), db_graph->kmer_size, tmp_seq)); //first kmer on second line of fasta
  CU_ASSERT(path_orientations[offset+60]==reverse);
  CU_ASSERT(path_labels[offset+59]==Thymine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+120]), db_graph->kmer_size, tmp_seq)); 
  CU_ASSERT(path_orientations[offset+120]==forward);
  CU_ASSERT(path_labels[offset+119]==Thymine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+180]), db_graph->kmer_size, tmp_seq)); 
  CU_ASSERT(path_orientations[offset+180]==forward);
  CU_ASSERT(path_labels[offset+179]==Adenine);

  CU_ASSERT_STRING_EQUAL("GGTTAGGGTTAGGGTTTAGGGTTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+240]), db_graph->kmer_size, tmp_seq)); 
  CU_ASSERT(path_orientations[offset+240]==reverse);
  CU_ASSERT(path_labels[offset+239]==Cytosine);
                          
  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+360]), db_graph->kmer_size, tmp_seq)); 
  CU_ASSERT(path_orientations[offset+360]==forward);
  CU_ASSERT(path_labels[offset+359]==Cytosine);


  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCTAACCC", binary_kmer_to_seq(element_get_kmer(path_nodes[offset+399]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[offset+399]==forward);
  CU_ASSERT(path_labels[offset+398]==Cytosine);


  //first batch correctly loaded. Now load another 400, and ensure the transition is OK. num_of_nodes_to_read is already = 400
  //from now on, each new base will be a new node. 

  
  //push everything along in the arrays:
  for (i=0; i<length_of_arrays-num_of_nodes_to_read; i++)
    {
      path_nodes[i]       = path_nodes[i+num_of_nodes_to_read];
      path_orientations[i]= path_orientations[i+num_of_nodes_to_read];
      path_labels[i]      = path_labels[i+num_of_nodes_to_read];
      path_string[i]      = path_string[i+num_of_nodes_to_read];
    }


  //move final kmer to the start of seq. 

  //Note there are two separate complications relating to how far we push seq to the left
  // 1. When you start, to get n nodes, you need k + n-1 bases. After that, you only need n
  // 2. After the first batch, you always want to leave a kmer's worth of bases in seq before loading the next batch, to allow calculation of the dge between the last node of first batch, and first of next batch.
  
  //this time, you don't have an extra kmer of bases in seq from the last time around.
 for (i=0; i<db_graph->kmer_size; i++)
   {
     seq->seq[i]=seq->seq[num_of_nodes_to_read-1+i];  //because num_nodes = num_bases-k+1; so num_bases=num_nodes+k-1; so last k bases are the num_nodes,....num_nodes+k-1 -th bases. But we count from 0... 
   }
 seq->seq[db_graph->kmer_size]='\0';

  
  //double check, compare with what used to be in path_nodes[offset+399] - see above
  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCTAACCC", seq->seq);




  //load in another 400
  expecting_new_fasta_entry=false;


  // *********** second batch load *************
  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);
  CU_ASSERT(retvalue==num_of_nodes_to_read);
  

  // First the obvious - these nodes that were between offset and offset+399 should now be between 0 and 399
  // this is just checking my test essentially,


  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[offset]), db_graph->kmer_size, tmp_seq)); //reverse complement of first kmer in fasta - TAACCCTAACCCTAACCCTAACCCTAACCCT
  CU_ASSERT(path_orientations[0]==reverse);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[1]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[1]==forward);
  CU_ASSERT(path_labels[1]==Adenine);

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[2]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[2]==forward);
  CU_ASSERT(path_labels[1]==Adenine);

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[3]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[3]==forward);
  CU_ASSERT(path_labels[2]==Cytosine);

  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[60]), db_graph->kmer_size, tmp_seq)); //first kmer on second line of fasta
  CU_ASSERT(path_orientations[60]==reverse);
  CU_ASSERT(path_labels[59]==Thymine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCCT", binary_kmer_to_seq(element_get_kmer(path_nodes[120]), db_graph->kmer_size, tmp_seq)); 
  CU_ASSERT(path_orientations[120]==forward);
  CU_ASSERT(path_labels[119]==Thymine);

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[180]), db_graph->kmer_size, tmp_seq)); 
  CU_ASSERT(path_orientations[180]==forward);
  CU_ASSERT(path_labels[179]==Adenine);

  CU_ASSERT_STRING_EQUAL("GGTTAGGGTTAGGGTTTAGGGTTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[240]), db_graph->kmer_size, tmp_seq)); 
  CU_ASSERT(path_orientations[240]==reverse);
  CU_ASSERT(path_labels[239]==Cytosine);
                          
  CU_ASSERT_STRING_EQUAL("CCTAACCCTAACCCTAACCCTAACCCTAACC", binary_kmer_to_seq(element_get_kmer(path_nodes[360]), db_graph->kmer_size, tmp_seq)); 
  CU_ASSERT(path_orientations[360]==forward);
  CU_ASSERT(path_labels[359]==Cytosine);

  CU_ASSERT_STRING_EQUAL("CTAACCCTAACCCTAACCCTAACCCTAACCC", binary_kmer_to_seq(element_get_kmer(path_nodes[399]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[399]==forward);
  CU_ASSERT(path_labels[398]==Cytosine);

  //now the new nodes :-)



  //with path_string test the whole lot in one go

  CU_ASSERT_STRING_EQUAL( "AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCTAACCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTCGCGGTACCCTCAGCCGGCCCGCCCGCCCGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGACACATGCTAGCGCGTCGGGGTGGAGGCGTGGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGACACATGCTACCGCGTCCAGGGGTGGAGGCGTGGCGCAGGCGCAGAGAGGCGCACCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCC",path_string);


  
  CU_ASSERT_STRING_EQUAL("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", binary_kmer_to_seq(element_get_kmer(path_nodes[400]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[400]==reverse);

  // ******** these next two asserts are VERY IMPORTANT **********************
  // ********  they test whether we get the correct edges across points where we load more nodes into the array

  CU_ASSERT(path_labels[399]==Thymine);   
  CU_ASSERT(path_string[399]=='T'); //iqbal

  CU_ASSERT_STRING_EQUAL("AACCCTAACCCTAACCCTAACCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[401]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[401]==forward);
  CU_ASSERT(path_labels[400]==Adenine);
  CU_ASSERT(path_string[400]=='A');

  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCTAACCCTAA", binary_kmer_to_seq(element_get_kmer(path_nodes[402]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[402]==forward);
  CU_ASSERT(path_labels[401]==Adenine);
  CU_ASSERT(path_string[401]=='A');

  CU_ASSERT_STRING_EQUAL("CCCTAACCCTAACCCTAACCCTAACCCTAAC", binary_kmer_to_seq(element_get_kmer(path_nodes[403]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[403]==forward);
  CU_ASSERT(path_labels[402]==Cytosine);
  CU_ASSERT(path_string[402]=='C');


  CU_ASSERT_STRING_EQUAL("ACCCTAACCCTAACCCTAACCCCTAACCCTA", binary_kmer_to_seq(element_get_kmer(path_nodes[420]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[420]==forward);
  CU_ASSERT(path_labels[419]==Adenine);
  CU_ASSERT(path_string[419]=='A');

  CU_ASSERT_STRING_EQUAL("GAGAGGCGCACCGCGCCGGCGCAGGCGCAGA", binary_kmer_to_seq(element_get_kmer(path_nodes[780]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[780]==forward);
  CU_ASSERT(path_labels[779]==Adenine);
  CU_ASSERT(path_string[779]=='A');


  CU_ASSERT_STRING_EQUAL("CGCAGGCGCAGAGACACATGCTAGCGCGTCC", binary_kmer_to_seq(element_get_kmer(path_nodes[799]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[799]==forward);
  CU_ASSERT(path_labels[798]==Cytosine);
  CU_ASSERT(path_string[798]=='C');


  //prepare to load next 400
  
  //push everything along in the arrays:
  for (i=0; i<length_of_arrays-num_of_nodes_to_read; i++)
    {
      path_nodes[i]       = path_nodes[i+num_of_nodes_to_read];
      path_orientations[i]= path_orientations[i+num_of_nodes_to_read];
      path_labels[i]      = path_labels[i+num_of_nodes_to_read];
      path_string[i]      = path_string[i+num_of_nodes_to_read];
    }


  //move final kmer to the start of seq. Note final means at the end of seq, which has got num_of_nodes_to_read PLUS ONE in it (you preloaded with the previous kmer)
  for (i=0; i<db_graph->kmer_size; i++)
   {
     seq->seq[i]=seq->seq[num_of_nodes_to_read+i];  //num_bases=num_nodes+k-1, but we have an extra base,  so last k bases are the num_nodes+1,....num_nodes+k -th bases. 
                                                                               //But we count from 0...
   }
  seq->seq[db_graph->kmer_size]='\0';


  //load the next 400 - or at least try to. This time there are only 1140-(800+31-1)=310 nodes before you run out of file.
  
  expecting_new_fasta_entry=false;


  // *************** third batch load *****************
  retvalue = load_seq_into_array(fptr, num_of_nodes_to_read, length_of_arrays, path_nodes, path_orientations, path_labels, path_string, seq, kmer_window, expecting_new_fasta_entry, db_graph);
  CU_ASSERT(retvalue==310);

  //check we have the transition correct, between previously loaded nodes and new ones

  //this node used to be at index 799:
  CU_ASSERT_STRING_EQUAL("CGCAGGCGCAGAGACACATGCTAGCGCGTCC", binary_kmer_to_seq(element_get_kmer(path_nodes[399]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[399]==forward);
  CU_ASSERT(path_labels[398]==Cytosine);
  CU_ASSERT(path_string[398]=='C');

  
  //finally, just check the last one. What is it's index? Well, we asked for 400, but only got 310. So we expect to find last one at 799-90=709.

  CU_ASSERT_STRING_EQUAL("ACCACTCTAAGCAAGAGAGCCCTGCAGTTGC", binary_kmer_to_seq(element_get_kmer(path_nodes[709]), db_graph->kmer_size, tmp_seq));
  CU_ASSERT(path_orientations[709]==reverse);
  CU_ASSERT(path_labels[708]==Thymine);
  CU_ASSERT(path_string[708]=='T');




  free_sequence(&seq);
  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  hash_table_free(&db_graph);



  free(kmer_window->kmer);
  free(kmer_window);


  //  ****************************************************************************************************
  // Example 12: Example with N at the end of the previous batch

  // ****************************************************************************************************
  // Example 13: Case where there is an N at start of next batch







}
