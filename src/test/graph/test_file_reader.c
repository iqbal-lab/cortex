#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <stdlib.h>
#include <limits.h>
#include <test_file_reader.h>

void test_dump_load_binary(){

  int kmer_size = 3;
  int number_of_bits_pre = 4; 
  int number_of_bits_post = 8;
  int bucket_size = 5;
  long long bad_reads = 0; 
  FILE * fout;
  int seq_length_pre,seq_length_post;
  dBGraph * db_graph_pre;
  dBGraph * db_graph_post;
  int max_retries = 10;
  int max_chunk_len_reading_from_fasta = 200;

  void print_node_binary(dBNode * node){
    db_node_print_binary(fout,node);
  }

  db_graph_pre = hash_table_new(number_of_bits_pre,bucket_size,max_retries,kmer_size);


  seq_length_pre = load_fasta_from_filename_into_graph("../data/test/graph/test_dB_graph.fasta", &bad_reads, max_chunk_len_reading_from_fasta, db_graph_pre);


  fout = fopen("../data/test/graph/dump_graph.bin", "w");

  if (fout == NULL){
    fprintf(stderr,"cannot open ../data/test/graph/dump_graph.bin");
    exit(1);
  }

  hash_table_traverse(&print_node_binary,db_graph_pre);
  
  hash_table_free(&db_graph_pre);
  CU_ASSERT(db_graph_pre==NULL);

  fclose(fout);

  
  db_graph_post = hash_table_new(number_of_bits_post,bucket_size,max_retries,kmer_size);

  seq_length_post = load_binary_from_filename_into_graph("../data/test/graph/dump_graph.bin", db_graph_post);

  CU_ASSERT_EQUAL(seq_length_post,15);
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),5);
  
  BinaryKmer tmp_seq, tmp_seq2;

  //all the kmers and their reverse complements from the reads
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);

  //kmers that should not be in the graph

  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
  dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);


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


  //check arrows
  
  Nucleotide base;

  // AAA -A-> AAA
  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base)==true);
  CU_ASSERT_EQUAL(base,Adenine);

  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base)==true);
  CU_ASSERT_EQUAL(base,Thymine);

  //GGC -T-> GCT
  CU_ASSERT(db_node_has_precisely_one_edge(test_element3, reverse,&base)==true);
  CU_ASSERT_EQUAL(base,Thymine);

  //TAG -G-> AGG
  CU_ASSERT(db_node_has_precisely_one_edge(test_element7, reverse,&base)==true);
  CU_ASSERT_EQUAL(base,Guanine);

  //AGC -C-> GCC
  CU_ASSERT(db_node_has_precisely_one_edge(test_element6, forward,&base)==true);
  CU_ASSERT_EQUAL(base,Cytosine);

  //CCT -A-> CTA
  CU_ASSERT(db_node_has_precisely_one_edge(test_element10, reverse,&base)==true);
  CU_ASSERT_EQUAL(base,Adenine);

  //add one extra arrows by had -- this breaks the graph - it is only to check that the arrows get set correctly

  db_node_set_edges(test_element1,0x02);
  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base)==false);

  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base)==true);
  CU_ASSERT_EQUAL(base,Thymine);

  CU_ASSERT(db_node_edge_exist(test_element1, Adenine, forward)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, forward)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Guanine, forward)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Thymine, forward)==false);
  
  db_node_set_edges(test_element1,0x20);
  CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base)==false);

  CU_ASSERT(db_node_edge_exist(test_element1, Adenine, forward)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, forward)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Guanine, forward)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Thymine, forward)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Adenine, reverse)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Cytosine, reverse)==true);
  CU_ASSERT(db_node_edge_exist(test_element1, Guanine, reverse)==false);
  CU_ASSERT(db_node_edge_exist(test_element1, Thymine, reverse)==true);
  
  

  hash_table_free(&db_graph_post);
  CU_ASSERT(db_graph_post == NULL);


  //Now try the same thing with big kmers

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {

      kmer_size = 33;
      number_of_bits_pre = 10; 
      number_of_bits_post =12;
      bucket_size = 10;
      bad_reads = 0; 
      
      db_graph_pre = hash_table_new(number_of_bits_pre,bucket_size,max_retries,kmer_size);


      seq_length_pre = load_fasta_from_filename_into_graph("../data/test/graph/person2.fasta", &bad_reads, max_chunk_len_reading_from_fasta, db_graph_pre);
      
      /*
	> 6 unique 33-mers
	TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
	> 13 unique 33-mers
	ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
	> 12 unique 33-mers
	GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
      */

      fout = fopen("../data/test/graph/dump_graph.bin", "w");
      
      if (fout == NULL){
	fprintf(stderr,"cannot open ../data/test/graph/dump_graph.bin");
	exit(1);
      }

      hash_table_traverse(&print_node_binary,db_graph_pre);
  
      hash_table_free(&db_graph_pre);
      CU_ASSERT(db_graph_pre==NULL);
      
      fclose(fout);
      
  
      db_graph_post = hash_table_new(number_of_bits_post,bucket_size,max_retries,kmer_size);

      seq_length_post = load_binary_from_filename_into_graph("../data/test/graph/dump_graph.bin", db_graph_post);

      CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),31);
  

      //some the kmers and their reverse complements from the reads
      test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TAACCCTAACCCTAACCCTAACCCTAACCCTAA", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("AACCCTAACCCTAACCCTAACCCTAACCCTAAC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCTAACCCTAACC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGT", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("CCCTAACCCTAACCCTAACCCTAACCCTAACCC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("CCTAACCCTAACCCTAACCCTAACCCTAACCCT", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("CTAACCCTAACCCTAACCCTAACCCTAACCCTA", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      
      //from read 2:
      test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCCTAACCCTAAC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTAGGGTTAGGGGTTAGGGTTAGGGTTAGGGT", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("CCCTAACCCTAACCCTAACCCCTAACCCTAACC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);
      test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("GGTTAGGGTTAGGGGTTAGGGTTAGGGTTAGGG", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);

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

      CU_ASSERT(test_element11 != NULL);
      CU_ASSERT(test_element12 != NULL);
      CU_ASSERT(test_element11 == test_element12);


      CU_ASSERT(test_element13 != NULL);
      CU_ASSERT(test_element14 != NULL);
      CU_ASSERT(test_element13 == test_element14);

      CU_ASSERT(test_element15 != NULL);
      CU_ASSERT(test_element16 != NULL);
      CU_ASSERT(test_element15 == test_element16);


      //now check arrows

      //Note we know that read1 forms a supernode that is a loop.

      // TAACCCTAACCCTAACCCTAACCCTAACCCTAA ----- C ----> AACCCTAACCCTAACCCTAACCCTAACCCTAAC
      CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base)==true);
      CU_ASSERT_EQUAL(base,Cytosine);
      
      CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base)==true);
      CU_ASSERT_EQUAL(base,Guanine);

      // AACCCTAACCCTAACCCTAACCCTAACCCTAAC ----C  ----> ACCCTAACCCTAACCCTAACCCTAACCCTAACC
      CU_ASSERT(db_node_has_precisely_one_edge(test_element3, forward,&base)==true);
      CU_ASSERT_EQUAL(base,Cytosine);
      
      CU_ASSERT(db_node_has_precisely_one_edge(test_element3, reverse,&base)==true);
      CU_ASSERT_EQUAL(base,Adenine);
      
      // ACCCTAACCCTAACCCTAACCCTAACCCTAACC ---C -----> CCCTAACCCTAACCCTAACCCTAACCCTAACCC
      CU_ASSERT(db_node_has_precisely_one_edge(test_element5, forward,&base)==true);
      CU_ASSERT_EQUAL(base,Cytosine);
      
      CU_ASSERT(db_node_has_precisely_one_edge(test_element5, reverse,&base)==true);
      CU_ASSERT_EQUAL(base,Thymine);
      
      //OK. Looks good.

      hash_table_free(&db_graph_post);

    }



  // Finally a test case which found a bug in binary read/write which no other test case found

      kmer_size = 17;
      number_of_bits_pre = 10; 
      number_of_bits_post =10;
      bucket_size = 30;
      bad_reads = 0; 
      
      db_graph_pre = hash_table_new(number_of_bits_pre,bucket_size,max_retries,kmer_size);


      seq_length_pre = load_fasta_from_filename_into_graph("../data/test/graph/person3.fasta", &bad_reads, max_chunk_len_reading_from_fasta, db_graph_pre);
      



      /*
	>read1 overlaps human chrom 1 
	TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
	> read 2 overlaps human chrom 1
	ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
	> read 3 does not
	GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
	> read 3 does not
	GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
	> read 3 does not
	GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
	> read 4 does not, but has too low coverage
	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
	
      */

      fout = fopen("../data/test/graph/dump_graph.bin", "w");
      
      if (fout == NULL){
	fprintf(stderr,"cannot open ../data/test/graph/dump_graph.bin");
	exit(1);
      }

      hash_table_traverse(&print_node_binary,db_graph_pre);//remember print_node_binary was defined locally, above, and prints to fout.
  
      hash_table_free(&db_graph_pre);
      CU_ASSERT(db_graph_pre==NULL);
      
      fclose(fout);
      
  
      db_graph_post = hash_table_new(number_of_bits_post,bucket_size,max_retries,kmer_size);

      seq_length_post = load_binary_from_filename_into_graph("../data/test/graph/dump_graph.bin", db_graph_post);

      //now try to traverse a supernode. This is effectively a regressiontest for a bug in graph/element.c: print_binary/read_binary
      test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TAACCCTAACCCTAACC", kmer_size, &tmp_seq),kmer_size, &tmp_seq2) ,db_graph_post);


      // this node is in the middle of this supernode: CCCTAACCCTAACCCTAACCC. You can't extend forward because there is one copy in the reads with a T afterwards
      // and one with a C afterwards. And you can' extend the other way because the first kmer CCCTAACCCTAACCCTA has two arrows in. ie is preceded by two different bases (A and T) in
      // different reads
      
      CU_ASSERT(test_element1 != NULL);


      dBNode * path_nodes[100];
      Orientation path_orientations[100];
      Nucleotide path_labels[100];
      char path_string[100];
      int limit=100;
      double avg_covg;
      int min_covg;
      int max_covg;
      boolean is_cycle;
      int len  = db_graph_supernode(test_element1, limit, &db_node_action_do_nothing, 
				    path_nodes, path_orientations, path_labels, path_string,
				    &avg_covg, &min_covg, &max_covg, &is_cycle, db_graph_post);

      CU_ASSERT(len==4); 
      CU_ASSERT(is_cycle==false);

      CU_ASSERT( (!strcmp(path_string, "AGGG") ) || (!strcmp(path_string, "ACCC"))) ;
	

      //free(path_nodes);
      //free(path_orientations);
      //free(path_labels);
      //free(path_string);
      //hash_table_free(&db_graph_post);
  
  
}

void test_coverage_is_correctly_counted_on_loading_from_file()
{


   //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits=4;
  int bucket_size   = 10;
  int seq_length;
  long long bad_reads = 0;


  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  
  int max_chunk_length=100;
  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/file_to_test_covg_of_reads.fasta", &bad_reads,max_chunk_length,  db_graph);

  /*
    >occurs once
    AATAATAATAATAATAATAATAATAAT
    >occurs twice
    GGGCAGTCTCT
    >occurs twice
    GGGCAGTCTCT
    >also occurs once
    TTTTTTTTTT
  */

  CU_ASSERT_EQUAL(seq_length,59);

  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),12);
  CU_ASSERT_EQUAL(bad_reads,0);

  BinaryKmer tmp_kmer;
  BinaryKmer tmp_kmer2;

  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==8);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==9);

  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GCA", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("CAG", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==2)
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AGT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTC", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TCT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(element_get_coverage(test_element1)==4);

  
  



  hash_table_free(&db_graph);




}



void test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph()
{

  int kmer_size = 17;
  int number_of_bits = 10; 
  int bucket_size = 30;
  long long bad_reads = 0; 
  int seq_length;
  dBGraph * db_graph;

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person3.fasta", &bad_reads, 200, db_graph);


  //OK - we have graph. Now test getting sliding windows from 
  // 1. a sequence that is all in the graph
  // 2. a sequence that has one bad base in the middle
  // 3. total garbage sequence

  FILE* fp = fopen("../data/test/graph/person3_with_errors_extended_file.fastq", "r");
  if (fp==NULL)
    {
      printf("Cannot open ../data/test/graph/person3_with_errors_extended_file.fastq in test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph");
      exit(1);
    }
  

  //allocations 
  int max_read_length=100;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);


  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);




  // GET READ 1 - this si full of sequencing errors, and it should be impossible to find any kmers that are in the graph
  int len = read_sequence_from_fastq(fp, seq, max_read_length);
  int quality_cutoff=0;
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==0);



  // GET READ 2 - this one lies entirely in the graph
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers=28);

  BinaryKmer test_kmer;
  seq_to_binary_kmer("ACCCTAACCCTAACCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("CCCTAACCCTAACCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("CCTAACCCTAACCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 

  seq_to_binary_kmer("CTAACCCTAACCCTAAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true ); 
  seq_to_binary_kmer("TAACCCTAACCCTAACC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true ); 
  seq_to_binary_kmer("AACCCTAACCCTAACCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[5],test_kmer)==true ); 
  seq_to_binary_kmer("ACCCTAACCCTAACCCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[6],test_kmer)==true ); 

  seq_to_binary_kmer("CCCTAACCCTAACCCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[7],test_kmer)==true ); 
  seq_to_binary_kmer("CCTAACCCTAACCCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[8],test_kmer)==true ); 
  seq_to_binary_kmer("CTAACCCTAACCCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[9],test_kmer)==true ); 


  seq_to_binary_kmer("TAACCCTAACCCCTAAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[10],test_kmer)==true ); 
  seq_to_binary_kmer("AACCCTAACCCCTAACC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[11],test_kmer)==true ); 
  seq_to_binary_kmer("ACCCTAACCCCTAACCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[12],test_kmer)==true ); 

  seq_to_binary_kmer("CCCTAACCCCTAACCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[13],test_kmer)==true ); 
  seq_to_binary_kmer("CCTAACCCCTAACCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[14],test_kmer)==true ); 
  seq_to_binary_kmer("CTAACCCCTAACCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[15],test_kmer)==true ); 

  // ..not going all the way to the end.
  // move on to next read



  // GET READ 3 - lots of errors again
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==0);


  // GET fourth read - this one has a single base in the middle which means kmers containing it won't be in thr graph. Otherwise the same as read3 in person3.fasta

  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==2);
  CU_ASSERT((windows->window[0]).nkmers==3);

  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 


  //and in the second window, after the interrupting T:

  CU_ASSERT((windows->window[1]).nkmers==8);
  seq_to_binary_kmer("GGGGCGGGGCCCCCTCA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("GGGCGGGGCCCCCTCAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("GGCGGGGCCCCCTCACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[2],test_kmer)==true ); 
  seq_to_binary_kmer("GCGGGGCCCCCTCACAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[3],test_kmer)==true ); 
  seq_to_binary_kmer("CGGGGCCCCCTCACACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[4],test_kmer)==true ); 
  seq_to_binary_kmer("GGGGCCCCCTCACACAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[5],test_kmer)==true ); 
  seq_to_binary_kmer("GGGCCCCCTCACACACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[6],test_kmer)==true ); 
  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[7],test_kmer)==true ); 
  

  // GET FIFTH READ - lies entirely in graph, so should just get one window
  
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==28);

  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 
  seq_to_binary_kmer("GCGGGGCGGGGCGGGGC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true ); 
  seq_to_binary_kmer("CGGGGCGGGGCGGGGCG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true ); 
  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[5],test_kmer)==true ); 
  
  //... etc ...

  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[27],test_kmer)==true ); 
  

  // READ 6 - entirely in graph
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==15);
  
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  int i;
  for (i=0; i<15; i++)
    {
      CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[i],test_kmer)==true ); 
    }

  // READ 7 - first character wrong. Vital test this one. Slightly different code path if the first kmer of all is bad.
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==27);

  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("GCGGGGCGGGGCGGGGC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 
  seq_to_binary_kmer("CGGGGCGGGGCGGGGCG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true ); 
  
  // etc

  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[26],test_kmer)==true ); 



  // READ 7 - errors spaces such that precisely two kmers can be pulled out, each in their own window

  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==2);
  CU_ASSERT((windows->window[0]).nkmers==1);

  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 

 
  CU_ASSERT((windows->window[1]).nkmers==1);
  seq_to_binary_kmer("GGGGCGGGGCCCCCTCA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[0],test_kmer)==true ); 



  // last read contains two kmers that are in the graph, but which have no edge between them. This function is not sensitive to that, and should just get one window

  printf("Start interesting\n");
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==5);

  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true ); 
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true ); 



}




void test_get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph()
{

  int kmer_size = 17;
  int number_of_bits = 10; 
  int bucket_size = 30;
  long long bad_reads = 0; 
  int seq_length;
  dBGraph * db_graph;

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person3.fasta", &bad_reads, 200, db_graph);


  //OK - we have graph. Now test getting sliding windows from 
  // 1. a sequence that is all in the graph
  // 2. a sequence that has one bad base in the middle
  // 3. total garbage sequence

  FILE* fp = fopen("../data/test/graph/person3_with_errors_extended_file.fastq", "r");
  if (fp==NULL)
    {
      printf("Cannot open ../data/test/graph/person3_with_errors_extended_file.fastq in test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph");
      exit(1);
    }
  

  //allocations 
  int max_read_length=100;
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);


  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);




  // GET READ 1 - this si full of sequencing errors, and it should be impossible to find any kmers that are in the graph
  int len = read_sequence_from_fastq(fp, seq, max_read_length);
  int quality_cutoff=0;
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==0);



  // GET READ 2 - this one lies entirely in the graph
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers=28);

  BinaryKmer test_kmer;
  seq_to_binary_kmer("ACCCTAACCCTAACCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("CCCTAACCCTAACCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("CCTAACCCTAACCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 

  seq_to_binary_kmer("CTAACCCTAACCCTAAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true ); 
  seq_to_binary_kmer("TAACCCTAACCCTAACC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true ); 
  seq_to_binary_kmer("AACCCTAACCCTAACCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[5],test_kmer)==true ); 
  seq_to_binary_kmer("ACCCTAACCCTAACCCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[6],test_kmer)==true ); 

  seq_to_binary_kmer("CCCTAACCCTAACCCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[7],test_kmer)==true ); 
  seq_to_binary_kmer("CCTAACCCTAACCCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[8],test_kmer)==true ); 
  seq_to_binary_kmer("CTAACCCTAACCCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[9],test_kmer)==true ); 


  seq_to_binary_kmer("TAACCCTAACCCCTAAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[10],test_kmer)==true ); 
  seq_to_binary_kmer("AACCCTAACCCCTAACC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[11],test_kmer)==true ); 
  seq_to_binary_kmer("ACCCTAACCCCTAACCC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[12],test_kmer)==true ); 

  seq_to_binary_kmer("CCCTAACCCCTAACCCT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[13],test_kmer)==true ); 
  seq_to_binary_kmer("CCTAACCCCTAACCCTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[14],test_kmer)==true ); 
  seq_to_binary_kmer("CTAACCCCTAACCCTAA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[15],test_kmer)==true ); 

  // ..not going all the way to the end.
  // move on to next read



  // GET READ 3 - lots of errors again
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==0);


  // GET fourth read - this one has a single base in the middle which means kmers containing it won't be in thr graph. Otherwise the same as read3 in person3.fasta

  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										   windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==2);
  CU_ASSERT((windows->window[0]).nkmers==3);

  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 

  //and in the second window, after the interrupting T:

  CU_ASSERT((windows->window[1]).nkmers==8);
  seq_to_binary_kmer("GGGGCGGGGCCCCCTCA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("GGGCGGGGCCCCCTCAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("GGCGGGGCCCCCTCACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[2],test_kmer)==true ); 
  seq_to_binary_kmer("GCGGGGCCCCCTCACAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[3],test_kmer)==true ); 
  seq_to_binary_kmer("CGGGGCCCCCTCACACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[4],test_kmer)==true ); 
  seq_to_binary_kmer("GGGGCCCCCTCACACAC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[5],test_kmer)==true ); 
  seq_to_binary_kmer("GGGCCCCCTCACACACA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[6],test_kmer)==true ); 
  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[7],test_kmer)==true ); 


  // GET FIFTH READ - lies entirely in graph, so should just get one window
  
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==28);

  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 
  seq_to_binary_kmer("GCGGGGCGGGGCGGGGC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true ); 
  seq_to_binary_kmer("CGGGGCGGGGCGGGGCG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[4],test_kmer)==true ); 
  seq_to_binary_kmer("GGGGCGGGGCGGGGCGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[5],test_kmer)==true ); 
  
  //... etc ...

  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[27],test_kmer)==true ); 
  

  // READ 6 - entirely in graph
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==15);
  
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  int i;
  for (i=0; i<15; i++)
    {
      CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[i],test_kmer)==true ); 
    }

  // READ 7 - first character wrong. Vital test this one. Slightly different code path if the first kmer of all is bad.
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==27);

  seq_to_binary_kmer("GGGCGGGGCGGGGCGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true ); 
  seq_to_binary_kmer("GGCGGGGCGGGGCGGGG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true ); 
  seq_to_binary_kmer("GCGGGGCGGGGCGGGGC", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true ); 
  seq_to_binary_kmer("CGGGGCGGGGCGGGGCG", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true ); 
  
  // etc

  seq_to_binary_kmer("GGCCCCCTCACACACAT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[26],test_kmer)==true ); 


  // READ 7 - errors spaced apart

  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==2);

  //last read - we shoud get two separate windows, as the kmers have no edge joining them
  len = read_sequence_from_fastq(fp, seq, max_read_length);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==2);

  CU_ASSERT((windows->window[0]).nkmers==4);
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[0],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[1],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[2],test_kmer)==true );
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[3],test_kmer)==true );


  CU_ASSERT((windows->window[1]).nkmers==1);
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTA", kmer_size, &test_kmer);
  CU_ASSERT(binary_kmer_comparison_operator( (windows->window[1]).kmer[0],test_kmer)==true );


}





//Assumption is that you use a bunch of fastq to build a graph, then clean it.
//You then want access to a set of fasta files that correspond to the good reads only.
void test_dumping_of_clean_fasta()
{
  
  int kmer_size = 17;
  int number_of_bits = 10; 
  int bucket_size = 30;
  long long bad_reads = 0; 
  int seq_length;
  dBGraph * db_graph;

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person3.fasta", &bad_reads, 200, db_graph);


  //OK, we now have a graph. Let's see if dumping a clean fasta gives the right answers. 

  FILE* fptr = fopen("../data/test/graph/person3_with_errors.fastq", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/graph/person3_with_errors.fastq  in test_dumping_of_clean_fasta. Exiting.\n");
      exit(1);
    }

  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      puts("new_entry has to be true for fastq\n");
      exit(1);
    }

    return read_sequence_from_fastq(fp,seq,max_read_length);
  }


  int max_read_length = 100;



  //alloc array to hold results
  char** array_of_reads= (char**) calloc(30,sizeof(char*));
  int i;
  for (i=0; i<30; i++)
    {
      array_of_reads[i]= (char*)calloc(100,sizeof(char));
    }
  int number_of_reads=0;
  read_fastq_and_print_reads_that_lie_in_graph(fptr, stdout, &file_reader, 
					       &bad_reads, max_read_length, db_graph, 
					       //false, NULL, NULL);
					       true, array_of_reads, &number_of_reads);

  CU_ASSERT_STRING_EQUAL(array_of_reads[0], "ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC");
  CU_ASSERT_STRING_EQUAL(array_of_reads[1], "GGGGCGGGGCGGGGCGGGG");
  CU_ASSERT_STRING_EQUAL(array_of_reads[2], "GGGGCGGGGCCCCCTCACACACAT");

  CU_ASSERT_STRING_EQUAL(array_of_reads[3], "GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT");
  CU_ASSERT_STRING_EQUAL(array_of_reads[4], "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

  for (i=0; i<29; i++)
    {
      free(array_of_reads[i]);
    }
  free(array_of_reads);
  hash_table_free(&db_graph);

  /*> read
ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
Examing this:  AGGGCGGGGCGGGGCAGGGCAGGGCGGAGCCCACTCACACACAT
Examing this:  GGGGCGGGGCGGGGCGGGGTGGGGCGGGGCCCCCTCACACACAT
Examing this:  GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
> read
GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
Examing this:  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
> read
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  */


}




void test_loading_of_paired_end_reads_removing_duplicates()
{

   //first set up the hash/graph
  int kmer_size = 21;
  int number_of_bits=10;
  int bucket_size   = 10;
  int seq_length;
  long long bad_reads = 0;
 
  char quality_cut_off=1; 



  // first a test where you do not remove duplicates - does all the data get loaded?
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  
  int max_read_length=100;
  seq_length = load_paired_fastq_from_filenames_into_graph("../data/test/graph/paired_end_file1_1.fastq", "../data/test/graph/paired_end_file1_2.fastq",
							   &bad_reads, quality_cut_off, max_read_length,  false, db_graph);
  CU_ASSERT(seq_length==720);
  CU_ASSERT(db_graph->unique_kmers == 243);

  hash_table_free(&db_graph);

  // then a test where you do remove duplicates from files that don't contain any

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  
  seq_length = load_paired_fastq_from_filenames_into_graph("../data/test/graph/paired_end_file1_1.fastq", "../data/test/graph/paired_end_file1_2.fastq",
							   &bad_reads, quality_cut_off, max_read_length,  true, db_graph);
  CU_ASSERT(seq_length==720);
  CU_ASSERT(db_graph->unique_kmers == 243);

  hash_table_free(&db_graph);


  //now try a pair of files containing one duplicate pair of reads (duplicated in both files) and also one read that is duplicated only in the _1 file
  // only the former pair of reads should be discarded


  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  
  seq_length = load_paired_fastq_from_filenames_into_graph("../data/test/graph/paired_end_file2_with_dup_1.fastq", "../data/test/graph/paired_end_file2_with_dup_2.fastq",
							   &bad_reads, quality_cut_off, max_read_length,  true, db_graph);

  //as before, but with four more 36bp reads
  CU_ASSERT(seq_length==864);
  CU_ASSERT(db_graph->unique_kmers == 245);

  hash_table_free(&db_graph);
  


}
