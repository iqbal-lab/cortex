/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */

#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <dB_graph_population.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <stdlib.h>
#include <test_file_reader.h>
#include <seq.h>
#include <graph_info.h>
#include <limits.h>

void test_dump_load_sv_trio_binary(){

  int kmer_size = 3;
  int number_of_bits_pre = 4; 
  int number_of_bits_post = 8;
  int bucket_size = 5;
  long long bad_reads = 0; 
  int max_retries=10;
  int max_chunk_len_reading_from_fasta = 200;

  //FILE * fout = fopen("../data/test/pop_graph/dump_cortex_var_graph.bin", "w");
  long long  seq_length_parsed_pre=0;
  long long seq_length_loaded_pre=0;


  dBGraph * db_graph_pre;
  dBGraph * db_graph_post;

  //void print_node_binary(dBNode * node){
  //  db_node_print_multicolour_binary(fout,node);
  //}

  db_graph_pre = hash_table_new(number_of_bits_pre,bucket_size,max_retries,kmer_size);

  //we need the following arguments for the API but we will not use them - for duplicate removal and homopolymer breaking
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  long long dup_reads=0;
  int homopolymer_cutoff=0;

  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/test_dB_graph.fasta", &seq_length_parsed_pre, &seq_length_loaded_pre, &bad_reads, &dup_reads, 20, 
										      remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
										      db_graph_pre, individual_edge_array,0);

  CU_ASSERT(seq_length_parsed_pre==16);
  CU_ASSERT(seq_length_loaded_pre==16);

  //  if (fout == NULL){
  //  fprintf(stderr,"cannot open ../data/test/pop_graph/dump_sv_trio_graph.bin");
  //  exit(1);
  //}

  //hash_table_traverse(&print_node_binary,db_graph_pre);
  GraphInfo ginfo;
  graph_info_initialise(&ginfo);
  graph_info_set_seq(&ginfo, 0, seq_length_parsed_pre);
  db_graph_dump_binary("../data/test/pop_graph/dump_cortex_var_graph.bin", &db_node_condition_always_true, db_graph_pre, &ginfo);


  hash_table_free(&db_graph_pre);
  CU_ASSERT(db_graph_pre==NULL);

  //fclose(fout);

  db_graph_post = hash_table_new(number_of_bits_post,bucket_size,10,kmer_size);
  int num_cols_in_binary=-1;

  int array_mean_readlens[NUMBER_OF_COLOURS];
  int* array_mean_readlens_ptrs[NUMBER_OF_COLOURS];
  long long array_total_seq[NUMBER_OF_COLOURS];
  long long* array_total_seq_ptrs[NUMBER_OF_COLOURS];

  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      array_mean_readlens[i]=0;
      array_mean_readlens_ptrs[i]=&array_mean_readlens[i];
      array_total_seq[i]=0;
      array_total_seq_ptrs[i]=&array_total_seq[i];
    }
  long long seq_length_post = load_multicolour_binary_from_filename_into_graph("../data/test/pop_graph/dump_cortex_var_graph.bin", db_graph_post, &num_cols_in_binary,
									       array_mean_readlens_ptrs, array_total_seq_ptrs);


  CU_ASSERT(num_cols_in_binary==NUMBER_OF_COLOURS);
  CU_ASSERT(array_mean_readlens[0]==0);
  CU_ASSERT(array_total_seq[0]==seq_length_parsed_pre);

  //load_multicolour_binary_data_from_filename_into_graph returns total number of unique kmers loaded, times kmer_length
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



  //Now try the same thing with big kmers

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {

      kmer_size = 33;
      number_of_bits_pre = 10; 
      number_of_bits_post =12;
      bucket_size = 10;
      bad_reads = 0; 
      
      db_graph_pre = hash_table_new(number_of_bits_pre,bucket_size,max_retries,kmer_size);


      load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/person2.fasta", &seq_length_parsed_pre, &seq_length_loaded_pre,
									 &bad_reads,&dup_reads, max_chunk_len_reading_from_fasta, 
									 remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
									 db_graph_pre, individual_edge_array, 0);
      
      
      //	> 6 unique 33-mers
      //	TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
      //	> 13 unique 33-mers
      //	ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
      //	> 12 unique 33-mers
      //	GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
      

      graph_info_initialise(&ginfo);
      graph_info_set_seq(&ginfo, 0, seq_length_parsed_pre);
      db_graph_dump_binary("../data/test/pop_graph/dump_cortex_var_graph_2.bin", &db_node_condition_always_true, db_graph_pre, &ginfo);
      //hash_table_traverse(&print_node_binary,db_graph_pre);
  
      hash_table_free(&db_graph_pre);
      CU_ASSERT(db_graph_pre==NULL);
      
      db_graph_post = hash_table_new(number_of_bits_post,bucket_size,max_retries,kmer_size);

      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  array_mean_readlens[i]=0;
	  array_total_seq[i]=0;
	}
      
      seq_length_post = load_multicolour_binary_from_filename_into_graph("../data/test/pop_graph/dump_cortex_var_graph_2.bin", db_graph_post, &num_cols_in_binary,
									 array_mean_readlens_ptrs, array_total_seq_ptrs);


      CU_ASSERT(num_cols_in_binary==NUMBER_OF_COLOURS);

      CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post),31);
  

      //some the kmers and their reverse complements from the reads
      test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TAACCCTAACCCTAACCCTAACCCTAACCCTAA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("AACCCTAACCCTAACCCTAACCCTAACCCTAAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCTAACCCTAACC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("CCCTAACCCTAACCCTAACCCTAACCCTAACCC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("CCTAACCCTAACCCTAACCCTAACCCTAACCCT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("CTAACCCTAACCCTAACCCTAACCCTAACCCTA", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      
      //from read 2:
      test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCCTAACCCTAAC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTAGGGTTAGGGGTTAGGGTTAGGGTTAGGGT", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("CCCTAACCCTAACCCTAACCCCTAACCCTAACC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
      test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("GGTTAGGGTTAGGGGTTAGGGTTAGGGTTAGGG", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);

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
      CU_ASSERT(db_node_has_precisely_one_edge(test_element1, forward,&base, individual_edge_array, 0)==true);
      CU_ASSERT_EQUAL(base,Cytosine);
      
      CU_ASSERT(db_node_has_precisely_one_edge(test_element1, reverse,&base, individual_edge_array, 0)==true);
      CU_ASSERT_EQUAL(base,Guanine);

      // AACCCTAACCCTAACCCTAACCCTAACCCTAAC ----C  ----> ACCCTAACCCTAACCCTAACCCTAACCCTAACC
      CU_ASSERT(db_node_has_precisely_one_edge(test_element3, forward,&base, individual_edge_array, 0)==true);
      CU_ASSERT_EQUAL(base,Cytosine);
      
      CU_ASSERT(db_node_has_precisely_one_edge(test_element3, reverse,&base, individual_edge_array, 0)==true);
      CU_ASSERT_EQUAL(base,Adenine);
      
      // ACCCTAACCCTAACCCTAACCCTAACCCTAACC ---C -----> CCCTAACCCTAACCCTAACCCTAACCCTAACCC
      CU_ASSERT(db_node_has_precisely_one_edge(test_element5, forward,&base, individual_edge_array, 0)==true);
      CU_ASSERT_EQUAL(base,Cytosine);
      
      CU_ASSERT(db_node_has_precisely_one_edge(test_element5, reverse,&base, individual_edge_array, 0)==true);
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
  
  
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/person3.fasta", &seq_length_parsed_pre, &seq_length_loaded_pre,
								     &bad_reads,&dup_reads, max_chunk_len_reading_from_fasta,
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
								     db_graph_pre, individual_edge_array,0);
  



//    >read1 overlaps human chrom 1 
 //   TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
//    > read 2 overlaps human chrom 1
//    ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
//    > read 3 does not
//    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
//    > read 3 does not
//    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
//    > read 3 does not
//    GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
//    > read 4 does not, but has too low coverage
//    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    

  graph_info_initialise(&ginfo);
  graph_info_set_seq(&ginfo, 0, seq_length_parsed_pre);
  db_graph_dump_binary("../data/test/pop_graph/dump_cortex_var_graph_3.bin", &db_node_condition_always_true, db_graph_pre, &ginfo);
  hash_table_free(&db_graph_pre);
  CU_ASSERT(db_graph_pre==NULL);
  
  
  db_graph_post = hash_table_new(number_of_bits_post,bucket_size,max_retries,kmer_size);

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      array_mean_readlens[i]=0;
      array_total_seq[i]=0;
    }
  seq_length_post = load_multicolour_binary_from_filename_into_graph("../data/test/pop_graph/dump_cortex_var_graph_3.bin", db_graph_post, &num_cols_in_binary,
								     array_mean_readlens_ptrs, array_total_seq_ptrs);


  CU_ASSERT(num_cols_in_binary==NUMBER_OF_COLOURS);

  //now try to traverse a supernode. This is effectively a regressiontest for a bug in graph/element.c: print_binary/read_binary
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TAACCCTAACCCTAACC", kmer_size, &tmp_kmer1),kmer_size, &tmp_kmer2) ,db_graph_post);
  
  
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
  int len  = db_graph_supernode_for_specific_person_or_pop(test_element1, limit, &db_node_action_do_nothing, 
							   path_nodes, path_orientations, path_labels, path_string,
							   &avg_covg, &min_covg, &max_covg, &is_cycle, db_graph_post, individual_edge_array, 0);
  
  CU_ASSERT(len==4); 
  CU_ASSERT(is_cycle==false);
  
  CU_ASSERT( (!strcmp(path_string, "AGGG") ) || (!strcmp(path_string, "ACCC"))) ;
  
  
  hash_table_free(&db_graph_post);
  


  
}



void test_load_singlecolour_binary()
{

  int kmer_size = 5;
  int number_of_bits = 7;
  int bucket_size = 5;
  // long long bad_reads = 0;
  //int seq_len=0;
  int max_retries=10;
  BinaryKmer tmp_kmer1, tmp_kmer2;
  
  dBGraph* db_graph_pre = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);



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




  //we need the following arguments for the API but we will not use them - for duplicate removal and homopolymer breaking
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  long long dup_reads=0;
  long long bad_reads=0;
  
  int homopolymer_cutoff=0;
  long long seq_length_parsed_pre=0;
  long long seq_length_loaded_pre=0;
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/pop_graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio.fasta", 
								     &seq_length_parsed_pre, &seq_length_loaded_pre,
								     &bad_reads, &dup_reads, 20, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
								     db_graph_pre, individual_edge_array,0);


  GraphInfo ginfo;
  graph_info_initialise(&ginfo);
  graph_info_set_seq(&ginfo, 0, seq_length_parsed_pre);
  db_graph_dump_single_colour_binary_of_colour0("../data/test/pop_graph/dump_single_colour_cortex_var_graph.bin", &db_node_condition_always_true, db_graph_pre, &ginfo);
  hash_table_free(&db_graph_pre);
  CU_ASSERT(db_graph_pre==NULL);

  dBGraph* db_graph_post = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  int mean_readlen=0;
  long long total_seq=0;
  int seq_length_post = load_single_colour_binary_data_from_filename_into_graph("../data/test/pop_graph/dump_single_colour_cortex_var_graph.bin", db_graph_post,
										&mean_readlen, &total_seq,
										true, individual_edge_array,0, false,0);



  CU_ASSERT(seq_length_post==25);//kmers loaded * length of kmer
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph_post), 5);


  //all the nodes and their rev complements from the graph
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AACGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACGTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGTTC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GAACG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("GGAAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("GTTCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("TGAAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("TTTTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);

  // nodes that should not be in the graph
  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("ATATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("TGGGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("AATAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("CTCTC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("GGCGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TACTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


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

  hash_table_free(&db_graph_post);

}



void test_load_individual_binaries_into_sv_trio()
{


  if (NUMBER_OF_COLOURS<3)
    {
      printf("This test is redundant unless the compile-time flag NUMBER_OF_COLOURS is set to a value >=3. This is set by typing make NUM_COLS=3 blah\n");
      return;
    }
  else
    {
      
      
      //prepare a hash table
      
      int kmer_size = 31;
      int number_of_bits = 0;
      int bucket_size = 41;
      //int seq_len=0;
      int max_retries=82;
      BinaryKmer tmp_kmer1, tmp_kmer2;
      
      dBGraph* db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
      
      


     //first dump 3 single-colour graphs
      
      //Person 0: created binary from this fasta
      //  >read1 taken from an Alu
      //  GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
      
      //Person 1:
      //  > read1 different line from same Alu as person1, so will be supernode of its own
      //  GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
      
      //Person2:
      //  > read1 matches first kmer of person 0, followed by a final A not G
      //  GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGA


     boolean remove_duplicates_single_endedly=false;
     boolean break_homopolymers=false;
     long long dup_reads=0;
     int homopolymer_cutoff=0;
     long long bad_reads=0;
     long long seq_read=0;
     long long seq_loaded1=0;
     long long seq_loaded2=0;
     long long seq_loaded3=0;
     load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_person0.fasta", &seq_read, &seq_loaded1,
									&bad_reads, &dup_reads, 200,
									remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
									db_graph, individual_edge_array,0);

     GraphInfo ginfo;
     graph_info_initialise(&ginfo);
     graph_info_set_seq(&ginfo, 0, seq_loaded1);
     db_graph_dump_single_colour_binary_of_colour0("../data/test/pop_graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_person0_kmer31.ctx", 
						   &db_node_condition_always_true, db_graph, &ginfo);
     hash_table_free(&db_graph);


     db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
     load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_person1.fasta", &seq_read,&seq_loaded2,
									&bad_reads, &dup_reads, 200,
									remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
									db_graph, individual_edge_array,0);
     graph_info_initialise(&ginfo);
     graph_info_set_seq(&ginfo, 0, seq_loaded2);
     db_graph_dump_single_colour_binary_of_colour0("../data/test/pop_graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_person1_kmer31.ctx", 
						   &db_node_condition_always_true, db_graph, &ginfo);
     hash_table_free(&db_graph);




     db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
     load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_person2.fasta", &seq_read,&seq_loaded3,
									&bad_reads, &dup_reads, 200,
									remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
									db_graph, individual_edge_array,0);
     graph_info_initialise(&ginfo);
     graph_info_set_seq(&ginfo, 0, seq_loaded3);
     db_graph_dump_single_colour_binary_of_colour0("../data/test/pop_graph/fasta_for_dumping_by_graph_and_reload_by_sv_trio_person2_kmer31.ctx", 
						   &db_node_condition_always_true, db_graph, &ginfo);
     hash_table_free(&db_graph);


     
     db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
     graph_info_initialise(&ginfo);
     int first_colour=0;
     load_population_as_binaries_from_graph("../data/test/pop_graph/trio_filelist_for_testing_loading_singlecolour_bins_into_multicol_bin", first_colour, true, db_graph, &ginfo,
					    false, 0);
     CU_ASSERT(ginfo.total_sequence[0]==seq_loaded1);
     CU_ASSERT(ginfo.total_sequence[1]==seq_loaded2);
     CU_ASSERT(ginfo.total_sequence[2]==seq_loaded3);
      
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

}






void test_coverage_is_correctly_counted_on_loading_from_file()
{


   //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits=4;
  int bucket_size   = 10;
  long long seq_read=0;
  long long seq_loaded=0;
  long long bad_reads = 0; long long dup_reads=0;
  boolean remove_duplicates_single_endedly=false; 
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;


  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  
  int max_chunk_length=100;
  
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/file_to_test_covg_of_reads.fasta", &seq_read, &seq_loaded,
								     &bad_reads,&dup_reads,max_chunk_length,  
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
								     db_graph, individual_edge_array, 0);

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

  CU_ASSERT_EQUAL(seq_read,59);
  CU_ASSERT_EQUAL(seq_loaded,59);

  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),12);
  CU_ASSERT_EQUAL(bad_reads,0);

  BinaryKmer tmp_kmer;
  BinaryKmer tmp_kmer2;

  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==8);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==9);

  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GCA", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("CAG", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==2)
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AGT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTC", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==2);
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TCT", kmer_size, &tmp_kmer), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(db_node_get_coverage(test_element1, individual_edge_array,0)==4);

  

  hash_table_free(&db_graph);




}


void test_getting_sliding_windows_where_you_break_at_kmers_not_in_db_graph()
{

  int kmer_size = 17;
  int number_of_bits = 10; 
  int bucket_size = 30;
  long long bad_reads = 0; long long dup_reads=0; 
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;

  long long seq_read=0;
  long long seq_loaded=0;
  dBGraph * db_graph;

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/person3.fasta", &seq_read, &seq_loaded,
								     &bad_reads,&dup_reads, 200, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
								     db_graph, individual_edge_array, 0);


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

  int fq_ascii_offset = 33;


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
  int len = read_sequence_from_fastq(fp, seq, max_read_length, fq_ascii_offset);
  int quality_cutoff=0;
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==0);



  // GET READ 2 - this one lies entirely in the graph
  len = read_sequence_from_fastq(fp, seq, max_read_length, fq_ascii_offset);
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
  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
  get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph);  
  CU_ASSERT(windows->nwindows==0);


  // GET fourth read - this one has a single base in the middle which means kmers containing it won't be in thr graph. Otherwise the same as read3 in person3.fasta

  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
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
  
  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
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
  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
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
  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
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

  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
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
  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
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


  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  hash_table_free(&db_graph);
}




void test_get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph()
{

  int kmer_size = 17;
  int number_of_bits = 10; 
  int bucket_size = 30;
  long long bad_reads = 0; long long dup_reads=0; 
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;

  long long seq_read=0;
  long long seq_loaded=0;
  dBGraph * db_graph;

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/person3.fasta", &seq_read, &seq_loaded,
								     &bad_reads,&dup_reads, 200, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
								     db_graph, individual_edge_array, 0);


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

  int fq_ascii_offset=33; 


  // GET READ 1 - this si full of sequencing errors, and it should be impossible to find any kmers that are in the graph
  int len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
  int quality_cutoff=0;
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
  CU_ASSERT(windows->nwindows==0);



  // GET READ 2 - this one lies entirely in the graph
  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
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
  len = read_sequence_from_fastq(fp, seq, max_read_length,fq_ascii_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
  CU_ASSERT(windows->nwindows==0);


  // GET fourth read - this one has a single base in the middle which means kmers containing it won't be in thr graph. Otherwise the same as read3 in person3.fasta

  len = read_sequence_from_fastq(fp, seq, max_read_length, fq_ascii_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										   windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
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
  
  len = read_sequence_from_fastq(fp, seq, max_read_length, fq_ascii_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
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
  len = read_sequence_from_fastq(fp, seq, max_read_length, fq_ascii_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
  CU_ASSERT(windows->nwindows==1);
  CU_ASSERT((windows->window[0]).nkmers==15);
  
  seq_to_binary_kmer("TTTTTTTTTTTTTTTTT", kmer_size, &test_kmer);
  int i;
  for (i=0; i<15; i++)
    {
      CU_ASSERT(binary_kmer_comparison_operator( (windows->window[0]).kmer[i],test_kmer)==true ); 
    }

  // READ 7 - first character wrong. Vital test this one. Slightly different code path if the first kmer of all is bad.
  len = read_sequence_from_fastq(fp, seq, max_read_length, fq_ascii_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
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

  len = read_sequence_from_fastq(fp, seq, max_read_length, fq_ascii_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
  CU_ASSERT(windows->nwindows==2);

  //last read - we shoud get two separate windows, as the kmers have no edge joining them
  len = read_sequence_from_fastq(fp, seq, max_read_length, fq_ascii_offset);
  get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, len, quality_cutoff, 
										windows, max_windows, max_kmers, db_graph, individual_edge_array, 0);  
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

  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  hash_table_free(&db_graph);
}





//Assumption is that you use a bunch of fastq to build a graph, then clean it.
//You then want access to a set of fasta files that correspond to the good reads only.
void test_dumping_of_clean_fasta()
{
  
  int kmer_size = 17;
  int number_of_bits = 10; 
  int bucket_size = 30;
  long long bad_reads = 0; long long dup_reads=0; 
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;

  long long seq_read=0;
  long long seq_loaded=0;
  dBGraph * db_graph;

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/person3.fasta", &seq_read, &seq_loaded,
								     &bad_reads,&dup_reads, 200, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
								     db_graph, individual_edge_array, 0);


  //OK, we now have a graph. Let's see if dumping a clean fasta gives the right answers. 

  FILE* fptr = fopen("../data/test/graph/person3_with_errors.fastq", "r");
  if (fptr==NULL)
    {
      printf("Cannot open ../data/test/graph/person3_with_errors.fastq  in test_dumping_of_clean_fasta. Exiting.\n");
      exit(1);
    }

  int fq_ascii_offset=33; 
  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      puts("new_entry has to be true for fastq\n");
      exit(1);
    }

    return read_sequence_from_fastq(fp,seq,max_read_length, fq_ascii_offset);
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

  for (i=0; i<30; i++)
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

  long long bad_reads = 0; long long dup_reads=0;
  boolean remove_duplicates=false;
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  int ascii_offset=33;
  char quality_cut_off=0; 

  long long seq_read=0;
  long long seq_loaded=0;

  // first a test where the file contains no duplicates and you do not try to remove duplicates - does all the data get loaded?
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

  FileFormat format=FASTQ;
  int max_read_length=100;
  long long readlen_array[101];
  long long* readlen_array_ptrs[101];
  int j;
  for (j=0; j<=100; j++)
    {
      readlen_array[j]=0;
      readlen_array_ptrs[j]=&readlen_array[j];
    }
  load_paired_end_data_from_filenames_into_graph_of_specific_person_or_pop("../data/test/graph/paired_end_file1_1.fastq", "../data/test/graph/paired_end_file1_2.fastq",
									   format,
									   &seq_read, &seq_loaded,readlen_array_ptrs,
									   &bad_reads, quality_cut_off, max_read_length,  
									   &dup_reads, remove_duplicates, break_homopolymers, homopolymer_cutoff, ascii_offset, 
									   db_graph, individual_edge_array, 0);
  CU_ASSERT(seq_read==720);
  CU_ASSERT(seq_loaded==665);//because there are various N's in the sequence, we lose some bases

  CU_ASSERT(db_graph->unique_kmers == 243);

  hash_table_free(&db_graph);

  // then a test where you try to  remove duplicates from files that don't contain any

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  remove_duplicates=true;
  seq_read=0;
  seq_loaded=0;
  for (j=0; j<=100; j++)
    {
      readlen_array[j]=0;
    }
  load_paired_end_data_from_filenames_into_graph_of_specific_person_or_pop("../data/test/graph/paired_end_file1_1.fastq", "../data/test/graph/paired_end_file1_2.fastq",
									   format,
									   &seq_read, &seq_loaded,readlen_array_ptrs,
									   &bad_reads, quality_cut_off, max_read_length,  
									   &dup_reads,remove_duplicates, break_homopolymers, homopolymer_cutoff, 
									   ascii_offset, db_graph, individual_edge_array, 0);
  CU_ASSERT(seq_read==720);
  CU_ASSERT(seq_loaded==665);
  CU_ASSERT(db_graph->unique_kmers == 243);
  CU_ASSERT(dup_reads==0);
  hash_table_free(&db_graph);


  //now try a pair of files containing one duplicate pair of reads (corresponding mates duplicated in both files) and also one read that is duplicated only in the _1 file
  // only the former pair of reads should be discarded

  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  seq_read=0;
  seq_loaded=0;
  for (j=0; j<=100; j++)
    {
      readlen_array[j]=0;
    }

  load_paired_end_data_from_filenames_into_graph_of_specific_person_or_pop("../data/test/graph/paired_end_file2_with_dup_1.fastq", 
									   "../data/test/graph/paired_end_file2_with_dup_2.fastq",
									   format,&seq_read, &seq_loaded,readlen_array_ptrs,
									   &bad_reads, quality_cut_off, max_read_length, 
									   &dup_reads, remove_duplicates, break_homopolymers, homopolymer_cutoff, ascii_offset, 
									   db_graph, individual_edge_array, 0);

  
  //as before, but with four more 36bp reads
  CU_ASSERT(seq_read==864);
  CU_ASSERT(seq_loaded==737);//all the sequence we have read is loaded except for 2 reads, plus the effect of the Ns

  CU_ASSERT(db_graph->unique_kmers == 245);
  CU_ASSERT(dup_reads==2);


  hash_table_free(&db_graph);
  dup_reads=0;
  
  // now take a pair of files where there is an original mate pair, one proper duplicate, and then 3 other mate-pairs that might be confused for duplicates
  // - one read is a dup of the original, and the other is reverse complement of the other.
  // Of these, only one pair is what we want to call a duplicate. However our filter cannot distinguish between the following cases
  //   1. a second mate pair that is identical to a previous mate pair
  //   2. Pair1_read1 is identical to Pair2_read1, and Pair2_read2 is identical to Pair3_read1.
  // the final mate pair in our file has one mate which is the same as a left_mate in one read pair, and a right mate which is the same
  // as the right mate in another read pair. So it isn't really a dup, but we discard it anyway as we cannot tell the difference.

  int count_file_pairs=0;
  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  seq_read=0;
  seq_loaded=0;
  for (j=0; j<=100; j++)
    {
      readlen_array[j]=0;
    }

  load_list_of_paired_end_files_into_graph_of_specific_person_or_pop( "../data/test/graph/list_paired_end_file3_left", 
								      "../data/test/graph/list_paired_end_file3_right", 
								      format,
								      quality_cut_off, max_read_length,
								      &seq_read, &seq_loaded,readlen_array_ptrs,
								      &bad_reads, &dup_reads, &count_file_pairs, 
								      remove_duplicates, break_homopolymers, homopolymer_cutoff, ascii_offset, 
								      db_graph, individual_edge_array, 0);
  

  CU_ASSERT(seq_read==360); //five 36bp reads, left and right
  CU_ASSERT(seq_loaded==216);
  CU_ASSERT(dup_reads==4);
  CU_ASSERT(count_file_pairs==1);

  hash_table_free(&db_graph);
}



void test_loading_of_single_ended_reads_removing_duplicates()
{

   //first set up the hash/graph
  int kmer_size = 21;
  int number_of_bits=10;
  int bucket_size   = 10;
  long long seq_read=0;
  long long seq_loaded=0;
  long long bad_reads = 0; long long dup_reads=0;
  boolean remove_duplicates=false;
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  
  char quality_cut_off=1; 



  // first a test where the file contains no duplicates and you do not try to remove duplicates - does all the data get loaded?
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  
  int max_read_length=100;
  int ascii_offset=33;
  long long readlen_array[101];
  long long* readlen_array_ptrs[101];
  int j;
  for (j=0; j<=100; j++)
    {
      readlen_array[j]=0;
      readlen_array_ptrs[j]=&readlen_array[j];
    }

  load_fastq_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/paired_end_file1_1.fastq",&seq_read, &seq_loaded,readlen_array_ptrs,
								     &bad_reads, quality_cut_off, &dup_reads, max_read_length, 
								     remove_duplicates, break_homopolymers, homopolymer_cutoff,
								     ascii_offset,
								     db_graph, individual_edge_array, 0);

  CU_ASSERT(seq_read==360);
  CU_ASSERT(seq_loaded==311);//just lose 1 read plus part of another due to Ns
  CU_ASSERT(db_graph->unique_kmers == 105);
  CU_ASSERT(dup_reads==0);
  hash_table_free(&db_graph);

  // then a test where you try to  remove duplicates from files that don't contain any

  dup_reads=0;
  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  remove_duplicates=true;
  seq_read=0;
  seq_loaded=0;
  for (j=0; j<=100; j++)
    {
      readlen_array[j]=0;
    }

  load_fastq_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/paired_end_file1_1.fastq",&seq_read, &seq_loaded,readlen_array_ptrs,
								     &bad_reads, quality_cut_off, &dup_reads, max_read_length, 
								     remove_duplicates, break_homopolymers, homopolymer_cutoff,ascii_offset,
								     db_graph, individual_edge_array, 0);



  CU_ASSERT(seq_read==360);
  CU_ASSERT(seq_loaded==311);
  CU_ASSERT(db_graph->unique_kmers == 105);
  CU_ASSERT(dup_reads==0);
  hash_table_free(&db_graph);
  dup_reads=0;


  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  seq_read=0;
  seq_loaded=0;
  for (j=0; j<=100; j++)
    {
      readlen_array[j]=0;
    }

  load_fastq_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/paired_end_file2_with_dup_1.fastq", &seq_read, &seq_loaded,readlen_array_ptrs,
								     &bad_reads, quality_cut_off, &dup_reads, max_read_length, 
								     remove_duplicates, break_homopolymers, homopolymer_cutoff, 
								     ascii_offset,
								     db_graph,  individual_edge_array, 0);
						   
  
  //as before, but with four more 36bp reads
  CU_ASSERT(seq_read==432);

  CU_ASSERT(seq_loaded==311);//basically same as previous fastq - the extra reads are ignored
  CU_ASSERT(dup_reads==2);

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
  long long seq_read=0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson1really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  
  //yes, only 4 bases in this fasta!
  CU_ASSERT(seq_read==4);
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

  seq_read=0;
  seq_loaded = 0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson2really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  CU_ASSERT(seq_read==4);
  CU_ASSERT(seq_loaded==0);
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


  seq_read=0;
  seq_loaded=0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson3really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  CU_ASSERT(seq_read==11);
  CU_ASSERT(seq_loaded==0);
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
  seq_read=0;
  seq_loaded=0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson4really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(seq_read==13);
  
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
  
  seq_read=0;
  seq_loaded=0;

  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson5really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  CU_ASSERT(seq_loaded==360);
  CU_ASSERT(seq_read==360);
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

  seq_read=0;
  seq_loaded=0;

  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson6really", &seq_read, &seq_loaded,&bad_reads, db_graph);
  CU_ASSERT(seq_read==68);
  CU_ASSERT(seq_loaded==67);//one N
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

  seq_read=0;
  seq_loaded=0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson7really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  CU_ASSERT(seq_read==72);
  CU_ASSERT(seq_loaded==67);//5 N's
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

  seq_read=0;
  seq_loaded=0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson8really", &seq_read, &seq_loaded,&bad_reads, db_graph);
  CU_ASSERT(seq_read==10);
  CU_ASSERT(seq_loaded==6);
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

  seq_read=0;
  seq_loaded=0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson9really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  CU_ASSERT(seq_read==22);
  CU_ASSERT(seq_loaded==5);//kmer is 5 and there are many Ns

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


  seq_read=0;
  seq_loaded = 0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson10really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  CU_ASSERT(seq_read==16);
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




  seq_read=0;
  seq_loaded = 0;
  load_population_as_fasta("../data/test/pop_graph/simplepeople_onlyperson11really", &seq_read, &seq_loaded, &bad_reads, db_graph);
  CU_ASSERT(seq_read==1140);
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




}



void test_align_next_read_to_graph_and_return_node_array()
{


  //set up db_graph 

  int kmer_size = 17;
  int number_of_bits = 10; 
  int bucket_size = 30;
  long long bad_reads = 0; long long dup_reads=0; 
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;

  long long seq_read=0;
  long long seq_loaded=0;
  dBGraph * db_graph;


  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/person3.fasta", &seq_read, &seq_loaded, &bad_reads,&dup_reads, 200, 
										  remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff, 
										  db_graph, individual_edge_array, 0);

  int max_read_length = 200;
  
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }


  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 



  //create file reader
   int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      printf("new_entry must be true in hsi test function");
      exit(1);
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }



  //Now we know person3.fasta is all in the graph in colour 0. Now let's just try to get our array of nodes:
   //person3 looks like this
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
>extra read
TTTTTTTTTTTTTTTTAAA

    */
  FILE* fp = fopen("../data/test/graph/person3.fasta", "r");
  if (fp==NULL)
    {
      printf("Cannot open ../data/test/graph/person3.fasta in test_align_next_read_to_graph_and_return_node_array");
      exit(1);
    }


  dBNode* array_nodes[200];//in fact there are 43 17-mers in the first line of the fasta
  Orientation array_or[200];
  int colour=0;
  int num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_read_length, array_nodes, array_or, true, file_reader,
								 seq, kmer_window, db_graph, colour);
  
  CU_ASSERT(num_kmers==43);
  BinaryKmer test_kmer, test_kmer_rev;
  seq_to_binary_kmer("TAACCCTAACCCTAACC", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[0]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[0]->kmer,test_kmer_rev)==true ); 

  seq_to_binary_kmer("AACCCTAACCCTAACCC", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[1]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[1]->kmer,test_kmer)==true ); 

  seq_to_binary_kmer("TAACCCTAACCCTAACC", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[42]==reverse);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[42]->kmer,test_kmer_rev)==true ); 
  

  //get the next read
  num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_read_length, array_nodes, array_or, true, file_reader,
								 seq, kmer_window, db_graph, colour);
  

  CU_ASSERT(num_kmers==45-17+1);//next read is 45 bases long

  seq_to_binary_kmer("ACCCTAACCCTAACCCT", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[0]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[0]->kmer,test_kmer)==true ); 

  seq_to_binary_kmer("CCCTAACCCTAACCCTA", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[1]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[1]->kmer,test_kmer)==true ); 

  seq_to_binary_kmer("CCTAACCCTAACCCTAA", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[2]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[2]->kmer,test_kmer)==true ); 

  seq_to_binary_kmer("CTAACCCTAACCCTAAC", 17, &test_kmer);
  binary_kmer_reverse_complement(&test_kmer, 17, &test_kmer_rev);
  CU_ASSERT(array_or[28]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(array_nodes[28]->kmer,test_kmer)==true ); 



  free(kmer_window->kmer);
  free(kmer_window);
  free_sequence(&seq);
  hash_table_free(&db_graph);
}



void test_read_next_variant_from_full_flank_file()
{


  //Load a single file containing a SNP between two reads, and use detect_vars and trusted SV caller to dump a pair of files.
  // then load these dumped files and check we get what we expect.

   //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits=6;
  int bucket_size   = 10;
  long long seq_read=0;
  long long seq_loaded=0;

  long long bad_reads = 0; 
  long long dup_reads=0;
  boolean remove_duplicates_single_endedly=false; 
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;


  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  
  int max_chunk_length=100;

  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/pop_graph/variations/one_person_with_SNP.fasta", &seq_read, &seq_loaded,
								     &bad_reads,&dup_reads,max_chunk_length,  
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
								     db_graph, individual_edge_array, 0);

  FILE* fout_bubble       = fopen("tmp_test.fff", "w");
  if (fout_bubble==NULL)
    {
      printf("Unable to open  tmp_test.fff ");
      exit(1);
    }

  int max_branch_len=10;
  void print_no_extra_info(VariantBranchesAndFlanks* var, FILE* fout)
  {
  }

  db_graph_detect_vars(fout_bubble, max_branch_len,db_graph, &detect_vars_condition_always_true, 
		       &db_node_action_set_status_visited, &db_node_action_set_status_visited, 
		       &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, &print_no_extra_info);
  fclose(fout_bubble);


  //Now check that we can read back in what we have printed out.
  dBNode* flank5p[20];
  dBNode* one_allele[20];
  dBNode* other_allele[20];
  dBNode* flank3p[20];
  Orientation flank5p_or[20];
  Orientation one_allele_or[20];
  Orientation other_allele_or[20];
  Orientation flank3p_or[20];
  int len_flank5p;
  int len_one_allele;
  int len_other_allele;
  int len_flank3p;
  int max_read_length = 50;

  FILE* var_fptr =  fopen("tmp_test_read_next_variant_from_full_flank_file_bubble.fff", "r");

  read_next_variant_from_full_flank_file(var_fptr, max_read_length,
					 flank5p, flank5p_or, &len_flank5p,
					 one_allele,  one_allele_or,  &len_one_allele,
					 other_allele,  other_allele_or,  &len_other_allele,
					 flank3p,     flank3p_or,     &len_flank3p,
					 db_graph, 0);
  
  CU_ASSERT( len_flank5p==3 );
  
  BinaryKmer tmp_kmer, tmp_kmer_rev;
  BinaryKmer tmp_kmer2, tmp_kmer_rev2;
  
  seq_to_binary_kmer("AGCTC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(flank5p[0]->kmer, tmp_kmer));
  CU_ASSERT(flank5p_or[0]==forward);
  
  seq_to_binary_kmer("GCTCA", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(flank5p[1]->kmer, tmp_kmer));
  CU_ASSERT(flank5p_or[1]==forward);
  
  seq_to_binary_kmer("CTCAT", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(flank5p[2]->kmer, tmp_kmer_rev));
  CU_ASSERT(flank5p_or[2]==reverse);
  
  //now the alleles
  CU_ASSERT(len_one_allele==2);
  CU_ASSERT(len_other_allele==2);
  //detect_vars will print alleles in alphabetical order of first base 
  seq_to_binary_kmer("AAGCG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("CAGCG", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(one_allele[0]->kmer, tmp_kmer) );
  CU_ASSERT(one_allele_or[0]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(other_allele[0]->kmer, tmp_kmer2) );
  CU_ASSERT(other_allele_or[0]==forward);
  
  seq_to_binary_kmer("AGCGC", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  seq_to_binary_kmer("AGCGC", 5, &tmp_kmer2);
  binary_kmer_reverse_complement(&tmp_kmer2, 5, &tmp_kmer_rev2);
  CU_ASSERT(binary_kmer_comparison_operator(one_allele[1]->kmer, tmp_kmer) );
  CU_ASSERT(one_allele_or[1]==forward);
  CU_ASSERT(binary_kmer_comparison_operator(other_allele[1]->kmer, tmp_kmer2) );
  CU_ASSERT(other_allele_or[1]==forward);
  
  //finally the 3p flank
  
  seq_to_binary_kmer("CTGCG", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(flank3p[0]->kmer, tmp_kmer_rev));
  CU_ASSERT(flank3p_or[0]==reverse);
  
  seq_to_binary_kmer("TGCGT", 5, &tmp_kmer);
  binary_kmer_reverse_complement(&tmp_kmer, 5, &tmp_kmer_rev);
  CU_ASSERT(binary_kmer_comparison_operator(flank3p[1]->kmer, tmp_kmer_rev));
  CU_ASSERT(flank3p_or[1]==reverse);
  



  hash_table_free(&db_graph);

}




void test_getting_readlength_distribution()
{
  //******************************************************************************************
  //1. Check that with single-ended parsing of fastq, when we get the read length distribution it correctly does so when there are
  //   low quality bases below threshold, homopolymers, Ns
  //******************************************************************************************

  

  int kmerlist[]={15,19,39};
  long long correct_answers_readlens[100][9];
  int p;
  for (p=0; p<=99; p++)
    {
      correct_answers_readlens[p][0]=0;
      correct_answers_readlens[p][1]=0;
      correct_answers_readlens[p][2]=0;
    }
  //for k=15, with no quality/homopol (ie jist due to Ns), we get reads of these lengths:41,  26, 17, 95, 41,39, 17
  correct_answers_readlens[41][0]=2;
  correct_answers_readlens[26][0]=1;
  correct_answers_readlens[17][0]=2;
  correct_answers_readlens[95][0]=1;
  correct_answers_readlens[39][0]=1;
  //for k=19, we get reads of tlength: 41, 11, 26, 17, 3, 95, 41,39, 17
  correct_answers_readlens[41][1]=2;
  correct_answers_readlens[26][1]=1;
  correct_answers_readlens[95][1]=1;
  correct_answers_readlens[39][1]=1;
  //for k=39, we get reads of tlength: 41,39,41, 95
  correct_answers_readlens[41][2]=2;
  correct_answers_readlens[95][2]=1;
  correct_answers_readlens[39][2]=1;


  // what about if we add a filter at q=10?

  for (p=0; p<=99; p++)
    {
      correct_answers_readlens[p][3]=0;
      correct_answers_readlens[p][4]=0;
      correct_answers_readlens[p][5]=0;
    }
  //for k=15, with q10 filte,r we get reads of lengths:  39, 18,17   41, 36    13, 16, 39 ,16
  correct_answers_readlens[39][3]=2;
  correct_answers_readlens[18][3]=1;
  correct_answers_readlens[17][3]=1;
  correct_answers_readlens[41][3]=1;
  correct_answers_readlens[36][3]=1;
  correct_answers_readlens[16][3]=2;
  //for k=19, with q10 filter, we get reads of lengths:  39, 41, 36 , 39
  correct_answers_readlens[39][4]=2;
  correct_answers_readlens[41][4]=1;
  correct_answers_readlens[36][4]=1;
  //for k=39, we get reads of lengths: 39,41,39
  correct_answers_readlens[39][5]=2;
  correct_answers_readlens[41][5]=1;


  // and now what if we add a homopolymer filter for homopolymers >=5 long?

  for (p=0; p<=99; p++)
    {
      correct_answers_readlens[p][6]=0;
      correct_answers_readlens[p][7]=0;
      correct_answers_readlens[p][8]=0;
    }
  //for k=15, with homopol 5 and q10 filter we get reads of lengths:  27, 11 , 18,17   12, 29, 36    13, 16, 39 ,16
  correct_answers_readlens[39][6]=1;
  correct_answers_readlens[27][6]=1;
  correct_answers_readlens[18][6]=1;
  correct_answers_readlens[17][6]=1;
  correct_answers_readlens[29][6]=1;
  correct_answers_readlens[36][6]=1;
  correct_answers_readlens[16][6]=2;
  //for k=19, with homopol5 and q10 filter, we get reads of lengths:  39, 41, 36 , 39
  correct_answers_readlens[27][7]=1;
  correct_answers_readlens[29][7]=1;
  correct_answers_readlens[36][7]=1;
  correct_answers_readlens[39][7]=1;
  //for k=39, with homopol5 and q10 filter we get reads of lengths: 39,41,39
  correct_answers_readlens[39][8]=1;






  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  int qfilter=0;
  
  void setup_hashtable_and_test(int j, int k_mer)
  {
    //prepare a hash table
      
    int number_of_bits = 10;
    int bucket_size = 100;
    int max_retries=82;

    dBGraph* db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,k_mer);
    
    long long dup_reads=0;
    long long bad_reads=0;
    long long seq_read, seq_loaded;
    int i;
    long long readlen_counts[100];
    long long* readlen_counts_ptrs[100];
    for (i=0; i<=99; i++)
      {
	readlen_counts[i]=0;
	readlen_counts_ptrs[i]=&readlen_counts[i];
      }
    if (j>2)
      {
	qfilter=10;
      }
    if (j>5)
      {
	break_homopolymers=true;
	homopolymer_cutoff=5;
      }
    //breakpoint:
    load_fastq_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/pop_graph/test_readlen_distrib.fastq", &seq_read, &seq_loaded, readlen_counts_ptrs,
									&bad_reads, qfilter, &dup_reads, 200, remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
								       33,db_graph, individual_edge_array, 0);
    //CU_ASSERT(seq_read==correct_answers_seq_read[j]);
    //    CU_ASSERT(seq_loaded==correct_answers_seq_loaded[j]);

    int n;

    for (n=0; n<=99; n++)
      {
	if (correct_answers_readlens[n][j] != *(readlen_counts_ptrs[n]) )
	  {
	    printf("n is %d and we expect %qd but get %qd\n", n, correct_answers_readlens[n][j], *(readlen_counts_ptrs[n]) );
	  }
	CU_ASSERT(correct_answers_readlens[n][j]==*(readlen_counts_ptrs[n]) );
      }
    

    hash_table_free(&db_graph);
  }
  


  setup_hashtable_and_test(0, kmerlist[0]);
  setup_hashtable_and_test(1,kmerlist[1]);
  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      setup_hashtable_and_test(2, kmerlist[2]);
    }
  setup_hashtable_and_test(3, kmerlist[0]);
  setup_hashtable_and_test(4, kmerlist[1]);
  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      setup_hashtable_and_test(5, kmerlist[2]);
    }

  setup_hashtable_and_test(6, kmerlist[0]);

  setup_hashtable_and_test(7, kmerlist[1]);
  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      setup_hashtable_and_test(8, kmerlist[2]);
    }

  


  //******************************************************************************************
  //2. Check that with PAIRED-END parsing of fastq, when we get the read length distribution it correctly does so when there are
  //   low quality bases below threshold, homopolymers, Ns
  //******************************************************************************************


  long long paired_end_answers[100][9];
  for (p=0; p<=99; p++)
    {
      paired_end_answers[p][0]=0;
      paired_end_answers[p][1]=0;
      paired_end_answers[p][2]=0;
      paired_end_answers[p][3]=0;
      paired_end_answers[p][4]=0;
      paired_end_answers[p][5]=0;
      paired_end_answers[p][6]=0;
      paired_end_answers[p][7]=0;
      paired_end_answers[p][8]=0;
    }

  // we are going to load 2 files. One is the same as the previous fastq we loaded, with an extra copy of the third read, and the other is its mate.
  // file1: read1, read2, read3, copy of read3
  // file2: read4, read5, copy of read5, copy of read5
  // so actually only the final pair will be removed by dup removal

  //first, what answers do you get with no quality or homopolymer filters, just the effects of Ns?
  //

  //for k=15, with no quality/homopol (ie jist due to Ns), we get reads of these lengths:41,  26, 17, 95, 41,39, 17 from lh read and 99,99,99 from the rh
  paired_end_answers[41][0]=3;
  paired_end_answers[26][0]=1;
  paired_end_answers[17][0]=3;
  paired_end_answers[95][0]=1;
  paired_end_answers[39][0]=2;
  paired_end_answers[99][0]=4;//from rh read
  //for k=19, we get reads of tlength: 41, 11, 26, 17, 3, 95, 41,39, 17
  paired_end_answers[41][1]=3;
  paired_end_answers[26][1]=1;
  paired_end_answers[95][1]=1;
  paired_end_answers[39][1]=2;
  paired_end_answers[99][1]=4;//from rh read
  //for k=39, we get reads of tlength: 41,39,41, 95
  paired_end_answers[41][2]=3;
  paired_end_answers[95][2]=1;
  paired_end_answers[39][2]=2;
  paired_end_answers[99][2]=4;//from rh read


  //for k=15, with q10 filte,r we get reads of lengths:  39, 18,17  |  41,36 |   39,16,16 |   39,16,16  and on rhs: 45,36,  99  99 99  
  paired_end_answers[39][3]=3;
  paired_end_answers[18][3]=1;
  paired_end_answers[16][3]=4;
  paired_end_answers[17][3]=1;
  paired_end_answers[41][3]=1;
  paired_end_answers[36][3]=2;//2nd one from rh read
  paired_end_answers[45][3]=1;//from rh read
  paired_end_answers[99][3]=3;//from rh read
  //for k=19, with q10 filter, we get reads of lengths:  39, 41, 36 , 39, 39 and on rhs: 45,36,99,99,99
  paired_end_answers[39][4]=3;
  paired_end_answers[41][4]=1;
  paired_end_answers[36][4]=2;//2nd one from rh read
  paired_end_answers[45][4]=1;//from rh read 
  paired_end_answers[99][4]=3;//from rh read 
  //for k=39, we get reads of lengths: 39,41,39,39  and on rhs: 45,99,99,99
  paired_end_answers[39][5]=3;
  paired_end_answers[41][5]=1;
  paired_end_answers[45][5]=1;//from rh read 
  paired_end_answers[99][5]=3;//from rh read 


  //and now with the homopolymer filter

  //for k=15, with homopol 5 and q10 filter we get reads of lengths:  27, 18,17 |    29, 36 |    16, 39 ,16  |  16,39,16, and on rhs: 36,45,99,99,99
  paired_end_answers[16][6]=4;
  paired_end_answers[17][6]=1;
  paired_end_answers[18][6]=1;
  paired_end_answers[27][6]=1;
  paired_end_answers[29][6]=1;
  paired_end_answers[36][6]=2;//2nd from rh read
  paired_end_answers[39][6]=2;

  paired_end_answers[45][6]=1;//from rh read
  paired_end_answers[99][6]=3;//from rh read   
  //for k=19, with homopol5 and q10 filter, we get reads of lengths:  39, 41, 36 , 39   and on rhs: 36,99
  paired_end_answers[27][7]=1;
  paired_end_answers[29][7]=1;
  paired_end_answers[36][7]=2;//2nd from rh read
  paired_end_answers[39][7]=2;
  paired_end_answers[45][7]=1;//from rh read 
  paired_end_answers[99][7]=3;//from rh read
  //for k=39, wuth homopol5 and q10 filter we get reads of lengths: 39,41,39   and on rhs: 45,99
  paired_end_answers[39][8]=2;
  paired_end_answers[45][8]=1;//from rh read 
  paired_end_answers[99][8]=3;//from rh read   


  boolean remove_duplicates_paired_endedly=false;
  break_homopolymers=false;
  homopolymer_cutoff=0;
  qfilter=0;
  
  void paired_setup_hashtable_and_test(int j, int k_mer)
  {
    //prepare a hash table
      
    int number_of_bits = 10;
    int bucket_size = 100;
    int max_retries=82;

    dBGraph* db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,k_mer);
    
    long long dup_reads=0;
    long long bad_reads=0;
    long long seq_read, seq_loaded;
    int i;
    long long readlen_counts[100];
    long long* readlen_counts_ptrs[100];
    for (i=0; i<=99; i++)
      {
	readlen_counts[i]=0;
	readlen_counts_ptrs[i]=&readlen_counts[i];
      }
    if (j>2)
      {
	qfilter=10;
      }
    if (j>5)
      {
	break_homopolymers=true;
	homopolymer_cutoff=5;
      }

    load_paired_end_data_from_filenames_into_graph_of_specific_person_or_pop("../data/test/pop_graph/test_readlen_distrib.pair1.fastq",
									     "../data/test/pop_graph/test_readlen_distrib.pair2.fastq",
									     FASTQ,
									     &seq_read, &seq_loaded, readlen_counts_ptrs,
									     &bad_reads, qfilter, 200, &dup_reads, remove_duplicates_paired_endedly, break_homopolymers, homopolymer_cutoff,
									     33,db_graph, individual_edge_array, 0);
    int n;
    for (n=0; n<=99; n++)
      {
	CU_ASSERT(paired_end_answers[n][j]==*(readlen_counts_ptrs[n]) );
	if (paired_end_answers[n][j] != *(readlen_counts_ptrs[n]) )
	  {
	    printf("n is %d and we expec %qd and we get %qd\n", n, paired_end_answers[n][j], *(readlen_counts_ptrs[n]) );
	  }
      }
     hash_table_free(&db_graph);
  }





  paired_setup_hashtable_and_test(0, kmerlist[0]);

  paired_setup_hashtable_and_test(1,kmerlist[1]);
  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      paired_setup_hashtable_and_test(2, kmerlist[2]);
    }
  
  paired_setup_hashtable_and_test(3, kmerlist[0]);
  paired_setup_hashtable_and_test(4, kmerlist[1]);
  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      paired_setup_hashtable_and_test(5, kmerlist[2]);
    }
     paired_setup_hashtable_and_test(6, kmerlist[0]);

  paired_setup_hashtable_and_test(7, kmerlist[1]);
  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      paired_setup_hashtable_and_test(8, kmerlist[2]);
    }
    



}


void test_loading_binary_data_iff_it_overlaps_a_fixed_colour()
{

  if (NUMBER_OF_COLOURS<=1)
    {
      printf("This test is redundant with only one colour\n");
      return;
    }

  int kmer_size = 3;
  int number_of_bits_pre = 4; 
  int number_of_bits_post = 8;
  int bucket_size = 5;
  long long bad_reads = 0; 
  int max_retries=10;
  int max_chunk_len_reading_from_fasta = 200;

  long long seq_length_parsed_pre=0;
  long long seq_length_loaded_pre=0;


  dBGraph * db_graph_pre;
  dBGraph * db_graph_post;

  db_graph_pre = hash_table_new(number_of_bits_pre,bucket_size,max_retries,kmer_size);

  //we need the following arguments for the API but we will not use them - for duplicate removal and homopolymer breaking
  boolean remove_duplicates_single_endedly=false;
  boolean break_homopolymers=false;
  long long dup_reads=0;
  int homopolymer_cutoff=0;


  //>read1
  //GATCGGGTGT
  //>read1 copy
  //GATCGGGTGT
  //>read2
  //GGCT
  //>read2 copy
  //GGCT
  //>read3
  //TAGG
  //>read3 copy
  //TAGG
  //>read4=read 1 with a single base error
  //GATCGGGAGT

  load_fasta_data_from_filename_into_graph_of_specific_person_or_pop("../data/test/graph/test_loading_binarynodes_if_overlap_a_colour.fasta", 
								     &seq_length_parsed_pre, &seq_length_loaded_pre, &bad_reads, &dup_reads, 20, 
								     remove_duplicates_single_endedly, break_homopolymers, homopolymer_cutoff,
								     db_graph_pre, individual_edge_array,0);

  CU_ASSERT(seq_length_parsed_pre==46);
  CU_ASSERT(seq_length_loaded_pre==46);

  GraphInfo ginfo;
  graph_info_initialise(&ginfo);
  graph_info_set_seq(&ginfo, 0, seq_length_parsed_pre);
  //dump a single colour binary just of this graph
  db_graph_dump_single_colour_binary_of_colour0("../data/test/pop_graph/dump_cortex_var_graph.singlecol.ctx", &db_node_condition_always_true, db_graph_pre, &ginfo);
  
  
  //OK, we now have dumped a binary corresponding to colour 0.
  //Now let's clean up, removing the bubble created by the single base error on the 2nd copy of read 4
  db_graph_remove_low_coverage_nodes_ignoring_colours(1, db_graph_pre);
  //and dump a clean graph,
  db_graph_dump_binary("../data/test/pop_graph/dump_cortex_var_graph.clean.ctx", &db_node_check_status_not_pruned, db_graph_pre, &ginfo);
  
  //cleanup, before starting all over
  hash_table_free(&db_graph_pre);
  CU_ASSERT(db_graph_pre==NULL);
  
  //Now load the clean graph, so the "dirty" nodes are not even there
  db_graph_post = hash_table_new(number_of_bits_post,bucket_size,10,kmer_size);
  int num_cols_in_binary=-1;
  int array_mean_readlens[NUMBER_OF_COLOURS];
  int* array_mean_readlens_ptrs[NUMBER_OF_COLOURS];
  long long array_total_seq[NUMBER_OF_COLOURS];
  long long* array_total_seq_ptrs[NUMBER_OF_COLOURS];

  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      array_mean_readlens[i]=0;
      array_mean_readlens_ptrs[i]=&array_mean_readlens[i];
      array_total_seq[i]=0;
      array_total_seq_ptrs[i]=&array_total_seq[i];
    }
  long long seq_length_post = load_multicolour_binary_from_filename_into_graph("../data/test/pop_graph/dump_cortex_var_graph.clean.ctx", db_graph_post, &num_cols_in_binary,
									       array_mean_readlens_ptrs, array_total_seq_ptrs);


  CU_ASSERT(num_cols_in_binary==NUMBER_OF_COLOURS);
  CU_ASSERT(array_mean_readlens[0]==0);
  CU_ASSERT(array_total_seq[0]==seq_length_parsed_pre);


  //OK, finally we are ready for our real test. Load the singlecolour uncleaned binary into colour1, but only load the bits that overlap the cleaned graph (ie colour 0)
  int mean_len;
  long long totseq;
  load_single_colour_binary_data_from_filename_into_graph("../data/test/pop_graph/dump_cortex_var_graph.singlecol.ctx",db_graph_post,&mean_len, &totseq,
							  false,individual_edge_array,1,true,0);





  BinaryKmer tmp_kmer1, tmp_kmer2;

  //all the kmers and their reverse complements from the cleaned graph - these should now be in colours 0 and 1
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GAT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ATC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("TCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("CGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("CGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("CCG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("GGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("ACC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("GGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ACC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("GTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("CAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element21 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element22 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element23 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element24 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element25 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element26 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


  //kmers that should not be in the graph, in either colour 0 or 1
  dBNode* test_element27 = hash_table_find(element_get_key(seq_to_binary_kmer("GGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element28 = hash_table_find(element_get_key(seq_to_binary_kmer("TCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element29 = hash_table_find(element_get_key(seq_to_binary_kmer("GAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element30 = hash_table_find(element_get_key(seq_to_binary_kmer("CTC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element31 = hash_table_find(element_get_key(seq_to_binary_kmer("AGT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);
  dBNode* test_element32 = hash_table_find(element_get_key(seq_to_binary_kmer("ACT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph_post);


  CU_ASSERT(test_element1 != NULL);
  CU_ASSERT(test_element2 != NULL);
  CU_ASSERT(test_element1 == test_element2);
  CU_ASSERT(db_node_get_coverage(test_element1,individual_edge_array,0)==6);
  CU_ASSERT(db_node_get_coverage(test_element1,individual_edge_array,1)==6);
  CU_ASSERT(get_edge_copy(*test_element1, individual_edge_array,0)==get_edge_copy(*test_element1, individual_edge_array,1));

  CU_ASSERT(test_element3 != NULL);
  CU_ASSERT(test_element4 != NULL);
  CU_ASSERT(test_element3 == test_element4);
  CU_ASSERT(db_node_get_coverage(test_element3,individual_edge_array,0)==3);
  CU_ASSERT(db_node_get_coverage(test_element3,individual_edge_array,1)==3);
  CU_ASSERT(get_edge_copy(*test_element3, individual_edge_array,0)==get_edge_copy(*test_element3, individual_edge_array,1));

  CU_ASSERT(test_element5 != NULL);
  CU_ASSERT(test_element6 != NULL);  
  CU_ASSERT(test_element5 == test_element6);
  CU_ASSERT(db_node_get_coverage(test_element5,individual_edge_array,0)==3);
  CU_ASSERT(db_node_get_coverage(test_element5,individual_edge_array,1)==3);
  CU_ASSERT(get_edge_copy(*test_element5, individual_edge_array,0)==get_edge_copy(*test_element5, individual_edge_array,1));

  CU_ASSERT(test_element7 != NULL);
  CU_ASSERT(test_element8 != NULL);
  CU_ASSERT(test_element7 == test_element8);
  CU_ASSERT(db_node_get_coverage(test_element7,individual_edge_array,0)==3);
  CU_ASSERT(db_node_get_coverage(test_element7,individual_edge_array,1)==3);
  CU_ASSERT(get_edge_copy(*test_element7, individual_edge_array,0)==get_edge_copy(*test_element7, individual_edge_array,1));

  CU_ASSERT(test_element9 != NULL);
  CU_ASSERT(test_element10 != NULL);
  CU_ASSERT(test_element9 == test_element10);  
  CU_ASSERT(db_node_get_coverage(test_element9,individual_edge_array,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element9,individual_edge_array,1)==2);
  CU_ASSERT(get_edge_copy(*test_element9, individual_edge_array,0)==get_edge_copy(*test_element9, individual_edge_array,1));

  CU_ASSERT(test_element11 != NULL);
  CU_ASSERT(test_element12 != NULL);
  CU_ASSERT(test_element11 == test_element12);  
  CU_ASSERT(db_node_get_coverage(test_element11,individual_edge_array,0)==3);
  CU_ASSERT(db_node_get_coverage(test_element11,individual_edge_array,1)==3);
  CU_ASSERT(get_edge_copy(*test_element11, individual_edge_array,0)==get_edge_copy(*test_element11, individual_edge_array,1));


  CU_ASSERT(test_element13 != NULL);
  CU_ASSERT(test_element14 != NULL);
  CU_ASSERT(test_element13 == test_element14);  
  CU_ASSERT(db_node_get_coverage(test_element13,individual_edge_array,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element13,individual_edge_array,1)==2);
  CU_ASSERT(get_edge_copy(*test_element13, individual_edge_array,0)==get_edge_copy(*test_element13, individual_edge_array,1));

  CU_ASSERT(test_element15 != NULL);
  CU_ASSERT(test_element16 != NULL);
  CU_ASSERT(test_element15 == test_element16);  
  CU_ASSERT(db_node_get_coverage(test_element15,individual_edge_array,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element15,individual_edge_array,1)==2);
  CU_ASSERT(get_edge_copy(*test_element15, individual_edge_array,0)==get_edge_copy(*test_element15, individual_edge_array,1));

  CU_ASSERT(test_element17 != NULL);
  CU_ASSERT(test_element18 != NULL);
  CU_ASSERT(test_element17 == test_element18);  
  CU_ASSERT(db_node_get_coverage(test_element17,individual_edge_array,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element17,individual_edge_array,1)==2);
  CU_ASSERT(get_edge_copy(*test_element17, individual_edge_array,0)==get_edge_copy(*test_element17, individual_edge_array,1));

  CU_ASSERT(test_element19 != NULL);
  CU_ASSERT(test_element20 != NULL);
  CU_ASSERT(test_element19 == test_element20);  
  CU_ASSERT(db_node_get_coverage(test_element19,individual_edge_array,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element19,individual_edge_array,1)==2);
  CU_ASSERT(get_edge_copy(*test_element19, individual_edge_array,0)==get_edge_copy(*test_element19, individual_edge_array,1));

  CU_ASSERT(test_element21 != NULL);
  CU_ASSERT(test_element22 != NULL);
  CU_ASSERT(test_element21 == test_element22);  
  CU_ASSERT(db_node_get_coverage(test_element21,individual_edge_array,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element21,individual_edge_array,1)==2);
  CU_ASSERT(get_edge_copy(*test_element21, individual_edge_array,0)==get_edge_copy(*test_element21, individual_edge_array,1));


  CU_ASSERT(test_element23 != NULL);
  CU_ASSERT(test_element24 != NULL);
  CU_ASSERT(test_element23 == test_element24);  
  CU_ASSERT(db_node_get_coverage(test_element23,individual_edge_array,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element23,individual_edge_array,1)==2);
  CU_ASSERT(get_edge_copy(*test_element23, individual_edge_array,0)==get_edge_copy(*test_element23, individual_edge_array,1));

  CU_ASSERT(test_element25 != NULL);
  CU_ASSERT(test_element26 != NULL);
  CU_ASSERT(test_element25 == test_element26);  
  CU_ASSERT(db_node_get_coverage(test_element25,individual_edge_array,0)==2);
  CU_ASSERT(db_node_get_coverage(test_element25,individual_edge_array,1)==2);
  CU_ASSERT(get_edge_copy(*test_element25, individual_edge_array,0)==get_edge_copy(*test_element25, individual_edge_array,1));



  //these nodes should just not be there
  CU_ASSERT(test_element27 == NULL);
  CU_ASSERT(test_element28 == NULL);
  CU_ASSERT(test_element29 == NULL);
  CU_ASSERT(test_element31 == NULL);
  CU_ASSERT(test_element32 == NULL);


  hash_table_free(&db_graph_post);
  CU_ASSERT(db_graph_post == NULL);
 
}
