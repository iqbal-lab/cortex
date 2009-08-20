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
#include <limits.h>

void test_db_graph_supernode_for_specific_person_or_pop()
{


 //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 8;
  int bucket_size    = 8;
  long long bad_reads = 0;
  int max_retries=10;

  
  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  int seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_one_long_supernode_with_conflict_at_end", &bad_reads, hash_table);
  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(bad_reads==0);

  //we have loaded a fasta containing these supernodes ACATT, TTC, TTG
  // reads are: ACATT, ATTC, ATTG
  //note we have loaded one person only
 

  int max_expected_supernode_length=100;
  dBNode * nodes_path[max_expected_supernode_length];
  Orientation orientations_path[max_expected_supernode_length];
  Nucleotide labels_path[max_expected_supernode_length];
  char seq[max_expected_supernode_length+hash_table->kmer_size+1];

  //get the supernodes one at a time. Remember length is the number of edges between nodes.

  dBNode* test_elem1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size),kmer_size), hash_table);
  dBNode* test_elem2 = hash_table_find(element_get_key(seq_to_binary_kmer("CAT", kmer_size),kmer_size), hash_table);
  dBNode* test_elem3 = hash_table_find(element_get_key(seq_to_binary_kmer("ATT", kmer_size),kmer_size), hash_table);
  CU_ASSERT(!(test_elem1 == NULL));
  CU_ASSERT(!(test_elem2 == NULL));
  CU_ASSERT(!(test_elem3 == NULL));


  CU_ASSERT(db_node_check_status(test_elem1,none));
  CU_ASSERT(db_node_check_status(test_elem2,none));
  CU_ASSERT(db_node_check_status(test_elem3,none));

  //coverage variables
  double avg_coverage=0;
  int max_coverage=0;
  int min_coverage=0;
  boolean is_cycle=false;




  int length = db_graph_supernode_for_specific_person_or_pop(test_elem1,max_expected_supernode_length, &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							     nodes_path,orientations_path, labels_path, seq,
							     &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							     hash_table, individual_edge_array, 0); //0 because we are looking at person 0 - the only person in the graph

  CU_ASSERT(length==2);
  CU_ASSERT(avg_coverage==1);
  CU_ASSERT(min_coverage==1);
  CU_ASSERT(max_coverage==1); //remember coverage only measured on internal nodes of a supernode.
  CU_ASSERT(is_cycle==false);

  CU_ASSERT_STRING_EQUAL(seq, "TT"); //does not put initial kmer "ACA" into the string, just the edges
  CU_ASSERT(nodes_path[0]==test_elem1);
  CU_ASSERT(nodes_path[1]==test_elem2);
  CU_ASSERT(nodes_path[2]==test_elem3);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==reverse);
  CU_ASSERT(orientations_path[2]==reverse);
  CU_ASSERT(labels_path[0]==Thymine);
  CU_ASSERT(labels_path[1]==Thymine);
	     
  CU_ASSERT(db_node_check_status(test_elem1, visited));
  CU_ASSERT(db_node_check_status(test_elem2, visited));
  CU_ASSERT(db_node_check_status(test_elem3, visited));







  hash_table_free(&hash_table);


  // ******** try another example ******************
  

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.txt", &bad_reads, hash_table);

  //person 1:
  //>person1_read1
  //AAAAAAAAAAAAAAA
  //>person1_read2
  //ACGTT     <<<<<<<<<<<<<<<<Note this generates supernode (AAC)GTT due to a hairpin
  //person 2:
  //>person2_read1
  //GGGGGGGGGGGGGGGGGG
  //>person2_read2
  //TTGACG   <<<<<<<<<<<<<<<< Also has a hairpin

  test_elem1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size),kmer_size), hash_table);
  test_elem2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACG", kmer_size),kmer_size), hash_table);
  test_elem3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGT", kmer_size),kmer_size), hash_table);
  dBNode* test_elem4 = hash_table_find(element_get_key(seq_to_binary_kmer("GTT", kmer_size),kmer_size), hash_table);
  dBNode* test_elem5 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size),kmer_size), hash_table);

  dBNode* test_elem6 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG", kmer_size),kmer_size), hash_table);

  dBNode* test_elem7 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA", kmer_size),kmer_size), hash_table);
  dBNode* test_elem8 = hash_table_find(element_get_key(seq_to_binary_kmer("GAC", kmer_size),kmer_size), hash_table);
  dBNode* test_elem9 = hash_table_find(element_get_key(seq_to_binary_kmer("ACG", kmer_size),kmer_size), hash_table);


  CU_ASSERT(!(test_elem1 == NULL));
  CU_ASSERT(!(test_elem2 == NULL));
  CU_ASSERT(!(test_elem3 == NULL));
  CU_ASSERT(!(test_elem4 == NULL));
  CU_ASSERT(!(test_elem5 == NULL));
  CU_ASSERT(!(test_elem6 == NULL));
  CU_ASSERT(!(test_elem7 == NULL));
  CU_ASSERT(!(test_elem8 == NULL));
  CU_ASSERT(!(test_elem9 == NULL));


  avg_coverage=0;
  min_coverage=0;
  max_coverage=0;
  is_cycle=false;




  length = db_graph_supernode_for_specific_person_or_pop(test_elem1,max_expected_supernode_length, &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, individual_edge_array, 0); //person 0 int his array is person in person1.fasta - sorry, offset by 1


  CU_ASSERT(length==1); //just one edge - self loop
  CU_ASSERT_STRING_EQUAL(seq, "A");
  CU_ASSERT(nodes_path[0]==test_elem1);
  CU_ASSERT(nodes_path[1]==test_elem1);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==forward);
  CU_ASSERT(labels_path[0]==Adenine);

  CU_ASSERT(avg_coverage==0);
  CU_ASSERT(min_coverage==0);
  CU_ASSERT(max_coverage==0);//there are not interior nodes of supernode
  CU_ASSERT(is_cycle==true);

  CU_ASSERT(db_node_check_status(test_elem1, visited)==true); 
  CU_ASSERT(db_node_check_status(test_elem2, none)==true);
  CU_ASSERT(db_node_check_status(test_elem3, none)==true);
  CU_ASSERT(db_node_check_status(test_elem4, none)==true);
  CU_ASSERT(db_node_check_status(test_elem5, none)==true);
  CU_ASSERT(db_node_check_status(test_elem6, none)==true);
  CU_ASSERT(db_node_check_status(test_elem7, none)==true);
  CU_ASSERT(db_node_check_status(test_elem8, none)==true);
  CU_ASSERT(db_node_check_status(test_elem9, none)==true);

  db_graph_set_all_visited_nodes_to_status_none(hash_table);
  CU_ASSERT(db_node_check_status(test_elem1, none)==true);


  avg_coverage=0;
  min_coverage=0;
  max_coverage=0;
  is_cycle=false;

  length = db_graph_supernode_for_specific_person_or_pop(test_elem2,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, individual_edge_array, 0); //person 0 int his array is person in person1.fasta - sorry, offset by 1


  CU_ASSERT(length==3);
  CU_ASSERT_STRING_EQUAL(seq, "GTT");
  CU_ASSERT(nodes_path[0]==test_elem4);
  CU_ASSERT(nodes_path[1]==test_elem2);
  CU_ASSERT(nodes_path[2]==test_elem2);//hairpin!
  CU_ASSERT(nodes_path[3]==test_elem4);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==forward);
  CU_ASSERT(orientations_path[2]==reverse);
  CU_ASSERT(orientations_path[3]==reverse);
  CU_ASSERT(labels_path[0]==Guanine);
  CU_ASSERT(labels_path[1]==Thymine);
  CU_ASSERT(labels_path[2]==Thymine);
  CU_ASSERT(avg_coverage==2);
  CU_ASSERT(min_coverage==2);
  CU_ASSERT(max_coverage==2);//note that the single interior node ALSO exists in the graph of person2, but that should not incrememtn this coverage to 3.
  CU_ASSERT(is_cycle==false);
  
  CU_ASSERT(db_node_check_status(test_elem1, none)==true);
  CU_ASSERT(db_node_check_status(test_elem2, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem3, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem4, visited)==true); 
  CU_ASSERT(db_node_check_status(test_elem5, none)==true);
  CU_ASSERT(db_node_check_status(test_elem6, none)==true);
  CU_ASSERT(db_node_check_status(test_elem7, none)==true);
  CU_ASSERT(db_node_check_status(test_elem8, none)==true);
  CU_ASSERT(db_node_check_status(test_elem9, visited)==true);

  
  //now try again with a different initial node within the same supernode. Should return the same as before.
  //In general it doesn't guarantee to trvaerse the supernode in the same direction; we know which was it will go however

  avg_coverage=0;
  min_coverage=0;
  max_coverage=0;
  is_cycle=false;

  length = db_graph_supernode_for_specific_person_or_pop(test_elem3,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, individual_edge_array, 0); //person 0 int his array is person in person1.fasta - sorry, offset by 1
  

  CU_ASSERT(length==3);
  CU_ASSERT_STRING_EQUAL(seq, "GTT");
  CU_ASSERT(nodes_path[0]==test_elem4);
  CU_ASSERT(nodes_path[1]==test_elem2);
  CU_ASSERT(nodes_path[2]==test_elem2);//hairpin!
  CU_ASSERT(nodes_path[3]==test_elem4);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==forward);
  CU_ASSERT(orientations_path[2]==reverse);
  CU_ASSERT(orientations_path[3]==reverse);
  CU_ASSERT(labels_path[0]==Guanine);
  CU_ASSERT(labels_path[1]==Thymine);
  CU_ASSERT(labels_path[2]==Thymine);
  CU_ASSERT(avg_coverage==2);
  CU_ASSERT(min_coverage==2);
  CU_ASSERT(max_coverage==2);//note that the single interior node ALSO exists in the graph of person2, but that should not incrememtn this coverage to 3.
  CU_ASSERT(is_cycle==false);
  
  CU_ASSERT(db_node_check_status(test_elem1, none)==true);
  CU_ASSERT(db_node_check_status(test_elem2, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem3, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem4, visited)==true); 
  CU_ASSERT(db_node_check_status(test_elem5, none)==true);
  CU_ASSERT(db_node_check_status(test_elem6, none)==true);
  CU_ASSERT(db_node_check_status(test_elem7, none)==true);
  CU_ASSERT(db_node_check_status(test_elem8, none)==true);
  CU_ASSERT(db_node_check_status(test_elem9, visited)==true);




  length = db_graph_supernode_for_specific_person_or_pop(test_elem4,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, individual_edge_array, 0); //person 0 int his array is person in person1.fasta - sorry, offset by 1
  
  
  CU_ASSERT(length==3);
  CU_ASSERT_STRING_EQUAL(seq, "GTT");
  CU_ASSERT(nodes_path[0]==test_elem4);
  CU_ASSERT(nodes_path[1]==test_elem2);
  CU_ASSERT(nodes_path[2]==test_elem2);//hairpin!
  CU_ASSERT(nodes_path[3]==test_elem4);
  CU_ASSERT(orientations_path[0]==forward);
  CU_ASSERT(orientations_path[1]==forward);
  CU_ASSERT(orientations_path[2]==reverse);
  CU_ASSERT(orientations_path[3]==reverse);
  CU_ASSERT(labels_path[0]==Guanine);
  CU_ASSERT(labels_path[1]==Thymine);
  CU_ASSERT(labels_path[2]==Thymine);
  CU_ASSERT(avg_coverage==2);
  CU_ASSERT(min_coverage==2);
  CU_ASSERT(max_coverage==2);//note that the single interior node ALSO exists in the graph of person2, but that should not incrememtn this coverage to 3.
  CU_ASSERT(is_cycle==false);
  
  CU_ASSERT(db_node_check_status(test_elem1, none)==true);
  CU_ASSERT(db_node_check_status(test_elem2, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem3, visited)==true);
  CU_ASSERT(db_node_check_status(test_elem4, visited)==true); 
  CU_ASSERT(db_node_check_status(test_elem5, none)==true);
  CU_ASSERT(db_node_check_status(test_elem6, none)==true);
  CU_ASSERT(db_node_check_status(test_elem7, none)==true);
  CU_ASSERT(db_node_check_status(test_elem8, none)==true);
  CU_ASSERT(db_node_check_status(test_elem9, visited)==true);



  //and now a trap! This node is not there in person 1, but is in person 2
  length = db_graph_supernode_for_specific_person_or_pop(test_elem6,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, individual_edge_array, 0);   
  CU_ASSERT(length==0);



  length = db_graph_supernode_for_specific_person_or_pop(test_elem6,max_expected_supernode_length,&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
							 nodes_path,orientations_path, labels_path, seq,
							 &avg_coverage, &min_coverage, &max_coverage, &is_cycle,
							 hash_table, individual_edge_array, 1); 
  CU_ASSERT(length==3);
  CU_ASSERT_STRING_EQUAL(seq, "CAA");  


  

  hash_table_free(&hash_table);
}


void test_is_supernode_end()
{

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 8;
  int bucket_size    = 8;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a graph/hash table.
  //         The test then picks an element in the graph, and calls get_seq_from_elem_to_end_of_supernode
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****

  int seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_two_self_loops",  &bad_reads, hash_table);
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
  number_of_bits = 4;
  bucket_size    = 4;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_one_long_supernode_with_conflict_at_end", &bad_reads, hash_table);

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
  number_of_bits = 4;
  bucket_size    = 4;
  bad_reads = 0;
  max_retries=10;


  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);


  seq_loaded=load_population_as_fasta("../data/test/pop_graph/supernode/one_person_one_long_supernode_with_inward_conflict_at_end",&bad_reads, hash_table);

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
  number_of_bits = 4;
  bucket_size    = 4;
  bad_reads = 0;
  max_retries=10;


  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_infiniteloop",  &bad_reads, hash_table);
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


  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.txt", &bad_reads, hash_table);
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


void test_get_min_and_max_covg_of_nodes_in_supernode()
{

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  long long seq_loaded=0;



  //start with an example with just oner person - note this is an example with the same kmer occuring twice in one read :) - CCG
  // so need to make sure that coverga eof that kmer is not double-incremented for the two times in that read

  // >r1
  //AACCGG
  //>r2
  //AAC
  //>r3
  //AAC
  //>r4
  //AAC
  //>r5
  //ACC
  //>r6
  //CCG
  //>r7
  //CCG

  //supernode will be
  // AAC --> ACC --> CCG ___|
  // TTG <-- TGG <-- GGC <--|
  // hairpin :-)

  //double coverness of hairpin must not confuse read coverage


  seq_loaded = load_population_as_fasta("../data/test/pop_graph/coverage/one_person", &bad_reads, hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 24);
  CU_ASSERT(bad_reads==0);
  

  dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_get_coverage(query_node, individual_edge_array, 0)==4);
  //printf("Expect covg aac is 4, but is %d", db_node_get_coverage(query_node, individual_edge_array, 0));
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("ACC",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_get_coverage(query_node, individual_edge_array, 0)==2);
  //printf("Expect covg acc is 2, but is %d", db_node_get_coverage(query_node, individual_edge_array, 0));
  query_node = hash_table_find(element_get_key(seq_to_binary_kmer("CCG",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node==NULL));
  CU_ASSERT(db_node_get_coverage(query_node, individual_edge_array, 0)==4); //NOTE this is only 4 because we cound kmer coverage. CCG occurs twice in read r1, once on each strand
  //printf("Expect covg ccg is 4, but is %d", db_node_get_coverage(query_node, individual_edge_array, 0));


  hash_table_free(&hash_table);

}


void test_db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta()
{

  int length_of_arrays=14;
  int number_of_nodes_to_load=7;


  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (db_graph==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/one_person_for_testing_array_loading", &bad_reads, db_graph);

  //>one read
  //AATAGACGCCCACACCTGATAGACCCCACACTCTAA

  
  CU_ASSERT(seq_loaded==36);
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/one_person.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ./data/test/pop_graph/one_person.fasta\n");
      exit(1);
    }

  //*************************************
  // malloc and initialising
  //*************************************
  dBNode**     chrom_path_array        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); //everything these dBNode*'s point to will be in the hash table - ie owned by the hash table
  Orientation* chrom_orientation_array = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
  Nucleotide*  chrom_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
  char*        chrom_string            = (char*) malloc(sizeof(char)*length_of_arrays+1); //+1 for \0

  int n;
  for (n=0; n<length_of_arrays; n++)
    {
      chrom_path_array[n]=NULL;
      chrom_orientation_array[n]=forward;
      chrom_labels[n]=Undefined;
      chrom_string[n]='N';
    }
  chrom_string[0]='\0';


  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,number_of_nodes_to_load+db_graph->kmer_size+1,LINE_MAX);
  seq->seq[0]='\0';


  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window . Exit.\n");
      exit(1);
    }
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*1000);    //*(number_of_nodes_to_load + db_graph->kmer_size));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer ");
      exit(1);
    }
  
  kmer_window->nkmers=0;

  // ***********************************************************
  //********** end of malloc and initialise ********************



  int ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, 0, 
													     length_of_arrays,
													     chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													     seq, kmer_window, 
													     true, false,
													     db_graph);
  CU_ASSERT(ret==number_of_nodes_to_load);

  char tmp_seqzam[db_graph->kmer_size+1];
  tmp_seqzam[db_graph->kmer_size]='\0';


  int i;
  //length of arrays is 14, number of nodes to load is 7
  for(i=0; i<length_of_arrays-number_of_nodes_to_load; i++)
    {
      CU_ASSERT(chrom_path_array[i]==NULL);
    }
  CU_ASSERT_STRING_EQUAL( "AATAGAC", binary_kmer_to_seq(chrom_path_array[7]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ATAGACG", binary_kmer_to_seq(chrom_path_array[8]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GCGTCTA", binary_kmer_to_seq(chrom_path_array[9]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGACGCC", binary_kmer_to_seq(chrom_path_array[10]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GACGCCC", binary_kmer_to_seq(chrom_path_array[11]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACGCCCA", binary_kmer_to_seq(chrom_path_array[12]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CGCCCAC", binary_kmer_to_seq(chrom_path_array[13]->kmer, db_graph->kmer_size, tmp_seqzam) );


  //one more batch, then array is full,
  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, false,
													 db_graph);
  CU_ASSERT(ret==number_of_nodes_to_load);



  CU_ASSERT_STRING_EQUAL( "AATAGAC", binary_kmer_to_seq(chrom_path_array[0]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ATAGACG", binary_kmer_to_seq(chrom_path_array[1]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GCGTCTA", binary_kmer_to_seq(chrom_path_array[2]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGACGCC", binary_kmer_to_seq(chrom_path_array[3]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GACGCCC", binary_kmer_to_seq(chrom_path_array[4]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACGCCCA", binary_kmer_to_seq(chrom_path_array[5]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CGCCCAC", binary_kmer_to_seq(chrom_path_array[6]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GCCCACA", binary_kmer_to_seq(chrom_path_array[7]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCACAC", binary_kmer_to_seq(chrom_path_array[8]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCACACC", binary_kmer_to_seq(chrom_path_array[9]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGGTGTG", binary_kmer_to_seq(chrom_path_array[10]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACACCTG", binary_kmer_to_seq(chrom_path_array[11]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CACCTGA", binary_kmer_to_seq(chrom_path_array[12]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACCTGAT", binary_kmer_to_seq(chrom_path_array[13]->kmer, db_graph->kmer_size, tmp_seqzam) );
  

  //from now on, it is always true that LAST time was not a new fasta entry, so penultimate argument is TRUE
  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, true,
                                                                                                         db_graph);

  CU_ASSERT(ret==number_of_nodes_to_load);


  CU_ASSERT_STRING_EQUAL( "GCCCACA", binary_kmer_to_seq(chrom_path_array[0]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCACAC", binary_kmer_to_seq(chrom_path_array[1]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCACACC", binary_kmer_to_seq(chrom_path_array[2]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGGTGTG", binary_kmer_to_seq(chrom_path_array[3]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACACCTG", binary_kmer_to_seq(chrom_path_array[4]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CACCTGA", binary_kmer_to_seq(chrom_path_array[5]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACCTGAT", binary_kmer_to_seq(chrom_path_array[6]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCTGATA", binary_kmer_to_seq(chrom_path_array[7]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CTATCAG", binary_kmer_to_seq(chrom_path_array[8]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "TCTATCA", binary_kmer_to_seq(chrom_path_array[9]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GATAGAC", binary_kmer_to_seq(chrom_path_array[10]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ATAGACC", binary_kmer_to_seq(chrom_path_array[11]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GGGTCTA", binary_kmer_to_seq(chrom_path_array[12]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGACCCC", binary_kmer_to_seq(chrom_path_array[13]->kmer, db_graph->kmer_size, tmp_seqzam) );
  

  //and again

  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, true,
                                                                                                         db_graph);

  CU_ASSERT(ret==number_of_nodes_to_load);


  CU_ASSERT_STRING_EQUAL( "CCTGATA", binary_kmer_to_seq(chrom_path_array[0]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CTATCAG", binary_kmer_to_seq(chrom_path_array[1]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "TCTATCA", binary_kmer_to_seq(chrom_path_array[2]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GATAGAC", binary_kmer_to_seq(chrom_path_array[3]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ATAGACC", binary_kmer_to_seq(chrom_path_array[4]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "GGGTCTA", binary_kmer_to_seq(chrom_path_array[5]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGACCCC", binary_kmer_to_seq(chrom_path_array[6]->kmer, db_graph->kmer_size, tmp_seqzam) );

  CU_ASSERT_STRING_EQUAL( "GACCCCA", binary_kmer_to_seq(chrom_path_array[7]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACCCCAC", binary_kmer_to_seq(chrom_path_array[8]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCCACA", binary_kmer_to_seq(chrom_path_array[9]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCACAC", binary_kmer_to_seq(chrom_path_array[10]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGTGTGG", binary_kmer_to_seq(chrom_path_array[11]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CACACTC", binary_kmer_to_seq(chrom_path_array[12]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACACTCT", binary_kmer_to_seq(chrom_path_array[13]->kmer, db_graph->kmer_size, tmp_seqzam) );


  //now - does it cope with hitting the end of the entry before getting the required number of nodes

  ret = db_graph_load_array_with_next_batch_of_nodes_corresponding_to_consecutive_bases_in_a_chrom_fasta(chrom_fptr, number_of_nodes_to_load, number_of_nodes_to_load, 
													 length_of_arrays,
													 chrom_path_array, chrom_orientation_array, chrom_labels, chrom_string,
													 seq, kmer_window, 
													 false, true,
                                                                                                         db_graph);

  CU_ASSERT(ret==2);

  CU_ASSERT_STRING_EQUAL( "GACCCCA", binary_kmer_to_seq(chrom_path_array[0]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACCCCAC", binary_kmer_to_seq(chrom_path_array[1]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCCACA", binary_kmer_to_seq(chrom_path_array[2]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CCCACAC", binary_kmer_to_seq(chrom_path_array[3]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "AGTGTGG", binary_kmer_to_seq(chrom_path_array[4]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "CACACTC", binary_kmer_to_seq(chrom_path_array[5]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACACTCT", binary_kmer_to_seq(chrom_path_array[6]->kmer, db_graph->kmer_size, tmp_seqzam) );

  CU_ASSERT_STRING_EQUAL( "CACTCTA", binary_kmer_to_seq(chrom_path_array[7]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT_STRING_EQUAL( "ACTCTAA", binary_kmer_to_seq(chrom_path_array[8]->kmer, db_graph->kmer_size, tmp_seqzam) );
  CU_ASSERT(chrom_path_array[9]==NULL);
  CU_ASSERT(chrom_path_array[10]==NULL);
  CU_ASSERT(chrom_path_array[11]==NULL);
  CU_ASSERT(chrom_path_array[12]==NULL);
  CU_ASSERT(chrom_path_array[13]==NULL);



  //cleanup

  free(kmer_window->kmer);
  free(kmer_window);
  free(chrom_path_array);
  free(chrom_orientation_array);
  free(chrom_string);
  free(chrom_labels);
  free_sequence(&seq);
  hash_table_free(&db_graph);
  fclose(chrom_fptr);

}



  //A fundamental assumption in our implementation is that the fasta whose path you are following is long, much longer than the longest supernode
  // you are going to find. This is fine when running using reference chromosomes. However in tests, this is a pain, and so I am using fasta's for these tests
  // that are padded at the end with N's.



/* this has been split into separate tests

void test_db_graph_make_reference_path_based_sv_calls()
{


  // ******************************************************************************************************************************
  // 1. NULL test. Load short fake reference, and load a person whose sequence is identical. This should find nothing.
  // ******************************************************************************************************************************


  // toy example  - kmer_size=7
  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;


  
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/one_person_with_Ns_on_end", &bad_reads, hash_table);

  //>one read
  //AATAGACGCCCACACCTGATAGACCCCACAC 


  // The de Bruijn graph of this in 7mers is as follows:
  // a line of 8 nodes, with a loop of length 16 off the 9th node. ie all nodes have 1-in and 1-out, except the 9th node, CCCACAC, which has 2 ins and 2 outs.
  // so it's a supernode of length 8, a loop-supernode of length 16,

  // our algorithm should look at the supernode starting with AATAGAC, then the supernode starting at CCCACAC, then stop


  
  CU_ASSERT(seq_loaded==1462);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/first_person_with_one_read_and_Ns_on_end.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ./data/test/pop_graph/first_person_with_one_read_and_Ns_on_end.fasta\n");
      exit(1);
    }

  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= 3;// I want to try to attach anchors, and the supernodes are quite short, so for this test only ask for (ridiculously) small anchors
  int max_anchor_span = 20;
  int length_of_arrays=40;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=20;



  
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
							individual_edge_array,1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
							0, NULL, NULL, NULL, NULL, NULL);

  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);





  // ******************************************************************************************************************************
  // 2. Harder NULL test. Reference=an ALU and load a person whose sequence is also that ALU. This should find nothing.
  //     Our implementation of db_graph_make_reference_path_based_sv_calls assumes we have a very large ref fasta, so we load big arrays of bases
  //     which are much longer than we expect any of the supernodes to be. In this case, we know the supernode is as long as the Alu.  
  // ******************************************************************************************************************************

  kmer_size = 31;
  number_of_bits = 8;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;

  //just load one person who's sequence is an Alu. Then take reference which is that same sequence and try to find variants
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/just_one_of_the_two_people.txt", &bad_reads, hash_table);

  //   >7SLRNA#SINE/Alu  plus GTTCAGAG at start and GTCAGCGTAG at end
  //   GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  //   AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG
  // plus loads of N's on end
   
  
  CU_ASSERT(seq_loaded==5158);
  
  chrom_fptr = fopen("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/person1.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/person1.fasta");
      exit(1);
    }

  min_fiveprime_flank_anchor = 10;
  min_threeprime_flank_anchor= 10;
  max_anchor_span = 370;
  length_of_arrays=740;
  min_covg =1;
  max_covg = 10;
  max_expected_size_of_supernode=370;
  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
						    0, NULL, NULL, NULL, NULL, NULL);
  
  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);



  // ******************************************************************************************************************************
  // 3. Reference = Alu-NNNNN- same Alu, and person is identical to reference. This should find nothing
  //    Note this really is testing something different. Since the reference is basically the same sequence repeated twice (with N's in between),
  //    there is a risk that we would match the 5' anchor of the (one and only) supernode at the start of the reference, and the 3' anchor at
  //    the end of the 2nd copy of the repeat.
  //    To be clear, look at this:
  //        Ref:  ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //        Indiv ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //      So we look at supernode ZamZamzZam. We could, if we implemented things wrong, do this:
  //
  //        Ref: ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //             Zam............................Zam
  //              ^                              ^
  //             5' anchor                    3' anchor    match start of supernode ZamZamZam with start of reference, and end of supernode with end of reference
  //     this would be a  bug - failing to realise that the entire supernode matches exactly the reference on the first instance of the repeat ZamZamZam
  //
  //    In fact, I found exactly this bug through this test.
  // ******************************************************************************************************************************


  kmer_size = 31;
  number_of_bits = 8;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/one_person_is_alu_Ns_then_same_alu", &bad_reads, hash_table);

  
  //   >7SLRNA#SINE/Alu  
  // GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  // GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  // GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  // CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  // ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  // GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  // AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAGNNNNNNNNNNNNNNNNNNN
  // GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  // GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  // GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  // CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  // ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  // GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  // AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG
    
  
  
  CU_ASSERT(seq_loaded==697);
  
  chrom_fptr = fopen("../data/test/pop_graph/variations/one_person_aluNsalu.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/one_person_aluNsalu.fasta");
      exit(1);
    }

  min_fiveprime_flank_anchor = 10;
  min_threeprime_flank_anchor= 10;
  max_anchor_span = 370;
  length_of_arrays= 740;
  min_covg =1;
  max_covg = 10;
  max_expected_size_of_supernode=370;
  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
						    0, NULL, NULL, NULL, NULL, NULL);

  
  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);



  // ******************************************************************************************************************************
  // 4. Reference = 10kb of chromosome 1. Individual is that same sequence. This should find nothing.
  // ******************************************************************************************************************************

  kmer_size = 31;
  number_of_bits = 13;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/one_person_is_10kb_of_chrom1", &bad_reads, hash_table);

  
  CU_ASSERT(seq_loaded==16320);
  
  chrom_fptr = fopen("../data/test/pop_graph/variations/person_1kb_chrom1.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_1kb_chrom1.fasta");
      exit(1);
    }

  min_fiveprime_flank_anchor = 20;
  min_threeprime_flank_anchor= 10;
  max_anchor_span = 10000;
  length_of_arrays=20000;
  min_covg =1;
  max_covg = 10;
  max_expected_size_of_supernode=10000;
  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
						    0, NULL, NULL, NULL, NULL, NULL);

  
  CU_ASSERT(ret==0);



  hash_table_free(&hash_table);
  fclose(chrom_fptr);


  // ***************************************************************************************************************************************************
  // 5. Reference  = short sequence, individual is idential except for one base change
  // ***************************************************************************************************************************************************



  printf("Start subtest 5 - short sequence for ref; individual differs by one base\n");

  kmer_size = 7;
  number_of_bits = 8;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_short_seq_with_one_base_difference", &bad_reads, hash_table);

  CU_ASSERT(seq_loaded==192);
  
  chrom_fptr = fopen("../data/test/pop_graph/variations/second_person_same_short_seq_one_base_diff.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/second_person_same_short_seq_one_base_diff.fasta");
      exit(1);
    }

  min_fiveprime_flank_anchor = 5;
  min_threeprime_flank_anchor= 5;
  max_anchor_span =40;
  length_of_arrays=80;
  min_covg =1;
  max_covg = 10;
  max_expected_size_of_supernode=40;


  //this time, we expect results. So  malloc some arrays to hold the results. Only use this option for testing, as in general use, we expect to find millions
  // of variants

  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 30);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 30);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 30);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 30);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';
  
  int return_variant_start_coords_array[2];
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  //now do the test!!!
  FILE* fp = fopen("../bin/temp_outputfile_for_testing_subtest5", "w");
  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
  						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);

  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("AATAGACGCCCACACCTGATAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCACA", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("AAGCCACA", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("CTGTACTTGTA", return_flank3p_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==23);



  //cleanup
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);






  // ***************************************************************************************************************************************************
  // 6. Reference = Alu-NNNNN- same Alu, and person is identical to reference except for a single base difference. This should find the single variant!
  //     It is also an example where we traverse the supernode in the reverse direction to that in which is appears in our array. This is 
  //      therefore a good unit test to ensure that our code works in this case
  //     Final comment - note that the coverage is double on a lot of these nodes because we have 2 copies of the Alu. But that doesn't in any way add credence
  //      to the variant call - it's paralogous sequence that doubles the coverage.
  // ***************************************************************************************************************************************************

  printf("Start subtest 6. Reference has an exact repeat of an ALU,and individual differs by one base in one of those copies\n");

  kmer_size = 31;
  number_of_bits = 8;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_both_alu_Ns_alu_with_one_base_difference", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/one_person_aluNsalu_PLUS_SINGLE_BASE_CHANGE.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/one_person_aluNsalu_PLUS_SINGLE_BASE_CHANGE.fasta");
      exit(1);
    }

  min_fiveprime_flank_anchor = 10;
  min_threeprime_flank_anchor= 10;
  max_anchor_span =350;
  length_of_arrays=700;
  min_covg =1;
  max_covg = 10;
  max_expected_size_of_supernode=350;


  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 400);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 400);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 400);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 400);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  fp = fopen("../bin/temp_outputfile_for_testing_subtest6", "w");

  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array,return_variant_start_coords_array_ptr );
  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("GCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("CTGGGAGGATCGCTTGAGTCCAGGAGTTCTGG", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGG", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==59);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);






  


  // ***************************************************************************************************************************************************
  // Reference is as in example 1, just a short sequence. Individual is missing 2 bases in the middle--> we should find this
  //     Note this test is also interesting because the reference and individual both contain a closed loop in their de Bruijn graph
  // ***************************************************************************************************************************************************

  printf("Start subtest 7 - deletion of 2 bases from a short sequence\n");

  kmer_size = 7;
  number_of_bits = 8;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_2_bases_missing", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/person_without_2_bases_missing.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_without_2_bases_missing.fasta");
      exit(1);
    }


  min_fiveprime_flank_anchor = 2;
  min_threeprime_flank_anchor= 2;
  max_anchor_span =40;
  length_of_arrays=80;
  min_covg =1;
  max_covg = 1000;
  max_expected_size_of_supernode=40;


  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 80);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 80);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  fp = fopen("../bin/temp_outputfile_for_testing_subtest7", "w");


  //individual has 2 missing bases, reference does not
  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr );
  
  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("CCACACCTGA", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("ACAC", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGACCCC", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("GACCCC", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==20);


  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);








  // ***************************************************************************************************************************************************
  // Reference is as in example 1, just a short sequence. Individual has an extra 2 bases in the middle--> we should find this
  // ***************************************************************************************************************************************************

  printf("Start subtest 8 - insertion of 2 bases in a short sequence\n");

  kmer_size = 7;
  number_of_bits = 8;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;


  //we use the same two people as last time, but swap their roles
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_2_bases_missing", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_2_bases_missing.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_with_2_bases_missing.fasta");
      exit(1);
    }


  min_fiveprime_flank_anchor = 2;
  min_threeprime_flank_anchor= 2;
  max_anchor_span =40;
  length_of_arrays=80;
  min_covg =1;
  max_covg = 1000;
  max_expected_size_of_supernode=40;


  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 80);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 80);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  fp = fopen("../bin/temp_outputfile_for_testing_subtest8", "w");



  //cf previous test - we are using individua with index 1, not 0 as last time. ie we have swapped which is individual and which is reference
  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr );

  


  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("CCACACCTGA", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("ACAC", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGACCCC", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("GACCCC", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==20);


  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);






  // ***************************************************************************************************************************************************
  // Reference is two copies of a single sequence, tandem repeat of about 36 bases. Individual is the same, with an Alu inserted between
  // Since the supernode in the individual has one copy of th repeat, and then the Alu, and then stops, you should be unable to find
  // an anchor at the Alu-end of the supernode. So this should find nothing.
  // ***************************************************************************************************************************************************


  printf("Start subtest 9 - should find nothing. Insertion between two tandem repeats (each of which is a full supernode)\n");

  kmer_size = 31;
  number_of_bits = 15;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;


  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_alu_insertion", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/person_without_alu.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_without_alu.fasta");
      exit(1);
    }


  min_fiveprime_flank_anchor = 3;
  min_threeprime_flank_anchor= 3;
  max_anchor_span =400;
  length_of_arrays=800;
  min_covg =1;
  max_covg = 100000;
  max_expected_size_of_supernode=400;



  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
						    individual_edge_array,0,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
						    0, NULL, NULL, NULL, NULL, NULL);
  
  CU_ASSERT(ret==0);


  hash_table_free(&hash_table);
  fclose(chrom_fptr);



  // ***************************************************************************************************************************************************
  // 10. Reference is a single sequence which is a single supernode. Individual is the same, with an Alu inserted in the middle
  // Should be able to find anchors at both ends this time, and find the insertion
  // Note that by missing the alu in the reference, it has a load of nodes that won't be seen in the individual, so for example here
  // we will ask for a min 3' flank of 8. This actually means 31+8 bases of agreement.
  // ***************************************************************************************************************************************************


  printf("Start subtest 10 - single Alu insertion\n");

  kmer_size = 31;
  number_of_bits = 15;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;



  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_alu_inserted_mid_supernode", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_one_supernode_and_without_alu.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_with_one_supernode_and_without_alu.fasta");
      exit(1);
    }


  min_fiveprime_flank_anchor = 3;
  min_threeprime_flank_anchor=  8;
  max_anchor_span =500;
  length_of_arrays=1000;
  min_covg =1;
  max_covg = 100000;
  max_expected_size_of_supernode=500;


  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  
  fp = fopen("../bin/temp_outputfile_for_testing_subtest10", "w");
  
  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);
  fclose(fp);




  CU_ASSERT(ret==1);

  //first base of Alu insert = first base of displaced sequence. ie we had xxxxGyyyy and then we insert Alu: xxxx[G-Alu]Gyyyy. So first base of Alu look slike is in flank5p
  CU_ASSERT_STRING_EQUAL("TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCAGG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("AAACTCCGTCTGTACTAATAGTACAAAAAT", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("CCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGAAACTCCGTCTGTACTAATAGTACAAAAAT", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==40); //note the insertion happens at coordinate 39, but the first inserted base is the same as what would be there anyway - a G. 
  
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);





  // ***************************************************************************************************************************************************
  // 11. Reverse the roles of previous test. Reference has an Alu, individual does not. 
  // Should be able to find anchors at both ends this time, and find the deletion
  // ***************************************************************************************************************************************************


  printf("Start subtest 11 - deletion of an Alu\n");

  kmer_size = 31;
  number_of_bits = 15;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;


  //we use the same two people as last time, but swap their roles
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_alu_inserted_mid_supernode", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_alu_in_middle_of_supernode.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_with_alu_in_middle_of_supernode.fasta");
      exit(1);
    }



  min_fiveprime_flank_anchor = 3;
  min_threeprime_flank_anchor= 8;
  max_anchor_span =500;
  length_of_arrays=1000;
  min_covg =1;
  max_covg = 100000;
  max_expected_size_of_supernode=500;


  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  
  fp = fopen("../bin/temp_outputfile_for_testing_subtest11", "w");


  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
						    individual_edge_array,0,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
                                                    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);
  fclose(fp);


  CU_ASSERT(ret==1);

  //first base of Alu insert = first base of displaced sequence. ie we had xxxxGyyyy and then we insert Alu: xxxx[G-Alu]Gyyyy. So first base of Alu look slike is in flank5p
  CU_ASSERT_STRING_EQUAL("TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCAGG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("AAACTCCGTCTGTACTAATAGTACAAAAAT", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("CCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGAAACTCCGTCTGTACTAATAGTACAAAAAT", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==40); //note the insertion happens at coordinate 39, but the first inserted base is the same as what would be there anyway - a G. 


  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);





  // ***************************************************************************************************************************************************
  // 12. Reference is an Alu with a different Alu inserted in the middle. Individual is just the first Alu - ie they have a deletion of an Alu from within an Alu.
  // 
  // ***************************************************************************************************************************************************


  printf("Start subtest 12 - deletion  of one Alu that had been inserted within another\n");

  kmer_size = 31;
  number_of_bits = 15;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;


  //we use the same two people as last time, but swap their roles
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_is_alu_other_has_2nd_alu_inserted", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_alu_in_middle_of_alu.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_with_alu_in_middle_of_alu.fasta");
      exit(1);
    }


  
 //    reference is
 //   > AluJo#SINE/Alu inserted in middle of 7SLRNA#SINE/Alu
 //   GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
 //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
 //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
 //   GCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGA
 //   GGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACAT
 //   AGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCGGGCGTGGTGG
 //   CGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCT
 //   TGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACT
 //   CCAGCCTGGGCGACAGAGCGAGACCCTGTCTCAAAAAAAAAAAAAAAAAA
 //   AAAAAAAAAAAA
 //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
 //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
 //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
 //   AAAAAAAAAAAAAAAAAAAA
 //   ..plus N's on the end
    
    
 //   and individual is
    
 //   >7SLRNA#SINE/Alu
 //   GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
 //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
 //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
 //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
 //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
 //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
 //   AAAAAAAAAAAAAAAAAAAA
    
    


  min_fiveprime_flank_anchor = 3;
  min_threeprime_flank_anchor= 3;
  max_anchor_span =500;
  length_of_arrays=1000;
  min_covg =1;
  max_covg = 100000;
  max_expected_size_of_supernode=500;


  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  fp = fopen("../bin/temp_outputfile_for_testing_subtest12", "w");



  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
						    individual_edge_array,0,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
                                                    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("GAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAA", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("CGGGGGACCACCAGGTTGCCTAAGGAGGGGT", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("GCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACATAGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCTTGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGGGGGACCACCAGGTTGCCTAAGGAGGGGT", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==151); 


  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);




  // ***************************************************************************************************************************************************
  // 13. Reference is 10 kb of chromosome 1 plus 1kb of sequence inserted mid-supernode, and person is identical, except for that 1kb is missing 
  // ***************************************************************************************************************************************************

  printf("Start subtest 13: deletion of 1kb within 10kb of sequence\n");

  kmer_size = 31;
  number_of_bits = 15;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_1kb_deletion", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/first_person_10kb_chrom1_plus_1kb_inserted_mid_supernode.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/first_person_10kb_chrom1_plus_1kb_inserted_mid_supernode.fasta");
      exit(1);
    }

  // Let's be clear about what this test looks like. 




  //       The individual is 10kb of chromosome 1, which has the following supernodes

//>node_0  cylcle:false length:4546 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:1 lstf:2
//TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACTCCGTCTGTACTAATAGTACAAAAATTAGCCAGGCGTGATGGTGTGCACCTGTAGTCCTTGCTACTCAGAAAGCTGAGGCAGGAGAATCGCTTGTACCCAGGAGGCAGAGGTTGCAGTGAGCAGAGATTGTGCCACTGCACTCCATCCTGGGTGACAGAGTGCTATGAGTCACCACACCTGGTATGAGCCACCGTGCCTGGCCCACAATGACTTTTACACATGTTGTTAAATCATCTTACAGATTTTATAATTTGGGGGAAGAAAAGTTTTACTAAATTGTCTTTTAATGGAAACTCTACAAGAACCAGAATCTTTGCTTTGTTCACTTATGTATCCATTCCTAGGCCTAGAAAAATGTCTGACACATAGCGGCAATTATTCATTGAATAAATGGACCCAGGGATAGTACATTAGCTATGCTATATGCATACATTAAAGATGTAGATTATCGACTTTCAAAAGATAATTAATGTAACTTCTTACTGCTTCTGAACATGTTTGTGAGTTATATTGCTGAGGGACCTTTATCTTCTCATTCTTTCATCTTAACCCAGTGTTATAAAATTGAAATCACCAATATTATTCCATATCTAAAATTAATATCTACCTTGTAAAAAATATCACTCTGCTGCATTTGAGAATAGACTTTTTAGGTAATAATGATGCAATCCATAGGGTTTTTTGGGGGCACAGAGGGATTCATGCTAACAGAACATTTTATTTTCTATTTTCCCAGAGCTGTAAAACATGAAATTACGGTAGTATAAGGCATATTTTTACTCTTTTTATAATTTTTTCTAAAAAAAATTAGTGTTTGTTCCCTATATAACTTTTAACTTTATAGGTAAATATTTGTTTCTTTCAGCTCCAGTTTTATGTGAAATAGAGTTTTCAGATTTATGTAGCATGGAAAGTTTTAATACGTCAGAGTTACTGATTTTTGCCAATCATTTTCTCAATTATTTCTTTTTTATCTTTAGTTGATTTTTTTGTAGTGACACATTTTGTTTCTAGTCTCATTTCCTTTTGTTTATATTCTATGTATATTTCATTTTTGGTTACTATGAGAATTACATATAACATCCTAGAGTTATAACATTTTAATTTGAATTTATTTCAACTTAAGTTCAATCACATACCAAAATTCTACTGCTATATATATAGCTCTACTCTTTTTATGTTATTGATGTGACAAATTATATCTTTATTCATTGTATACCAGCTAACAGATTTACAATTACATTTTATGCATTTGCCTTTTAAATTATGTAGAAAATAAAAAGCAGAGTTACAAACCAAAATTACAATAGGACTGTTTTTATGTTTGTTTATGTATTTACCTTTACCAGAGAGCTTTGTATATTCATACAGCTTGCTTATTTACTTATATAGTTATTGCCTAGAGTTCATTTATTTCAACCTGAAGGACTTAACACTTCCTGAATGTCAAATTCAGGGATAAATGGATTTTTTTCAGTTTTAAAAAAAAATCCGGAAATGTCTTAATTTCTCCTTCATTTTTGAAGGATAAGTTTTCCAGCTATATATTTCTCAATTGACAGGTTTCTTCATTATTTTAAATATATAATCCACTGCCTACTGGCCTTCAAGGTTTCTGCTGAGAAATCAGCTGCTAATGTTATCTGGATCCCTATCTGTGAGAGTTGCTCTTCTCTCTGAGTTTTCAACATTCTCCCATTATCTTTTTTTTGTTTGTTTTTGAGACAAATAATTGTACATATTCATGGGATACAGAGTGATATTTTGATACATGTATACAATGCCCAGTGATCAAATAAGGATAATTAGCGTATCCATCACCTCAAATAGTTGTCATTTATTTGTATTGTGAACAGTCAACATTCTTTCTTCTAGTTTTTTAAATTTATAAACATTTAAATTTTATTACAGAAATTTAAATTTTTTGATTCTGAAAAAGTCATATATGTATGCAACATTTTTTATCGTTTATTTATATATTTATGCATCTTTCCTTTTAGTTTTGACAGAGATTTTCTATTTTATCATTATTTCAAAAGAACTCTTACCTGTATTTATTTATCAAGTATATTTCCCTTGTTTTTTCCTAGTATATTAATTTATTTACTTATCTTCTAAAAATCCTCCATATAATCTGTTTATTTTGTTTCCTTTCTATAATTTCTTCAATAATTAGTTCTGTTCTATTTTCCATTAAAATATTTAAATCTTGTATGAATTTTTGTCAGATTAGAAATTTAGGGCGTTTCTTAATTTCTCTATACTCTAGCTTTTGACTTTTTTTTTCTGACCTAAGAGGTATTTAGAGCACATTTTAGATTTTTTATTTTGACTAATCATTTAAAATGTATACTAATCTTCAATTTAAATAAAAAACTGGTCTATAGTGACAAAAATTACAAATGAGCCTAACTAATAAATTATCAGCTGTGTTTATATGTATAAGCATGCACAGATTTTGGTAAATATGTACATAGTATATTGGTGAGCTTATTTTTATCATTCTTAACTCATTGTGTAGTCTAAACGTTGGGGAAAAAATAACATACAATAATCAGATGGTGTGAATAAGAAAATTGTTCTAATGTTTGTAAACCAAGCAACTGTTTTAACTGCTCCCCTCTTCCTGATTGACTTCTAAAAGGGATTGATCCATATTGGGTCCTATCATGTACGTCACGGTATAACATCTCCAGCTATAAAATGGAAATTTGAGAATAACTTTGCTGCTACTCAGATATATTTTATTTCAAAAACATACACTAAGGTGTTGCTGTTGGATCTTTCCAAAAACATATTCACACAGAACTTTCAATCACACTGAGCCATATTTGAACAATCTTTCAAGGTCAGCTCTGGCATAAGCTAACATTATACCATTTAACTCAGAAATTTCTTTAGTATTTGATTAATGGGTTTATGTTTGATATGTAATGTAATTTTCTAATACTAAATCAAGTGGTAATTTTGTTAGTCAAGTTGATTTAGTGGCTTGGGAAGAAAGCTTTTAATGTTCCCCTAATTTTTCTTTCCTTTGACATGATCCTTCACATGTCTTATTTTGCTTAGTGATTTTTCTTTTTTTTTTTTTTTTTTTGAGACAGGGTCTTACTCTACCACCCAGGCTTGAGTGCAGTGGTGCAATCACAGCTCATTGCAGCCTTGACCTCCCAGACTCAAGCTATTCTTCCACCTCAGCCTCCCAAGTAGCTGGTACTACAGGCACATGCCACCAAACTTGGCTAATTTTTGTATTTTTTGTAGAGACAGAGTTTTGCCAAATTCTCAGGCTGGTCTGGAATTTCTGGGCTCAAGTAATCCTGCCTTGGCCTCCCAACATGCTGATATTACAGACATAAGCCACAGTACCTGGCCAGTTTTCTTTTTAAAAAAATCTATTGGTTATTAATTTGAAGCCTTCCTTTTCATAGCTGTGCTCCTTAATTGGGAGCAAACATGAATGGACCACAACTTAGCCAATTTTCTATATACGATCTTTGCCATCCTAATTTAAAGGAATATTAATTCTTTCTTTTCCTCTTTCATTCCACAAACCTGTATTGACTACATCTAAGTTCTAAATGGTGCACTGGATGTTGAAAAAGTTGATGATGAGCAAGAACAAAATTCCTGCTTTCAGGAGACTTACAGTTCAATATGGGAAATATAATTTGTTAAAATATAAAAGTGCAATTGTGTTACATGCTGTACGAAGTACATGTTGACATGTGAGCATATAATAAATGGGCTGGAGGCCAGAGGATTGCCAAAGAGAATGGGCCTCCTGCTGAGATGAAAAGTTGAGCAGGGATTAGTTGGCAAAAGTGGAGGGACGATCCTTTCTAGGCAGGAGGAAGAACATGTACAGAATCTCTGAGGTGTGATGCGACAAAGTCTATATAAAAAACTGAAGAAAGGTCTAATATGGCTTAAATACAGAAGCTAGTAGGAGAGGAGTCGAAAAGAGGCTGGAGAAGTAGAAAGTGTCTGCATTCTGCAGGAACTTATATTGTATAAAATATATTGTATATTGTATAAAAAGAATTTCTCTTTATTCTAAGTGCAATGTGAAGCCAATGAAGTGCTTTAAACAGGTGATGTGATTTGATTGAATTTATTACTTCACTTAACAAATATTCATTACATGCCCACTGTTTGTCAGATATTGCTCTAGCCCCTGGTGATACAGTAGGGAATAAAACAGGCAAAAATCCCTGTCCTCTTGCAGCTTATAATGGACTGCAATGTTTAATATGTCAGAGGAGGTCCACGGAGGAGTGACTTCTAAGCAAGAATCTGAAAAAAATGAGGATATCTAAGGAGGGAACAAATGGTTCAAAAGCCCTATAATTGCAAGCAGGCATGATGAAGCAATTGCAGTTGTCCTGACTCTCAACACCGTGGAACTCAAAGGAGATGGAAAGATTCCTTCTCTCCCTCATATATTTTCTCTCTTTCTGTCTATATATATAGAATATGAGACATTTCCCTAATCATTATGT
//>node_1  cylcle:false length:1696 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//ACATAATGATTAGGGAAATGTCTCATATTCTTCTACTCAGAAATAAGCAATATAGCAATTACTGTTTTTTACATTTTACAGTTACAGTTTCAGAGAAAGTTTGATATTTATCTAAAATTTTTCAATGTATGAACTTTTTCATTTGACAAACCATAATTGTACATATTCTTGGGATACAGAGTGATATTTTCTTTACATGTATAGAATGTGTAGTGATCAAATCAGGGTAATTTCCACTAATTTAAAATGCCACCTTTATGTTATTGTAATTTATATATATACTATATATATACACACACACATATATATATATACATGTCCACATACAGTGTGTGTGTGCACATGTACACACATGCATATGTGTATATAATGCCCAGTATAAGCAATGTGCACAAATAAAATTAGCTAACAGAGATAGTATAGAGTGAGAGGAGAGGCAGATTAATCTTTGAGGAAAAGCACAATTTTATAGCTGAATGGAGAAAGCTGAGGTGGTTTCTAAGATGGAGAATAAGACGAAAAATGTAAGTACGTTGTCTGACTGAATTCAAGAAAGAAGGGTAAAAGAGAAGAAAGTAGTGGTCTTATCATTAAATGCCACAGAGAGGTAAAGATAAAGACAACATATTGTTTTGGGTTTAGTAATTTAAGGGTTACCAAATTCCGTTTTGGAGGAGGAACAGATTCCATGTCCACTAGAATGGAATGAACAAGAAATGGAGGAGGAAAATAGGTAGTTTTTCAAAAGTTTTCAAAAATATGAAAAGAAGAAATGAAGTGGTACTTGGAAGAGATTGTTGAAATGGGAGAGACTATGGTGGCTTGTTTAGAAGCAGTTGAGATAGATCCAATTGAGATAGAGATATTGACTATATAAACAAAAGAATGACAAATTAATAGTGTAATGGATAACTTGACTTTGGCAAATATTGTGAATTTTTGTGAAAGTACAACTAAAAGGCAATGTCACTCCAATAATCACCAGAGTAATCAATTTGCTTATTGCTGTCCCTTTAAATATAGTTCTCTGGTATCAACTAACATGTTTTTAACTAATGATGCTTCTTAAAGAAAAGGGAAAAGACCTTTTTCTTTCTTTCAGTCTTCAATGATTCACTGCTTCATCTCGCTCCACCAAAGATAAATGAAATCTACATCTCTTATACATTAACAATGCATGACAATTTACAAATAGCTAAATTTTTGGAGCTAACTTTAAGTACCTGAATGGAATTTAATCAACCCACTAATCTCCTTCTCACTTCTCAGTTATTTATCAAGTTTATGTCAAGGGACAAGGAAAAATTATCCAAACATTGTTTAAAACAATCATCATTAATTAGTAACACTTATCCAGGGGGGTTTTTAACCTTTCCCCCACTCAAGGATTATTCTAATGTCAGAGTAGAATAAAAAATAAGTGCAGTGATGCTGACTCTTCCAAGCTTAACATTTCTCACAAGTCAATTAGCTTTGTACTGGGAGGAGGGCGTGAAGGGCTGCTTGCGGTAGTTGTGTAGCAGCAGCACAATGGCCGCAGACAAGGAAAACAGTTTCTAGGAATTCCTCGTATATAATTTTATATTTTTGACAAGATTAATGACCCATGCTCCCTTCCTCTCCATTTCTTTTTTTGGAATTCTGTTGGTATGTAGTTACTATATTTTATTAAAGGAAATTAGCCTTATCTCTTATT
//>node_2  cylcle:false length:566 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//AGGAAATTAGCCTTATCTCTTATTATATTTTTTATGACCTTCAAAGTAGTGTCTCTGCTTAAAAGTGTACCCTGGCCGGGCGTGGTGGCTCACACCTGTAATTCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTGTACTAAAAATACAAAAAATTAGCAGGGCATAGTGGCGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTCAGGCAGGAGAATGGCGTGAACCCGGGAGACGGAGCTTGCGGTGAGCTGAGATCGCACCGCTGCACTCCAGCCTGGGCGACAGAGCAAGACTCCGTCTCAAAAAAAAAAAAAAAGTATACCCTGAGGCACACATCAAGCGACATGTAGAGTTCATAAATTCTGGCCAAATGGTCATACCTCAAACCTCATCAGCAGTAAGGCTCTTTACTTGCACTGACAAATATGAACGCTGGGGAATTTGGAAATGATATATAATATATAATATTATATATATAATAGATATATAATATATAATATATATAATACATAT
//>node_3  cylcle:false length:2989 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:1 fst:0 lst_coverage:2 lstr:1 lstf:2
//AATTTGTTACAGCAGTAATAGAACAAGTGGTTATCCATATGAGGCAAATTAGATTGGATACCTATCTCCAATAGAAATCAATTCAAGGTGAATTCCAGGAAAATACTTAAAACATTTAGATTAAAAATAAATGAGAATTTTTGTTACTTTTGGTAGGTCATAGAACCAAGAAAAACAAACATTAAGGAGGAAAAATGAACATATGACTACATCAAAATATAAAGCTTCTCTATTTGGAAGATATCATAAGGTGACAAATCATAAACTGTAATATTTACAACATATATATAAGTGAATAAATATACATTTAGAATATATATGAACTCCCAAAAATCAACAGGAAAAATAAGACATAGAACAAGCAAAATGCATAAACAAAAGAAGGCAAAACAAAAATAATGACTCATAATTATATGAAAAGAAGCTCATCTTCATAGATGAGCAGATAAATGCAAATTAAAACCACCCTGAGATGCTTTTTACATCCATGAGCCTGATAAAAGTTAGAGTCTAAAAGTAATAACAAAGATGGGAAGTAATAGAAAATCTTGTCCATTACTGGTTAAAGTATAAACTGATACAGCTACTTTATAGAATATTACATTATAGAATAAAGTTGTGAGTATGTATATGCAGTGACTCAGCATATTCATTGCTAGTATGTACTCAAGAGAAACTTACAGGAGTGGACTAGGAAGTAAATACAAAATGATTACAACATTGTTTGTTATATCAAAAAATAAAAAAGACACCCAATTTTCCAGCAAAAAAAATAAGTAAAAATAAATCCTGGTGTATTCTAACAATGGAATAATATATAGCCATTAAAATAAATCAACTATTACTGTACATATGAATGTAAATATCAGCAAAACATATTGTTTAGTGAAAAACTAAGAAGCTGAAGAAGAATATATACAATATGGTTACATTTATATGAAGTCCAAAAACTTGCAAAATAAAGAAATGTATTTAGAAATAGATTCACATGTGAGAAAACTAGAAGAAAATTAATGAAAGGATAAGAGGGATAGCAGTAATTCTGAGTAGTTGAGGGAATTTCAATTGGAAAAAAATAATATCATATTCTTTAAGTCAGGTAGTGGGTATTAGCATTTGTTTTACCATCGTTCTTTATTCTTATAGCTACACTATATATTTTCAATGTATTTAATGTATTTTTTGCATAATTAAATATTATGCAATAAAAATGAGAAAACAAAAAAGTAGAAAATGATAAATTACAATAAAGAAATGGAGAAAAAATTATAATCTAGTTGAGTAATGGTATATTACATAGCTATTTTCTTAAGTAGATGTATGTACATGATGTATGCACGATTGTACATACATGTTCTTAATTATATATAAATATATATGTACATATTTTTAATATAAAATACTAAACAAAGTACACCAAAATATTAGCTCCTATGTTAGTGAGATAATGTTTTGTTTTTTTGTATTTTAAGTTTTACATAGTAGGTGTATTTTTCTGTTTTCATACTGCTATAAAGAACTGCCCAAGACTGGGTAATTTATAAAGGAAAGAAGTTTAATTGGCTCACCGTTCAGCACAGCTTGGGAGGCCTCAGGAAATCTACAATCATGGCGGAAGACAAAGAGGAAGCAAGCCAGCTTCTTCGCAAGGCAGCATGAAGAAGTGCCGAGCAAAGGGGAAAGAATCCCTTATAAAACCATCAAATCTCGTGAGAACTCACTATCACAAGAACAGCACAGGGGAAACTGCCCCCATGATTCAATTACCTCCACCTGGTCTCTCCCTTGACCTGTGGGGATTATGGGGACTATGGGGATTACAATTCAAGATGAGATTCAGGTGGGGATACAAAGCCTAACCATATCAGTAGGCATGTATTGAATTTTAAACTCAGAGAAAAATACTAGTGTTTTTATAGGATTCTTACTAAAGAAAAACCAGAAAGTAATAAACCATCTACGCTAAGACATAAAATTCAGTTGTTTAGTTACAAGATAGAATGTGGCCTTGTAAGAAAGCAAATTAACTTCTAACATACAAAGCCTTAGAGAAGATTCAAGTGACTGACGGATCTTAAACAGAGCTATTATTACAACTCGAACTGCAGTAAAATATCCTCAGCAACATAGATGTGTGTGTTTCACTAGTCAGAGCAATACAAATTTAATGAAACTCCACTGGTGGTGTTTTTAATCAGACAATTTCTGAAGATGTCCTGGCTTATTCACAGATGCAAGCCAAATCTCTAGAAGAGTACCATAATAAGAAAAAAAAGAATACAGGCAATTGAGAGCTGTTCCAAAGTTTAGGGAGTTTTTGTAAGGAATTAATAAATAAAAATGTTCTTGAAAGACAGAAATTAATATGCAGTTCATACTGTCAGAATTGCAGGCAATTTATCAAAGTCCCCTAATCCTCCAAAATCGCTATTTTTTTTTTGACACACACTTTACAGTACAGAAGAAAATGTCTCCGGCAATAAATCACAAAGTTAAAATTACCTAGTCTACAATTAACTACACAGTGATGGTAAATCATTTTCTACCAAAAGAAAGAAATGTCTTGTCTATTCAGGTTCTGCTCTACTTAAAAGTTTTCCTTGTTGGCGAGCAAGTGGTTAGAAAATTATATTTTATACGTACATTCAGCTTAACTATCATTCAGCTCAGGAAGATGACTCAGGGCCTTATCCATACCTTCAAGTTTGCTCTTAGCAAGTAATTGTTTCAGTATCTATATCAAAAATGGCTTAAGCCTGCAACATGTTTCTGAATGATTAACAAGGTGATAGTCAGTTCTTCATTGAATCCTGGATGCTTTATTTTTCTTAATAAGAGGAATTCATATGGATCAGCTAGAAAAAAATTAAGAGGAAAATCACATGGAAAGTTATATATTATATATCTATTATATATAATATATATCTATTACATATTATATATTGTATATCTATTACATATATATTATATATGTATTATATATATTATA
//>node_4  cylcle:true length:1691 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:2
//TGGCCAGGCTGGTCTCGAACTCCTGACCTCAAATGATCCATCTACCTCGGCCTCCCAAAGTGCTGGAATTACAGATGTGAGCCACAATGCCCGGCCTTATTTTCTACAACTTTGGTAACTTTAGCATATACCCCAAATCTGTAAGACATAATATTATAATTCAAATGCAACTCATGGCTTCTCTTTGTACTCTTTCTCTAGCTTTTGAATTATTTATTCTAATACCAGTTTTAATTCTGACACAAAATCATGGGAGTTCTAATCAAAATCCAACCTTTTATCATAAAAACTATGAAGAAATTATGAGTAGAATTTAAAAAGGAAAATAGGCCTATTAATTAGATTTGTCTTTGTAGCATTTAACTCTATAATAAATAATATTTTATGCCTATGAGTCCCCAACAAAGCCTCCAGCTTCTATTTAGATATAAACTGTAAAAGTCACTACTGGATCCACAAGCAAGACTATGGTAAATAAATTTCTCCACCTAACCAGCTTCTTTTACATGATGTTACATGTTTCTTTTGTTTTTTCATTTTGGCAAATATTGATTGTCATCTTCGTGTTTGTCTATGTCCTAAGTGCTGGGATACAGAATCTGAAAAGATGGACACAGGACCTGCCTTCAAGTTCACCCCCTTTTTTTTTTTTTTTTGAGATGCAGTTTTGCTCTTGTCACCCAGGCTGGAGTGTACTGGTGAGATCTCTGCTCACTGCAACCTCCACCTTCAGGGTTCAAGTGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGTGATTACAGGTCCCAGCCACCACGCCTAGCTAATTTTTGTATTTTTAGTAGAGACAGCGTTTCATCATGCTGGTCAGGCTGGTCTCGAACTCCTAACCTCAGGTAGTCGACCCACCTCGGCCTCCCACAGTGCTGAGATTACAGGCATGAGCCACCACGCCCTGCTAGGAGTTCACGCTTTAGTTGGGGAAAATATACAATAAGCAAGCCAGTTTTTAAAATGAGAACTGCAATTAGAGTTAAATGCTACAAAGACAAACTCACAGGAAGATGGGATGTAGAATGATAAGGCTCTCAGAATAGTAAGAGAAACTATTGCTTCTTACGATGTTTGTCTTTCTTTGTATCGGTGCTCAGCTGAGTCTGCAGTGCTTCAGAGGCAGCTTTCATTTTATAAAAATCTATGATTTCTCCTTCCAGTTTTTTTTTCTCTTCCTCGAGCTTCCTTATCTCCTCCTGTTGAATCATTTTAAGATGCTCGAACTTGTCCTGCAGCTGTGAAACCAATGTGCAGTTGTGACACCAAAGCAGTGTGGCTGAACACCTAAAAGAATACGCTTTTTTTCTGATTATCAAACAAACCCAAATCATCACAGTAGACCACGATCTTAATAACAATCTCAAAAACTCAGGAGTAAACACTCAGATATGGAATTTTTCTTTTCTTTCTTTTTTCCTTTTATAAGATGGAGTCTCACTCTGTTGCCCAGGCTGGAGTGCACTGGTGCGATCTCAGCTCACTGCAACCTCCATCTCCCAGTTCAAGTGATTCTCCTGCCTCAGCCTCTTGAGTAGCTGGGACTATAGGCATGCACCACCACTACAGGCGTGTGCCACCACACCTGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTTGCCATGATGGCCAGGCTGGTCTCGAACTCCTGACCTCA
//>node_5  cylcle:false length:462 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:1 lstr:0 lstf:1
//TGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCCTCCCGCTTTGGCCTCCCAAAGACTTTTTTTTTTTTTTTAATATAGAGACAAGTTCTCAGTACGTTGCCCAGGCTGGTCTCAAACTCCTGAGCTCAAGTGATCCTCCCACCTCAGCTTCCCAAAGTGCTGGGACTGACTGGATGCAGTGGCTCATGCTTGTAAACTCAGCACTTTGGGAGGCCAAGGTGGGAGGATCGCTTGAGCCCAGGAGTTCAAGACCAGACTGGGTGATATAACACAATAGTCAACTTCAACAGGAGAGAGAATCTGTAAACTTGAATATAGATCTTCCGAAATTATCCAGTCAGTGGACAGAGAAAAAAAGAATAAAAGAGAGAAAAGAAGGCTGGGTGTGGTGGCTCAAGCCTGTAATCCCAACACTTTGGGAGGCCGAGGCAGGCAGATTAAGAGGTCAGGAGTTCA
//>node_6  cylcle:false length:116 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//AATAAGAGATAAGGCTAATTTCCTTTAATAATAATAAAATCCTTTAATAAAAATATAAAGGAATAATATAATAATTTTCTTTAATAAAATATAATAAGAGATAAGGCTAATTTCCT
//>node_7  cylcle:false length:41 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//ATATATAATATATAATATATATAATACATATATAATATATA
//>node_8  cylcle:false length:38 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//AAAATATAATAAGAGATAAGGCTAATTTCCTTTAATAA
//>node_9  cylcle:false length:48 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//TATAATATATATAATACATATATAATATATAATATATATAATACATAT
//>node_10  cylcle:false length:100 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//AGAATATGAGACATTTCCCTAATCATTATGTGTAATTACAATTACATATATATATGTAATTGTAATTACACATAATGATTAGGGAAATGTCTCATATTCT


//The reference is the same as the individual, except it has about 1kb inserted in node 0. To be clearer, the reverse complement of node 0 is what is in the 10kb without insertion:

//GTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATCTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCAATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAATGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA  TTGCACCACTGCACTCAAGCCTGGGTGGTAGAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCTTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTAAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGC

//  and the insert happens just after TTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA, where I have put a space above


//The inserted sequence is

// ACCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACAT
// which I took from chromosome 2.




//NOTE - the inserted sequence starts with an A. The sequence that follows the inserted sequence also starts with an A. So our algorithm will add that A to the flanking region.
// We allow for this is checking results match expectations



  min_fiveprime_flank_anchor = 21;
  min_threeprime_flank_anchor= 21;
  max_anchor_span =8000;
  length_of_arrays=16000;
  min_covg =1;
  max_covg = 100000000;
  max_expected_size_of_supernode=8000;


  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *5600);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  fp = fopen("../bin/temp_outputfile_for_testing_subtest13", "w");





  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table,  fp,
                                                    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("CATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGCAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTATTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA", return_flank5p_array[0]);

  CU_ASSERT_STRING_EQUAL("GAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGTGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAAAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCACTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGGTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCA", 
			 return_flank3p_array[0]);

  CU_ASSERT_STRING_EQUAL("TTGCACCACTGCACTCAAGCCTGGGTGGTA", return_branch2_array[0]);

  CU_ASSERT_STRING_EQUAL("CCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACATATTGCACCACTGCACTCAAGCCTGGGTGGTA", return_trusted_branch_array[0]);

  CU_ASSERT(return_variant_start_coords_array[0]==6722); 



  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);





  // ***************************************************************************************************************************************************
  // 14. Identical to previous test, but each person has identical 600 lines of chromosome 12 added (roughly 30kb) before the sequence that was there in previous tes
  //     This is purely to test that the start coordinate of the variant is found correctly, even when you have to load more and more sequence into the array.
  // ***************************************************************************************************************************************************
  

  printf("Start subtest 14: same deletion of 1kb within 10kb of sequence as in subtest13, but with some extra sequence bunged on front of fastas, so check whether we get start coord of variant right\n");

  kmer_size = 31;
  number_of_bits = 15;
  bucket_size    = 10;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_1kb_deletion_both_with_600lineschrom12_beforehand", &bad_reads, hash_table);

  chrom_fptr = fopen("../data/test/pop_graph/variations/first_person_600lineschrom12_then_10kb_chrom1_plus_1kb_inserted_mid_supernode.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/first_person_600lineschrom12_then_10kb_chrom1_plus_1kb_inserted_mid_supernode.fasta");
      exit(1);
    }


  min_fiveprime_flank_anchor = 21;
  min_threeprime_flank_anchor= 21;
  max_anchor_span =8000;
  length_of_arrays=16000;
  min_covg =1;
  max_covg = 100000000;
  max_expected_size_of_supernode=8000;


  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *5600);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  fp = fopen("../bin/temp_outputfile_for_testing_subtest14", "w");





  ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table,  fp,
                                                    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("CATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGCAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTATTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA", return_flank5p_array[0]);

  CU_ASSERT_STRING_EQUAL("GAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGTGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAAAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCACTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGGTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCA", 
			 return_flank3p_array[0]);

  CU_ASSERT_STRING_EQUAL("TTGCACCACTGCACTCAAGCCTGGGTGGTA", return_branch2_array[0]);

  CU_ASSERT_STRING_EQUAL("CCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACATATTGCACCACTGCACTCAAGCCTGGGTGGTA", return_trusted_branch_array[0]);


  CU_ASSERT(return_variant_start_coords_array[0]==42662); 
  
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);



}



*/









void test_db_graph_make_reference_path_based_sv_calls_null_test_1()
{


  // ******************************************************************************************************************************
  // 1. NULL test. Load short fake reference, and load a person whose sequence is identical. This should find nothing.
  // ******************************************************************************************************************************


  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;


  
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/one_person_with_Ns_on_end", &bad_reads, hash_table);

  //>one read
  //AATAGACGCCCACACCTGATAGACCCCACAC 


  // The de Bruijn graph of this in 7mers is as follows:
  // a line of 8 nodes, with a loop of length 16 off the 9th node. ie all nodes have 1-in and 1-out, except the 9th node, CCCACAC, which has 2 ins and 2 outs.
  // so it's a supernode of length 8, a loop-supernode of length 16,

  // our algorithm should look at the supernode starting with AATAGAC, then the supernode starting at CCCACAC, then stop


  
  CU_ASSERT(seq_loaded==1462);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/first_person_with_one_read_and_Ns_on_end.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ./data/test/pop_graph/first_person_with_one_read_and_Ns_on_end.fasta\n");
      exit(1);
    }

  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= 3;// I want to try to attach anchors, and the supernodes are quite short, so for this test only ask for (ridiculously) small anchors
  int max_anchor_span = 20;
  int length_of_arrays=40;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=20;



  
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
							individual_edge_array,1,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
							0, NULL, NULL, NULL, NULL, NULL);

  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);


}





void test_db_graph_make_reference_path_based_sv_calls_null_test_2()
{

  // ******************************************************************************************************************************
  // 2. Harder NULL test. Reference=an ALU and load a person whose sequence is also that ALU. This should find nothing.
  //     Our implementation of db_graph_make_reference_path_based_sv_calls assumes we have a very large ref fasta, so we load big arrays of bases
  //     which are much longer than we expect any of the supernodes to be. In this case, we know the supernode is as long as the Alu.  
  // ******************************************************************************************************************************

  int kmer_size = 31;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  //just load one person who's sequence is an Alu. Then take reference which is that same sequence and try to find variants
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/just_one_of_the_two_people.txt", &bad_reads, hash_table);

  //   >7SLRNA#SINE/Alu  plus GTTCAGAG at start and GTCAGCGTAG at end
  //   GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  //   AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG
  // plus loads of N's on end
   
  
  CU_ASSERT(seq_loaded==5158);
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/person1.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/person1.fasta");
      exit(1);
    }

  int min_fiveprime_flank_anchor = 10;
  int min_threeprime_flank_anchor= 10;
  int max_anchor_span = 370;
  int length_of_arrays=740;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=370;
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
						    0, NULL, NULL, NULL, NULL, NULL);
  
  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);

}




void test_db_graph_make_reference_path_based_sv_calls_null_test_3()
{


  // ******************************************************************************************************************************
  // 3. Reference = Alu-NNNNN- same Alu, and person is identical to reference. This should find nothing
  //    Note this really is testing something different. Since the reference is basically the same sequence repeated twice (with N's in between),
  //    there is a risk that we would match the 5' anchor of the (one and only) supernode at the start of the reference, and the 3' anchor at
  //    the end of the 2nd copy of the repeat.
  //    To be clear, look at this:
  //        Ref:  ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //        Indiv ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //      So we look at supernode ZamZamzZam. We could, if we implemented things wrong, do this:
  //
  //        Ref: ZamZamZamNNNNNNNNNNNNNNNNZamZamZam
  //             Zam............................Zam
  //              ^                              ^
  //             5' anchor                    3' anchor    match start of supernode ZamZamZam with start of reference, and end of supernode with end of reference
  //     this would be a  bug - failing to realise that the entire supernode matches exactly the reference on the first instance of the repeat ZamZamZam
  //
  //    In fact, I found exactly this bug through this test.
  // ******************************************************************************************************************************


  int kmer_size = 31;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/one_person_is_alu_Ns_then_same_alu", &bad_reads, hash_table);

  
  //   >7SLRNA#SINE/Alu  
  // GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  // GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  // GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  // CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  // ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  // GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  // AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAGNNNNNNNNNNNNNNNNNNN
  // GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
  // GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
  // GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
  // CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
  // ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
  // GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
  // AAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG
    
  
  
  CU_ASSERT(seq_loaded==697);
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/one_person_aluNsalu.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/one_person_aluNsalu.fasta");
      exit(1);
    }

  int min_fiveprime_flank_anchor = 10;
  int min_threeprime_flank_anchor= 10;
  int max_anchor_span = 370;
  int length_of_arrays= 740;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=370;
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
						    0, NULL, NULL, NULL, NULL, NULL);

  
  CU_ASSERT(ret==0);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);


}




void test_db_graph_make_reference_path_based_sv_calls_null_test_4()
{

  // ******************************************************************************************************************************
  // 4. Reference = 10kb of chromosome 1. Individual is that same sequence. This should find nothing.
  // ******************************************************************************************************************************

  int kmer_size = 31;
  int number_of_bits = 13;
  int bucket_size    = 10;
  long long  bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/one_person_is_10kb_of_chrom1", &bad_reads, hash_table);

  
  CU_ASSERT(seq_loaded==16320);
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_1kb_chrom1.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_1kb_chrom1.fasta");
      exit(1);
    }

  int min_fiveprime_flank_anchor = 20;
  int min_threeprime_flank_anchor= 10;
  int max_anchor_span = 10000;
  int length_of_arrays=20000;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=10000;
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
						    0, NULL, NULL, NULL, NULL, NULL);

  
  CU_ASSERT(ret==0);



  hash_table_free(&hash_table);
  fclose(chrom_fptr);


}



void test_db_graph_make_reference_path_based_sv_calls_null_test_5()
{

  // ***************************************************************************************************************************************************
  // Reference is two copies of a single sequence, tandem repeat of about 36 bases. Individual is the same, with an Alu inserted between
  // Since the supernode in the individual has one copy of th repeat, and then the Alu, and then stops, you should be unable to find
  // an anchor at the Alu-end of the supernode. So this should find nothing.
  // ***************************************************************************************************************************************************


  int kmer_size = 31;
  int number_of_bits = 15;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;


  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_alu_insertion", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_without_alu.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_without_alu.fasta");
      exit(1);
    }


  int min_fiveprime_flank_anchor = 3;
  int min_threeprime_flank_anchor= 3;
  int max_anchor_span =400;
  int length_of_arrays=800;
  int min_covg =1;
  int max_covg = 100000;
  int max_expected_size_of_supernode=400;



  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
						    individual_edge_array,0,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, NULL,
						    0, NULL, NULL, NULL, NULL, NULL);
  
  CU_ASSERT(ret==0);


  hash_table_free(&hash_table);
  fclose(chrom_fptr);


}

//numbering of tests goes back to 1 - these unti tests are for cases where we DO expect to find a variant
void test_db_graph_make_reference_path_based_sv_calls_test_1()
{


  // ***************************************************************************************************************************************************
  // 1. Reference  = short sequence, individual is idential except for one base change
  // ***************************************************************************************************************************************************



  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_short_seq_with_one_base_difference", &bad_reads, hash_table);

  CU_ASSERT(seq_loaded==192);
  
  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/second_person_same_short_seq_one_base_diff.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/second_person_same_short_seq_one_base_diff.fasta");
      exit(1);
    }

  int min_fiveprime_flank_anchor = 5;
  int min_threeprime_flank_anchor= 5;
  int max_anchor_span =40;
  int length_of_arrays=80;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=40;


  //this time, we expect results. So  malloc some arrays to hold the results. Only use this option for testing, as in general use, we expect to find millions
  // of variants

  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 30);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 30);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 30);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 30);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';
  
  int return_variant_start_coords_array[2];
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  //now do the test!!!
  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test1", "w");
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
  						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);

  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("AATAGACGCCCACACCTGATAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCACA", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("AAGCCACA", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("CTGTACTTGTA", return_flank3p_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==23);



  //cleanup
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}







void test_db_graph_make_reference_path_based_sv_calls_test_2()
{


  // ***************************************************************************************************************************************************
  // 2. Reference = Alu-NNNNN- same Alu, and person is identical to reference except for a single base difference. This should find the single variant!
  //     It is also an example where we traverse the supernode in the reverse direction to that in which is appears in our array. This is 
  //      therefore a good unit test to ensure that our code works in this case
  //     Final comment - note that the coverage is double on a lot of these nodes because we have 2 copies of the Alu. But that doesn't in any way add credence
  //      to the variant call - it's paralogous sequence that doubles the coverage.
  // ***************************************************************************************************************************************************

  int kmer_size = 31;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_both_alu_Ns_alu_with_one_base_difference", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/one_person_aluNsalu_PLUS_SINGLE_BASE_CHANGE.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/one_person_aluNsalu_PLUS_SINGLE_BASE_CHANGE.fasta");
      exit(1);
    }

  int min_fiveprime_flank_anchor = 10;
  int min_threeprime_flank_anchor= 10;
  int max_anchor_span =350;
  int length_of_arrays=700;
  int min_covg =1;
  int max_covg = 10;
  int max_expected_size_of_supernode=350;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 400);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 400);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 400);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 400);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];  
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test2", "w");

  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array,return_variant_start_coords_array_ptr );
  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("GCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("CTGGGAGGATCGCTTGAGTCCAGGAGTTCTGG", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGG", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==59);

  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}




void test_db_graph_make_reference_path_based_sv_calls_test_3()
{
  // ***************************************************************************************************************************************************
  // 3. Reference is just a short sequence. Individual is missing 2 bases in the middle--> we should find this
  //     Note this test is also interesting because the reference and individual both contain a closed loop in their de Bruijn graph
  // ***************************************************************************************************************************************************


  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_2_bases_missing", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_without_2_bases_missing.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_without_2_bases_missing.fasta");
      exit(1);
    }


  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= 2;
  int max_anchor_span =40;
  int length_of_arrays=80;
  int min_covg =1;
  int max_covg = 1000;
  int max_expected_size_of_supernode=40;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 80);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 80);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];    
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test3", "w");


  //individual has 2 missing bases, reference does not
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr );
  
  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("CCACACCTGA", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("ACAC", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGACCCC", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("GACCCC", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==20);


  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);


  
}




void test_db_graph_make_reference_path_based_sv_calls_test_4()
{
  // ***************************************************************************************************************************************************
  // Reference is just a short sequence. Individual has an extra 2 bases in the middle--> we should find this
  // ***************************************************************************************************************************************************

  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;


  //we use the same two people as last time, but swap their roles
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_2_bases_missing", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_2_bases_missing.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_with_2_bases_missing.fasta");
      exit(1);
    }


  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= 2;
  int max_anchor_span =40;
  int length_of_arrays=80;
  int min_covg =1;
  int max_covg = 1000;
  int max_expected_size_of_supernode=40;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 80);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 80);
  return_branch2_array[0] = (char*) malloc( sizeof(char) * 80);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];      
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test4", "w");



  //cf previous test - we are using individua with index 1, not 0 as last time. ie we have swapped which is individual and which is reference
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr );

  


  fclose(fp);

  CU_ASSERT(ret==1);
  CU_ASSERT_STRING_EQUAL("CCACACCTGA", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("ACAC", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGACCCC", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("GACCCC", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==20);


  hash_table_free(&hash_table);
  fclose(chrom_fptr);

  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);


}



void test_db_graph_make_reference_path_based_sv_calls_test_5()
{


  // ***************************************************************************************************************************************************
  // 5. Reference is a single sequence which is a single supernode. Individual is the same, with an Alu inserted in the middle
  // Should be able to find anchors at both ends this time, and find the insertion
  // Note that by missing the alu in the reference, it has a load of nodes that won't be seen in the individual, so for example here
  // we will ask for a min 3' flank of 8. This actually means 31+8 bases of agreement.
  // ***************************************************************************************************************************************************


  int kmer_size = 31;
  int number_of_bits = 15;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;



  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_alu_inserted_mid_supernode", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_one_supernode_and_without_alu.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_with_one_supernode_and_without_alu.fasta");
      exit(1);
    }


  int   min_fiveprime_flank_anchor = 3;
  int   min_threeprime_flank_anchor=  8;
  int   max_anchor_span =500;
  int   length_of_arrays=1000;
  int   min_covg =1;
  int   max_covg = 100000;
  int   max_expected_size_of_supernode=500;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];        
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];  
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);


  
  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test5", "w");
  
  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
						    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);
  fclose(fp);




  CU_ASSERT(ret==1);

  //first base of Alu insert = first base of displaced sequence. ie we had xxxxGyyyy and then we insert Alu: xxxx[G-Alu]Gyyyy. So first base of Alu look slike is in flank5p
  CU_ASSERT_STRING_EQUAL("TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCAGG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("AAACTCCGTCTGTACTAATAGTACAAAAAT", return_trusted_branch_array[0]);
  CU_ASSERT_STRING_EQUAL("CCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGAAACTCCGTCTGTACTAATAGTACAAAAAT", return_branch2_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==40); //note the insertion happens at coordinate 39, but the first inserted base is the same as what would be there anyway - a G. 
  
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);




}



void test_db_graph_make_reference_path_based_sv_calls_test_6()
{


  // ***************************************************************************************************************************************************
  // 6. Reverse the roles of previous test. Reference has an Alu, individual does not. 
  // Should be able to find anchors at both ends this time, and find the deletion
  // ***************************************************************************************************************************************************


  int kmer_size = 31;
  int number_of_bits = 15;
  int bucket_size    = 10;
  long long  bad_reads = 0;
  int max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;


  //we use the same two people as last time, but swap their roles
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_alu_inserted_mid_supernode", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_alu_in_middle_of_supernode.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_with_alu_in_middle_of_supernode.fasta");
      exit(1);
    }



  int  min_fiveprime_flank_anchor = 3;
  int  min_threeprime_flank_anchor= 8;
  int  max_anchor_span =500;
  int  length_of_arrays=1000;
  int  min_covg =1;
  int  max_covg = 100000;
  int  max_expected_size_of_supernode=500;


  char**  return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char**  return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char**  return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char**  return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';

  int return_variant_start_coords_array[2];          
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];    
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  
  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test6", "w");


  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
						    individual_edge_array,0,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
                                                    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);
  fclose(fp);


  CU_ASSERT(ret==1);

  //first base of Alu insert = first base of displaced sequence. ie we had xxxxGyyyy and then we insert Alu: xxxx[G-Alu]Gyyyy. So first base of Alu look slike is in flank5p
  CU_ASSERT_STRING_EQUAL("TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("TAGCCAGG", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("AAACTCCGTCTGTACTAATAGTACAAAAAT", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("CCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGAAACTCCGTCTGTACTAATAGTACAAAAAT", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==40); //note the insertion happens at coordinate 39, but the first inserted base is the same as what would be there anyway - a G. 


  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}





void test_db_graph_make_reference_path_based_sv_calls_test_7()
{

  // ***************************************************************************************************************************************************
  // 7. Reference is an Alu with a different Alu inserted in the middle. Individual is just the first Alu - ie they have a deletion of an Alu from within an Alu.
  // 
  // ***************************************************************************************************************************************************


  int  kmer_size = 31;
  int  number_of_bits = 15;
  int  bucket_size    = 10;
  long long bad_reads = 0;
  int  max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;


  //we use the same two people as last time, but swap their roles
  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_is_alu_other_has_2nd_alu_inserted", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/person_with_alu_in_middle_of_alu.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/person_with_alu_in_middle_of_alu.fasta");
      exit(1);
    }


  
 //    reference is
 //   > AluJo#SINE/Alu inserted in middle of 7SLRNA#SINE/Alu
 //   GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
 //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
 //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
 //   GCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGA
 //   GGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACAT
 //   AGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCGGGCGTGGTGG
 //   CGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCT
 //   TGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACT
 //   CCAGCCTGGGCGACAGAGCGAGACCCTGTCTCAAAAAAAAAAAAAAAAAA
 //   AAAAAAAAAAAA
 //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
 //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
 //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
 //   AAAAAAAAAAAAAAAAAAAA
 //   ..plus N's on the end
    
    
 //   and individual is
    
 //   >7SLRNA#SINE/Alu
 //   GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG
 //   GTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCC
 //   GATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG
 //   CGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAA
 //   ACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGT
 //   GAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
 //   AAAAAAAAAAAAAAAAAAAA
    
    


  int  min_fiveprime_flank_anchor = 3;
  int  min_threeprime_flank_anchor= 3;
  int  max_anchor_span =500;
  int  length_of_arrays=1000;
  int  min_covg =1;
  int  max_covg = 100000;
  int  max_expected_size_of_supernode=500;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 500);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 500);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *500);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];          
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];    
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test7", "w");



  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
						    individual_edge_array,0,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table, fp,
                                                    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAG", return_flank5p_array[0]);
  CU_ASSERT_STRING_EQUAL("GAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAA", return_flank3p_array[0]);
  CU_ASSERT_STRING_EQUAL("CGGGGGACCACCAGGTTGCCTAAGGAGGGGT", return_branch2_array[0]);
  CU_ASSERT_STRING_EQUAL("GCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGAGGATCGCTTGAGCCCAGGAGTTCGAGACCAGCCTGGGCAACATAGCGAGACCCCGTCTCTACAAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGGATCGCTTGAGCCCAGGAGTTCGAGGCTGCAGTGAGCTATGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGGGGGACCACCAGGTTGCCTAAGGAGGGGT", return_trusted_branch_array[0]);
  CU_ASSERT(return_variant_start_coords_array[0]==151); 


  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);

}



void test_db_graph_make_reference_path_based_sv_calls_test_8()
{


  // ***************************************************************************************************************************************************
  // 8. Reference is 10 kb of chromosome 1 plus 1kb of sequence inserted mid-supernode, and person is identical, except for that 1kb is missing 
  // ***************************************************************************************************************************************************


  int  kmer_size = 31;
  int  number_of_bits = 15;
  int  bucket_size    = 10;
  long long bad_reads = 0;
  int  max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_1kb_deletion", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/first_person_10kb_chrom1_plus_1kb_inserted_mid_supernode.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/first_person_10kb_chrom1_plus_1kb_inserted_mid_supernode.fasta");
      exit(1);
    }

  // Let's be clear about what this test looks like. 

  //       The individual is 10kb of chromosome 1, which has the following supernodes

//>node_0  cylcle:false length:4546 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:1 lstf:2
//TGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACTCCGTCTGTACTAATAGTACAAAAATTAGCCAGGCGTGATGGTGTGCACCTGTAGTCCTTGCTACTCAGAAAGCTGAGGCAGGAGAATCGCTTGTACCCAGGAGGCAGAGGTTGCAGTGAGCAGAGATTGTGCCACTGCACTCCATCCTGGGTGACAGAGTGCTATGAGTCACCACACCTGGTATGAGCCACCGTGCCTGGCCCACAATGACTTTTACACATGTTGTTAAATCATCTTACAGATTTTATAATTTGGGGGAAGAAAAGTTTTACTAAATTGTCTTTTAATGGAAACTCTACAAGAACCAGAATCTTTGCTTTGTTCACTTATGTATCCATTCCTAGGCCTAGAAAAATGTCTGACACATAGCGGCAATTATTCATTGAATAAATGGACCCAGGGATAGTACATTAGCTATGCTATATGCATACATTAAAGATGTAGATTATCGACTTTCAAAAGATAATTAATGTAACTTCTTACTGCTTCTGAACATGTTTGTGAGTTATATTGCTGAGGGACCTTTATCTTCTCATTCTTTCATCTTAACCCAGTGTTATAAAATTGAAATCACCAATATTATTCCATATCTAAAATTAATATCTACCTTGTAAAAAATATCACTCTGCTGCATTTGAGAATAGACTTTTTAGGTAATAATGATGCAATCCATAGGGTTTTTTGGGGGCACAGAGGGATTCATGCTAACAGAACATTTTATTTTCTATTTTCCCAGAGCTGTAAAACATGAAATTACGGTAGTATAAGGCATATTTTTACTCTTTTTATAATTTTTTCTAAAAAAAATTAGTGTTTGTTCCCTATATAACTTTTAACTTTATAGGTAAATATTTGTTTCTTTCAGCTCCAGTTTTATGTGAAATAGAGTTTTCAGATTTATGTAGCATGGAAAGTTTTAATACGTCAGAGTTACTGATTTTTGCCAATCATTTTCTCAATTATTTCTTTTTTATCTTTAGTTGATTTTTTTGTAGTGACACATTTTGTTTCTAGTCTCATTTCCTTTTGTTTATATTCTATGTATATTTCATTTTTGGTTACTATGAGAATTACATATAACATCCTAGAGTTATAACATTTTAATTTGAATTTATTTCAACTTAAGTTCAATCACATACCAAAATTCTACTGCTATATATATAGCTCTACTCTTTTTATGTTATTGATGTGACAAATTATATCTTTATTCATTGTATACCAGCTAACAGATTTACAATTACATTTTATGCATTTGCCTTTTAAATTATGTAGAAAATAAAAAGCAGAGTTACAAACCAAAATTACAATAGGACTGTTTTTATGTTTGTTTATGTATTTACCTTTACCAGAGAGCTTTGTATATTCATACAGCTTGCTTATTTACTTATATAGTTATTGCCTAGAGTTCATTTATTTCAACCTGAAGGACTTAACACTTCCTGAATGTCAAATTCAGGGATAAATGGATTTTTTTCAGTTTTAAAAAAAAATCCGGAAATGTCTTAATTTCTCCTTCATTTTTGAAGGATAAGTTTTCCAGCTATATATTTCTCAATTGACAGGTTTCTTCATTATTTTAAATATATAATCCACTGCCTACTGGCCTTCAAGGTTTCTGCTGAGAAATCAGCTGCTAATGTTATCTGGATCCCTATCTGTGAGAGTTGCTCTTCTCTCTGAGTTTTCAACATTCTCCCATTATCTTTTTTTTGTTTGTTTTTGAGACAAATAATTGTACATATTCATGGGATACAGAGTGATATTTTGATACATGTATACAATGCCCAGTGATCAAATAAGGATAATTAGCGTATCCATCACCTCAAATAGTTGTCATTTATTTGTATTGTGAACAGTCAACATTCTTTCTTCTAGTTTTTTAAATTTATAAACATTTAAATTTTATTACAGAAATTTAAATTTTTTGATTCTGAAAAAGTCATATATGTATGCAACATTTTTTATCGTTTATTTATATATTTATGCATCTTTCCTTTTAGTTTTGACAGAGATTTTCTATTTTATCATTATTTCAAAAGAACTCTTACCTGTATTTATTTATCAAGTATATTTCCCTTGTTTTTTCCTAGTATATTAATTTATTTACTTATCTTCTAAAAATCCTCCATATAATCTGTTTATTTTGTTTCCTTTCTATAATTTCTTCAATAATTAGTTCTGTTCTATTTTCCATTAAAATATTTAAATCTTGTATGAATTTTTGTCAGATTAGAAATTTAGGGCGTTTCTTAATTTCTCTATACTCTAGCTTTTGACTTTTTTTTTCTGACCTAAGAGGTATTTAGAGCACATTTTAGATTTTTTATTTTGACTAATCATTTAAAATGTATACTAATCTTCAATTTAAATAAAAAACTGGTCTATAGTGACAAAAATTACAAATGAGCCTAACTAATAAATTATCAGCTGTGTTTATATGTATAAGCATGCACAGATTTTGGTAAATATGTACATAGTATATTGGTGAGCTTATTTTTATCATTCTTAACTCATTGTGTAGTCTAAACGTTGGGGAAAAAATAACATACAATAATCAGATGGTGTGAATAAGAAAATTGTTCTAATGTTTGTAAACCAAGCAACTGTTTTAACTGCTCCCCTCTTCCTGATTGACTTCTAAAAGGGATTGATCCATATTGGGTCCTATCATGTACGTCACGGTATAACATCTCCAGCTATAAAATGGAAATTTGAGAATAACTTTGCTGCTACTCAGATATATTTTATTTCAAAAACATACACTAAGGTGTTGCTGTTGGATCTTTCCAAAAACATATTCACACAGAACTTTCAATCACACTGAGCCATATTTGAACAATCTTTCAAGGTCAGCTCTGGCATAAGCTAACATTATACCATTTAACTCAGAAATTTCTTTAGTATTTGATTAATGGGTTTATGTTTGATATGTAATGTAATTTTCTAATACTAAATCAAGTGGTAATTTTGTTAGTCAAGTTGATTTAGTGGCTTGGGAAGAAAGCTTTTAATGTTCCCCTAATTTTTCTTTCCTTTGACATGATCCTTCACATGTCTTATTTTGCTTAGTGATTTTTCTTTTTTTTTTTTTTTTTTTGAGACAGGGTCTTACTCTACCACCCAGGCTTGAGTGCAGTGGTGCAATCACAGCTCATTGCAGCCTTGACCTCCCAGACTCAAGCTATTCTTCCACCTCAGCCTCCCAAGTAGCTGGTACTACAGGCACATGCCACCAAACTTGGCTAATTTTTGTATTTTTTGTAGAGACAGAGTTTTGCCAAATTCTCAGGCTGGTCTGGAATTTCTGGGCTCAAGTAATCCTGCCTTGGCCTCCCAACATGCTGATATTACAGACATAAGCCACAGTACCTGGCCAGTTTTCTTTTTAAAAAAATCTATTGGTTATTAATTTGAAGCCTTCCTTTTCATAGCTGTGCTCCTTAATTGGGAGCAAACATGAATGGACCACAACTTAGCCAATTTTCTATATACGATCTTTGCCATCCTAATTTAAAGGAATATTAATTCTTTCTTTTCCTCTTTCATTCCACAAACCTGTATTGACTACATCTAAGTTCTAAATGGTGCACTGGATGTTGAAAAAGTTGATGATGAGCAAGAACAAAATTCCTGCTTTCAGGAGACTTACAGTTCAATATGGGAAATATAATTTGTTAAAATATAAAAGTGCAATTGTGTTACATGCTGTACGAAGTACATGTTGACATGTGAGCATATAATAAATGGGCTGGAGGCCAGAGGATTGCCAAAGAGAATGGGCCTCCTGCTGAGATGAAAAGTTGAGCAGGGATTAGTTGGCAAAAGTGGAGGGACGATCCTTTCTAGGCAGGAGGAAGAACATGTACAGAATCTCTGAGGTGTGATGCGACAAAGTCTATATAAAAAACTGAAGAAAGGTCTAATATGGCTTAAATACAGAAGCTAGTAGGAGAGGAGTCGAAAAGAGGCTGGAGAAGTAGAAAGTGTCTGCATTCTGCAGGAACTTATATTGTATAAAATATATTGTATATTGTATAAAAAGAATTTCTCTTTATTCTAAGTGCAATGTGAAGCCAATGAAGTGCTTTAAACAGGTGATGTGATTTGATTGAATTTATTACTTCACTTAACAAATATTCATTACATGCCCACTGTTTGTCAGATATTGCTCTAGCCCCTGGTGATACAGTAGGGAATAAAACAGGCAAAAATCCCTGTCCTCTTGCAGCTTATAATGGACTGCAATGTTTAATATGTCAGAGGAGGTCCACGGAGGAGTGACTTCTAAGCAAGAATCTGAAAAAAATGAGGATATCTAAGGAGGGAACAAATGGTTCAAAAGCCCTATAATTGCAAGCAGGCATGATGAAGCAATTGCAGTTGTCCTGACTCTCAACACCGTGGAACTCAAAGGAGATGGAAAGATTCCTTCTCTCCCTCATATATTTTCTCTCTTTCTGTCTATATATATAGAATATGAGACATTTCCCTAATCATTATGT
//>node_1  cylcle:false length:1696 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//ACATAATGATTAGGGAAATGTCTCATATTCTTCTACTCAGAAATAAGCAATATAGCAATTACTGTTTTTTACATTTTACAGTTACAGTTTCAGAGAAAGTTTGATATTTATCTAAAATTTTTCAATGTATGAACTTTTTCATTTGACAAACCATAATTGTACATATTCTTGGGATACAGAGTGATATTTTCTTTACATGTATAGAATGTGTAGTGATCAAATCAGGGTAATTTCCACTAATTTAAAATGCCACCTTTATGTTATTGTAATTTATATATATACTATATATATACACACACACATATATATATATACATGTCCACATACAGTGTGTGTGTGCACATGTACACACATGCATATGTGTATATAATGCCCAGTATAAGCAATGTGCACAAATAAAATTAGCTAACAGAGATAGTATAGAGTGAGAGGAGAGGCAGATTAATCTTTGAGGAAAAGCACAATTTTATAGCTGAATGGAGAAAGCTGAGGTGGTTTCTAAGATGGAGAATAAGACGAAAAATGTAAGTACGTTGTCTGACTGAATTCAAGAAAGAAGGGTAAAAGAGAAGAAAGTAGTGGTCTTATCATTAAATGCCACAGAGAGGTAAAGATAAAGACAACATATTGTTTTGGGTTTAGTAATTTAAGGGTTACCAAATTCCGTTTTGGAGGAGGAACAGATTCCATGTCCACTAGAATGGAATGAACAAGAAATGGAGGAGGAAAATAGGTAGTTTTTCAAAAGTTTTCAAAAATATGAAAAGAAGAAATGAAGTGGTACTTGGAAGAGATTGTTGAAATGGGAGAGACTATGGTGGCTTGTTTAGAAGCAGTTGAGATAGATCCAATTGAGATAGAGATATTGACTATATAAACAAAAGAATGACAAATTAATAGTGTAATGGATAACTTGACTTTGGCAAATATTGTGAATTTTTGTGAAAGTACAACTAAAAGGCAATGTCACTCCAATAATCACCAGAGTAATCAATTTGCTTATTGCTGTCCCTTTAAATATAGTTCTCTGGTATCAACTAACATGTTTTTAACTAATGATGCTTCTTAAAGAAAAGGGAAAAGACCTTTTTCTTTCTTTCAGTCTTCAATGATTCACTGCTTCATCTCGCTCCACCAAAGATAAATGAAATCTACATCTCTTATACATTAACAATGCATGACAATTTACAAATAGCTAAATTTTTGGAGCTAACTTTAAGTACCTGAATGGAATTTAATCAACCCACTAATCTCCTTCTCACTTCTCAGTTATTTATCAAGTTTATGTCAAGGGACAAGGAAAAATTATCCAAACATTGTTTAAAACAATCATCATTAATTAGTAACACTTATCCAGGGGGGTTTTTAACCTTTCCCCCACTCAAGGATTATTCTAATGTCAGAGTAGAATAAAAAATAAGTGCAGTGATGCTGACTCTTCCAAGCTTAACATTTCTCACAAGTCAATTAGCTTTGTACTGGGAGGAGGGCGTGAAGGGCTGCTTGCGGTAGTTGTGTAGCAGCAGCACAATGGCCGCAGACAAGGAAAACAGTTTCTAGGAATTCCTCGTATATAATTTTATATTTTTGACAAGATTAATGACCCATGCTCCCTTCCTCTCCATTTCTTTTTTTGGAATTCTGTTGGTATGTAGTTACTATATTTTATTAAAGGAAATTAGCCTTATCTCTTATT
//>node_2  cylcle:false length:566 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//AGGAAATTAGCCTTATCTCTTATTATATTTTTTATGACCTTCAAAGTAGTGTCTCTGCTTAAAAGTGTACCCTGGCCGGGCGTGGTGGCTCACACCTGTAATTCCAGCACTTTGGGAGGCCGAGGCGGGTGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTGTACTAAAAATACAAAAAATTAGCAGGGCATAGTGGCGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTCAGGCAGGAGAATGGCGTGAACCCGGGAGACGGAGCTTGCGGTGAGCTGAGATCGCACCGCTGCACTCCAGCCTGGGCGACAGAGCAAGACTCCGTCTCAAAAAAAAAAAAAAAGTATACCCTGAGGCACACATCAAGCGACATGTAGAGTTCATAAATTCTGGCCAAATGGTCATACCTCAAACCTCATCAGCAGTAAGGCTCTTTACTTGCACTGACAAATATGAACGCTGGGGAATTTGGAAATGATATATAATATATAATATTATATATATAATAGATATATAATATATAATATATATAATACATAT
//>node_3  cylcle:false length:2989 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:1 fst:0 lst_coverage:2 lstr:1 lstf:2
//AATTTGTTACAGCAGTAATAGAACAAGTGGTTATCCATATGAGGCAAATTAGATTGGATACCTATCTCCAATAGAAATCAATTCAAGGTGAATTCCAGGAAAATACTTAAAACATTTAGATTAAAAATAAATGAGAATTTTTGTTACTTTTGGTAGGTCATAGAACCAAGAAAAACAAACATTAAGGAGGAAAAATGAACATATGACTACATCAAAATATAAAGCTTCTCTATTTGGAAGATATCATAAGGTGACAAATCATAAACTGTAATATTTACAACATATATATAAGTGAATAAATATACATTTAGAATATATATGAACTCCCAAAAATCAACAGGAAAAATAAGACATAGAACAAGCAAAATGCATAAACAAAAGAAGGCAAAACAAAAATAATGACTCATAATTATATGAAAAGAAGCTCATCTTCATAGATGAGCAGATAAATGCAAATTAAAACCACCCTGAGATGCTTTTTACATCCATGAGCCTGATAAAAGTTAGAGTCTAAAAGTAATAACAAAGATGGGAAGTAATAGAAAATCTTGTCCATTACTGGTTAAAGTATAAACTGATACAGCTACTTTATAGAATATTACATTATAGAATAAAGTTGTGAGTATGTATATGCAGTGACTCAGCATATTCATTGCTAGTATGTACTCAAGAGAAACTTACAGGAGTGGACTAGGAAGTAAATACAAAATGATTACAACATTGTTTGTTATATCAAAAAATAAAAAAGACACCCAATTTTCCAGCAAAAAAAATAAGTAAAAATAAATCCTGGTGTATTCTAACAATGGAATAATATATAGCCATTAAAATAAATCAACTATTACTGTACATATGAATGTAAATATCAGCAAAACATATTGTTTAGTGAAAAACTAAGAAGCTGAAGAAGAATATATACAATATGGTTACATTTATATGAAGTCCAAAAACTTGCAAAATAAAGAAATGTATTTAGAAATAGATTCACATGTGAGAAAACTAGAAGAAAATTAATGAAAGGATAAGAGGGATAGCAGTAATTCTGAGTAGTTGAGGGAATTTCAATTGGAAAAAAATAATATCATATTCTTTAAGTCAGGTAGTGGGTATTAGCATTTGTTTTACCATCGTTCTTTATTCTTATAGCTACACTATATATTTTCAATGTATTTAATGTATTTTTTGCATAATTAAATATTATGCAATAAAAATGAGAAAACAAAAAAGTAGAAAATGATAAATTACAATAAAGAAATGGAGAAAAAATTATAATCTAGTTGAGTAATGGTATATTACATAGCTATTTTCTTAAGTAGATGTATGTACATGATGTATGCACGATTGTACATACATGTTCTTAATTATATATAAATATATATGTACATATTTTTAATATAAAATACTAAACAAAGTACACCAAAATATTAGCTCCTATGTTAGTGAGATAATGTTTTGTTTTTTTGTATTTTAAGTTTTACATAGTAGGTGTATTTTTCTGTTTTCATACTGCTATAAAGAACTGCCCAAGACTGGGTAATTTATAAAGGAAAGAAGTTTAATTGGCTCACCGTTCAGCACAGCTTGGGAGGCCTCAGGAAATCTACAATCATGGCGGAAGACAAAGAGGAAGCAAGCCAGCTTCTTCGCAAGGCAGCATGAAGAAGTGCCGAGCAAAGGGGAAAGAATCCCTTATAAAACCATCAAATCTCGTGAGAACTCACTATCACAAGAACAGCACAGGGGAAACTGCCCCCATGATTCAATTACCTCCACCTGGTCTCTCCCTTGACCTGTGGGGATTATGGGGACTATGGGGATTACAATTCAAGATGAGATTCAGGTGGGGATACAAAGCCTAACCATATCAGTAGGCATGTATTGAATTTTAAACTCAGAGAAAAATACTAGTGTTTTTATAGGATTCTTACTAAAGAAAAACCAGAAAGTAATAAACCATCTACGCTAAGACATAAAATTCAGTTGTTTAGTTACAAGATAGAATGTGGCCTTGTAAGAAAGCAAATTAACTTCTAACATACAAAGCCTTAGAGAAGATTCAAGTGACTGACGGATCTTAAACAGAGCTATTATTACAACTCGAACTGCAGTAAAATATCCTCAGCAACATAGATGTGTGTGTTTCACTAGTCAGAGCAATACAAATTTAATGAAACTCCACTGGTGGTGTTTTTAATCAGACAATTTCTGAAGATGTCCTGGCTTATTCACAGATGCAAGCCAAATCTCTAGAAGAGTACCATAATAAGAAAAAAAAGAATACAGGCAATTGAGAGCTGTTCCAAAGTTTAGGGAGTTTTTGTAAGGAATTAATAAATAAAAATGTTCTTGAAAGACAGAAATTAATATGCAGTTCATACTGTCAGAATTGCAGGCAATTTATCAAAGTCCCCTAATCCTCCAAAATCGCTATTTTTTTTTTGACACACACTTTACAGTACAGAAGAAAATGTCTCCGGCAATAAATCACAAAGTTAAAATTACCTAGTCTACAATTAACTACACAGTGATGGTAAATCATTTTCTACCAAAAGAAAGAAATGTCTTGTCTATTCAGGTTCTGCTCTACTTAAAAGTTTTCCTTGTTGGCGAGCAAGTGGTTAGAAAATTATATTTTATACGTACATTCAGCTTAACTATCATTCAGCTCAGGAAGATGACTCAGGGCCTTATCCATACCTTCAAGTTTGCTCTTAGCAAGTAATTGTTTCAGTATCTATATCAAAAATGGCTTAAGCCTGCAACATGTTTCTGAATGATTAACAAGGTGATAGTCAGTTCTTCATTGAATCCTGGATGCTTTATTTTTCTTAATAAGAGGAATTCATATGGATCAGCTAGAAAAAAATTAAGAGGAAAATCACATGGAAAGTTATATATTATATATCTATTATATATAATATATATCTATTACATATTATATATTGTATATCTATTACATATATATTATATATGTATTATATATATTATA
//>node_4  cylcle:true length:1691 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:2
//TGGCCAGGCTGGTCTCGAACTCCTGACCTCAAATGATCCATCTACCTCGGCCTCCCAAAGTGCTGGAATTACAGATGTGAGCCACAATGCCCGGCCTTATTTTCTACAACTTTGGTAACTTTAGCATATACCCCAAATCTGTAAGACATAATATTATAATTCAAATGCAACTCATGGCTTCTCTTTGTACTCTTTCTCTAGCTTTTGAATTATTTATTCTAATACCAGTTTTAATTCTGACACAAAATCATGGGAGTTCTAATCAAAATCCAACCTTTTATCATAAAAACTATGAAGAAATTATGAGTAGAATTTAAAAAGGAAAATAGGCCTATTAATTAGATTTGTCTTTGTAGCATTTAACTCTATAATAAATAATATTTTATGCCTATGAGTCCCCAACAAAGCCTCCAGCTTCTATTTAGATATAAACTGTAAAAGTCACTACTGGATCCACAAGCAAGACTATGGTAAATAAATTTCTCCACCTAACCAGCTTCTTTTACATGATGTTACATGTTTCTTTTGTTTTTTCATTTTGGCAAATATTGATTGTCATCTTCGTGTTTGTCTATGTCCTAAGTGCTGGGATACAGAATCTGAAAAGATGGACACAGGACCTGCCTTCAAGTTCACCCCCTTTTTTTTTTTTTTTTGAGATGCAGTTTTGCTCTTGTCACCCAGGCTGGAGTGTACTGGTGAGATCTCTGCTCACTGCAACCTCCACCTTCAGGGTTCAAGTGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGTGATTACAGGTCCCAGCCACCACGCCTAGCTAATTTTTGTATTTTTAGTAGAGACAGCGTTTCATCATGCTGGTCAGGCTGGTCTCGAACTCCTAACCTCAGGTAGTCGACCCACCTCGGCCTCCCACAGTGCTGAGATTACAGGCATGAGCCACCACGCCCTGCTAGGAGTTCACGCTTTAGTTGGGGAAAATATACAATAAGCAAGCCAGTTTTTAAAATGAGAACTGCAATTAGAGTTAAATGCTACAAAGACAAACTCACAGGAAGATGGGATGTAGAATGATAAGGCTCTCAGAATAGTAAGAGAAACTATTGCTTCTTACGATGTTTGTCTTTCTTTGTATCGGTGCTCAGCTGAGTCTGCAGTGCTTCAGAGGCAGCTTTCATTTTATAAAAATCTATGATTTCTCCTTCCAGTTTTTTTTTCTCTTCCTCGAGCTTCCTTATCTCCTCCTGTTGAATCATTTTAAGATGCTCGAACTTGTCCTGCAGCTGTGAAACCAATGTGCAGTTGTGACACCAAAGCAGTGTGGCTGAACACCTAAAAGAATACGCTTTTTTTCTGATTATCAAACAAACCCAAATCATCACAGTAGACCACGATCTTAATAACAATCTCAAAAACTCAGGAGTAAACACTCAGATATGGAATTTTTCTTTTCTTTCTTTTTTCCTTTTATAAGATGGAGTCTCACTCTGTTGCCCAGGCTGGAGTGCACTGGTGCGATCTCAGCTCACTGCAACCTCCATCTCCCAGTTCAAGTGATTCTCCTGCCTCAGCCTCTTGAGTAGCTGGGACTATAGGCATGCACCACCACTACAGGCGTGTGCCACCACACCTGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTTGCCATGATGGCCAGGCTGGTCTCGAACTCCTGACCTCA
//>node_5  cylcle:false length:462 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:2 lst_coverage:1 lstr:0 lstf:1
//TGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCCTCCCGCTTTGGCCTCCCAAAGACTTTTTTTTTTTTTTTAATATAGAGACAAGTTCTCAGTACGTTGCCCAGGCTGGTCTCAAACTCCTGAGCTCAAGTGATCCTCCCACCTCAGCTTCCCAAAGTGCTGGGACTGACTGGATGCAGTGGCTCATGCTTGTAAACTCAGCACTTTGGGAGGCCAAGGTGGGAGGATCGCTTGAGCCCAGGAGTTCAAGACCAGACTGGGTGATATAACACAATAGTCAACTTCAACAGGAGAGAGAATCTGTAAACTTGAATATAGATCTTCCGAAATTATCCAGTCAGTGGACAGAGAAAAAAAGAATAAAAGAGAGAAAAGAAGGCTGGGTGTGGTGGCTCAAGCCTGTAATCCCAACACTTTGGGAGGCCGAGGCAGGCAGATTAAGAGGTCAGGAGTTCA
//>node_6  cylcle:false length:116 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//AATAAGAGATAAGGCTAATTTCCTTTAATAATAATAAAATCCTTTAATAAAAATATAAAGGAATAATATAATAATTTTCTTTAATAAAATATAATAAGAGATAAGGCTAATTTCCT
//>node_7  cylcle:false length:41 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//ATATATAATATATAATATATATAATACATATATAATATATA
//>node_8  cylcle:false length:38 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//AAAATATAATAAGAGATAAGGCTAATTTCCTTTAATAA
//>node_9  cylcle:false length:48 min_coverage:1 max_coverage:1 average_coverage: 1.00 fst_coverage:2 fst:1 lst_coverage:2 lstr:1 lstf:2
//TATAATATATATAATACATATATAATATATAATATATATAATACATAT
//>node_10  cylcle:false length:100 min_coverage:2 max_coverage:2 average_coverage: 2.00 fst_coverage:2 fst:2 lst_coverage:2 lstr:2 lstf:1
//AGAATATGAGACATTTCCCTAATCATTATGTGTAATTACAATTACATATATATATGTAATTGTAATTACACATAATGATTAGGGAAATGTCTCATATTCT


//The reference is the same as the individual, except it has about 1kb inserted in node 0. To be clearer, the reverse complement of node 0 is what is in the 10kb without insertion:

//GTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATCTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCAATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAATGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA  TTGCACCACTGCACTCAAGCCTGGGTGGTAGAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCTTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTAAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGC

//  and the insert happens just after TTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA, where I have put a space above


//The inserted sequence is

// ACCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACAT
// which I took from chromosome 2.






//NOTE - the inserted sequence starts with an A. The sequence that follows the inserted sequence also starts with an A. So our algorithm will add that A to the flanking region.
// We allow for this is checking results match expectations



  int  min_fiveprime_flank_anchor = 21;
  int  min_threeprime_flank_anchor= 21;
  int  max_anchor_span =8000;
  int  length_of_arrays=16000;
  int  min_covg =1;
  int  max_covg = 100000000;
  int  max_expected_size_of_supernode=8000;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *5600);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];            
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];      
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test8", "w");





  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table,  fp,
                                                    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("CATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGCAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTATTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA", return_flank5p_array[0]);

  CU_ASSERT_STRING_EQUAL("GAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGTGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAAAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCACTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGGTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCA", 
			 return_flank3p_array[0]);

  CU_ASSERT_STRING_EQUAL("TTGCACCACTGCACTCAAGCCTGGGTGGTA", return_branch2_array[0]);

  CU_ASSERT_STRING_EQUAL("CCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACATATTGCACCACTGCACTCAAGCCTGGGTGGTA", return_trusted_branch_array[0]);

  CU_ASSERT(return_variant_start_coords_array[0]==6722); 



  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);


}

void test_db_graph_make_reference_path_based_sv_calls_test_9()
{

  // ***************************************************************************************************************************************************
  // 9. Identical to previous test, but each person has identical 600 lines of chromosome 12 added before the sequence that was there in previous tes
  //     This is purely to test that the start coordinate of the variant is found correctly, even when you have to load more and more sequence into the array.
  // ***************************************************************************************************************************************************
  

  int  kmer_size = 31;
  int  number_of_bits = 15;
  int  bucket_size    = 10;
  long long   bad_reads = 0;
  int  max_retries=10;

  dBGraph* hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }
  long long seq_loaded=0;

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/variations/two_people_one_with_1kb_deletion_both_with_600lineschrom12_beforehand", &bad_reads, hash_table);

  FILE* chrom_fptr = fopen("../data/test/pop_graph/variations/first_person_600lineschrom12_then_10kb_chrom1_plus_1kb_inserted_mid_supernode.fasta", "r");
  if (chrom_fptr==NULL)
    {
      printf("Cannot open ../data/test/pop_graph/variations/first_person_600lineschrom12_then_10kb_chrom1_plus_1kb_inserted_mid_supernode.fasta");
      exit(1);
    }


  int min_fiveprime_flank_anchor = 21;
  int  min_threeprime_flank_anchor= 21;
  int  max_anchor_span =8000;
  int  length_of_arrays=16000;
  int  min_covg =1;
  int  max_covg = 100000000;
  int  max_expected_size_of_supernode=8000;


  char** return_flank5p_array = (char**) malloc( sizeof(char*) *2);
  char** return_flank3p_array = (char**) malloc( sizeof(char*) *2);
  char** return_trusted_branch_array = (char**) malloc( sizeof(char*) *2);
  char** return_branch2_array = (char**) malloc( sizeof(char*) *2);
  
  if ( (return_flank5p_array==NULL) || (return_flank3p_array==NULL) || (return_trusted_branch_array==NULL) || (return_branch2_array==NULL) )
    {
      printf("Failed to alloc return_something_array - cannot start test\n");
      exit(1);
    }

  return_flank5p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_flank3p_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_trusted_branch_array[0] = (char*) malloc( sizeof(char) * 5600);
  return_branch2_array[0] = (char*) malloc( sizeof(char) *5600);

  if ( (return_flank5p_array[0]==NULL )||(return_flank3p_array[0]==NULL )|| (return_trusted_branch_array[0] ==NULL) || (return_branch2_array[0]==NULL))
     {
      printf("Failed to malloc the [0] entry of one of the return arrays. Cannot start test.");
      exit(1);
    }
  
  return_flank5p_array[0][0]='\0';
  return_flank3p_array[0][0]='\0';
  return_trusted_branch_array[0][0]='\0';
  return_branch2_array[0][0]='\0';


  int return_variant_start_coords_array[2];              
  return_variant_start_coords_array[0]=0;
  return_variant_start_coords_array[1]=0;

  int* return_variant_start_coords_array_ptr[2];        
  return_variant_start_coords_array_ptr[0]=&(return_variant_start_coords_array[0]);
  return_variant_start_coords_array_ptr[1]=&(return_variant_start_coords_array[1]);

  FILE* fp = fopen("../bin/temp_outputfile_trustedpath_sv_caller_test9", "w");





  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
						    individual_edge_array,1,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, hash_table,  fp,
                                                    1, return_flank5p_array, return_trusted_branch_array, return_branch2_array, return_flank3p_array, return_variant_start_coords_array_ptr);


  fclose(fp);


  CU_ASSERT(ret==1);

  CU_ASSERT_STRING_EQUAL("CATAATGATTAGGGAAATGTCTCATATTCTATATATATAGACAGAAAGAGAGAAAATATATGAGGGAGAGAAGGAATCTTTCCATCTCCTTTGAGTTCCACGGTGTTGAGAGTCAGGACAACTGCAATTGCTTCATCATGCCTGCTTGCAATTATAGGGCTTTTGAACCATTTGTTCCCTCCTTAGATATCCTCATTTTTTTCAGATTCTTGCTTAGAAGTCACTCCTCCGTGGACCTCCTCTGACATATTAAACATTGCAGTCCATTATAAGCTGCAAGAGGACAGGGATTTTTGCCTGTTTTATTCCCTACTGTATCACCAGGGGCTAGAGCAATATCTGACAAACAGTGGGCATGTAATGAATATTTGTTAAGTGAAGTAATAAATTCAATCAAATCACATCACCTGTTTAAAGCACTTCATTGGCTTCACATTGCACTTAGAATAAAGAGAAATTCTTTTTATACAATATACAATATATTTTATACAATATAAGTTCCTGCAGAATGCAGACACTTTCTACTTCTCCAGCCTCTTTTCGACTCCTCTCCTACTAGCTTCTGTATTTAAGCCATATTAGACCTTTCTTCAGTTTTTTATATAGACTTTGTCGCATCACACCTCAGAGATTCTGTACATGTTCTTCCTCCTGCCTAGAAAGGATCGTCCCTCCACTTTTGCCAACTAATCCCTGCTCAACTTTTCATCTCAGCAGGAGGCCCATTCTCTTTGGCAATCCTCTGGCCTCCAGCCCATTTATTATATGCTCACATGTCAACATGTACTTCGTACAGCATGTAACACAATTGCACTTTTATATTTTAACAAATTATATTTCCCATATTGAACTGTAAGTCTCCTGAAAGCAGGAATTTTGTTCTTGCTCATCATCAACTTTTTCAACATCCAGTGCACCATTTAGAACTTAGATGTAGTCAATACAGGTTTGTGGAATGAAAGAGGAAAAGAAAGAATTAATATTCCTTTAAATTAGGATGGCAAAGATCGTATATAGAAAATTGGCTAAGTTGTGGTCCATTCATGTTTGCTCCCAATTAAGGAGCACAGCTATGAAAAGGAAGGCTTCAAATTAATAACCAATAGATTTTTTTAAAAAGAAAACTGGCCAGGTACTGTGGCTTATGTCTGTAATATCAGCATGTTGGGAGGCCAAGGCAGGATTACTTGAGCCCAGAAATTCCAGACCAGCCTGAGAATTTGGCAAAACTCTGTCTCTACAAAAAATACAAAAATTAGCCAAGTTTGGTGGCATGTGCCTGTAGTACCAGCTACTTGGGAGGCTGAGGTGGAAGAATAGCTTGAGTCTGGGAGGTCAAGGCTGCAATGAGCTGTGA", return_flank5p_array[0]);

  CU_ASSERT_STRING_EQUAL("GAGTAAGACCCTGTCTCAAAAAAAAAAAAAAAAAAAGAAAAATCACTAAGCAAAATAAGACATGTGAAGGATCATGTCAAAGGAAAGAAAAATTAGGGGAACATTAAAAGCTTTCTTCCCAAGCCACTAAATCAACTTGACTAACAAAATTACCACTTGATTTAGTATTAGAAAATTACATTACATATCAAACATAAACCCATTAATCAAATACTAAAGAAATTTCTGAGTTAAATGGTATAATGTTAGCTTATGCCAGAGCTGACCTTGAAAGATTGTTCAAATATGGCTCAGTGTGATTGAAAGTTCTGTGTGAATATGTTTTTGGAAAGATCCAACAGCAACACCTTAGTGTATGTTTTTGAAATAAAATATATCTGAGTAGCAGCAAAGTTATTCTCAAATTTCCATTTTATAGCTGGAGATGTTATACCGTGACGTACATGATAGGACCCAATATGGATCAATCCCTTTTAGAAGTCAATCAGGAAGAGGGGAGCAGTTAAAACAGTTGCTTGGTTTACAAACATTAGAACAATTTTCTTATTCACACCATCTGATTATTGTATGTTATTTTTTCCCCAACGTTTAGACTACACAATGAGTTAAGAATGATAAAAATAAGCTCACCAATATACTATGTACATATTTACCAAAATCTGTGCATGCTTATACATATAAACACAGCTGATAATTTATTAGTTAGGCTCATTTGTAATTTTTGTCACTATAGACCAGTTTTTTATTTAAATTGAAGATTAGTATACATTTTAAATGATTAGTCAAAATAAAAAATCTAAAATGTGCTCTAAATACCTCTTAGGTCAGAAAAAAAAAGTCAAAAGCTAGAGTATAGAGAAATTAAGAAACGCCCTAAATTTCTAATCTGACAAAAATTCATACAAGATTTAAATATTTTAATGGAAAATAGAACAGAACTAATTATTGAAGAAATTATAGAAAGGAAACAAAATAAACAGATTATATGGAGGATTTTTAGAAGATAAGTAAATAAATTAATATACTAGGAAAAAACAAGGGAAATATACTTGATAAATAAATACAGGTAAGAGTTCTTTTGAAATAATGATAAAATAGAAAATCTCTGTCAAAACTAAAAGGAAAGATGCATAAATATATAAATAAACGATAAAAAATGTTGCATACATATATGACTTTTTCAGAATCAAAAAATTTAAATTTCTGTAATAAAATTTAAATGTTTATAAATTTAAAAAACTAGAAGAAAGAATGTTGACTGTTCACAATACAAATAAATGACAACTATTTGAGGTGATGGATACGCTAATTATCCTTATTTGATCACTGGGCATTGTATACATGTATCAAAATATCACTCTGTATCCCATGAATATGTACAATTATTTGTCTCAAAAACAAACAAAAAAAAGATAATGGGAGAATGTTGAAAACTCAGAGAGAAGAGCAACTCTCACAGATAGGGATCCAGATAACATTAGCAGCTGATTTCTCAGCAGAAACCTTGAAGGCCAGTAGGCAGTGGATTATATATTTAAAATAATGAAGAAACCTGTCAATTGAGAAATATATAGCTGGAAAACTTATCCTTCAAAAATGAAGGAGAAATTAAGACATTTCCGGATTTTTTTTTAAAACTGAAAAAAATCCATTTATCCCTGAATTTGACATTCAGGAAGTGTTAAGTCCTTCAGGTTGAAATAAATGAACTCTAGGCAATAACTATATAAGTAAATAAGCAAGCTGTATGAATATACAAAGCTCTCTGGTAAAGGTAAATACATAAACAAACATAAAAACAGTCCTATTGTAATTTTGGTTTGTAACTCTGCTTTTTATTTTCTACATAATTTAAAAGGCAAATGCATAAAATGTAATTGTAAATCTGTTAGCTGGTATACAATGAATAAAGATATAATTTGTCACATCAATAACATAAAAAGAGTAGAGCTATATATATAGCAGTAGAATTTTGGTATGTGATTGAACTTAAGTTGAAATAAATTCAAATTAAAATGTTATAACTCTAGGATGTTATATGTAATTCTCATAGTAACCAAAAATGAAATATACATAGAATATAAACAAAAGGAAATGAGACTAGAAACAAAATGTGTCACTACAAAAAAATCAACTAAAGATAAAAAAGAAATAATTGAGAAAATGATTGGCAAAAATCAGTAACTCTGACGTATTAAAACTTTCCATGCTACATAAATCTGAAAACTCTATTTCACATAAAACTGGAGCTGAAAGAAACAAATATTTACCTATAAAGTTAAAAGTTATATAGGGAACAAACACTAATTTTTTTTAGAAAAAATTATAAAAAGAGTAAAAATATGCCTTATACTACCGTAATTTCATGTTTTACAGCTCTGGGAAAATAGAAAATAAAATGTTCTGTTAGCATGAATCCCTCTGTGCCCCCAAAAAACCCTATGGATTGCATCATTATTACCTAAAAAGTCTATTCTCAAATGCAGCAGAGTGATATTTTTTACAAGGTAGATATTAATTTTAGATATGGAATAATATTGGTGATTTCAATTTTATAACACTGGGTTAAGATGAAAGAATGAGAAGATAAAGGTCCCTCAGCAATATAACTCACAAACATGTTCAGAAGCAGTAAGAAGTTACATTAATTATCTTTTGAAAGTCGATAATCTACATCTTTAATGTATGCATATAGCATAGCTAATGTACTATCCCTGGGTCCATTTATTCAATGAATAATTGCCGCTATGTGTCAGACATTTTTCTAGGCCTAGGAATGGATACATAAGTGAACAAAGCAAAGATTCTGGTTCTTGTAGAGTTTCCATTAAAAGACAATTTAGTAAAACTTTTCTTCCCCCAAATTATAAAATCTGTAAGATGATTTAACAACATGTGTAAAAGTCATTGTGGGCCAGGCACGGTGGCTCATACCAGGTGTGGTGACTCATAGCACTCTGTCACCCAGGATGGAGTGCAGTGGCACAATCTCTGCTCACTGCAACCTCTGCCTCCTGGGTACAAGCGATTCTCCTGCCTCAGCTTTCTGAGTAGCAAGGACTACAGGTGCACACCATCACGCCTGGCTAATTTTTGTACTATTAGTACAGACGGAGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCA", 
			 return_flank3p_array[0]);

  CU_ASSERT_STRING_EQUAL("TTGCACCACTGCACTCAAGCCTGGGTGGTA", return_branch2_array[0]);

  CU_ASSERT_STRING_EQUAL("CCTCAATCTATAACCACTACCTTCTGGGTCCTTTCTAAAAATTGACAAATAATAATCATATATAATTAATGTACAATGTTATGTTTCAATACATGTTTGCATTGTGGAATAATTAAATCAAGCTACTTGGCATGTCAGTAACATCACATGCTTAGCATTTTTGTAGTGAAAATATTTAAAATCTACTCTTTTAGCAATTTTGAAATATACAATACAGTACTTACTCACTTAACATCATTGATTGGTCCTCAGAAACTGCAACTTGGAGTGAAAAGATGTATAAAGAAACCAATTTTCCCATAGGCTAATAGCTATAAATAAGAGTTAGATTCTTACAGCATATTTTTGGTCACAAAATATCACCAAACTTCTAAATAAAGACCAAAACACTTCAAATATTAAACATTGAAATAAATATGAGCTTTGCATACATTTAAGAAAGATTAATAAAAACAAGTAAGATAATTATTTGCCCAATTATTTCATTCAGGGTTGGGGAGACTGGAGTCTGTGCTGGAAGCTCAGGGCTCAAGCTGGGCAACAGCCCTGGACAGGATGCCATCCCACTGCAGGATGGCTCACACATGCCCACAGCCACTCAGCCTGGGACCATTTGGACACAGCAATTAACCTTACCTGCATGTCTTTGTGGGGAGGAAACCAGAGTGCTTAGAAAAACCCATGCAGACAGACACAGAGCAAACATGCAAACCTCACAAAGATATTGTTTCTTCTGTCACCTGTGCTTTTGGGTCATATTCAAGAAATCATTAACCAAATAAAAGTCGTGGAGCTTTTCCCTATGTTTTCTTTTAGTAGTTTTATAGTTTCAGGTCTTACATTTAACTCCTTAATCCATTTTGATTTTTGCATATGGTGTGAGATAAGCTTCTGGTTTCATTCTCCCACATGTGGATATCCAGTTCTCTGAACACCATATATGGAAGAGACTGTCATTTCCTCATGATATGTTCCTGGCACTTTTGTTGAAATCAATTGACCATAGATCTGTGGGTTTATTTCTGGCTTTTTATTCTGTTCCATTGGCCAATGTACCTGTGTTTATGCTTGTGCCTTGCTGTTTTGATTATTATAGCTTTATAATATGTTTTGAAATCAGGTAGTGTGATGCCTCCATCTTTGCTTTTTATGCTCAAGATAGTTTGGATATTCAGAGTGTTTTATGGTTCCATATACATATTGCACCACTGCACTCAAGCCTGGGTGGTA", return_trusted_branch_array[0]);


  CU_ASSERT(return_variant_start_coords_array[0]==42662); 
  
  hash_table_free(&hash_table);
  fclose(chrom_fptr);
  
  free(return_flank5p_array[0]);
  free(return_flank3p_array[0]);
  free(return_trusted_branch_array[0]);
  free(return_branch2_array[0]);
  free(return_flank5p_array);
  free(return_flank3p_array);
  free(return_trusted_branch_array);
  free(return_branch2_array);





}




void test_get_covg_of_nodes_in_one_but_not_other_of_two_arrays()
{


 //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 8;
  int bucket_size    = 8;
  long long bad_reads = 0;
  int max_retries=10;

  
  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
  int seq_loaded = load_population_as_fasta("../data/test/pop_graph/one_person_two_reads", &bad_reads, hash_table);

  //>read1
  //AAAACGAAAAAATTCGAG
  //>read2
  //TCGAGAAAACGAGA

  //We will make two arrays of nodes, one corresponding to each read, and compare these two arrays. There is no reason to be associated with a read, this is just to make the test easier.
  // Where the nodes come from is irrelevant

  dBNode* n1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAC", kmer_size),kmer_size), hash_table); //seen in both reads
  dBNode* n2 = hash_table_find(element_get_key(seq_to_binary_kmer("AAACG", kmer_size),kmer_size), hash_table); //covg 2
  dBNode* n3 = hash_table_find(element_get_key(seq_to_binary_kmer("AACGA", kmer_size),kmer_size), hash_table); //covg 2
  dBNode* n4 = hash_table_find(element_get_key(seq_to_binary_kmer("ACGAA", kmer_size),kmer_size), hash_table); //covg 1
  dBNode* n5 = hash_table_find(element_get_key(seq_to_binary_kmer("CGAAA", kmer_size),kmer_size), hash_table); //covg 1
  dBNode* n6 = hash_table_find(element_get_key(seq_to_binary_kmer("GAAAA", kmer_size),kmer_size), hash_table); //seen in both reads
  dBNode* n7 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA", kmer_size),kmer_size), hash_table); //covg 2
  dBNode* n8 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAA", kmer_size),kmer_size), hash_table); //covg 2
  dBNode* n9 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAT", kmer_size),kmer_size), hash_table); //covg1 
  dBNode* n10 = hash_table_find(element_get_key(seq_to_binary_kmer("AAATT", kmer_size),kmer_size), hash_table); //covg 1
  dBNode* n11 = hash_table_find(element_get_key(seq_to_binary_kmer("AATTC", kmer_size),kmer_size), hash_table); //covg 1
  dBNode* n12 = hash_table_find(element_get_key(seq_to_binary_kmer("ATTCG", kmer_size),kmer_size), hash_table); //covg 1
  dBNode* n13 = hash_table_find(element_get_key(seq_to_binary_kmer("TTCGA", kmer_size),kmer_size), hash_table); //covg1
  dBNode* n14 = hash_table_find(element_get_key(seq_to_binary_kmer("TCGAG", kmer_size),kmer_size), hash_table); // seen in both reads

  dBNode* e1 = hash_table_find(element_get_key(seq_to_binary_kmer("TCGAG", kmer_size),kmer_size), hash_table); //seen in both reads
  dBNode* e2 = hash_table_find(element_get_key(seq_to_binary_kmer("CGAGA", kmer_size),kmer_size), hash_table); //covg 2
  dBNode* e3 = hash_table_find(element_get_key(seq_to_binary_kmer("GAGAA", kmer_size),kmer_size), hash_table); //covg 1
  dBNode* e4 = hash_table_find(element_get_key(seq_to_binary_kmer("AGAAA", kmer_size),kmer_size), hash_table); //covg 1
  dBNode* e5 = hash_table_find(element_get_key(seq_to_binary_kmer("GAAAA", kmer_size),kmer_size), hash_table); //seen in both reads
  dBNode* e6 = hash_table_find(element_get_key(seq_to_binary_kmer("AAAAC", kmer_size),kmer_size), hash_table); //seen in both reads
  //note we haven't gone all the way to the end of the second read



  CU_ASSERT(n1 != NULL);
  CU_ASSERT(n2 != NULL);
  CU_ASSERT(n3 != NULL);
  CU_ASSERT(n4 != NULL);
  CU_ASSERT(n5 != NULL);
  CU_ASSERT(n6 != NULL);
  CU_ASSERT(n7 != NULL);
  CU_ASSERT(n8 != NULL);
  CU_ASSERT(n9 != NULL);
  CU_ASSERT(n10 != NULL);
  CU_ASSERT(n11 != NULL);
  CU_ASSERT(n12 != NULL);
  CU_ASSERT(n13 != NULL);
  CU_ASSERT(n14 != NULL);

  CU_ASSERT(e1 != NULL);
  CU_ASSERT(e2 != NULL);
  CU_ASSERT(e3 != NULL);
  CU_ASSERT(e4 != NULL);
  CU_ASSERT(e5 != NULL);
  CU_ASSERT(e6 != NULL);


  

  dBNode* array1[] = {n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14};
  dBNode* array2[] = {e1,e2,e3,e4,e5,e6};


  int num_nodes_in_array1_not_array2 = 0;
  int num_nodes_in_array2_not_array1 = 0;

  int coverages_of_nodes_in_array1_not_array2[14];
  int coverages_of_nodes_in_array2_not_array1[6];
  
  int* ptrs_to_coverages_of_nodes_in_array1_not_array2[14];
  int* ptrs_to_coverages_of_nodes_in_array2_not_array1[6];

  dBNode* working_array1[20];
  dBNode* working_array2[20];

  int i;
  for (i=0; i<14; i++)
    {
      coverages_of_nodes_in_array1_not_array2[i]=0;
      ptrs_to_coverages_of_nodes_in_array1_not_array2[i]=&coverages_of_nodes_in_array1_not_array2[i];
      
    }
  for (i=0; i<6; i++)
    {
      coverages_of_nodes_in_array2_not_array1[i]=0;
      ptrs_to_coverages_of_nodes_in_array2_not_array1[i]=&coverages_of_nodes_in_array2_not_array1[i];

    }

  get_covg_of_nodes_in_one_but_not_other_of_two_arrays(array1, array2, 14,6,&num_nodes_in_array1_not_array2 , &num_nodes_in_array2_not_array1, 
						       ptrs_to_coverages_of_nodes_in_array1_not_array2, ptrs_to_coverages_of_nodes_in_array2_not_array1, 
						       working_array1, working_array2, individual_edge_array,0);

  CU_ASSERT(num_nodes_in_array1_not_array2==10 );
  CU_ASSERT(num_nodes_in_array2_not_array1==3 );

  //we also get the coverages of these nodes, in some random order (depends on the mem addresses of the nodes)
  //Note that in advance, we don't know how many such nodes there will be, if any
  qsort(coverages_of_nodes_in_array1_not_array2, num_nodes_in_array1_not_array2, sizeof(int), int_cmp); 
  qsort(coverages_of_nodes_in_array2_not_array1, num_nodes_in_array2_not_array1, sizeof(int), int_cmp); 


  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[0]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[1]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[2]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[3]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[4]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[5]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[6]==1);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[7]==2);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[8]==2);
  CU_ASSERT(coverages_of_nodes_in_array1_not_array2[9]==2);

  CU_ASSERT(coverages_of_nodes_in_array2_not_array1[0]==1);
  CU_ASSERT(coverages_of_nodes_in_array2_not_array1[1]==1);
  CU_ASSERT(coverages_of_nodes_in_array2_not_array1[2]==2);


  hash_table_free(&hash_table);
  
}
