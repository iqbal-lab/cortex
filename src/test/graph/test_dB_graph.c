#include <CUnit.h>
#include <Basic.h>
#include <dB_graph.h>
#include <element.h>
#include <binary_kmer.h>
#include <file_reader.h>
#include <test_dB_graph.h>
#include <stdlib.h>


void test_hash_table_find()
{

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0; 

  int max_retries=10;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  //Load the following fasta:
  //    >read1
  //    AAAAAAAA
  //    >read2
  //    GGCT
  //    >read3
  //    TAGG


  int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/test_dB_graph.fasta",&bad_reads, 20, db_graph);

  //length of total sequence
  CU_ASSERT(seq_length == 16);
  
  //number of kmers
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 5);

  //bad reads
  CU_ASSERT_EQUAL(bad_reads, 0);

  //all the kmers and their reverse complements from the reads
  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;
  
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);

  //some kmers that should not be in the graph

  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);

  //kmers in the graph
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

  //kmers not in the graph
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


void test_tip_clipping()
{

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0;

  int min_coverage, max_coverage;
  double avg_coverage;
  boolean is_cycle;
  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  

  //Load the following fasta:
  
  //>main_trunk
  //GCGTCCCAT
  //>tip
  //CGTTT

  int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_tip.fasta",&bad_reads, 20, db_graph);


  CU_ASSERT_EQUAL(seq_length,14);

  dBNode* node1 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);

  CU_ASSERT(node1 != NULL);
  
  int tip_length = db_graph_db_node_clip_tip(node1, 10,&db_node_action_set_status_pruned,db_graph);
  
  CU_ASSERT_EQUAL(tip_length,2);

  dBNode* node2 = hash_table_find(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* node3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* node4 = hash_table_find(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);


  //check the change in the status
  CU_ASSERT(db_node_check_status(node1,pruned));
  CU_ASSERT(db_node_check_status(node2,pruned));
  
  //check status didn't change
  CU_ASSERT(db_node_check_status(node3,none));
  CU_ASSERT(db_node_check_status(node4,none));

  //check tip clip works well with db_graph_supernode (ie that the tip was really removed)
  dBNode * nodes_path[100];
  Orientation orientations_path[100];
  Nucleotide labels_path[100];
  char tmp_seq[100+1];

  int length_supernode = db_graph_supernode(node3,100,&db_node_action_set_status_visited,
					    nodes_path,orientations_path,labels_path,
					    tmp_seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
					    db_graph);

  CU_ASSERT_EQUAL(length_supernode,6);
  CU_ASSERT_STRING_EQUAL("GGACGC",tmp_seq);
  

  //check ends
  dBNode* node5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* node6 = hash_table_find(element_get_key(seq_to_binary_kmer("CAT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);

  //check the status are correct
  
  CU_ASSERT(db_node_check_status(node1,pruned));
  CU_ASSERT(db_node_check_status(node2,pruned));

  CU_ASSERT(db_node_check_status(node3,visited));  
  CU_ASSERT(db_node_check_status(node4,visited));
  CU_ASSERT(db_node_check_status(node5,visited));
  CU_ASSERT(db_node_check_status(node6,visited));
   

  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
  

}

void test_node_prunning_low_coverage()
{

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0;
  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  

  //Load the following fasta:
  
  //>main_trunk
  //GCGTCCCAT
  //>tip
  //CGTTT


  int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_tip.fasta",&bad_reads, 20, db_graph);

  CU_ASSERT_EQUAL(seq_length,14);

  dBNode* node1 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* node2 = hash_table_find(element_get_key(seq_to_binary_kmer("TCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  dBNode* node3 = hash_table_find(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
  

  boolean node_pruned = db_graph_db_node_prune_low_coverage(node1,1,&db_node_action_set_status_pruned,db_graph);

  CU_ASSERT(node_pruned);
  CU_ASSERT(db_node_check_status(node1,pruned));

  CU_ASSERT(!db_node_edge_exist(node3,Guanine,reverse));
  CU_ASSERT(db_node_edge_exist(node3,Thymine,forward));

  CU_ASSERT(!db_node_edge_exist(node2,Cytosine,reverse));
  CU_ASSERT(db_node_edge_exist(node2,Cytosine,forward));

  //check all arrows were removed from node1
  void nucleotide_action(Nucleotide n){
    CU_ASSERT(!db_node_edge_exist(node1,n,reverse));
    CU_ASSERT(!db_node_edge_exist(node1,n,forward));
  }
  
  nucleotide_iterator(&nucleotide_action);
  hash_table_free(&db_graph);

  CU_ASSERT(db_graph == NULL);
}


void test_get_perfect_path() //test db_graph_get_perfect_path
 {
 
  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits=4;
  int bucket_size   = 10;
  int seq_length;
  long long bad_reads = 0;
  dBNode * path_nodes[100];
  Orientation path_orientations[100];
  Nucleotide path_labels[100];
  char tmp_seq[100];

  double  avg_coverage;
  int min_coverage, max_coverage;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a dB_graph.
  //         The test then picks an element in the graph, and calls get_perfect_path
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****


  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta", &bad_reads, 20,  db_graph);

   CU_ASSERT_EQUAL(seq_length,6);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),2);
   CU_ASSERT_EQUAL(bad_reads,0);
   

  // now start at GTA and get all the sequence from there to the end of the supernode, and see
  // if that is right.

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size,&tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
   CU_ASSERT(test_element1!=NULL);
   
   boolean is_cycle=false;

   
   int test1_length = db_graph_get_perfect_path(test_element1,forward,100,
						&db_node_action_do_nothing,
						path_nodes,path_orientations,path_labels,
						tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
						&is_cycle,db_graph);

    CU_ASSERT(is_cycle);
    CU_ASSERT_EQUAL(test1_length,4);
    
   CU_ASSERT_STRING_EQUAL(tmp_seq,"CGTA");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);
 
 
   /*   // **** */
   /*   //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end */
   /*   //   caused by two outward/exiting edges */
   /*   // **** */
   
   //first set up the hash/graph
   kmer_size = 3;
   number_of_bits= 4;
   bucket_size   = 10;

   bad_reads = 0;
   
   db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
   
   seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta",&bad_reads,20,db_graph);

   CU_ASSERT_EQUAL(seq_length,13);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),5);
   CU_ASSERT_EQUAL(bad_reads,0);

   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
   CU_ASSERT(test_element1!=NULL);

   //ACA < TGT so forward gives TT
   is_cycle=false;
   

   //go forward
   int test2_length = db_graph_get_perfect_path(test_element1,forward,100,
						&db_node_action_do_nothing,
						path_nodes,path_orientations,path_labels,
						tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
						&is_cycle,db_graph);

   CU_ASSERT(!is_cycle);
   CU_ASSERT_EQUAL(test2_length,2);
   CU_ASSERT_STRING_EQUAL(tmp_seq,"TT");


   //ACA < TGT so backward gives ""
   
   //go backward
   test2_length = db_graph_get_perfect_path(test_element1,reverse,100,
					    &db_node_action_do_nothing,
					    path_nodes,path_orientations,path_labels,
					    tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
					    &is_cycle,db_graph);

  
   CU_ASSERT(!is_cycle);
   CU_ASSERT_EQUAL(test2_length,0);

   
   CU_ASSERT_STRING_EQUAL(tmp_seq,"");
   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


   // ****
   //1.3 Fasta file that generate a graph with one long supernode, with a conflict at the end
   //   caused by two INWARD edges in the opposite direction
   // ****
   
   //first set up the hash/graph
   kmer_size = 3;
   number_of_bits = 4;
   bucket_size    = 5;
   bad_reads = 0;

   db_graph = hash_table_new(number_of_bits, bucket_size,10,kmer_size);
   

   seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta", &bad_reads, 20,  db_graph);

   CU_ASSERT_EQUAL(seq_length,13);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),5);
   CU_ASSERT_EQUAL(bad_reads,0);
   

   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
   CU_ASSERT(test_element1!=NULL);

   //ACA < TGT so forward gives TT
   is_cycle=false;

   //go forward
   int test3_length = db_graph_get_perfect_path(test_element1,forward,100,
						&db_node_action_do_nothing,
						path_nodes,path_orientations,path_labels,
						tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
						&is_cycle,db_graph);

   CU_ASSERT(!is_cycle);

   CU_ASSERT_EQUAL(test3_length,2);
   CU_ASSERT_STRING_EQUAL(tmp_seq,"TT");

   
   //ACA < TGT so backwar gives ""
   is_cycle=false;

   //go reverse
   test3_length = db_graph_get_perfect_path(test_element1,reverse,100,
					    &db_node_action_do_nothing,
					    path_nodes,path_orientations,path_labels,
					    tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
					    &is_cycle,db_graph);

   CU_ASSERT(!is_cycle);

   CU_ASSERT_EQUAL(test3_length,0);
   CU_ASSERT_STRING_EQUAL(tmp_seq,"");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


   // ****
   //1.4 Fasta file that generate a graph with an infinite loop at a single kmer
   //
   // ****

  
   //first set up the hash/graph
   kmer_size = 3;
   number_of_bits=8;
   bucket_size   =4;

   bad_reads = 0;

   db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

   seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_infinite_loop.fasta",&bad_reads,30,db_graph);

   CU_ASSERT_EQUAL(seq_length,25);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),1);
   CU_ASSERT_EQUAL(bad_reads,0);

   test_element1 = hash_table_find(seq_to_binary_kmer("AAA", kmer_size, &tmp_kmer1) ,db_graph);
   CU_ASSERT(test_element1!=NULL);

   //forward
   is_cycle=false;

   int test4_length = db_graph_get_perfect_path(test_element1,forward,100,
						&db_node_action_do_nothing,
						path_nodes,path_orientations,path_labels,
						tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
						&is_cycle,db_graph);
   CU_ASSERT(is_cycle);//this time it should find a cycle
   CU_ASSERT_STRING_EQUAL(tmp_seq,"A");
   CU_ASSERT_EQUAL(test4_length,1);

   //backward
   is_cycle=false;

   test4_length = db_graph_get_perfect_path(test_element1,reverse,100,
					    &db_node_action_do_nothing,
					    path_nodes,path_orientations,path_labels,
					    tmp_seq, &avg_coverage, &min_coverage, &max_coverage,
					    &is_cycle,db_graph);

   CU_ASSERT(is_cycle);//this time it should find a cycle
   CU_ASSERT_EQUAL(test4_length,1);
   CU_ASSERT_STRING_EQUAL(tmp_seq,"T");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


   
   // ****
   // 1.5 check parameters (path nodes,labels,etc) for get_perfect_path
   //
   // ****
   
   kmer_size = 3;
   number_of_bits = 4;
   bucket_size = 4;
   db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
   
   boolean found1, found2, found3;
   dBNode * node1;
   dBNode * node2;
   dBNode * node3;
   dBNode * node4;

   
   node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CGT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found1, db_graph);
   node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("GTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found2, db_graph);
   db_node_add_edge(node1, node2, reverse,reverse, db_graph->kmer_size);

   node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found3, db_graph);
   db_node_add_edge(node2, node3, reverse,reverse, db_graph->kmer_size);
 
   node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TAG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found3, db_graph);
   db_node_add_edge(node3, node4, reverse,reverse, db_graph->kmer_size);
  

   //add coverage
   element_update_coverage(node1,10);
   element_update_coverage(node2,2);
   element_update_coverage(node3,1);
   element_update_coverage(node4,20);
  


   dBNode * nodes[10];
   Orientation orientations[10];
   Nucleotide bases[10];
   int test5_length = 0;
 

   test5_length = db_graph_get_perfect_path(node1,reverse, 10,
					    &db_node_action_set_status_visited,
					    nodes, orientations, bases,
					    tmp_seq,&avg_coverage,&min_coverage,&max_coverage,
					    &is_cycle,db_graph);

  CU_ASSERT_EQUAL(test5_length,3);
  
  //check nodes
  CU_ASSERT_EQUAL(node1,nodes[0]);
  CU_ASSERT_EQUAL(node2,nodes[1]);
  CU_ASSERT_EQUAL(node3,nodes[2]);
  CU_ASSERT_EQUAL(node4,nodes[3]);

  //check labels
  CU_ASSERT_EQUAL(bases[0],Thymine);
  CU_ASSERT_EQUAL(bases[1],Adenine);
  CU_ASSERT_EQUAL(bases[2],Guanine);
  
  //check orientations
  CU_ASSERT_EQUAL(orientations[0],reverse);
  CU_ASSERT_EQUAL(orientations[1],reverse);
  CU_ASSERT_EQUAL(orientations[2],reverse);

  //check statuses
  CU_ASSERT(db_node_check_status(nodes[0], none));
  CU_ASSERT(db_node_check_status(nodes[1], visited));
  CU_ASSERT(db_node_check_status(nodes[2], visited));
  CU_ASSERT(db_node_check_status(nodes[3], none));

  //check coverage
  CU_ASSERT_EQUAL(avg_coverage,1.5);
  CU_ASSERT_EQUAL(min_coverage,1);
  CU_ASSERT_EQUAL(max_coverage,2);
  
  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
}


void test_writing_reading_binary(){

  dBNode node1, node2;
  int kmer_size = 3;
  boolean test_read;

  char seq[kmer_size];
  
  BinaryKmer tmp_kmer;

  element_initialise(&node1,seq_to_binary_kmer("AAA",kmer_size, &tmp_kmer),kmer_size);
  node1.edges = 'a'; // ie 64 = (0110 0100)2
  node1.coverage = 10;

  FILE* fp1 = fopen("../data/test/graph/dump_element.ctx", "w");
  db_node_print_binary(fp1,&node1);
  fclose(fp1);

  FILE* fp2 = fopen("../data/test/graph/dump_element.ctx", "r");
  
  test_read = db_node_read_binary(fp2,kmer_size,&node2);

  CU_ASSERT_EQUAL(test_read, true);

  CU_ASSERT(binary_kmer_comparison_operator(node1.kmer, node2.kmer));
  CU_ASSERT_EQUAL(node1.edges, node2.edges);

  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(&(node1.kmer),kmer_size,seq));
  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(&(node2.kmer),kmer_size,seq));
  
  CU_ASSERT_EQUAL(node1.edges, node2.edges);
  CU_ASSERT_EQUAL(node1.edges, 'a');
  CU_ASSERT_EQUAL(node2.edges, 'a');
 
  CU_ASSERT_EQUAL(node2.coverage, 10);

  test_read = db_node_read_binary(fp2,kmer_size,&node2);
  CU_ASSERT_EQUAL(test_read, false);

  fclose(fp2);
}
   


void test_detect_and_smoothe_bubble(){
  int kmer_size = 3;
  int number_of_buckets = 4;
  int bucket_size = 3;
  dBGraph * db_graph = hash_table_new(number_of_buckets,bucket_size,10,kmer_size);
  boolean found;
  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

 
  dBNode * node1, * node2, * node3, * node4, * node5, * node6, * node7, * node8, * node9;
  Orientation orientations1[100],orientations2[100];
  Nucleotide labels1[100], labels2[100];
  dBNode * path_nodes1[100], * path_nodes2[100];

  int length1, length2;
  char seq1[101], seq2[101];

  boolean bubble;

  double avg_coverage1,avg_coverage2;
  int min_coverage1,max_coverage1, min_coverage2,max_coverage2;
  
  //check a perfect bubble 

  //start point
  node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);

  //branch1
  node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);


  //modify coverage -- used below to clip one branch
  element_update_coverage(node2,5);
  CU_ASSERT_EQUAL(element_get_coverage(node2),5);


  node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CAC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
  node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("ACT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);

  //end point (branch merge here)
  node5 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);

  //branch 2
  node6 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);

  element_update_coverage(node6,1);
  CU_ASSERT_EQUAL(element_get_coverage(node6),1);

  node7 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
  node8 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);

  //add 3p extension
  node9 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),  &found, db_graph);
 
  //branch1
  db_node_add_edge(node1, node2, forward,forward, db_graph->kmer_size);
  db_node_add_edge(node2, node3, forward,forward, db_graph->kmer_size);
  db_node_add_edge(node3, node4, forward,forward, db_graph->kmer_size);
  db_node_add_edge(node4, node5, forward,reverse, db_graph->kmer_size);

  
  //branch2
  db_node_add_edge(node1, node6, forward,reverse, db_graph->kmer_size);
  db_node_add_edge(node6, node7, reverse,forward, db_graph->kmer_size);
  db_node_add_edge(node7, node8, forward,reverse, db_graph->kmer_size);
  db_node_add_edge(node8, node5, reverse,reverse, db_graph->kmer_size);

  CU_ASSERT(db_node_edge_exist(node1,Thymine,forward));
  CU_ASSERT(db_node_edge_exist(node1,Adenine,forward));

  //add 3p extension
  db_node_add_edge(node5, node9, reverse,reverse, db_graph->kmer_size);


  bubble = db_graph_detect_bubble(node1,forward,db_graph->kmer_size+1,
				  &db_node_action_set_status_visited,
				  &length1,path_nodes1,orientations1,labels1,
				  seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
				  &length2,path_nodes2,orientations2,labels2,
				  seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
				  db_graph);


  CU_ASSERT(bubble);

  CU_ASSERT(((labels1[0] == Adenine) && (labels2[0] == Thymine)) || ((labels2[0] == Adenine) && (labels1[0] == Thymine)));

  CU_ASSERT_EQUAL(labels1[1],Cytosine);
  CU_ASSERT_EQUAL(labels1[2],Thymine);
  CU_ASSERT_EQUAL(labels1[3],Thymine);

  CU_ASSERT_EQUAL(labels2[1],Cytosine);
  CU_ASSERT_EQUAL(labels2[2],Thymine);
  CU_ASSERT_EQUAL(labels2[3],Thymine);



  //check extensions
  
  Orientation orientations5p[100];
  Orientation orientations3p[100];
  Nucleotide labels_flank5p[100];
  Nucleotide labels_flank3p[100];
  boolean is_cycle5p, is_cycle3p;
  int length_flank5p,length_flank3p;
  dBNode * nodes5p[100];
  dBNode * nodes3p[100];
  char tmp_seq[101];

  double avg_coverage;
  int min_coverage, max_coverage;
  
  length_flank5p = db_graph_get_perfect_path(node1, reverse,100,
					     &db_node_action_set_status_visited,
					     nodes5p,orientations5p, labels_flank5p,
					     tmp_seq,&avg_coverage,&min_coverage,&max_coverage,
  					     &is_cycle5p,db_graph);

  CU_ASSERT_EQUAL(length_flank5p,0);
  CU_ASSERT_EQUAL(nodes5p[0],node1);


  length_flank3p = db_graph_get_perfect_path(node5, reverse,100,
					     &db_node_action_set_status_visited,
					     nodes3p,orientations3p, labels_flank3p,
					     tmp_seq,&avg_coverage,&min_coverage,&max_coverage,
					     &is_cycle3p,db_graph);
 
  
  CU_ASSERT_EQUAL(length_flank3p,1);
  CU_ASSERT_EQUAL(nodes3p[0],node5);
  CU_ASSERT_EQUAL(nodes3p[1],node9);
  CU_ASSERT_EQUAL(labels_flank3p[0],Adenine);

  //check statuses
  CU_ASSERT(db_node_check_status(node1,none));
  CU_ASSERT(db_node_check_status(node2,visited));
  CU_ASSERT(db_node_check_status(node3,visited));
  CU_ASSERT(db_node_check_status(node4,visited));
  CU_ASSERT(db_node_check_status(node6,visited));
  CU_ASSERT(db_node_check_status(node7,visited));
  CU_ASSERT(db_node_check_status(node8,visited));
  //this last node shouldn't be marked as visited
  CU_ASSERT(db_node_check_status(node5,none));


  //check the bubble from node5 perspective
  bubble = db_graph_detect_bubble(node5,forward,db_graph->kmer_size+1,
				  &db_node_action_set_status_visited,
				  &length1,path_nodes1,orientations1,labels1,
				  seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
				  &length2,path_nodes2,orientations2,labels2,
				  seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
				  db_graph);
  

  CU_ASSERT(bubble);
  CU_ASSERT(((labels1[0] == Adenine) && (labels2[0] == Thymine)) || (((labels2[0] == Adenine) && (labels1[0] == Thymine))));
  

  //test smooth bubble

  //remove a branch
  //this one should fail -- because of limit 0 on min coverage
  boolean removed = db_graph_db_node_smooth_bubble(node1,forward,db_graph->kmer_size+1,1,
						   &db_node_action_set_status_pruned,db_graph);
  CU_ASSERT(!removed);


  //this one should work
  removed = db_graph_db_node_smooth_bubble(node1,forward,db_graph->kmer_size+1,10,
					   &db_node_action_set_status_pruned,db_graph);
  
  CU_ASSERT(removed);
  //check that correct path was removed

  //the first node shouldn't be marked
  CU_ASSERT(db_node_check_status(node1,none));
  CU_ASSERT(db_node_check_status(node2,visited));
  CU_ASSERT(db_node_check_status(node3,visited));
  CU_ASSERT(db_node_check_status(node4,visited));
  CU_ASSERT(db_node_check_status(node6,pruned));
  CU_ASSERT(db_node_check_status(node7,pruned));
  CU_ASSERT(db_node_check_status(node8,pruned));
  //the last node shouldn't be marked 
  CU_ASSERT(db_node_check_status(node5,none));

  CU_ASSERT(!db_node_edge_exist(node1,Thymine,forward));
  CU_ASSERT(db_node_edge_exist(node1,Adenine,forward));

  //check arrows by getting a supernode

  dBNode * nodes_path[100];
  Orientation orientations_path[100];
  Nucleotide labels_path[100];
  boolean is_cycle;

  int length_supernode = db_graph_supernode(node2,100,&db_node_action_set_status_visited,
  					    nodes_path,orientations_path,labels_path,
  					    tmp_seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
  					    db_graph);
  
  CU_ASSERT_STRING_EQUAL(tmp_seq,"ACTTA");
  
  
  length_supernode = db_graph_supernode(node7,100,&db_node_action_set_status_visited,
					nodes_path,orientations_path,labels_path,
					tmp_seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
					db_graph);
  
  CU_ASSERT_STRING_EQUAL(tmp_seq,"");
  

  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);

  //test bubbles with unequal branch sizes
   
  kmer_size = 21;
  number_of_buckets = 4; //2**4 buckets
  bucket_size = 40;
  db_graph = hash_table_new(number_of_buckets,bucket_size,10,kmer_size);
  
  long long bad_reads;
 
  int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generate_bubble_with_unequal_branch_sizes.fa",&bad_reads, 200,  db_graph);
  
  CU_ASSERT_EQUAL(seq_length,343);

  //fetch the kmer where the path splits (bubble appears) GCCAACCATGCCTGTTAAGGG
  node1 = hash_table_find(element_get_key(seq_to_binary_kmer("GCCAACCATGCCTGTTAAGGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
 
  CU_ASSERT(node1!=NULL);

  bubble = db_graph_detect_bubble(node1,reverse,100,
				  &db_node_action_set_status_visited,
				  &length1,path_nodes1,orientations1,labels1,
				  seq1,&avg_coverage1,&min_coverage1,&max_coverage1,
				  &length2,path_nodes2,orientations2,labels2,
				  seq2,&avg_coverage2,&min_coverage2,&max_coverage2,
				  db_graph);


  CU_ASSERT(bubble);
  
  CU_ASSERT((length1 == 21 && length2==22) || (length2 == 21 && length1==22));

  if (length1==21)
    {
      CU_ASSERT_EQUAL(labels1[0],Guanine);
      CU_ASSERT_STRING_EQUAL(seq1,"GGGTTATTTTTCTAGAGAGTT");
      CU_ASSERT_EQUAL(labels2[0],Thymine);
      CU_ASSERT_STRING_EQUAL(seq2,"TGGGTTATTTTTCTAGAGAGTT");
    }
  else{
    CU_ASSERT_EQUAL(labels2[0],Guanine);
    CU_ASSERT_STRING_EQUAL(seq2,"GGGTTATTTTTCTAGAGAGTT");
    CU_ASSERT_EQUAL(labels1[0],Thymine);
    CU_ASSERT_STRING_EQUAL(seq1,"TGGGTTATTTTTCTAGAGAGTT");
  }
    
  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);

}


void test_db_graph_db_node_has_precisely_n_edges_with_status(){ 
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

  boolean found;
  dBNode * node1, * node2, * node3, * node4;
  dBNode  * next_node[4];
  Orientation next_orientation[4];
  Nucleotide next_base[4];

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
  node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
  node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);

 

  db_node_add_edge(node1, node2, forward,reverse, db_graph->kmer_size);
  db_node_add_edge(node1, node3, forward,forward, db_graph->kmer_size);

  boolean one_edge1 = db_graph_db_node_has_precisely_n_edges_with_status(node1,forward,none,1,
									 next_node,next_orientation,next_base,db_graph);


  CU_ASSERT_EQUAL(one_edge1,false);

  db_node_set_status(node2,visited);

  boolean one_edge2 = db_graph_db_node_has_precisely_n_edges_with_status(node1,forward,none,1,
									 next_node,next_orientation,next_base,db_graph);

 
  CU_ASSERT_EQUAL(one_edge2,true);
  CU_ASSERT_EQUAL(node3,next_node[0]);
  CU_ASSERT_EQUAL(next_orientation[0],forward);
  CU_ASSERT_EQUAL(next_base[0],Guanine);
  

  node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2), &found, db_graph);
  db_node_add_edge(node1, node4, forward, forward, db_graph->kmer_size);
   
  boolean one_edge3 = db_graph_db_node_has_precisely_n_edges_with_status(node1,forward,none,2,
									 next_node,next_orientation,next_base,db_graph);
  

  
  CU_ASSERT_EQUAL(one_edge3,true);

  //observe that the tests below require to know in wich order the edges are visited

  CU_ASSERT_EQUAL(node4,next_node[0]);
  CU_ASSERT_EQUAL(node3,next_node[1]);

  CU_ASSERT_EQUAL(next_base[0],Adenine);
  CU_ASSERT_EQUAL(next_base[1],Guanine);
  
  CU_ASSERT_EQUAL(next_orientation[0],forward);
  CU_ASSERT_EQUAL(next_orientation[1],forward);

  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
}




void test_get_N50()
{

  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 5;
  int bucket_size    = 4;
  int max_retries = 10;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  int seq_length=0;
  //long long count_kmers = 0;
  long long bad_reads = 0;


  //1. One fasta file containing two reads:
  //  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA which becomes just a one node supernode
  //  CACGGTATTATTTACCT which is a 13 node supernode.
  // so the N50 is 13
  
   int max_read_length=70;
   seq_length = load_fasta_from_filename_into_graph("../data/test/graph/n50_example1.fasta",&bad_reads, max_read_length,  db_graph);


   CU_ASSERT_EQUAL(seq_length,58);
   //CU_ASSERT_EQUAL(count_kmers,14);
   CU_ASSERT_EQUAL(bad_reads,0);

   CU_ASSERT(db_graph_get_N50_of_supernodes(db_graph)==13);

   //clean-up
   hash_table_free(&db_graph);



   //***************************************
   // example 2:

   // >read1
   // AAAAAA
   // >read2
   // CCCCCC
   // >read3
   // ACGTA
   // >read4
   // CTATG
   // >read5
   // TTTAT
   // >read6
   // GCTTA
   // >read7
   // AGGCT
   // >read8
   // CACGGTATT
   // 7 singleton super´nodes, and one 5-node supernode. So N50 is 1
   
   //first set up the hash/graph
   db_graph = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);


   //count_kmers=0;
   bad_reads=0;
   seq_length=0;
   seq_length = load_fasta_from_filename_into_graph("../data/test/graph/n50_example2.fasta",&bad_reads, max_read_length,  db_graph);

   CU_ASSERT_EQUAL(seq_length,46);
   //  CU_ASSERT_EQUAL(count_kmers,12);
   CU_ASSERT_EQUAL(bad_reads,0);


   CU_ASSERT(db_graph_get_N50_of_supernodes(db_graph)==1); //remember this leaves the nodes all visited

   //clean-up
   hash_table_free(&db_graph);
   

   



}


void test_is_condition_true_for_all_nodes_in_supernode()
{
  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits=4;
  int bucket_size   = 10;
  int seq_length;
  long long bad_reads = 0;

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a dB_graph.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****


  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta", &bad_reads, 20,  db_graph);
  
  CU_ASSERT_EQUAL(seq_length,6);
  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),2);
  CU_ASSERT_EQUAL(bad_reads,0);


  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element1!=NULL);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("ACG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element2!=NULL);

  int limit =50;
  dBNode * nodes_path[limit];
  Orientation orientations_path[limit];
  Nucleotide labels_path[limit];
  char seq[limit+1];
  int length_path;
  double avg_coverage;
  boolean is_cycle;
  int min_coverage, max_coverage;

  //  printf("Check if all nodes have status none - this should be true\n");
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,						   
								  &db_node_check_status_none,
								  &db_node_action_do_nothing, 
								  nodes_path, orientations_path, labels_path, &length_path,
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph));
  CU_ASSERT(db_node_check_status(test_element1, none));
  //  printf("set status of one node to visited\n");
  db_node_set_status(test_element1, visited);
  CU_ASSERT(db_node_check_status(test_element1, visited));
  CU_ASSERT(db_node_check_status(test_element2, none));
  // printf("Checkif all nodes have status none - this should not be true - confirm this\n");
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								   &db_node_check_status_none,
								   &db_node_action_do_nothing, 
								   nodes_path, orientations_path, labels_path, &length_path, 
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph));


  //  printf("double check statuses of nodes unchanged\n");
  CU_ASSERT(db_node_check_status(test_element1, visited));
  CU_ASSERT(db_node_check_status(test_element2, none));

  // printf("Check if all nodes have status visited - this should not be true - confirm this\n");
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								   &db_node_check_status_visited,
								   &db_node_action_do_nothing, 
								   nodes_path, orientations_path, labels_path, &length_path,
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph));


  //check again, but this time, set all nodes to visited in the process
  //printf("Set all nodes to visited while checking them all\n");
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50, 
								   &db_node_check_status_visited,
								   &db_node_action_set_status_visited_or_visited_and_exists_in_reference, 
								   nodes_path, orientations_path, labels_path, &length_path, 
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph));
  //and now this time it SHOULD BE TRUE
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50, 
								  &db_node_check_status_visited,
								  &db_node_action_do_nothing, 
								  nodes_path, orientations_path, labels_path, &length_path, 
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph));

  // and should still be true
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								  &db_node_check_status_visited,
								  &db_node_action_do_nothing, 
								  nodes_path, orientations_path, labels_path, &length_path, 
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph));
  

  //printf("Nodes currently all visited. Set all nodes to none while checking them all to confirm that they are currently all visited (before doing the set-to-none)\n");
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								  &db_node_check_status_visited,
								  &db_node_action_set_status_none, 
								  nodes_path, orientations_path, labels_path, &length_path,
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph));

  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								  &db_node_check_status_none,
								  &db_node_action_do_nothing, 
								  nodes_path, orientations_path, labels_path, &length_path, 
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph));
  
  db_node_set_status(test_element2, pruned);
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,
								   &db_node_check_status_none,
								   &db_node_action_do_nothing, 
								   nodes_path, orientations_path, labels_path, &length_path, 
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph));

  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element2, 50,
								   &db_node_check_status_none,
								   &db_node_action_do_nothing, 
								   nodes_path, orientations_path, labels_path, &length_path,
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph));

   
  hash_table_free(&db_graph);
  
  
}



void test_read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference()
{

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits=4;
  int bucket_size   = 10;
  long long bad_reads=0;
  int max_read_length=30;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

  int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person.fasta", &bad_reads, max_read_length, db_graph);
  
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/chrom1.fasta", db_graph);
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/chrom2.fasta", db_graph);
  
  //Now see if it correctly gets the supernode that does not intersect a chromosome
  dBNode * nodes_path[100];
  Orientation orientations_path[100];
  Nucleotide labels_path[100];
  char seq[100];
  int length_path;
  double avg_coverage;
  boolean is_cycle;
  int min_coverage, max_coverage;

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;


  //element on supernode we know intersects chromosomes
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element1!=NULL);
  //elemtn on node that does not intersect chromosomes - is "novel"
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element2!=NULL);


  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50, 
								   &db_node_check_status_is_not_exists_in_reference,
								   &db_node_action_set_status_visited, 
								   nodes_path, orientations_path, labels_path, &length_path, 
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph));

  CU_ASSERT(length_path==3);
  CU_ASSERT_STRING_EQUAL(seq, "TTC");

  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element2, 50,
								  &db_node_check_status_is_not_exists_in_reference,
								  &db_node_action_set_status_visited,
								  nodes_path, orientations_path, labels_path, &length_path, 
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph));
  CU_ASSERT(length_path==2);
  CU_ASSERT_STRING_EQUAL(seq, "CC");

  hash_table_free(&db_graph);



  // now a harder example.
  //first set up the hash/graph
  kmer_size = 31;
  number_of_bits=10;
  bucket_size   = 10;
  bad_reads=0;
  max_read_length=2000;
  int max_rehash_tries=10;

  db_graph = hash_table_new(number_of_bits,bucket_size,max_rehash_tries,kmer_size);

  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person2.fasta", &bad_reads, max_read_length, db_graph);
  
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", db_graph);


  //Now see if it correctly gets the supernode that does not intersect a chromosome

  //element on supernode we know intersects chromosomes
  test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AACCCTAACCCTAACCCTAACCCTAACCCTA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element1!=NULL);
  //elemtn on node that does not intersect chromosomes - is "novel"
  test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("GCGGGGCGGGGCGGGGCGGGGCGGGGCCCCC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element2!=NULL);


  length_path=0;
  CU_ASSERT(!db_graph_is_condition_true_for_all_nodes_in_supernode(test_element1, 50,  
								   &db_node_check_status_is_not_exists_in_reference,
								   &db_node_action_set_status_visited, 
								   nodes_path, orientations_path, labels_path, &length_path, 
								   seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								   db_graph));

  CU_ASSERT(length_path==6);

  CU_ASSERT_STRING_EQUAL(seq,"ACCCTA");
  
  CU_ASSERT(db_graph_is_condition_true_for_all_nodes_in_supernode(test_element2, 50,
								  &db_node_check_status_is_not_exists_in_reference,
								  &db_node_action_set_status_visited, 
								  nodes_path, orientations_path, labels_path, &length_path,
								  seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
								  db_graph));

  CU_ASSERT(length_path==13);
  CU_ASSERT_STRING_EQUAL(seq, "CCCTCACACACAT");
	    

  hash_table_free(&db_graph);

  
  // and now another

  //first set up the hash/graph
  kmer_size = 31;
  number_of_bits=20;
  bucket_size   = 10;
  bad_reads=0;
  max_read_length=2000;
  max_rehash_tries=10;
  seq_length=0;
  db_graph = hash_table_new(number_of_bits,bucket_size,max_rehash_tries,kmer_size);

  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person3.fasta", &bad_reads, max_read_length, db_graph);


  // fasta looks like this:
  /*
    >read1 overlaps human chrom 1 - it is 9 copies of TAACCC and then TAACC on the end.  
    TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC
    > read 2 overlaps human chrom 1 - is has a CCCC in the middle which means any 31-mer will contain CCCC and so not overlap the above
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
  

  // This person has the following supernodes:

  // >this supernode is the entrety of read2 and is entirely contained in chrom1
  // ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAAC
  // >node_1 - This supernode lies entirely in chrom1 also. However note the supernode is a loop, so you could start printing at various places....
  // TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
  // > node 3 - overlaps chrom but has low covg
  // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  // >node 4 - no overlap with chromosome
  // GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCCCCCTCACACACAT
  

  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", db_graph);

  char** array_of_supernodes_for_person3= (char**) calloc(10,sizeof(char*));
  array_of_supernodes_for_person3[0]= (char*)calloc(100,sizeof(char));
  array_of_supernodes_for_person3[1]= (char*)calloc(100,sizeof(char));
  array_of_supernodes_for_person3[2]= (char*)calloc(100,sizeof(char));

  int number_of_supernodes=0;
  int min_covg_required = 2;


  //element on supernode we know intersects chromosomes
   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACCCTAACCCTAACCCTAACCCTAACCCTAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element1!=NULL);
  //elemtn on node that does not intersect chromosomes - is "novel" - and has coverage 3
  test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("GGGCGGGGCGGGGCGGGGCGGGGCCCCCTCA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element2!=NULL);
  //element does not intersect chrom but has covg only 1
  dBNode* test_element3 =  hash_table_find(element_get_key(seq_to_binary_kmer("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element3!=NULL);





  db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(db_graph, &db_node_check_status_is_not_exists_in_reference, min_covg_required, NULL,
									       true, array_of_supernodes_for_person3, &number_of_supernodes);


  CU_ASSERT(number_of_supernodes==1);
  CU_ASSERT( !strcmp(array_of_supernodes_for_person3[0], "CCCGCCCCGCCCC")
	     || !strcmp(array_of_supernodes_for_person3[0], "CCCTCACACACAT"));
  
  
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);


  number_of_supernodes=0;
  //here just a quick check of the _at_least_ version of the function
  min_covg_required=1;
  db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(db_graph,
										       &db_node_check_status_exists_in_reference,  min_covg_required, 
										       NULL, true, array_of_supernodes_for_person3, 
										       &number_of_supernodes);


  CU_ASSERT(number_of_supernodes==2);


  //some of read 1 is in chrom1, and the supernode is printed above, just after we loaded the fasta
  //the whole of read2 is a supernode, and is in chrom 1. Note this print function prints only the edges, not the first kmer in the path


  //one of the supernodes is unambiguous
  CU_ASSERT(    (!strcmp(array_of_supernodes_for_person3[0],"ACCCTAACCCTAAC")) 
		|| (!strcmp(array_of_supernodes_for_person3[1],"ACCCTAACCCTAAC"))
		|| (!strcmp(array_of_supernodes_for_person3[0],"GTTAGGGTTAGGGT")) 
		|| (!strcmp(array_of_supernodes_for_person3[1],"GTTAGGGTTAGGGT")) 
		);
	     

  //the other one is ambiguous because you could start printing at one of 5 places (is a loop of 5 kmers)

  /* - this is what we test for if the print prints the whole supernode including the first 31 bases
  CU_ASSERT( 
            !strcmp(array_of_supernodes_for_person3[0],"TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT") || !strcmp(array_of_supernodes_for_person3[1],"TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA") || !strcmp(array_of_supernodes_for_person3[1],"AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA") || !strcmp(array_of_supernodes_for_person3[1],"AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT") || !strcmp(array_of_supernodes_for_person3[1],"TAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA") || !strcmp(array_of_supernodes_for_person3[1],"ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGT") || !strcmp(array_of_supernodes_for_person3[1],"TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGT")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC") || !strcmp(array_of_supernodes_for_person3[1],"CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG") || !strcmp(array_of_supernodes_for_person3[1],"GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC") || !strcmp(array_of_supernodes_for_person3[1],"CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG") || !strcmp(array_of_supernodes_for_person3[1],"GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC") || !strcmp(array_of_supernodes_for_person3[1],"CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG") || !strcmp(array_of_supernodes_for_person3[1],"GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG")
	    
	     );
  */



 CU_ASSERT( !strcmp(array_of_supernodes_for_person3[0],"AACCCT") || !strcmp(array_of_supernodes_for_person3[1],"AACCCT")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"AGGGTT") || !strcmp(array_of_supernodes_for_person3[1],"AGGGTT")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"ACCCTA") || !strcmp(array_of_supernodes_for_person3[1],"ACCCTA")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"TAGGGT") || !strcmp(array_of_supernodes_for_person3[1],"TAGGGT")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"CCCTAA") || !strcmp(array_of_supernodes_for_person3[1],"CCCTAA")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"TTAGGG") || !strcmp(array_of_supernodes_for_person3[1],"TTAGGG")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"CCTAAC") || !strcmp(array_of_supernodes_for_person3[1],"CCTAAC")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"GTTAGG") || !strcmp(array_of_supernodes_for_person3[1],"GTTAGG")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"CTAACC") || !strcmp(array_of_supernodes_for_person3[1],"CTAACC")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"GGTTAG") || !strcmp(array_of_supernodes_for_person3[1],"GGTTAG")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"TAACCC") || !strcmp(array_of_supernodes_for_person3[1],"TAACCC")
	    ||
	    !strcmp(array_of_supernodes_for_person3[0],"GGGTTA") || !strcmp(array_of_supernodes_for_person3[1],"GGGTTA")

	    
	     );


  free(array_of_supernodes_for_person3[0]) ;
  free(array_of_supernodes_for_person3[1]) ;
  free(array_of_supernodes_for_person3[2]) ;
  free(array_of_supernodes_for_person3) ;

  hash_table_free(&db_graph);

  

  
}





void test_indel_discovery_simple_test_1()
{


  //first set up the hash/graph
  int kmer_size = 31;
  int number_of_bits=10;
  int bucket_size   = 5;
  long long bad_reads=0;
  int max_read_length=600;

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

  int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person_with_sv.fasta", &bad_reads, max_read_length, db_graph);
  CU_ASSERT(seq_length==612);

  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", db_graph);

  BinaryKmer tmp_kmer1;
  BinaryKmer tmp_kmer2;

  
  //element on supernode we know intersects chromosome entirely
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("TGTGCAGAGGACAACGCAGCTCCGCCCTCGC", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element1!=NULL);

  //element on supernode that overlaps at start and end but not middle, with chromosome
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("CCGGCGCAGGCGCAGTTGTTGTAGAGGCGCG", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2),db_graph);
  CU_ASSERT(test_element2!=NULL);


  dBNode * nodes_path[100];
  Orientation orientations_path[100];
  Nucleotide labels_path[100];
  char seq[100];
  int length_path;
  int min_diff=1; //enought o have one node different in the middle
  double avg_coverage;
  int min_coverage, max_coverage;
  boolean is_cycle;


  int num_nodes_we_demand_overlap_with_reference_at_start=2;
  int num_nodes_we_demand_overlap_with_reference_at_end=22;
  
  CU_ASSERT(!db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(test_element1, 200,
											 &db_node_check_status_exists_in_reference,
											 &db_node_action_set_status_visited_or_visited_and_exists_in_reference,
											 num_nodes_we_demand_overlap_with_reference_at_start,
											 num_nodes_we_demand_overlap_with_reference_at_end, min_diff, 
											 nodes_path, orientations_path, labels_path, &length_path, 
											 seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
											 db_graph));
  
  CU_ASSERT(db_graph_is_condition_true_for_start_and_end_but_not_all_nodes_in_supernode(test_element2, 200, 
											&db_node_check_status_exists_in_reference,
											
											&db_node_action_set_status_visited_or_visited_and_exists_in_reference,
											num_nodes_we_demand_overlap_with_reference_at_start,
											num_nodes_we_demand_overlap_with_reference_at_end, min_diff, 
											nodes_path, orientations_path, labels_path, &length_path, 
											seq,&avg_coverage,&min_coverage,&max_coverage,&is_cycle,
											db_graph));
  

  //remove the visited markings
  hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);
  
  
  char** array_of_supernodes= (char**) calloc(10,sizeof(char*));
  array_of_supernodes[0]= (char*)calloc(100,sizeof(char));
  array_of_supernodes[1]= (char*)calloc(100,sizeof(char));
  array_of_supernodes[2]= (char*)calloc(100,sizeof(char));

  int number_of_supernodes=0;
  int min_covg_required = 2;
  int min_start = 2;
  int min_end = 2;
  
  db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
												    min_start, min_end, min_diff, NULL,
												    true, array_of_supernodes, &number_of_supernodes);
												    



  printf("We get %d supernodes\n\n\n\n", number_of_supernodes);
  CU_ASSERT(number_of_supernodes==1);

  CU_ASSERT( !strcmp(array_of_supernodes[0], "AGTTGTTGTAGAGGCGCGCCGCGCCGGCGCAGGCGCAGACACATGCTAGCGCGTCGGGGTGGAGGCGT") || !strcmp(array_of_supernodes[0], "TGCGCCTGCGCCGGCGCGGCGCGCCTCTACAACAACTGCGCCTGCGCCGGCGGAGTTGCGTTCTCCTC"));
    //  CU_ASSERT( !strcmp(array_of_supernodes[0], "GAGGAGAACGCAACTCCGCCGGCGCAGGCGCAGTTGTTGTAGAGGCGCGCCGCGCCGGCGCAGGCGCAGACACATGCTAGCGCGTCGGGGTGGAGGCGT") || !strcmp(array_of_supernodes[0], "ACGCCTCCACCCCGACGCGCTAGCATGTGTCTGCGCCTGCGCCGGCGCGGCGCGCCTCTACAACAACTGCGCCTGCGCCGGCGGAGTTGCGTTCTCCTC"));


  free(array_of_supernodes[0]) ;
  free(array_of_supernodes[1]) ;
  free(array_of_supernodes[2]) ;
  free(array_of_supernodes) ;
  
  

  hash_table_free(&db_graph);
}




//suppose we are given a pair of start-end coords within which we think there is a deletion.
// Load your fasta, then load precisely that section of the chromosome as reference, and print supernodes that match ref at start and end
// (WARNING: if you delete the middle of a sequence, you create new kmers in the middle that might not have been there before.
//           so the nodes in the new supernode are NOT necessatrily all in the reference)
// Then build a new graph just out of the section of the chromosome, and load the supernodes above as reference. Now look for supernodes
// in the chromosome graph that match the "reference" supernodes at the start and end only. These are potential deletions.
void test_deletion_validation()
{

  //first set up the hash/graph
  int kmer_size = 31;
  int number_of_bits=10;
  int bucket_size   = 5;
  long long bad_reads=0;
  int max_read_length=2000;



  //STEP 1: get supernodes from our person which match reference at start and end.
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  int seq_length = load_fasta_from_filename_into_graph("../data/test/graph/person_with_deletion_in_chrom.fasta",  &bad_reads, max_read_length, db_graph);
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", db_graph);


  FILE* intermediate_output = fopen("../data/test/graph/test_db_graph_intermediate_output_file", "w");
  int min_covg_required = 1;
  int min_start = 1;
  int min_end = 31; //iei 31 bases at start and end
  int min_diff = 2;
  db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
                                                                                                    intermediate_output, false, NULL, 0);

  //  db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
  //                                                                                                  min_start, min_end, min_diff, intermediate_output,
  //												    false, NULL, 0);
  fclose(intermediate_output);
  hash_table_free(&db_graph);

  //STEP 2: Load the chromosome as a person, and the previosu supernodes as reference
  db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  seq_length = load_fasta_from_filename_into_graph("../data/test/graph/Homo_sapiens.NCBI36.52.dna.chromosome.1.first_20_lines.fasta", &bad_reads, max_read_length,db_graph);
  read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("../data/test/graph/test_d_graph_intermediate_output_file", db_graph);



  char** array_of_supernodes= (char**) calloc(10,sizeof(char*));
  array_of_supernodes[0]= (char*)calloc(100,sizeof(char));
  array_of_supernodes[1]= (char*)calloc(100,sizeof(char));
  array_of_supernodes[2]= (char*)calloc(100,sizeof(char));

  int number_of_supernodes=0;

  //STEP 3: Print supernodes in chromosome that match our person's supernodes (here the reference) at start and end - deletions
   db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
                                                                                                    min_start, min_end, min_diff, NULL,
                                                                                                    true, array_of_supernodes, &number_of_supernodes);
  //db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(db_graph, &db_node_check_status_exists_in_reference, min_covg_required,
   //										       NULL,true, array_of_supernodes, &number_of_supernodes);


   
  printf("deletion %s and number fo supernodes is %d\n", array_of_supernodes[0], number_of_supernodes);
  CU_ASSERT_STRING_EQUAL("GAGAGGCGCGCCGCGCCGGCGC", array_of_supernodes[0]);

  free(array_of_supernodes[0]);
  free(array_of_supernodes[1]);
  free(array_of_supernodes[2]);
  free(array_of_supernodes);
  hash_table_free(&db_graph);

}

