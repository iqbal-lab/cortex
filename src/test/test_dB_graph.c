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
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0; 
  


  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  
  //Load the following fasta:
  //    >read1
  //    AAAAAAAA
  //    >read2
  //    GGCT
  //    >read3
  //    TAGG

  int seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/test_dB_graph.fasta",&bad_reads, 20, db_graph);
  
  //length of total sequence
  CU_ASSERT(seq_length == 16);
  
  //number of kmers

  CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph), 5);

  //bad reads
  CU_ASSERT_EQUAL(bad_reads, 0);

  //all the kmers and their reverse complements from the reads
  dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("AAA", kmer_size),kmer_size) ,db_graph);
  dBNode* test_element2 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element3 = hash_table_find(element_get_key(seq_to_binary_kmer("GGC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element5 = hash_table_find(element_get_key(seq_to_binary_kmer("GCT",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element6 = hash_table_find(element_get_key(seq_to_binary_kmer("AGC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element7 = hash_table_find(element_get_key(seq_to_binary_kmer("TAG",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element8 = hash_table_find(element_get_key(seq_to_binary_kmer("CTA",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element9 = hash_table_find(element_get_key(seq_to_binary_kmer("AGG",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element10 = hash_table_find(element_get_key(seq_to_binary_kmer("CCT",  kmer_size),kmer_size) ,db_graph);

  //kmers that should not be in the graph

  dBNode* test_element11 = hash_table_find(element_get_key(seq_to_binary_kmer("GGG", kmer_size),kmer_size) ,db_graph);
  dBNode* test_element12 = hash_table_find(element_get_key(seq_to_binary_kmer("CCC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element13 = hash_table_find(element_get_key(seq_to_binary_kmer("TAT",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element14 = hash_table_find(element_get_key(seq_to_binary_kmer("ATA",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element15 = hash_table_find(element_get_key(seq_to_binary_kmer("TAC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element16 = hash_table_find(element_get_key(seq_to_binary_kmer("ATG",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element17 = hash_table_find(element_get_key(seq_to_binary_kmer("TTG",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element18 = hash_table_find(element_get_key(seq_to_binary_kmer("AAC",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element19 = hash_table_find(element_get_key(seq_to_binary_kmer("TGA",  kmer_size),kmer_size) ,db_graph);
  dBNode* test_element20 = hash_table_find(element_get_key(seq_to_binary_kmer("TCA",  kmer_size),kmer_size) ,db_graph);


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



void test_tip_clipping()
{

  //first set up the hash/graph
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size    = 4;
  long long bad_reads = 0; 

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  

  //Load the following fasta:
  
  //>main_trunk
  //GCGTCCCAT
  //>tip
  //CGTTT

  int seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_tip.fasta",&bad_reads, 20, db_graph);


 

  CU_ASSERT_EQUAL(seq_length,14);

  dBNode* node1 = hash_table_find(element_get_key(seq_to_binary_kmer("TTT", kmer_size),kmer_size) ,db_graph);

  CU_ASSERT(node1 != NULL);
  
  int tip_length = db_graph_db_node_clip_tip(node1, 10 , &db_node_check_status_none,&db_node_action_set_status_pruned,db_graph);
  
  CU_ASSERT_EQUAL(tip_length,2);

  dBNode* node2 = hash_table_find(element_get_key(seq_to_binary_kmer("GTT", kmer_size),kmer_size) ,db_graph);
  dBNode* node3 = hash_table_find(element_get_key(seq_to_binary_kmer("CGT", kmer_size),kmer_size) ,db_graph);



  CU_ASSERT(db_node_check_status(node1,pruned));
  CU_ASSERT(db_node_check_status(node2,pruned));
  
  CU_ASSERT(db_node_check_status(node3,none));


  dBNode * nodes_path[100];
  Orientation orientations_path[100];
  Nucleotide labels_path[100];
  char tmp_seq[100+db_graph->kmer_size+1];


  int length_supernode = db_graph_supernode(node3,100,&db_node_check_status_none,&db_node_action_set_status_visited,tmp_seq,nodes_path,orientations_path,labels_path,db_graph);

  CU_ASSERT_EQUAL(length_supernode,6);

  
  CU_ASSERT_STRING_EQUAL("ATGGGACGC",tmp_seq);


  //check ends
  dBNode* node4 = hash_table_find(element_get_key(seq_to_binary_kmer("GCG", kmer_size),kmer_size) ,db_graph);
  dBNode* node5 = hash_table_find(element_get_key(seq_to_binary_kmer("CAT", kmer_size),kmer_size) ,db_graph);

  
  CU_ASSERT(db_node_check_status(node3,visited));
  CU_ASSERT(db_node_check_status(node4,visited));
  CU_ASSERT(db_node_check_status(node5,visited));
   

  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);
  

}
  

void test_supernode_walking() //test db_graph_get_perfect_path 
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

  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);

 
  //1. Sequence of tests as follows
  //         Each test loads a single specifically designed fasta file into a dB_graph.
  //         The test then picks an element in the graph, and calls get_seq_from_elem_to_end_of_supernode
  //         and checks that it gets the right sequence.
  

  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****


  seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_two_self_loops.fasta", &bad_reads, 20,  db_graph);

   CU_ASSERT_EQUAL(seq_length,6);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),2);
   CU_ASSERT_EQUAL(bad_reads,0);
   

  // now start at GTA and get all the sequence from there to the end of the supernode, and see
  // if that is right. Remember that once encoded GTA, might be in a kmer where to read GTA you must go
  // in the reverse orientation.

   dBNode* test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA", kmer_size), kmer_size),db_graph);
   CU_ASSERT(test_element1!=NULL);
   


   boolean is_cycle=false;
   
   int test1_length = db_graph_get_perfect_path(test_element1,forward,100,&db_node_action_do_nothing,path_nodes,path_orientations,path_labels,&is_cycle,db_graph);

    CU_ASSERT(is_cycle);


   //depending on what orientation GTA is with respect to it's kmer, the answer you get will be either CGT or GT
   //Zam we know GTA<TAC!
   
   CU_ASSERT_STRING_EQUAL(nucleotides_to_string(path_labels,test1_length,tmp_seq),"CGTA");
   

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
   seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_conflict_at_end.fasta",&bad_reads,20,db_graph);

   CU_ASSERT_EQUAL(seq_length,13);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),5);
   CU_ASSERT_EQUAL(bad_reads,0);


   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size), kmer_size),db_graph);
   CU_ASSERT(test_element1!=NULL);

   //ACA < TGT so forward gives TT
   is_cycle=false;
   
   int test2_length = db_graph_get_perfect_path(test_element1,forward,100,&db_node_action_do_nothing,path_nodes,path_orientations,path_labels,&is_cycle,db_graph);

   CU_ASSERT(!is_cycle);
   
   CU_ASSERT_STRING_EQUAL(nucleotides_to_string(path_labels,test2_length,tmp_seq),"TT");

   //ACA < TGT so backward gives ""

   test2_length = db_graph_get_perfect_path(test_element1,reverse,100,db_node_action_do_nothing,path_nodes,path_orientations,path_labels,&is_cycle,db_graph);
  
   CU_ASSERT(!is_cycle);
   CU_ASSERT_EQUAL(test2_length,0);
   CU_ASSERT_STRING_EQUAL(nucleotides_to_string(path_labels,test2_length,tmp_seq),"");


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
   
   seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_one_long_supernode_with_inward_conflict_at_end.fasta", &bad_reads, 20,  db_graph);

   CU_ASSERT_EQUAL(seq_length,13);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),5);
   CU_ASSERT_EQUAL(bad_reads,0);
   

   test_element1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA", kmer_size), kmer_size) ,db_graph);
   CU_ASSERT(test_element1!=NULL);

   //ACA < TGT so forward gives TT
   is_cycle=false;

   int test3_length = db_graph_get_perfect_path(test_element1,forward,100,db_node_action_do_nothing,path_nodes,path_orientations,path_labels,&is_cycle,db_graph);

   CU_ASSERT(!is_cycle);
   CU_ASSERT_STRING_EQUAL(nucleotides_to_string(path_labels,test3_length,tmp_seq),"TT");

   
   //ACA < TGT so backwar gives ""
   is_cycle=false;
   test3_length = db_graph_get_perfect_path(test_element1,reverse,100,db_node_action_do_nothing,path_nodes,path_orientations,path_labels,&is_cycle,db_graph);

   CU_ASSERT(!is_cycle);
   CU_ASSERT_STRING_EQUAL(nucleotides_to_string(path_labels,test3_length,tmp_seq),"");

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
   seq_length = load_fasta_data_from_filename_into_graph("../data/test/graph/generates_graph_with_infinite_loop.fasta",&bad_reads,30,db_graph);

   CU_ASSERT_EQUAL(seq_length,25);
   CU_ASSERT_EQUAL(hash_table_get_unique_kmers(db_graph),1);
   CU_ASSERT_EQUAL(bad_reads,0);

   test_element1 = hash_table_find(seq_to_binary_kmer("AAA", kmer_size) ,db_graph);
   CU_ASSERT(test_element1!=NULL);

   //forward
   is_cycle=false;
   int test4_length = db_graph_get_perfect_path(test_element1,forward,100,db_node_action_do_nothing,path_nodes,path_orientations,path_labels,&is_cycle,db_graph);
   
   CU_ASSERT(is_cycle);//this time it should find a cycle
   CU_ASSERT_STRING_EQUAL(nucleotides_to_string(path_labels,test4_length,tmp_seq),"A");

   //backward
   is_cycle=false;
   test4_length = db_graph_get_perfect_path(test_element1,reverse,100,db_node_action_do_nothing,path_nodes,path_orientations,path_labels,&is_cycle,db_graph);

   CU_ASSERT(is_cycle);//this time it should find a cycle
   CU_ASSERT_STRING_EQUAL(nucleotides_to_string(path_labels,test4_length,tmp_seq),"T");

   hash_table_free(&db_graph);
   CU_ASSERT(db_graph == NULL);


}


void test_writing_reading_graph(){

  dBNode node1, node2;
  int kmer_size = 3;
  boolean test_read;

  char seq[kmer_size];

  element_initialise(&node1,seq_to_binary_kmer("AAA",kmer_size),kmer_size);
  node1.edges = 'a'; // ie 64 = (0110 0100)2
  node1.coverage = 10;

  FILE* fp1 = fopen("../data/test/graph/dump_element.bin", "w");
  db_node_print_binary(fp1,&node1);
  fclose(fp1);

  FILE* fp2 = fopen("../data/test/graph/dump_element.bin", "r");
  
  test_read = db_node_read_binary(fp2,kmer_size,&node2);

  CU_ASSERT_EQUAL(test_read, true);

  CU_ASSERT_EQUAL(node1.kmer, node2.kmer);
  CU_ASSERT_EQUAL(node1.edges, node2.edges);

  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(node1.kmer,kmer_size,seq));
  CU_ASSERT_STRING_EQUAL("AAA", binary_kmer_to_seq(node2.kmer,kmer_size,seq));
  
  CU_ASSERT_EQUAL(node1.edges, node2.edges);
  CU_ASSERT_EQUAL(node1.edges, 'a');
  CU_ASSERT_EQUAL(node2.edges, 'a');
 
  CU_ASSERT_EQUAL(node2.coverage, 10);

  test_read = db_node_read_binary(fp2,kmer_size,&node2);
  CU_ASSERT_EQUAL(test_read, false);

  fclose(fp2);
}
   
void test_get_perfect_path(){
  int kmer_size = 3;
  int number_of_bits = 4;
  int bucket_size = 4;
  dBGraph * db_graph = hash_table_new(number_of_bits,bucket_size,10,kmer_size);
  boolean found1, found2, found3;
  dBNode * node1;
  dBNode * node2;
  dBNode * node3;
  dBNode * node4;
 
 
  node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CGT", kmer_size),kmer_size), &found1, db_graph);
  node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("GTT", kmer_size),kmer_size),  &found2, db_graph);
  db_node_add_edge(node1, node2, reverse,reverse, db_graph->kmer_size);	  

  node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size),kmer_size),  &found3, db_graph);
  db_node_add_edge(node2, node3, reverse,reverse, db_graph->kmer_size);	 
 
  node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TAG", kmer_size),kmer_size),  &found3, db_graph);
  db_node_add_edge(node3, node4, reverse,reverse, db_graph->kmer_size);	  
  
  dBNode * nodes[10];
  Orientation orientations[10];
  Nucleotide bases[10];
  boolean is_cycle;
  int length = 0;

  length = db_graph_get_perfect_path(node1,reverse, 10,db_node_action_set_status_visited,nodes, orientations, bases, &is_cycle,db_graph);

  CU_ASSERT_EQUAL(length,3);
  

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

  hash_table_free(&db_graph);
  CU_ASSERT(db_graph == NULL);

}

void test_get_perfect_bubble(){
  int kmer_size = 3;
  int number_of_buckets = 4;
  int bucket_size = 3;
  dBGraph * db_graph = hash_table_new(number_of_buckets,bucket_size,10,kmer_size);
  boolean found;
 
  dBNode * node1, * node2, * node3, * node4, * node5, * node6, * node7, * node8, * node9;
  Orientation orientation;
  Nucleotide labels[kmer_size], base1, base2;
  dBNode * end_node;
  Orientation end_orientation;

  boolean bubble;

  //start point
  node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size),kmer_size), &found, db_graph);

  //branch1
  node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size),kmer_size),  &found, db_graph);
  node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CAC", kmer_size),kmer_size),  &found, db_graph);
  node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("ACT", kmer_size),kmer_size),  &found, db_graph);

  //end point
  node5 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTT", kmer_size),kmer_size),  &found, db_graph);

  //branch 2
  node6 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size),kmer_size),  &found, db_graph);
  node7 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CTC", kmer_size),kmer_size),  &found, db_graph);
  node8 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TCT", kmer_size),kmer_size),  &found, db_graph);

  //add 3p extension
  node9 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("TTA", kmer_size),kmer_size),  &found, db_graph);
 
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

  //add 3p extension
  db_node_add_edge(node5, node9, reverse,reverse, db_graph->kmer_size);


  bubble = db_graph_detect_perfect_bubble(node1,&orientation,&db_node_check_status_none,&db_node_action_set_status_visited,
&base1,&base2,labels,&end_node,&end_orientation,db_graph);

  CU_ASSERT(bubble);
  CU_ASSERT_EQUAL(orientation,forward);
  CU_ASSERT(base1 == Adenine || base1 == Thymine);
  CU_ASSERT(base2 == Adenine || base2 == Thymine);

  CU_ASSERT_EQUAL(labels[0],Cytosine);
  CU_ASSERT_EQUAL(labels[1],Thymine);
  CU_ASSERT_EQUAL(labels[2],Thymine);

  CU_ASSERT_EQUAL(end_node,node5);
  CU_ASSERT_EQUAL(end_orientation,reverse);
  
  Orientation orientations5p[100];
  Orientation orientations3p[100];
  Nucleotide labels_flank5p[100]; 
  Nucleotide labels_flank3p[100];
  boolean is_cycle5p, is_cycle3p;
  int length_flank5p,length_flank3p;
  dBNode * nodes5p[100];
  dBNode * nodes3p[100];

  length_flank5p = db_graph_get_perfect_path(node1, opposite_orientation(orientation),100,&db_node_action_set_status_visited,nodes5p,orientations5p, labels_flank5p,
  					     &is_cycle5p,db_graph); 

  CU_ASSERT_EQUAL(length_flank5p,0);
  CU_ASSERT_EQUAL(nodes5p[0],node1);

  length_flank3p = db_graph_get_perfect_path(end_node, end_orientation,100,&db_node_action_set_status_visited,nodes3p,orientations3p, labels_flank3p,
					     &is_cycle3p,db_graph);    
 
  
  CU_ASSERT_EQUAL(length_flank3p,1);
  CU_ASSERT_EQUAL(nodes3p[0],end_node);
  CU_ASSERT_EQUAL(nodes3p[1],node9);  
  CU_ASSERT_EQUAL(labels_flank3p[0],Adenine);

  //if searched again should return false -- becuase nodes are marked as visited
  bubble = db_graph_detect_perfect_bubble(node1,&orientation,&db_node_check_status_none,&db_node_action_set_status_visited,&base1,&base2,labels,&end_node,&end_orientation,db_graph);
  CU_ASSERT(!bubble);

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

  node1 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCC", kmer_size),kmer_size), &found, db_graph);
  node2 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCT", kmer_size),kmer_size), &found, db_graph);
  node3 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCG", kmer_size),kmer_size), &found, db_graph);

 

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
  

  node4 = hash_table_find_or_insert(element_get_key(seq_to_binary_kmer("CCA", kmer_size),kmer_size), &found, db_graph);
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
