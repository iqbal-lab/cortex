#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <dB_graph_population.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <stdlib.h>
#include "supernode_cmp.h"

void test_load_two_people_in_same_populations_and_print_separately_their_supernodes()
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

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_individuals_simple.txt", &bad_reads, hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 44);
  CU_ASSERT(bad_reads==0);

  char** array_of_supernodes_for_person1= (char**) calloc(10,sizeof(char*));
  char** array_of_supernodes_for_person2= (char**) calloc(10,sizeof(char*));

  if ( (array_of_supernodes_for_person1==NULL) || (array_of_supernodes_for_person2==NULL))
    {
      printf("cant start - OOM");
      exit(1);
    }

  //these counters are used to make sure none of the files of printed out supernodes get too big. In fact
  //they are completely unused in test code, as we don't print anything
  long supernode_count_person1=0;
  long supernode_count_person2=0;
  
  //this on the other hand is used in testing.
  int number_of_supernodes_in_person_1=0;

  //print_supernode will, in debug mode, alloc memory for you in your array, and put the supernode in it


  int i;


  //USE NEW PRINT FUNCTION and then uncomment the below 
  /*

  //  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person1, individual_edge_array, 0, 
  //								  true, array_of_supernodes_for_person1,&number_of_supernodes_in_person_1);

  //printf("PERSON 1 has %d supernodes\n", number_of_supernodes_in_person_1);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  int number_of_supernodes_in_person_2=0;
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person2, individual_edge_array, 1, 
  					   true, array_of_supernodes_for_person2,&number_of_supernodes_in_person_2);
  //printf("PERSON 2 has %d supernodes\n", number_of_supernodes_in_person_2);

  CU_ASSERT(number_of_supernodes_in_person_1==2);
  CU_ASSERT(number_of_supernodes_in_person_2==2);


  //for (i=0; i<number_of_supernodes_in_person_1; i++)
  // {
  //   printf("\nPerson 1 node %d : %s\n",i,array_of_supernodes_for_person1[i]);
  // }
  //for (i=0; i<number_of_supernodes_in_person_2; i++)
  // {
  //   printf("\nPerson 2 node %d : %s\n",i,array_of_supernodes_for_person2[i]);
  // }


  //Quicksort these results:
  qsort((void*) array_of_supernodes_for_person1, number_of_supernodes_in_person_1 , sizeof(char*),  &supernode_cmp);
  qsort((void*) array_of_supernodes_for_person2, number_of_supernodes_in_person_2 , sizeof(char*),  &supernode_cmp);
  // for (i=0; i<number_of_supernodes_in_person_1; i++)
  // {
  //   printf("\nAfter sort Person 1 node %d : %s\n",i,array_of_supernodes_for_person1[i]);
  // }
  //for (i=0; i<number_of_supernodes_in_person_2; i++)
  // {
  //   printf("\nafter sort Person 2 node %d : %s\n",i,array_of_supernodes_for_person2[i]);
  // }


  //Expected results are

  char* correct_answer_person_1[] ={"AAAA", "AACGTT"};
  char* correct_answer_person_2[] ={"CCCC", "CGTCAA"};


  int j;

  CU_ASSERT(2==number_of_supernodes_in_person_1);
  for (j=0; j<number_of_supernodes_in_person_1; j++)
    {
  //   printf("\nj is %d, Person 1 should have %s and we see %s\n", j,correct_answer_person_1[j], array_of_supernodes_for_person1[j]);
      CU_ASSERT_STRING_EQUAL(correct_answer_person_1[j], array_of_supernodes_for_person1[j]);
    }


  
  CU_ASSERT(2==number_of_supernodes_in_person_2);
  for (j=0; j<number_of_supernodes_in_person_2; j++)
    {
      // printf("Person 2 should have %s and we see%s\n", correct_answer_person_2[j], array_of_supernodes_for_person2[j]);

     CU_ASSERT_STRING_EQUAL(correct_answer_person_2[j], array_of_supernodes_for_person2[j]);
    }




  */


  //cleanup

  for (i=0; i<2; i++)
    {
      free(array_of_supernodes_for_person1[i]);
    }
 for (i=0; i<2; i++)
    {
      free(array_of_supernodes_for_person2[i]);
    }

 free(array_of_supernodes_for_person1);
 free(array_of_supernodes_for_person2);


  hash_table_free(&hash_table);
}


// Three people, with slight variants at one of two loci
// This is a good test set for seeing if we find the right shared variants
// However for this test case, just check that we get the right supernodes for each person,
// and that we can find the subsets of the supernodes that they ALL share

void test_take_three_people_each_with_one_read_and_find_variants()
{

  // FIX TO USE NEW PRINT

  /*


  //first set up the hash/graph
  int kmer_size = 5;
  int number_of_bits = 5;
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

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/three_indiv_simple/three_individuals_simple.txt", &bad_reads, hash_table);

  //printf("take 3 people test Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 55);
  CU_ASSERT(bad_reads==0);

  char** array_of_supernodes_for_person1= (char**) calloc(20,sizeof(char*));
  char** array_of_supernodes_for_person2= (char**) calloc(20,sizeof(char*));
  char** array_of_supernodes_for_person3= (char**) calloc(20,sizeof(char*));
  
  if ( (array_of_supernodes_for_person1==NULL) || (array_of_supernodes_for_person2==NULL) || (array_of_supernodes_for_person3==NULL) )
    {
      printf("cant start - OOM");
      exit(1);
    }

  
  //these counters are used to make sure none of the files of printed out supernodes get too big. In fact
  //they are completely unused in test code, as we don't print anything
  long supernode_count_person1=0;
  long supernode_count_person2=0;
  long supernode_count_person3=0;


  int number_of_supernodes_in_person_1=0;
  //print_supernode will, in debug mode, alloc memory for you in your array, and put the supernode in it
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person1, individual_edge_array, 0, 
					   true, array_of_supernodes_for_person1,&number_of_supernodes_in_person_1);
  //printf("PERSON 1 has %d supernodes\n", number_of_supernodes_in_person_1);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  int i;
  //    for (i=0; i<number_of_supernodes_in_person_1; i++)
  //{
  //  printf("SUPERNODE %s\n", array_of_supernodes_for_person1[i]);
  // }

  
  int number_of_supernodes_in_person_2=0;
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table,&supernode_count_person2, individual_edge_array, 1, 
					   true, array_of_supernodes_for_person2,&number_of_supernodes_in_person_2);
  // printf("PERSON 2 has %d supernodes\n", number_of_supernodes_in_person_2);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  // for (i=0; i<number_of_supernodes_in_person_2; i++)
  //{
  //  printf("SUPERNODE %s\n", array_of_supernodes_for_person2[i]);
  //  }


  int number_of_supernodes_in_person_3=0;
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person3, individual_edge_array, 2, 
					   true, array_of_supernodes_for_person3,&number_of_supernodes_in_person_3);
  printf("PERSON 3 has %d supernodes\n", number_of_supernodes_in_person_3);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  //for (i=0; i<number_of_supernodes_in_person_3; i++)
  // {
  //  printf("SUPERNODE %s\n", array_of_supernodes_for_person3[i]);
  // }


  //qsort results
  qsort((void*) array_of_supernodes_for_person1, number_of_supernodes_in_person_1 , sizeof(char*),  &supernode_cmp);
  qsort((void*) array_of_supernodes_for_person2, number_of_supernodes_in_person_2 , sizeof(char*),  &supernode_cmp);
  qsort((void*) array_of_supernodes_for_person3, number_of_supernodes_in_person_3 , sizeof(char*),  &supernode_cmp);


  //Expected results are

  char* correct_answer_person_1[] ={"AAGCCTCGACAGCCATGC"};
  char* correct_answer_person_2[]={"AAGCCTCGTTCGGCCATGC"};
  char* correct_answer_person_3[]={"AAGCCTCGCTAG","GCATGGCTAG"};
  char* rev_correct_answer_person_1[]={"GCATGGCTGTCGAGGCTT"};
  char* rev_correct_answer_person_2[]={"GCATGGCCGAACGAGGCTT"};
  char* rev_correct_answer_person_3[]={"CTAGCGAGGCTT","CTAGCCATGC"};


  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      //  printf("\ni is %d, person 1, compare %s and %s\n", i, correct_answer_person_1[i], array_of_supernodes_for_person1[i]);
      CU_ASSERT(!strcmp(correct_answer_person_1[i], array_of_supernodes_for_person1[i]) || !strcmp(rev_correct_answer_person_1[i], array_of_supernodes_for_person1[i]) );
    }
  for (i=0; i<number_of_supernodes_in_person_2; i++)
    {
      // printf("\n i is %d, person 2, compare %s and %s\n", i, correct_answer_person_2[i], array_of_supernodes_for_person2[i]);
      CU_ASSERT(!strcmp(correct_answer_person_2[i], array_of_supernodes_for_person2[i]) || !strcmp(rev_correct_answer_person_2[i], array_of_supernodes_for_person2[i]) );
    }
  for (i=0; i<number_of_supernodes_in_person_3; i++)
    {
      //printf("\n i is %d, person 3, compare %s and %s\n", i, correct_answer_person_3[i], array_of_supernodes_for_person3[i]);
      CU_ASSERT(!strcmp(correct_answer_person_3[i], array_of_supernodes_for_person3[i]) || !strcmp(rev_correct_answer_person_3[i], array_of_supernodes_for_person3[i]) );
    }



  //cleanup


  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      free(array_of_supernodes_for_person1[i]);
    }
  for (i=0; i<number_of_supernodes_in_person_2; i++)
    {
      free(array_of_supernodes_for_person2[i]);
    }
  for (i=0; i<number_of_supernodes_in_person_3; i++)
    {
      free(array_of_supernodes_for_person3[i]);
    }

  free(array_of_supernodes_for_person1);
  free(array_of_supernodes_for_person2);
  free(array_of_supernodes_for_person3);


  hash_table_free(&hash_table);
  */

}



void test_take_two_people_sharing_an_alu_and_find_supernodes()
{

  // FIX to use new print

  /*

  //first set up the hash/graph
  int kmer_size = 31;
  int number_of_bits = 20;
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

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/two_people.txt",  &bad_reads,hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  CU_ASSERT(seq_loaded == 677);
  CU_ASSERT(bad_reads ==0);

  char** array_of_supernodes_for_person1= (char**) calloc(316,sizeof(char*));
  char** array_of_supernodes_for_person2= (char**) calloc(255,sizeof(char*));
  
  if ( (array_of_supernodes_for_person1==NULL) || (array_of_supernodes_for_person2==NULL) )
    {
      printf("cant start - OOM");
      exit(1);
    }

  
  //dummy counter used for making sure don't print too many supernodes to a file
  long supernode_count_person1=0;
  long supernode_count_person2=0;
  
  int  number_of_supernodes_in_person_1=0;
  //print_supernode will, in debug mode, alloc memory for you in your array, and put the supernode in it
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person1, individual_edge_array, 0, 
					   true, array_of_supernodes_for_person1,&number_of_supernodes_in_person_1);
  printf("PERSON 1 has %d supernodes\n", number_of_supernodes_in_person_1);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  int i;
  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      printf("SUPERNODE %s\n", array_of_supernodes_for_person1[i]);
    }

  
  int number_of_supernodes_in_person_2=0;
  db_graph_traverse_specific_person_or_pop_for_supernode_printing(&db_graph_choose_output_filename_and_print_supernode_for_specific_person_or_pop, hash_table, &supernode_count_person2, individual_edge_array, 1, 
					   true, array_of_supernodes_for_person2,&number_of_supernodes_in_person_2);
 printf("PERSON 2 has %d supernodes\n", number_of_supernodes_in_person_2);
 db_graph_set_all_visited_nodes_to_status_none(hash_table);
 
 for (i=0; i<number_of_supernodes_in_person_2; i++)
   {
     printf("SUPERNODE %s\n", array_of_supernodes_for_person2[i]);
   }


   printf("\n******   TODO ********Check these supernodes are correct\n");

   

   //Now see which bits of the two supernodes the two individuals have in common



  //cleanup


  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      free(array_of_supernodes_for_person1[i]);
    }
  for (i=0; i<number_of_supernodes_in_person_2; i++)
    {
      free(array_of_supernodes_for_person2[i]);
    }


  free(array_of_supernodes_for_person1);
  free(array_of_supernodes_for_person2);

  hash_table_free(&hash_table);

  */

}


/*   NEED TO UPDATE THESE TESTS TO USE UP-TO-DATE FRAMEWORK - CURRENTLY USING THE CHROMOSOME OVERLAP FRAMEWORK, WHICH IS OBSOLETE
     BUT THE TEST CASES ARE USEFUL ..

*/







// If we want to trust our SV call, then we want to be sure that this section of the genome does not atch multiple chromosomes
void test_db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome()
{

  printf("Update this test to use new framework. Currentluy null test. Checks for section that matches two different chromosomes\n");

  /*

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

  
  // ****
  //1.1 Fasta file that generate a graph with two hairpins, and a single edge (in each rorientation) joining them.
  //  Sequence is :  ACGTAC
  // ****


  int seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_two_self_loops", &bad_reads, hash_table);
  CU_ASSERT(seq_loaded==6);
  CU_ASSERT(bad_reads==0);


  //GTA is not the end of a supernode in either direction - will overlap with chrom 2
  dBNode* query_node1 = hash_table_find(element_get_key(seq_to_binary_kmer("GTA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node1==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node1, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node1, reverse, individual_edge_array, 0, hash_table));

  //CGT is not the end of a supernode in either direction - will overlap with chrom1
  dBNode* query_node2 = hash_table_find(element_get_key(seq_to_binary_kmer("CGT",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node2==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node2, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(!db_node_is_supernode_end(query_node2, reverse, individual_edge_array, 0, hash_table));
  
  //now load a few dummy chromosomes, each of which intersects precisely one different kmer in the graph, so that no node in the graph overlaps >1 chromosome
  //need to take care because we have hairpins

  load_chromosome_overlap_data("../data/test/pop_graph/dummy_chromosomes/simple2/chrom1.fasta", hash_table, 1);
  load_chromosome_overlap_data("../data/test/pop_graph/dummy_chromosomes/simple2/chrom2.fasta", hash_table, 2);

  int number_of_chromosomes_overlapped=-100;
  CU_ASSERT(db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome(query_node1, individual_edge_array, 0, hash_table, &number_of_chromosomes_overlapped));
  //printf("Number of overlapped chroms is %d\n", number_of_chromosomes_overlapped);
  CU_ASSERT(number_of_chromosomes_overlapped==2);
  CU_ASSERT(db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome(query_node2, individual_edge_array, 0, hash_table, &number_of_chromosomes_overlapped))
    //printf("Number of overlapped chroms is %d\n", number_of_chromosomes_overlapped);
  CU_ASSERT(number_of_chromosomes_overlapped==2);

  hash_table_free(&hash_table);



  // ****
  //1.2 Fasta file that generate a graph with one long supernode, with a conflict at the end
  //   caused by two outward/exiting edges 
  //
  // >read1
  // ACATT
  // >read2
  // ATTC
  // >read3
  // ATTG
  // ****


  //first set up the hash/graph
  kmer_size = 3;
  number_of_bits = 4;
  bucket_size    = 4;
  bad_reads = 0;
  max_retries=10;

  hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  seq_loaded = load_population_as_fasta("../data/test/pop_graph/supernode/one_person_one_long_supernode_with_conflict_at_end",&bad_reads, hash_table);

  CU_ASSERT(seq_loaded==13);
  CU_ASSERT(bad_reads==0);


  load_chromosome_overlap_data("../data/test/pop_graph/dummy_chromosomes/simple2/chrom1.fasta", hash_table, 1);//does not overlap the graph
  load_chromosome_overlap_data("../data/test/pop_graph/dummy_chromosomes/simple2/chrom2.fasta", hash_table, 2);//does not overlap the graph
  load_chromosome_overlap_data("../data/test/pop_graph/dummy_chromosomes/simple2/chrom3.fasta", hash_table, 3);//overlaps node ATT
  load_chromosome_overlap_data("../data/test/pop_graph/dummy_chromosomes/simple2/chrom4.fasta", hash_table, 4);//overlaps node ATT

  //ACA IS  the end of a supernode in the reverse direction
  query_node1 = hash_table_find(element_get_key(seq_to_binary_kmer("ACA",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node1==NULL));
  CU_ASSERT(!db_node_is_supernode_end(query_node1, forward, individual_edge_array, 0, hash_table));
  CU_ASSERT(db_node_is_supernode_end(query_node1, reverse, individual_edge_array, 0, hash_table));

  number_of_chromosomes_overlapped=-100;
  CU_ASSERT(!db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome(query_node1, individual_edge_array, 0, hash_table, &number_of_chromosomes_overlapped));
  

  //ATT is not a supernode end
  query_node1 = hash_table_find(element_get_key(seq_to_binary_kmer("ATT",hash_table->kmer_size), hash_table->kmer_size), hash_table);
  CU_ASSERT(!(query_node1==NULL));
  CU_ASSERT(!db_node_has_at_most_one_intersecting_chromosome(query_node1, &number_of_chromosomes_overlapped));

  CU_ASSERT(number_of_chromosomes_overlapped==-1);

  
  hash_table_free(&hash_table);


  */


}





//void test_printing_supernode_with_chromosome_intersections_simple_alu_example()
// this deleted test was as folows. take a person whose geneome = an Alu. Take a chromosome which = same Alu. Make sure we do not call any SV.
//This should be covered in the new tests for ref based SV calling with the new path/supernode algorithm





void test_printing_supernode_with_chromosome_intersections_simple_alu_example_2()
{

  printf("Update this test to use new framework. Currentluy null test, but will be a useful example when we extend to cover translocations. Take one person, whose fasta is just an Alu.Take 2 chromosomes, one matches the start of the Alu, and the other which matches the end.So our one supernode (the alu) intersects 2 chromosomes, with no node intersecting >1 chromosome => is a  potential sv (translocation) locus\n");

  /*



 // ********************************
 // Another example. Take one person, whose fasta is just an Alu. 
 // Take 2 chromosomes, one matches the start of the Alu, and the other which matches the end.
 // So our one supernode (the alu) intersects 2 chromosomes, with no node intersecting >1 chromosome => is a  potential sv (translocation) locus.
 // **********************************

  //first set up the hash/graph
  int kmer_size = 31;
  int number_of_bits = 15;
  int bucket_size    = 10;
  long long bad_reads = 0;
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);



  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }


  load_population_as_fasta("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/just_one_of_the_two_people.txt",  &bad_reads,hash_table);
  
  CU_ASSERT(bad_reads ==0);

  char** array_of_supernodes_for_person1= (char**) calloc(316,sizeof(char*));
  if (array_of_supernodes_for_person1==NULL)
    {
      printf("unable to alloc the supernode array. dead before we even started. OOM");
      exit(1);
    }

  int i;
  for (i=0; i<316; i++)
    {
      array_of_supernodes_for_person1[i]="rubbish also";
    }

  char** array_of_chrom_overlaps_for_person1= (char**) calloc(316,sizeof(char*));
  if (array_of_chrom_overlaps_for_person1==NULL)
    {
      printf("unable to alloc the chrom overlaps. dead before we even started. OOM");
      exit(1);
    }

  for (i=0; i<316; i++)
    {
      array_of_chrom_overlaps_for_person1[i]="rubbish also";
    }
  

  //now load two chromosomes; one matches the start of the supernode, and the other matches the end of the supernode
  load_chromosome_overlap_data("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/matches_start_of_person1.fasta",hash_table, 1);
  load_chromosome_overlap_data("../data/test/pop_graph/test_pop_load_and_print/two_people_sharing_alu/matches_end_of_person1.fasta",hash_table, 2);

  //this counters are used to make sure none of the files of printed out supernodes get too big. In fact
  //it is completely unused in test code, as we don't print anything
  long supernode_count_person1=0;


  //this on the other hand is used in testing.
  int number_of_supernodes_that_are_potential_sv_loci=0;
  int number_of_chrom_overlaps_to_print_in_potential_sv_loci=0;//these two should end up being the same

  //print_supernode will, in debug mode, alloc memory for you in your array, and put the supernode in it
  int min_covg=0;
  int max_covg=1000;

  db_graph_traverse_specific_person_or_pop_for_supernode_and_chromosome_overlap_printing(&db_graph_choose_output_filename_and_print_potential_transloc_for_specific_person_or_pop, hash_table,                    
											 &supernode_count_person1, individual_edge_array, 0, min_covg, max_covg,
											 //false, NULL, NULL, &number_of_supernodes_that_are_potential_sv_loci,
											 //&number_of_chrom_overlaps_to_print_in_potential_sv_loci);
											 true, array_of_supernodes_for_person1,
  											 array_of_chrom_overlaps_for_person1, &number_of_supernodes_that_are_potential_sv_loci, 
											 &number_of_chrom_overlaps_to_print_in_potential_sv_loci);


  //there is one supernode and it should be seen as a potential SV locus
 CU_ASSERT((number_of_supernodes_that_are_potential_sv_loci==1));
 CU_ASSERT((number_of_chrom_overlaps_to_print_in_potential_sv_loci==1));

 // printf("SUPERNODE %s\n", array_of_supernodes_for_person1[0]);
 //printf("CHROMS start%sstop\n", array_of_chrom_overlaps_for_person1[0]);

 //should be the smae as or reverse complement of the sequence in person1, which is all one supernode
 CU_ASSERT(  !strcmp(array_of_supernodes_for_person1[0],"CTACGGCTGACTTTTTTTTTTTTTTTTTTTTAAGAGACGGGGTCTCGCTATGTTGCTCAGGCTGGAGTGCAGTGGCTATTCACAGGCGCGATCCCACTACTGATCAGCACGGGAGTTTTGACCTGCTCCGTTTCCGACCTGGGCCGGTTCACCCCTCCTTAGGCAACCTGGTGGTCCCCCGCTCCCGGGAGGTCACCATATTGATGCCGAACTTAGTGCGGACACCCGATCGGCATAGCGCACTACAGCCCAGAACTCCTGGACTCAAGCGATCCTCCCACCTCAGCCTCCCGAGTAGCTGGGACTACAGGCACGCGCCACCGCGCCCGGCCTCTGAAC") || !strcmp(array_of_supernodes_for_person1[0], "GTTCAGAGGCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGATCGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAGTTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGGTGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGATCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTTAAAAAAAAAAAAAAAAAAAAGTCAGCCGTAG") );

 // first 101 nodes intesect chom2, the last 78 nodes intersect chrom 1, and the middle (309-78-101)=130 don't intersect anything
 // all directions should be the same

 CU_ASSERT(!strcmp(array_of_chrom_overlaps_for_person1[0],"2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 2R 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R 1R ") || !strcmp(array_of_chrom_overlaps_for_person1[0], "1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 1F 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F 2F ") );


 free(array_of_supernodes_for_person1[0]);
 free(array_of_chrom_overlaps_for_person1[0]);
 free(array_of_supernodes_for_person1);
 free(array_of_chrom_overlaps_for_person1);
 hash_table_free(&hash_table);

  */


  
}



void test_printing_of_supernode_that_might_be_an_inversion_simple()
{


  printf("Update this test to use new framework. Currentluy null test, but is a simple example of an inversion - useful\n");


  /*

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

  
  // ****
  // fasta file will have ACAGATGCT
  // chrom will be        ACAATCGCT
  //                         *** <---- inversion marked with *
  //  If you draw it the overlaps are F 0 0 R 0 0 F
  // ****


  int seq_loaded = load_population_as_fasta("../data/test/pop_graph/sv/inversion/one_person", &bad_reads, hash_table);
  CU_ASSERT(seq_loaded==9);
  CU_ASSERT(bad_reads==0);


  load_chromosome_overlap_data("../data/test/pop_graph/sv/inversion/chromosome", hash_table, 1);


  char** array_of_supernodes_for_person1= (char**) calloc(10,sizeof(char*));

  if (array_of_supernodes_for_person1==NULL)
    {
      printf("cant start - OOM");
      exit(1);
    }
 int i;
 for (i=0; i<10; i++)
   {
     array_of_supernodes_for_person1[i]="rubbish";
   }

  char** array_of_chrom_overlaps_for_person1= (char**) calloc(10,sizeof(char*));

  if (array_of_chrom_overlaps_for_person1==NULL)
    {
      printf("cant start - OOM");
      exit(1);
    }
  for (i=0; i<10; i++)
    {
      array_of_chrom_overlaps_for_person1[i]="rubbish_chrom";
    }
  //this counters are used to make sure none of the files of printed out supernodes get too big. In fact
  //it is completely unused in test code, as we don't print anything
  long supernode_count_person1=0;

  
  //this on the other hand is used in testing.
  int number_of_supernodes_in_person_1=0;
  int number_of_chrom_overlaps_lists_in_person_1=0;//these two should end up being the same

  //print_supernode will, in debug mode, alloc memory for you in your array, and put the supernode in it
  int min_required_covg=0;
  int max_required_covg=1000;
  db_graph_traverse_specific_person_or_pop_for_supernode_and_chromosome_overlap_printing(&db_graph_choose_output_filename_and_print_potential_inversion_for_specific_person_or_pop, hash_table, 
											 &supernode_count_person1, individual_edge_array, 0, min_required_covg, max_required_covg,
											 //false, NULL, NULL, &number_of_supernodes_in_person_1, 
											//&number_of_chrom_overlaps_lists_in_person_1);
											 true, array_of_supernodes_for_person1,
											array_of_chrom_overlaps_for_person1, &number_of_supernodes_in_person_1, &number_of_chrom_overlaps_lists_in_person_1);



 CU_ASSERT((number_of_supernodes_in_person_1==1));
 CU_ASSERT_EQUAL(number_of_chrom_overlaps_lists_in_person_1,1);

 
 //  for (i=0; i< number_of_supernodes_in_person_1; i++)
 // {
 //  printf("SUPERNODE %s\n", array_of_supernodes_for_person1[i]);
 // }

 // for (i=0; i< number_of_chrom_overlaps_lists_in_person_1; i++)
 // {
 //   printf("CHROM XS %s\n", array_of_chrom_overlaps_for_person1[i]);
 // }



 
 CU_ASSERT_STRING_EQUAL(array_of_supernodes_for_person1[0],"ACAGATGCT");
 CU_ASSERT_STRING_EQUAL(array_of_chrom_overlaps_for_person1[0],"1F 00 00 1R 00 00 1F ");

 //cleanup
 for (i=0; i< number_of_supernodes_in_person_1; i++)
   {
     free(array_of_supernodes_for_person1[i]);
   }
 for (i=0; i< number_of_chrom_overlaps_lists_in_person_1; i++)
    {
    free(array_of_chrom_overlaps_for_person1[i]);
  }

 free(array_of_supernodes_for_person1);
 free(array_of_chrom_overlaps_for_person1);
 hash_table_free(&hash_table);

  */  
}







