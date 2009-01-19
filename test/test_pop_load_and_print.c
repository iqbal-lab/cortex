#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <dB_graph_population.h>
#include <element.h>
#include <hash_table.h>
#include <stdlib.h>
#include "supernode_cmp.h"

void test_load_two_people_in_same_populations_and_print_separately_their_supernodes()
{

  int kmer_size = 3;
  int number_of_buckets=5;
  HashTable* hash_table = hash_table_new(number_of_buckets,kmer_size);
  
  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  int seq_loaded = load_population_as_fasta("../test/data/test_pop_load_and_print/two_individuals_simple.txt", hash_table);
  //printf("Number of bases loaded is %d",seq_loaded);
  //CU_ASSERT(seq_loaded == 33);

  char** array_of_supernodes_for_person1= (char**) calloc(10,sizeof(char*));
  char** array_of_supernodes_for_person2= (char**) calloc(10,sizeof(char*));

  if ( (array_of_supernodes_for_person1==NULL) || (array_of_supernodes_for_person2==NULL))
    {
      printf("cant start - OOM");
      exit(1);
    }


  int number_of_supernodes_in_person_1=0;
  db_graph_traverse_specific_person_or_pop(&print_supernode_for_specific_person_or_pop, hash_table, individual_edge_array, 0, true, array_of_supernodes_for_person1,&number_of_supernodes_in_person_1);
  //printf("PERSON 1 has %d supernodes\n", number_of_supernodes_in_person_1);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  int number_of_supernodes_in_person_2=0;
  db_graph_traverse_specific_person_or_pop(&print_supernode_for_specific_person_or_pop, hash_table, individual_edge_array, 1, true, array_of_supernodes_for_person2,&number_of_supernodes_in_person_2);
  //printf("PERSON 2 has %d supernodes\n", number_of_supernodes_in_person_2);

  CU_ASSERT(number_of_supernodes_in_person_1==2);
  CU_ASSERT(number_of_supernodes_in_person_2==2);

  int i;
  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      printf("Person 1 has supernode %s\n", array_of_supernodes_for_person1[i]);
    }
  for (i=0; i<number_of_supernodes_in_person_2; i++)
    {
      printf("Person 2 has supernode %s\n", array_of_supernodes_for_person2[i]);
    }


  //Quicksort these results:
  qsort((void*) array_of_supernodes_for_person1, sizeof(array_of_supernodes_for_person1)/sizeof(char*), sizeof(char*),  supernode_cmp);
  qsort((void*) array_of_supernodes_for_person2, sizeof(array_of_supernodes_for_person2)/sizeof(char*), sizeof(char*),  supernode_cmp);



  //Expected results are
  char** correct_answer_person_1 = (char**) malloc(30*sizeof(char*));

  correct_answer_person_1[0]= (char*) malloc(4*sizeof(char));
  correct_answer_person_1[1]= (char*) malloc(7*sizeof(char));
  if ((correct_answer_person_1==NULL) || (correct_answer_person_1[0]==NULL) || (correct_answer_person_1[1]==NULL) )
    {
      printf("OOM - malloc answers");
      exit(1);
    }

  correct_answer_person_1[0]="AAA";  
  correct_answer_person_1[1]="AACGTT";

  char** correct_answer_person_2 = (char**) malloc(3*sizeof(char*));
  correct_answer_person_2[0]= (char*) malloc(4*sizeof(char));
  correct_answer_person_2[1]= (char*) malloc(7*sizeof(char));
  if ((correct_answer_person_2==NULL) || (correct_answer_person_2[0]==NULL) || (correct_answer_person_2[1]==NULL) )
    {
      printf("OOM - malloc answers");
      exit(1);
    }

  correct_answer_person_2[0]="CCC";
  correct_answer_person_2[1]="CGTCAA";


  printf("person1 should have supernodes %s and %s\n", correct_answer_person_1[0], correct_answer_person_1[1]);
  printf("person2 should have supernodes %s and %s\n", correct_answer_person_2[0], correct_answer_person_2[1]);
  //Quicksort the expected answers:
  
  size_t corr1_size = sizeof(correct_answer_person_1)/sizeof(char*);
  size_t corr2_size = sizeof(correct_answer_person_2)/sizeof(char*);

  qsort((void*) correct_answer_person_1, corr1_size, sizeof(char*),  supernode_cmp);
  qsort((void*) correct_answer_person_2, corr2_size, sizeof(char*),  supernode_cmp);



  int j;

  CU_ASSERT(2==number_of_supernodes_in_person_1);
  //  for (j=0; j<2; j++)
  // {
  //   CU_ASSERT_STRING_EQUAL(correct_answer_person_1[j], array_of_supernodes_for_person1[j]);
  // }


  CU_ASSERT_STRING_EQUAL(correct_answer_person_1[0], array_of_supernodes_for_person1[0]);

  printf(" person 2 first corr %s and got %s\n", correct_answer_person_2[0],  array_of_supernodes_for_person2[0] );
  printf(" person 2 second corr %s and got %s\n", correct_answer_person_2[1],  array_of_supernodes_for_person2[1] );

  //  CU_ASSERT_STRING_EQUAL(correct_answer_person_2[0], array_of_supernodes_for_person2[0]);


  CU_ASSERT(2==number_of_supernodes_in_person_2);
  //for (j=0; j<2; j++)
  // {
  //   CU_ASSERT_STRING_EQUAL(correct_answer_person_2[j], array_of_supernodes_for_person2[j]);
  // }




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
}
