#include <CUnit.h>
#include <Basic.h>
#include <file_reader.h>
#include <dB_graph_population.h>
#include <element.h>
#include <hash_table.h>
#include <stdlib.h>

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

  char** array_of_supernodes_for_person1= malloc(10*sizeof(char*));
  char** array_of_supernodes_for_person2= malloc(10*sizeof(char*));

  if ( (array_of_supernodes_for_person1==NULL) || (array_of_supernodes_for_person2==NULL))
    {
      printf("cant start - OOM");
      exit(1);
    }


  int number_of_supernodes_in_person_1=0;
  db_graph_traverse_specific_person_or_pop(&print_supernode_for_specific_person_or_pop, hash_table, individual_edge_array, 0, true, array_of_supernodes_for_person1,&number_of_supernodes_in_person_1);
  printf("PERSON 1 has %d supernodes\n", number_of_supernodes_in_person_1);
  db_graph_set_all_visited_nodes_to_status_none(hash_table);

  int number_of_supernodes_in_person_2=0;
  db_graph_traverse_specific_person_or_pop(&print_supernode_for_specific_person_or_pop, hash_table, individual_edge_array, 1, true, array_of_supernodes_for_person2,&number_of_supernodes_in_person_2);
  printf("PERSON 2 has %d supernodes\n", number_of_supernodes_in_person_2);

  CU_ASSERT(number_of_supernodes_in_person_1==number_of_supernodes_in_person_2);

  int i;
  for (i=0; i<number_of_supernodes_in_person_1; i++)
    {
      printf("Person 1 has supernode %s\n", array_of_supernodes_for_person1[i]);
    }
  for (i=0; i<number_of_supernodes_in_person_2; i++)
    {
      printf("Person 2 has supernode %s\n", array_of_supernodes_for_person2[i]);
    }


  //Expected results are
  char** correct_answer_person_1 = malloc(2*sizeof(char*));
  correct_answer_person_1[0]= malloc(4*sizeof(char));
  correct_answer_person_1[1]= malloc(6*sizeof(char));
  correct_answer_person_1[0]="AAA";

  char** correct_answer_person_2 = malloc(2*sizeof(char*));
  correct_answer_person_2[0]= malloc(*sizeof(char));
  correct_answer_person_2[1]= malloc(*sizeof(char));






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
