
#include <pqueue_pop.h>

void pqueue_traverse_specific_person_or_pop_for_supernode_printing(void (*f)(HashTable*, Element *, long*, EdgeArrayType, int, boolean, char**, int*),HashTable* hash_table,  PQueue * pqueue, long* supernode_count, 
					    EdgeArrayType type, int index, boolean is_for_testing, char** for_test, int* index_for_test)
{
  int i;
  for(i=0;i<pqueue->number_entries;i++)
    {
      f(hash_table, &(pqueue->elements[i]), supernode_count, type, index, is_for_testing, for_test, index_for_test);
    }
}


void pqueue_traverse_specific_person_or_pop_for_supernode_and_chromosome_overlap_printing(void (*f)(HashTable*, Element *, long*, EdgeArrayType, int, int, int, boolean, char**, char**, int*,  int*),
											  HashTable* hash_table,  PQueue * pqueue,  long* supernode_count,  EdgeArrayType type, int index, 
											  int min_covg, int max_covg,
											  boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2)
{
  int i;
  for(i=0;i<pqueue->number_entries;i++)
    {
      f(hash_table, &(pqueue->elements[i]), supernode_count, type, index, min_covg, max_covg, is_for_testing, for_test1, for_test2, index_for_test1, index_for_test2);
    }
}

void pqueue_traverse_to_gather_statistics_about_people(void (*f)(HashTable*, Element *, int**, int), PQueue * pqueue, HashTable * hash_table, int** array, int num_people)
{
  int i;
  for(i=0;i<pqueue->number_entries;i++){
    f(hash_table, &(pqueue->elements[i]), array, num_people);
  }

}
