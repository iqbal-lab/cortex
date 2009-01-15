
#include <pqueue_pop.h>

void pqueue_traverse_specific_person_or_pop(void (*f)(HashTable*, Element *, EdgeArrayType, int, char**, int*),HashTable* hash_table, PQueue * pqueue, EdgeArrayType type, int index, char** for_test, int* index_for_test)
{
  int i;
  for(i=0;i<pqueue->number_entries;i++)
    {
      f(hash_table, &(pqueue->elements[i]), type, index, for_test, index_for_test);
    }
}
