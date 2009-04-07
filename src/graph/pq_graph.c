#include <pq_graph.h>
#include <open_hash/hash_table.h>
#include <element.h>
#include <priority_queue.h>


void pqueue_traverse_with_array(void (*f)(HashTable*, Element *, int**, int),HashTable* hash_table, PQueue * pqueue, int** array_of_ints, int length_of_array)
{
  int i;
  for(i=0;i<pqueue->number_entries;i++){
    f(hash_table, &(pqueue->elements[i]), array_of_ints, length_of_array);
  }
}
