
/*
  pq_graph.h
  
Wrapper for the pq

*/

#ifndef PQ_GRAPH_H_
#define PQ_GRAPH_H_


#include <hash_table.h>
#include <priority_queue.h>
#include <element.h>

void pqueue_traverse_with_array(void (*f)(HashTable*, Element *, int**, int), HashTable* hash_table, PQueue * pqueue, int** array_of_ints, int length_of_array);

#endif /* PQ_GRAPH_H_ */
