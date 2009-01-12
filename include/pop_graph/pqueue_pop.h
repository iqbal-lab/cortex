#ifndef PQ_POP_H_
#define PQ_POP_H_

#include<priority_queue.h>
#include <element.h>
#include <hash_table.h>

void pqueue_traverse_specific_person_or_pop(void (*f)(HashTable*, Element*, EdgeArrayType, int ),HashTable* hash_table, PQueue * pqueue, EdgeArrayType type, int index);


#endif
