#ifndef PQ_POP_H_
#define PQ_POP_H_

#include<priority_queue.h>
#include <element.h>
#include <hash_table.h>

//void pqueue_traverse_specific_person_or_pop(void (*f)(HashTable*, Element*, EdgeArrayType, int ),HashTable* hash_table, PQueue * pqueue, EdgeArrayType type, int index);

void pqueue_traverse_specific_person_or_pop_for_supernode_printing(void (*f)(HashTable*, Element *, long*, EdgeArrayType, int, boolean, char**, int*),
								   HashTable* hash_table,  PQueue * pqueue, long* supernode_count, EdgeArrayType type, int index, 
								   boolean is_for_testing, char** for_test, int* index_for_test);

void pqueue_traverse_specific_person_or_pop_for_supernode_and_chromosome_overlap_printing(void (*f)(HashTable*, Element *, long*, EdgeArrayType, int, int, int, boolean, char**, char**, int*,  int*),
                                                                                          HashTable* hash_table,  PQueue * pqueue,  long* supernode_count,  EdgeArrayType type, int index,
											  int min_covg, int max_covg,
                                                                                          boolean is_for_testing, char** for_test1, char** for_test2, int* index_for_test1, int* index_for_test2);

void pqueue_traverse_to_gather_statistics_about_people(void (*f)(HashTable*, Element *, int**, int), PQueue *, HashTable *, int**, int);


#endif
