/*
 * priority_queue.h
 
 routines are prefixed with priority_queue_

 */



#ifndef PQUEUE_H_
#define PQUEUE_H_

#include <element.h>


typedef struct{
  int number_entries;
  Element * elements;
} PQueue;


PQueue * pqueue_new();


boolean pqueue_apply_or_insert(Key, void (*f)(Element*), PQueue *, short kmer_size);

Element * pqueue_find_or_insert(Key key,boolean * found, PQueue *pqueue, short kmer_size);

void pqueue_traverse(void (*f)(Element *),PQueue *);

Element * pqueue_find(Key,PQueue *, short kmer_size);

void pqueue_free(PQueue** pqueue);
void pqueue_free_elements(PQueue* pqueue);

#endif /* PQUEUE_H_ */
