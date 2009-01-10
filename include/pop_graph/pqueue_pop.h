#ifndef PQ_POP_H_
#define PQ_POP_H_

#include<priority_queue.h>
#include <element.h>


void pqueue_traverse_specific_person_or_pop(void (*f)(Element *, EdgeArrayType, int ),PQueue * pqueue, EdgeArrayType type, int index);


#endif
