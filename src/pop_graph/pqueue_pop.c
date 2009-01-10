
#include <pqueue_pop.h>

void pqueue_traverse_specific_person_or_pop(void (*f)(Element *, EdgeArrayType, int),PQueue * pqueue, EdgeArrayType type, int index)
{
  int i;
  for(i=0;i<pqueue->number_entries;i++){
    f(&(pqueue->elements[i]), type, index);
  }
}
