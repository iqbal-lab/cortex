#include <stdlib.h>
#include <stdio.h>
#include <priority_queue.h>


//assumes the element is not already stored -- check this before calling this function
void pqueue_bubble_up(int current_index, PQueue *pqueue){
  
  while (current_index > 0){

    int parent_index = (current_index-1)/2;
    if (element_smaller(pqueue->elements[parent_index],pqueue->elements[current_index])){
      Element tmp;
      tmp = pqueue->elements[current_index];
      pqueue->elements[current_index] = pqueue->elements[parent_index];
      pqueue->elements[parent_index] = tmp;
      current_index = parent_index;
    }
    else{
      break;
    }
  }
}

boolean pqueue_apply_or_insert(Key key, void (*f)(Element*),PQueue *pqueue, short kmer_size){

  int i =0;
  boolean found = false;
  Element element;

  for(i=0;i<pqueue->number_entries;i++){
    if (element_is_key(key,pqueue->elements[i], kmer_size)){
      f(&pqueue->elements[i]);
      pqueue_bubble_up(i,pqueue);
      found = true;
      break;
    }
  }

  if (! found){
    int current_index = pqueue->number_entries;
    pqueue->number_entries++;
    
    Element * new_ptr;

    new_ptr = realloc(pqueue->elements,(current_index+1) * sizeof(Element));
    
    if (new_ptr == NULL){
      puts("priority queue: cannot allocate memory");
      exit(1);//so no need to worry about orphaned pointer at this stage.
    }
    
    pqueue->elements = new_ptr;

    element_initialise(&element,key, kmer_size);
    
    pqueue->elements[current_index] = element; //structure assignment

  }
  return found;
}

void pqueue_traverse(void (*f)(Element *),PQueue * pqueue)
{
  int i;
  for(i=0;i<pqueue->number_entries;i++){
    f(&(pqueue->elements[i]));
  }
}


Element * pqueue_find_or_insert(Key key,boolean * found, PQueue * pqueue, short kmer_size){  
  int i;
  Element element;

  *found = false;
  for(i=0;i<pqueue->number_entries;i++){ 
    if (element_is_key(key,pqueue->elements[i], kmer_size)){
      *found = true;
      return &(pqueue->elements[i]);
    }
  }
  
  int current_index = pqueue->number_entries;
  pqueue->number_entries++;
  
  Element * new_ptr;

  new_ptr = realloc(pqueue->elements,(current_index+1) * sizeof(Element));
    
  if (new_ptr == NULL){
    puts("priority queue: cannot allocate memory");
    exit(1);//so no need to worry about orphaned pointer at this stage.
  }
  
  pqueue->elements = new_ptr;

  element_initialise(&element,key, kmer_size);

  pqueue->elements[current_index] = element; //structure assignment!
  

  return &(pqueue->elements[i]);  
}

Element * pqueue_find(Key key,PQueue * pqueue, short kmer_size){
  
  int i;

  for(i=0;i<pqueue->number_entries;i++){ 
    if (element_is_key(key,pqueue->elements[i], kmer_size)){
      return &(pqueue->elements[i]);
    }
  }
  
  return NULL;
}

void pqueue_free(PQueue** pqueue){
  
  free((*pqueue)->elements);
  
  free(*pqueue);
  *pqueue = NULL;
  
}

void pqueue_free_elements(PQueue* pqueue){
  
  free(pqueue->elements);
  pqueue->elements = NULL;  
}





