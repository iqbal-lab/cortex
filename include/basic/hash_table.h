/*
  hash_table.h 

  all the routines as prefixed with hash_table
*/


#ifndef HASH_H_
#define HASH_H_

#include <element.h>
#include <priority_queue.h>


typedef struct{
  short kmer_size;
  unsigned int number_buckets;
  PQueue * table; //every bucket is implemented as a priority queue
} HashTable;


HashTable * hash_table_new(int number_buckets, short kmer_size);

void hash_table_free(HashTable * * hash_table);

//if the key is present applies f otherwise adds a new element for kmer
boolean hash_table_apply_or_insert(Key key, void (*f)(Element*), HashTable *);

//applies f to every element of the table
void hash_table_traverse(void (*f)(Element *),HashTable *);

//if the element is not in table create an element with key and adds it
Element * hash_table_find_or_insert(Key key, boolean * found, HashTable * hash_table);

//return entry for kmer
Element * hash_table_find(Key key, HashTable * hash_table);


#endif /* HASH_H_ */
