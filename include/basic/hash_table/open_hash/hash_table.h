/*
  hash_table.h 
  
  open hash table implementation - ie every bucket has a predifined size 
  overloads results in rehashing
  all the routines as prefixed with hash_table
*/


#ifndef HASH_H_
#define HASH_H_

#include <element.h>



typedef struct{
  short kmer_size;
  long long number_buckets;
  int bucket_size;
  Element * table; 
  long long * collisions;
  long long unique_kmers;
  int max_rehash_tries;
} HashTable;


HashTable * hash_table_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size);

void hash_table_free(HashTable * * hash_table);

//if the key is present applies f otherwise adds a new element for kmer
boolean hash_table_apply_or_insert(Key key, void (*f)(Element*), HashTable *);

//applies f to every element of the table
void hash_table_traverse(void (*f)(Element *),HashTable *);

//if the element is not in table create an element with key and adds it
Element * hash_table_find_or_insert(Key key, boolean * found, HashTable * hash_table);

void hash_table_print_stats(HashTable *);

long long hash_table_get_unique_kmers(HashTable *);

//return entry for kmer
Element * hash_table_find(Key key, HashTable * hash_table);

#endif /* HASH_H_ */
