/*
  hash_table.c -- implementation
 */

#include <stdlib.h>
#include <stdio.h>
#include <hash_table.h>
#include <priority_queue.h>
#include <hash_value.h>

HashTable * hash_table_new(int number_buckets, short kmer_size){ 
  
  HashTable *hash_table = malloc(sizeof(HashTable));

  hash_table->table = calloc(number_buckets, sizeof(PQueue));
  
  if (hash_table == NULL || hash_table->table == NULL) {
    printf("could not allocate hash table of size %i\n",number_buckets);
    exit(1);
  }
  
  hash_table->number_buckets = number_buckets;
  hash_table->kmer_size = kmer_size;
  return hash_table;
}

void hash_table_free(HashTable ** hash_table)
{ 
  int i;

  for(i=0;i<(*hash_table)->number_buckets;i++){
    pqueue_free_elements(&(*hash_table)->table[i]);
  }
  
  free((*hash_table)->table);
  free(*hash_table);
  *hash_table = NULL;
}


boolean hash_table_apply_or_insert(Key key,void (*f)(Element *), HashTable * hash_table){

  int hashval = hash_value(key,hash_table->number_buckets);
 
  if (hash_table == NULL) {
    puts("NULL table!");
    exit(1);
  }
  
  return pqueue_apply_or_insert(key,f,&(hash_table->table[hashval]),hash_table->kmer_size);
}


void hash_table_traverse(void (*f)(Element *),HashTable * hash_table){
  int i;
  for(i=0;i<hash_table->number_buckets;i++){
    pqueue_traverse(f,&(hash_table->table[i]));
  }
}



Element * hash_table_find(Key key, HashTable * hash_table){
  
  int hashval = hash_value(key,hash_table->number_buckets);

  return pqueue_find(key,&(hash_table->table[hashval]), hash_table->kmer_size);
}




Element * hash_table_find_or_insert(Key key, boolean * found, HashTable * hash_table){
  
  int hashval = hash_value(key,hash_table->number_buckets);

  return pqueue_find_or_insert(key,found, &(hash_table->table[hashval]), hash_table->kmer_size);
}




