/*
  hash_table.c -- implementation
 */

#include <stdlib.h>
#include <stdio.h>
#include <open_hash/hash_table.h>
#include <hash_value.h>

HashTable * hash_table_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size){ 
  
  HashTable *hash_table = malloc(sizeof(HashTable));

  if (hash_table == NULL) {
    fprintf(stderr,"could not allocate hash table of size %qd\n", (long long) bucket_size * (1 << number_bits));
    exit(1);
  }
  
  hash_table->collisions = calloc(max_rehash_tries, sizeof(long long));
  if (hash_table->collisions == NULL) {
    fprintf(stderr,"could not allocate meme\n");
    exit(1);
  }
  
  hash_table->unique_kmers = 0;

  hash_table->max_rehash_tries = max_rehash_tries;

  hash_table->number_buckets = (long long) 1 << number_bits;
  hash_table->bucket_size   = bucket_size;
  
  hash_table->table = calloc(hash_table->number_buckets * hash_table->bucket_size, sizeof(Element));
  
 
  if (hash_table->table == NULL) {
    fprintf(stderr,"could not allocate hash table of size %qd\n",hash_table->number_buckets * hash_table->bucket_size);
    exit(1);
  }
  
  hash_table->kmer_size      = kmer_size;
  return hash_table;
}

void hash_table_free(HashTable ** hash_table)
{ 
  free((*hash_table)->table);
  free(*hash_table);
  *hash_table = NULL;
}

boolean hash_table_find_with_position(Key key, long long * current_pos, boolean * overflow, HashTable * hash_table, int rehash){

  int hashval = hash_value(key+rehash,hash_table->number_buckets);

  boolean found = false;
  int i=0;

  *overflow    = false;
  *current_pos   = (long long) hashval * hash_table->bucket_size;
  
   
  while(i<hash_table->bucket_size){
    
    printf("in the hash table with postion while loop, hashval is %qd and status is %d\n", hashval, (hash_table->table[*current_pos]).status);
    //sanity check -- to avoid out of boundery access
    if (*current_pos >= hash_table->number_buckets * hash_table->bucket_size || *current_pos<0){
      printf("problem found\n");
      exit(1);
    }

    //reach an empty space
    if (db_node_check_status(&hash_table->table[*current_pos],unassigned)){
      break;
    }
    else
      {
	if (hash_table->table[*current_pos].status == none)
	  {
	    printf("Did not break - status was NONE\n");
	  }
	else if (hash_table->table[*current_pos].status == visited)
          {
            printf("Did not break - status was VISITED\n");
          }
	else if (hash_table->table[*current_pos].status == unassigned)
          {
            printf("Did not break - status was UNASSIGNED\n");
          }


      }

    //element found
    if (element_is_key(key,hash_table->table[*current_pos], hash_table->kmer_size)){
       found = true;
       break;
     }

    (*current_pos)++;
    i++;
  }

  if (i == hash_table->bucket_size){
    *overflow = true;
  }
  
  //assert(!found || !overflow);

  return found;
}


boolean hash_table_apply_or_insert(Key key, void (*f)(Element *), HashTable * hash_table){
  if (hash_table == NULL) {
    puts("NULL table!");
    exit(1);
  }
  
  long long current_pos;
  Element element;
  boolean overflow; 
  boolean found = hash_table_find_with_position(key,&current_pos,&overflow, hash_table,0);
  
  if (found){
    f(&hash_table->table[current_pos]);
  }
  else{
    if (overflow){
      //sanity check
      if (!db_node_check_status(&hash_table->table[current_pos],unassigned)){

	exit(1);
      }

      element_initialise(&element,key, hash_table->kmer_size);
      hash_table->table[current_pos] = element; //structure assignment
    }
    else{
      fprintf(stderr,"bucket overflow!");
    }
  }
  
  return found;
     
}


void hash_table_traverse(void (*f)(Element *),HashTable * hash_table){
  long long i;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      f(&hash_table->table[i]);
    }
  }
}




Element * hash_table_find(Key key, HashTable * hash_table){
  if (hash_table == NULL) {
    puts("NULL table!");
    exit(1);
  }
  Element * ret = NULL;
  long long current_pos;
  boolean overflow;
  int rehash = 0;
  boolean found; 

  do{
    found = hash_table_find_with_position(key,&current_pos, &overflow, hash_table,rehash);
    
    if (found){
      ret =  &hash_table->table[current_pos];
    }

    if (overflow){ //rehash
       rehash++; 
    }
  } while(overflow);
  
  return ret;
}
  

Element * hash_table_find_or_insert(Key key, boolean * found,  HashTable * hash_table){
  
  if (hash_table == NULL) {
    puts("NULL table!");
    exit(1);
  }
  
  Element element;
  Element * ret = NULL;
  int rehash = 0;
  boolean overflow; 

  long long current_pos;

  do{

    *found = hash_table_find_with_position(key,&current_pos,&overflow,hash_table,rehash);

    if (! *found){
      if (!overflow){
	//sanity check
	if (!db_node_check_status(&hash_table->table[current_pos],unassigned)){
	  fprintf(stderr,"error trying to write on an occupied element\n");
	  exit(1);
	}
      
	//insert element
	element_initialise(&element,key, hash_table->kmer_size);
	printf("after elem initialise the elemtn status is %d\n", element.status);
	hash_table->table[current_pos] = element; //structure assignment
	ret = &hash_table->table[current_pos];
	hash_table->unique_kmers++;

      }
      else{ //overflow -> rehashing
	
	rehash++;
	if (rehash>hash_table->max_rehash_tries){
	  fprintf(stderr,"too much rehashing!!\n");
	  exit(1);
	}
      }
    }
    else{
      ret = &hash_table->table[current_pos];
    }
  } while (overflow);
  
  hash_table->collisions[rehash]++;
  return ret;
}


void hash_table_print_stats(HashTable * hash_table){
  int k;
  fprintf(stderr,"Collisions:\n");
  for(k=0;k<10;k++){
    if (hash_table->collisions[k] != 0){
      fprintf(stderr,"\t tries %i: %qd\n",k,hash_table->collisions[k]);
    }
  }
}



long long hash_table_get_unique_kmers(HashTable * hash_table){
  return hash_table->unique_kmers;
}



