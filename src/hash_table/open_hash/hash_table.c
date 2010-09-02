/*
  hash_table.c -- implementation
 */

#include <stdlib.h>
#include <stdio.h>
#include <open_hash/hash_table.h>
#include <hash_value.h>
#include <assert.h>


HashTable * hash_table_new(int number_bits, int bucket_size, int max_rehash_tries, short kmer_size){ 
  
  HashTable *hash_table = malloc(sizeof(HashTable));

  if (hash_table == NULL) {
    fprintf(stderr,"could not allocate hash table of size %qd\n", (long long) bucket_size * (1 << number_bits));
    return NULL;
    //exit(1);
  }
  
  hash_table->collisions = calloc(max_rehash_tries, sizeof(long long));
  if (hash_table->collisions == NULL) {
    fprintf(stderr,"could not allocate memory\n");
    return NULL;
    //exit(1);
  }
  
  hash_table->unique_kmers = 0;
  hash_table->max_rehash_tries = max_rehash_tries;
  hash_table->number_buckets = (long long) 1 << number_bits;
  hash_table->bucket_size   = bucket_size;

  //calloc is vital - we want to make sure initialised to zero
  hash_table->table = calloc(hash_table->number_buckets * hash_table->bucket_size, sizeof(Element));
 
  if (hash_table->table == NULL) {
    fprintf(stderr,"could not allocate hash table of size %qd\n",hash_table->number_buckets * hash_table->bucket_size);
    return NULL;
    //exit(1);
  }
  
  hash_table->kmer_size      = kmer_size;
  return hash_table;
}

void hash_table_free(HashTable ** hash_table)
{ 
  free((*hash_table)->table);
  free((*hash_table)->collisions);
  free(*hash_table);
  *hash_table = NULL;
}


// Lookup for key in bucket defined by the hash value. 
// If key is in bucket, returns true and the position of the key/element in current_pos.
// If key is not in bucket, and bucket is not full, returns the next available position in current_pos (and overflow is returned as false)
// If key is not in bucket, and bucket is full, returns overflow=true
boolean hash_table_find_in_bucket(Key key, long long * current_pos, boolean * overflow, HashTable * hash_table, int rehash){


  //add the rehash to the final bitfield in the BinaryKmer
  BinaryKmer bkmer_with_rehash_added;
  binary_kmer_initialise_to_zero(&bkmer_with_rehash_added);
  binary_kmer_assignment_operator(bkmer_with_rehash_added, *key);
  bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] =   bkmer_with_rehash_added[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]+ (bitfield_of_64bits) rehash;

  int hashval = hash_value(&bkmer_with_rehash_added,hash_table->number_buckets);


  boolean found = false;
  int i=0;                     //position in bucket
  *overflow    = false;
  *current_pos   = (long long) hashval * hash_table->bucket_size;   //position in hash table

  while( (i<hash_table->bucket_size) &&   // still within the bucket
	 (!db_node_check_status(&hash_table->table[*current_pos],unassigned) )  && // not yet reached an empty space
	 (!found)
	 )
    {

    //sanity check -- to avoid out of boundary access
    if (*current_pos >= hash_table->number_buckets * hash_table->bucket_size || *current_pos<0)
      {
	printf("out of bounds problem found in hash table_find_with_position\n");
	exit(1);
      }
    
    //element found

   
    if (element_is_key(key,hash_table->table[*current_pos], hash_table->kmer_size))
      {
	found = true;
      }
    else
      {
	(*current_pos)++;
	i++;
      }
    
    }
  

  if (i == hash_table->bucket_size)
    {
      *overflow = true;
    }
  

  assert(!found || !(*overflow));
  return found;
}


//currently not used, and must add a test
boolean hash_table_apply_or_insert(Key key, void (*f)(Element *), HashTable * hash_table){
  if (hash_table == NULL) {
    puts("NULL table!");
    exit(1);
  }
  
  long long current_pos;
  Element element;
  boolean overflow;
  int rehash=0;
  boolean found;
  do
    {
      found = hash_table_find_in_bucket(key,&current_pos,&overflow, hash_table,rehash);
      
      if (!found)
	{
	  if (!overflow)
	    {
	      //sanity check
	      if (!db_node_check_status(&hash_table->table[current_pos],unassigned)){
		printf("Out of bounds - trying to insert new node beyond end of bucket\n");
		exit(1);
	      }

	      element_initialise(&element,key, hash_table->kmer_size);
	      element_assign( &(hash_table->table[current_pos]),  &element); 
	      hash_table->unique_kmers++;
	    }
	  else//overflow
	    {
	      rehash++;
	      if (rehash>hash_table->max_rehash_tries)
		{
		  printf("too much rehashing!! Rehash=%d\n", rehash);
		  exit(1);
		}
	    }
	}
      else
	{
	  f(&hash_table->table[current_pos]);
	}


    }while(overflow);

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

long long hash_table_traverse_returning_sum(long long (*f)(Element *),HashTable * hash_table){
  long long i;
  long long ret=0;
  for(i=0;i<hash_table->number_buckets * hash_table->bucket_size;i++){
    if (!db_node_check_status(&hash_table->table[i],unassigned)){
      ret += f(&hash_table->table[i]);
    }
  }
  return ret;
}



Element * hash_table_find(Key key, HashTable * hash_table)
{
  if (hash_table == NULL) 
    {
      puts("hash_table_find has been called with a NULL table! Exiting");
      exit(1);
    }

  Element * ret = NULL;
  long long current_pos;
  boolean overflow;
  int rehash = 0;
  boolean found; 

  do
    {
      found = hash_table_find_in_bucket(key,&current_pos, &overflow, hash_table,rehash);
      
      if (found) //then we know overflow is false - this is checked in find_in_bucket
	{
	  ret =  &hash_table->table[current_pos];
	}
      else if (overflow)
	{ //rehash
	  rehash++; 
	  if (rehash>hash_table->max_rehash_tries)
	    {
	      printf("too much rehashing!! Rehash=%d\n", rehash);
	      exit(1);
	    }
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

    *found = hash_table_find_in_bucket(key,&current_pos,&overflow,hash_table,rehash);

    if (! *found)
      {
	if (!overflow) //it is definitely nowhere in the hashtable, so free to insert
	  {
	    //sanity check
	    if (!db_node_check_status(&hash_table->table[current_pos],unassigned))
	      {
		fprintf(stderr,"error trying to write on an occupied element\n");
		exit(1);
	      }
      
	    //insert element
	    //printf("Inserting element at position %qd in bucket \n", current_pos);
	    element_initialise(&element,key, hash_table->kmer_size);

	    //hash_table->table[current_pos] = element; //structure assignment
	    element_assign(&(hash_table->table[current_pos]) , &element);
	    
	    ret = &hash_table->table[current_pos];
	    hash_table->unique_kmers++;
	    
	  }
	else
	  { //overflow -> rehashing
	    
	    rehash++;
	    if (rehash>hash_table->max_rehash_tries)
	      {
		printf("too much rehashing!! Rehash=%d\n", rehash);
		exit(1);
	      }
	  }
      }
    else //it is found
      {
	ret = &hash_table->table[current_pos];
      }
  } while (overflow);
  
  hash_table->collisions[rehash]++;
  return ret;
}


void hash_table_print_stats(HashTable * hash_table)
{
  int k;
  fprintf(stderr,"Collisions:\n");
  for(k=0;k<10;k++)
    {
      if (hash_table->collisions[k] != 0){
	fprintf(stderr,"\t tries %i: %qd\n",k,hash_table->collisions[k]);
      }
    }
}



long long hash_table_get_unique_kmers(HashTable * hash_table)
{
  return hash_table->unique_kmers;
}



