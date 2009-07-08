#include <CUnit.h>
#include <Basic.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <stdlib.h>


void test_hash_table_find_or_insert()
{

  short kmer_size;
  long long max_key_given_kmer_size;

      

  //adds binary kmers to a hash, then tries to find them
  //will go through all possible bin kmers given the kmer_size, stepping 
  // with granulatiry step. Since for kmer_size large (eg 21) it takes to long to try ALL
  // possible 21-mers, the last argument allows you to specify the largest one to test

  
  void test(short kmer_size, int num_bits, int bucket, int max_tries, int step, long long max_kmer_to_test)
    {
      
      int number_of_bits      = num_bits;
      int bucket_size         = bucket;
      long long bad_reads     = 0; 
      int max_retries         = max_tries;
      boolean found           = false;
      
      HashTable* hash_table  = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
      
      
      long long i;
          
      for (i=0; i< max_kmer_to_test; i++)
	{
	  hash_table_find_or_insert(element_get_key(i, hash_table->kmer_size),&found, hash_table);
	  if (found==true)
	    {
	      CU_ASSERT(i>element_get_key(i, hash_table->kmer_size)); //the reverse complement has already been added
	    }
	  else
	    {
	      CU_ASSERT(i==element_get_key(i, hash_table->kmer_size)); 
	    }
	}
      
      Element* e=NULL;
      
      
      for (i=0; i< max_kmer_to_test; i=i+step)
	{
	  e = hash_table_find(element_get_key(i,hash_table->kmer_size), hash_table);
	  CU_ASSERT(e!=NULL);
	  if (e !=NULL)
	    {
	      CU_ASSERT(element_get_key(e->kmer, hash_table->kmer_size)==element_get_key(i, hash_table->kmer_size));
	    }
	  else
	    {
	      printf("e is NULL for i=%lld - unable to find\n",i);
	      //exit(1);
	    }
	}
      
      hash_table_free(&hash_table);
      CU_ASSERT(hash_table == NULL);
      
    }

  int num_bits;
  int bucket;
  int max_tries;
  int step;

  for (kmer_size=3; kmer_size<10; kmer_size+=2)
    {
      num_bits=15;
      bucket=100;
      max_tries=10;
      step=1;
      max_key_given_kmer_size = (BinaryKmer) 1 <<(kmer_size*2);
      test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size); 
    }


  /*  
  kmer_size=17;
  num_bits=16;
  bucket=400;
  max_tries=20;
  step=10000000;
  max_key_given_kmer_size = (BinaryKmer) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);
  */

  //now test with just one bucket
  kmer_size=3;
  num_bits=20;
  bucket=1;
  max_tries=10;
  step=1;
  max_key_given_kmer_size = (BinaryKmer) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);


  //now test with many very shallow buckets
  kmer_size=3;
  num_bits=2;
  bucket=10000000;
  max_tries=10;
  step=1;
  max_key_given_kmer_size = (BinaryKmer) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);


  /* this takes ages:  
  //and with large kmer but still shallow buckets
  kmer_size=31;
  step=200000;
  long long max_kmer_to_test=100000;
  test(kmer_size, num_bits, bucket, max_tries, step, max_kmer_to_test);
  */
}


void test_hash_table_apply_or_insert()
{
  short kmer_size;
  long long max_key_given_kmer_size;
  

  //adds kmers to hash, then goes through and tries applying to them all, and checks if the apply was applied
  void test(short kmer_size, int num_bits, int bucket, int max_tries, int step, long long max_kmer_to_test)
    {
      
      int number_of_bits      = num_bits;
      int bucket_size         = bucket;
      long long bad_reads     = 0; 
      int max_retries         = max_tries;
      boolean found           = false;
      
      HashTable* hash_table  = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);
      
      
      long long i;
      
      //insert all possible binary kmers, prior to start of test
      for (i=0; i< max_kmer_to_test; i++)
	{
	  hash_table_find_or_insert(element_get_key(i, hash_table->kmer_size),&found, hash_table);
	}
      
      
      boolean applied=false;
      
      //now test the "apply"
      for (i=0; i< max_kmer_to_test; i=i+step)
	{
	  applied = hash_table_apply_or_insert(element_get_key(i,hash_table->kmer_size),&db_node_action_set_status_pruned , hash_table);
	  CU_ASSERT(applied==true);
	  CU_ASSERT(db_node_check_status(hash_table_find(element_get_key(i,hash_table->kmer_size),hash_table), pruned)==true);
	}
      
      hash_table_free(&hash_table);
      CU_ASSERT(hash_table == NULL);
      
    }

  int num_bits;
  int bucket;
  int max_tries;
  int step;

  for (kmer_size=3; kmer_size<10; kmer_size+=2)
    {
      max_key_given_kmer_size = (BinaryKmer) 1 <<(kmer_size*2);
      num_bits=15;
      bucket=100;
      max_tries=10;
      step=1;
      test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);
    }

  //now test with just one bucket
  kmer_size=3;
  num_bits=20;
  bucket=1;
  max_tries=10;
  step=1;
  max_key_given_kmer_size = (BinaryKmer) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);


  //now test with many very shallow buckets
  kmer_size=3;
  num_bits=2;
  bucket=10000000;
  max_tries=10;
  step=1;
  max_key_given_kmer_size = (BinaryKmer) 1 <<(kmer_size*2);
  test(kmer_size, num_bits, bucket, max_tries, step, max_key_given_kmer_size);

  //and with large kmer but still shallow buckets
  kmer_size=31;
  step=200000;
  long long max_kmer_to_test=100000;
  test(kmer_size, num_bits, bucket, max_tries, step, max_kmer_to_test);


  
}
