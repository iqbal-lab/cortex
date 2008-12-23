/*
 * test_binary_kmer.c
 *
 *  Created on: Sep 5, 2008
 *      Author: zi
 */

#include <CUnit.h>
#include <Basic.h>
#include <binary_kmer.h>
#include <stdio.h>
#include<stdlib.h>
#include <assert.h>
#include "test_binary_kmer.h"

void test_seq_to_binary_kmer()
{

  CU_ASSERT_EQUAL(seq_to_binary_kmer("A",1),Adenine);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("C",1),Cytosine);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("G",1),Guanine);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("T",1),Thymine);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("a",1),Adenine);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("c",1),Cytosine);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("g",1),Guanine);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("t",1),Thymine);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("X",1),Undefined);
  CU_ASSERT_EQUAL(seq_to_binary_kmer("N",1),Undefined);



  //if behaving properly, this function call should lead to an assertion, which will be caught by CU_PASS
  // i.e. we expect it to hit the assert, ad that counts as a pass.
  //CU_PASS(seq_to_binary_kmer('N'));
  //CU_PASS(seq_to_binary_kmer('X'));

}

void test_binary_kmer_to_seq()
{
  //fors soem reason i get a "bus error" if I use CU_ASSERT_STRING_EQUAL
  char* test_A = binary_kmer_to_seq(Adenine,1);
  char* test_C = binary_kmer_to_seq(Cytosine,1);
  char* test_G = binary_kmer_to_seq(Guanine,1);
  char* test_T = binary_kmer_to_seq(Thymine,1);

  CU_ASSERT_STRING_EQUAL("A",test_A);
  CU_ASSERT_STRING_EQUAL("C",test_C);
  CU_ASSERT_STRING_EQUAL("G",test_G);
  CU_ASSERT_STRING_EQUAL("T",test_T);

  free(test_A);
  free(test_C);
  free(test_G);
  free(test_T);
  
  //CU_PASS(binary_kmer_to_seq(100));
}


void test_binary_kmer_creation_and_deletion()
{
  //create a KmerArray object
  char* seq = "AAAA";
  KmerArray* k = get_binary_kmers_from_sequence(seq,4,2);
  if (k==NULL)
    {
      assert(0);
    }
  binary_kmer_free_kmers(&k);
  
}

