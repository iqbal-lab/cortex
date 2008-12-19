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
#include <assert.h>
#include "test_binary_kmer.h"

void test_seq_to_bin()
{

  CU_ASSERT_EQUAL(seq_to_bin("A",1),Adenine);
  CU_ASSERT_EQUAL(seq_to_bin("C",1),Cytosine);
  CU_ASSERT_EQUAL(seq_to_bin("G",1),Guanine);
  CU_ASSERT_EQUAL(seq_to_bin("T",1),Thymine);
  CU_ASSERT_EQUAL(seq_to_bin("a",1),Adenine);
  CU_ASSERT_EQUAL(seq_to_bin("c",1),Cytosine);
  CU_ASSERT_EQUAL(seq_to_bin("g",1),Guanine);
  CU_ASSERT_EQUAL(seq_to_bin("t",1),Thymine);
  CU_ASSERT_EQUAL(seq_to_bin("X",1),Undefined);
  CU_ASSERT_EQUAL(seq_to_bin("N",1),Undefined);



  //if behaving properly, this function call should lead to an assertion, which will be caught by CU_PASS
  // i.e. we expect it to hit the assert, ad that counts as a pass.
  //CU_PASS(seq_to_bin('N'));
  //CU_PASS(seq_to_bin('X'));

}

void test_bin_to_seq()
{
  //fors soem reason i get a "bus error" if I use CU_ASSERT_STRING_EQUAL
  CU_ASSERT_STRING_EQUAL("A",bin_to_seq(Adenine,1));
  CU_ASSERT_STRING_EQUAL("C",bin_to_seq(Cytosine,1));
  CU_ASSERT_STRING_EQUAL("G",bin_to_seq(Guanine,1));
  CU_ASSERT_STRING_EQUAL("T",bin_to_seq(Thymine,1));

  CU_PASS(bin_to_seq(100));
}


void test_binary_kmer_creation_and_deletion()
{
  //create a Kmers object
  char* seq = "AAAA";
  Kmers* k = get_binary_kmers(seq,4,2);
  if (k==NULL)
    {
      assert(0);
    }
  free_kmers(&k);
  
}

