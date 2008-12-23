/*
 * test_seq.c
 *
 *  Created on: Sep 5, 2008
 *      Author: zi
 */

#include <CUnit.h>
#include <Basic.h>
#include <seq.h>
#include <stdio.h>
#include<stdlib.h>
#include <assert.h>
#include "test_seq.h"

void test_read_sequence_from_fastq()
{

  // 1. Read from simple fastq:
  //     @Zam
  //     ACGT
  //     +
  //     &&&&


  Sequence* seq;
  int seq_length=0;

  FILE* fp = fopen("../test/test_seq_dir/one_entry.fastq", "r");
  seq=read_sequence_from_fastq(fp);
  
  CU_ASSERT_STRING_EQUAL(seq->seq, "ACGT");
  CU_ASSERT_STRING_EQUAL(seq->qual, "&&&&");

  free_sequence(&seq);
  fclose(fp);

// 2. Read from longer fastq:
  //     @Zam
  //     ACGT
  //     +
  //     &&&&
  //     @iqbal
  //     AAAAAAAA
  //     +
  //     !@£$%^&*
  //     @zamzamzam
  //     ATATATAT
  //     -
  //     @@@@@@@@

 
  fp = fopen("../test/test_seq_dir/three_entries.fastq", "r");
 
 
  seq=read_sequence_from_fastq(fp);
  CU_ASSERT_STRING_EQUAL(seq->seq, "ACGT");
  CU_ASSERT_STRING_EQUAL(seq->qual, "&&&&");
  free_sequence(&seq);
 
  seq=read_sequence_from_fastq(fp);
  CU_ASSERT_STRING_EQUAL(seq->seq, "AAAAAAAA");
  CU_ASSERT_STRING_EQUAL(seq->qual, "!@@$%^&*");
  free_sequence(&seq);

 
  seq=read_sequence_from_fastq(fp);
  CU_ASSERT_STRING_EQUAL(seq->seq, "ATATATAT");
  CU_ASSERT_STRING_EQUAL(seq->qual, "@@@@@@@@");
  free_sequence(&seq);



  fclose(fp);




}
