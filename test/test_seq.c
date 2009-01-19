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

void test_read_sequence_from_fasta()
{

  // 1. Read from fasta with one short and one long read

  Sequence* seq;
  //  int  seq_length=0;

  FILE* fp = fopen("../test/test_seq_dir/two_entry_oneofwhich_very_long.fasta", "r");
  seq=read_sequence_from_fasta(fp);
  CU_ASSERT_STRING_EQUAL(seq->seq, "ACGT");
  free_sequence(&seq);

  seq=read_sequence_from_fasta(fp);
  CU_ASSERT_STRING_EQUAL(seq->seq, "ATAATTTCACTATTTAAGGTGTAAAGAAGTACTTTAGTCAAAGTAGGAGAGTGACAATAAAGTAAAAGTAAAATAGAAAGGTTATACAGGTGTTATACCTGGAGTGGACGTAGTTAGGGTTTACCGTTTTTAAAAAACGGGGTAAAAGGGTTTTAAAAACTTAAAACCCGTAAAACGTTTAAACGGGGGAAAAGGGACCGTTAGGGGGTAAGGGGGGTAAGGGGGGGGGGTTTGGGTTTAGGGGTTTAGGGGGGTAAAAAGGGGGGGGGGTAAAAAAAAACGTTTAAAAAAAAAGGGAAAACCGAAAAAAGGGTAAAACGGGTAAAAAACCCGAAAACCCGTAAAACCCCCCCGTTTTACCGGGGTTAACCCCCCCCGTTTTTAAAACCCCCCCCCGTTTTTACCCCCCGTTTTTTAACCCCCCCCCCCCTTTTTTTTAAAAAAAACCCCCCGTTTTTTAAAAAAACCCCCCGTTTAAAAAAAACGGGTTTTTTAAAAAACGGGGGGGGGGTTTTTTACCGGGGGTTTTTTAACCCCCCCCCCGGGTTTTTTACGGGGGGGGGGTTTTTTAACGGGGGGGGGGTTTTAAAAAACGGGGGTTTTAAAAAACGGTTTTTTAAGGGGGTTTTTAAAAAGGGTTTAAAACTTTTTTTTTTTTTAACCCCGTTTAAAAAAACGGGGTTTTTTTTAAAAAAAAAAAAAAAAACGGGGTTAAAAAAAAAAAAAAAAAAAAAAAAAAACCGTTTTAACCCCGGTTAAAAACGGGGGTTTTTTTTAAACCCGGGGGGGGTTAAAAAAAACCCCCCCGGGGGGTTTAAAAAAAAAACGGGGGGGGGGTTTTTTTTTTAAAGGGGTTTAAAAAAACCGTTTTTTTAAACCCCCCGTTTAAAAAAAAAACCCCCCCCCCCGGTTTTTTTTTTTTAACCCCCCGGGGGGTTTTTTTAAACCCCCCCCGTTAAAAAAACCCCCCGGGGTTAACCCCCCCGTTAAAAGGTTTTTAACGGGGGTTTTAAAACGGGGGG");
  free_sequence(&seq);
  fclose(fp);


}

void test_read_sequence_from_fasta_efficient()
{

  // 1. Read from fasta with one short and one long read

  int longest_expected_read_length = 3000;

  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL)
    {
      printf("Out of memory allocating a Sequence object");
      exit(1);
    }

  seq->seq = malloc(sizeof(char)*longest_expected_read_length ) ;
  if (seq->seq == NULL)
    {
      printf("Out of memory allocating a seq->seq object");
      exit(1);
    }

  seq->max=longest_expected_read_length;
  seq->qual=NULL;
  seq->length=0;

  int did_realloc_happen=0;

  FILE* fp = fopen("../test/test_seq_dir/two_entry_oneofwhich_very_long.fasta", "r");
  int result =read_sequence_from_fasta_efficient(fp, seq, &did_realloc_happen);
  CU_ASSERT_STRING_EQUAL(seq->seq, "ACGT");
  CU_ASSERT(did_realloc_happen==0);

  result=read_sequence_from_fasta_efficient(fp, seq, &did_realloc_happen);
  CU_ASSERT_STRING_EQUAL(seq->seq, "ATAATTTCACTATTTAAGGTGTAAAGAAGTACTTTAGTCAAAGTAGGAGAGTGACAATAAAGTAAAAGTAAAATAGAAAGGTTATACAGGTGTTATACCTGGAGTGGACGTAGTTAGGGTTTACCGTTTTTAAAAAACGGGGTAAAAGGGTTTTAAAAACTTAAAACCCGTAAAACGTTTAAACGGGGGAAAAGGGACCGTTAGGGGGTAAGGGGGGTAAGGGGGGGGGGTTTGGGTTTAGGGGTTTAGGGGGGTAAAAAGGGGGGGGGGTAAAAAAAAACGTTTAAAAAAAAAGGGAAAACCGAAAAAAGGGTAAAACGGGTAAAAAACCCGAAAACCCGTAAAACCCCCCCGTTTTACCGGGGTTAACCCCCCCCGTTTTTAAAACCCCCCCCCGTTTTTACCCCCCGTTTTTTAACCCCCCCCCCCCTTTTTTTTAAAAAAAACCCCCCGTTTTTTAAAAAAACCCCCCGTTTAAAAAAAACGGGTTTTTTAAAAAACGGGGGGGGGGTTTTTTACCGGGGGTTTTTTAACCCCCCCCCCGGGTTTTTTACGGGGGGGGGGTTTTTTAACGGGGGGGGGGTTTTAAAAAACGGGGGTTTTAAAAAACGGTTTTTTAAGGGGGTTTTTAAAAAGGGTTTAAAACTTTTTTTTTTTTTAACCCCGTTTAAAAAAACGGGGTTTTTTTTAAAAAAAAAAAAAAAAACGGGGTTAAAAAAAAAAAAAAAAAAAAAAAAAAACCGTTTTAACCCCGGTTAAAAACGGGGGTTTTTTTTAAACCCGGGGGGGGTTAAAAAAAACCCCCCCGGGGGGTTTAAAAAAAAAACGGGGGGGGGGTTTTTTTTTTAAAGGGGTTTAAAAAAACCGTTTTTTTAAACCCCCCGTTTAAAAAAAAAACCCCCCCCCCCGGTTTTTTTTTTTTAACCCCCCGGGGGGTTTTTTTAAACCCCCCCCGTTAAAAAAACCCCCCGGGGTTAACCCCCCCGTTAAAAGGTTTTTAACGGGGGTTTTAAAACGGGGGG");
  CU_ASSERT(did_realloc_happen==0);

  free_sequence(&seq);
  fclose(fp);


}



void test_read_sequence_from_fastq()
{

  // 1. Read from simple fastq:
  //     @Zam
  //     ACGT
  //     +
  //     &&&&


  Sequence* seq;
  //  int seq_length=0;

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
