/*
 * test_seq.c
 *
 */

#include <CUnit.h>
#include <Basic.h>
#include <seq.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <test_seq.h>

void test_read_sequence_from_fasta(){

  Sequence * seq = malloc(sizeof(Sequence));
 
  if (seq == NULL){							
    fputs("Out of memory trying to allocate Sequence\n",stderr);	
    exit(1);								
  }
  //pre-allocate space where to read the sequences
  alloc_sequence(seq,200,LINE_MAX);

  int length_seq;
  FILE* fp1 = fopen("../data/test/basic/one_entry.fasta", "r");

  // 1. Read from simple fasta:
  // >Zam
  // ACGT
  // ACGTACGTACGT
  length_seq = read_sequence_from_fasta(fp1,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 16);
  CU_ASSERT_STRING_EQUAL("Zam",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGT",seq->seq);

  fclose(fp1);

  FILE* fp2 = fopen("../data/test/basic/three_entries.fasta", "r");

  // 2. Read from fasta:
  //>Zam1
  //ACGT
  //ACGTACGTACGT
  //>Zam2
  //ACGT
  //ACGTACGTACGT
  //TTTTTTTT
  //>Zam3
  //ACGTNNACGTACGTACGT

  length_seq = read_sequence_from_fasta(fp2,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 16);
  CU_ASSERT_STRING_EQUAL("Zam1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGT",seq->seq);

  length_seq = read_sequence_from_fasta(fp2,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 24);
  CU_ASSERT_STRING_EQUAL("Zam2",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGTTTTTTTTt",seq->seq);

  length_seq = read_sequence_from_fasta(fp2,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 18);
  CU_ASSERT_STRING_EQUAL("Zam3",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGTNNACGTACGTACGT",seq->seq);

  length_seq = read_sequence_from_fasta(fp2,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 0);

  fclose(fp2);
  free_sequence(&seq);

}

void test_read_sequence_from_fastq(){

  //pre-allocate space where to read the sequences
  Sequence* seq = malloc(sizeof(Sequence));
  if (seq==NULL){
    fputs("Out of memory trying to allocate a Sequence",stderr);
      exit(1);
  } 
  
  alloc_sequence(seq,200,LINE_MAX);
  
  int length_seq;
  FILE* fp1 = fopen("../data/test/basic/one_entry.fastq", "r");

  // 1. Read from simple fasta:
  // >Zam
  // ACGT
  // +
  // &&&&

  length_seq = read_sequence_from_fastq(fp1,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("Zam",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  CU_ASSERT_STRING_EQUAL("&&&&",seq->qual);

  FILE* fp2 = fopen("../data/test/basic/three_entries.fastq", "r");
  
  //2. Read from fastq:
  
  // @Zam1
  // ACGT
  // +
  // &&&&
  // @Zam2
  // AAAAAAAA
  // +
  // !@@$%^&*
  // @Zam3
  // ATATATAT
  // TTTTTTTTTT
  // -
  // @@@@@@@&AAAAAABAAA

  length_seq = read_sequence_from_fastq(fp2,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("Zam1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  CU_ASSERT_STRING_EQUAL("&&&&",seq->qual);

  length_seq = read_sequence_from_fastq(fp2,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 8);
  CU_ASSERT_STRING_EQUAL("Zam2",seq->name);
  CU_ASSERT_STRING_EQUAL("AAAAAAAA",seq->seq);
  CU_ASSERT_STRING_EQUAL("!@@$%^&*",seq->qual);

  length_seq = read_sequence_from_fastq(fp2,seq,1000);
  
  CU_ASSERT_EQUAL(length_seq, 18);
  CU_ASSERT_STRING_EQUAL("Zam3",seq->name);
  CU_ASSERT_STRING_EQUAL("ATATATATTTTTTTTTTT",seq->seq);
  CU_ASSERT_STRING_EQUAL("@@@@@@@&AAAAAABAAA",seq->qual);

  length_seq = read_sequence_from_fastq(fp2,seq,1000);

  CU_ASSERT_EQUAL(length_seq, 0);

  fclose(fp2);
  free_sequence(&seq);
}
