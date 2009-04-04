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


void test_read_sequence_from_fasta_when_file_has_long_and_bad_reads()
{
  Sequence * seq = malloc(sizeof(Sequence));
 
  if (seq == NULL){							
    fputs("Out of memory trying to allocate Sequence\n",stderr);	
    exit(1);								
  }
  //pre-allocate space where to read the sequences
  int max_read_length=50;
  alloc_sequence(seq,50,LINE_MAX);

  int length_seq;
  FILE* fp1 = fopen("../data/test/basic/includes_two_reads_that_are_too_long.fasta", "r");

  // >read1 100 bases
  // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  // >read2 5 bases
  // CCCCC
  // >read3 110 bases
  // GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
  // >read4 4 bases
  // TTTT

  length_seq = read_sequence_from_fasta(fp1,seq,max_read_length);
  
  CU_ASSERT_EQUAL(length_seq, 5);
  CU_ASSERT_STRING_EQUAL("read2",seq->name);
  CU_ASSERT_STRING_EQUAL("CCCCC",seq->seq);

  length_seq = read_sequence_from_fasta(fp1,seq,max_read_length);
  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read4",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTT",seq->seq);



  fclose(fp1);

  FILE* fp2= fopen("../data/test/basic/includes_reads_that_have_bad_characters.fasta", "r");

  // >read1
  // AAAAAAAAAAAA¢
  // >read2
  // ¡€#¢∞§¶#¶•#•#•#ª#ª#ª#ªº#º#º#º––––
  // >read3 4 c's
  // CCCC
  // >read4 10 Ts
  // TTTTTTTTTT
  // >read5
  // $
  // >read6
  // AAAAAAAAAAAAAAAAAA#A
  // >read7
  // AAA

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length);
  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT_STRING_EQUAL("CCCC",seq->seq);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length);
  
  CU_ASSERT_EQUAL(length_seq, 10);
  CU_ASSERT_STRING_EQUAL("read4",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTTTTTTTT",seq->seq);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length);
  
  CU_ASSERT_EQUAL(length_seq, 3);
  CU_ASSERT_STRING_EQUAL("read7",seq->name);
  CU_ASSERT_STRING_EQUAL("AAA",seq->seq);

  length_seq = read_sequence_from_fasta(fp2,seq,max_read_length);
  
  CU_ASSERT_EQUAL(length_seq, 0);
  fclose(fp2);


  //now make sure we do not get trapped in an infinite loop if the last read of a file is bad

  FILE* fp3= fopen("../data/test/basic/includes_final_read_that_has_bad_characters.fasta", "r");

  // >read1
  // AAAAAAAAAAAA¢
  // >read2
  // ¡€#¢∞§¶#¶•#•#•#ª#ª#ª#ªº#º#º#º––––
  // >read3 4 c's
  // CCCC
  // >read4 10 Ts
  // TTTTTTTTTT
  // >read5
  // $
  // >read6
  // AAAAAAAAAAAAAAAAAA#A


  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length);

  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT_STRING_EQUAL("CCCC",seq->seq);

  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length);
  
  CU_ASSERT_EQUAL(length_seq, 10);
  CU_ASSERT_STRING_EQUAL("read4",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTTTTTTTT",seq->seq);

  length_seq = read_sequence_from_fasta(fp3,seq,max_read_length);
  
  CU_ASSERT_EQUAL(length_seq, 0);
  fclose(fp3);



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
  // @@@$%^&*
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
  CU_ASSERT_STRING_EQUAL("@@@$%^&*",seq->qual);

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




void test_read_sequence_from_fastq_with_bad_reads_and_long_reads()
{

  //pre-allocate space where to read the sequences
  Sequence* seq = malloc(sizeof(Sequence));
  if (seq==NULL){
    fputs("Out of memory trying to allocate a Sequence",stderr);
      exit(1);
  } 

  int max_read_length=50;
  alloc_sequence(seq,max_read_length,LINE_MAX);
  

  
  int length_seq;
  /*
  FILE* fp1 = fopen("../data/test/basic/includes_one_read_that_is_too_long.fastq", "r");

  // @read1
  // ACGT
  // +
  // @@@@
  // @read2
  // CCCC
  // +
  // 5555
  // @read3
  // AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  // -
  // 0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
  // @read4
  // ACGT
  // +
  // 3333



  length_seq = read_sequence_from_fastq(fp1,seq,max_read_length);  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  CU_ASSERT_STRING_EQUAL("@@@@",seq->qual);

  length_seq = read_sequence_from_fastq(fp1,seq,max_read_length);  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read2",seq->name);
  CU_ASSERT_STRING_EQUAL("CCCC",seq->seq);
  CU_ASSERT_STRING_EQUAL("5555",seq->qual);

  length_seq = read_sequence_from_fastq(fp1,seq,max_read_length);  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read4",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGT",seq->seq);
  CU_ASSERT_STRING_EQUAL("3333",seq->qual);


  fclose(fp1);

  */


  FILE* fp2 = fopen("../data/test/basic/includes_reads_with_bad_characters.fastq", "r");

  //@read1
  //ACGTACGTACGTACGT
  //+
  //WEW£WEW£WEW£WEWA
  //@read2
  //AAAA#¢A
  //+
  //@@@@@@@
  //@read3
  //TTTT
  //+
  //3333


  length_seq = read_sequence_from_fastq(fp2,seq,max_read_length);  
  CU_ASSERT_EQUAL(length_seq, 16);
  CU_ASSERT_STRING_EQUAL("read1",seq->name);
  CU_ASSERT_STRING_EQUAL("ACGTACGTACGTACGT",seq->seq);
  CU_ASSERT_STRING_EQUAL("WEW£WEW£WEW£WEWA",seq->qual);

  length_seq = read_sequence_from_fastq(fp2,seq,max_read_length);  
  CU_ASSERT_EQUAL(length_seq, 4);
  CU_ASSERT_STRING_EQUAL("read3",seq->name);
  CU_ASSERT_STRING_EQUAL("TTTT",seq->seq);
  CU_ASSERT_STRING_EQUAL("3333",seq->qual);

  length_seq = read_sequence_from_fastq(fp2,seq,max_read_length);  
  CU_ASSERT_EQUAL(length_seq, 0);
  
  fclose(fp2);


  free_sequence(&seq);


}
