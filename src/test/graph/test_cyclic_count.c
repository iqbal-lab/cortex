#include <CUnit.h>
#include <Basic.h>
#include <cyclic_count.h>
#include <dB_graph.h>
#include <element.h>
#include <binary_kmer.h>
#include <seq.h>
#include <test_cyclic_count.h>


void test_rotate()
{
  char* seq1 = "ACG";
  char* seq1_r1 = "GAC";
  char* seq1_r2 = "CGA";


  char* seq2 = "AAAAAAAC";
  char* seq2_r1 = "CAAAAAAA";
  char* seq2_r2 = "ACAAAAAA";
  char* seq2_r3 = "AACAAAAA";
  char* seq2_r4 = "AAACAAAA";
  char* seq2_r5 = "AAAACAAA";
  char* seq2_r6 = "AAAAACAA";
  char* seq2_r7 = "AAAAAACA";
  char* seq2_r8 = "AAAAAAAC";

  char* seq3 =    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC";
  char* seq3_r1 = "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";


  BinaryKmer bk1 = seq_to_binary_kmer(seq1, 3);
  BinaryKmer bk1_r1 = seq_to_binary_kmer(seq1_r1, 3);
  BinaryKmer bk1_r2 = seq_to_binary_kmer(seq1_r2, 3);

  BinaryKmer bk2 = seq_to_binary_kmer(seq2, 8);
  BinaryKmer bk2_r1 = seq_to_binary_kmer(seq2_r1, 8);
  BinaryKmer bk2_r2 = seq_to_binary_kmer(seq2_r2, 8);
  BinaryKmer bk2_r3 = seq_to_binary_kmer(seq2_r3, 8);
  BinaryKmer bk2_r4 = seq_to_binary_kmer(seq2_r4, 8);
  BinaryKmer bk2_r5 = seq_to_binary_kmer(seq2_r5, 8);
  BinaryKmer bk2_r6 = seq_to_binary_kmer(seq2_r6, 8);
  BinaryKmer bk2_r7 = seq_to_binary_kmer(seq2_r7, 8);
  BinaryKmer bk2_r8 = seq_to_binary_kmer(seq2_r8, 8);

  BinaryKmer bk3 = seq_to_binary_kmer(seq3, 31);
  BinaryKmer bk3_r1 = seq_to_binary_kmer(seq3_r1, 31);


  BinaryKmer tmp;

  tmp = rotate_least_sig_2k_bits(bk1,3);
  CU_ASSERT(tmp == bk1_r1);
  tmp = rotate_least_sig_2k_bits(bk1_r1,3);
  CU_ASSERT(tmp == bk1_r2);
  tmp = rotate_least_sig_2k_bits(bk1_r2,3);
  CU_ASSERT(tmp == bk1);



  tmp = rotate_least_sig_2k_bits(bk2,8);
  CU_ASSERT(tmp == bk2_r1);
  tmp = rotate_least_sig_2k_bits(bk2_r1,8);
  CU_ASSERT(tmp == bk2_r2);
  tmp = rotate_least_sig_2k_bits(bk2_r2,8);
  CU_ASSERT(tmp == bk2_r3);
  tmp = rotate_least_sig_2k_bits(bk2_r3,8);
  CU_ASSERT(tmp == bk2_r4);
  tmp = rotate_least_sig_2k_bits(bk2_r4,8);
  CU_ASSERT(tmp == bk2_r5);
  tmp = rotate_least_sig_2k_bits(bk2_r5,8);
  CU_ASSERT(tmp == bk2_r6);
  tmp = rotate_least_sig_2k_bits(bk2_r6,8);
  CU_ASSERT(tmp == bk2_r7);
  tmp = rotate_least_sig_2k_bits(bk2_r7,8);
  CU_ASSERT(tmp == bk2_r8);


  tmp = rotate_least_sig_2k_bits(bk3,8);
  CU_ASSERT(tmp == bk3_r1);


  
  

}


void test_count_how_many_cyclic_perms_of_this_node_are_in_graph()
{

}
