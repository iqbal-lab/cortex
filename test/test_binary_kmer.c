/*
 * test_binary_kmer.c
 *
 */

#include <CUnit.h>
#include <Basic.h>
#include <binary_kmer.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "test_binary_kmer.h"


void test_seq_to_binary_kmer_and_binary_kmer_to_seq(){
  
  BinaryKmer kmer;
  char seq[7];
  
  kmer = seq_to_binary_kmer("ATCGCGC",7);
 
  CU_ASSERT_STRING_EQUAL("ATCGCGC",binary_kmer_to_seq(kmer,7,seq));
  
}

void test_binary_kmer_reverse_complement(){
  
  BinaryKmer kmer, kmer_reverse;
  char seq7[7];
  char seq1[1];
  
  kmer = seq_to_binary_kmer("ATCGCGC",7);
  kmer_reverse =  binary_kmer_reverse_complement(kmer, 7);
  
  CU_ASSERT_STRING_EQUAL("GCGCGAT",binary_kmer_to_seq(kmer_reverse,7,seq7));

  kmer = seq_to_binary_kmer("A",1);
  kmer_reverse =  binary_kmer_reverse_complement(kmer, 1);
  
  CU_ASSERT_STRING_EQUAL("T",binary_kmer_to_seq(kmer_reverse,1,seq1));
  
}



void test_get_sliding_windows_from_sequence(){
  
   //----------------------------------
    KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
    if (windows == NULL){
      fputs("Out of memory trying to allocate a KmerSlidingWindowSet",stderr);
      exit(1);
    } 
    
    //allocate memory for the sliding windows         
    windows->window = malloc(sizeof(KmerSlidingWindow) * 20);       
    if (windows->window== NULL){
      fputs("binary_kmer: Out of memory trying to allocate an array of KmerSlidingWindow",stderr);
      exit(1);
    }

    //allocate memory for every every sliding window
    int w;
    for(w=0;w<20;w++){
      KmerSlidingWindow * current_window =&(windows->window[w]);

      current_window->kmer = malloc(sizeof(BinaryKmer) * 20);
      if (current_window->kmer == NULL){
	fputs("binary_kmer: Out of memory trying to allocate an array of BinaryKmer",stderr);
	exit(1);
      }      
    }      
    //----------------------------------
  
    char * seq = "AAAAANNTTTTGGGG";
    char kmer_seq3[3];

    int nkmers1 = get_sliding_windows_from_sequence(seq,NULL,strlen(seq),0,3,windows,20,20);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 3);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[0].kmer[0],3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[0].kmer[1],3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[0].kmer[2],3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 6);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[0],3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[1],3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[2],3,kmer_seq3),"TTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[3],3,kmer_seq3),"TGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[4],3,kmer_seq3),"GGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[5],3,kmer_seq3),"GGG");

    CU_ASSERT_EQUAL(nkmers1, 9);
    
    char qual1[15] = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 10, 20, 20, 20};

    int nkmers2 = get_sliding_windows_from_sequence(seq,qual1,strlen(seq),15,3,windows,20,20);

    CU_ASSERT_EQUAL(windows->nwindows,3);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 3);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[0].kmer[0],3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[0].kmer[1],3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[0].kmer[2],3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 2);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[0],3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[1],3,kmer_seq3),"TTT");

    CU_ASSERT_EQUAL((windows->window[2]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[2].kmer[0],3,kmer_seq3),"GGG");

    CU_ASSERT_EQUAL(nkmers2, 6);

    char qual2[15] = { 5, 10, 20, 20, 20, 20, 0, 0, 0, 20, 20, 30, 20, 20, 20};

    int nkmers3 = get_sliding_windows_from_sequence(seq,qual2,strlen(seq),15,3,windows,20,20);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[0].kmer[0],3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[0],3,kmer_seq3),"TTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[1],3,kmer_seq3),"TGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[2],3,kmer_seq3),"GGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[3],3,kmer_seq3),"GGG");
    
    CU_ASSERT_EQUAL(nkmers3, 5);

    char qual3[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int nkmers4 = get_sliding_windows_from_sequence(seq,qual3,strlen(seq),15,3,windows,20,20);

    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers4, 0);

    char kmer_seq5[5];
    int nkmers5 = get_sliding_windows_from_sequence(seq,NULL,strlen(seq),0,5,windows,20,20);
    
    CU_ASSERT_EQUAL(windows->nwindows,2);
  
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[0].kmer[0],5,kmer_seq5),"AAAAA");

  
    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[0],5,kmer_seq5),"TTTTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[1],5,kmer_seq5),"TTTGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[2],5,kmer_seq5),"TTGGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(windows->window[1].kmer[3],5,kmer_seq5),"TGGGG");

    CU_ASSERT_EQUAL(nkmers5, 5);

    
  }
