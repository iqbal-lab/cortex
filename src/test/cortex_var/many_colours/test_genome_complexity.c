#include <stdio.h>
#include <time.h>
#include <CUnit.h>
#include <genome_complexity.h>
#include <limits.h>

void test_count_reads_where_snp_makes_clean_bubble1()
{

  //  simple test. Use this refernece genome which I load into colour 0:

  /// file: data/test/genome_complexity/test_allele_clean_file1.fa

  //  >ref contains CAAGTTC CACGTTC and CAGGTTC
  //  CAAGTTCAGAGTTACTCACACCCGATCGATAAGCGGTACAGAGCACGTTCAGAAAAAAACAGGTTCAGA
  //  >ref + SNP
  //  TATCCATGTTCAGAGTTACTGACACCCGATCGATAAGCG


  //first set up the hash/graph
  int kmer_size = 7;
  int number_of_bits = 8;
  int bucket_size    = 10;
  long long bad_reads = 0; 
  int max_retries=10;

  dBGraph * hash_table = hash_table_new(number_of_bits,bucket_size,max_retries,kmer_size);

  if (hash_table==NULL)
    {
      printf("unable to alloc the hash table. dead before we even started. OOM");
      exit(1);
    }

  long long seq_read=0;
  long long seq_loaded=0;

  load_population_as_fasta("../data/test/genome_complexity/pop_first_test", &seq_read, &seq_loaded, &bad_reads, hash_table, NULL);



  //and use this file of reads: /data/test/genome_complexity/test_allele_clean_file2.fa
  //  >read lies entirely in graph defined by test_allele_clean_file1.fa, and is clean (forms supernode at k=7)
  //  GTTCAGAGTTACT
  //  >read lies in graph defined by test_allele_clean_file1.fa, but overlaps a junction so is not clean
  //  TTACTGACACCCGATCG
  //  >read lies in ref but any change to the 7th base results in an overlap with the ref, so
  //  GTTCAGAGTTACTCACA

  int col_genome=0;
  int reads_tested=0;
  int reads_where_snp_makes_clean_bubble = 0;
  dBNode* array_nodes[500];
  Orientation array_or[500];
  
    //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  int max_read_length=2000;
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in test. Exit.\n");
      exit(1);
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-hash_table->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in test. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 
	  
	  
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = hash_table->kmer_size;
      //printf("new_entry must be true in hsi test function");
      //exit(1);
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  count_reads_where_snp_makes_clean_bubble(hash_table, "../data/test/genome_complexity/test_allele_clean_file2.fa", true,
					   col_genome, &reads_tested, &reads_where_snp_makes_clean_bubble, 
					   file_reader_fasta,
					   array_nodes, array_or, seq, kmer_window, max_read_length);


  
  CU_ASSERT(reads_tested==3);
  CU_ASSERT(reads_where_snp_makes_clean_bubble==2);

    
}
