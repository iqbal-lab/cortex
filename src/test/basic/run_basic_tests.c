
#include <test_binary_kmer.h>
#include <test_seq.h>
#include <CUnit.h>
#include <Basic.h>

int  main()
{

  CU_pSuite pSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();
  
  /* add a suite to the registry */
  pSuite = CU_add_suite("Suite_1", NULL, NULL);
  if (NULL == pSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suite */

    
  if (NULL == CU_add_test(pSuite, "test reading of fasta file",  test_read_sequence_from_fasta)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "test reading of fasta file when some reads that have bad characters",  test_read_sequence_from_fasta_when_file_has_bad_reads)){ 
    CU_cleanup_registry(); 
     return CU_get_error();
  }

  
   if (NULL == CU_add_test(pSuite, "test reading of fastq file",  test_read_sequence_from_fastq)){ 
     CU_cleanup_registry(); 
     return CU_get_error(); 
   } 
   

  if (NULL == CU_add_test(pSuite, "test reading of long fasta file",  test_read_sequence_from_long_fasta)){
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pSuite, "test shift last kmer to start",  test_shift_last_kmer_to_start_of_sequence)){
    CU_cleanup_registry();
    return CU_get_error();
  }

  
  if (NULL == CU_add_test(pSuite, "test reading of fastq file when some reads are too long or have bad characters",  test_read_sequence_from_fastq_with_bad_reads_and_long_reads)){
    
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "test conversion from binary nucleotide to C string", test_seq_to_binary_kmer_and_binary_kmer_to_seq)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

 if (NULL == CU_add_test(pSuite, "test binary kmer reverse complement", test_binary_kmer_reverse_complement)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

 if (NULL == CU_add_test(pSuite, "test seq reverse complement", test_seq_reverse_complement)) {
   CU_cleanup_registry();
   return CU_get_error();
 }

 if (NULL == CU_add_test(pSuite, "test nucleotide iterator", test_binary_kmer_nucleotide_iterator)) {
    CU_cleanup_registry();
    return CU_get_error();
  }


  if (NULL == CU_add_test(pSuite, "test creation of binary kmers from sequence - sliding window", test_get_sliding_windows_from_sequence)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}





