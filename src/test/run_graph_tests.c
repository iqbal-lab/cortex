#include "test_binary_kmer.h"
#include "test_count_kmers.h"
#include "test_hash_table.h"
#include "test_dB_graph.h"
#include "test_dB_graph_node.h"
#include "test_seq.h"
#include <CUnit.h>
#include <Basic.h>

int  main()
{

  CU_pSuite pGraphSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();

  /* add a suite to the registry */
  pGraphSuite = CU_add_suite("Test graph", NULL, NULL);
  if (NULL == pGraphSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suites */

  if (NULL == CU_add_test(pGraphSuite, "test conversion from char to binary nucleotide", test_seq_to_binary_kmer)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pGraphSuite, "test conversion from binary nucleotide to char", test_binary_kmer_to_seq)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pGraphSuite, "test creation and deletion of binary kmers", test_binary_kmer_creation_and_deletion)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  //  if (NULL == CU_add_test(pGraphSuite, "test kmer counting ", test_count_kmers)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
   // if (NULL == CU_add_test(pGraphSuite, "test CRC returns positive values", test_element_hashval)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
  //  if (NULL == CU_add_test(pGraphSuite, "test element_equal works", test_element_equal)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}


  if (NULL == CU_add_test(pGraphSuite, "test read_sequence_from_fasta with simple file with two reads, one of which is very long", test_read_sequence_from_fasta)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pGraphSuite, "test read_sequence_from_fasta_efficient, checking for unnecessary realloc", test_read_sequence_from_fasta_efficient)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pGraphSuite, "test read_sequence_from_fastq", test_read_sequence_from_fastq)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pGraphSuite, "Load a simple fasta and confirm that hash_table_find will find the kmers and reverse complements of the kmers in those reads, and will not find any others.", test_hash_table_find)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pGraphSuite, "test supernode walking", test_supernode_walking)) {
    CU_cleanup_registry();
    return CU_get_error();
  }


  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}





