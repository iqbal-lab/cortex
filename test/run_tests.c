#include "test_binary_kmer.h"
#include "test_count_kmers.h"
#include "test_dB_graph.h"
#include "test_dB_graph_node.h"
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

  if (NULL == CU_add_test(pSuite, "test conversion from char to binary nucleotide", test_seq_to_binary_kmer)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pSuite, "test conversion from binary nucleotide to char", test_binary_kmer_to_seq)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pSuite, "test creation and deletion of binary kmers", test_binary_kmer_creation_and_deletion)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  //  if (NULL == CU_add_test(pSuite, "test kmer counting ", test_count_kmers)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
  if (NULL == CU_add_test(pSuite, "test supernode walking", test_supernode_walking)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  // if (NULL == CU_add_test(pSuite, "test CRC returns positive values", test_element_hashval)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
  //  if (NULL == CU_add_test(pSuite, "test element_equal works", test_element_equal)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}





