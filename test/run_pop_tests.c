#include "test_binary_kmer.h"
#include "test_count_kmers.h"
#include "test_dB_graph.h"
#include "test_dB_graph_node.h"
#include "test_seq.h"
#include "test_pop_load_and_print.h"
#include <CUnit.h>
#include <Basic.h>

int  main()
{

  CU_pSuite pPopGraphSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();

  /* add a suite to the registry */
  pPopGraphSuite = CU_add_suite("Test pop graph", NULL, NULL);
  if (NULL == pPopGraphSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suites */

  if (NULL == CU_add_test(pPopGraphSuite, "test conversion from char to binary nucleotide", test_seq_to_binary_kmer)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "test conversion from binary nucleotide to char", test_binary_kmer_to_seq)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "test creation and deletion of binary kmers", test_binary_kmer_creation_and_deletion)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  //  if (NULL == CU_add_test(pPopGraphSuite, "test kmer counting ", test_count_kmers)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
  //if (NULL == CU_add_test(pPopGraphSuite, "test supernode walking", test_supernode_walking)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
  // if (NULL == CU_add_test(pPopGraphSuite, "test CRC returns positive values", test_element_hashval)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
  //  if (NULL == CU_add_test(pPopGraphSuite, "test element_equal works", test_element_equal)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
  //if (NULL == CU_add_test(pPopGraphSuite, "test read_sequence_from_fastq", test_read_sequence_from_fastq)) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}
  if (NULL == CU_add_test(pPopGraphSuite, "test loading two people, one fasta each, in same population, and print supernodes", test_load_two_people_in_same_populations_and_print_separately_their_supernodes)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}

