#include "test_binary_kmer.h"
#include "test_count_kmers.h"
#include "test_dB_graph.h"
#include "test_dB_graph_node.h"
#include "test_seq.h"
#include "test_pop_element.h"
#include "test_pop_load_and_print.h"
#include "test_pop_supernode_consensus.h"
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


  if (NULL == CU_add_test(pPopGraphSuite, "Test element - get edge", test_get_edge)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test element - get edge copy", test_get_edge_copy)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Loading two people, one fasta each, in same population, and print supernodes. \n\tEach person has 2 reads, and we have manually determined what the supernodes should be", test_load_two_people_in_same_populations_and_print_separately_their_supernodes)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Load four people in one population, each consists of one toy chromosome.\n\tThree differ at one locus, and one at another.\n\t Find their individual supernodes correctly", test_take_four_people_each_with_one_read_and_find_variants))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }
  if (NULL == CU_add_test(pPopGraphSuite, "Load two people sharing an Alu, differing by sequence before and after the Alu. Find supernodes.", test_take_two_people_sharing_an_alu_and_find_supernodes))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }
  if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that given a kmer, we can find the first node in the supernode in which that kmer lies.", test_find_first_node_in_supernode))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }
  if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that can pull out required subsection of a supernode.", test_correctly_find_subsection_of_supernode))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }



 

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}

