//#include "test_dB_graph.h"
#include "test_dB_graph_node.h"
#include "test_dB_graph_population.h"
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
  pPopGraphSuite = CU_add_suite("Test SV trio -specific aspects of the population/people aware dB graph", NULL, NULL);
  if (NULL == pPopGraphSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suites */

  if (NULL == CU_add_test(pPopGraphSuite, "Test element - get edge", test_get_edge)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test element - get edge copy", test_get_edge_copy)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test element - mark chromosome overlap", test_mark_chromosome_overlap)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test element -check if element overlaps <=1 chromosome", test_has_at_most_one_chromosome_intersection)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Test that can identify supernode ends",   test_is_supernode_end)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Loading two people, one fasta each, in same population, and print supernodes. \n\tEach person has 2 reads, and we have manually determined what the supernodes should be", test_load_two_people_in_same_populations_and_print_separately_their_supernodes)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, "Load three people in one population, each consists of one toy chromosome.\n\t Find their individual supernodes correctly", test_take_three_people_each_with_one_read_and_find_variants))
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

  if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that given a node in a person's graoh, you can find the next node in the supernode.", test_find_next_node_in_supernode))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }

  if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that can pull out required subsection of a supernode.", test_correctly_find_subsection_of_supernode))
    {
      CU_cleanup_registry();
      return CU_get_error();
    }
    if (NULL == CU_add_test(pPopGraphSuite, "Load two people in one population, and test that for each separately, can find the best sub_supernode, requiring certain people_coverage", test_find_best_subsection_of_supernode_with_just_two_people))
    {
     CU_cleanup_registry();
     return CU_get_error();
    }
    //   if (NULL == CU_add_test(pPopGraphSuite, "Load five people in one population, sharing different amounts of sequence. Test finding population consensus supernode.", test_get_population_consensus_supernode))
    //{
    //CU_cleanup_registry();
    //return CU_get_error();
    //}
    if (NULL == CU_add_test(pPopGraphSuite, "Check that correctly get stats on how many kmers are shared by 1,2,3,... people in a population", test_getting_stats_of_how_many_indivduals_share_a_node))
    {
     CU_cleanup_registry();
     return CU_get_error();
    }
    if (NULL == CU_add_test(pPopGraphSuite, "Load 24 fake (tiny) chromosomes and count intersections with a single person", test_loading_simple_fasta_and_getting_chromosome_intersections))
    {
     CU_cleanup_registry();
     return CU_get_error();
    }
    if (NULL == CU_add_test(pPopGraphSuite, "Check with simple examples that can correctly check supernode for whether each node intersects at most one chromosome, and if so, how many chromosomes altogether are overlapped", test_db_graph_do_all_nodes_in_supernode_intersect_at_most_one_chromosome))
    {
     CU_cleanup_registry();
     return CU_get_error();
    }
    //    if (NULL == CU_add_test(pPopGraphSuite, "Simple test that you can correctly print supernodes with chromosome intersections", test_printing_supernode_with_chromosome_intersections_simple))
    // {
    //CU_cleanup_registry();
    // return CU_get_error();
    // }
    if (NULL == CU_add_test(pPopGraphSuite, "Test that you can correctly print supernodes with chromosome intersections using simple alu example", test_printing_supernode_with_chromosome_intersections_alus))
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

