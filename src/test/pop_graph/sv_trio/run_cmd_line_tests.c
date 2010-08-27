//#include "test_dB_graph.h"
#include <cmd_line.h>
#include <test_cmd_line.h>

#include <CUnit.h>
#include <Basic.h>

int  main()
{

  CU_pSuite pPopGraphSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();

  /* add a suite to the registry */
  pPopGraphSuite = CU_add_suite("Test Command-Line of cortex_var", NULL, NULL);
  if (NULL == pPopGraphSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suites */

  if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for getting colours from input format", test_get_numbers_from_comma_sep_list)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  
  if (NULL == CU_add_test(pPopGraphSuite, "Test utility function for parsing args (lists of colours separated by /) for detect_bubbles", test_parse_colourinfo_argument)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Test that when given good input, internal variables are correctly set", test_parse_cmdline_inner_loop_are_basic_variables_correctly_set )) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  



 

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}

