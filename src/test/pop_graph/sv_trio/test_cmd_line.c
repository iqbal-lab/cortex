#include <CUnit.h>
#include <Basic.h>
#include <cmd_line.h>
#include <string.h>
#include <element.h>

void test_get_numbers_from_comma_sep_list()
{

  //simple list
  char list[]="1,2,3";
  int len = strlen(list);
  int return_list[10];
  int max=10;
  int ans =  get_numbers_from_comma_sep_list(list, len, return_list, max);

  CU_ASSERT(return_list[0]==1);
  CU_ASSERT(return_list[1]==2);
  CU_ASSERT(return_list[2]==3);
  CU_ASSERT(ans==3);

  //list of just one number
  char list2[]="100";
  int len2=1;
  int return_list2[10];
  int max2=10;
  int ans2=get_numbers_from_comma_sep_list(list2, len2, return_list2, max2);
  
  CU_ASSERT(ans2==1);
  printf("ans is %d \n", ans);
  CU_ASSERT(return_list2[0]==100);


  //list containing negative number
  char list3[]="-17";
  int len3=1;
  int return_list3[10];
  int max3=10;
  int ans3=get_numbers_from_comma_sep_list(list3, len3, return_list3, max3);
  
  CU_ASSERT(ans3==-1);


  //list containing negative number and other numbers
  char list4[]="10,1,2,3,4,-17";
  int len4=1;
  int return_list4[10];
  int max4=10;
  int ans4=get_numbers_from_comma_sep_list(list4, len4, return_list4, max4);
  
  CU_ASSERT(ans4==-1);

  
  //list containing too many numbers
  char list5[]="2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3";
  int len5=1;
  int return_list5[10];
  int max5=10;
  int ans5=get_numbers_from_comma_sep_list(list5, len5, return_list5, max5);
  
  CU_ASSERT(ans5==-1);



  //list containing a non-number
  char list6[]="10,1,2,3,6,a,2,132";
  int len6=1;
  int return_list6[10];
  int max6=10;
  int ans6=get_numbers_from_comma_sep_list(list6, len6, return_list6, max6);
  
  CU_ASSERT(ans6==-1);



}



void test_parse_colourinfo_argument()
{

  CmdLine cmd_line;

  char arg[]="1,2/3,4";
  int len_arg = strlen(arg);
  char text[]="dummy text";
  int which=1;

  parse_colourinfo_argument(&cmd_line, arg, len_arg, text, which);

  //printf("Num us  %d, ",  cmd_line.num_colours_in_detect_bubbles1_first_colour_list);
  CU_ASSERT(cmd_line.num_colours_in_detect_bubbles1_first_colour_list ==2);
  CU_ASSERT(cmd_line.detect_bubbles1_first_colour_list[0]==1);
  CU_ASSERT(cmd_line.detect_bubbles1_first_colour_list[1]==2);

  CU_ASSERT(cmd_line.num_colours_in_detect_bubbles1_second_colour_list==2);
  CU_ASSERT(cmd_line.detect_bubbles1_second_colour_list[0]==3);
  CU_ASSERT(cmd_line.detect_bubbles1_second_colour_list[1]==4);


  CmdLine cmd_line2;
  char arg2[]="1/234,17";
  int len_arg2=strlen(arg2);
  which=2;
  
  parse_colourinfo_argument(&cmd_line2, arg2, len_arg2, text, which);
  

  CU_ASSERT(cmd_line2.num_colours_in_detect_bubbles2_first_colour_list ==1);
  CU_ASSERT(cmd_line2.detect_bubbles2_first_colour_list[0]==1);

  CU_ASSERT(cmd_line2.num_colours_in_detect_bubbles2_second_colour_list==2);
  CU_ASSERT(cmd_line2.detect_bubbles2_second_colour_list[0]==234);
  CU_ASSERT(cmd_line2.detect_bubbles2_second_colour_list[1]==17);



  CmdLine cmd_line3;
  char arg3[]="1,28,19,278,8888/234,17,8";
  int len_arg3=strlen(arg3);
  which=2;
  parse_colourinfo_argument(&cmd_line3, arg3, len_arg3, text, which);


  CU_ASSERT(cmd_line3.num_colours_in_detect_bubbles2_first_colour_list ==5);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[0]==1);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[1]==28);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[2]==19);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[3]==278);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[4]==8888);

  CU_ASSERT(cmd_line3.num_colours_in_detect_bubbles2_second_colour_list==3);
  CU_ASSERT(cmd_line3.detect_bubbles2_second_colour_list[0]==234);
  CU_ASSERT(cmd_line3.detect_bubbles2_second_colour_list[1]==17);
  CU_ASSERT(cmd_line3.detect_bubbles2_second_colour_list[2]==8);



}


void test_parse_cmdline_inner_loop_are_basic_variables_correctly_set()
{
  

  // --COLOUR_LIST
  char* argv1[] = {"cortex", "--colour_list=../data/test/pop_graph/cmd_line/colour_list_1"};
  
  //set up a default cmdline
  CmdLine test1;
  default_opts(&test1);
  char error_string1[LEN_ERROR_STRING];
  int err = parse_cmdline_inner_loop(2, argv1, sizeof(Element), &test1, error_string1);
  CU_ASSERT_STRING_EQUAL(test1.colour_list, "../data/test/pop_graph/cmd_line/colour_list_1");
  

  //--MULTICOLOUR_BIN

  char* argv2[] = {"cortex", "--multicolour_bin=../data/test/pop_graph/cmd_line/multicolour_bin_1.ctx"};  
  //set up a default cmdline
  CmdLine test2;
  default_opts(&test2);
  char error_string2[LEN_ERROR_STRING];
  printf("start the parse we care about\n");
  int err2 = parse_cmdline_inner_loop(2, argv2, sizeof(Element), &test2, error_string2);
  printf("end the parse we care about\n");
  printf("We expect ../data/test/pop_graph/cmd_line/multicolour_bin_1.ctx, In fact it is %s\n", test2.multicolour_bin);
  CU_ASSERT_STRING_EQUAL(test2.multicolour_bin, "../data/test/pop_graph/cmd_line/multicolour_bin_1.ctx");
  


  // SE LIST

  char* argv3[] = {"cortex", "--multicolour_bin=../data/test/pop_graph/cmd_line/se_list1"};
  
  //set up a default cmdline
  CmdLine test3;
  default_opts(&test3);

  char error_string3[LEN_ERROR_STRING];
  int err3 = parse_cmdline_inner_loop(2, argv3, sizeof(Element), &test3, error_string3);

  CU_ASSERT_STRING_EQUAL(test3.multicolour_bin, "../data/test/pop_graph/cmd_line/se_list1");


  //PE_LIST

  char* argv4[] = {"cortex", "--pe_list=../data/test/pop_graph/cmd_line/pe_list1,../data/test/pop_graph/cmd_line/pe_list2"};
  
  //set up a default cmdline
  CmdLine test4;
  default_opts(&test4);

  char error_string4[LEN_ERROR_STRING];
  int err4 = parse_cmdline_inner_loop(2, argv4, sizeof(Element), &test4, error_string4);

  CU_ASSERT_STRING_EQUAL(test4.pe_list_lh_mates, "../data/test/pop_graph/cmd_line/pe_list1");
  CU_ASSERT_STRING_EQUAL(test4.pe_list_rh_mates, "../data/test/pop_graph/cmd_line/pe_list2");
  

  




}

