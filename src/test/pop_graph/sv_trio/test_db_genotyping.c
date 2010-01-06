#include <CUnit.h>
#include <Basic.h>
#include <test_db_genotyping.h>
#include <db_genotyping.h>

void test_get_node_multiplicities()
{
  //create a bunch of nodes, all o which have non-zero edges in colour 0.
  //don't make any effort to create them properly - this function we are testing
  // only looks at the memory address of the node, and whether it has an edge in the right colour

  dBNode* e1 = new_element();
  e1->individual_edges[0]=1;
  dBNode* e2 = new_element();
  e2->individual_edges[0]=1;
  dBNode* e3 = new_element();
  e3->individual_edges[0]=1;
  dBNode* e4 = new_element();
  e4->individual_edges[0]=1;
  dBNode* e5 = new_element();
  e5->individual_edges[0]=1;
  dBNode* e6 = new_element();
  e6->individual_edges[0]=1;
  dBNode* e7 = new_element();
  e7->individual_edges[0]=1;
  dBNode* e8 = new_element();
  e8->individual_edges[0]=1;
  dBNode* e9 = new_element();
  e9->individual_edges[0]=1;
  dBNode* e10 = new_element();
  e10->individual_edges[0]=1;


  //first test - two branches of equl length, each with globally unique nodes
  dBNode* array1[4]={e1,e2,e3,e4};
  dBNode* array2[4]={e5,e6,e7,e8};

  int mult1[4];
  int mult2[4];
  int mult1_in_2[4];
  int mult2_in_1[4];
  int* mult_ptr1[4]={&mult1[0], &mult1[1], &mult1[2], &mult1[3]};
  int* mult_ptr2[4]={&mult2[0], &mult2[1], &mult2[2], &mult2[3]};
  int* mult1_in_2_ptr[4]={&mult1_in_2[0], &mult1_in_2[1], &mult1_in_2[2], &mult1_in_2[3]};
  int* mult2_in_1_ptr[4]={&mult2_in_1[0], &mult2_in_1[1], &mult2_in_1[2], &mult2_in_1[3]};

  get_node_multiplicities(array1, 4, array2, 4, mult_ptr1, mult_ptr2, mult1_in_2_ptr, mult2_in_1_ptr, 
			  &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);

  int p;
  for (p=0; p<3; p++)
    {
      CU_ASSERT(mult1[p]==1);
      if (mult1[p] != 1)
	{
	  printf("Failed for p=%d. mult1[p] =%d\n", p, mult1[p]);
	}
      CU_ASSERT(mult2[p]==1);
      CU_ASSERT(mult1_in_2[p]==0);
      CU_ASSERT(mult2_in_1[p]==0);
    }


  //second test - two branches of equal length, one branch containing a repeated node
  array1[1]=e1;//so array1[0]==array1[1]
  get_node_multiplicities(array1, 4, array2, 4, mult_ptr1, mult_ptr2, mult1_in_2_ptr, mult2_in_1_ptr,
                          &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);

  CU_ASSERT(mult1[0]==2);
  CU_ASSERT(mult1[1]==2);
  CU_ASSERT(mult2[0]==1);
  CU_ASSERT(mult2[1]==1);
  CU_ASSERT(mult1_in_2[0]==0);
  CU_ASSERT(mult1_in_2[1]==0);
  CU_ASSERT(mult2_in_1[0]==0);
  CU_ASSERT(mult2_in_1[1]==0);

  for (p=2; p<4; p++)
    {
      CU_ASSERT(mult1[p]==1);
      if (mult1[p] != 1)
	{
	  printf("Failed for p=%d. mult1[p] =%d\n", p, mult1[p]);
	}
      CU_ASSERT(mult2[p]==1);
      CU_ASSERT(mult1_in_2[p]==0);
      CU_ASSERT(mult2_in_1[p]==0);
    }
  

  //third test - two branches, first shorter than second, with some nodes repeated within and between branches
  array1[1]=e2;//reset the second element in array1 to what it was before the last test
  //recall array1[4]={e1,e2,e3,e4}
  dBNode* array3[7]={e4,e1,e9,e10,e9,e9,e2};
  int mult3[7];
  int mult1_in_3[4];
  int mult3_in_1[7];
  int* mult_ptr3[7]={&mult3[0], &mult3[1], &mult3[2], &mult3[3], &mult3[4], &mult3[5], &mult3[6]};
  int* mult1_in_3_ptr[4]={&mult1_in_3[0], &mult1_in_3[1], &mult1_in_3[2], &mult1_in_3[3]};
  int* mult3_in_1_ptr[7]={&mult3_in_1[0], &mult3_in_1[1], &mult3_in_1[2], &mult3_in_1[3], &mult3_in_1[4], &mult3_in_1[5], &mult3_in_1[6]};

  get_node_multiplicities(array1, 4, array3, 7, mult_ptr1, mult_ptr3, mult1_in_3_ptr, mult3_in_1_ptr,
                          &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);


  CU_ASSERT(mult1[0]==1);
 
  CU_ASSERT(mult1[1]==1);
  CU_ASSERT(mult1[2]==1);
  CU_ASSERT(mult1[3]==1);

  CU_ASSERT(mult3[0]==1);
  CU_ASSERT(mult3[1]==1);
  CU_ASSERT(mult3[2]==3);
  CU_ASSERT(mult3[3]==1);
  CU_ASSERT(mult3[4]==3);
  CU_ASSERT(mult3[5]==3);
  CU_ASSERT(mult3[6]==1);

  CU_ASSERT(mult1_in_3[0]==1);
  CU_ASSERT(mult1_in_3[1]==1);
  CU_ASSERT(mult1_in_3[2]==0);
  CU_ASSERT(mult1_in_3[3]==1);

  CU_ASSERT(mult3_in_1[0]==1);
  CU_ASSERT(mult3_in_1[1]==1);
  CU_ASSERT(mult3_in_1[2]==0);
  CU_ASSERT(mult3_in_1[3]==0);
  CU_ASSERT(mult3_in_1[4]==0);
  CU_ASSERT(mult3_in_1[5]==0);
  CU_ASSERT(mult3_in_1[6]==1);



}
