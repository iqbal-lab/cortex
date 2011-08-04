#include <CUnit.h>
#include <Basic.h>
#include <test_db_complex_genotyping.h>
#include <db_complex_genotyping.h>
#include <stdlib.h>
#include <element.h>

void test_initialise_multiplicities_of_allele_nodes_wrt_both_alleles()
{
  //create a bunch of nodes, all o which have non-zero edges in colour 0.
  //don't make any effort to create them properly - this function we are testing
  //only looks at the memory address of the node, and whether it has an edge in the right colour

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


  //first test - two branches of equal length, each with globally unique nodes
  dBNode* array1[4]={e1,e2,e3,e4};
  dBNode* array2[4]={e5,e6,e7,e8};
  VariantBranchesAndFlanks var;
  set_variant_branches_and_flanks(&var,
				  NULL, NULL, 0,
				  array1, NULL, 4,
				  array2, NULL, 4,
				  NULL, NULL, 0,
				  unknown);

  
  int mult11[4];
  int mult22[4];
  int mult12[4];
  int mult21[4];
  MultiplicitiesAndOverlapsOfBiallelicVariant mobv;
  mobv.mult11=mult11;
  mobv.mult22=mult22;
  mobv.mult12=mult12;
  mobv.mult21=mult21;
  initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);

  int p;
  for (p=0; p<3; p++)
    {
      CU_ASSERT(mobv.mult11[p]==1);
      if (mobv.mult11[p] != 1)
	{
	  printf("Failed for p=%d. mult11[p] =%d\n", p, mobv.mult11[p]);
	}
      CU_ASSERT(mobv.mult22[p]==1);
      CU_ASSERT(mobv.mult12[p]==0);
      CU_ASSERT(mobv.mult21[p]==0);
    }


  //second test - two branches of equal length, one branch containing a repeated node
  array1[1]=e1;//so array1[0]==array1[1]
  initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);

  CU_ASSERT(mobv.mult11[0]==2);
  CU_ASSERT(mobv.mult11[1]==2);
  CU_ASSERT(mobv.mult22[0]==1);
  CU_ASSERT(mobv.mult22[1]==1);
  CU_ASSERT(mobv.mult12[0]==0);
  CU_ASSERT(mobv.mult12[1]==0);
  CU_ASSERT(mobv.mult21[0]==0);
  CU_ASSERT(mobv.mult21[1]==0);

  for (p=2; p<4; p++)
    {
      CU_ASSERT(mobv.mult11[p]==1);
      if (mobv.mult11[p] != 1)
	{
	  printf("Failed for p=%d. mult11[p] =%d\n", p, mobv.mult11[p]);
	}
      CU_ASSERT(mobv.mult22[p]==1);
      CU_ASSERT(mobv.mult12[p]==0);
      CU_ASSERT(mobv.mult21[p]==0);
    }
  

  //third test - two branches, first shorter than second, with some nodes repeated within and between branches
  array1[1]=e2;//reset the second element in array1 to what it was before the last test
  //recall array1[4]={e1,e2,e3,e4}
  dBNode* array3[7]={e4,e1,e9,e10,e9,e9,e2};


  set_variant_branches_and_flanks(&var,
				  NULL, NULL, 0,
				  array1, NULL, 4,
				  array3, NULL, 7,
				  NULL, NULL, 0,
				  unknown);

  
  int mult3_11[4];
  int mult3_22[7];
  int mult3_12[4];
  int mult3_21[7];
  mobv.mult11=mult3_11;
  mobv.mult22=mult3_22;
  mobv.mult12=mult3_12;
  mobv.mult21=mult3_21;
  initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, &mobv, true, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);



  CU_ASSERT(mobv.mult11[0]==1);
  CU_ASSERT(mobv.mult11[1]==1);
  CU_ASSERT(mobv.mult11[2]==1);
  CU_ASSERT(mobv.mult11[3]==1);

  CU_ASSERT(mobv.mult22[0]==1);
  CU_ASSERT(mobv.mult22[1]==1);
  CU_ASSERT(mobv.mult22[2]==3);
  CU_ASSERT(mobv.mult22[3]==1);
  CU_ASSERT(mobv.mult22[4]==3);
  CU_ASSERT(mobv.mult22[5]==3);
  CU_ASSERT(mobv.mult22[6]==1);

  CU_ASSERT(mobv.mult12[0]==1);
  CU_ASSERT(mobv.mult12[1]==1);
  CU_ASSERT(mobv.mult12[2]==0);
  CU_ASSERT(mobv.mult12[3]==1);

  CU_ASSERT(mobv.mult21[0]==1);
  CU_ASSERT(mobv.mult21[1]==1);
  CU_ASSERT(mobv.mult21[2]==0);
  CU_ASSERT(mobv.mult21[3]==0);
  CU_ASSERT(mobv.mult21[4]==0);
  CU_ASSERT(mobv.mult21[5]==0);
  CU_ASSERT(mobv.mult21[6]==1);

  free(e1);
  free(e2);
  free(e3);
  free(e4);
  free(e5);
  free(e6);
  free(e7);
  free(e8);
  free(e9);
  free(e10);

}


void test_calc_log_likelihood_of_genotype_with_complex_alleles()
{

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
  
}
