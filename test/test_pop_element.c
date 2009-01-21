#include <CUnit.h>
#include <Basic.h>
#include <element.h>
#include <hash_table.h>
#include <stdlib.h>


void test_get_edge()
{

  Element* my_element=new_element();

  my_element->individual_edges[0]=1;
  my_element->individual_edges[1]=3;

  Edges* edges=get_edge( *my_element, individual_edge_array, 0);
  CU_ASSERT(*edges==1);
  Edges* edges=get_edge( *my_element, individual_edge_array, 1);
  CU_ASSERT(*edges==3);

}

