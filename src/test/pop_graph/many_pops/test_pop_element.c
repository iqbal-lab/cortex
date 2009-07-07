#include <CUnit.h>
#include <Basic.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include "test_pop_element.h"
#include <stdlib.h>


void test_get_edge()
{

  Element* my_element=new_element();

  my_element->individual_edges[0]=1;
  my_element->individual_edges[1]=3;

  Edges* edges=get_edge( *my_element, individual_edge_array, 0);
  CU_ASSERT(*edges==1);
  edges=get_edge( *my_element, individual_edge_array, 1);
  CU_ASSERT(*edges==3);

  free_element(&my_element);




  my_element=new_element();
 
  my_element->individual_edges[0]=2;
  my_element->individual_edges[1]=4;

  edges=get_edge( *my_element, individual_edge_array, 0);
  CU_ASSERT(*edges==2);
  edges=get_edge( *my_element, individual_edge_array, 1);
  CU_ASSERT(*edges==4);

  free_element(&my_element);


}



void test_get_edge_copy()
{

  Element* my_element=new_element();

  my_element->individual_edges[0]=1;
  my_element->individual_edges[1]=3;

  Edges edges=get_edge_copy( *my_element, individual_edge_array, 0);
  CU_ASSERT(edges==1);
  edges=get_edge_copy( *my_element, individual_edge_array, 1);
  CU_ASSERT(edges==3);

  free_element(&my_element);




  my_element=new_element();
 
  my_element->individual_edges[0]=2;
  my_element->individual_edges[1]=4;

  edges=get_edge_copy( *my_element, individual_edge_array, 0);
  CU_ASSERT(edges==2);
  edges=get_edge_copy( *my_element, individual_edge_array, 1);
  CU_ASSERT(edges==4);

  free_element(&my_element);


}


