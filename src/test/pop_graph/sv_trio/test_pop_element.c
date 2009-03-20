#include <CUnit.h>
#include <Basic.h>
#include <element.h>
#include <hash_table.h>
#include "test_pop_element.h"
#include <stdlib.h>
#include <stdio.h>

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


void test_mark_chromosome_overlap()
{
  dBNode* e = new_element();

  //mark chromosome 1 in forward direction
  db_node_mark_chromosome_overlap(e,1,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==1 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 1 in reverse direction also
  db_node_mark_chromosome_overlap(e,1,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==3 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 2:
  e = new_element();

  //mark chromosome 2 in forward direction
  db_node_mark_chromosome_overlap(e,2,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==4 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 2 in reverse direction also
  db_node_mark_chromosome_overlap(e,2,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==12 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);



  //now chromosome 3:
  e = new_element();

  //mark chromosome 3 in forward direction
  db_node_mark_chromosome_overlap(e,3,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==16 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 3 in reverse direction also
  db_node_mark_chromosome_overlap(e,3,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==48 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 4:
  e = new_element();

  //mark chromosome 4 in forward direction
  db_node_mark_chromosome_overlap(e,4,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==64 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 4 in reverse direction also
  db_node_mark_chromosome_overlap(e,4,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==192 ); // = 64 + 128
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 5:
  e = new_element();

  //mark chromosome 5 in forward direction
  db_node_mark_chromosome_overlap(e,5,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==1 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 5 in reverse direction also
  db_node_mark_chromosome_overlap(e,5,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==3 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 6:
  e = new_element();

  //mark chromosome 6 in forward direction
  db_node_mark_chromosome_overlap(e,6,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==4 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 6 in reverse direction also
  db_node_mark_chromosome_overlap(e,6,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==12 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 7:
  e = new_element();

  //mark chromosome 7 in forward direction
  db_node_mark_chromosome_overlap(e,7,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==16 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 7 in reverse direction also
  db_node_mark_chromosome_overlap(e,7,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==48 );//16+32
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 8:
  e = new_element();

  //mark chromosome 8 in forward direction
  db_node_mark_chromosome_overlap(e,8,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==64 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 8 in reverse direction also
  db_node_mark_chromosome_overlap(e,8,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==192 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);



  //now chromosome 9:
  e = new_element();

  //mark chromosome 9 in forward direction
  db_node_mark_chromosome_overlap(e,9,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==1 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 9 in reverse direction also
  db_node_mark_chromosome_overlap(e,9,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==3 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 10:
  e = new_element();

  //mark chromosome 10 in forward direction
  db_node_mark_chromosome_overlap(e,10,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==4 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 10 in reverse direction also
  db_node_mark_chromosome_overlap(e,10,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==12 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 11:
  e = new_element();

  //mark chromosome 11 in forward direction
  db_node_mark_chromosome_overlap(e,11,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==16 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 11 in reverse direction also
  db_node_mark_chromosome_overlap(e,11,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==48 );//16+32
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 12:
  e = new_element();

  //mark chromosome 12 in forward direction
  db_node_mark_chromosome_overlap(e,12,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==64 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 12 in reverse direction also
  db_node_mark_chromosome_overlap(e,12,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==192 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 13:
  e = new_element();

  //mark chromosome 13 in forward direction
  db_node_mark_chromosome_overlap(e,13,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==1 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 13 in reverse direction also
  db_node_mark_chromosome_overlap(e,13,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==3 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 14:
  e = new_element();

  //mark chromosome 13 in forward direction
  db_node_mark_chromosome_overlap(e,14,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==4 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 14 in reverse direction also
  db_node_mark_chromosome_overlap(e,14,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==12 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 15:
  e = new_element();

  //mark chromosome 15 in forward direction
  db_node_mark_chromosome_overlap(e,15,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==16 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 15 in reverse direction also
  db_node_mark_chromosome_overlap(e,15,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==48 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 16:
  e = new_element();

  //mark chromosome 16 in forward direction
  db_node_mark_chromosome_overlap(e,16,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==64 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 16 in reverse direction also
  db_node_mark_chromosome_overlap(e,16,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==192 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 17:
  e = new_element();

  //mark chromosome 17 in forward direction
  db_node_mark_chromosome_overlap(e,17,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==1 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 17 in reverse direction also
  db_node_mark_chromosome_overlap(e,17,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==3 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 18:
  e = new_element();

  //mark chromosome 18 in forward direction
  db_node_mark_chromosome_overlap(e,18,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==4 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 18 in reverse direction also
  db_node_mark_chromosome_overlap(e,18,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==12 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);

  //now chromosome 19:
  e = new_element();

  //mark chromosome 19 in forward direction
  db_node_mark_chromosome_overlap(e,19,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==16 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 19 in reverse direction also
  db_node_mark_chromosome_overlap(e,19,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==48 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 20:
  e = new_element();

  //mark chromosome 20 in forward direction
  db_node_mark_chromosome_overlap(e,20,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==64 );
  CU_ASSERT(  e->chrom_xs[5]==0 );
  
  //makr chromosome 20 in reverse direction also
  db_node_mark_chromosome_overlap(e,20,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==192 );
  CU_ASSERT(  e->chrom_xs[5]==0 );

  free_element(&e);


  //now chromosome 21
  e = new_element();

  //mark chromosome 21 in forward direction
  db_node_mark_chromosome_overlap(e,21,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==1 );
  
  //makr chromosome 21 in reverse direction also
  db_node_mark_chromosome_overlap(e,21,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==3 );

  free_element(&e);

  //now chromosome 22:
  e = new_element();

  //mark chromosome 22 in forward direction
  db_node_mark_chromosome_overlap(e,22,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==4 );
  
  //makr chromosome 22  in reverse direction also
  db_node_mark_chromosome_overlap(e,22,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==12 );

  free_element(&e);

  //now chromosome 23 - ie X:
  e = new_element();

  //mark chromosome 23 in forward direction
  db_node_mark_chromosome_overlap(e,23,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==16 );
  
  //makr chromosome 23 in reverse direction also
  db_node_mark_chromosome_overlap(e,23,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==48 );

  free_element(&e);


  //now chromosome 24 - ie Y:
  e = new_element();

  //mark chromosome 24 in forward direction
  db_node_mark_chromosome_overlap(e,24,forward);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==64 );
  
  //makr chromosome 24 in reverse direction also
  db_node_mark_chromosome_overlap(e,24,reverse);
  //now check it
  CU_ASSERT(  e->chrom_xs[0]==0 );
  CU_ASSERT(  e->chrom_xs[1]==0 );
  CU_ASSERT(  e->chrom_xs[2]==0 );
  CU_ASSERT(  e->chrom_xs[3]==0 );
  CU_ASSERT(  e->chrom_xs[4]==0 );
  CU_ASSERT(  e->chrom_xs[5]==192 );

  free_element(&e);

}



void test_has_at_most_one_chromosome_intersection()
{
  int answer=-99;
  dBNode* e = new_element();


  //first check that if there really is only 1 intersection, it finds it:
  
  e->chrom_xs[0]=1;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==1);
  answer=-99;
  e->chrom_xs[0]=2;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==1);
  answer=-99;
  e->chrom_xs[0]=3;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==1);
  answer=-99;

  e->chrom_xs[0]=4;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==2);
  answer=-99;
  e->chrom_xs[0]=8;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==2);
  answer=-99;
  e->chrom_xs[0]=12;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==2);
  answer=-99;


  e->chrom_xs[0]=16;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==3);
  answer=-99;
  e->chrom_xs[0]=32;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==3);
  answer=-99;
  e->chrom_xs[0]=48;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==3);
  answer=-99;


  e->chrom_xs[0]=64;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==4);
  answer=-99;
  e->chrom_xs[0]=128;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==4);
  answer=-99;
  e->chrom_xs[0]=192;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==4);
  answer=-99;

  free_element(&e);
  e=new_element();
  
  e->chrom_xs[1]=1;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==5);
  answer=-99;
  e->chrom_xs[1]=2;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==5);
  answer=-99;
  e->chrom_xs[1]=3;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==5);
  answer=-99;

  e->chrom_xs[1]=4;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==6);
  answer=-99;
  e->chrom_xs[1]=8;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==6);
  answer=-99;
  e->chrom_xs[1]=12;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==6);
  answer=-99;


  e->chrom_xs[1]=16;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==7);
  answer=-99;
  e->chrom_xs[1]=32;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==7);
  answer=-99;
  e->chrom_xs[1]=48;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==7);
  answer=-99;


  e->chrom_xs[1]=64;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==8);
  answer=-99;
  e->chrom_xs[1]=128;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==8);
  answer=-99;
  e->chrom_xs[1]=192;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==8);
  answer=-99;


  free_element(&e);
  e=new_element();


  e->chrom_xs[2]=1;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==9);
  answer=-99;
  e->chrom_xs[2]=2;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==9);
  answer=-99;
  e->chrom_xs[2]=3;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==9);
  answer=-99;

  e->chrom_xs[2]=4;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==10);
  answer=-99;
  e->chrom_xs[2]=8;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==10);
  answer=-99;
  e->chrom_xs[2]=12;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==10);
  answer=-99;


  e->chrom_xs[2]=16;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==11);
  answer=-99;
  e->chrom_xs[2]=32;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==11);
  answer=-99;
  e->chrom_xs[2]=48;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==11);
  answer=-99;


  e->chrom_xs[2]=64;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==12);
  answer=-99;
  e->chrom_xs[2]=128;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==12);
  answer=-99;
  e->chrom_xs[2]=192;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==12);
  answer=-99;

  free_element(&e);
  e=new_element();


  e->chrom_xs[3]=1;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==13);
  answer=-99;
  e->chrom_xs[3]=2;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==13);
  answer=-99;
  e->chrom_xs[3]=3;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==13);
  answer=-99;

  e->chrom_xs[3]=4;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==14);
  answer=-99;
  e->chrom_xs[3]=8;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==14);
  answer=-99;
  e->chrom_xs[3]=12;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==14);
  answer=-99;


  e->chrom_xs[3]=16;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==15);
  answer=-99;
  e->chrom_xs[3]=32;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==15);
  answer=-99;
  e->chrom_xs[3]=48;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==15);
  answer=-99;


  e->chrom_xs[3]=64;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==16);
  answer=-99;
  e->chrom_xs[3]=128;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==16);
  answer=-99;
  e->chrom_xs[3]=192;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==16);
  answer=-99;

  free_element(&e);
  e=new_element();


  e->chrom_xs[4]=1;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==17);
  answer=-99;
  e->chrom_xs[4]=2;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==17);
  answer=-99;
  e->chrom_xs[4]=3;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==17);
  answer=-99;

  e->chrom_xs[4]=4;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==18);
  answer=-99;
  e->chrom_xs[4]=8;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==18);
  answer=-99;
  e->chrom_xs[4]=12;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==18);
  answer=-99;


  e->chrom_xs[4]=16;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==19);
  answer=-99;
  e->chrom_xs[4]=32;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==19);
  answer=-99;
  e->chrom_xs[4]=48;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==19);
  answer=-99;


  e->chrom_xs[4]=64;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==20);
  answer=-99;
  e->chrom_xs[4]=128;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==20);
  answer=-99;
  e->chrom_xs[4]=192;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==20);
  answer=-99;

  free_element(&e);


  //zam

  e=new_element();


  e->chrom_xs[5]=1;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==21);
  answer=-99;
  e->chrom_xs[5]=2;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==21);
  answer=-99;
  e->chrom_xs[5]=3;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==21);
  answer=-99;

  e->chrom_xs[5]=4;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==22);
  answer=-99;
  e->chrom_xs[5]=8;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==22);
  answer=-99;
  e->chrom_xs[5]=12;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==22);
  answer=-99;


  e->chrom_xs[5]=16;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==23);
  answer=-99;
  e->chrom_xs[5]=32;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==23);
  answer=-99;
  e->chrom_xs[5]=48;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==23);
  answer=-99;


  e->chrom_xs[5]=64;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==24);
  answer=-99;
  e->chrom_xs[5]=128;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==24);
  answer=-99;
  e->chrom_xs[5]=192;
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==24);
  answer=-99;

  free_element(&e);


  //now check that when you have NONE, you get the right answer

  e=new_element();
  CU_ASSERT(db_node_has_at_most_one_intersecting_chromosome(e, &answer));
  CU_ASSERT(answer==0);
  answer=-99;
  free_element(&e);

  //now try a few examples where there are more than one intersection

  e=new_element();
  e->chrom_xs[0]=64;
  e->chrom_xs[1]=1;
  CU_ASSERT(! (db_node_has_at_most_one_intersecting_chromosome(e, &answer)));
  CU_ASSERT(answer==-1);
  answer=-99;
  free_element(&e);

  e=new_element();
  e->chrom_xs[1]=64;
  e->chrom_xs[3]=64;
  CU_ASSERT(! (db_node_has_at_most_one_intersecting_chromosome(e, &answer)));
  CU_ASSERT(answer==-1);
  answer=-99;
  free_element(&e);

  e=new_element();
  e->chrom_xs[4]=123;
  e->chrom_xs[5]=5;
  CU_ASSERT(! (db_node_has_at_most_one_intersecting_chromosome(e, &answer)));
  CU_ASSERT(answer==-1);
  answer=-99;
  free_element(&e);

}


void test_get_chromosome_overlap_direction()
{
  dBNode* e = new_element();

  e->chrom_xs[0]=1;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,9)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,10)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,11)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,12)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,13)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,14)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,15)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,16)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,17)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,18)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,19)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,20)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,21)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,22)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,23)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,24)==does_not_overlap     );

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==overlaps_forwards_only     );
  e->chrom_xs[0]=2;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==overlaps_reverse_only     );
  e->chrom_xs[0]=3;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==overlaps_both_directions     );
  e->chrom_xs[0]=0;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==does_not_overlap     );



  e->chrom_xs[0]=4;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,9)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,10)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,11)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,12)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,13)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,14)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,15)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,16)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,17)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,18)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,19)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,20)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,21)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,22)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,23)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,24)==does_not_overlap     );

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==overlaps_forwards_only     );
  e->chrom_xs[0]=8;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==overlaps_reverse_only     );
  e->chrom_xs[0]=12;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==overlaps_both_directions     );
  e->chrom_xs[0]=0;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==does_not_overlap     );


  e->chrom_xs[0]=16;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,9)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,10)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,11)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,12)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,13)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,14)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,15)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,16)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,17)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,18)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,19)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,20)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,21)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,22)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,23)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,24)==does_not_overlap     );

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==overlaps_forwards_only     );
  e->chrom_xs[0]=32;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==overlaps_reverse_only     );
  e->chrom_xs[0]=48;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==overlaps_both_directions     );
  e->chrom_xs[0]=0;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==does_not_overlap     );


  e->chrom_xs[0]=64;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,9)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,10)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,11)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,12)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,13)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,14)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,15)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,16)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,17)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,18)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,19)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,20)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,21)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,22)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,23)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,24)==does_not_overlap     );

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==overlaps_forwards_only     );
  e->chrom_xs[0]=128;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==overlaps_reverse_only     );
  e->chrom_xs[0]=192;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==overlaps_both_directions     );
  e->chrom_xs[0]=0;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==does_not_overlap     );


  //zamzam

  e->chrom_xs[1]=1;

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,9)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,10)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,11)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,12)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,13)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,14)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,15)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,16)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,17)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,18)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,19)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,20)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,21)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,22)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,23)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,24)==does_not_overlap     );

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==overlaps_forwards_only     );
  e->chrom_xs[1]=2;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==overlaps_reverse_only     );
  e->chrom_xs[1]=3;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==overlaps_both_directions     );
  e->chrom_xs[1]=0;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==does_not_overlap     );



  e->chrom_xs[1]=4;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,9)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,10)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,11)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,12)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,13)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,14)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,15)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,16)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,17)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,18)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,19)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,20)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,21)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,22)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,23)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,24)==does_not_overlap     );

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==overlaps_forwards_only     );
  e->chrom_xs[1]=8;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==overlaps_reverse_only     );
  e->chrom_xs[1]=12;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==overlaps_both_directions     );
  e->chrom_xs[1]=0;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==does_not_overlap     );


  e->chrom_xs[1]=16;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,9)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,10)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,11)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,12)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,13)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,14)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,15)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,16)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,17)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,18)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,19)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,20)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,21)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,22)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,23)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,24)==does_not_overlap     );

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==overlaps_forwards_only     );
  e->chrom_xs[1]=32;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==overlaps_reverse_only     );
  e->chrom_xs[1]=48;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==overlaps_both_directions     );
  e->chrom_xs[1]=0;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==does_not_overlap     );


  e->chrom_xs[1]=64;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,1)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,2)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,3)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,4)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,5)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,6)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,7)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,9)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,10)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,11)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,12)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,13)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,14)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,15)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,16)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,17)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,18)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,19)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,20)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,21)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,22)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,23)==does_not_overlap     );
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,24)==does_not_overlap     );

  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==overlaps_forwards_only     );
  e->chrom_xs[1]=128;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==overlaps_reverse_only     );
  e->chrom_xs[1]=192;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==overlaps_both_directions     );
  e->chrom_xs[1]=0;
  CU_ASSERT( db_node_get_direction_through_node_in_which_chromosome_passes(e,8)==does_not_overlap     );




  free_element(&e);
}
