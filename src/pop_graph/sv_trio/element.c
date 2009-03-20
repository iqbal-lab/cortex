
/*
 Elemenx.c -- implements the nodes of the dBruijn graph
 */

#include <element.h>
#include <stdlib.h>
#include <global.h>
#include <stdio.h>
#include <string.h>


//const int NUMBER_OF_INDIVIDUALS_PER_POPULATION = 5;
//const int NUMBER_OF_POPULATIONS = 2;

const char mask1 = 255-3; //11111100
const char mask2 = 255-12;//11110011
const char mask3 = 255-48;//11001111
const char mask4 = 255-192;//00111111


//currently noone calls this in normal use
// In normal use, the priority queue allocates space to put the eloement directly within,
// and calls element_initialise
Element* new_element()
{


  Element* e = malloc(sizeof(Element));

  if (e==NULL)
    {
      printf("Unable to allocate a new element");
      exit(1);
    }
  
  e->kmer=0;

  int i;
  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
      e->coverage[i]=0;
    }

  //chromosome intersections struct has 6 chars.
  for (i=0; i<6; i++)
    {
      e->chrom_xs[i]=0;
    }

  e->status=none;
  e->kmer=0;

  return e;
}


void free_element(Element** element)
{
  free(*element);
  *element=NULL;
}

//gets you a pointer to the edge you are referring to
Edges* get_edge(Element e, EdgeArrayType type,int index)
{

  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}
      return &e.individual_edges[index];
    }
  else 
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array");
      exit(1);
    }

  exit(1);
}


//return a copy of the edge you are referring to
Edges get_edge_copy(const Element e, EdgeArrayType type,int index)
{

  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}
      return e.individual_edges[index];
    }
  else 
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array");
      exit(1);
    }

  exit(1);
}



Edges get_union_of_edges(Element e)
{

  int i;
  Edges edges=0;

  for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      edges |= e.individual_edges[i];
    }

  return edges;
}








//adds edges from edge_char to the appropriate person/population edgeset, without removing existing edges
void add_edges(Element* e, EdgeArrayType type, int index, Edges edge_char)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}
      e->individual_edges[index] |= edge_char;
    }

  else
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array");
      exit(1);
    }
  
}


void set_edges(Element* e, EdgeArrayType type, int index, Edges edge_char)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}
      e->individual_edges[index] = edge_char;
    }

  else
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array");
      exit(1);
    }
  
}


void db_node_reset_all_edges_for_all_people_and_pops_to_zero(Element* e)
{
  int i;

    for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
    }

}

void reset_one_edge(Element* e, Orientation orientation, Nucleotide nucleotide, EdgeArrayType type, int index)
{
  if (type == individual_edge_array)
    {
      if (index>=NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  exit(1);
	}

      char edge = 1 << nucleotide;      
      if (orientation == reverse){
	edge <<= 4;
      }
      //toggle 1->0 0->1
      edge ^= (unsigned char) 0xFF; //xor with all 1's, ie 00010000 -> 11101111
      
      e->individual_edges[index] &= edge; //reset one edge
      

    }
  else
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array");
      exit(1);
    }
  
}


int element_get_number_of_people_or_pops_containing_this_element(Element* e, EdgeArrayType type, int index)
{
  int i;
  int count=0;
  if (type == individual_edge_array)
    {
      for (i=0; i< NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
	{
	  if ( (e->individual_edges)[i] != 0)
	    {
	      count++;
	    }
	}
    }
  else
    {
      printf("Coding error. Only expecting enum of edge array types to contain one type - individual_edge_array");
      exit(1);
    }

  return count;
}

boolean element_smaller(Element  e1, Element e2){
 
  return get_union_of_edges(e1)  <  get_union_of_edges(e2);
  
}



BinaryKmer element_get_kmer(Element * e){
  return e->kmer;
}

boolean element_is_key(Key key, Element e, short kmer_size){
  return key == e.kmer;
}

Key element_get_key(BinaryKmer kmer, short kmer_size){
  
  BinaryKmer rev_kmer = binary_kmer_reverse_complement(kmer,kmer_size);
  
  if (rev_kmer < kmer){
    kmer = rev_kmer;
  }

  return kmer;

}


void element_initialise(Element * e, Key kmer, short kmer_size){

  e->kmer = element_get_key(kmer, kmer_size);

  int i;
  for (i=0; i<NUMBER_OF_INDIVIDUALS_PER_POPULATION; i++)
    {
      e->individual_edges[i]=0;
      e->coverage[i]=0;
    }

  for (i=0; i<6; i++)
    {
      e->chrom_xs[i]=0;
    }
  e->status = none;
}


void db_node_increment_coverage(dBNode* e, EdgeArrayType type, int index)
{
  //do nothing for now. Need to check about chars
}



Overlap db_node_get_direction_through_node_in_which_chromosome_passes(dBNode* node, int which_chromosome)
{

  int which_char;//there are 6 chars in which the chromosome xs are encoded. numbered 0 to 5. 0 carries chromosomes 1-4, etc.
  int which_probe;//two 1's in the two bits that encode this chromsomes's info within the char. So for example chromosome 2 is in the 3rd and 4th bots from the right in the 0th char. so probe is 12 (4+8).
  int which_pair_of_bits; //0 means the first two from the right, 1 means the next two, etc

  if (which_chromosome==0)
    {
      return does_not_overlap;
    }
  else if ( (which_chromosome==1)|| (which_chromosome==2) || (which_chromosome==3) || (which_chromosome==4) )
    {
      which_char=0;
    }
  else if  ( (which_chromosome==5)|| (which_chromosome==6) ||(which_chromosome==7) || (which_chromosome==8) )
    {
      which_char=1;
    }
  else if ( (which_chromosome==9)|| (which_chromosome==10) ||(which_chromosome==11) || (which_chromosome==12) )
    {
      which_char=2;
    }
  else if ( (which_chromosome==13)|| (which_chromosome==14) ||(which_chromosome==15) || (which_chromosome==16) )
    {
      which_char=3;
    }
  else if ( (which_chromosome==17)|| (which_chromosome==18) ||(which_chromosome==19) || (which_chromosome==20) )
    {
      which_char=4;
    }
  else if ( (which_chromosome==21)|| (which_chromosome==22) ||(which_chromosome==23) || (which_chromosome==24) )
    {
      which_char=5;
    }
  else
    {
      printf("Do not call db_node_get_chromosome_overlap_direction with second argument not between 1 and 24. You used %d", which_chromosome);
    }

  if (  (which_chromosome==1) || (which_chromosome==5) || (which_chromosome==9) || (which_chromosome==13) || (which_chromosome==17) || (which_chromosome==21) )
    {
      which_probe=3; //1+2
      which_pair_of_bits=0;
    }
  else  if (  (which_chromosome==2) || (which_chromosome==6) || (which_chromosome==10) || (which_chromosome==14) || (which_chromosome==18) || (which_chromosome==22) )
    {
      which_probe=12; //4+8
      which_pair_of_bits=1;
    }
  else  if (  (which_chromosome==3) || (which_chromosome==7) || (which_chromosome==11) || (which_chromosome==15) || (which_chromosome==19) || (which_chromosome==23) )
    {
      which_probe=48; //16+32
      which_pair_of_bits=2;
    }
  else if (  (which_chromosome==4) || (which_chromosome==8) || (which_chromosome==12) || (which_chromosome==16) || (which_chromosome==20) || (which_chromosome==24) )
    {
      which_probe=192; //64+128
      which_pair_of_bits=3;
    }
  else
    {
      printf("programming error. this should not be possible.");
      exit(1);
    }

  int fw=1<<(which_pair_of_bits*2);
  int rev=2<<(which_pair_of_bits*2);
  int both=3<<(which_pair_of_bits*2);
  int none=0;
  
  if ( ((node->chrom_xs[which_char]) &  which_probe)==fw)
    {
      return overlaps_forwards_only;
    }
  else if ( ((node->chrom_xs[which_char]) &  which_probe)==rev)
    {
      return overlaps_reverse_only;
    }
  else if ( ((node->chrom_xs[which_char]) &  which_probe)==both)
    {
      return overlaps_both_directions;
    }
  else if ( ((node->chrom_xs[which_char]) &  which_probe)==none)
    {
      return does_not_overlap;
    }
  else
    {
      printf("something wrong with finding direction of overlap");
      exit(1);
    }
  

	
}

char* overlap_to_char(Overlap ov, char* pre_alloced_string)
{
  pre_alloced_string[0]='\0';

  if (ov==overlaps_forwards_only)
    {
      strcat(pre_alloced_string, "F");
    }
  else if (ov==overlaps_reverse_only)
    {
      strcat(pre_alloced_string, "R");
    }
  else if (ov==overlaps_both_directions)
    {
      strcat(pre_alloced_string, "B");
    }
  else if (ov==does_not_overlap)
    {
      //want to return null string in this case
    }

  return pre_alloced_string;
}


char* compare_chrom_overlap_and_supernode_direction(Overlap ov, Orientation o, char* pre_alloced_string)
{
  pre_alloced_string[0]='\0';

  if (ov==overlaps_both_directions)
    {
      strcat(pre_alloced_string,"B");
    }
  else if (ov==does_not_overlap)
    {
      //do nothing - return null string
    }
  else if ( (ov==overlaps_forwards_only) && (o==forward) )
    {
      strcat(pre_alloced_string,"F");
    }
  else if ( (ov==overlaps_reverse_only) && (o==reverse) )
    {
      strcat(pre_alloced_string,"F");
    }
  else if ( (ov==overlaps_forwards_only) && (o==reverse) )
    {
      strcat(pre_alloced_string,"R");
    }
  else if ( (ov==overlaps_reverse_only) && (o==forward) )
    {
      strcat(pre_alloced_string,"R");
    }
  else
    {
      printf("problem with compairng overlap and orientation");
      exit(1);
    }
  return pre_alloced_string;

}



//of each 2 bit pair, least sig bit referes to forward and more sig to reverse. +1 means is present in chromosome, 0 means is not.
void db_node_mark_chromosome_overlap(dBNode* node, int which_chromosome, Orientation orientation)
{
  int which_char;

  if (which_chromosome==1)
    {
      which_char=0;

      if (orientation==forward)
	{
	  node->chrom_xs[which_char] |=  1;
	}
      else
	{
	  node->chrom_xs[which_char] |= 2;
	}
    }
  else if (which_chromosome==2)
    {
      which_char=0;
      if (orientation==forward)
	{
	  node->chrom_xs[which_char] |= 4;
	}
      else
	{
	  node->chrom_xs[which_char] |= 8;
	}
    }
  else if (which_chromosome==3)
    {
      which_char=0;
      if (orientation==forward)
	{
	  node->chrom_xs[which_char] |=  16;
	}
      else
	{
	  node->chrom_xs[which_char] |= 32;
	}
    }
  else if (which_chromosome==4)
    {
      which_char=0;
      if (orientation==forward)
	{
	  node->chrom_xs[which_char] |=  64;
	}
      else
	{
	  node->chrom_xs[which_char] |= 128;
	}
    }
  else if (which_chromosome==5)
    {
      which_char=1;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  1;
        }
      else
        {
          node->chrom_xs[which_char] |= 2;
        }

    }
  else if (which_chromosome==6)
    {
      which_char=1;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |= 4;
        }
      else
        {
          node->chrom_xs[which_char] |= 8;
        }

    }
  else if (which_chromosome==7)
    {
      which_char=1;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  16;
        }
      else
        {
          node->chrom_xs[which_char] |= 32;
        }

    }
  else if (which_chromosome==8)
    {
      which_char=1;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |= 64;
        }
      else
        {
          node->chrom_xs[which_char] |= 128;
        }

    }
  else if (which_chromosome==9)
    {
      which_char=2;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  1;
        }
      else
        {
          node->chrom_xs[which_char] |= 2;
        }

    }
  else if (which_chromosome==10)
    {
      which_char=2;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  4;
        }
      else
        {
          node->chrom_xs[which_char] |= 8;
        }

    }
  else if (which_chromosome==11)
    {
      which_char=2;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  16;
        }
      else
        {
          node->chrom_xs[which_char] |= 32;
        }

    }
  else if (which_chromosome==12)
    {
      which_char=2;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  64;
        }
      else
        {
          node->chrom_xs[which_char] |= 128;
        }

    }
  else if (which_chromosome==13)
    {
      which_char=3;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  1;
        }
      else
        {
          node->chrom_xs[which_char] |= 2;
        }

    }
  else if (which_chromosome==14)
    {
      which_char=3;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |= 4;
        }
      else
        {
          node->chrom_xs[which_char] |= 8;
        }

    }
  else if (which_chromosome==15)
    {
      which_char=3;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  16;
        }
      else
        {
          node->chrom_xs[which_char] |= 32;
        }

    }
  else if (which_chromosome==16)
    {
      which_char=3;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  64;
        }
      else
        {
          node->chrom_xs[which_char] |= 128;
        }

    }
  else if (which_chromosome==17)
    {
      which_char=4;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  1;
        }
      else
        {
          node->chrom_xs[which_char] |= 2;
        }

    }
  else if (which_chromosome==18)
    {
      which_char=4;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |= 4;
        }
      else
        {
          node->chrom_xs[which_char] |= 8;
        }

    }
  else if (which_chromosome==19)
    {
      which_char=4;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  16;
        }
      else
        {
          node->chrom_xs[which_char] |= 32;
        }

    }
  else if (which_chromosome==20)
    {
      which_char=4;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  64;
        }
      else
        {
          node->chrom_xs[which_char] |= 128;
        }

    }
  else if (which_chromosome==21)
    {
      which_char=5;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  1;
        }
      else
        {
          node->chrom_xs[which_char] |= 2;
        }

    }
  else if (which_chromosome==22)
    {
      which_char=5;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |= 4;
        }
      else
        {
          node->chrom_xs[which_char] |= 8;
        }

    }
  else if (which_chromosome==23) //X chromosome is number 23
    {
      which_char=5;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |=  16;
        }
      else
        {
          node->chrom_xs[which_char] |= 32;
        }

    }
  else if (which_chromosome==24)//Y chromosome is number 24
    {
      which_char=5;
      if (orientation==forward)
        {
          node->chrom_xs[which_char] |= 64;
        }
      else
        {
          node->chrom_xs[which_char] |= 128;
        }

    }

  
  


}


boolean db_node_has_at_most_one_intersecting_chromosome(dBNode* node, int* which_chromosome)
{

  //  if (  (node->chrom_xs[0]==0) && (node->chrom_xs[1]==0) && (node->chrom_xs[2]==0) && (node->chrom_xs[3]==0) && (node->chrom_xs[4]==0) && (node->chrom_xs[5]==0) )
  // {
  //   //intersects none at all
  //   *which_chromosome=0;
  //   return true;
  // } 
  //else if ((node->chrom_xs[0]==1) || (node->chrom_xs[0]==2) || (node->chrom_xs[0]==3)) //overlaps chrom 1, and none of 2,3,4
  // {
  //   if ( (node->chrom_xs[1]==0) && (node->chrom_xs[2]==0) && (node->chrom_xs[3]==0) && (node->chrom_xs[4]==0) && (node->chrom_xs[5]==0) ) //overlaps no other chromosome
  //	{
  //	  *which_chromosome=1;
  //	  return true;
  //	}
  //   else
  ////	{
  //	  *which_chromosome=-1;//overlaps more than 1
  //	  return false;
  //	}
  // }

  if (node==NULL)
    {
      printf("Dont call db_node_has_at_most_one_intersecting_chromosome wit null node");
      exit(1);
      return 0;
    }
  char chrom1_to_4 = node->chrom_xs[0];
  char chrom5_to_8 = node->chrom_xs[1];
  char chrom9_to_12 = node->chrom_xs[2];
  char chrom13_to_16 = node->chrom_xs[3];
  char chrom17_to_20 = node->chrom_xs[4];
  char chrom21_to_24 = node->chrom_xs[5];

  char all_except_1_to_4 = (chrom5_to_8 || chrom9_to_12 || chrom13_to_16 || chrom17_to_20 || chrom21_to_24);
  char all_except_5_to_8 = (chrom1_to_4 || chrom9_to_12 || chrom13_to_16 || chrom17_to_20 || chrom21_to_24);
  char all_except_9_to_12= (chrom1_to_4 || chrom5_to_8 || chrom13_to_16 || chrom17_to_20 || chrom21_to_24);
  char all_except_13_to_16=(chrom1_to_4 || chrom5_to_8 || chrom9_to_12 || chrom17_to_20 || chrom21_to_24);
  char all_except_17_to_20=(chrom1_to_4 || chrom5_to_8 || chrom9_to_12 || chrom13_to_16 || chrom21_to_24);
  char all_except_21_to_24=(chrom1_to_4 || chrom5_to_8 || chrom9_to_12 || chrom13_to_16 || chrom17_to_20);

  if ( (chrom1_to_4==0) && (chrom5_to_8==0) && (chrom9_to_12==0) && (chrom13_to_16==0) && (chrom17_to_20==0) && (chrom21_to_24==0) )
    {
      *which_chromosome=0;
      return true;
    }
  else if ( (all_except_1_to_4==0) && ( (chrom1_to_4 &  mask1) ==0) )
  {
    *which_chromosome=1;
    return true;
  }
  else if ( (all_except_1_to_4==0) && ( (chrom1_to_4 &  mask2) ==0) )
    {
      *which_chromosome=2;
      return true;
    }
  else if ( (all_except_1_to_4==0) && ( (chrom1_to_4 &  mask3) ==0) )
    {
      *which_chromosome=3;
      return true;
    }
  else if ( (all_except_1_to_4==0) && ((chrom1_to_4 &  mask4) ==0) )
    {
      *which_chromosome=4;
      return true;
    }
    else if ( (all_except_5_to_8==0) && ( (chrom5_to_8 &  mask1) ==0) )
    {
      *which_chromosome=5;
      return true;
    }
    else if ( (all_except_5_to_8==0) && ( (chrom5_to_8 &  mask2) ==0) )
    {
      *which_chromosome=6;
      return true;
    }
    else if ( (all_except_5_to_8==0) && ( (chrom5_to_8 &  mask3) ==0) )
    {
      *which_chromosome=7;
      return true;
    }
    else if ( (all_except_5_to_8==0) && ( (chrom5_to_8 &  mask4) ==0) )
    {
      *which_chromosome=8;
      return true;
    }
    else if ( (all_except_9_to_12==0) && ( (chrom9_to_12 &  mask1) ==0) )
    {
      *which_chromosome=9;
      return true;
    }
    else if ( (all_except_9_to_12==0) && ( (chrom9_to_12 &  mask2) ==0) )
    {
      *which_chromosome=10;
      return true;
    }
    else if ( (all_except_9_to_12==0) && ( (chrom9_to_12 &  mask3) ==0) )
    {
      *which_chromosome=11;
      return true;
    }
    else if ( (all_except_9_to_12==0) && ( (chrom9_to_12 &  mask4) ==0) )
    {
      *which_chromosome=12;
      return true;
    }
    else if ( (all_except_13_to_16==0) && ( (chrom13_to_16 &  mask1) ==0) )
    {
      *which_chromosome=13;
      return true;
    }
    else if ( (all_except_13_to_16==0) && ( (chrom13_to_16 &  mask2) ==0) )
    {
      *which_chromosome=14;
      return true;
    }
    else if ( (all_except_13_to_16==0) && ( (chrom13_to_16 &  mask3) ==0) )
    {
      *which_chromosome=15;
      return true;
    }
    else if ( (all_except_13_to_16==0) && ( (chrom13_to_16 &  mask4) ==0) )
    {
      *which_chromosome=16;
      return true;
    }
    else if ( (all_except_17_to_20==0) && ( (chrom17_to_20 &  mask1) ==0) )
    {
      *which_chromosome=17;
      return true;
    }
    else if ( (all_except_17_to_20==0) && ( (chrom17_to_20 &  mask2) ==0) )
    {
      *which_chromosome=18;
      return true;
    }
    else if ( (all_except_17_to_20==0) && ( (chrom17_to_20 &  mask3) ==0) )
    {
      *which_chromosome=19;
      return true;
    }
    else if ( (all_except_17_to_20==0) && ( (chrom17_to_20 &  mask4) ==0) )
    {
      *which_chromosome=20;
      return true;
    }
    else if ( (all_except_21_to_24==0) && ( (chrom21_to_24 &  mask1) ==0) )
    {
      *which_chromosome=21;
      return true;
    }
    else if ( (all_except_21_to_24==0) && ( (chrom21_to_24 &  mask2) ==0) )
    {
      *which_chromosome=22;
      return true;
    }
    else if ( (all_except_21_to_24==0) && ( (chrom21_to_24 &  mask3) ==0) )
    {
      *which_chromosome=23;
      return true;
    }
    else if ( (all_except_21_to_24==0) && ( (chrom21_to_24 &  mask4) ==0) )
    {
      *which_chromosome=24;
      return true;
    }

  //must be overlapping more than 1, so
  *which_chromosome=-1;
  return false;



}


Orientation opposite_orientation(Orientation o){
  return o ^ 1;
  
}

Orientation db_node_get_orientation(BinaryKmer k, dBNode * e, short kmer_size){

  if (e->kmer == k){
    return forward;
  }
  
  if (e->kmer == binary_kmer_reverse_complement(k,kmer_size)){
    return reverse;
  }

  exit(1);
  
}





//After specifying which individual or population you are talking about, this function
//adds one edge ("arrow") to the appropriate edge in the appropriate array in the element -- basically sets a bit in the correct edges char
void db_node_add_labeled_edge(dBNode * e, Orientation o, Nucleotide base, EdgeArrayType edge_type, int edge_index){

  //set edge 
  char edge = 1 << base; // A (0) -> 0001, C (1) -> 0010, G (2) -> 0100, T (3) -> 1000
  
  if (o == reverse){
    edge <<= 4; //move to next nibble 
  }

  //update node
  add_edges(e, edge_type, edge_index, edge);
  
}


//adding an edge between two nodes implies adding two labeled edges (one in each direction)
//be aware that in the case of self-loops in palindromes the two labeled edges collapse in one

boolean db_node_add_edge(dBNode * src_e, dBNode * tgt_e, Orientation src_o, Orientation tgt_o, short kmer_size, EdgeArrayType edge_type, int edge_index){

  BinaryKmer src_k, tgt_k; 

  src_k = src_e->kmer;
  tgt_k = tgt_e->kmer;
  char tmp_seq[kmer_size];
 
  if (src_o == reverse){
    src_k = binary_kmer_reverse_complement(src_k,kmer_size);
  }
    
  if (tgt_o == reverse){
    tgt_k = binary_kmer_reverse_complement(tgt_k,kmer_size);
  }
    
  
  if (DEBUG){
    printf("add edge %s -%c-> %s to edge type %d, and edge index %d\n",binary_kmer_to_seq(src_k,kmer_size, tmp_seq),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(tgt_k)),
	   binary_kmer_to_seq(tgt_k,kmer_size, tmp_seq), edge_type, edge_index);
  }

  db_node_add_labeled_edge(src_e,src_o,binary_kmer_get_last_nucleotide(tgt_k), edge_type, edge_index);

  if (DEBUG){

    printf("add edge %s -%c-> %s to edge type %d, and edge index %d\n",binary_kmer_to_seq(tgt_k,kmer_size,tmp_seq),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(src_k,kmer_size))),binary_kmer_to_seq(src_k,kmer_size, tmp_seq),  edge_type, edge_index);
  }

  db_node_add_labeled_edge(tgt_e,opposite_orientation(tgt_o),binary_kmer_get_last_nucleotide(binary_kmer_reverse_complement(src_k,kmer_size)), edge_type, edge_index );

  return true;
}



boolean db_node_edge_exist(dBNode * element,Nucleotide base,Orientation orientation, EdgeArrayType edge_type, int edge_index){

  //get the edge char for this specific person or pop:
  char edge = get_edge_copy(*element, edge_type, edge_index);


  edge >>= base;
  if (orientation == reverse){
    edge >>= 4;
  }
  
  edge &= 1;
  
  if (edge == 1){
    return true;
  }
  else{
    return false;
  }
}




void db_node_reset_edges(dBNode * node,EdgeArrayType edge_type, int edge_index){
  set_edges(node, edge_type, edge_index, 0);
}

void db_node_reset_edge(dBNode * node, Orientation orientation, Nucleotide nucleotide, EdgeArrayType edge_type, int edge_index){
  reset_one_edge(node, orientation, nucleotide, edge_type, edge_index);
}




boolean db_node_edges_reset(dBNode * node, EdgeArrayType edge_type, int edge_index){
  return get_edge_copy(*node,edge_type,edge_index) == 0;
}


boolean db_node_has_precisely_one_edge(dBNode * node, Orientation orientation, Nucleotide * nucleotide, EdgeArrayType edge_type, int edge_index){
  
  Nucleotide n;
  Edges edges = get_edge_copy(*node,edge_type,edge_index);
  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
 
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      *nucleotide = n;
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 1);
  
}


boolean db_node_has_precisely_one_edge_in_union_graph_over_all_people(dBNode * node, Orientation orientation, Nucleotide * nucleotide){
  
  Nucleotide n;
  Edges edges = get_union_of_edges(*node);
  short edges_count = 0;

  if (orientation == reverse){
    edges >>= 4;
  }
 
  
  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      *nucleotide = n;
      edges_count++;
    }
    
    edges >>= 1;    
  }
  
  return (edges_count == 1);
  
}


boolean db_node_is_blunt_end(dBNode * node, Orientation orientation, EdgeArrayType edge_type, int edge_index){
  
  Edges edges = get_edge_copy(*node, edge_type, edge_index);


  if (orientation == reverse){
    edges >>= 4;
  }
  
  edges &= 15; // AND with 00001111 so that we only look at the 4 least significant bits
  
  return edges == 0;
}

boolean db_node_check_status(dBNode * node, NodeStatus status){
  return node->status == status;
}
boolean db_node_check_status_not_pruned(dBNode * node){
  if ( db_node_check_status(node, none) || db_node_check_status(node,visited))
    {
      return true;
    }
  return false;
}


void db_node_set_status(dBNode * node,NodeStatus status){
  node->status = status;
}
void db_node_set_status_to_none(dBNode * node){
  node->status = none;
}

//assumes index 0 = NA12878, index 1 = NA12891, index 2 = NA12892
//semantics - call this when pruning node from person defined by index
void db_node_trio_aware_set_pruned_status(dBNode * node, int index)
{
  if (db_node_check_status(node,none) || db_node_check_status(node,visited))
    {

  	  if (index==0)//NA12878
	    {
	      db_node_set_status(node, pruned_from_NA12878);
	    }
	  else if (index==1) //NA12891
	    {
	      db_node_set_status(node, pruned_from_NA12891);
	    }
	  else if (index==2)//NA12892
	    {
	      db_node_set_status(node, pruned_from_NA12892);
	    }
    }

  if (db_node_check_status(node, pruned_from_NA12878))
    {
      if (index==0)//NA12878
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==1) //NA12891
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891);
	}
      else if (index==2)//NA12892
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12892);
	}
    }
  
  else if (db_node_check_status(node, pruned_from_NA12891))
    {
      if (index==0)//NA12878
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891);
	  
	}
      else if (index==1) //NA12891
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==2)//NA12892
	{
	  db_node_set_status(node, pruned_from_NA12891_and_NA12892);
	}
    }
  else if (db_node_check_status(node, pruned_from_NA12892))
    {
      if (index==0)//NA12878
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12892);
	  
	}
      else if (index==1) //NA12891
	{
	  db_node_set_status(node, pruned_from_NA12891_and_NA12892);
	}
      else if (index==2)//NA12892
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
    }
  else if (db_node_check_status(node, pruned_from_NA12878_and_NA12891))
    {
      if (index==0)//NA12878
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==1) //NA12891
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==2)//NA12892
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891_and_NA12892);
	}
    }
  
  else if (db_node_check_status(node, pruned_from_NA12878_and_NA12892))
    {
      if (index==0)//NA12878
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==1) //NA12891
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891_and_NA12892);
	}
      else if (index==2)//NA12892
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
    }
  else if (db_node_check_status(node, pruned_from_NA12891_and_NA12892))
    {
      if (index==0)//NA12878
	{
	  db_node_set_status(node, pruned_from_NA12878_and_NA12891_and_NA12892);
	}
      else if (index==1) //NA12891
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
      else if (index==2)//NA12892
	{
	  printf("WARNING. Pruning a node that is already pruned");
	}
    }
}



boolean db_node_is_this_node_in_this_person_or_populations_graph(dBNode* node, EdgeArrayType type, int index)
{

  if (node ==NULL)
    {
      return false;
    }

  Edges edge_for_this_person_or_pop = get_edge_copy(*node, type, index);

  if (edge_for_this_person_or_pop == 0)
    {

      return false;
    }
  else
    {
      return true;
    }
 
}


