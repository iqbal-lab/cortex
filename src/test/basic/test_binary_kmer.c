/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */

/*
 * test_binary_kmer.c
 *
 */

#include <CUnit.h>
#include <Basic.h>
#include <binary_kmer.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <test_binary_kmer.h>


void test_that_bitfield_really_is_64bits()
{
  CU_ASSERT(sizeof(bitfield_of_64bits)==8);
}


void test_binary_kmer_assignment_operator()
{
  int i;

  bitfield_of_64bits  array[10] = {(bitfield_of_64bits) 1,(bitfield_of_64bits)2,(bitfield_of_64bits)7,(bitfield_of_64bits)25,
				   (bitfield_of_64bits)25656464646,(bitfield_of_64bits)43,(bitfield_of_64bits)4444444,(bitfield_of_64bits)101010,(bitfield_of_64bits)999,
				   (bitfield_of_64bits)32};

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>10)
    {
      printf("Had assumed you would never try going higher than 10 long longs in a BinaryKmer. Fix the test - very simple fix");
      exit(1);
    }

  //set up your test kmer
  BinaryKmer test_kmer;
  binary_kmer_initialise_to_zero(&test_kmer);

  for (i=0 ; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      test_kmer[i]=array[i];

    }
  
  BinaryKmer assignee;
  binary_kmer_assignment_operator(assignee, test_kmer);

  for (i=0 ; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      CU_ASSERT(assignee[i]==test_kmer[i]);
    }



}

void test_binary_kmer_comparison_operator()
{

  BinaryKmer bk1;
  BinaryKmer bk2;
  BinaryKmer bk3;
  BinaryKmer bk4;
  BinaryKmer bk5;
  BinaryKmer bk6;
  BinaryKmer bk7;
  BinaryKmer bk8;

  int i;
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      bk1[i]=0;
      bk2[i]=1;
      bk3[i]=0;
      bk4[i]=i+2;
      bk5[i]=i+1;
      bk6[i]=~0;
      bk7[i]=1<<i;
      bk8[i]=1<<i;
    }
  CU_ASSERT(binary_kmer_comparison_operator(bk1,bk2)==false);
  CU_ASSERT(binary_kmer_comparison_operator(bk1,bk3)==true);
  CU_ASSERT(binary_kmer_comparison_operator(bk3,bk4)==false);
  CU_ASSERT(binary_kmer_comparison_operator(bk4,bk5)==false);
  CU_ASSERT(binary_kmer_comparison_operator(bk2,bk6)==false);
  CU_ASSERT(binary_kmer_comparison_operator(bk7,bk8)==true);

}

void test_binary_kmer_less_than()
{
  
  BinaryKmer bk1;
  BinaryKmer bk2;

  int i;
  
  bitfield_of_64bits j;
  short kmer_size;
  
  for (kmer_size=31; kmer_size < NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32; kmer_size=kmer_size+32)
    {
      int number_of_bitfields_fully_used = kmer_size/32;
      int number_of_bits_in_most_sig_bitfield = 2* (kmer_size-(32*number_of_bitfields_fully_used));

      for (j=0; j<10; j++)
	{
	  
	  for (i=NUMBER_OF_BITFIELDS_IN_BINARY_KMER-number_of_bitfields_fully_used-1; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	    {
	      bk1[i]=j;
	      bk2[i]=j;
	    }
	  
	  //Now make a single difference between the two
	  for (i=NUMBER_OF_BITFIELDS_IN_BINARY_KMER-number_of_bitfields_fully_used-1; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	    {
	      bk1[i]=(bitfield_of_64bits) 10;
	      bk2[i]=(bitfield_of_64bits) 11;
	      CU_ASSERT(binary_kmer_less_than(bk1,bk2, kmer_size)==true);
	
	if (binary_kmer_less_than(bk1,bk2, kmer_size)!=true)
	  {
	    printf("Breaks at kmer = %d, i=%d, j=%llu\n",kmer_size, i,j);
	    printf("num of bitfields fully used is %d, and number of bits in most sig is %d\n", number_of_bitfields_fully_used, number_of_bits_in_most_sig_bitfield);
	    exit(1);
	  }
	
	CU_ASSERT(binary_kmer_less_than(bk2,bk1, kmer_size)==false);
	
	//reset values, so now both arrays are identical
	bk1[i]=j;
	bk2[i]=j;
	
	    }
	  
	}
    }
  
}

void test_binary_kmer_right_shift()
{
  BinaryKmer test;
  int i;


  // First a simple test. Start with a 1 at the far left, and keep shifting it one bit to the right
  
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      test[i]=0;
    }

  bitfield_of_64bits init = ((bitfield_of_64bits) 1) <<63;
  test[0]= init ;//put a 1 at the far left of the left-most bitfield

  

  for (i=1; i<64; i++)
    {
      binary_kmer_right_shift(&test, 1);
      CU_ASSERT(test[0]==(init>>i));
      int j;
      for (j=1; j< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; j++)
	{
	  CU_ASSERT(test[j]==0);
	}
    }


  //so we should now have all bitfields=0, exceptfar left, which is 1
  CU_ASSERT(test[0]==1);
  int j;
  for (j=1; j< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; j++)
    {
      CU_ASSERT(test[j]==0);
    }


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      binary_kmer_right_shift(&test, 1);
      CU_ASSERT(test[0]==0);
      CU_ASSERT(test[1]==init);

      for (i=1; i<64; i++)
	{
	  binary_kmer_right_shift(&test, 1);
	  
	  CU_ASSERT(test[0]==0);
	  CU_ASSERT(test[1]==init>>i);
	  int j;
	  for (j=2; j< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; j++)
	    {
	      test[j]=0;
	    }
	}

    }


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>=2)
    {
      //cleanup
      for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	{
	  test[i]=0;
	}
      
      
      test[0]= ((bitfield_of_64bits) 1) <<63;
      
      binary_kmer_right_shift(&test,2);
      CU_ASSERT(test[0]== ((bitfield_of_64bits) 1)<<61);
      CU_ASSERT(test[1]==0);
      binary_kmer_right_shift(&test,20);
      CU_ASSERT(test[0]== ((bitfield_of_64bits) 1)<<41);
      CU_ASSERT(test[1]==0);
      binary_kmer_right_shift(&test,50);
      CU_ASSERT(test[0]==0);
      CU_ASSERT(test[1]== ((bitfield_of_64bits) 1)<<55) ;
      binary_kmer_right_shift(&test,50);
      CU_ASSERT(test[0]==0);
      CU_ASSERT(test[1]== ((bitfield_of_64bits) 1)<<5) ;
    }


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>=3)
    {
      for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
        {
          test[i]=0;
	}

      test[0]=~0;
      binary_kmer_right_shift(&test,1);
      CU_ASSERT(test[0]==((((bitfield_of_64bits) 1)<<63) -1) );
      CU_ASSERT(test[1]== ((bitfield_of_64bits) 1)<<63) ;
      binary_kmer_right_shift(&test,1);
      CU_ASSERT(test[0]==(( ((bitfield_of_64bits) 1)<<62) -1) );
      CU_ASSERT(test[1]==((bitfield_of_64bits) 3)<<62);
      binary_kmer_right_shift(&test,63);
      CU_ASSERT(test[0]==0);
      CU_ASSERT(test[1]==(( ((bitfield_of_64bits) 1)<<63) -1) );
      CU_ASSERT(test[2]==((bitfield_of_64bits) 1)<<63);

    }

}

void test_binary_kmer_left_shift()
{


  BinaryKmer test;
  int i;

  // First a simple test. Start with a 1 at the far right, and keep shifting it one bit to the left

  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      test[i]=0;
    }
  test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=1;

  int kmer_size = NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32-1;
  int number_of_bitfields_fully_used = kmer_size/32;
  int number_of_bits_in_most_sig_bitfield = 2* (kmer_size-(32*number_of_bitfields_fully_used));


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==1)
    {
      for (i=0; i<61; i++)
	{
	  binary_kmer_left_shift(&test, 1, kmer_size);
	  
	  //printf("Have shifted 1 left by %d, and furthest rh bitfield is %lld\n", i+1, test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]);
	  
	  CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==( ((bitfield_of_64bits) 1)<<(i+1)));
	  int j;
	  
	  for (j=0; j<NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; j++)
	    {
	      CU_ASSERT(test[j]==0);
	    }
	  
	}
      binary_kmer_left_shift(&test, 1, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==0)
      binary_kmer_left_shift(&test, 1, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==0)

    }
  else
    {
      for (i=0; i<63; i++)
	{
	  binary_kmer_left_shift(&test, 1, kmer_size);
	  
	  //printf("Have shifted 1 left by %d, and furthest rh bitfield is %lld\n", i+1, test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]);
	  
	  CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==( ((bitfield_of_64bits) 1)<<(i+1)));
	  int j;
	  
	  for (j=0; j<NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; j++)
	    {
	      CU_ASSERT(test[j]==0);
	    }
	  
	}
    }

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      binary_kmer_left_shift(&test, 1, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==0);

      if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER==2)
	{
	  for (i=0; i<61; i++)//ie go to the end of the kmer, last 2 bits are masked to zero
	    {
	      binary_kmer_left_shift(&test, 1, kmer_size);
	      
	      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==0);
	      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]==( ((bitfield_of_64bits) 1)<<(i+1)));

	      int j;
	      for (j=0; j< NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2; j++)
		{
		  test[j]=0;
		}
	    }
	}
	else
	  {
	  for (i=0; i<63; i++)
	    {
	      binary_kmer_left_shift(&test, 1, kmer_size);
	      
	      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==0);
	      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]==( ((bitfield_of_64bits) 1)<<(i+1)));

	      int j;
	      for (j=0; j< NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2; j++)
		{
		  test[j]=0;
		}
	    }

	  }
      
    }


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>=2)
    {
      for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	{
	  test[i]=0;
	}
      

      test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=1;
      
      binary_kmer_left_shift(&test,2, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] == ((bitfield_of_64bits) 1)<<2);
      binary_kmer_left_shift(&test,20,kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]== ((bitfield_of_64bits) 1)<<22);
      binary_kmer_left_shift(&test,50, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==0);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]== ((bitfield_of_64bits) 1)<<8) ;
      binary_kmer_left_shift(&test,50, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==0);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]== ((bitfield_of_64bits) 1)<<58) ;
    }
  

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>=3)
    {
      for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	{
          test[i]=0;
	}

      test[0]=~0;
      binary_kmer_left_shift(&test,63, kmer_size);
      binary_kmer_left_shift(&test,1, kmer_size);

      for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	{
          CU_ASSERT(test[i]==0);
	  if (test[i] !=0)
	    {
	      printf("i is %d when test[i] is not zero - it is %llu\n", i, test[0]);
	      exit(1);
	    }
	}
      
      test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=~0;
      binary_kmer_left_shift(&test,1, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]  ==~1 );
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]==1) ;
      binary_kmer_left_shift(&test,1, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==~3 );
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]==3)  ;
      binary_kmer_left_shift(&test,63, kmer_size);
      binary_kmer_left_shift(&test,1, kmer_size);

      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]  == 0);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]== ~3);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-3]== 3);

    }


  // Now run some checks with different kmer_sizes, so that we for example do not use all of the bitfields completely.
  
  for (kmer_size=1; kmer_size<=NUMBER_OF_BITFIELDS_IN_BINARY_KMER*32; kmer_size=kmer_size+1)
    {
      //cleanup to start with
      for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
        {
          test[i]=0;
	}

      int number_of_bitfields_fully_used = kmer_size/32;
      int number_of_bits_in_most_sig_bitfield = 2* (kmer_size-(32*number_of_bitfields_fully_used));

      //simple example - 1s in every bit up to size of kmer
      //Set all the bitfields in test except the leftmost that holds any kmer bits


      for (i=NUMBER_OF_BITFIELDS_IN_BINARY_KMER - number_of_bitfields_fully_used-1; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	{
	  test[i]=~0;
	}

      bitfield_of_64bits init = (((bitfield_of_64bits)1)<<number_of_bits_in_most_sig_bitfield) -1; //all 1's as far as can be fitted within the kmer

      if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER - number_of_bitfields_fully_used>=1)
	{
	  //set the final bitfield of our test BinaryKmer
	  test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER - number_of_bitfields_fully_used-1]=  init;
	      
	  //at first, every time you shift left, you still have 1's in the leftmost places, as you're pushing 1's in,
	  // and the 1's that go beyond the kmer-length are masked off to zero.
	  for (i=1; i<= number_of_bitfields_fully_used*64; i++)
	    {
	      binary_kmer_left_shift(&test,1, kmer_size);
	      
	      int j;
	      for (j=0; j< NUMBER_OF_BITFIELDS_IN_BINARY_KMER - number_of_bitfields_fully_used-1; j++)
		{
		  CU_ASSERT(test[j]==0); //zero's to the left of the bits used to hold the kmer
		}
	      
	      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER - number_of_bitfields_fully_used-1]== init);
	    }

	  //now all that is left is the 1's in the leftmost of the bitfields that held any kmer bits, and each shift left will push one off the end
	  for (i=0; i<number_of_bits_in_most_sig_bitfield; i++)
	    {
	      binary_kmer_left_shift(&test,1, kmer_size);
	      int j;
              for (j=0; j< NUMBER_OF_BITFIELDS_IN_BINARY_KMER - number_of_bitfields_fully_used-1; j++)
                {
                  CU_ASSERT(test[j]==0); //zero's to the left of the bits used to hold the kmer
                }
	      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER - number_of_bitfields_fully_used-1]== init - ((((bitfield_of_64bits)1)<<(i+1))-1) );
	    }


	}


    }


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>=3)
    {
      //cleanup to start afresh
      for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	{
	  test[i]=0;
	}
      
      kmer_size=67;
      test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=32; // 000...000100000
      test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]=1;
      test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-3]=31; // 00..011111

      binary_kmer_left_shift(&test,1, kmer_size);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]==64);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-2]==2);
      CU_ASSERT(test[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-3]==62);//00...0111110
      
      for (i=0; i<NUMBER_OF_BITFIELDS_IN_BINARY_KMER-3; i++)
	{
	  CU_ASSERT(test[i]==0);
	}
      
      
    }

}


void test_seq_to_binary_kmer_and_binary_kmer_to_seq(){
  
  BinaryKmer kmer;
  binary_kmer_initialise_to_zero(&kmer);
  char seq[7];
  
  seq_to_binary_kmer("ATCGCGC",7, &kmer);
 
  CU_ASSERT_STRING_EQUAL("ATCGCGC",binary_kmer_to_seq(&kmer,7,seq));

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
      char seq2[37];
      binary_kmer_initialise_to_zero(&kmer);
      seq_to_binary_kmer("CATCAGTGGGACATAAACCACACAGATGACACACACA",37, &kmer);
      CU_ASSERT_STRING_EQUAL("CATCAGTGGGACATAAACCACACAGATGACACACACA",binary_kmer_to_seq(&kmer,37,seq2));


      char seq3[63];
      binary_kmer_initialise_to_zero(&kmer);
      seq_to_binary_kmer("CATCAGTGGGACATAAACCACACAGATGACACACACACATCAGTGGGACATAAACCACACAGA",63, &kmer);
      CU_ASSERT_STRING_EQUAL("CATCAGTGGGACATAAACCACACAGATGACACACACACATCAGTGGGACATAAACCACACAGA",binary_kmer_to_seq(&kmer,63,seq3));
    }

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>2)
    {
      char seq2[80];
      binary_kmer_initialise_to_zero(&kmer);
      seq_to_binary_kmer("CATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTCATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGT",80, &kmer);
      CU_ASSERT_STRING_EQUAL("CATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTCATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGT",binary_kmer_to_seq(&kmer,80,seq2));

      char seq3[95];
      binary_kmer_initialise_to_zero(&kmer);
      seq_to_binary_kmer("CATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTCATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTGGGGGGGGTTTTCAC",95, &kmer);
      CU_ASSERT_STRING_EQUAL("CATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTCATCAGGGGTCAGTCAGTCACACAAAAACGTCAGTCGTGTGGGGGGGGTTTTCAC",binary_kmer_to_seq(&kmer,95,seq3));

    }



  
}

void test_binary_kmer_reverse_complement(){
  
  BinaryKmer kmer, kmer_reverse;
  char seq7[7];
  char seq1[1];
  char seq31[31];
  char seq33[33];
  char seq37[37];
  char seq57[57];
  char seq63[63];
  char seq65[65];
  char seq81[81];
  char seq97[97];
  char seq127[127];
  char seq159[159];
  char seq255[255];


  //test with various different kmer sizes


  seq_to_binary_kmer("ATCGCGC",7, &kmer);
  binary_kmer_reverse_complement(&kmer, 7, &kmer_reverse);  
  CU_ASSERT_STRING_EQUAL("GCGCGAT",binary_kmer_to_seq(&kmer_reverse,7,seq7));

  seq_to_binary_kmer("A",1, &kmer);
  binary_kmer_reverse_complement(&kmer, 1, &kmer_reverse);
  CU_ASSERT_STRING_EQUAL("T",binary_kmer_to_seq(&kmer_reverse,1,seq1));


  seq_to_binary_kmer("GGCCCCGCCCCGCCCCGCCCCGCCCCGCCCC",31, &kmer);
  binary_kmer_reverse_complement(&kmer, 31, &kmer_reverse);
  CU_ASSERT_STRING_EQUAL("GGGGCGGGGCGGGGCGGGGCGGGGCGGGGCC", binary_kmer_to_seq(&kmer_reverse,31,seq31));


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>1)
    {
  
      seq_to_binary_kmer("AAACGTAACGTAACGTAACGTTTTTTCATGGCA",33, &kmer);
      binary_kmer_reverse_complement(&kmer, 33, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("TGCCATGAAAAAACGTTACGTTACGTTACGTTT", binary_kmer_to_seq(&kmer_reverse,33,seq33));
  
      seq_to_binary_kmer("AAACGTAACGTAACGTAACGTTTTTTCATGGCAACGT",37, &kmer);
      binary_kmer_reverse_complement(&kmer, 37, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("ACGTTGCCATGAAAAAACGTTACGTTACGTTACGTTT", binary_kmer_to_seq(&kmer_reverse,37,seq37));
  
      seq_to_binary_kmer("ACGTTGCCATGAAAAAACGTTACGTTACGTTACGTTTAAAAAAAAAACCCCCCCCCC",57, &kmer);
      binary_kmer_reverse_complement(&kmer, 57, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("GGGGGGGGGGTTTTTTTTTTAAACGTAACGTAACGTAACGTTTTTTCATGGCAACGT", binary_kmer_to_seq(&kmer_reverse,57,seq57));

      seq_to_binary_kmer("ACGTTGCCATGAAAAAACGTTACGTTACGTTACGTTTAAAAAAAAAACCCCCCCCCCGGGTAC",63, &kmer);
      binary_kmer_reverse_complement(&kmer, 63, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("GTACCCGGGGGGGGGGTTTTTTTTTTAAACGTAACGTAACGTAACGTTTTTTCATGGCAACGT", binary_kmer_to_seq(&kmer_reverse,63,seq63));
    }

  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>2)
    {

      seq_to_binary_kmer("TACGTTGCCATGAAAAAACGTTACGTTACGTTACGTTTAAAAAAAAAACCCCCCCCCCGGGTACC",65, &kmer);
      binary_kmer_reverse_complement(&kmer, 65, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("GGTACCCGGGGGGGGGGTTTTTTTTTTAAACGTAACGTAACGTAACGTTTTTTCATGGCAACGTA", binary_kmer_to_seq(&kmer_reverse,65,seq65));
      
      seq_to_binary_kmer("AACGTGTGTGCCCATACAGGAACGTGTGTGCCCATACAGGAACGTGTGTGCCCATACAGGAACGTGTGTGCCCATACAGGG", 81, &kmer);
      binary_kmer_reverse_complement(&kmer, 81, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("CCCTGTATGGGCACACACGTTCCTGTATGGGCACACACGTTCCTGTATGGGCACACACGTTCCTGTATGGGCACACACGTT", binary_kmer_to_seq(&kmer_reverse,81,seq81));
    }
  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>3)
    {
      seq_to_binary_kmer("CCCCCCCCCTATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCA",97, &kmer);
      binary_kmer_reverse_complement(&kmer, 97, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("TGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATAGGGGGGGGG", binary_kmer_to_seq(&kmer_reverse,97,seq97));

      seq_to_binary_kmer("CCCCCCCCTATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCAC", 127, &kmer);
      binary_kmer_reverse_complement(&kmer, 127, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("GTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATAGGGGGGGG", binary_kmer_to_seq(&kmer_reverse,127,seq127));

    }




  //if you think all these tests are paranoid - this next one found a bug in how the far left bitfield was treated.


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>4)
    {

      seq_to_binary_kmer("TGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATAGGGGGGGGG", 159, &kmer);
      binary_kmer_reverse_complement(&kmer, 159, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("CCCCCCCCCTATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACA",  
			     binary_kmer_to_seq(&kmer_reverse,159,seq159));
    }


  if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>7)
    {
      
      seq_to_binary_kmer("CCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCCCATATGTGCGGGGG", 255, &kmer);
      binary_kmer_reverse_complement(&kmer, 255, &kmer_reverse);
      CU_ASSERT_STRING_EQUAL("CCCCCGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGGGCACATATGG",  binary_kmer_to_seq(&kmer_reverse,255,seq255));
      
  
    }  



  
}


void test_seq_reverse_complement(){
  
  char out[100];
  char * seq = "AAAAAA";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq,6,out),"TTTTTT");
  
  seq = "ATAAAA";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq,6,out),"TTTTAT");
  
  seq = "CGATAAAA";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq,8,out),"TTTTATCG");
 
  seq = "CGATAAAAGG";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq,10,out),"CCTTTTATCG");
   
  seq = "";
  CU_ASSERT_STRING_EQUAL(seq_reverse_complement(seq,0,out),"");
  
}



void test_binary_kmer_nucleotide_iterator(){
  
  int count = 0;
  
  void check_nucleotides(Nucleotide base){
    
    if (base == Adenine || base == Guanine || base == Cytosine || base == Thymine){
      count ++;
    }
    else{
      count --;
    }
  }
  
  nucleotide_iterator(check_nucleotides);
  
  CU_ASSERT_EQUAL(count,4);
}



void test_get_sliding_windows_from_sequence(){

    
   //----------------------------------
    KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
    if (windows == NULL){
      fputs("Out of memory trying to allocate a KmerSlidingWindowSet",stderr);
      exit(1);
    } 

    binary_kmer_alloc_kmers_set(windows,20,30);
    
    CU_ASSERT_EQUAL(windows->max_nwindows,20);

    //----------------------------------


    // ****************************************************
    // First tests do not involve any homopolymer breaking:
    // ****************************************************

  
    char * seq  = "AAAAANNTTTTGGGG";
    char qual0[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    char kmer_seq3[3];

    int nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, false, 0);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 3);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 6);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),3,kmer_seq3),"TTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),3,kmer_seq3),"TGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[4]),3,kmer_seq3),"GGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[5]),3,kmer_seq3),"GGG");

    CU_ASSERT_EQUAL(nkmers1, 9);
    
    char qual1[15] = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 10, 20, 20, 20};

    int nkmers2 = get_sliding_windows_from_sequence(seq,qual1,strlen(seq),15,3,windows,20,20, false,0);

    CU_ASSERT_EQUAL(windows->nwindows,3);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 3);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),3,kmer_seq3),"AAA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 2);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),3,kmer_seq3),"TTT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),3,kmer_seq3),"TTT");

    CU_ASSERT_EQUAL((windows->window[2]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[0]),3,kmer_seq3),"GGG");

    CU_ASSERT_EQUAL(nkmers2, 6);

    char qual2[15] = { 5, 10, 20, 20, 20, 20, 0, 0, 0, 20, 20, 30, 20, 20, 20};

    int nkmers3 = get_sliding_windows_from_sequence(seq,qual2,strlen(seq),15,3,windows,20,20, false,0);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"AAA");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),3,kmer_seq3),"TTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),3,kmer_seq3),"TGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),3,kmer_seq3),"GGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),3,kmer_seq3),"GGG");
    
    CU_ASSERT_EQUAL(nkmers3, 5);

    char qual3[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int nkmers4 = get_sliding_windows_from_sequence(seq,qual3,strlen(seq),15,3,windows,20,20, false,0);

    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers4, 0);

    char kmer_seq5[5];
    int nkmers5 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,5,windows,20,20, false,0);
    
    CU_ASSERT_EQUAL(windows->nwindows,2);
  
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),5,kmer_seq5),"AAAAA");

  
    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),5,kmer_seq5),"TTTTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),5,kmer_seq5),"TTTGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),5,kmer_seq5),"TTGGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),5,kmer_seq5),"TGGGG");

    CU_ASSERT_EQUAL(nkmers5, 5);
    


    // Now try some examples with large kmers

    
    if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>2)
      {
	
	char* seq_contains_only_one_good_33mer="CCCCCCCCTANTGGGCACATATGGGCACATATGGNGCACATATGGGCACATATGGGCACATATGGNGCACATATGGNNNNNNNNNNGCACATATGGGNCANCATATGNNGCACATATGGGCACATATGGGCACATATGGGCANC";
        char qual_144_zeroes[144] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	char qual_144_30s[144] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30};
	char qual_all_30s_except_for_the_one_33mer[144] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30, 30};


	char kmer_seq33[33];
	int quality_cutoff=0;
	short kmer_size=33;
	int nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_144_zeroes,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, false,0);

	CU_ASSERT_EQUAL(windows->nwindows,1);
	CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
	CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),33,kmer_seq33),"GCACATATGGGCACATATGGGCACATATGGGCA");
	CU_ASSERT(nkmers_33==1);

	quality_cutoff=20;
	nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_144_30s,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, false,0);
	CU_ASSERT_EQUAL(windows->nwindows,1);
	CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
	CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),33,kmer_seq33),"GCACATATGGGCACATATGGGCACATATGGGCA");
	CU_ASSERT(nkmers_33==1);

	nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_all_30s_except_for_the_one_33mer,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, false,0);
	CU_ASSERT_EQUAL(windows->nwindows,0);
	CU_ASSERT(nkmers_33==0);


      }



    binary_kmer_free_kmers_set(&windows);
    
    CU_ASSERT(windows == NULL);


    
}


void test_breaking_homopolymers_in_get_sliding_windows ()
{
   //----------------------------------
    KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
    if (windows == NULL){
      fputs("Out of memory trying to allocate a KmerSlidingWindowSet",stderr);
      exit(1);
    } 

    binary_kmer_alloc_kmers_set(windows,20,30);
    
    CU_ASSERT_EQUAL(windows->max_nwindows,20);

    //----------------------------------


    char * seq  = "AAAAANNTCAGAT";
    char qual0[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    char kmer_seq3[3];

    //break at homopolymers of length 4. ie you should never load any homopolymer of length >3, and when you reach one of length >=4, you don't start your next window
    //until after the end of the homopolymer
    int nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 4);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"AAA");
 
    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),3,kmer_seq3),"TCA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),3,kmer_seq3),"CAG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),3,kmer_seq3),"AGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),3,kmer_seq3),"GAT");

    CU_ASSERT_EQUAL(nkmers1, 5);

    //now try length 3, which is also the kmer_length. Should only get one window
    nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 3);

    CU_ASSERT_EQUAL(windows->nwindows,1);
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"TCA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),3,kmer_seq3),"CAG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),3,kmer_seq3),"AGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[3]),3,kmer_seq3),"GAT");

    CU_ASSERT_EQUAL(nkmers1, 4);

    //does it break if you enter hom cutoff of 0 or 1 or 2?

    //cutoff is 2
    nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 2);

    CU_ASSERT_EQUAL(windows->nwindows,1);
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 4);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),3,kmer_seq3),"TCA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),3,kmer_seq3),"CAG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),3,kmer_seq3),"AGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[3]),3,kmer_seq3),"GAT");

    CU_ASSERT_EQUAL(nkmers1, 4);

    //cutoff is 1
    nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 1);

    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers1, 0);

    //cutoff is zero
    nkmers1 = get_sliding_windows_from_sequence(seq,qual0,strlen(seq),0,3,windows,20,20, true, 0);

    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers1, 0);



    // another example:

    seq  = "AAAAANNTCTCCCCCTAGATGGGGGGGGGGGGGCCACCCNCCCCGTGATAT";
    char qual_other[51] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    char kmer_seq7[7];

    //first of all - not breaking homopolymers, only the N's cause problems:
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, false, 0);

    CU_ASSERT_EQUAL(windows->nwindows,2);

    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 26);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"TCTCCCC");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),7,kmer_seq7),"CTCCCCC");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[2]),7,kmer_seq7),"TCCCCCT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[3]),7,kmer_seq7),"CCCCCTA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[4]),7,kmer_seq7),"CCCCTAG");
    // ... not doing them all
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[25]),7,kmer_seq7),"GCCACCC");



    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 5);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),7,kmer_seq3),"CCCCGTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),7,kmer_seq3),"CCCGTGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),7,kmer_seq3),"CCGTGAT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[3]),7,kmer_seq3),"CGTGATA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[4]),7,kmer_seq3),"GTGATAT");

    CU_ASSERT_EQUAL(nkmers1, 31);


    //then, breaking homopolymers - no homopolymer of length 1 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 1);
    CU_ASSERT_EQUAL(windows->nwindows,0);
    CU_ASSERT_EQUAL(nkmers1, 0);
    //then, breaking homopolymers - no homopolymer of length 2 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 2);
    CU_ASSERT_EQUAL(windows->nwindows,1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"GTGATAT");
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_EQUAL(nkmers1, 1);
    //then, breaking homopolymers - no homopolymer of length 3 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 3);
    CU_ASSERT_EQUAL(windows->nwindows,2);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"TAGATGG");
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),7,kmer_seq7),"GTGATAT");
    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 1);
    CU_ASSERT_EQUAL(nkmers1, 2);
    //then, breaking homopolymers - no homopolymer of length 4 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 4);
    CU_ASSERT_EQUAL(windows->nwindows,2);
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 2);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"TAGATGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[1]),7,kmer_seq7),"AGATGGG");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),7,kmer_seq7),"GTGATAT");
    CU_ASSERT_EQUAL(nkmers1, 3);

    //"AAAAANNTCTCCCCCTAGATGGGGGGGGGGGGGCCACCCNCCCCGTGATAT"
    //then, breaking homopolymers - no homopolymer of length 5 or greater
    nkmers1 = get_sliding_windows_from_sequence(seq,qual_other,strlen(seq),0,7,windows,40,40, true, 5);
    CU_ASSERT_EQUAL(windows->nwindows,3);
    CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),7,kmer_seq7),"TCTCCCC");

    CU_ASSERT_EQUAL((windows->window[1]).nkmers, 3);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[0]),7,kmer_seq7),"TAGATGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[1]),7,kmer_seq7),"AGATGGG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[1].kmer[2]),7,kmer_seq7),"GATGGGG");

    CU_ASSERT_EQUAL((windows->window[2]).nkmers, 5);
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[0]),7,kmer_seq7),"CCCCGTG");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[1]),7,kmer_seq7),"CCCGTGA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[2]),7,kmer_seq7),"CCGTGAT");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[3]),7,kmer_seq7),"CGTGATA");
    CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[2].kmer[4]),7,kmer_seq7),"GTGATAT");

    CU_ASSERT_EQUAL(nkmers1, 9);

    



    



    // Now try some examples with large kmers

    
    if (NUMBER_OF_BITFIELDS_IN_BINARY_KMER>2)
      {
	
	char* seq_contains_only_one_good_33mer="CCCCCCCCTANTGGGCACATATGGGCACATATGGNGCACATATGGGCACATATGGGCACATATGGNGCACATATGGNNNNNNNNNNGCACATATGGGNCANCATATGNNGCACATATGGGCACATATGGGCACATATGGGCANC";
        char qual_144_zeroes[144] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	char kmer_seq33[33];
	int quality_cutoff=0;
	short kmer_size=33;
	//let's put a homopolymer cutoff of 4 - should have no effect, as the one good 33mer contains none longer than 3
	int nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_144_zeroes,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, true,4);

	CU_ASSERT_EQUAL(windows->nwindows,1);
	CU_ASSERT_EQUAL((windows->window[0]).nkmers, 1);
	CU_ASSERT_STRING_EQUAL(binary_kmer_to_seq(&(windows->window[0].kmer[0]),33,kmer_seq33),"GCACATATGGGCACATATGGGCACATATGGGCA");
	CU_ASSERT(nkmers_33==1);


	//let's put a homopolymer cutoff of 3 - should find nothing
	nkmers_33 = get_sliding_windows_from_sequence(seq_contains_only_one_good_33mer,qual_144_zeroes,strlen(seq_contains_only_one_good_33mer),quality_cutoff,kmer_size,windows,40,40, true,3);

	CU_ASSERT_EQUAL(windows->nwindows,0);
	CU_ASSERT(nkmers_33==0);



      }



    binary_kmer_free_kmers_set(&windows);
    
    CU_ASSERT(windows == NULL);

  
}
