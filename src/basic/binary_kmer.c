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
   binary_kmer.c - routines to 
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <binary_kmer.h>
#include <event_encoding.h>
#include <global.h>
#include <string.h>

void binary_kmer_initialise_to_zero(BinaryKmer* bkmer)
{
  int i;
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      ((*bkmer)[i])=(bitfield_of_64bits) 0;
    }
}


void binary_kmer_assignment_operator(BinaryKmer left, BinaryKmer right)
{
  int i;

  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      left[i]=right[i];
    }
}


void binary_kmer_set_all_bitfields(BinaryKmer assignee, bitfield_of_64bits val)
{
  int i;
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      assignee[i]=val;
    }
}


//returns true if they are the same
boolean binary_kmer_comparison_operator(const BinaryKmer const left, const BinaryKmer const right)
{

  boolean they_are_the_same=true;

  int i;
  for (i=0; i< NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      if (left[i]!=right[i])
	{
	  they_are_the_same=false;
	  break;                             //sorry Mario, I know you hate breaks
	}
      
    }

  return they_are_the_same;
}


//TODO - this wrongly says left<right when they are the same! IS really a <= operator
boolean binary_kmer_less_than(const BinaryKmer const left, const BinaryKmer const right, short kmer_size)
{
  boolean left_is_less_than_right=false;

  //need the following to work out which bits to ignore
  int number_of_bitfields_fully_used = kmer_size/32;
  //int number_of_bits_in_most_sig_bitfield = 2* (kmer_size-(32*number_of_bitfields_fully_used));

  int i;


  //start at most significant end
  // this would break if we had number_of_bitfields_fully_used==NUMBER_OF_BITFIELDS_IN_BINARY_KMER. But we can never have that as k is always odd.
  for (i=NUMBER_OF_BITFIELDS_IN_BINARY_KMER-number_of_bitfields_fully_used-1; i<NUMBER_OF_BITFIELDS_IN_BINARY_KMER ; i++)
    {
      if (left[i]<right[i])
	{
	  left_is_less_than_right=true;
	  break;                             
	}
      else if (left[i]>right[i])
	{
	  left_is_less_than_right=false;
	  break;
	}

      
    }

  return left_is_less_than_right;
  
}



//implicit in this is the idea that you shift left, and mask to 0 the bits that fall off the left hand end
void binary_kmer_right_shift_one_base(BinaryKmer kmer)
{
  int i;
  for(i = NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; i > 0; i--)
  {
    kmer[i] >>= 2;
    kmer[i] |= (kmer[i-1] << 62); // & 0x3
  }

  kmer[0] >>= 2;
}

void binary_kmer_left_shift_one_base(BinaryKmer kmer, short kmer_size)
{
  int top_word = NUMBER_OF_BITFIELDS_IN_BINARY_KMER - (kmer_size+31)/32;

  int i;
  for(i = top_word; i < NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1; i++)
  {
    kmer[i] <<= 2;
    kmer[i] |= (kmer[i+1] >> 62); // & 0x3
  }

  kmer[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] <<= 2;

  // Mask top word
  short top_bits = 2 * (kmer_size % 32); // bits in top word
  kmer[top_word] &= (~(unsigned long long)0 >> (64 - top_bits));
}


void binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(BinaryKmer* bkmer, Nucleotide n, short kmer_size)
{

  //shift left by one base,
  binary_kmer_left_shift_one_base(*bkmer, kmer_size);

  // add new base at right hand end 
  (*bkmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] |= n;

}




char reverse_char_nucleotide(char c)
{
  switch (c)
    {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    case 'a':
      return 't';
    case 'c':
      return 'g';
    case 'g':
      return 'c';
    case 't':
      return 'a';
    default:
      printf("Non-existent nucleotide %c\n",c);
      assert(0);
      return 'N';
      //return Adenine; 
    }
}


//length is the length in number of bases; the char* should have one MORE base than that allocated, to hold '\0'
char * seq_reverse_complement(char * in, int length, char * out){

  int k;
  for(k=0;k<length;k++){
    out[k] = reverse_char_nucleotide(in[length-k-1]);
  }
  out[length] = '\0';
  return out;
}


/*
Nucleotide reverse_binary_nucleotide(Nucleotide n)
{
  switch (n)
    {
    case Adenine:
      return Thymine;
    case Cytosine:
      return Guanine;
    case Guanine:
      return Cytosine;
    case Thymine:
      return Adenine;
    default:
      printf("Calling reverse_binary_nucleotide on non-existent nucleotide %i\n",n);
      exit(1);
    }
}
*/


char * nucleotides_to_string(Nucleotide * nucleotides, int length, char * string){
  
  if (string == NULL){
    fputs("seq argument cannot be NULL",stderr);
    exit(1);
  }

  int i;
  for(i=0;i<length;i++){
    string[i] = binary_nucleotide_to_char(nucleotides[i]);
  }
  
  string[length] = '\0';
  return string;
}



//The first argument - seq - is a C string in A,C,G,T format
//The second argument - quality - is a string of qualities for the sequence, one byte per base.
//quality cutoff argument defines the threshold for quality
//return total number of kmers read
//The third argument - length - is the length in bases of the sequence.
//Final argument - homopolymer_cutoff - allows you to break a sliding window at the end of a homopolymer.
//                 If homopolymer_cutoff= n > 0, then as soon as the latest base is the n-th in a homopolymeric sequence, the window is broken, and the next
//                  window only starts when there is  a new base.
//return total number of kmers read in - ie good kmers that go into windows
int get_sliding_windows_from_sequence(char * sequence,  char * qualities, int length, char quality_cut_off, short kmer_size, KmerSlidingWindowSet * windows, 
				      int max_windows, int max_kmers, boolean break_homopolymers, int homopolymer_cutoff){  

  char first_kmer[kmer_size+1];
  first_kmer[kmer_size]='\0';

  int i=0; //current index
  int count_kmers = 0;

  if (sequence == NULL){
    fputs("in get_sliding_windows_from_sequence, sequence is NULL\n",stderr);    
    exit(1);
  }

  if (length < kmer_size || max_windows == 0 || max_kmers == 0){
    return 0;
  }

  int index_windows = 0;
  
  //loop over the bases in the sequence
  //index i is the current position in input sequence -- it nevers decreases. 
  
  int hom_ct; //count how long current homopolymer run is

  do{

    //built first kmer, ie a stretch of kmer_size good qualities bases
    int j = 0; //count how many good bases
    hom_ct=0; 
    while ((i<length) && (j<kmer_size)){

      //collects the bases in the first kmer
      first_kmer[j] = sequence[i];

      //what is current homopolymer length
      if ( (j>0) && (first_kmer[j]==first_kmer[j-1]) )
	{
	  hom_ct++;
	}
      else if (j>=0)
	{
	  hom_ct=1;
	}


      if ((char_to_binary_nucleotide(sequence[i]) == Undefined) || 
	  (quality_cut_off>0 && qualities[i]<= quality_cut_off)){
	j=0; //restart the first kmer 
      }
      else if ( (break_homopolymers==true) && (hom_ct>=homopolymer_cutoff) )
	{
	  
	  //now we may be in the middle of a very long hompoppolymer run.So we want to increment i sufficiently to hit the next base
	  int first_base_after_homopolymer=i;
	  while ( (first_base_after_homopolymer<length) && (sequence[first_base_after_homopolymer]==first_kmer[j]) )
	    {
	      first_base_after_homopolymer++;
	    }

	  i=first_base_after_homopolymer-1; //we are going to add one at the end of the loop, just below
	  j=0; //restart the first kmer

	}
      
      else{
	j++;
      }

      i++; 
    }

    if (j==kmer_size){ //ie we did not parse the entire sequence looking for a single good kmer, the first kmer
      
      count_kmers++;

      //new sliding window
      if (index_windows>=max_windows){
	  fputs("number of windows is bigger than max_windows in get_sliding_windows_from_sequence",stderr);
	  exit(1);
	}

      KmerSlidingWindow * current_window =&(windows->window[index_windows]);

      int index_kmers = 0;
      //do first kmer
      BinaryKmer tmp_bin_kmer;
      seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer);
      binary_kmer_assignment_operator(current_window->kmer[index_kmers] , tmp_bin_kmer);

      //do the rest --
      index_kmers++;
    
      while(i<length){
	
	if (index_kmers>=max_kmers){
	  fputs("number of kmers is bigger than max_kmers in get_sliding_windows_from_sequence - second check\n",stderr);
	  exit(1);
	}

	Nucleotide current_base = char_to_binary_nucleotide(sequence[i]);

	if ( (i>0) && (sequence[i] == sequence[i-1]) )
	  {
	    hom_ct++;
	  }
	else 
	  {
	    hom_ct=1;
	  }

	if ((current_base == Undefined) ||
	    (quality_cut_off!=0 && qualities[i]<= quality_cut_off)){
	  i++; 
	  break;
	}
	else if ( (break_homopolymers==true) && (hom_ct>=homopolymer_cutoff) )
	  {
	    //now we may be in the middle of a very long homopolymer run.So we want to increment i sufficiently to go beyond         
	    int first_base_after_homopolymer=i;
	    while ( (first_base_after_homopolymer<length) && (sequence[first_base_after_homopolymer]==sequence[i]) )
	      {
		first_base_after_homopolymer++;
	      }

	    i=first_base_after_homopolymer; 
	    break;
	  }


	//set the kmer to previous
	binary_kmer_assignment_operator(current_window->kmer[index_kmers], current_window->kmer[index_kmers-1]);
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&(current_window->kmer[index_kmers]), current_base, kmer_size);

	index_kmers++;
	count_kmers++;
	i++;
      }

      current_window->nkmers = index_kmers; 
    
      index_windows++;
            
    }
  } while (i<length);
 
  windows->nwindows = index_windows;

  return count_kmers;
}






//caller passes in preallocated BinaryKmer, which is also returned in the return value
BinaryKmer* seq_to_binary_kmer(char * seq, short kmer_size, BinaryKmer* prealloced_kmer){

  //sanity checks
  if (seq==NULL)
    {
      printf("DO not passs null ptr to seq_to_binary_kmer. Exiting..\n");
      exit(1);
    }
  if (strlen(seq) != kmer_size)
    {
      printf("Calling seq_to_binary_kmer with  a sequence %s of length %d, but kmer size %d, which is different. Exiting", seq, (int) strlen(seq),  kmer_size);
      exit(1);
    }
  
  int j;
  binary_kmer_initialise_to_zero(prealloced_kmer);

  for(j=0;j<kmer_size;j++){

    if (char_to_binary_nucleotide(seq[j]) == Undefined){
      fputs("seq contains an undefined char\n",stderr);
      exit(1);
    }

    binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(prealloced_kmer, char_to_binary_nucleotide(seq[j]), kmer_size ); 
    
  }

  return prealloced_kmer;

}


  
  



//caller passes in allocated char*. This is returned and also set in 3rd argument.
//user of this method is responsible for deallocating the returned sequence
//note that the allocated space has to be kmer_size+1;
char * binary_kmer_to_seq(BinaryKmer* bkmer, short kmer_size, char * seq){
 
  BinaryKmer local_bkmer;
  binary_kmer_assignment_operator(local_bkmer, *bkmer);

  if (seq == NULL){
      fputs("seq argument cannot be NULL",stderr);
      exit(1);
    }

  int mask = 3; // 0000011 mask used to extract the two least significative bits
  int j;
  
  for(j=kmer_size-1; j>=0; j--){ //start from the back of the sequence
    
    //get translation for the two least significant bits
    seq[j] =  binary_nucleotide_to_char(local_bkmer[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] & mask);
    //shift right zam
    binary_kmer_right_shift_one_base(local_bkmer); //note this is a local copy internal to this function - not altering the original BinaryKmer

  }
  
  seq[kmer_size] = '\0';
  
  return seq;
}


//does not affect the kmer passed in as argument 1
BinaryKmer* binary_kmer_reverse_complement(BinaryKmer* kmer, short kmer_size, BinaryKmer* prealloc_reverse_kmer){
  binary_kmer_initialise_to_zero(prealloc_reverse_kmer);
  BinaryKmer local_copy_of_input_kmer;
  binary_kmer_assignment_operator(local_copy_of_input_kmer, *kmer);


  bitfield_of_64bits  mask = 3; //000..0011
  int j;

  //first complement the original kmer - xor with all 1's  
  for (j=0; j<NUMBER_OF_BITFIELDS_IN_BINARY_KMER; j++)
    {
      local_copy_of_input_kmer[j] ^= ~0;           
    }


  //then reverse
  for(j=0;j<kmer_size;j++){

    //make space for new base  
    binary_kmer_left_shift_one_base(*prealloc_reverse_kmer, kmer_size);

    //add base
    (*prealloc_reverse_kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] 
      = (*prealloc_reverse_kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] | (local_copy_of_input_kmer[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] & mask);

    binary_kmer_right_shift_one_base(local_copy_of_input_kmer);

  }

  return prealloc_reverse_kmer;
}


Nucleotide binary_kmer_get_last_nucleotide(BinaryKmer* kmer){
  
  bitfield_of_64bits bf = (*kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1] & 3; // mask against (11)base 2

  return (Nucleotide) bf;
  
}

Nucleotide binary_kmer_get_first_nucleotide(BinaryKmer* kmer,short kmer_size){

  int number_of_bitfields_fully_used = kmer_size/32;
  int number_of_bits_in_most_sig_bitfield = 2* (kmer_size-(32*number_of_bitfields_fully_used));

  bitfield_of_64bits bf = (*kmer)[NUMBER_OF_BITFIELDS_IN_BINARY_KMER - number_of_bitfields_fully_used -1] & (3 <<  (number_of_bits_in_most_sig_bitfield-2));
  
  bf >>= (number_of_bits_in_most_sig_bitfield-2);
  return bf;
  
}





void binary_kmer_alloc_kmers_set(KmerSlidingWindowSet * windows, int max_windows, int max_kmers){

  if (windows == NULL){
    fputs("Cannot pass a NULL window to alloc",stderr);
    exit(1);
  } 
  
  //allocate memory for the sliding windows         
  windows->max_nwindows = max_windows;

  windows->window = malloc(sizeof(KmerSlidingWindow) * max_windows);       
  if (windows->window== NULL){
    fputs("Out of memory trying to allocate an array of KmerSlidingWindow",stderr);
    exit(1);
  }
  windows->nwindows = 0;
  
  //allocate memory for every every sliding window
  int w;
  for(w=0;w<max_windows;w++){
    windows->window[w].nkmers=0;
    //windows->window[w].kmer = malloc(sizeof(BinaryKmer) * max_kmers);
    windows->window[w].kmer = (BinaryKmer*) malloc(sizeof(bitfield_of_64bits)*NUMBER_OF_BITFIELDS_IN_BINARY_KMER * max_kmers);
    if (windows->window[w].kmer == NULL){
      fputs("binary_kmer: Out of memory trying to allocate an array of BinaryKmer",stderr);
      exit(1);
    }      
  }      
  
}

void binary_kmer_free_kmers(KmerSlidingWindow * * kmers)
{

  free((*kmers)->kmer);
  free(*kmers); 
  *kmers = NULL;
}

void binary_kmer_free_kmers_set(KmerSlidingWindowSet * * kmers_set)
{
  int w;
  for(w=0;w<(*kmers_set)->max_nwindows;w++){
    free((*kmers_set)->window[w].kmer);
  }
  
  free((*kmers_set)->window);
  free(*kmers_set);
  *kmers_set = NULL;
}

void nucleotide_iterator(void (*f)(Nucleotide)){
  
  int i;
  for (i=0;i<4;i++){
    f(i);
  }
  
}

void nucleotide_iterator_orientation(void (*f)(Nucleotide n,Orientation o)){
  
  int i;
  for (i=0;i<4;i++){
    f(i,forward);
    f(i,reverse);
  }
  
}
