/* 
   binary_kmer.c - routines to 
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <binary_kmer.h>
#include <global.h>
#include <string.h>

//returns Undefined if given non AGCT character
Nucleotide char_to_binary_nucleotide(char c)
{
	switch (c)
	{
	case 'A':
	  return Adenine;
	case 'C':
	  return Cytosine;
	case 'G':
	  return Guanine;
	case 'T':
	  return Thymine;
	case 'a':
	  return Adenine;
	case 'c':
	  return Cytosine;
	case 'g':
	  return Guanine;
	case 't':
	  return Thymine;
	default:
	  return Undefined;
	}
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


char * seq_reverse_complement(char * in, int length, char * out){

  int k;
  for(k=0;k<length;k++){
    out[k] = reverse_char_nucleotide(in[length-k-1]);
  }
  out[length] = '\0';
  return out;
}



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
      printf("Non-existent nucleotide %i\n",n);
      exit(1);
    }
}

char binary_nucleotide_to_char(Nucleotide n)
{
	switch (n) {
	case Adenine:
	  return 'A';
	case Cytosine:
	  return 'C';
	case Guanine:
	  return 'G';
	case Thymine:
	  return 'T';
	default:
	  printf("Non existent binary nucleotide %d\n",n);
	  assert(0); 
	  return 'N'; //Don't really return this, must fail before this point. But stops compiler warning.
	}
}


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
//return total number of kmers read

int get_sliding_windows_from_sequence(char * seq,  char * qualities, int length, char quality_cut_off, short kmer_size, KmerSlidingWindowSet * windows, int max_windows, int max_kmers){  

  char first_kmer[kmer_size];
  int i=0; //current index
  int count_kmers = 0;

  BinaryKmer mask = (( (BinaryKmer) 1 << (2*kmer_size)) - 1); // mask binary 00..0011..11 as many 1's as kmer_size * 2 (every base takes 2 bits)

  if (seq == NULL){
    fputs("seq is NULL\n",stderr);    
    exit(1);
  }

  if (length < kmer_size || max_windows == 0 || max_kmers == 0){
    return 0;
  }


  int index_windows = 0;
  

  //loop over the bases in the sequence
  //index i is the current position in input sequence -- it nevers decreases. 
  
  do{

    //built first kmer, ie a stretch of kmer_size good qualities bases
    int j = 0; //count how many good bases
    while ((i<length) && (j<kmer_size)){

      //collects the bases in the first kmer
      first_kmer[j] = seq[i];

      if ((char_to_binary_nucleotide(seq[i]) == Undefined) || 
	  (quality_cut_off!=0 && qualities[i]<= quality_cut_off)){
	j=0; //restart the first kmer 
	
      }
      else{
	j++;
      }

      i++; 
    }    

    if (j==kmer_size){ //ie we did not parse the entire sequence looking for a single good kmer, the first kmer
      
      count_kmers++;
      first_kmer[kmer_size]='\0';

      //new sliding window
      if (index_windows>=max_windows){
	  fputs("number of windows is bigger than max_windows",stderr);
	  exit(1);
	}

      KmerSlidingWindow * current_window =&(windows->window[index_windows]);

      int index_kmers = 0;
      //do first kmer
      current_window->kmer[index_kmers]= seq_to_binary_kmer(first_kmer,kmer_size);

      //do the rest --
      index_kmers++;
    
      while(i<length){
	
	if (index_kmers>=max_kmers){
	  fputs("number of kmers is bigger than max_kmers\n",stderr);
	  exit(1);
	}

	Nucleotide current_base = char_to_binary_nucleotide(seq[i]);
	if ((current_base == Undefined) ||
	    (quality_cut_off!=0 && qualities[i]<= quality_cut_off)){
	  break;
	}
	//set the kmer to previous
	current_window->kmer[index_kmers]= current_window->kmer[index_kmers-1];
	//shift left - one base (ie 2 bits)
	current_window->kmer[index_kmers] <<= 2;
	//remove most significant base (using the mask x000.0011..11)
	current_window->kmer[index_kmers] &= mask;
	//add new base
	current_window->kmer[index_kmers] |= current_base;
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


//The first argument - seq - is a C string in A,C,G,T,N format. (Function handles bad characters)
//The second argument - length - is the length in bases of the sequence.
//we want a single sliding window, with all A's where the kmer would include an N, or Undefined nucleotide
//this is not ideal - but the caller is presumably going to break the kmer at N's, and here we force them to break at AAAAAA also.
int get_single_kmer_sliding_window_from_sequence(char * seq, int length, short kmer_size, KmerSlidingWindow* kmer_window)
{  

  if ( (kmer_window==NULL) || (seq==NULL))
    {
      printf("Do not pass NULL pointer to get_single_kmer_sliding_window_from_sequence\n");
      exit(1);
    }

  int number_of_steps_before_current_kmer_is_good=0; //good means free of bad characters. 
  int latest_base_we_have_read=0;
  int num_kmers=0;
  char first_kmer[kmer_size+1]; //as string
  first_kmer[kmer_size]='\0';
  Nucleotide current_base;

  long long current_good_kmer=~0;
  
  BinaryKmer mask = (( (BinaryKmer) 1 << (2*kmer_size)) - 1); // mask binary 00..0011..11 as many 1's as kmer_size * 2 (every base takes 2 bits)

  if (length < kmer_size )
    {
      return 0;
    }


  //set up first kmer
  int i;
  for (i=0; i<kmer_size; i++)
    {

      first_kmer[i]=seq[latest_base_we_have_read];
      current_base = char_to_binary_nucleotide(seq[latest_base_we_have_read]);

      //shift left - one base (ie 2 bits)
      current_good_kmer <<= 2;
      //remove most significant base (using the mask x000.0011..11)
      current_good_kmer &= mask;

      if (current_base==Undefined)
	{
	  //we will ignore contents of the string  first_kmer as it contains a bad character
	  number_of_steps_before_current_kmer_is_good=i+1;
	  current_good_kmer |= Adenine;
	}      
      else
	{
	  //add new base
	  current_good_kmer |= current_base;
	}

      latest_base_we_have_read++;
    }


  //add first kmer to window
  num_kmers++;

  if (number_of_steps_before_current_kmer_is_good==0)
    {
      kmer_window->kmer[num_kmers-1]=seq_to_binary_kmer(first_kmer,kmer_size);	  
    }
  else
    {
      kmer_window->kmer[num_kmers-1]=~0;
      number_of_steps_before_current_kmer_is_good--;
    }


  while (latest_base_we_have_read<length)
    {

      while ( (latest_base_we_have_read<length) && (number_of_steps_before_current_kmer_is_good>0))
	{
	  current_base = char_to_binary_nucleotide(seq[latest_base_we_have_read]);

	  //shift left - one base (ie 2 bits)
	  current_good_kmer <<= 2;
	  //remove most significant base (using the mask x000.0011..11)
	  current_good_kmer &= mask;
	  
	  if (current_base==Undefined)
	    {
	      current_good_kmer |= Adenine;
	      number_of_steps_before_current_kmer_is_good=kmer_size;
	    }
	  else
	    {

	      //add new base
	      current_good_kmer |= current_base;

	    }

	  num_kmers++;

	  //add a marked kmer to the window

	  kmer_window->kmer[num_kmers-1]=~0;
	  number_of_steps_before_current_kmer_is_good--;
	  latest_base_we_have_read++; 

	}

      //as long as previous kmer was good, you loop through this while loop
      while (latest_base_we_have_read<length)
	{
	  current_base = char_to_binary_nucleotide(seq[latest_base_we_have_read]);
	  //shift left - one base (ie 2 bits)
	  current_good_kmer <<= 2;
	  //remove most significant base (using the mask x000.0011..11)
	  current_good_kmer &= mask;

	  if (current_base==Undefined)
	    {
	      current_good_kmer |= Adenine;
	      number_of_steps_before_current_kmer_is_good=kmer_size;

	      //add a marked kmer to the window 
	      num_kmers++;
	      kmer_window->kmer[num_kmers-1]=~0; 
	      number_of_steps_before_current_kmer_is_good--; 
	      latest_base_we_have_read++;
	      break;
	    }
	  else
	    {
	  
	      //add new base
	      current_good_kmer |= current_base;
	      num_kmers++;
	      kmer_window->kmer[num_kmers-1]=current_good_kmer;      
	      latest_base_we_have_read++;	      
	    }

	}
  
    }
    
  kmer_window->nkmers=num_kmers;
  return num_kmers;

}



//Mario - worth noting this does direct translation, and does not return the smaller of kmer and rev_comp(kmer)
// answer: this is the correct behaviour, isn't it?

BinaryKmer seq_to_binary_kmer(char * seq, short kmer_size){
  
  int j;
  BinaryKmer kmer = 0;
  
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
  

  for(j=0;j<kmer_size;j++){
    //shift left,
    kmer <<= 2;
    //and then insert next nucleotide (in binary form) at the end.
    kmer |= char_to_binary_nucleotide(seq[j]);
    
    if (char_to_binary_nucleotide(seq[j]) == Undefined){
      fputs("seq contains an undefined char\n",stderr);
      exit(1);
    }
    
  }
  return kmer;
}



//caller passes in allocated char*. This is returned and also set in 3rd argument.
//user of this method is responsible for deallocating the returned sequence
//note that the allocated space has to be kmer_size+1;

char * binary_kmer_to_seq(BinaryKmer kmer, short kmer_size, char * seq){
 
  if (seq == NULL){
      fputs("seq argument cannot be NULL",stderr);
      exit(1);
    }

  int mask = 3; // 0000011 mask used to extract the two less significative bits
  int j;
  
  for(j=kmer_size-1; j>=0; j--){ //start from the back of the sequence
    
    //get translation for the two less significative bits
    seq[j] =  binary_nucleotide_to_char(kmer & mask);
    //shift right
    kmer >>=2;
  }
  
  seq[kmer_size] = '\0';
  
  return seq;
}


BinaryKmer binary_kmer_reverse_complement(BinaryKmer kmer, short kmer_size){

  BinaryKmer reverse_kmer = 0;
  BinaryKmer mask = 3; //000..0011
  int j;

  //xor with the mask first
  kmer ^= (( (BinaryKmer) 1 << (2*kmer_size)) - 1);// mask binary 00..0011..11 as many 1's as kmer_size * 2 (every base takes 2 bits)

  //reverse
  for(j=0;j<kmer_size;j++){
    //make space for new base
    reverse_kmer<<=2;
    //add base
    reverse_kmer |= (kmer & mask);
    kmer >>=2;
  }
  return reverse_kmer;
}


Nucleotide binary_kmer_get_last_nucleotide(BinaryKmer kmer){
  
  kmer &= 3; // mask against (11)base 2

  return kmer;
  
}

Nucleotide binary_kmer_get_first_nucleotide(BinaryKmer kmer,short kmer_size){
  
  kmer &= 3 << 2*(kmer_size-1); // mask against (11)base 2

  kmer >>= 2*(kmer_size-1);
  return kmer;
  
}


//this routines adds nucleotide add the end of kmer and removes the first base -- generating a new kmer that overlaps
//with kmer passed as argument

BinaryKmer binary_kmer_add_nucleotide_shift(BinaryKmer kmer,Nucleotide nucleotide, short kmer_size){
 
  kmer &= (((BinaryKmer) 1 << 2*(kmer_size-1))-1); // remove the last 2 bits - one base
  kmer <<= 2;
  kmer |= nucleotide; //add the new base

  return kmer;

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
    windows->window[w].kmer = malloc(sizeof(BinaryKmer) * max_kmers);
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
