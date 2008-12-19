/*
 * element.c
 */

#include <element.h>
#include <stdlib.h>

boolean element_is_key(Key key, Element e, short kmer_size){
  return key == e.kmer;
}


boolean element_smaller(Element  e1, Element e2){
	return e1.count < e2.count;
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
  e->count = 1;
}

void element_increment_count(Element * e){
  e->count++;
}


void element_print(FILE * file, Element * e,short kmer_size, char * prefix){

  char * string = binary_kmer_to_seq(e->kmer,kmer_size);

  if (prefix == NULL){
    fprintf(file,"%s %qd\n",string,e->count);
  }
  else{
    fprintf(file,"%s %s %qd\n",prefix,string,e->count);
  }
  
  free(string);
}

