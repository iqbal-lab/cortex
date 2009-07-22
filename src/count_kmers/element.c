/*
 * element.c
 */

#include <element.h>
#include <stdlib.h>

boolean element_is_key(Key key, Element e, short kmer_size){

  char tmp_seq[kmer_size+1];

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

boolean db_node_check_status(Element * element, Status status){
  return element->status == status;
}

void element_initialise(Element * e, Key kmer, short kmer_size){

  char tmp_seq[kmer_size+1];

  e->kmer = element_get_key(kmer, kmer_size);
  e->count = 0;
  e->status = none;

}

void element_increment_count(Element * e){
  e->count++;
}


void element_print(FILE * file, Element * e,short kmer_size, char * prefix){

  char tmp_seq[kmer_size+1];

  binary_kmer_to_seq(e->kmer,kmer_size,tmp_seq);

  if (prefix == NULL){
    fprintf(file,"%s %qd\n",tmp_seq,e->count);
  }
  else{
    fprintf(file,"%s %s %qd\n",prefix,tmp_seq,e->count);
  }
  
}

