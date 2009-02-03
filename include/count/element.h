/*
  element.h defines the interface to count kmers in a set of reads
*/

#ifndef ELEMENT_H_
#define ELEMENT_H_
#include <binary_kmer.h>
#include <global.h>
#include <stdio.h>


typedef struct{
	BinaryKmer kmer;
	long long count;
} Element;

typedef BinaryKmer Key;

boolean element_smaller(Element,Element);

boolean element_is_key(Key,Element,short kmer_size);

Key element_get_key(BinaryKmer,short kmer_size);

void element_initialise(Element *,Key, short kmer_size);

void element_increment_count(Element*);

void element_print(FILE * file, Element* e,short kmer_size, char * prefix);
#endif /* ELEMENT_H_ */
