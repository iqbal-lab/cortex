/*
  binary_kmer.h - basic manipulation of binary kmers
  
*/

#ifndef EVENT_ENCODING_H_
#define EVENT_ENCODING_H_
#include <global.h>
#include <event_encoding.h>

typedef enum
 {
    Zero      = 0,
    One       = 1,
    Two       = 2,
    Three     = 3,
    Undefined = 4,
  } Nucleotide ;




Nucleotide char_to_binary_nucleotide(char c);
char binary_nucleotide_to_char(Nucleotide n);
Nucleotide reverse_binary_nucleotide(Nucleotide n);
boolean good_symbol(char c);


#endif /* EVENT_ENCODING_H_ */
