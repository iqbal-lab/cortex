/*
  file_format.h
  
*/

#ifndef FILE_FORMAT_H_
#define FILE_FORMAT_H_

typedef enum
 {
   FASTA = 0,
   FASTQ = 1,
   CTX   = 2,
   UNSPECIFIED = 3,
   //VAR   = 3,
 } FileFormat ;

#endif //FILE_FORMAT_H_
