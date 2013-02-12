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
  error_correction.c 
*/

#ifndef ERROR_CORRECTION_H_
#define ERROR_CORRECTION_H_


// system libraries
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <errno.h>


 
// third party libraries
#include <seq_file.h>
#include <string_buffer.h>

// cortex_var headers
#include "binary_kmer.h"
#include "seq.h"
#include "dB_graph.h"


typedef enum {
  PrintUncorrected =0,
  PrintCorrected   =1,
  Discard          =2
} ReadCorrectionDecison;

typedef enum {
  DiscardReadIfLowQualBaseUnCorrectable  = 0,
  DontWorryAboutLowQualBaseUnCorrectable = 1,
} HandleLowQualUncorrectable;

typedef enum {
  Left = 0,
  Right= 1,
}WhichEndOfKmer;

inline void error_correct_file_against_graph(char* fastq_file, char quality_cutoff, char ascii_qual_offset,
					     dBGraph *db_graph, 
					     unsigned long long  *bases_modified,//total bases modified in dataset
					     int *bases_modified_count_array,//distribution across reads; how many of the read_length bases are fixed
					     int *posn_modified_count_array,//where in the read are we making corrections?
					     int bases_modified_count_array_size,
					     HandleLowQualUncorrectable policy);

ReadCorrectionDecison check_kmers_good(StrBuf* seq, StrBuf* qual, 
				       int num_kmers, int read_len,
				       int* kmers_in_graph, int* quals_good,
				       char quality_cutoff, int* first_good_kmer, 
				       dBGraph* db_graph);
  

boolean fix_end_if_unambiguous(WhichEndOfKmer which_end, StrBuf* read_buffer, int pos, 
			       StrBuf* kmer_buf, char* kmer_str,
			       dBGraph* dbg);

boolean mutate_base(StrBuf* strbuf, int which_base, int which_mutant);

#endif
