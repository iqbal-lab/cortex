// system libraries
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <libgen.h> // dirname
#include <errno.h>
#include <ctype.h> // tolower
#include <sys/types.h>
#include <dirent.h>

// third party libraries
#include <seq_file.h>
#include <string_buffer.h>


#include "binary_kmer.h"
#include "dB_graph.h"



//typedef signed char boolean; //aded for f_reader test only
//next two added for f_reader only
//#define MAX(x,y) ((x) >= (y) ? (x) : (y))
//#define MIN(x,y) ((x) <= (y) ? (x) : (y))

// cortex_var headers


#define is_base_char(x) ((x) == 'a' || (x) == 'A' || \
                         (x) == 'c' || (x) == 'C' || \
                         (x) == 'g' || (x) == 'G' || \
                         (x) == 't' || (x) == 'T')

// Note: parsing of homopolymers means input is different if reads are read
//       forward vs reverse-complemented!


//
// Sequence file loading
//

// cut-offs:
//  > quality_cutoff valid
//  < homopolymer_cutoff valid

void invalid_base_warning(char b)
{
  printf("WARNING");
  //warn("FYI (no cause for great concern) - Invalid sequence [%c] in read\n", b);
}

short _kmer_errors(char *kmer_str, char* qual_str,
                   short kmer_size, boolean is_fastq,
                   char quality_cutoff, int homopolymer_cutoff,
                   char *skip_homorun_base) // set to != 0 if need to skip a base
{
  short num_bp_to_skip = 0;
  int i;

  for(i = kmer_size-1; i >= 0; i--)
  {
    if(!is_base_char(kmer_str[i]))
    {
      if(kmer_str[i] != 'n' && kmer_str[i] != 'N')
      {
        // Invalid base
        invalid_base_warning(kmer_str[i]);
      }

      // skip
      num_bp_to_skip = i+1;
      break;
    }
  }

  if(is_fastq && quality_cutoff>0)
  {
    // Find latest poor-quality base position
    for(i = kmer_size-1; i >= 0; i--)
    {
      // Dev: we don't do any validation on the quality scores
      if(qual_str[i] <= quality_cutoff)
      {
        num_bp_to_skip = i+1;
        break;
      }
    }
  }

  if(homopolymer_cutoff != 0)
  {
    int run_length = 1;

    for(i = 1; i < kmer_size; i++)
    {
      if(kmer_str[i] == kmer_str[i-1])
        run_length++;
      else
        run_length = 1;

      if(run_length >= homopolymer_cutoff)
        num_bp_to_skip = MAX(i+1, num_bp_to_skip);
    }

    *skip_homorun_base = (run_length >= homopolymer_cutoff ? kmer_str[kmer_size-1]
                                                           : 0);
  }

  return num_bp_to_skip;
}

//pass in (0-based) position of previously read base (pos)
inline char _read_base(read_t *read_ob, int* pos, char *b, char *q, char read_qual)
{
  if (*pos+1 >= read_ob->seq.end)
  {
    return 0;
  }
  b=&(read_ob->seq.b[*pos+1]);
  if(!is_base_char(*b) && *b != 'n' && *b != 'N')
  {
    // Invalid base
    invalid_base_warning(*b);
  
    // Default to 'N'
    *b = 'N';
  }

  // Convert to upper case
  //  *b = toupper(*b);

  *pos = *pos+1;
  return 1;
}

//was _read_k_bases
inline char _next_k_bases(read_t *read_ob, int* pos, char *bases, char *quals, int k,
                          boolean is_fastq)
{
  printf("5\n");
  if(*pos+k>=read_ob->seq.end)
    {
      return 0;
    }
  
  if (is_fastq && (*pos+k>=read_ob->qual.end)   )
    {
      printf("WARNING\n");
    }

  memmove(bases, (read_ob->seq.b) + *pos+1, k);
  memmove(quals, (read_ob->qual.b) + *pos+1, k);
      // bases[i] = toupper(bases[i]);   // Convert to upper case

  return 1;
}

/*
inline char _read_all_bases(SeqFile *sf, StrBuf *bases, StrBuf *quals,
                            char read_qual)
{
  if(!seq_read_all_bases(sf, bases))
    return 0;

  if(read_qual && !seq_read_all_quals(sf, quals))
  {
    warn("%s:%d: Couldn't read quality scores [read: %s; path: %s; line: %lu]",
         __FILE__, __LINE__, seq_get_read_name(sf), seq_get_path(sf),
         seq_curr_line_number(sf));
  }

  // Convert to upper case
  strbuf_to_uppercase(bases);

  return 1;
}
*/



//was _read_first_kmer
char _get_first_kmer(read_t* read_ob, int* pos, 
		      char *kmer_str, char* qual_str,
                      short kmer_size, boolean is_fastq,
                      char quality_cutoff, int homopolymer_cutoff,
                      char first_base, char first_qual)
{
  short bp_to_skip;
  printf("6\n");
  if(first_base != 0)
  {
  printf("7\n");
    // Already got one base
    kmer_str[0] = first_base;
    qual_str[0] = first_qual;

    if(!_next_k_bases(read_ob, pos, kmer_str+1, qual_str+1, kmer_size-1, is_fastq))
      return 0;
  }
  else if(!_next_k_bases(read_ob, pos, kmer_str, qual_str, kmer_size, is_fastq))
  {
    return 0;
  }

  char skip_homorun_base = 0;

  while((bp_to_skip = _kmer_errors(kmer_str, qual_str, kmer_size, is_fastq,
                                   quality_cutoff, homopolymer_cutoff,
                                   &skip_homorun_base)) > 0)
  {
    //printf("proposed kmer: '%s' bases to skip: %i\n", kmer_str, bp_to_skip);
    
    while(bp_to_skip == kmer_size)
    {
      // Re read whole kmer
      if(!_next_k_bases(read_ob, pos, kmer_str, qual_str, kmer_size, is_fastq))
      {
        return 0;
      }

      bp_to_skip = 0;
      while(bp_to_skip < kmer_size && kmer_str[bp_to_skip] == skip_homorun_base)
      {
        bp_to_skip++;
      }
    }
    
    if(bp_to_skip > 0)
    {
      // Skip bp_to_skip bases
      int bp_to_keep = kmer_size - bp_to_skip;
      printf("Before move kmer_str is %s\n", kmer_str);
      memmove(kmer_str, kmer_str+bp_to_skip, bp_to_keep);
      memmove(qual_str, qual_str+bp_to_skip, bp_to_keep);

      if(!_next_k_bases(read_ob, pos, kmer_str+bp_to_keep, qual_str+bp_to_keep,
                        bp_to_skip, is_fastq))
      {
        return 0;
      }
    }
  }

  return 1;
}

inline void _process_read(seq_file_t* sf, read_t* read_obj, short kmer_size,
			  char* kmer_str, char* qual_str, 
                          char quality_cutoff, int homopolymer_cutoff,
                          int colour_index,
			  unsigned long long *bases_loaded)
{
  boolean is_fastq = seq_is_fastq(sf);

  char base, qual, prev_base;

  int homopol_length = 1;
  char keep_reading = 1;
  char is_first_kmer = 1;

  //read an entire read into buffer
  if (seq_read(sf, read_obj)>0)
  {
	memmove(kmer_str, (read_ob->seq.b), k);
	prev_base = kmer_str[kmer_size-1];

	// Set homopol_length for the first kmer in this contig
	if(homopolymer_cutoff != 0)
	  {
	    homopol_length = 1;
	    
	    int i;
	    for(i = kmer_size-2; i >= 0 && kmer_str[i] == prev_base; i--)
	      homopol_length++;
	  }
	
	int i=0;
	for(i=0; i<read_obj->seq.end; i++)
	  {
	    printf("2\n");
	    base=read_obj->seq.b[i];
	    printf("Base is %c\n", base);
	    if (is_fastq)
	      {
		qual = read_obj->qual.b[i];
	      }
	    
	    // Check for Ns and low quality scores
	    if(base == 'N' || (is_fastq && qual <= quality_cutoff))
	      {
		//pass in i by reference, skip beyond this and find first good kmer
		_get_first_kmer(read_obj, &i,  
				kmer_str, qual_str,
				kmer_size, is_fastq,
				quality_cutoff, homopolymer_cutoff,
				0, 0);
		printf("Now got first kmer %s\n", kmer_str);
		break;
	      }
	    
	    // Check homopolymer run length
	    if(homopolymer_cutoff != 0)
	      {
		if(prev_base == base)
		  {
		    homopol_length++;
		    
		    if(homopol_length >= homopolymer_cutoff)
		      {
			// Skip the rest of these bases
			while((keep_reading = _read_base(read_obj, &i, &base, &qual, is_fastq)) &&
			      base == prev_base);
			
			// Pass on first base that is != prev_base
			keep_reading = _get_first_kmer(read_obj, &i, 
						       kmer_str, qual_str,
						       kmer_size, is_fastq,
						       quality_cutoff, homopolymer_cutoff,
						       base, qual);
			printf("After homopol, first kmer is %s\n", kmer_str);
			break;
		      }
		  }
		else
		  {
		    // Reset homopolymer length
		    homopol_length = 1;
		  }
	      }
	
	    prev_base = base;
	    
	  }
	printf("Kmer is %s\n", kmer_str);
	// Store bases that made it into the graph
	(*bases_loaded) += contig_length;
	

  }

}

int main()
{
    // Open file
  seq_file_t *sf = seq_open("test.fa");
  if(sf == NULL)
  {
    printf("Cant open\n");
    exit(1);
  }
  short kmer_size=5;
    // Vars for reading in
  char kmer_str[kmer_size+1];
  kmer_str[0]='\0';
  char qual_str[kmer_size+1];
  qual_str[0]='\0';
  char quality_cutoff=0;
  int homopolymer_cutoff=0;
  int colour_index=0;
  unsigned long long bases_loaded=0;
  read_t* sr = seq_read_new();
  _process_read(sf, sr, kmer_size,
		kmer_str, qual_str,
		quality_cutoff, homopolymer_cutoff,
		colour_index,
		&bases_loaded);

  return 1;
}
