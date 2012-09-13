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
  routines to load files into dB Graph 
 */

// System libraries
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <libgen.h> // dirname
#include <errno.h>
#include <ctype.h> // tolower

// Third party libraries
#include <seq_file.h>

// Our headers
#include <binary_kmer.h>
#include <dB_graph.h>
#include <seq.h>
#include <file_reader.h>
#include <global.h>
#include <dB_graph_supernode.h>
#include <dB_graph_population.h>

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

#define is_base_char(x) ((x) == 'a' || (x) == 'A' || \
                         (x) == 'c' || (x) == 'C' || \
                         (x) == 'g' || (x) == 'G' || \
                         (x) == 't' || (x) == 'T' || \
                         (x) == 'n' || (x) == 'N')

// Uncomment this to print contigs+read names as they are loaded
// #define DEBUG_CONTIGS 1

// Note: parsing of homopolymers means input is different if reads are read
//       forward vs reverse-complemented!

int MAX_FILENAME_LENGTH=500;
int MAX_READ_LENGTH=10000;//should ONLY be used by test code

// Returns 1 on success, 0 on failure
// Sets errno to ENOTDIR if already exists but is not directory
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
char _ensure_dir_exists(const char *path, mode_t mode)
{
  struct stat st;

  if(stat(path, &st) != 0)
  {
    // Directory does not exist   
    return mkdir(path, mode) == 0 ? 1 : 0;
  }
  else if(!S_ISDIR(st.st_mode))
  {
    errno = ENOTDIR;
    return 0;
  }

  return 1;
}

// mkpath - ensure all directories in path exist
// Returns 1 on success, 0 on failure
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
char mkpath(const char *path, mode_t mode)
{
  char *copypath = strdup(path);

  size_t i = 0, j = 0;
  char status = 1;

  while(1)
  {
    while(path[i] == '.' || path[i] == '/')
      i++;

    j = i;

    while(path[j] != '.' && path[j] != '/' && path[j] != '\0')
      j++;

    if(i == j)
      break;

    char tmp = copypath[j];

    copypath[j] = '\0';

    if(!(status = _ensure_dir_exists(copypath, mode)))
      break;

    if(tmp == '\0')
      break;

    copypath[j] = tmp;
    i = j + 1;
  }

  free(copypath);

  return status;
}

//
// Isaac's re-write
//

// cut-offs:
//  > quality_cutoff valid
//  < homopolymer_cutoff valid

short _kmer_errors(SeqFile *sf, char *kmer_str, char* qual_str,
                   short kmer_size, char read_qual,
                   char quality_cutoff, int homopolymer_cutoff,
                   char *skip_homorun_base) // set to != 0 if need to skip a base
{
  short num_bp_to_skip = 0;
  int i;

  for(i = kmer_size-1; i >= 0; i--)
  {
    if(!is_base_char(kmer_str[i]))
    {
      // die -- invalid base
      fprintf(stderr, "%s:%d: Invalid base character '%i' "
                      "[path: %s; line: %lu]\n",
              __FILE__, __LINE__, (int)kmer_str[i],
              seq_get_path(sf), seq_curr_line_number(sf));
      exit(EXIT_FAILURE);
    }
    else if(kmer_str[i] == 'n' || kmer_str[i] == 'N')
    {
      // skip N
      num_bp_to_skip = i+1;
      break;
    }
  }

  if(read_qual)
  {
    // Find latest poor-quality base position
    for(i = kmer_size-1; i >= 0; i--)
    {
      // lt/le
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

      // lt/le
      if(run_length >= homopolymer_cutoff)
        num_bp_to_skip = MAX(i+1, num_bp_to_skip);
    }

    *skip_homorun_base = (run_length >= homopolymer_cutoff ? kmer_str[kmer_size-1]
                                                           : 0);
  }

  return num_bp_to_skip;
}

inline char _read_base(SeqFile *sf, char *b, char *q, char read_qual)
{
  if(!seq_read_base(sf, b))
  {
    return 0;
  }

  if(!is_base_char(*b))
  {
    // die -- invalid base
    fprintf(stderr, "%s:%d: Invalid base character %i [path: %s; line: %lu]\n",
            __FILE__, __LINE__, (int)*b,
            seq_get_path(sf), seq_curr_line_number(sf));
    exit(EXIT_FAILURE);
  }

  if(read_qual && !seq_read_qual(sf, q))
  {
    fprintf(stderr, "%s:%d: Couldn't read quality scores "
                    "[read: %s; path: %s; line: %lu]",
            __FILE__, __LINE__, seq_get_read_name(sf), seq_get_path(sf),
            seq_curr_line_number(sf));
  }

  // Convert to upper case
  *b = toupper(*b);

  return 1;
}

inline char _read_k_bases(SeqFile *sf, char *bases, char *quals, int k,
                          char read_qual)
{
  if(!seq_read_k_bases(sf, bases, k))
  {
    return 0;
  }

  if(read_qual && !seq_read_k_quals(sf, quals, k))
  {
    fprintf(stderr, "%s:%d: Couldn't read quality scores "
                    "[read: %s; path: %s; line: %lu]",
            __FILE__, __LINE__, seq_get_read_name(sf), seq_get_path(sf),
            seq_curr_line_number(sf));
  }

  // Convert to upper case
  int i;
  for(i = 0; i < k; i++)
    bases[i] = toupper(bases[i]);

  return 1;
}

inline char _read_all_bases(SeqFile *sf, StrBuf *bases, StrBuf *quals,
                            char read_qual)
{
  if(!seq_read_all_bases(sf, bases))
    return 0;

  if(read_qual && !seq_read_all_quals(sf, quals))
  {
    fprintf(stderr, "%s:%d: Couldn't read quality scores "
                    "[read: %s; path: %s; line: %lu]",
            __FILE__, __LINE__, seq_get_read_name(sf), seq_get_path(sf),
            seq_curr_line_number(sf));
  }

  // Convert to upper case
  strbuf_to_uppercase(bases);

  return 1;
}

void _print_kmer(BinaryKmer bkmer, short kmer_size)
{
  char str[kmer_size+1];
  binary_kmer_to_seq((BinaryKmer*)bkmer, kmer_size, str);
  str[kmer_size] = '\0';
  printf("%s", str);
}

char _read_first_kmer(SeqFile *sf, char *kmer_str, char* qual_str,
                      short kmer_size, char read_qual,
                      char quality_cutoff, int homopolymer_cutoff,
                      char first_base, char first_qual)
{
  short bp_to_skip;

  if(first_base != 0)
  {
    // Already got one base
    kmer_str[0] = first_base;
    qual_str[0] = first_qual;

    if(!_read_k_bases(sf, kmer_str+1, qual_str+1, kmer_size-1, read_qual))
      return 0;
  }
  else if(!_read_k_bases(sf, kmer_str, qual_str, kmer_size, read_qual))
  {
    return 0;
  }

  char skip_homorun_base = 0;

  while((bp_to_skip = _kmer_errors(sf, kmer_str, qual_str, kmer_size, read_qual,
                                   quality_cutoff, homopolymer_cutoff,
                                   &skip_homorun_base)) > 0)
  {
    //printf("proposed kmer: '%s' bases to skip: %i\n", kmer_str, bp_to_skip);
    
    while(bp_to_skip == kmer_size)
    {
      // Re read whole kmer
      if(!_read_k_bases(sf, kmer_str, qual_str, kmer_size, read_qual))
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

      memmove(kmer_str, kmer_str+bp_to_skip, bp_to_keep);
      memmove(qual_str, qual_str+bp_to_skip, bp_to_keep);

      if(!_read_k_bases(sf, kmer_str+bp_to_keep, qual_str+bp_to_keep,
                        bp_to_skip, read_qual))
      {
        return 0;
      }
    }
  }

  return 1;
}

inline void _process_read(SeqFile *sf, char* kmer_str, char* qual_str,
                          char quality_cutoff, int homopolymer_cutoff,
                          dBGraph *db_graph, int colour_index,
                          BinaryKmer curr_kmer, Element *curr_node,
                          Orientation curr_orient,
                          unsigned long long *bases_loaded,
                          unsigned long *readlen_count_array,
                          unsigned long readlen_count_array_size)
{
  // Hash table stuff
  Element *prev_node = NULL; // Element is a hash table entry
  Orientation prev_orient = forward;
  boolean curr_found;
  BinaryKmer tmp_key;

  short kmer_size = db_graph->kmer_size;
  char read_qual = seq_has_quality_scores(sf);

  char base, qual, prev_base;
  unsigned long contig_length;
  int homopol_length = 1;

  char keep_reading = 1;
  char is_first_kmer = 1;

  #ifdef DEBUG_CONTIGS
  printf(">%s\n", seq_get_read_name(sf));
  #endif

  while(keep_reading)
  {
    if(!is_first_kmer)
    {
      seq_to_binary_kmer(kmer_str, kmer_size, (BinaryKmer*)curr_kmer);

      element_get_key((BinaryKmer*)curr_kmer, kmer_size, &tmp_key);
      curr_node = hash_table_find_or_insert(&tmp_key, &curr_found, db_graph);
      curr_orient = db_node_get_orientation((BinaryKmer*)curr_kmer, curr_node,
                                             kmer_size);

      // Update coverage
      db_node_update_coverage(curr_node, 0, colour_index, 1);
      //db_node_update_coverage(curr_node, colour_index, 1); // EdgeArrayType removed
    }

    #ifdef DEBUG_CONTIGS
    _print_kmer(curr_kmer, kmer_size);
    #endif

    is_first_kmer = 0;

    contig_length = kmer_size;
    prev_base = kmer_str[kmer_size-1];

    // Set homopol_length for the first kmer in this contig
    if(homopolymer_cutoff != 0)
    {
      homopol_length = 1;

      int i;
      for(i = kmer_size-2; i >= 0 && kmer_str[i] == prev_base; i--)
        homopol_length++;
    }

    while((keep_reading = _read_base(sf, &base, &qual, read_qual)))
    {
      // Check for Ns and low quality scores
      // lt/le
      if(base == 'n' || base == 'N' || (read_qual && qual <= quality_cutoff))
      {
        keep_reading = _read_first_kmer(sf, kmer_str, qual_str,
                                        kmer_size, read_qual,
                                        quality_cutoff, homopolymer_cutoff,
                                        0, 0);
        break;
      }

      // Check homopolymer run length
      if(homopolymer_cutoff != 0)
      {
        if(prev_base == base)
        {
          homopol_length++;

          // lt/le
          if(homopol_length >= homopolymer_cutoff)
          {
            // Skip the rest of these bases
            while((keep_reading = _read_base(sf, &base, &qual, read_qual)) &&
                  base == prev_base);

            // Pass on first base that is != prev_base
            keep_reading = _read_first_kmer(sf, kmer_str, qual_str,
                                            kmer_size, read_qual,
                                            quality_cutoff, homopolymer_cutoff,
                                            base, qual);
            break;
          }
        }
        else
        {
          // Reset homopolymer length
          homopol_length = 1;
        }
      }

      // base accepted
      contig_length++;

      // Update curr -> prev hash table entry
      prev_node = curr_node;
      prev_orient = curr_orient;
      prev_base = base;

      // Construct new kmer
      //binary_kmer_assignment_operator(curr_kmer, prev_kmer);
      Nucleotide nuc = char_to_binary_nucleotide(base);
      binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end((BinaryKmer*)curr_kmer,
                                                                       nuc,
                                                                       kmer_size);

      #ifdef DEBUG_CONTIGS
      printf("%c", base);
      #endif

      // Lookup in db
      element_get_key((BinaryKmer*)curr_kmer, kmer_size, &tmp_key);
      curr_node = hash_table_find_or_insert(&tmp_key, &curr_found, db_graph);
      curr_orient = db_node_get_orientation((BinaryKmer*)curr_kmer, curr_node, kmer_size);

      // Update covg
      db_node_update_coverage(curr_node, 0, colour_index, 1);
      //db_node_update_coverage(curr_node, colour_index, 1); // EdgeArrayType removed

      // Add edge
      db_node_add_edge(prev_node, curr_node,
                       prev_orient, curr_orient,
                       kmer_size, 0, colour_index);
      // db_node_add_edge(prev_node, curr_node,
      //                  prev_orient, curr_orient,
      //                  kmer_size, colour_index); // EdgeArrayType removed
    }

    // Store bases that made it into the graph
    (*bases_loaded) += contig_length;

    #ifdef DEBUG_CONTIGS
    printf(" [%i]\n", contig_length);
    #endif

    // Store contig length
    if(readlen_count_array != NULL)
    {
      contig_length = MIN(contig_length, readlen_count_array_size-1);
      readlen_count_array[contig_length]++;
    }
  }
}

void load_se_seq_data_into_graph_colour(
  const char *file_path,
  char quality_cutoff, int homopolymer_cutoff, boolean remove_dups_se,
  char ascii_fq_offset, int colour_index, dBGraph *db_graph,
  unsigned long long *bad_reads, // number of reads that have no good kmers
  unsigned long long *dup_reads, // number of reads that pcr duplicates
  unsigned long long *bases_read, // total bases in file
  unsigned long long *bases_loaded, // bases that make it into the graph
  unsigned long *readlen_count_array, // histogram of contigs lengths
  unsigned long readlen_count_array_size) // contigs bigger go in final bin
{
  short kmer_size = db_graph->kmer_size;

  // DEV:
  // First check if this is a cortex binary
  // -> Load binary

  // Open file
  SeqFile *sf = seq_file_open(file_path);

  if(sf == NULL)
  {
    fprintf(stderr, "Couldn't open single-end sequence file '%s'\n", file_path);
    exit(EXIT_FAILURE);
  }

  seq_set_fastq_ascii_offset(sf, ascii_fq_offset);

  // Are we using quality scores
  char read_qual = seq_has_quality_scores(sf);

  // Vars for reading in
  char kmer_str[kmer_size+1];
  char qual_str[kmer_size+1];

  //
  // Hash table variables
  BinaryKmer curr_kmer;

  // Element and dBNode are the same thing
  Element *curr_node = NULL;
  Orientation curr_orient = forward;
  BinaryKmer tmp_key;
  boolean curr_found;

  while(seq_next_read(sf))
  {
    //printf("Started seq read: %s\n", seq_get_read_name(sf));

    if(_read_first_kmer(sf, kmer_str, qual_str, kmer_size, read_qual,
                        quality_cutoff, homopolymer_cutoff, 0, 0))
    {
      // Check if we want this read
      char is_dupe = 0;

      // Look up in table
      curr_found = false;

      seq_to_binary_kmer(kmer_str, kmer_size, &curr_kmer);

      element_get_key(&curr_kmer, kmer_size, &tmp_key);
      curr_node = hash_table_find_or_insert(&tmp_key, &curr_found, db_graph);
      curr_orient = db_node_get_orientation(&curr_kmer, curr_node, kmer_size);

      if(remove_dups_se == true)
      {
        // We are assuming that if you want to remove dups, you deal with one sample
        // at a time, in a single coloured graph.  I don't want to comtempate
        // multiple colours but only one type of read_start label. So I'm happy to
        // use hash_table_find here and not check if it is in this person's/colour
        // of the graph

        if(curr_found == true &&
           db_node_check_read_start(curr_node, curr_orient) == true)
        {
          (*dup_reads)++;
          is_dupe = 1;
        }
        else
        {
          db_node_set_read_start_status(curr_node, curr_orient);
        }
      }

      //debug
      //char tmpstr[32];
      //binary_kmer_to_seq(curr_kmer, kmer_size, tmpstr);
      //printf("First kmer: %s\n", tmpstr);

      if(!is_dupe)
      {
        // Update coverage
        db_node_update_coverage(curr_node, 0, colour_index, 1);
        //db_node_update_coverage(curr_node, colour_index, 1); // EdgeArrayType removed
      
        _process_read(sf, kmer_str, qual_str,
                      quality_cutoff, homopolymer_cutoff,
                      db_graph, colour_index,
                      curr_kmer, curr_node, curr_orient,
                      bases_loaded,
                      readlen_count_array, readlen_count_array_size);
      }
    }
    else
    {
      // Couldn't get a single kmer from read
      (*bad_reads)++;
    }
  }

  // Update with bases read in
  (*bases_read) += seq_total_bases_passed(sf) + seq_total_bases_skipped(sf);

  seq_file_close(sf);
}

void load_pe_seq_data_into_graph_colour(
  const char *file_path1, const char *file_path2,
  char quality_cutoff, int homopolymer_cutoff, boolean remove_dups_pe,
  char ascii_fq_offset, int colour_index, dBGraph *db_graph,
  unsigned long long *bad_reads, // number of reads that have no good kmers
  unsigned long long *dup_reads, // number of reads that are pcr duplicates
  unsigned long long *bases_read, // total bases in file
  unsigned long long *bases_loaded, // bases that make it into the graph
  unsigned long *readlen_count_array, // length of contigs loaded
  unsigned long readlen_count_array_size) // contigs bigger go in final bin
{
  short kmer_size = db_graph->kmer_size;

  // Open files
  SeqFile *sf1 = seq_file_open(file_path1);
  SeqFile *sf2 = seq_file_open(file_path2);

  if(sf1 == NULL)
  {
    fprintf(stderr, "Couldn't open paired-end sequence file '%s'\n", file_path1);
    exit(EXIT_FAILURE);
  }
  else if(sf2 == NULL)
  {
    fprintf(stderr, "Couldn't open paired-end sequence file '%s'\n", file_path2);
    exit(EXIT_FAILURE);
  }

  seq_set_fastq_ascii_offset(sf1, ascii_fq_offset);
  seq_set_fastq_ascii_offset(sf2, ascii_fq_offset);

  // Are we using quality scores
  char read_qual1 = seq_has_quality_scores(sf1);
  char read_qual2 = seq_has_quality_scores(sf2);

  // Vars for reading in
  char kmer_str1[kmer_size+1], kmer_str2[kmer_size+1];
  char qual_str1[kmer_size+1], qual_str2[kmer_size+1];

  //
  // Hash table variables
  BinaryKmer curr_kmer1, curr_kmer2;

  // Element and dBNode are the same thing
  Element *curr_node1 = NULL, *curr_node2 = NULL;
  Orientation curr_orient1 = forward, curr_orient2 = forward;
  BinaryKmer tmp_key;
  boolean curr_found1, curr_found2;

  while(1)
  {
    char read1 = seq_next_read(sf1);
    char read2 = seq_next_read(sf2);

    if(read1 != read2)
    {
      fprintf(stderr, "Paired-end files don't have the same number of reads\n");
      fprintf(stderr, "file1: %s\n", file_path1);
      fprintf(stderr, "file2: %s\n", file_path2);
      exit(EXIT_FAILURE);
    }
    else if(!read1)
    {
      // Neither have read
      break;
    }

    read1 = _read_first_kmer(sf1, kmer_str1, qual_str1, kmer_size, read_qual1,
                             quality_cutoff, homopolymer_cutoff, 0, 0);

    read2 = _read_first_kmer(sf2, kmer_str2, qual_str2, kmer_size, read_qual2,
                             quality_cutoff, homopolymer_cutoff, 0, 0);

    if(read1)
    {
      curr_found1 = false;

      // Get binary kmer
      seq_to_binary_kmer(kmer_str1, kmer_size, &curr_kmer1);

      // Look up first kmer
      element_get_key(&curr_kmer1, kmer_size, &tmp_key);
      curr_node1 = hash_table_find_or_insert(&tmp_key, &curr_found1, db_graph);
      curr_orient1 = db_node_get_orientation(&curr_kmer1, curr_node1, kmer_size);
    }
    else
    {
      // Couldn't get a single kmer from read
      (*bad_reads)++;
    }

    if(read2)
    {
      curr_found2 = false;

      // Get binary kmer
      seq_to_binary_kmer(kmer_str2, kmer_size, &curr_kmer2);

      // Look up second kmer
      element_get_key(&curr_kmer2, kmer_size, &tmp_key);
      curr_node2 = hash_table_find_or_insert(&tmp_key, &curr_found2, db_graph);
      curr_orient2 = db_node_get_orientation(&curr_kmer2, curr_node2, kmer_size);
    }
    else
    {
      // Couldn't get a single kmer from read
      (*bad_reads)++;
    }

    char is_dupe = 0;

    if(remove_dups_pe == true && (read1 || read2))
    {
      // We are assuming that if you want to remove dups, you deal with one sample
      // at a time, in a single coloured graph.  I don't want to comtempate
      // multiple colours but only one type of read_start label. So I'm happy to
      // use hash_table_find here and not check if it is in this person's/colour
      // of the graph

      // Check if all reads that were read in are marked as dupe => ie.
      //   A) Both reads read in and both are marked as read starts OR
      //   B) Just one read in and it is marked as read start

      if((!read1 || (curr_found1 == true &&
           db_node_check_read_start(curr_node1, curr_orient1) == true)) &&
         (!read2 || (curr_found2 == true &&
           db_node_check_read_start(curr_node2, curr_orient2) == true)))
      {
        // Both reads are dupes or neither are
        is_dupe = 1;
        (*dup_reads) += 2;
      }
      else
      {
        if(read1)
          db_node_set_read_start_status(curr_node1, curr_orient1);

        if(read2)
          db_node_set_read_start_status(curr_node2, curr_orient2);
      }
    }

    /*
    // For debugging
    char tmpstr1[kmer_size+1], tmpstr2[kmer_size+1];
    binary_kmer_to_seq((BinaryKmer*)curr_kmer1, kmer_size, tmpstr1);
    binary_kmer_to_seq((BinaryKmer*)curr_kmer2, kmer_size, tmpstr2);
    printf("Paired first kmers: %s %s\n", read1 ? tmpstr1 : "", read2 ? tmpstr2 : "");
    */

    if(!is_dupe)
    {
      if(read1)
      {
        // Update coverage
        db_node_update_coverage(curr_node1, 0, colour_index, 1);
        //db_node_update_coverage(curr_node1, colour_index, 1); // EdgeArrayType removed

        _process_read(sf1, kmer_str1, qual_str1, 
                      quality_cutoff, homopolymer_cutoff,
                      db_graph, colour_index,
                      curr_kmer1, curr_node1, curr_orient1,
                      bases_loaded,
                      readlen_count_array, readlen_count_array_size);
      }
      
      if(read2)
      {
        // Update coverage
        db_node_update_coverage(curr_node2, 0, colour_index, 1);
        //db_node_update_coverage(curr_node2, colour_index, 1); // EdgeArrayType removed

        _process_read(sf2, kmer_str2, qual_str2, 
                      quality_cutoff, homopolymer_cutoff,
                      db_graph, colour_index,
                      curr_kmer2, curr_node2, curr_orient2,
                      bases_loaded,
                      readlen_count_array, readlen_count_array_size);
      }
    }
  }

  // Update with bases read in
  (*bases_read) += seq_total_bases_passed(sf1) + seq_total_bases_skipped(sf1) +
                   seq_total_bases_passed(sf2) + seq_total_bases_skipped(sf2);

  seq_file_close(sf1);
  seq_file_close(sf2);
}

StrBuf* file_reader_get_strbuf_of_dir_path(char* path)
{
  char *tmp = strdup(path);
  StrBuf *dir = strbuf_create(dirname(tmp));
  strbuf_append_char(dir, '/');
  free(tmp);

  return dir;
}

// Go through all the files, loading data into the graph
// pass in bases_read to track amount of sequence read in, and
// bases_loaded to see how much passed filters and
// got into the graph
void load_se_filelist_into_graph_colour(
  char* se_filelist_path,
  int qual_thresh, int homopol_limit, boolean remove_dups_se,
  char ascii_fq_offset, int colour, dBGraph* db_graph, char is_colour_list,
  unsigned int *total_files_loaded,
  unsigned long long *total_bad_reads, unsigned long long *total_dup_reads,
  unsigned long long *total_bases_read, unsigned long long *total_bases_loaded,
  unsigned long *readlen_count_array, unsigned long readlen_count_array_size)
{
  qual_thresh += ascii_fq_offset;

  printf(is_colour_list ? "Load single-ended sequence colour list\n"
                        : "Load single-ended files\n");

  // Get absolute path
  char absolute_path[PATH_MAX+1];
  char* se_filelist_abs_path = realpath(se_filelist_path, absolute_path);

  if(se_filelist_abs_path == NULL)
  {
    fprintf(stderr, "Cannot get absolute path to filelist of SE files: %s\n",
            se_filelist_path);
    exit(EXIT_FAILURE);
  }

  printf("  path: %s\n", se_filelist_abs_path);

  FILE* se_list_file = fopen(se_filelist_abs_path, "r");

  if(se_list_file == NULL)
  {
    fprintf(stderr, "Cannot open filelist of SE files: %s\n",
            se_filelist_abs_path);
    exit(EXIT_FAILURE);
  }

  // Get directory path
  StrBuf *dir = file_reader_get_strbuf_of_dir_path(se_filelist_abs_path);

  // Stats
  unsigned int se_files_loaded = 0;

  unsigned long long se_bad_reads = 0;
  unsigned long long se_dup_reads = 0;

  unsigned long long se_bases_read = 0;
  unsigned long long se_bases_loaded = 0;

  StrBuf *line = strbuf_new();

  while(strbuf_reset_readline(line, se_list_file))
  {
    strbuf_chomp(line);

    if(strbuf_len(line) > 0)
    {
      // Get paths relative to filelist dir
      if(strbuf_get_char(line, 0) != '/')
        strbuf_insert(line, 0, dir, 0, strbuf_len(dir));

      // Get absolute paths
      char* path_ptr = realpath(line->buff, absolute_path);

      if(path_ptr == NULL)
      {
        fprintf(stderr, is_colour_list ? "Cannot find filelist of SE files: %s\n"
                                       : "Cannot find sequence file: %s\n",
                line->buff);
        exit(EXIT_FAILURE);
      }

      if(is_colour_list)
      {
        load_se_filelist_into_graph_colour(path_ptr,
          qual_thresh, homopol_limit, remove_dups_se,
          ascii_fq_offset, colour, db_graph, 0,
          &se_files_loaded,
          &se_bad_reads, &se_dup_reads,
          &se_bases_read, &se_bases_loaded,
          readlen_count_array, readlen_count_array_size);

        colour++;
      }
      else
      {
        load_se_seq_data_into_graph_colour(path_ptr,
          qual_thresh, homopol_limit, remove_dups_se,
          ascii_fq_offset, colour, db_graph,
          &se_bad_reads, &se_dup_reads,
          &se_bases_read, &se_bases_loaded,
          readlen_count_array, readlen_count_array_size);
      }

      se_files_loaded++;
    }
  }

  strbuf_free(line);
  strbuf_free(dir);

  // Finished reading single-ended file list
  fclose(se_list_file);

  // Update cumulative stats
  *total_files_loaded += se_files_loaded;
  *total_bad_reads += se_bad_reads;
  *total_dup_reads += se_dup_reads;
  *total_bases_read += se_bases_read;
  *total_bases_loaded += se_bases_loaded;

  // Print SE stats for this set of files
  printf("\nNum SE files loaded:%u\n", se_files_loaded);
  printf("\tKmers:%llu\n", hash_table_get_unique_kmers(db_graph));
  printf("\tNumber of bad reads:%llu\n", se_bad_reads);
  printf("\tNumber of dupe reads:%llu\n", se_dup_reads);
  printf("\tSE sequence parsed:%llu\n", se_bases_read);
  printf("\tTotal SE sequence that passed filters:%llu\n", se_bases_loaded);
}

// Go through all the files, loading data into the graph
void load_pe_filelists_into_graph_colour(
  char* pe_filelist_path1, char* pe_filelist_path2,
  int qual_thresh, int homopol_limit, boolean remove_dups_pe,
  char ascii_fq_offset, int colour, dBGraph* db_graph, char is_colour_lists,
  unsigned int *total_file_pairs_loaded,
  unsigned long long *total_bad_reads, unsigned long long *total_dup_reads,
  unsigned long long *total_bases_read, unsigned long long *total_bases_loaded,
  unsigned long *readlen_count_array, unsigned long readlen_count_array_size)
{
  qual_thresh += ascii_fq_offset;

  printf(is_colour_lists ? "Load paired-ended sequence colour list\n"
                         : "Load paired-ended files\n");

  // Get absolute paths
  char absolute_path1[PATH_MAX+1], absolute_path2[PATH_MAX+1];

  char* pe_filelist_abs_path1 = realpath(pe_filelist_path1, absolute_path1);
  char* pe_filelist_abs_path2 = realpath(pe_filelist_path2, absolute_path2);

  if(pe_filelist_abs_path1 == NULL)
  {
    fprintf(stderr, "Cannot get absolute path to filelist of PE files: %s\n",
            pe_filelist_path1);
    exit(EXIT_FAILURE);
  }
  else if(pe_filelist_abs_path2 == NULL)
  {
    fprintf(stderr, "Cannot get absolute path to filelist of PE files: %s\n",
            pe_filelist_path2);
    exit(EXIT_FAILURE);
  }

  // Open files
  FILE* pe_list_file1 = fopen(pe_filelist_abs_path1, "r");
  FILE* pe_list_file2 = fopen(pe_filelist_abs_path2, "r");

  if(pe_list_file1 == NULL)
  {
    fprintf(stderr, "Cannot open filelist of PE files: %s\n",
            pe_filelist_abs_path1);
    exit(EXIT_FAILURE);
  }
  else if(pe_list_file2 == NULL)
  {
    fprintf(stderr, "Cannot open filelist of PE files: %s\n",
            pe_filelist_abs_path2);
    exit(EXIT_FAILURE);
  }

  // Get directory paths for filelist files
  StrBuf *dir1 = file_reader_get_strbuf_of_dir_path(pe_filelist_abs_path1);
  StrBuf *dir2 = file_reader_get_strbuf_of_dir_path(pe_filelist_abs_path2);

  // Stats
  unsigned int pe_file_pairs_loaded = 0;

  unsigned long long pe_bad_reads = 0;
  unsigned long long pe_dup_reads = 0;

  unsigned long long pe_bases_read = 0;
  unsigned long long pe_bases_loaded = 0;

  StrBuf *line1 = strbuf_new();
  StrBuf *line2 = strbuf_new();

  while(1)
  {
    char read1 = strbuf_reset_readline(line1, pe_list_file1) > 0;
    char read2 = strbuf_reset_readline(line2, pe_list_file2) > 0;

    if(!read1 && !read2)
      break;

    strbuf_chomp(line1);
    strbuf_chomp(line2);

    if((strbuf_len(line1) == 0) != (strbuf_len(line2) == 0))
    {
      fprintf(stderr, "Paired-end files don't have the same number of lines:\n");
      fprintf(stderr, "File 1: %s", pe_filelist_path1);
      fprintf(stderr, "File 2: %s", pe_filelist_path2);
      exit(EXIT_FAILURE);
    }
    else if(strbuf_len(line1) > 0)
    {
      // Get paths relative to filelist dir
      if(strbuf_get_char(line1, 0) != '/')
        strbuf_insert(line1, 0, dir1, 0, strbuf_len(dir1));

      if(strbuf_get_char(line2, 0) != '/')
        strbuf_insert(line2, 0, dir2, 0, strbuf_len(dir2));

      // Get absolute paths
      char* path_ptr1 = realpath(line1->buff, absolute_path1);
      char* path_ptr2 = realpath(line2->buff, absolute_path2);

      if(path_ptr1 == NULL)
      {
        fprintf(stderr, is_colour_lists ? "Cannot find filelist of PE files: %s\n"
                                        : "Cannot find sequence file: %s\n",
                line1->buff);
        exit(EXIT_FAILURE);
      }
      else if(path_ptr2 == NULL)
      {
        fprintf(stderr, is_colour_lists ? "Cannot find filelist of PE files: %s\n"
                                        : "Cannot find sequence file: %s\n",
                line2->buff);
        exit(EXIT_FAILURE);
      }

      // Print file paths
      printf(is_colour_lists ? "Paired-end colourlists:\n"
                             : "Paired-end seq files:\n");

      printf("  File 1: %s\n", path_ptr1);
      printf("  File 2: %s\n", path_ptr2);

      // Read
      if(is_colour_lists)
      {
        load_pe_filelists_into_graph_colour(path_ptr1, path_ptr2,
          qual_thresh, homopol_limit, remove_dups_pe,
          ascii_fq_offset, colour, db_graph, 0,
          &pe_file_pairs_loaded, &pe_bad_reads, &pe_dup_reads,
          &pe_bases_read, &pe_bases_loaded,
          readlen_count_array, readlen_count_array_size);

        colour++;
      }
      else
      {
        load_pe_seq_data_into_graph_colour(path_ptr1, path_ptr2,
          qual_thresh, homopol_limit, remove_dups_pe,
          ascii_fq_offset, colour, db_graph,
          &pe_bad_reads, &pe_dup_reads,
          &pe_bases_read, &pe_bases_loaded,
          readlen_count_array, readlen_count_array_size);
      }

      pe_file_pairs_loaded++;
    }
  }

  /// Free line buffers
  strbuf_free(line1);
  strbuf_free(line2);

  // Free dir paths
  strbuf_free(dir1);
  strbuf_free(dir2);

  // Finished reading paired-ended file list
  fclose(pe_list_file1);
  fclose(pe_list_file2);

  // Update cumulative stats
  *total_file_pairs_loaded += pe_file_pairs_loaded;
  *total_bad_reads += pe_bad_reads;
  *total_dup_reads += pe_dup_reads;
  *total_bases_read += pe_bases_read;
  *total_bases_loaded += pe_bases_loaded;

  // Print SE stats for this set of files
  printf("\nNum PE file pairs loaded:%u\n", pe_file_pairs_loaded);
  printf("\tKmers:%llu\n", hash_table_get_unique_kmers(db_graph));
  printf("\tNumber of bad reads:%llu\n", pe_bad_reads);
  printf("\tNumber of dupe reads:%llu\n", pe_dup_reads);
  printf("\tPE sequence parsed:%llu\n", pe_bases_read);
  printf("\tTotal PE sequence that passed filters:%llu\n", pe_bases_loaded);
}


//
// End of loading sequence data
//


void initialise_binary_header_info(BinaryHeaderInfo* binfo, GraphInfo* ginfo)
{
  binfo->version=0;
  binfo->kmer_size=0;
  binfo->number_of_bitfields=0;
  binfo->number_of_colours=0;
  binfo->ginfo=ginfo;
}


//do not export the folloiwing internal function
// this is called repeatedly in a loop.  read_len_count_array is used to collect statistics on the length of reads,
// and needs to allow for the fact that we CUT reads at N's, low quality bases
void  load_kmers_from_sliding_window_into_graph_marking_read_starts_of_specific_person_or_pop(KmerSlidingWindowSet * windows, boolean* prev_full_ent, 
											      boolean* full_ent, long long* bases_loaded, boolean mark_read_starts, 
											      dBGraph* db_graph, EdgeArrayType type, int index, long long** read_len_count_array)
											      
{
  long long total_bases_loaded=0;

  Element * current_node  = NULL;
  Element * previous_node  = NULL;
  Orientation current_orientation=forward;
  Orientation previous_orientation=forward;
  BinaryKmer tmp_kmer;
  int i,j;

      for(i=0;i<windows->nwindows;i++){ //for each window
	KmerSlidingWindow * current_window = &(windows->window[i]);
	
	//update total bases loaded
	long long length_this_window = (long long) (current_window->nkmers+db_graph->kmer_size-1);
	total_bases_loaded+=length_this_window;
	
	if (read_len_count_array !=NULL)
	  {
	    //log that we have loaded a "read" of this length, provided it is not the last window - the last window length we will return, so  it can be cumulated, in case this is a long read
	    // this should never happen with fastq,as we have set the file_reader not to allow this
	    (*(read_len_count_array[length_this_window]))++;
	  }
	
	
	for(j=0;j<current_window->nkmers;j++){ //for each kmer in window
	  boolean found = false;
	  current_node = hash_table_find_or_insert(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),&found,db_graph);	  
	  if (current_node == NULL){
	    fputs("file_reader: problem - current kmer not found\n",stderr);
	    exit(1);
	  }
	  
	  if (! (i==0 && j==0 && *prev_full_ent == false && current_node == previous_node)){ //otherwise is the same old last entry
	    db_node_update_coverage(current_node,type, index, 1);
	  }
	  
	  current_orientation = db_node_get_orientation(&(current_window->kmer[j]),current_node, db_graph->kmer_size);
	  
	  
	  if (mark_read_starts==true)
	    {
	      if ( (i==0) && (j==0) && (current_orientation==forward))
		{
		  db_node_set_read_start_status(current_node, forward);
		}
	      else if ( (i==0) && (j==0) && (current_orientation==reverse))
		{
		  db_node_set_read_start_status(current_node, reverse);
		}
	    }
	  
	  
	  if (DEBUG){
	    char kmer_seq[db_graph->kmer_size+1];
	    char kmer_seq2[db_graph->kmer_size+1];
	    
	    printf("kmer i:%i j:%i:  %s %s %i\n",i,j,binary_kmer_to_seq(&(current_window->kmer[j]),db_graph->kmer_size,kmer_seq),
		   binary_kmer_to_seq(binary_kmer_reverse_complement(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph->kmer_size,kmer_seq2),
		   db_node_get_coverage(current_node, type, index));
	  }
	  
	  if (j>0){
	    
	    if (previous_node == NULL){
	      printf("PROBLEM i:%i j:%i bases loaded:%qd nkmers:%i prev_full_entry:%s\n",i,j,*(bases_loaded),current_window->nkmers,*prev_full_ent == true ? "true" : "false");
	      fputs("file_reader: problem - prev kmer not found\n",stderr);
	      exit(1);
	    }
	    else{ //this is the point at which the new element/node gets associated with the specific person
	      db_node_add_edge(previous_node,current_node,previous_orientation,current_orientation, db_graph->kmer_size, type, index);
	    }
	  }
	  previous_node = current_node;
	  previous_orientation = current_orientation;
	  
	}
      }
      
      *bases_loaded = *bases_loaded + total_bases_loaded;
    

}


//pass in a single kmer sliding window and the Sequence* it was derived from. Will find the nodes correspinding to this seqeunce
//and put them in array. Also will check that edges exist as expected from the Sequence*
void load_kmers_from_sliding_window_into_array(KmerSlidingWindow* kmer_window, Sequence* seq, dBGraph* db_graph, 
					       dBNode** array_nodes, Orientation* array_orientations, 
					       int max_array_size, 
					       boolean require_nodes_to_lie_in_given_colour, int colour)

{

      Element * current_node  = NULL;
      Element * previous_node  = NULL;
      Orientation current_orientation=forward;
      Orientation previous_orientation=forward;
      BinaryKmer tmp_kmer;
      int j;
	
      if (kmer_window->nkmers>max_array_size)
	{
	  printf("Cannot load_kmers_from_sliding_window_into_array as max_array_size %d < number f kmers in the window %d\n", max_array_size, kmer_window->nkmers);
	  exit(1);
	}

      for(j=0;j<kmer_window->nkmers;j++){ //for each kmer in window
	current_node = hash_table_find(element_get_key(&(kmer_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  
	if (current_node == NULL){
	  //printf("load_kmers_from_sliding_window_into_array: problem - current kmer not found\n");
	  //exit(1);
	}
	if ( (require_nodes_to_lie_in_given_colour==true) && 
	     (db_node_is_this_node_in_this_person_or_populations_graph(current_node, individual_edge_array, colour)==false) )
	  {
	    printf("This current node does not exist in colour %d\n", colour);
	    exit(1);
	  }

	  
	if (current_node !=NULL)
	  {
	    current_orientation = db_node_get_orientation(&(kmer_window->kmer[j]),current_node, db_graph->kmer_size);
	  }
	else
	  {
	    current_orientation = forward;
	  }
	//add to array
	array_nodes[j]        = current_node;
	array_orientations[j] = current_orientation;
	    
	previous_node = current_node;
	previous_orientation = current_orientation;
	  
      }
}





void set_binary_kmer_to_something_not_in_the_hash_table(BinaryKmer* bkmer, dBGraph* db_graph)
{
  BinaryKmer b;
  binary_kmer_initialise_to_zero(&b);
  boolean found=true;
  int count=0;
  int i=1;
  while (found==true)
    {
      i++;
      b[NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1]=(bitfield_of_64bits) i;

      dBNode* e = hash_table_find(element_get_key(&b, db_graph->kmer_size, &b), db_graph);
      if (e==NULL)
	{
	  found=false;
	}
      count++;
      if ((count>10000) && (found==true) )
	{
	  printf("Cortex needs, for non-obvious reasons, to find a kmer which is NOT in your graph, but after %d random tries, has failed to find one. Still looking.\nIf you are using a very small k, so the graph is saturating the kmer-space, this will go on forever....\n", count);
	}
    }
  
  binary_kmer_assignment_operator(*bkmer, b);
}

//The first argument - seq - is a C string in A,C,G,T,N format. (Function handles bad characters)
//The second argument - length - is the length in bases of the sequence.
//we want a single sliding window, using a kmer that does not exists in the hash where the kmer would include an N, or Undefined nucleotide
int get_single_kmer_sliding_window_from_sequence(char * seq, int length, short kmer_size, KmerSlidingWindow* kmer_window, dBGraph* db_graph)
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

  BinaryKmer marked_kmer; 
  set_binary_kmer_to_something_not_in_the_hash_table(&marked_kmer, db_graph);
  int i;
  /*

  for (i=0; i<NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
    {
      marked_kmer[i]=~0;
    }
  */
  BinaryKmer current_good_kmer;
  //binary_kmer_assignment_operator(current_good_kmer, marked_kmer); //initialisation
  binary_kmer_initialise_to_zero(&current_good_kmer); //zam DEBUG

  //long long current_good_kmer=~0;
  
  // don't think need this given new API --> BinaryKmer mask = (( (BinaryKmer) 1 << (2*kmer_size)) - 1); // mask binary 00..0011..11 as many 1's as kmer_size * 2 (every base takes 2 bits)

  if (length < kmer_size )
    {
      return 0;
    }


  //set up first kmer
  for (i=0; i<kmer_size; i++)
    {

      first_kmer[i]=seq[latest_base_we_have_read];
      current_base = char_to_binary_nucleotide(seq[latest_base_we_have_read]);

      if (current_base==Undefined)
	{
	  //we will ignore contents of the string  first_kmer as it contains a bad character
	  binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&current_good_kmer, Adenine, kmer_size );
	  number_of_steps_before_current_kmer_is_good=i+1;
	}      
      else
	{
	  binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&current_good_kmer, current_base, kmer_size );
	}

      latest_base_we_have_read++;
    }


  //add first kmer to window
  num_kmers++;

  BinaryKmer tmp_bin_kmer;
  binary_kmer_assignment_operator(tmp_bin_kmer, marked_kmer);
  //binary_kmer_initialise_to_zero(&tmp_bin_kmer);//zam debug

  if (number_of_steps_before_current_kmer_is_good==0)
    {
      seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer);
    }
  else
    {
      number_of_steps_before_current_kmer_is_good--;
    }
  binary_kmer_assignment_operator(kmer_window->kmer[num_kmers-1], tmp_bin_kmer);



  while (latest_base_we_have_read<length)
    {

      while ( (latest_base_we_have_read<length) && (number_of_steps_before_current_kmer_is_good>0))
	{
	  current_base = char_to_binary_nucleotide(seq[latest_base_we_have_read]);

	  if (current_base==Undefined)
	    {
	      binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&current_good_kmer, Adenine, kmer_size );
	      number_of_steps_before_current_kmer_is_good=kmer_size;
	    }
	  else
	    {
	      //add new base
	      binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&current_good_kmer, current_base, kmer_size );
	    }

	  num_kmers++;

	  //add a marked kmer to the window
	  binary_kmer_assignment_operator(kmer_window->kmer[num_kmers-1], marked_kmer);

	  number_of_steps_before_current_kmer_is_good--;
	  latest_base_we_have_read++; 

	}

      //as long as previous kmer was good, you loop through this while loop
      while (latest_base_we_have_read<length)
	{
	  current_base = char_to_binary_nucleotide(seq[latest_base_we_have_read]);

	  if (current_base==Undefined)
	    {
	      binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&current_good_kmer, Adenine, kmer_size );
	      number_of_steps_before_current_kmer_is_good=kmer_size;

	      //add a marked kmer to the window 
	      num_kmers++;
	      binary_kmer_assignment_operator(kmer_window->kmer[num_kmers-1], marked_kmer);

	      number_of_steps_before_current_kmer_is_good--; 
	      latest_base_we_have_read++;
	      break;
	    }
	  else
	    {
	  
	      //add new base
	      binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&current_good_kmer, current_base, kmer_size );
	      num_kmers++;
	      binary_kmer_assignment_operator(kmer_window->kmer[num_kmers-1],current_good_kmer);      
	      latest_base_we_have_read++;	      
	    }

	}
  
    }
    
  kmer_window->nkmers=num_kmers;
  return num_kmers;

}





// gets the next number_of_bases_to_load bases from fasta file, and returns them in the array of nodes.
// assumes this file has already been loaded into the graph.
// returns the number of nodes loaded. If this is less than what you asked for, you know it has hit the end of the file.
// We expect this to be used as follows:
// repeated calls of this function load etc into the LAST number_of_bases_to_load places of the relevant arrays

int load_seq_into_array(FILE* chrom_fptr, int number_of_nodes_to_load, int length_of_arrays, 
			dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
			Sequence* seq, KmerSlidingWindow* kmer_window, boolean expecting_new_fasta_entry, dBGraph * db_graph)
{
  
  int offset; //nodes are placed in the array from offset to offset+number_of_nodes_to_load-1, offset is calculated here to ensure they go at the end of the array
              //path_labels and path_string will contain only the edges - ie not include the first kmer.

  int offset_for_filereader=0;//this is the length of seq->seq that is preserved from last time

  if (expecting_new_fasta_entry==false)
    {
      offset = length_of_arrays-number_of_nodes_to_load-1;//the last node of the previous chunk will be found again, and put back in same place in array , hence the -1
      offset_for_filereader=db_graph->kmer_size;

      if (offset<0)
	{
	  printf("Error, offset<0. NOT expecting new fasta entry, length of arrays is %d, num of bases to load is %d, and offset is %d, their difference\n", length_of_arrays, number_of_nodes_to_load, offset);
	  exit(1);
	}
   
    }
  else //expecting a new entry
    {
      offset = length_of_arrays-number_of_nodes_to_load;

      if (offset<0)
        {
          printf("Error, offset<0. YES, expecting new fasta entry, length of arrays is %d, num of bases to load is %d, and offset is %d, defined by length - (bases to load _ kmer_size -1)\n", 
		 length_of_arrays, number_of_nodes_to_load, offset);
	  exit(1);
        }
    }

  long long seq_length=0;
  short kmer_size = db_graph->kmer_size;
  
  boolean full_entry=false;

  int chunk_length;
  int j;

  BinaryKmer marked_kmer; 
  set_binary_kmer_to_something_not_in_the_hash_table(&marked_kmer, db_graph);
  //  binary_kmer_set_all_bitfields(marked_kmer, ~( (bitfield_of_64bits) 0) );//will have all longlongs in array being ~0
  


  if (expecting_new_fasta_entry==false)
    {
      //3rd argument is limit set on number of bases in seq before read_sequence_from_fasta returns. We want number_of_nodes_to_load new bases, plus the kmer's woorth of bases already in seq
      chunk_length = read_sequence_from_fasta(chrom_fptr,seq,number_of_nodes_to_load+db_graph->kmer_size, expecting_new_fasta_entry, &full_entry, offset_for_filereader);
    }
  else
    {
      chunk_length = read_sequence_from_fasta(chrom_fptr,seq,number_of_nodes_to_load+db_graph->kmer_size-1, expecting_new_fasta_entry, &full_entry, offset_for_filereader);
    }




  //doesn't matter whether start of fasta entry or not. If YES, then we ignore the first kmer_size bases, as we are only interested in edges between nodes.
  // If NO, then seq has been preloaded with the last k bases from the previous time, which we want to ignore.
  for (j=0; j < number_of_nodes_to_load; j++)
    {
      path_string[offset+j]=seq->seq[db_graph->kmer_size+j];
    }
  path_string[offset+number_of_nodes_to_load]='\0';

    if (DEBUG){
    printf ("\n sequence returned from read_sequence_from_fasta to load_seq_into_array is %s - kmer size: %i - number of bases loaded inc the preassigned ones at start f seq  length: %i \n",seq->seq,db_graph->kmer_size, chunk_length);
    }
  
  j=0;
  seq_length += (long long) chunk_length;
  
  //number of nodes may be less than what we asked for, if we hit the end of the file
  int num_nodes = get_single_kmer_sliding_window_from_sequence(seq->seq,chunk_length,db_graph->kmer_size, kmer_window, db_graph);
  
  //sanity
  if (expecting_new_fasta_entry==true)
    {      
      if (num_nodes > number_of_nodes_to_load) 
	{
	  printf("Returned more nodes than asked for in load_seq_into_array. Exit/\n");
	  exit(1);
	}
    }
  else
    {
      //we expect number of nodes = num_of_nodes_to_load , except we had to preload seq with kmer_size bases, so that will give us an extra node at the start, which we need to get the first edge 
      if (num_nodes > number_of_nodes_to_load+1) 
	{
          printf("Returned morenodes than asked for inload_seq_into_array. Exit/\n");
          exit(1);

	}
    }
  

  Element * current_node  = NULL;
  Element * previous_node = NULL;
  
  Orientation current_orientation=forward;
  Orientation previous_orientation=forward;
  BinaryKmer tmp_kmer;

  //if num_kmers=0 (ie have hit end of file), then kmer_window->nkmers=0, so will skip this next "for" loop
  for(j=0;j<kmer_window->nkmers;j++)
    { //for each kmer in window

      if ( binary_kmer_comparison_operator(kmer_window->kmer[j], marked_kmer) ) //a non-kmer - ie anything that would have had an N in it
	                                                                        
	{
	  //corresponds to a kmer that contains an N
	  path_nodes[offset+j]        =NULL;
	  path_orientations[offset+j] =forward;

	  //iqbal
	  //if (j>1)
	  //  {  //edges to/from this node don't exist  remember path_nodes[m] is edge from node m to node m+1.
	  //    path_labels[offset+j-1] = Undefined;
	  //    path_labels[offset+j]   = Undefined;
	  //  }
	  previous_node=NULL;
	  continue;
	}
      
      //boolean found = false;
      current_node = hash_table_find(element_get_key(&(kmer_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  

      if (current_node == NULL){
	BinaryKmer tmp_dbg_kmer;
	char tmp_dbg_seq[db_graph->kmer_size+1];

	printf("Problem in load_seq_into_array - current kmer not found %s\n",
	       binary_kmer_to_seq(element_get_key(&(kmer_window->kmer[j]),db_graph->kmer_size, &tmp_dbg_kmer), db_graph->kmer_size,tmp_dbg_seq));
	exit(1);
	/* potential bug, and fix/more robust solution:
	previous_node=NULL; 
	continue;
	*/
      }

      current_orientation = db_node_get_orientation(&(kmer_window->kmer[j]),current_node, db_graph->kmer_size);

      path_nodes[offset+j]        = current_node;
      path_orientations[offset+j] = current_orientation;


      if (DEBUG){
	char kmer_seq[db_graph->kmer_size+1];
	printf("j=%d, Current node kmer is  %s\n",j, binary_kmer_to_seq(&(kmer_window->kmer[j]),db_graph->kmer_size,kmer_seq));
	if (current_orientation==forward)
	  printf("Current orientation is forward\n");
	else
	  printf("Current orientation is reverse");
      }
      
      if (j>0)
	{
	  if (previous_node == NULL)
	    {
	      if (j>0)
		{
		  path_labels[offset+j-1]=Undefined;
		}
	    }
	  else
	    {
	      BinaryKmer previous_k, current_k; 
	      char seq1[kmer_size+1];
	      char seq2[kmer_size+1];
	      
	      binary_kmer_assignment_operator(previous_k, previous_node->kmer);
	      binary_kmer_assignment_operator(current_k, current_node->kmer);
	      
	      if (previous_orientation == reverse){
		binary_kmer_assignment_operator(previous_k, *(binary_kmer_reverse_complement(&previous_k,kmer_size, &tmp_kmer)) );
	      }
	      
	      if (current_orientation == reverse){
		binary_kmer_assignment_operator(current_k, *(binary_kmer_reverse_complement(&current_k,kmer_size, &tmp_kmer)) ) ;
	      }
    
	      
	      if (DEBUG){
		printf("Found edge %s -%c-> %s\n",binary_kmer_to_seq(&previous_k,kmer_size,seq1),binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&current_k)),binary_kmer_to_seq(&current_k,kmer_size,seq2));
	      }
	      if (j>0)
		{
		  path_labels[offset+j-1]=binary_kmer_get_last_nucleotide(&current_k); //iqbal zam added if j>0 and added a -1
		}
	      
	    }
	}


	  
      previous_node = current_node;
      previous_orientation = current_orientation;
      
    }


  if ( (expecting_new_fasta_entry==true) || (num_nodes==0) )
    {
      return num_nodes;
    }
  else
    {
      return num_nodes-1;// remember the first node is just the same again as the last node of the last batch
    }

}


//returns number of kmers loaded*kmer_length
 //array_mean_readlens and array_total_seqs are arrays of length NUMBER_OF_COLOURS, so they can hold the mean read length+total seq in every colour
long long load_multicolour_binary_from_filename_into_graph(char* filename,  dBGraph* db_graph, GraphInfo* ginfo, int* num_cols_in_loaded_binary) 
//							    int** array_mean_readlens, long long** array_total_seqs)
{

  //printf("Load this binary - %s\n", filename);
  FILE* fp_bin = fopen(filename, "r");
  long long  seq_length = 0;
  dBNode node_from_file;
  element_initialise_kmer_covgs_edges_and_status_to_zero(&node_from_file);


  int count=0;

  if (fp_bin == NULL){
    printf("load_multicolour_binary_from_filename_into_graph cannot open file:%s\n",filename);
    exit(1); 
  }

  BinaryHeaderErrorCode ecode = EValid;
  BinaryHeaderInfo binfo;
  initialise_binary_header_info(&binfo, ginfo);//pass in the main graph info

  if (!(check_binary_signature_NEW(fp_bin, db_graph->kmer_size, &binfo, &ecode)))
    {
      printf("Cannot load this binary(%s) - signature check fails. Wrong max kmer, number of colours, or binary version. Exiting, error code %d\n", 
	     filename, ecode);
      exit(1);
    }
  else
    {
      *num_cols_in_loaded_binary = binfo.number_of_colours;
    }

  //always reads the multicol binary into successive colours starting from 0 - assumes the hash table is empty prior to this
  while (db_node_read_multicolour_binary(fp_bin,db_graph->kmer_size,&node_from_file, *num_cols_in_loaded_binary, binfo.version)){
    count++;
    
    dBNode * current_node  = NULL;
    BinaryKmer tmp_kmer;
    //current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size, &tmp_kmer),&found,db_graph);
    current_node = hash_table_insert(element_get_key(element_get_kmer(&node_from_file),db_graph->kmer_size, &tmp_kmer),db_graph);
    
    seq_length+=db_graph->kmer_size;
   
    int i;
    for (i=0; i<(*num_cols_in_loaded_binary) ; i++)
      {
	add_edges(current_node,individual_edge_array, i, get_edge_copy(node_from_file, individual_edge_array, i));
	db_node_update_coverage(current_node, individual_edge_array, i, db_node_get_coverage(&node_from_file, individual_edge_array,i));
      }

  }
  
  fclose(fp_bin);
  return seq_length;
}



// In some special cases (eg if you have pooled individuals and cleaned the graph, dumped a binary, and then reloaded  that in colour 0,
// and now want to load binaries into colour 1 only if the kmer is in the (cleaned) colour 0),
// then we set only_load_kmers_already_in_hash==true. In this case, we only load edges in that overlap with the cleaned graph, in colour_clean
long long load_single_colour_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph, 
								  int* mean_readlen, long long* total_seq,
								  boolean all_entries_are_unique, EdgeArrayType type, int index,
								  boolean only_load_kmers_already_in_hash, int colour_clean,
								  char* sample_identity,
								  boolean load_all_kmers_but_only_increment_covg_on_new_ones)
								  //last arg is to load the "union" of two graphs.)
{

  if ( (only_load_kmers_already_in_hash==true) && (colour_clean>=NUMBER_OF_COLOURS) )
    {
      printf("Called load_single_colour_binary_data_from_filename_into_graph and specified as clean-colour, colour %d, when this executable is compiled for %d colours only. Exit.\n", colour_clean, NUMBER_OF_COLOURS);
      exit(1);
    }


  printf("Open single colour binary: %s\n", filename);

  FILE* fp_bin = fopen(filename, "r");
  long long  seq_length = 0;
  dBNode tmp_node;
  element_initialise_kmer_covgs_edges_and_status_to_zero(&tmp_node);

  boolean found;
  int count=0;

  if (fp_bin == NULL)
    {
      //  exit(1); //TODO - prefer to print warning and skip file and return an error code?
      printf("Unable to open this binary %s\n", filename);
      exit(1);
    }


  //this will be inefficient once we have thousands of colours, but is only run once.
  GraphInfo* ginfo=graph_info_alloc_and_init();//no need to check return code
  BinaryHeaderErrorCode ecode=EValid;
  BinaryHeaderInfo binfo;
  initialise_binary_header_info(&binfo, ginfo);


  if (!(check_binary_signature_NEW(fp_bin, db_graph->kmer_size, &binfo, &ecode)))
    {
      printf("Cannot load this binary - fails signature check with error code %d. Exiting.\n", ecode);
      exit(1);
    }
  else
    {
      *mean_readlen = ginfo->mean_read_length[0];
      *total_seq    = ginfo->total_sequence[0];
      if  (strcmp(sample_identity, "undefined")==0)//we don't already have a name for this colour
	{
	  if (strcmp(ginfo->sample_ids[0], "undefined")!=0)//we've got  a name from this binary
	    {
	      sample_identity[0]='\0';
	      strcat(sample_identity, ginfo->sample_ids[0]);
	    }
	}
      else//we already have a name for this colour in our main GraphInfo
	{
	  if ( (strcmp(ginfo->sample_ids[0], "undefined")!=0) //we got a proper name from loading this binary
	       &&
	       (strcmp(ginfo->sample_ids[0], sample_identity)!=0) )
	    {
	      //then we have two different names for this colour. It must be a pool!
	      sample_identity[0]='\0';
	      strcat(sample_identity, "pool");
	    }
	}
    }

  if (binfo.number_of_colours!=1)
    {
      printf("Expecting a single colour binary, but instead this one has %d colours\n. Exiting.\n", binfo.number_of_colours);
      exit(1);
    }
  graph_info_free(ginfo);
  
  //Go through all the entries in the binary file
  // each time you load the info into a temporary node, and load them *** into colour number index ***
  while (db_node_read_single_colour_binary(fp_bin,db_graph->kmer_size,&tmp_node, type, index, binfo.version))
    {
      count++;
      found=false;//zam just added this
      dBNode * current_node  = NULL;
      BinaryKmer tmp_kmer;
      if (only_load_kmers_already_in_hash==false) //normal case
	{
	  if (!all_entries_are_unique)
	    {
	      current_node = hash_table_find_or_insert(element_get_key(element_get_kmer(&tmp_node),db_graph->kmer_size, &tmp_kmer),&found,db_graph);
	    }
	  else
	    {
	      current_node = hash_table_insert(element_get_key(element_get_kmer(&tmp_node),db_graph->kmer_size, &tmp_kmer), db_graph);
	    }
	  seq_length+=db_graph->kmer_size;
	  add_edges(current_node,individual_edge_array, index, get_edge_copy(tmp_node, individual_edge_array, index));
	  if ( (load_all_kmers_but_only_increment_covg_on_new_ones==false)//usual case
	       ||
	       ( (load_all_kmers_but_only_increment_covg_on_new_ones==true) && (found==false)) //loading union, and this is new
	       
	       )
	    {
	      db_node_update_coverage(current_node, individual_edge_array, index, db_node_get_coverage(&tmp_node, individual_edge_array,index) );
	    }
	}
      else
	{//check if node exists in hash already. If yes, then load edge and covg info into the appropriate colour
	  current_node = hash_table_find(element_get_key(element_get_kmer(&tmp_node),db_graph->kmer_size, &tmp_kmer), db_graph);
	  if (current_node !=NULL)
	    {
	      Edges pre_existing_edge = get_edge_copy(*current_node, individual_edge_array, colour_clean);
	      Edges edge_from_binary  = get_edge_copy(tmp_node, individual_edge_array, index);
	      Edges edge_to_load = pre_existing_edge & edge_from_binary; //only load edge from binary if is in the cleaned colour also.
	      add_edges(current_node,individual_edge_array, index,edge_to_load);
	      db_node_update_coverage(current_node, individual_edge_array, index, db_node_get_coverage(&tmp_node, individual_edge_array,index));
	    }
	}

    }

  fclose(fp_bin);
  return seq_length;

}



//ordinarily, only_load_kmers_already_in_hash==false, and colour_clean is ignored.
// If you have a clean graph in colour 0, and you only want to load nodes from the binaries that overlap with this,
// then set only_load_kmers_already_in_hash==true, and specify colour_clean to be that clean graph colour. Usually this is zero,
// we compile for 2 colours only, and we are loading into colour 1.
long long load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(char* filename, //_plus_maybe_samplename,  
											   dBGraph* db_graph, 
											   GraphInfo* db_graph_info, 
											   boolean all_entries_are_unique, 
											   EdgeArrayType type, 
											   int index,
											   boolean only_load_kmers_already_in_hash, 
											   int colour_clean,
											   boolean load_all_kmers_but_only_increment_covg_on_new_ones)
{
  // Get absolute path
  char absolute_path[PATH_MAX+1];
  char* filename_abs_path = realpath(filename, absolute_path);

  if(filename_abs_path == NULL)
  {
    fprintf(stderr, "Cannot get absolute path to ctxlist: %s\n", filename);
    exit(EXIT_FAILURE);
  }

  printf("loading ctxlist: %s\n", filename_abs_path);

  FILE* fptr = fopen(filename, "r");
  if (fptr == NULL)
  {
    printf("cannot open %s which is supposed to list all .ctx files for person "
           "with index %d\n", filename, index);
    exit(1); 
  }

  // Get directory path
  StrBuf *dir = file_reader_get_strbuf_of_dir_path(filename_abs_path);

  //file contains a list of .ctx filenames
  StrBuf *line = strbuf_new();

  int total_seq_loaded = 0;

  while(strbuf_reset_readline(line, fptr))
  {
    strbuf_chomp(line);

    if(strbuf_len(line) > 0)
    {
      // Get paths relative to filelist dir
      if(strbuf_get_char(line, 0) != '/')
        strbuf_insert(line, 0, dir, 0, strbuf_len(dir));

      // Get absolute paths
      char* path_ptr = realpath(line->buff, absolute_path);
      
      if(path_ptr == NULL)
      {
        fprintf(stderr, "Cannot find .ctx binary: %s\n", line->buff);
        exit(EXIT_FAILURE);
      }

      //printf("Load this binary: %s, into this colour : %d\n", line, index);
      int mean_read_len_in_this_binary = 0;
      long long total_seq_in_this_binary = 0;

      total_seq_loaded += 
        load_single_colour_binary_data_from_filename_into_graph(path_ptr, db_graph,
          &mean_read_len_in_this_binary,&total_seq_in_this_binary,
          all_entries_are_unique, individual_edge_array, index,
          only_load_kmers_already_in_hash, colour_clean,
          db_graph_info->sample_ids[index],
          load_all_kmers_but_only_increment_covg_on_new_ones);
      
      all_entries_are_unique = false;
      
      graph_info_update_mean_readlen_and_total_seq(db_graph_info, index,
        mean_read_len_in_this_binary, total_seq_in_this_binary);

      printf("Loaded next binary; total kmers in graph is now %qd\n",
             hash_table_get_unique_kmers(db_graph));
    }
  }

  strbuf_free(line);
  strbuf_free(dir);

  fclose(fptr);
  //  free(temp1);
  return total_seq_loaded;  
}






//takes a filename 
// this file contains a list of filenames, each of these represents an individual (and contains a list of binaries for that individual).
// these go into successive colours, starting with first_colour
long long load_population_as_binaries_from_graph(char* filename, int first_colour,boolean about_to_load_first_binary_into_empty_graph, 
						 dBGraph* db_graph, GraphInfo* db_graph_info, boolean only_load_kmers_already_in_hash, int colour_clean,
						 boolean load_all_kmers_but_only_increment_covg_on_new_ones)
{
  if(about_to_load_first_binary_into_empty_graph == true &&
     only_load_kmers_already_in_hash == true)
  {
    printf("You are trying to load binaries into an empty hash table, but are "
           "specifying that they should be compared with the existing hash "
           "table. User error\n");
    exit(EXIT_FAILURE);
  }

  // Get absolute path
  char absolute_path[PATH_MAX+1];
  char* filename_abs_path = realpath(filename, absolute_path);

  if(filename_abs_path == NULL)
  {
    fprintf(stderr, "Cannot get absolute path to colours: %s\n", filename_abs_path);
    exit(EXIT_FAILURE);
  }

  // Get directory path
  StrBuf *dir = file_reader_get_strbuf_of_dir_path(filename_abs_path);

  printf("Open this list of colours: %s\n", filename_abs_path);

  FILE* fp = fopen(filename, "r");
  if (fp == NULL)
  {
    //TODO - prefer to print warning and skip file and reutnr an error code?
    printf("load_population_as_binaries_from_graph cannot open file:%s\n", filename);
    exit(EXIT_FAILURE);
  }

  StrBuf *line = strbuf_new();

  int total_seq_loaded = 0;
  int which_colour = first_colour;

  while(strbuf_reset_readline(line, fp))
  {
    strbuf_chomp(line);

    if(strbuf_len(line) > 0)
    {
      if(which_colour >= NUMBER_OF_COLOURS)
      {
        printf("This filelist contains too many people, remember we have set a "
               "population limit of %d in variable NUMBER_OF_COLOURS. Cannot "
               "load into colour %d\n", NUMBER_OF_COLOURS, which_colour);
        exit(EXIT_FAILURE);
      }

      // Get paths relative to filelist dir
      if(strbuf_get_char(line, 0) != '/')
        strbuf_insert(line, 0, dir, 0, strbuf_len(dir));

      // Replace the first '\t' with '\0'
      strtok(line->buff, "\t");

      // Get absolute paths
      char* path_ptr = realpath(line->buff, absolute_path);

      if(path_ptr == NULL)
      {
        fprintf(stderr, "Cannot find ctxlist: %s\n",
                line->buff);
        exit(EXIT_FAILURE);
      }

      //printf("Open this filelist of binaries, %s,  all corresponding to the same colour:%d\n",
      //	     line, which_colour-1);

      total_seq_loaded = total_seq_loaded + 
        load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(
          path_ptr, db_graph,db_graph_info, 
          about_to_load_first_binary_into_empty_graph, 
          individual_edge_array, which_colour,
          only_load_kmers_already_in_hash, colour_clean,
          load_all_kmers_but_only_increment_covg_on_new_ones);

      about_to_load_first_binary_into_empty_graph = false;

      //printf("Loaded person %d, total kmers in graph %qd\n", which_colour,
      //       hash_table_get_unique_kmers(db_graph));

      which_colour++;
    }
  }

  strbuf_free(line);
  strbuf_free(dir);

  fclose(fp);

  //printf("Finished loading population, with total seq loaded %d\n",total_seq_loaded); 
  return total_seq_loaded;
}


// takes a filename (a colour_list)
// this file contains a list of filenames (one per colour), each of these represents an individual (and contains a list of single-colour binaries for that individual).
// Also takes two colour numbers. clean_colour is the clean colour, and in_colour is the colour ALL of these individuals are loaded into, thus:
// First take person 0's list of binaries and load them all into colour in_colour, BUT only load nodes that are already in the hash,
//  and only load those edges that are in the clean colour. Then dump a single-colour binary of colour in_colour,
// with filename = colour name PLUS a suffix added on the end.
void dump_successive_cleaned_binaries(char* filename, int in_colour,
  int clean_colour, char* suffix, dBGraph* db_graph, GraphInfo* db_graph_info)
{
  if(in_colour == clean_colour)
  {
    printf("In dump_successive_cleaned_binaries You cannot specify the same "
           "colour as both clean_colour and in_colour\n");
    exit(EXIT_FAILURE);
  }

  // Get absolute path
  char absolute_path[PATH_MAX+1];
  char* filename_abs_path = realpath(filename, absolute_path);

  if(filename_abs_path == NULL)
  {
    fprintf(stderr, "Cannot get absolute path to .colours file: %s\n", filename);
    exit(EXIT_FAILURE);
  }

  //printf("Open this list of colours: %s\n", filename);

  FILE* fp = fopen(filename, "r");
  if(fp == NULL)
  {
    printf("dump_successive_cleaned_binaries cannot open file:%s\n", filename);
    exit(EXIT_FAILURE);
  }

  int total_seq_loaded = 0;

  // Get directory path
  StrBuf *dir = file_reader_get_strbuf_of_dir_path(filename_abs_path);
  // Create buffer for reading in lines
  StrBuf *line = strbuf_new();

  while(strbuf_reset_readline(line, fp))
  {
    strbuf_chomp(line);

    if(strbuf_len(line) > 0)
    {
      // Get paths relative to filelist dir
      if(strbuf_get_char(line, 0) != '/')
        strbuf_insert(line, 0, dir, 0, strbuf_len(dir));

      // Get absolute paths
      char* path_ptr = realpath(line->buff, absolute_path);

      if(path_ptr == NULL)
      {
        fprintf(stderr, "Cannot find ctxlist: %s\n", line->buff);
        exit(EXIT_FAILURE);
      }

      total_seq_loaded +=
        load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(
          path_ptr, db_graph,db_graph_info, false, individual_edge_array,
          in_colour, true, clean_colour, false);
      
      char outfile[1000];
      outfile[0]='\0';
      sprintf(outfile, "%s_%s.ctx", line->buff, suffix);

      db_graph_dump_single_colour_binary_of_specified_colour(
        outfile, &db_node_check_status_not_pruned, db_graph, db_graph_info,
        in_colour, BINVERSION);

      //reset that colour:
      db_graph_wipe_colour(in_colour, db_graph);
      graph_info_set_seq(db_graph_info, in_colour, 0);
      graph_info_set_mean_readlen(db_graph_info, in_colour, 0);
    }
  }

  strbuf_free(line);
  strbuf_free(dir);

  fclose(fp);
}





void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
									    int max_read_length, dBGraph * db_graph)
{
    //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
  int seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  //we will be calling functions that apply to fastq as well as fasta. Since we are loading a ref fasta, we do not want to cut homopolymers,
  // and since we have fastq, the qualit-cutfoff os redundant
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  char quality_cut_off=0;
  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);

  BinaryKmer tmp_kmer;

  boolean full_entry = true;
  boolean prev_full_entry = true;
  
  int entry_length;

  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

    
    //last two arguments mean not to break at homopolymers
    int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,quality_cut_off,db_graph->kmer_size,windows,max_windows, 
						   max_kmers,break_homopolymers, homopolymer_cutoff);
    
    if (nkmers == 0) {
      //(*bad_reads)++;
    }
    else {
      Element * current_node  = NULL;

      
      for(i=0;i<windows->nwindows;i++){ //for each window
	KmerSlidingWindow * current_window = &(windows->window[i]);
	
	for(j=0;j<current_window->nkmers;j++){ //for each kmer in window

	  	  current_node = hash_table_find(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  	 

	  if (current_node){
	    db_node_set_status(current_node, exists_in_reference);
	  }
	  
	}
	
      }
    }
   
    if (full_entry == false){
      shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
    }

    prev_full_entry = full_entry;

  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);

}





// This is just like get_sliding_windows_from_sequence is seq.c, but this one breaks the window also when a kmer is not in the graph.
//The first argument - seq - is a C string in A,C,G,T format
//The second argument - quality - is a string of qualities for the sequence, one byte per base.
//quality cutoff argument defines the threshold for quality
//return total number of kmers read
//The third argument - length - is the length in bases of the sequence.
//return total number of kmers read
int get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										  KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph)  
{


  short kmer_size = db_graph->kmer_size;

  char first_kmer[kmer_size+1];
  first_kmer[kmer_size]='\0';

  BinaryKmer tmp_bin_kmer;
  BinaryKmer tmp_bin_kmer2;

      


  int i=0; //current index
  int count_kmers = 0;

  if (seq == NULL){
    fputs("in get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph, seq is NULL\n",stderr);    
    exit(1);
  }

  if (length < kmer_size || max_windows == 0 || max_kmers == 0){
    return 0;
  }

  int index_windows = 0;
  
  //loop over the bases in the sequence
  //index i is the current position in input sequence -- it nevers decreases. 
  
  do{

    //build first kmer, ie a stretch of kmer_size good qualities bases
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
      
      if (  hash_table_find(element_get_key(seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer), kmer_size, &tmp_bin_kmer2), db_graph) == NULL )
	{
	  j=0;
	  i=i-(kmer_size-1); // first kmer may be bad because of the first base. Start again from just after that
	  continue; //want to restart building a first kmer
	}

      count_kmers++;

      //new sliding window
      if (index_windows>=max_windows){
	  fputs("number of windows is bigger than max_windows",stderr);
	  exit(1);
	}

      KmerSlidingWindow * current_window =&(windows->window[index_windows]);

      int index_kmers = 0;
      //do first kmer
      seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer);
      binary_kmer_assignment_operator(current_window->kmer[index_kmers] , tmp_bin_kmer);

      //do the rest --
      index_kmers++;
    
      while(i<length){
	
	if (index_kmers>=max_kmers){
	  fputs("number of kmers is bigger than max_kmers\n",stderr);
	  exit(1);
	}

	Nucleotide current_base = char_to_binary_nucleotide(seq[i]);
	BinaryKmer tmp_next_kmer;

	//set to previous kmer, then shift an add new base
	binary_kmer_assignment_operator(tmp_next_kmer, current_window->kmer[index_kmers-1]);
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&tmp_next_kmer, current_base, kmer_size);

	if (  (current_base == Undefined) 
	      ||
	      (quality_cut_off!=0 && qualities[i]<= quality_cut_off)
	      )
	  {
	    break;
	  }
	else if (hash_table_find(element_get_key(&tmp_next_kmer, kmer_size, &tmp_bin_kmer2), db_graph) == NULL)
	  {
	    i++;
	    break;
	  } 
	binary_kmer_assignment_operator(current_window->kmer[index_kmers], tmp_next_kmer);
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


// This is just like get_sliding_windows_from_sequence is seq.c, but this one breaks a window if any of the kmers is not in the graph
// OR if any of the edges is not in the graph. ie entire window and edges must lie in graph

//The first argument - seq - is a C string in A,C,G,T format
//The second argument - quality - is a string of qualities for the sequence, one byte per base.
//quality cutoff argument defines the threshold for quality
//return total number of kmers read
//The third argument - length - is the length in bases of the sequence.
//return total number of kmers read
//note - this does not take arguments for homopolymer cutting - this in an interna function for aligning fastq to the graph
int get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										     KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph,
										     EdgeArrayType type, int index)  
{
  short kmer_size = db_graph->kmer_size;

  char first_kmer[kmer_size+1];
  first_kmer[kmer_size]='\0';

  BinaryKmer tmp_bin_kmer;
  BinaryKmer tmp_bin_kmer2;
  BinaryKmer tmp_bin_kmer3;
      


  int i=0; //current index
  int count_kmers = 0;

  if (seq == NULL){
    fputs("in get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph, seq is NULL\n",stderr);    
    exit(1);
  }

  if (length < kmer_size || max_windows == 0 || max_kmers == 0){
    return 0;
  }

  int index_windows = 0;
  
  //loop over the bases in the sequence
  //index i is the current position in input sequence -- it nevers decreases. 
  
  do{ //first do

    //build first kmer, ie a stretch of kmer_size good qualities bases
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
      
      if (  hash_table_find(element_get_key(seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer), kmer_size, &tmp_bin_kmer2), db_graph) == NULL )
	{
	  j=0;
	  i=i-(kmer_size-1); // first kmer may be bad because of the first base. Start again from just after that
	  continue; //want to restart building a first kmer       

	  // give up
	  //i=length;
	  //index_windows=0;
	  //break;
	  
	}

      count_kmers++;

      //new sliding window
      if (index_windows>=max_windows){
	  fputs("number of windows is bigger than max_windows",stderr);
	  exit(1);
	}

      KmerSlidingWindow * current_window =&(windows->window[index_windows]);

      int index_kmers = 0;
      //do first kmer
      seq_to_binary_kmer(first_kmer,kmer_size, &tmp_bin_kmer);
      binary_kmer_assignment_operator(current_window->kmer[index_kmers] , tmp_bin_kmer);

      //do the rest --
      index_kmers++;
    
      while(i<length){
	
	if (index_kmers>=max_kmers){
	  fputs("number of kmers is bigger than max_kmers\n",stderr);
	  exit(1);
	}

	Nucleotide current_base = char_to_binary_nucleotide(seq[i]);
	BinaryKmer tmp_curr_kmer;

	//set to previous kmer, then shift an add new base
	binary_kmer_assignment_operator(tmp_curr_kmer, current_window->kmer[index_kmers-1]);
	binary_kmer_left_shift_one_base_and_insert_new_base_at_right_end(&tmp_curr_kmer, current_base, kmer_size);

	dBNode* test_prev_node = hash_table_find(element_get_key(&(current_window->kmer[index_kmers-1]) , kmer_size, &tmp_bin_kmer2), db_graph);
	dBNode* test_curr_node = hash_table_find(element_get_key(&tmp_curr_kmer,                          kmer_size, &tmp_bin_kmer3), db_graph);


	if (  (current_base == Undefined) 
	      ||
	      (quality_cut_off!=0 && qualities[i]<= quality_cut_off)
	      )
	  {
	    //give up
	    //i=length;
	    //index_windows=0;
	    break;
	  }
	else if  (test_curr_node == NULL)
	  {
	    //give up
	    //i=length;
	    //index_windows=0;
	    i++;
	    break;
	  } 
	else //check edge exists between prev kmer and this one   
	  {
	    Orientation prev_orientation = db_node_get_orientation(&(current_window->kmer[index_kmers-1]),test_prev_node, kmer_size);
	    
	    if (! (db_node_edge_exist(test_prev_node, current_base, prev_orientation, type, index) ) )
	      {
		//give up
		//i=length;
		i=i-(kmer_size-1);
		//i++;
		//index_windows=0;
		break;
	      }

	  }
	binary_kmer_assignment_operator(current_window->kmer[index_kmers], tmp_curr_kmer);
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




void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
						  long long * bad_reads, int max_read_length, dBGraph * db_graph,
						  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index)
{
    //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
  int seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);

  char tmpseq[db_graph->kmer_size+1];
  tmpseq[db_graph->kmer_size]='\0';
  tmpseq[0]='\0';
  

  BinaryKmer tmp_kmer;
  boolean full_entry = true;
  boolean prev_full_entry = true;

  int entry_length;

  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

    
    //use quality cutoff of 0, arg 4 below
    // int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,0,db_graph->kmer_size,windows,max_windows, max_kmers);
    int nkmers = get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq,seq->qual,entry_length,0,
											       windows,max_windows, max_kmers, db_graph);
    
    if (nkmers == 0) 
      {
	(*bad_reads)++;
      }
    else 
      {
	Element * current_node  = NULL;
	//we will throw out any window that does not have ALL kmers in the graph
	for(i=0;i<windows->nwindows;i++)
	  { //for each window
	    KmerSlidingWindow * current_window = &(windows->window[i]);
	    
	    boolean all_kmers_in_this_window_are_in_graph=true;
	    

	    for(j=0;j<current_window->nkmers;j++)
	      { //for each kmer in window
		
		current_node = hash_table_find(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  	 
		if (current_node==NULL)
		  {
		    all_kmers_in_this_window_are_in_graph=false;
		  }
	      }

	    
	    if (all_kmers_in_this_window_are_in_graph==true)
	      {
		//print out this window as a "read". If this read has many windows, we will print each as a separate read (provided they lie in the graph)
		if (is_for_testing == false)
		  {
		    fprintf(fout, "> %s part %d \n", seq->name, i );
		    fprintf(fout, "%s", binary_kmer_to_seq(&(current_window->kmer[0]), db_graph->kmer_size, tmpseq) );
		    for(j=1;j<current_window->nkmers;j++){ 
		      fprintf(fout, "%c", binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&(current_window->kmer[j]))) );
		    }
		    fprintf(fout, "\n");
		  }
		else
		  {
		    for_test_array_of_clean_reads[*for_test_index][0]='\0';
		    strcat(for_test_array_of_clean_reads[*for_test_index], binary_kmer_to_seq(&(current_window->kmer[0]), db_graph->kmer_size, tmpseq) );
		    for(j=1;j<current_window->nkmers;j++){ 
		      char tmp_ch[2];
		      tmp_ch[0]=binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&(current_window->kmer[j])));
		      tmp_ch[1]='\0';
		      strcat(for_test_array_of_clean_reads[*for_test_index], tmp_ch );
		    }
		    *for_test_index=*for_test_index+1;
		  }
	      }
	  //else
	  //  {
	  //  }
	  
	  }
      }
    
    if (full_entry == false){
      shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
    }

    prev_full_entry = full_entry;

  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);

}


void read_fastq_and_print_subreads_that_lie_in_graph_breaking_at_edges_or_kmers_not_in_graph(FILE* fp, FILE* fout,
											     int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, 
														 boolean * full_entry), 
											     long long * bad_reads, int max_read_length, dBGraph * db_graph, 
											     EdgeArrayType type, int index,
											     boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index)
{
    //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
  int seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);

  char tmpseq[db_graph->kmer_size+1];
  tmpseq[db_graph->kmer_size]='\0';
  tmpseq[0]='\0';
  

  //BinaryKmer tmp_kmer;
  boolean full_entry = true;
  boolean prev_full_entry = true;

  int entry_length;

  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

    
    //use quality cutoff of 0, arg 4 below
    // int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,0,db_graph->kmer_size,windows,max_windows, max_kmers);
    //    int nkmers = get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(seq->seq,seq->qual,entry_length,0,
    //											       windows,max_windows, max_kmers, db_graph);
    
    int nkmers = get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(seq->seq, seq->qual, entry_length, 0,
												  windows, max_windows, max_kmers, db_graph, type, index);
    

    if (nkmers == 0) 
      {
	(*bad_reads)++;
      }
    else 
      {
	//Element * current_node  = NULL;
	
	for(i=0;i<windows->nwindows;i++)
	  { //for each window
	    KmerSlidingWindow * current_window = &(windows->window[i]);
	    
	    //print out this window as a "read". If this read has many windows, we will print each as a separate read (provided they lie in the graph)
	    if (is_for_testing == false)
	      {
		fprintf(fout, "> %s part %d \n", seq->name, i );
		fprintf(fout, "%s", binary_kmer_to_seq(&(current_window->kmer[0]), db_graph->kmer_size, tmpseq) );
		for(j=1;j<current_window->nkmers;j++){ 
		  fprintf(fout, "%c", binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&(current_window->kmer[j]))) );
		}
		fprintf(fout, "\n");
	      }
	    else
	      {
		for_test_array_of_clean_reads[*for_test_index][0]='\0';
		strcat(for_test_array_of_clean_reads[*for_test_index], binary_kmer_to_seq(&(current_window->kmer[0]), db_graph->kmer_size, tmpseq) );
		for(j=1;j<current_window->nkmers;j++){ 
		  char tmp_ch[2];
		  tmp_ch[0]=binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&(current_window->kmer[j])));
		  tmp_ch[1]='\0';
		  strcat(for_test_array_of_clean_reads[*for_test_index], tmp_ch );
		}
		*for_test_index=*for_test_index+1;
	      }
	    
		    
	  }
      }
	
	
    if (full_entry == false){
      shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
    }
    
    prev_full_entry = full_entry;
    
  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);
  
}



void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph)
{

  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry)
  {
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  int max_read_length = 3000;


  FILE* fptr = fopen(f_name, "r");
  if (fptr == NULL)
    {
      printf("cannot open chromosome fasta file:%s\n",f_name);
      exit(1);
    }

  read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(fptr, &file_reader, max_read_length, db_graph);


}


//returns the number of kmers loaded
int align_next_read_to_graph_and_return_node_array(FILE* fp, int max_read_length, dBNode** array_nodes, Orientation* array_orientations, 
						   boolean require_nodes_to_lie_in_given_colour,
						   boolean* full_entry,
						   int (* file_reader)(FILE * fp, Sequence * seq, 
								       int max_read_length,boolean new_entry, 
								       boolean * full_entry), 
						   Sequence* seq, KmerSlidingWindow* kmer_window,dBGraph * db_graph, int colour)
{
  

  //get next read as a C string and put it in seq. Entry_length is the length of the read.
  int entry_length = file_reader(fp,seq,max_read_length,*full_entry,full_entry);
  if (entry_length>0)
    {
      //turn it into a sliding window 
      int nkmers = get_single_kmer_sliding_window_from_sequence(seq->seq,entry_length, db_graph->kmer_size, kmer_window, db_graph);
      //work through the sliding window and put nodes into the array you pass in. Note this may find NULL nodes if the kmer is not in the graph
      load_kmers_from_sliding_window_into_array(kmer_window, seq, db_graph, array_nodes, array_orientations, 
						max_read_length-db_graph->kmer_size+1, require_nodes_to_lie_in_given_colour, colour);

      return nkmers;
    }
  else
    {
      return 0;
    }
}




// suppose 5p flank is : xxxxxxxxxxxxxxxxx and branch1 (as printed out) is yyyyyyyyyyy
// then the 5p flank in VariantFlanksAndBranches end at xxxxx,
// AND branch1 in VariantFlanksAndBranches starts    at xxxxx
// and branch2, however it starts, must end with the same kmer as branch1

// I want to read in the yyyyyy read, and get the kmers xxxxx, xxxxy, xxxyy, xxyyy, xyyyy, yyyyy
// so I want to be able to pass in xxxxx to this function and get all the nodes for the branch
//returns the number of kmers loaded
// MAKE SURE you pass in kmer_window capable of holding read-length + KMER  bases.
int given_prev_kmer_align_next_read_to_graph_and_return_node_array_including_overlap(char* prev_kmer, FILE* fp, int max_read_length, 
										     dBNode** array_nodes, Orientation* array_orientations, 
										     boolean require_nodes_to_lie_in_given_colour,
										     boolean* full_entry,
										     int (* file_reader)(FILE * fp, Sequence * seq, 
													 int max_read_length,boolean new_entry, 
													 boolean * full_entry), 
										     Sequence* seq, Sequence* seq_inc_prev_kmer, 
										     KmerSlidingWindow* kmer_window,dBGraph * db_graph, int colour)

{
  

  //get next read as a C string and put it in seq. Entry_length is the length of the read.
  int entry_length = file_reader(fp,seq,max_read_length,*full_entry,full_entry);

  seq_inc_prev_kmer->seq[0]='\0';
  seq_inc_prev_kmer->seq[entry_length+db_graph->kmer_size]='\0';
  // char prev_kmer_and_this_read[entry_length+db_graph->kmer_size+1];
  //prev_kmer_and_this_read[0]='\0';
  //prev_kmer_and_this_read[entry_length+db_graph->kmer_size]='\0';
  //strcpy(prev_kmer_and_this_read, prev_kmer);
  //strcat(prev_kmer_and_this_read, seq->seq);
  strcpy(seq_inc_prev_kmer->seq, prev_kmer);
  strcat(seq_inc_prev_kmer->seq, seq->seq);
  //turn it into a sliding window 
  int nkmers = get_single_kmer_sliding_window_from_sequence(seq_inc_prev_kmer->seq,entry_length+db_graph->kmer_size, db_graph->kmer_size, kmer_window, db_graph);
  //work through the sliding window and put nodes into the array you pass in. Note this may find NULL nodes if the kmer is not in the graph
  load_kmers_from_sliding_window_into_array(kmer_window, seq, db_graph, array_nodes, array_orientations, 
					    max_read_length+1, require_nodes_to_lie_in_given_colour, colour);

  return nkmers;
}



// returns 0 when hits the end of the file
// returns -1 when this variant is bad
// otherwise returns 1 (good variant)
// MAKE SURE your kmer_window is malloced to allow max_read_length PLUS KMER bases in a "read", as we want the transitions between
// flank and branches etc handled properly
// MAKE SURE seq_inc_prev_kmer also has space for an extra k bases at the start
int read_next_variant_from_full_flank_file(FILE* fptr, int max_read_length,
					   VariantBranchesAndFlanks* var, dBGraph* db_graph, 
					   int (file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
					   Sequence* seq, Sequence* seq_inc_prev_kmer, KmerSlidingWindow* kmer_window)
{
  int colour=-1; //ignored.
  
  boolean f_entry=true;
  var->len_flank5p = -1 + align_next_read_to_graph_and_return_node_array(fptr, max_read_length, var->flank5p, var->flank5p_or,  false, &f_entry, file_reader,
									 seq, kmer_window, db_graph, colour);
  if (!f_entry)
    {
      printf("One of these reads (5p flank) is longer than specified max read length. Last chunk we get was %s\n", seq->seq);
      exit(1);
    }

  if (var->len_flank5p==-1)
    {
      return 0; //end of file (should never have a zero length read as 5prime flank btw - should at least be one kmer
    }


  //use the read-id of the 5prime flank to get the variant name, and to double check that this IS a 5prime flank
  char* sub_ptr = strstr(seq->name, "_5p_flank");
  if (sub_ptr==NULL)
    {
      printf("Abort. Mandatory format for read-id of 5p flanks is >..some text.._5p_flank - but in this read, names %s, I cannot find the text \"_5p_flank\"\n", seq->name);
      exit(1);
    }
  else
    {
      size_t len = sub_ptr-seq->name;
      strncpy(var->var_name, seq->name, (int)len);
      var->var_name[(int)len]='\0';
      //printf("Found var name %s\n", var->var_name);
    }

  //save the sequence we have read:
  strncpy(var->seq5p, seq->seq, (int) strlen(seq->seq));
  var->seq5p[(int) strlen(seq->seq)]='\0';
  //so we have got the 5prime flank. Now we need to get all the kmers joining it to the branches
  char last_kmer_5p[db_graph->kmer_size+1];
  last_kmer_5p[0]='\0';
  last_kmer_5p[db_graph->kmer_size]='\0';
  strncpy(last_kmer_5p, seq->seq+ (int)strlen(seq->seq)-db_graph->kmer_size, db_graph->kmer_size);
  //printf("We think this %s is the last kmer in the 5p flank %s\n", last_kmer_5p, seq->seq);

  // we will also need the last kmer of either branch, to prepend in front of the 3p flank. 
  // At first sight, this seems complicated by the fact that sometimes one branch or even both branches are very short (<kmer)
  // however we have helpfully passed in the last kmer of the 50 flank, so we definitely have >k bases available to us


  char last_kmer_of_branch1[db_graph->kmer_size+1];
  last_kmer_of_branch1[0]='\0';
  last_kmer_of_branch1[db_graph->kmer_size]='\0';

  var->len_one_allele = -1 + 
    given_prev_kmer_align_next_read_to_graph_and_return_node_array_including_overlap(last_kmer_5p, fptr, max_read_length, 
										     var->one_allele, var->one_allele_or, 
										     false, &f_entry, file_reader,
										     seq, seq_inc_prev_kmer,kmer_window, db_graph, colour);
  if (!f_entry)
    {
      printf("One of these reads (branch1) is longer than specified max read length. The last chunk we got was %s\n", seq->seq);
      exit(1);
    }

  strncpy(last_kmer_of_branch1, seq_inc_prev_kmer->seq + (int)strlen(seq_inc_prev_kmer->seq)-db_graph->kmer_size, db_graph->kmer_size);

  //save the sequence we have read:
  strncpy(var->seq_one, seq->seq, (int) strlen(seq->seq));
  var->seq_one[(int) strlen(seq->seq)]='\0';
  
  var->len_other_allele = -1 + 
    given_prev_kmer_align_next_read_to_graph_and_return_node_array_including_overlap(last_kmer_5p, fptr, max_read_length, 
										     var->other_allele, var->other_allele_or, 
										     false, &f_entry, file_reader,
										     seq, seq_inc_prev_kmer, kmer_window, db_graph, colour);
  //printf("alt allele: %s, length %d\n", seq->seq, *len_branch_other );
  if (!f_entry)
    {
      printf("One of these reads (branch2) is longer than specified max read length Last chunk we got was %s\n\n", seq->seq);
      exit(1);
    }


  //save the sequence we have read:
  strncpy(var->seq_other, seq->seq, (int) strlen(seq->seq));
  var->seq_other[(int) strlen(seq->seq)]='\0';

  var->len_flank3p = -1 + 
    given_prev_kmer_align_next_read_to_graph_and_return_node_array_including_overlap(last_kmer_of_branch1, fptr, max_read_length, 
										     var->flank3p, var->flank3p_or, 
										     false, &f_entry, file_reader,
										     seq, seq_inc_prev_kmer, kmer_window, db_graph, colour);

  if (!f_entry)
    {
      printf("One of these reads (3p flank) is longer than specified max read length\n");
      exit(1);
    }

  //save the sequence we have read:
  strncpy(var->seq3p, seq->seq, (int) strlen(seq->seq));
  var->seq3p[(int) strlen(seq->seq)]='\0';
  
  return 1;

}
					   

/*

//array_mean_readlens is an array of length num_cols, giving the mean read length of data loaded into each colour
//array_total_seq is an array of length num_cols, giving the total amount of sequence in each colour (ie sequence loaded, AFTER filtering out by quality, PCR dups, homopolymers etc)
void print_binary_signature(FILE * fp,int kmer_size, int num_cols, int* array_mean_readlens, long long* array_total_seq){
  char magic_number[6];
  int version = BINVERSION;
  
  magic_number[0]='C';
  magic_number[1]='O';
  magic_number[2]='R';
  magic_number[3]='T';
  magic_number[4]='E';
  magic_number[5]='X';
  

  int num_bitfields = NUMBER_OF_BITFIELDS_IN_BINARY_KMER;

  fwrite(magic_number,sizeof(char),6,fp);
  fwrite(&version,sizeof(int),1,fp);
  fwrite(&kmer_size,sizeof(int),1,fp);
  fwrite(&num_bitfields, sizeof(int),1,fp);
  fwrite(&num_cols, sizeof(int), 1, fp);

  int i;
  for (i=0; i<num_cols; i++)
    {
      fwrite(&(array_mean_readlens[i]), sizeof(int), 1, fp);
    }
  for (i=0; i<num_cols; i++)
    {
      fwrite(&(array_total_seq[i]), sizeof(long long), 1, fp);
    }
  fwrite(magic_number,sizeof(char),6,fp);

}
*/


//we assume we are dumping N consecutive colours from the graph (the UI only supports
//1 colour (colour0) or all, but also I have a function for dumping one specific colour.
//let's say we want to dump from first_col to first_col + num_cols-1. These colours
//also correspond to those in the graph_info of course
//the final argument, version, is always BINVERSION in normal use, but in testing it might be an old version
// so I can test backward compatbility
void print_binary_signature_NEW(FILE * fp,int kmer_size, int num_cols, GraphInfo* ginfo, int first_col, int version)
{
  char magic_number[6];
  //  int version = BINVERSION;
  
  magic_number[0]='C';
  magic_number[1]='O';
  magic_number[2]='R';
  magic_number[3]='T';
  magic_number[4]='E';
  magic_number[5]='X';
  

  int num_bitfields = NUMBER_OF_BITFIELDS_IN_BINARY_KMER;

  fwrite(magic_number,sizeof(char),6,fp);
  fwrite(&version,sizeof(int),1,fp);
  fwrite(&kmer_size,sizeof(int),1,fp);
  fwrite(&num_bitfields, sizeof(int),1,fp);
  fwrite(&num_cols, sizeof(int), 1, fp);

  int i;
  for (i=first_col; i<first_col+num_cols; i++)
    {
      fwrite(&(ginfo->mean_read_length[i]), sizeof(int), 1, fp);
    }
  for (i=first_col; i<first_col+num_cols; i++)
    {
      fwrite(&(ginfo->total_sequence[i]), sizeof(long long), 1, fp);
    }

  if (version>5)
    {
      for (i=first_col; i<first_col+num_cols; i++)
	{
	  fwrite(&(ginfo->sample_id_lens[i]), sizeof(int), 1, fp);
	  fwrite(ginfo->sample_ids[i], sizeof(char), ginfo->sample_id_lens[i], fp);
	}
      for (i=first_col; i<first_col+num_cols; i++)
	{
	  fwrite(&(ginfo->seq_err[i]), sizeof(long double), 1, fp);
	}
      for (i=first_col; i<first_col+num_cols; i++)
	{
	  print_error_cleaning_object(fp, ginfo, i);
	}
    }

  fwrite(magic_number,sizeof(char),6,fp);

}

void print_error_cleaning_object(FILE* fp, GraphInfo* ginfo, int colour)
{
  fwrite(&(ginfo->cleaning[colour]->tip_clipping), sizeof(boolean), 1, fp);
  fwrite(&(ginfo->cleaning[colour]->remv_low_cov_sups), sizeof(boolean), 1, fp);
  fwrite(&(ginfo->cleaning[colour]->remv_low_cov_nodes), sizeof(boolean), 1, fp);
  fwrite(&(ginfo->cleaning[colour]->cleaned_against_another_graph), sizeof(boolean), 1, fp);
  fwrite(&(ginfo->cleaning[colour]->remv_low_cov_sups_thresh), sizeof(int), 1, fp);
  fwrite(&(ginfo->cleaning[colour]->remv_low_cov_nodes_thresh), sizeof(int), 1, fp);
  fwrite(&(ginfo->cleaning[colour]->len_name_of_graph_against_which_was_cleaned), sizeof(int), 1, fp);
  fwrite(ginfo->cleaning[colour]->name_of_graph_against_which_was_cleaned, sizeof(char), 
	 ginfo->cleaning[colour]->len_name_of_graph_against_which_was_cleaned, fp);
  
}

/*
//return yes if signature is consistent
boolean check_binary_signature(FILE * fp,int kmer_size, int bin_version, 
			       int* number_of_colours_in_binary, 
			       int** array_mean_readlens, 
			       long long** array_total_seqs,
			       int *return_binversion
			       )
{
  int read;
  char magic_number[6];
  boolean ret = false;

  read = fread(magic_number,sizeof(char),6,fp);
  if (read>0 &&
      magic_number[0]=='C' &&
      magic_number[1]=='O' &&
      magic_number[2]=='R' &&
      magic_number[3]=='T' &&
      magic_number[4]=='E' &&
      magic_number[5]=='X' ){

    int version;
    read = fread(&version,sizeof(int),1,fp);
    *return_binversion = version;

    //version 4 was the version I had for ages, up until I moved to using Uint32 in the binary
    if (read>0 && ( (version==BINVERSION)|| (version==4)) )
      {

	int kmer_size2;
	read = fread(&kmer_size2,sizeof(int),1,fp);
	if ((read>0) && (kmer_size2 == kmer_size) )
	  {

	    int num_bitfields;
	    read = fread(&num_bitfields,sizeof(int),1,fp);

	    if ( (read>0) && (num_bitfields==NUMBER_OF_BITFIELDS_IN_BINARY_KMER) )
	    
	      {
		int num_cols;
		read = fread(&num_cols,sizeof(int),1,fp);
		
		if ( (read>0) && (num_cols<=NUMBER_OF_COLOURS) && (num_cols>0)  )
		  { 
		    *number_of_colours_in_binary = num_cols;

		    int i;
		    boolean problem=false;
		    for (i=0; (i<num_cols) && (problem==false); i++)
		      {
			int mean_read_len;
			read = fread(&mean_read_len,sizeof(int),1,fp);
			if (read==0)
			  {
			    problem=true;
			  }
			else
			  {
			    *(array_mean_readlens[i]) = mean_read_len;
			  }
		      }
		    if (problem==false)
		      {
			for (i=0; (i<num_cols) && (problem==false); i++)
			  {
			    long long total_seq;
			    read = fread(&total_seq,sizeof(long long),1,fp);
			    if (read==0)
			      {
				problem=true;
			      }
			    else
			      {
				*(array_total_seqs[i]) = total_seq;
			      }
			  }
			
			magic_number[0]='\0';
			magic_number[1]='\0';
			magic_number[2]='\0';
			magic_number[3]='\0';
			magic_number[4]='\0';
			magic_number[5]='\0';
			read = fread(magic_number,sizeof(char),6,fp);
			if ( (problem==false)&&
			     (read>0)&&
			     magic_number[0]=='C' &&
			     magic_number[1]=='O' &&
			     magic_number[2]=='R' &&
			     magic_number[3]=='T' &&
			     magic_number[4]=='E' &&
			     magic_number[5]=='X' )
			  {
			    ret = true;
			  }
			else
			  {
			    printf("Binary is missing read-length/sequence info, or the end-of-header magic number.\n");
			  }
		      }
		    else
		      {
			printf("Problem reading chunk from binary. \n");
		      }
		  }
		else
		  {
		    printf("You are loading  binary with %d colours into a graph with %d colours - incompatible\n",
			   num_cols, NUMBER_OF_COLOURS);
		  }
	      }
	    else
	      {
		printf("Kmer of binary matches the current graph. However this binary was dumped with a different max_kmer size to that of the current graph\n");
		printf("This binary uses %d bitfields, and current graph uses %d\n", num_bitfields, NUMBER_OF_BITFIELDS_IN_BINARY_KMER);

	      }
	  }
	else
	  {
	printf("You are loading a binary with kmer=%d into a graph with kmer=%d - incompatible\n", kmer_size2, kmer_size);
	  }
      }
    else
      {

      }
    
  }
  else
    {
      printf("Binary fails signature check - was build for different kmer, or with different binary version. (For debug purposes: Magic number is %s and read is %d)\n", magic_number, read);
    }


  return ret;

}
*/




//return yes if signature is consistent
boolean check_binary_signature_NEW(FILE * fp,int kmer_size, 
				   BinaryHeaderInfo* binfo, BinaryHeaderErrorCode* ecode)
				   
{
  boolean bin_header_ok = query_binary_NEW(fp, binfo, ecode);

  if (bin_header_ok==true)
    {
      //just need to check the kmer is ok
      if (kmer_size==binfo->kmer_size)
	{
	  return true;
	}
      else
	{
	  return false;
	}
    }
  return false;
}


//return true if signature is readable, checks binversion, number of bitfields, magic number.
//does not check kmer is compatible with number of bitfields, leaves that to caller.
boolean query_binary_NEW(FILE * fp, BinaryHeaderInfo* binfo, BinaryHeaderErrorCode* ecode)
{
  int read;
  char magic_number[6];
  
  *ecode = EValid;
  read = fread(magic_number,sizeof(char),6,fp);
  if (read>0)
    {
      if (       magic_number[0]=='C' &&
		 magic_number[1]=='O' &&
		 magic_number[2]=='R' &&
		 magic_number[3]=='T' &&
		 magic_number[4]=='E' &&
		 magic_number[5]=='X' )
	{
	  
	  read = fread(&(binfo->version),sizeof(int),1,fp);
	  if (read>0)
	    {//can read version
	      if ((binfo->version >=4) && (binfo->version<=BINVERSION) )
		{//version is good
		  read = fread(&(binfo->kmer_size),sizeof(int),1,fp);
		  if (read>0)
		    {//can read bitfields
		      read = fread(&(binfo->number_of_bitfields),sizeof(int),1,fp);
		     
		      if (binfo->number_of_bitfields==NUMBER_OF_BITFIELDS_IN_BINARY_KMER)
			{//bitfuields are good

			  read = fread(&(binfo->number_of_colours),sizeof(int),1,fp);
			  
			  if ( read>0  )
			    {//can read colours

			      if (binfo->number_of_colours<=NUMBER_OF_COLOURS)
				{//colours are good
				  
				  //ok, the basic information looks OK
				  // get extra information. In all cases, return false and an error code if anything looks bad.
				  return get_extra_data_from_header(fp, binfo, ecode);
				  //also checks for the magic number at the end of the header
				}
			      else
				{//colours bad
				  *ecode = EBadColours;
				  return false;				      
				}
			    }
			  else
			    {//cant read colours
			      *ecode = ECannotReadNumColours;
			      return false;
			    }
			}//bitfields are good
		      else
			{//bitfields are bad
			  *ecode =  EWrongNumberBitfields;
			  return false;
			}
		    }//can read bitfields
		  else
		    {//cannot read bitfields
		      *ecode = ECannotReadNumBitfields;
		      return false;
		    }
		}//version is good
	      else
		{//version is bad
		  *ecode = EInvalidBinversion;
		  return false;
		}
	    }//can read version
	  else
	    {//cannot read version
	      *ecode =  ECannotReadBinversion;
	      return false;
	    }
	}//magic number good
      else
	{
	  *ecode=  ECanReadMagicNumberButIsWrong;
	  return false;
	}
    }
  else
    {
      *ecode = ECannotReadMagicNumber;
      return false;
    }

  
}


//assume num colours, number of bitfields, etc were fine. So now
//
// Binary Version 5: mean read lengths (NUMBER_OF_COLOURS of them) and then total sequences
// Binary Version 6: as 5, and then the Sample Id's for each colour, then the sequencing error rates
//                         and then the ErrorCleanings
boolean get_extra_data_from_header(FILE * fp, BinaryHeaderInfo* binfo, BinaryHeaderErrorCode* ecode)
{
  boolean no_problem=true;
  if ( (binfo->version==5) || (binfo->version==4) )//legacy
    {
      no_problem = get_read_lengths_and_total_seqs_from_header(fp, binfo, ecode);
    }
  else if (binfo->version==6)
    {
      no_problem = get_read_lengths_and_total_seqs_from_header(fp, binfo, ecode);
      if (no_problem==true)
	{
	  //get the extra stuff for binary version 6
	  no_problem = get_binversion6_extra_data(fp, binfo, ecode);
	}
    }
  else
    {
      *ecode = EInvalidBinversion;
      no_problem=false;
    }

  if (no_problem==true)
    {
      //only thing remaining to check is the end of header magic number
      int read;
      char magic_number[6];
      magic_number[0]='\0';
      magic_number[1]='\0';
      magic_number[2]='\0';
      magic_number[3]='\0';
      magic_number[4]='\0';
      magic_number[5]='\0';
      read = fread(magic_number,sizeof(char),6,fp);
      if (read==0)
	{
	  *ecode = ECannotReadEndOfHeaderMagicNumber;
	  no_problem=false;
	}
      else if (
	   magic_number[0]=='C' &&
	   magic_number[1]=='O' &&
	   magic_number[2]=='R' &&
	   magic_number[3]=='T' &&
	   magic_number[4]=='E' &&
	   magic_number[5]=='X' )
	{
	  //all good.
	}
      else
	{
	  no_problem=false;
	  *ecode = ECanReadEndOfHeaderMagicNumberButIsWrong;
	}
    }

  return no_problem;

}


boolean get_read_lengths_and_total_seqs_from_header(FILE * fp, BinaryHeaderInfo* binfo, BinaryHeaderErrorCode* ecode)
{
  int read;
  int i;
  boolean no_problem=true;
  
  for (i=0; (i<binfo->number_of_colours) && (no_problem==true); i++)
    {
      int mean_read_len=0;
      read = fread(&mean_read_len, sizeof(int),1,fp);
      if (read==0)
	{
	  no_problem=false;
	  *ecode= EFailedToReadReadLensAndCovgs;
	}
      else
	{
	  binfo->ginfo->mean_read_length[i] = mean_read_len;
	}
    }
  if (no_problem==true)
    {
      for (i=0; (i<binfo->number_of_colours) && (no_problem==true); i++)
	{
	  long long tot=0;
	  read = fread(&tot, sizeof(long long),1,fp);
	  if (read==0)
	    {
	      no_problem=false;
	      *ecode= EFailedToReadReadLensAndCovgs;
	    }
	  else
	    {
	      binfo->ginfo->total_sequence[i] =tot;
	    }
	}
    }

  return no_problem;
}

//Binary header version 6 includes sample id's, and Seq Error rate and Error Cleaning Info.
boolean  get_binversion6_extra_data(FILE * fp, BinaryHeaderInfo* binfo, BinaryHeaderErrorCode* ecode)
{
  int read;
  int i;
  boolean no_problem=true;

  //first get sample information
  for (i=0; (i<binfo->number_of_colours) && (no_problem==true); i++)
    {
      //get the length of the sample id name
      read = fread(&(binfo->ginfo->sample_id_lens[i]),sizeof(int),1,fp);
      if (read==0)
	{
	  no_problem=false;
	  *ecode = EFailedToReadSampleIds;
	}
      else
	{
	  //now get the actual sample id for colour i
	  if (binfo->ginfo->sample_id_lens[i]<MAX_LEN_SAMPLE_NAME)
	    {
	      read = fread(binfo->ginfo->sample_ids[i],sizeof(char),binfo->ginfo->sample_id_lens[i],fp);
	      if (read==0)
		{
		  no_problem=false;
		  *ecode = EFailedToReadSampleIds;
		}
	    }
	  else
	    {
	      no_problem=false;
	      *ecode = EFailedToReadSampleIdsSeemsTooLong;

	    }
	}
    }

  //now get the sequencing error rate
  for (i=0; (i<binfo->number_of_colours) && (no_problem==true); i++)
    {
      read = fread(&(binfo->ginfo->seq_err[i]),sizeof(long double),1,fp);
      if (read==0)
	{
	  no_problem=false;
	  printf("i is %d and num colours is %d, but read is zero  - problem reading binary header. Contact Zam\n", i, binfo->number_of_colours);
	  *ecode = EFailedToReadSeqErrRates;
	}
    }


  //now get the error cleaning information for each colour
  for (i=0; (i<binfo->number_of_colours) && (no_problem==true); i++)
    {
      no_problem = read_next_error_cleaning_object(fp, (binfo->ginfo->cleaning[i]) );
      if (no_problem==false)
	{
	  *ecode = EFailedToReadErrorCleaningInfo;
	}
    }

  return no_problem;

}

boolean read_next_error_cleaning_object(FILE* fp, ErrorCleaning* cl)
{
  int read;
  boolean no_problem=true;

  read = fread(&(cl->tip_clipping),sizeof(boolean),1,fp);
  if (read==0)
    {
      no_problem=false;
    }
  if (no_problem==true)
    {
      read = fread(&(cl->remv_low_cov_sups),sizeof(boolean),1,fp);
      if (read==0)
	{
	  no_problem=false;
	}
    }

  if (no_problem==true)
    {
      read = fread(&(cl->remv_low_cov_nodes),sizeof(boolean),1,fp);
      if (read==0)
	{
	  no_problem=false;
	}
    }
  if (no_problem==true)
    {
      read = fread(&(cl->cleaned_against_another_graph),sizeof(boolean),1,fp);
      if (read==0)
	{
	  no_problem=false;
	}
    }
  if (no_problem==true)
    {
      read = fread(&(cl->remv_low_cov_sups_thresh),sizeof(int),1,fp);
      if (read==0)
	{
	  no_problem=false;
	}
    }
  if (no_problem==true)
    {
      read = fread(&(cl->remv_low_cov_nodes_thresh),sizeof(int),1,fp);
      if (read==0)
	{
	  no_problem=false;
	}
    }

  if (no_problem==true)
    {
      read = fread(&(cl->len_name_of_graph_against_which_was_cleaned),sizeof(int),1,fp);
      if (read==0)
	{
	  no_problem=false;
	}
    }
  if (no_problem==true)
    {
      read = fread(cl->name_of_graph_against_which_was_cleaned,sizeof(char),cl->len_name_of_graph_against_which_was_cleaned,fp);
      if (read==0)
	{
	  no_problem=false;
	}
    }

  return no_problem;

}

/*
//return true if signature is readable
boolean query_binary(FILE * fp,int* binary_version, int* kmer_size, int* number_of_bitfields, int* number_of_colours_in_binary){
  int read;
  char magic_number[6];
  

  read = fread(magic_number,sizeof(char),6,fp);
  if (read>0 &&
      magic_number[0]=='C' &&
      magic_number[1]=='O' &&
      magic_number[2]=='R' &&
      magic_number[3]=='T' &&
      magic_number[4]=='E' &&
      magic_number[5]=='X' ){

    int version;
    read = fread(&version,sizeof(int),1,fp);
    if (read>0)
      {
	int kmer_size2;
	read = fread(&kmer_size2,sizeof(int),1,fp);
	if (read>0)
	  {
	    int num_bitfields;
	    read = fread(&num_bitfields,sizeof(int),1,fp);

	    if ( read>0 )
	      {
		int num_cols;
		read = fread(&num_cols,sizeof(int),1,fp);
		
		if ( read>0  )
		  { 
		    //int* binary_version, int* kmer_size, int* num_bitfields, int* number_of_colours_in_binar
		    *binary_version = version;
		    *kmer_size = kmer_size2;
		    *number_of_bitfields = num_bitfields;
		    *number_of_colours_in_binary = num_cols;
		  }
		else
		  {
		    return false;
		  }
	      }
	    else
	      {
		return false;
	      }
	  }
	else
	  {
	    return false;
	  }
      }
    else
      {
	return false;
      }
  }
  else
    {
      return false;
    }

  return true;
}
*/


//given a list of filenames, check they all exist, and return the number of them
int get_number_of_files_and_check_existence_from_filelist(char* filelist)
{
  
  int count=0;
  

  FILE* fp = fopen(filelist, "r");
  if (fp==NULL)
    {
      printf("Cannot open %s\n", filelist);
      exit(1);
    }


  char file[MAX_FILENAME_LENGTH+MAX_LEN_SAMPLE_NAME+1];
  file[0]='\0';
  
  while (feof(fp)==0)
    {
      if (fgets(file, MAX_FILENAME_LENGTH+MAX_LEN_SAMPLE_NAME, fp) != NULL)
	{
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(file, '\n')) != NULL)
	    {
	      *p = '\0';
	    }
	  
	  char temp1[MAX_FILENAME_LENGTH];
	  temp1[0]='\0';
	  strcpy(temp1, file);
	  char delims[] = "\t";
	  char* proper_filename = strtok(temp1, delims);

	  if (access(proper_filename, R_OK)==-1)
	    {
	      printf("Cannot access file %s listed in %s\n", proper_filename, filelist);
	      exit(1);
	    }
	  else
	    {
	      count++;
	    }

	}
    }
  fclose(fp);

  return count;
}

//assumes we already know how many files in the list, and have preallocated an array to hold the filenames
void get_filenames_from_list(char* filelist, char** array, int len)
{
  int count=0;

  FILE* fp = fopen(filelist, "r");
  if (fp==NULL)
    {
      printf("Cannot open %s\n", filelist);
      exit(1);
    }

  char file[MAX_FILENAME_LENGTH+MAX_LEN_SAMPLE_NAME+1];
  file[0]='\0';
  
  while ( (feof(fp)==0) && (count<len) )
    {
      if (fgets(file, MAX_FILENAME_LENGTH, fp) != NULL)
	{
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(file, '\n')) != NULL)
	    {
	      *p = '\0';
	    }


	  char temp1[MAX_FILENAME_LENGTH];
	  temp1[0]='\0';
	  strcpy(temp1, file);
	  char delims[] = "\t";
	  char* proper_filename = strtok(temp1, delims);

	  strcpy(array[count], proper_filename);
	  count++;
	}
    }
  fclose(fp);


  
}

// filename is a list of files, one for each colour (with optional second column of sample-ids). Check they all exists, there are not too many,
// ad that each of them contains a alist of valid binaries.
boolean check_colour_list(char* filename, int kmer)
{
  int num_cols_in_list = get_number_of_files_and_check_existence_from_filelist(filename);


  char** list_cols = malloc( sizeof(char*) * num_cols_in_list);
  if (list_cols==NULL)
    {
      printf("OOM. Give up can't even allocate space for the names of colourss\n");
      exit(1);
    }
  int i;
  for (i=0; i< num_cols_in_list; i++)
    {
      list_cols[i] = malloc(sizeof(char)*500);
      if (list_cols[i]==NULL)
	{
	  printf("OOM. Giveup can't even allocate space for the names of the colour file i = %d\n",i);
	  exit(1);
	}
    }

  get_filenames_from_list(filename, list_cols,num_cols_in_list);

  //first - none of these files should be a cortex binary - typical mistake, and worth checking for
  for (i=0; i< num_cols_in_list; i++)
    {
      FILE* fptr = fopen(list_cols[i], "r");
      if (fptr==NULL)
	{
	  printf("Cannot open %s from list %s\n", list_cols[i], filename);
	  exit(1);
	}
      GraphInfo* ginfo=graph_info_alloc_and_init();
      BinaryHeaderErrorCode ecode=EValid;
      BinaryHeaderInfo binfo;
      initialise_binary_header_info(&binfo, ginfo);
      int check = query_binary_NEW(fptr, &binfo, &ecode);
      if (check==true)
	{
	  printf("Error with input arguments. Summary: you have passed in a list of binaries, not a list of lists of binaries.\n--colour_list requires a list of files, each of which represent a colour.\n Inside each of those ");
	  printf("should be a list of cortex binaries. \nHowever your colour list %s contains this file %s which is itself a cortex binary not a list of binaries\n", filename, list_cols[i]);
	  exit(1);
	}
      else
	{
	  //check it contains only cortex binaries.

	  int count=0;

	  FILE* fp = fopen(list_cols[i], "r");
	  if (fp==NULL)
	    {
	      printf("Cannot open %s\n", list_cols[i]);
	      exit(1);
	    }

	  char file[MAX_FILENAME_LENGTH+1];
	  file[0]='\0';
	  
	  while ( feof(fp)==0 )
	    {
	      if (fgets(file, MAX_FILENAME_LENGTH, fp) != NULL)
		{
		  //remove newline from end of line - replace with \0
		  char* p;
		  if ((p = strchr(file, '\n')) != NULL)
		    {
		      *p = '\0';
		    }
		  FILE* fp_putative_binary=fopen(file, "r");
		  if (fp_putative_binary==NULL)
		    {
		      printf("Cannot open this file %s\n", file);
		      exit(1);
		    }
		  
		  GraphInfo* ginfo2=graph_info_alloc_and_init();
		  BinaryHeaderErrorCode ecode2=EValid;
		  BinaryHeaderInfo binfo2;
		  initialise_binary_header_info(&binfo2, ginfo2);
		  int check2 = query_binary_NEW(fp_putative_binary, &binfo2, &ecode2);
		  fclose(fp_putative_binary);
		  if (check2==false)
		    {
		      printf("Error with input arguments. --colour_list requires a list of files, each of which represent a colour. Inside each of those ");
		      printf("should be a list of cortex binaries. Therefore your colour list %s contains this file %s which should be a list of cortex binaries\n", filename, list_cols[i]);
		      printf("However it contains %s, which is not a cortex binary\n", file);
		      exit(1);
		    }
		  else if (binfo2.kmer_size != kmer)
		    {
		      printf("This binary %s is a kmer %d binary, and cannot be loaded into the main graph which has kmer %d\n", file,binfo.kmer_size, kmer);
		      exit(1);
		    }
		  else if (binfo2.number_of_bitfields != NUMBER_OF_BITFIELDS_IN_BINARY_KMER)
		    {
		      printf("Cannot load this binary %s, as it was built with max_kmer_size=%d, whereas the current graph has this set to %d. Incompatible binary.\n",
			     file, binfo2.number_of_bitfields, NUMBER_OF_BITFIELDS_IN_BINARY_KMER);
		      exit(1);
		    }
		  else if (binfo2.number_of_colours !=1)
		    {
		      printf("Input error. --colour_list requires a list of colour files, each one containing a list of single-colour binaries\n");
		      printf("This binary %s is not a single colour binary - it has %d colours.\n", file, binfo2.number_of_colours);
		      exit(1);
		    }
		  count++;
		  graph_info_free(ginfo2);
		}
	    }
	  fclose(fp);




	}


      fclose(fptr);
      graph_info_free(ginfo);
    }


  //cleanup
  for (i=0; i<num_cols_in_list; i++)
    {
      free(list_cols[i]);
    }
  free(list_cols);

  return true;
}


