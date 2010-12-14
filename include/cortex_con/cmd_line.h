/*
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
  cmd_line.h - manipulation of command line
  
*/

#ifndef CMD_LINE_H_
#define CMD_LINE_H_

#include <stdio.h>
#include <global.h>

typedef enum
 {
   FASTA = 0,
   FASTQ = 1,
   CTX   = 2,
   HASH  = 3,
 } FileFormat ;

typedef struct
{
  //core parameters 
  boolean verbose;
  int kmer_size;
  int bucket_size;
  int number_of_buckets_bits;
  

  //actions on/off
  boolean tip_clip;
  boolean remove_bubbles;
  boolean dump_binary;
  boolean dump_hash_table;
  boolean check_binary;

  boolean print_supernodes_fasta;
  boolean print_paths_fasta;
  boolean low_coverage_node_clip;
  boolean detect_bubbles;
  boolean detect_dirty_bubbles;
  boolean detect_forks;

  boolean print_coverages;
  boolean input_file; //if it is present
  boolean input_file_format_known; 
  boolean remove_low_coverage_edges; 
  boolean remove_low_coverage_supernodes; 
  boolean health_check_binary; 
  boolean remove_weak_edges;

  //parameters 
  FileFormat input_file_format;
  int quality_score_threshold;
  int node_coverage_threshold;
  char input_filename[LENGTH_FILENAME];
  char output_ctx_filename[LENGTH_FILENAME];
  char output_supernodes_fasta_filename[LENGTH_FILENAME];
  char output_paths_fasta_filename[LENGTH_FILENAME];
  char output_hash_filename[LENGTH_FILENAME];
  int tip_length;
  int detect_vars_delta;
  int detect_vars_branch_length;
  int quality_score_offset;
  int coverage_threshold;
  int remove_low_coverage_supernodes_threshold;
  int max_length_low_coverage_supernode;
  int min_coverage_threshold_remove_weak_edges;
  int max_read_len;
  int binary_version;
} CmdLine;

CmdLine parse_cmdline( int argc, char* argv[],int unit_size); 
int default_opts(CmdLine *);

#endif /* CMD_LINE_H_ */
