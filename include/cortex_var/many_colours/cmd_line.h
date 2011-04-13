
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
  cmd_line.h - manipulation of command line
  
*/


#ifndef CMD_LINE_H_
#define CMD_LINE_H_

#include <stdio.h>
#include <global.h>
#include <file_format.h>

#define MAX_FILENAME_LEN 1000
#define MAX_SUFFIX_LEN 100
#define  MAX_LEN_DETECT_BUB_COLOURINFO 500 //will be info of form 1,2,3/5,6,7,8,9 specifying how you call vars between colours
#define  MAX_COLOURS_ALLOWED_TO_MERGE 100 //arbitrary limit, can be increased
#define LEN_ERROR_STRING 200

//typedef enum
// {
//   FASTA = 0,
//   FASTQ = 1,
//   CTX   = 2,
//   UNSPECIFIED = 3,
   //VAR   = 3,
// } FileFormat ;

typedef struct
{

  int kmer_size;
  int bucket_size;
  int number_of_buckets_bits;
  int ref_colour;
  int clean_colour;
  int homopolymer_limit;
  int quality_score_threshold;
  int node_coverage_threshold;
  int quality_score_offset;
  int max_read_length;
  int max_var_len;

  int detect_bubbles1_first_colour_list[MAX_COLOURS_ALLOWED_TO_MERGE];//these are a set of colours to be considered as merged, for purposes of bubble calling
  int num_colours_in_detect_bubbles1_first_colour_list;
  int detect_bubbles1_second_colour_list[MAX_COLOURS_ALLOWED_TO_MERGE];// these are another set of colours to be considered as merged, for purposes of bubble calling
  int num_colours_in_detect_bubbles1_second_colour_list;

  int detect_bubbles2_first_colour_list[MAX_COLOURS_ALLOWED_TO_MERGE];
  int num_colours_in_detect_bubbles2_first_colour_list;
  int detect_bubbles2_second_colour_list[MAX_COLOURS_ALLOWED_TO_MERGE];
  int num_colours_in_detect_bubbles2_second_colour_list;

  int colour_overlap_first_colour_list[MAX_COLOURS_ALLOWED_TO_MERGE];
  int num_colours_in_colour_overlap_first_list;
  int colour_overlap_second_colour_list[MAX_COLOURS_ALLOWED_TO_MERGE];
  int num_colours_in_colour_overlap_second_list;


  int pd_colour_list[MAX_COLOURS_ALLOWED_TO_MERGE];
  int num_colours_in_pd_colour_list;

  char colour_list[MAX_FILENAME_LEN];
  int num_colours_in_input_colour_list;
  char multicolour_bin[MAX_FILENAME_LEN];
  int num_colours_in_multicol_bin;

  char se_list[MAX_FILENAME_LEN];
  char pe_list_lh_mates[MAX_FILENAME_LEN];
  char pe_list_rh_mates[MAX_FILENAME_LEN];
  char output_binary_filename[MAX_FILENAME_LEN];
  char output_supernodes[MAX_FILENAME_LEN];
  char output_detect_bubbles1[MAX_FILENAME_LEN];
  char output_detect_bubbles2[MAX_FILENAME_LEN];
  char path_divergence_caller_output_stub[MAX_FILENAME_LEN];
  char ref_chrom_fasta_list[MAX_FILENAME_LEN];
  char config[MAX_FILENAME_LEN];
  char covg_distrib_outfile[MAX_FILENAME_LEN];
  char readlen_distrib_outfile[MAX_FILENAME_LEN];
  char successively_dump_cleaned_colours_suffix[MAX_SUFFIX_LEN];
  char list_fastaq_to_align[MAX_FILENAME_LEN];

  boolean cut_homopolymers;
  boolean remove_pcr_dups;
  //boolean clip_tips;
  boolean remove_seq_errors;
  boolean print_colour_coverages;
  boolean load_colours_only_where_overlap_clean_colour;
  boolean successively_dump_cleaned_colours;
  boolean dump_binary;
  boolean print_supernode_fasta;
  boolean remove_low_coverage_nodes;
  boolean detect_bubbles1;
  boolean detect_bubbles2;
  boolean make_pd_calls;
  boolean using_ref;  
  boolean seq_file_format_known; 
  boolean input_colours;
  boolean input_multicol_bin;//flag - has it been input
  boolean input_seq; //flag - has it been input
  boolean dump_covg_distrib; //flag - shoud we dump covg distrib
  boolean dump_readlen_distrib; //flag - shoud we dump distrib of filtered read lengths
  boolean health_check;
  boolean align_given_list;
  boolean print_colour_overlap_matrix;
  //int detect_vars_delta;
  //int detect_vars_branch_length;
  //int quality_score_offset;
  FileFormat format_of_input_seq;
  FileFormat format_of_files_to_align;
} CmdLine;


int parse_cmdline_inner_loop(int argc, char* argv[], int unit_size, CmdLine* cmdline_ptr, char* error_string);
int check_cmdline(CmdLine* cmd_ptr, char* error_string);
CmdLine parse_cmdline( int argc, char* argv[],int unit_size); 
int default_opts(CmdLine *);
int get_numbers_from_comma_sep_list(char* list,  int* return_list, int max_len_return_list);
int parse_colourinfo_argument(CmdLine* cmd, char* arg, int len_arg, char* text_for_error_describing_which_option_this_is, int which_detect_bubbles);
int parse_commasep_list(CmdLine* cmd, char* arg, int len_arg, char* text_for_error_describing_which_option_this_is);

#endif /* CMD_LINE_H_ */
