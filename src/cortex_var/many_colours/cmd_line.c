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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <cmd_line.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <file_reader.h>
#include <file_format.h>
#include <err.h>

#define MAX_NUM_DIGITS_IN_COLOUR_ENTERED_ON_CMDLINE_COMMASEP 4
#define LEN_ERROR_STRING 200
#define MAX_LINE 500


boolean more_than_one_colour_in_list(char* file)
{

  if (file[0]=='\0')
    {
      return false;
    }

  FILE* fp = fopen(file, "r");
  if (fp==NULL)
    {
      printf("Cannot open %s\n", file);
      exit(1);
    }

  char filename[MAX_LINE];
  int count=0;
  while (feof(fp) ==0)
    {
      if (fgets(filename,MAX_LINE, fp) !=NULL)
        {
	  count++;
	}
    }
  fclose(fp);

  if (count>1)
    {
      return true;
    }
  else
    {
      return false;
    }
}



boolean more_than_one_colour_in_multicol_binary(char* file, int kmer_size)
{

  if (file[0]=='\0')
    {
      return false;
    }

  FILE* fp = fopen(file, "r");
  if (fp==NULL)
    {
      printf("Cannot open %s\n", file);
      exit(1);
    }

  int num_cols=0;
  int mean_readlens[NUMBER_OF_COLOURS];
  long long total_seqs[NUMBER_OF_COLOURS];
  int* mean_readlens_ptrs[NUMBER_OF_COLOURS];
  long long* total_seqs_ptrs[NUMBER_OF_COLOURS];
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      mean_readlens[i]=0;
      mean_readlens_ptrs[i]=&(mean_readlens[i]);
      total_seqs[i]=0;
      total_seqs_ptrs[i]=&(total_seqs[i]);
    }
  check_binary_signature(fp, kmer_size, BINVERSION, &num_cols, mean_readlens_ptrs, total_seqs_ptrs);
  fclose(fp);

  if (num_cols>1)
    {
      return true;
    }
  else
    {
      return false;
    }

}


const char* usage=
"Cortex, multicoloured target, cortex_var by Z. Iqbal (zam@well.ox.ac.uk) (primary contact for cortex_var) and M. Caccamo (mario.caccamo@bbsrc.ac.uk)\n" \
"\nusage: to build a single-colour binary:\ncortex_var --se_list <filename> --pe_list <filename> --format FASTQ --quality_score_threshold 5 --remove_pcr_duplicates --remove_seq_errors --dump_binary some_name.ctx\n" \
"\nusage: to build a multicolour graph from single-colour graphs and call variants between colours 1 and 2:\ncortex_var --colour_list <filename> --detect_bubbles1 1/2 --output_bubbles1 vars_between_cols1_and_2\n" \
"\nusage: to load a multicolour graph from single-colour graphs and call heterozygous variants in colour 0:\ncortex_var --colour_list <filename> --detect_bubbles1 0/0 --output_bubbles1 hets_in_colour_0\n" \
"\n" \
"   [--help] \t\t\t\t\t\t\t=\t This help screen.\n" \
"   [--colour_list FILENAME] \t\t\t\t\t=\t File of filenames, one per colour. n-th file is a list of\n\t\t\t\t\t\t\t\t\t single-colour binaries to be loaded into colour n.\n\t\t\t\t\t\t\t\t\t Cannot be used with --se_list or --pe_list \n" \
"   [--multicolour_bin FILENAME] \t\t\t\t=\t Filename of a multicolour binary, will be loaded first, into colours 0..n.\n\t\t\t\t\t\t\t\t\t If using --colour_list also, those will be loaded into subsequent colours, after this.\n" \
"   [--se_list FILENAME] \t\t\t\t\t=\t List of single-end fasta/q to be loaded into a single-colour graph.\n\t\t\t\t\t\t\t\t\t Cannot be used with --colour_list\n" \
"   [--pe_list FILENAME] \t\t\t\t\t=\t Two filenames, comma-separated: each is a list of paired-end fasta/q to be \n\t\t\t\t\t\t\t\t\t loaded into a single-colour graph. Lists are assumed to ordered so that \n\t\t\t\t\t\t\t\t\t corresponding paired-end fasta/q files are at the same positions in their lists.\n\t\t\t\t\t\t\t\t\t Currently Cortex only use paired-end information to remove\n\t\t\t\t\t\t\t\t\t PCR duplicate reads (if that flag is set).\n\t\t\t\t\t\t\t\t\t Cannot be used with --colour_list\n" \
"   [--kmer_size INT] \t\t\t\t\t\t=\t Kmer size (default 21). Must be an odd number.\n" \
"   [--mem_width INT] \t\t\t\t\t\t=\t Size of hash table buckets (default 100).\n" \
  //-g 
"   [--mem_height INT] \t\t\t\t\t\t=\t Number of buckets in hash table in bits (default 10). \n\t\t\t\t\t\t\t\t\t Actual number of buckets will be 2^(the number you enter)\n" \
  // -i
"   [--ref_colour INT] \t\t\t\t\t\t=\t Colour of reference genome.\n" \
  // -j
"   [--remove_pcr_duplicates] \t\t\t\t\t=\t Removes PCR duplicate reads by ignoring read pairs if both \n\t\t\t\t\t\t\t\t\t reads start at the same k-mer as a previous read,\n\t\t\t\t\t\t\t\t\t and single-ended reads if they start at the same k-mer as a previous read\n" \
  // -k
"   [--cut_homopolymers INT] \t\t\t\t\t=\t Breaks reads at homopolymers of length >= this threshold.\n\t\t\t\t\t\t (i.e. max homopolymer in filtered read==threshold-1, and New read starts after homopolymer)\n" \
  // -l
"   [--path_divergence_caller COMMA_SEP_COLOURS] \t\t=\t Make Path Divergence variant calls.\n\t\t\t\t\t\t\t\t\t Must specify colour of sample in which you want to find\n\t\t\t\t\t\t\t\t\t variants compared with the reference.\n\t\t\t\t\t\t\t\t\t This sample colour can be a union of colours (comma-separated list). \n\t\t\t\t\t\t\t\t\t Must also specify --ref_colour and --list_ref_fasta\n" \
  // -m
"   [--quality_score_threshold INT] \t\t\t\t=\t Filter for quality scores in the input file (default 0).\n" \
  // -n 
"   [--fastq_offset INT] \t\t\t\t\t=\t Default 33, for standard fastq.\n\t\t\t\t\t\t\t\t\t Some fastq directly from different versions of Illumina machines require different offsets.\n" \
  // -o
"   [--remove_seq_errors] \t\t\t\t\t=\t Remove tips + remove nodes if their coverage is more likely\n\t\t\t\t\t\t\t\t\t to be due to single-base sequencing error than sampling.\n" \
  // -p
"   [--dump_binary FILENAME] \t\t\t\t\t=\t Dump a binary file, with this name (after applying error-cleaning, if specified).\n" \
  // -q
"   [--output_supernodes FILENAME] \t\t\t\t\t=\t Dump a fasta file of all the supernodes (after applying all specified actions on graph).\n" \
  // -r
"   [--detect_bubbles1 COMMA_SEP_COLOURS/COMMA_SEP_COLOURS] \t=\t Find all the bubbles in the graph where the two branches lie in the specified colours\n\t\t\t\t\t\t\t\t\t (after applying all specified actions on graph).\n\t\t\t\t\t\t\t\t\t Typical use would be --detect_bubbles1 1/1 to find hets in colour 1,\n\t\t\t\t\t\t\t\t\t or --detect_bubbles1 0/1 to find homozygous non-reference bubbles where one branch is in colour 0 (and not colour1)\n\t\t\t\t\t\t\t\t\t and the other branch is in colour1 (but not colour 0).\n\t\t\t\t\t\t\t\t\t However, one can do more complex things:\n\t\t\t\t\t\t\t\t\t e.g.  --detect_bubbles1 1,2,3/4,5,6 to find bubbles where one branch is in 1,2 or 3 (and not 4,5 or 6)\n\t\t\t\t\t\t\t\t\t and the other branch in colour 4,5 or 6 (but not 1,2, or 3).\n" \
  // -s
"   [--output_bubbles1 FILENAME]\t\t\t\t\t=\t Bubbles called in detect_bubbles1 are dumped to this file.\n" \
  // -t
"   [--detect_bubbles2 COMMA_SEP_COLOURS/COMMA_SEP_COLOURS] \t=\t Exactly the same as detect_bubbles1, but allows you to make\n\t\t\t\t\t\t\t\t\t a second set of bubble calls immediately afterwards.\n\t\t\t\t\t\t\t\t\t This is to accomodate the common use-case where one loads a reference\n\t\t\t\t\t\t\t\t\t and an individual, and then wants to call homs, and hets.\n" \
  // -u
"   [--output_bubbles2 FILENAME]\t\t\t\t\t=\t Bubbles called in detect_bubbles2 are dumped to this file.\n" \
  // -v
"   [--format TYPE] \t\t\t\t\t\t=\t File format for input in se_list and pe_list. All files assumed to be of the same format.\n\t\t\t\t\t\t\t\t\t Type must be FASTQ, FASTA or CTX\n" \
  // -w
"   [--max_read_len] \t\t\t\t\t\t=\t Maximum read length over all input files. (Mandatory if fastq or fasta files are input.)\n" \
  // -x
"   [--print_colour_coverages]\t\t\t\t\t=\t Print coverages in all colours for supernodes and variants.\n" \
  // -y
"   [--max_var_len INT] \t\t\t\t\t\t=\t Maximum variant size searched for. Default 10kb. \n" \
 // -z
"   [--list_ref_fasta FILENAME] \t\t\t\t\t=\t File listing reference chromosome fasta file(s); needed for path-divergence calls. \n" \
  // -A
"   [--dump_covg_distribution FILENAME] \t\t\t\t=\t Print k-mer coverage distribution to the file specified\n" \
  // -B
"   [--remove_low_coverage_kmers INT] \t\t\t\t=\t Filter for kmers with coverage less than or equal to  threshold.\n"  \
  // -D
"   [--dump_filtered_readlen_distribution FILENAME] \t\t=\t Dump to file the distribution of \"effective\" read lengths after quality/homopolymer/PCR dup filters \n"  \
  // -E
"   [--load_colours_only_where_overlap_clean_colour INT] \t=\t Only load nodes from binary files in the colour-list when they overlap a\n\t\t\t\t\t\t\t\t\t specific colour (e.g. that contains a cleaned pooled graph);\n\t\t\t\t\t\t\t\t\t requires you to specify this particular colour. You must have loaded that colour beforehand, using --multicolour_bin\n"  \
  // -F
"   [--successively_dump_cleaned_colours SUFFIX] \t\t\t=\t Only to be used when also using --load_colours_only_where_overlap_clean_colour and --multicolour_bin\n\t\t\t\t\t\t\t\t\t Used to allow error-correction of low-coverage data on large numbers of individuals with large genomes.\n\t\t\t\t\t\t\t\t\t Requires the user specify a suffix which will be added to the names of cleaned binaries. See manual for details.\n"  \
  // -G
"   [--align FILENAME] \t\t\t\t\t\t=\t Aligns a list of fasta/q files to the graph, and prints coverage of each kmer in each read in each colour.\n\t\t\t\t\t\t\t\t\t Must also specify --align_input_format, and --max_read_len\n"  \
  // -H
"   [--align_input_format TYPE] \t\t\t\t\t=\t --align requires a list of fasta or fastq. This option specifies for format as LIST_OF_FASTQ or LIST_OF_FASTA\n"  \
  // -I
"   [--path_divergence_caller_output PATH_STUB]\t\t\t=\t Specifies the path and beginning of filenames of Path Divergence caller output files.\n\t\t\t\t\t\t\t\t\t One output file will be created per reference fasta listed in --list_ref_fasta\n" \
  // -J
"   [--colour_overlaps COMMA_SEP_COLOURS/COMMA_SEP_COLOURS]\t=\t Compares each coloured subgraph in the first list with all of the \n\t\t\t\t\t\t\t\t\t coloured subgraphs in the second list. Outputs a matrix to stdout;\n\t\t\t\t\t\t\t\t\t (i,j)-element is the number of nodes in both colour-i (on first list)\n\t\t\t\t\t\t\t\t\t and colour-j (on second list).\n" \
  "\n";






//assumes you have already checked the lenghts of the arrays  
void copy_list(int* from_list, int len_from_list, int* to_list, int len_to_list)
{
  int i;
  for (i=0; i<len_from_list; i++)
    {
      to_list[i]=from_list[i];
    }
}





void set_string_to_null(char* str, int len)
{
  int i;
  for (i=0; i<len; i++)
    {
      str[i]='\0';
    }
}


void initialise_int_list(int* list, int len)
{
  int i;
  for (i=0; i<len; i++)
    {
      list[i]=-1;
    }
}


void initialise_longlong_list(long long* list, int len)
{
  int i;
  for (i=0; i<len; i++)
    {
      list[i]=-1;
    }
}





int default_opts(CmdLine * c)
{

  c->kmer_size = 21;
  c->bucket_size = 100;
  c->number_of_buckets_bits = 10;
  c->ref_colour=-1;
  c->homopolymer_limit=-1;
  c->quality_score_threshold=0;
  c->node_coverage_threshold=0;
  c->quality_score_offset = 33;//standard fastq, not illumina v-whatever fastq  
  c->max_read_length = 0;
  c->max_var_len = 10000;
  c->clean_colour=NUMBER_OF_COLOURS+1;//default to an impossible value
  c->num_colours_in_detect_bubbles1_first_colour_list=0;
  initialise_int_list(c->detect_bubbles1_first_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);
  c->num_colours_in_detect_bubbles1_second_colour_list=0;
  initialise_int_list(c->detect_bubbles1_second_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);
  
  c->num_colours_in_detect_bubbles2_first_colour_list=0;
  initialise_int_list(c->detect_bubbles2_first_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);  
  c->num_colours_in_detect_bubbles2_second_colour_list=0;
  initialise_int_list(c->detect_bubbles2_second_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);

  c->num_colours_in_colour_overlap_first_list=0;
  initialise_int_list(c->colour_overlap_first_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);
  c->num_colours_in_colour_overlap_second_list=0;
  initialise_int_list(c->colour_overlap_second_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);

  c->num_colours_in_pd_colour_list=0;
  initialise_int_list(c->pd_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);
  

  c->num_colours_in_input_colour_list=0;

  //filenames/strings
  set_string_to_null(c->colour_list, MAX_FILENAME_LEN);
  set_string_to_null(c->multicolour_bin,MAX_FILENAME_LEN);
  set_string_to_null(c->se_list,MAX_FILENAME_LEN);
  set_string_to_null(c->pe_list_lh_mates,MAX_FILENAME_LEN);
  set_string_to_null(c->pe_list_rh_mates,MAX_FILENAME_LEN);
  set_string_to_null(c->output_binary_filename,MAX_FILENAME_LEN);
  set_string_to_null(c->output_supernodes,MAX_FILENAME_LEN);
  set_string_to_null(c->output_detect_bubbles1,MAX_FILENAME_LEN);
  set_string_to_null(c->output_detect_bubbles2,MAX_FILENAME_LEN);
  set_string_to_null(c->path_divergence_caller_output_stub, MAX_FILENAME_LEN);
  set_string_to_null(c->ref_chrom_fasta_list,MAX_FILENAME_LEN);
  set_string_to_null(c->config,MAX_FILENAME_LEN);
  set_string_to_null(c->covg_distrib_outfile,MAX_FILENAME_LEN);
  set_string_to_null(c->readlen_distrib_outfile,MAX_FILENAME_LEN);
  set_string_to_null(c->successively_dump_cleaned_colours_suffix, MAX_SUFFIX_LEN);
  set_string_to_null(c->list_fastaq_to_align, MAX_FILENAME_LEN);

  //booleans
  c->cut_homopolymers = false;
  c->remove_pcr_dups=false;
  //c->clip_tips=false;
  c->remove_seq_errors=false;
  c->print_colour_coverages=false;
  c->dump_binary=false;
  c->print_supernode_fasta=false;
  c->remove_low_coverage_nodes=false;
  c->detect_bubbles1=false;
  c->detect_bubbles2=false;
  c->make_pd_calls=false;
  c->using_ref=false;
  c->seq_file_format_known=false;
  c->input_colours=false;
  c->input_multicol_bin = false;
  c->input_seq = false;
  c->format_of_input_seq=UNSPECIFIED;
  c->dump_covg_distrib=false;
  c->dump_readlen_distrib=false;
  c->health_check=false;
  c->load_colours_only_where_overlap_clean_colour=false;
  c->successively_dump_cleaned_colours=false;
  c->align_given_list=false;
  c->print_colour_overlap_matrix=false;
  c->format_of_files_to_align=UNSPECIFIED;
  return 1;
}


//inner loop called by the parse_cmdline function, which returns various error codes.
int parse_cmdline_inner_loop(int argc, char* argv[], int unit_size, CmdLine* cmdline_ptr, char* error_string)
{
  int opt;
  int longopt_index;
  
  static struct option long_options[] = {
    {"help",no_argument,NULL,'h'},
    {"colour_list",required_argument, NULL, 'a'},
    {"multicolour_bin", required_argument, NULL, 'b'},
    {"se_list",required_argument, NULL, 'c'},
    {"pe_list",required_argument, NULL, 'd'},
    {"kmer_size", required_argument, NULL, 'e'},
    {"mem_width", required_argument, NULL, 'f'},
    {"mem_height",required_argument, NULL, 'g'},
    {"ref_colour",required_argument, NULL, 'i'},
    {"remove_pcr_duplicates",no_argument,NULL,'j'},
    {"cut_homopolymers",required_argument, NULL, 'k'},
    {"path_divergence_caller",required_argument,NULL,'l'},
    {"quality_score_threshold",required_argument, NULL, 'm'},
    {"fastq_offset",required_argument, NULL, 'n'},
    {"remove_seq_errors",no_argument,NULL,'o'},
    {"dump_binary",required_argument,NULL,'p'},
    {"output_supernodes",required_argument,NULL,'q'},
    {"detect_bubbles1",required_argument, NULL, 'r'},
    {"output_bubbles1",required_argument, NULL, 's'},    
    {"detect_bubbles2",required_argument, NULL, 't'},
    {"output_bubbles2",required_argument, NULL, 'u'},    
    {"format",required_argument,NULL,'v'},
    {"max_read_len",required_argument,NULL,'w'},
    {"print_colour_coverages",no_argument,NULL,'x'},
    {"max_var_len",required_argument,NULL,'y'},
    {"list_ref_fasta",required_argument,NULL,'z'},
    {"dump_covg_distribution",required_argument,NULL,'A'},
    {"remove_low_coverage_kmers",required_argument,NULL,'B'},
    {"health_check",no_argument,NULL,'C'},//hidden
    {"dump_filtered_readlen_distribution",required_argument,NULL,'D'},
    {"load_colours_only_where_overlap_clean_colour",required_argument,NULL,'E'},
    {"successively_dump_cleaned_colours",required_argument,NULL,'F'},
    {"align",required_argument,NULL,'G'},
    {"align_input_format",required_argument,NULL,'H'},
    {"path_divergence_caller_output",required_argument,NULL,'I'},
    {"colour_overlaps",required_argument,NULL,'J'},
    {0,0,0,0}	
  };
  

  //do not change this! Only really matters for testing, but getopt_long uses 
  //variables which are not local to this function, so when we run tests, and call this function
  // repeatedly, those variable values are carried across. The following line resets this.
  optind=1;
  
 
  opt = getopt_long(argc, argv, "ha:b:c:d:e:f:g:i:jk:l:m:n:op:q:r:s:t:u:v:w:xy:z:A:B:C:D:E:F:G:H:I:J:", long_options, &longopt_index);

  while ((opt) > 0) {
	       
    //Parse the default options
    switch(opt) {

    case 'h':
      {
	printf("***********************\n");
	printf("Cortex version %d.%d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION,SUBSUBSUBVERSION);
	printf("Compiled with support for %d colours\n", NUMBER_OF_COLOURS);
	printf("***********************\n");

	printf("%s",usage);
	exit(0);
	break;
      }

    case 'a'://colour list
      {
	if (optarg==NULL)
	  errx(1,"[--colour_list] option requires a filename [file of filenames, one for each colour]");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    //how any colours in this filelist?
	    int num_cols_in_input_list=get_number_of_files_and_check_existence_from_filelist(optarg);
	    if (num_cols_in_input_list==0)
	      {
		errx(1, "[--colour_list] filename %s contains nothing", optarg);
	      }
	    cmdline_ptr ->num_colours_in_input_colour_list=num_cols_in_input_list;
	    strcpy(cmdline_ptr->colour_list,optarg);
	    printf("colour list is set to %s\n", cmdline_ptr->colour_list);
	  }
	else
	  {
	    errx(1,"[--colour_list] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,R_OK)==-1){
	  errx(1,"[--colour_list] filename [%s] cannot be accessed",optarg);
	}

	cmdline_ptr->input_colours = true;
	break; 
      }


    case 'b'://multicolour bin
      {
	if (optarg==NULL)
	  errx(1,"[--multicolour_bin] option requires a filename (of a multicolour binary)");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->multicolour_bin,optarg);
	  }
	else
	  {
	    errx(1,"[--multicolour_bin] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,R_OK)==-1){
	  errx(1,"[--multicolour_bin] filename [%s] cannot be accessed",optarg);
	}
	cmdline_ptr->input_multicol_bin = true;
	//we set num_colours_in_multicol_bin later on - need to open the file and check signature, and cant do that til we know what kmer, and cant be sure in this case 
	//that we have already parsed the --kmer_size case already
	break; 
      }





    case 'c': //se_list
      {
	if (optarg==NULL)
	  errx(1,"[--se_list] option requires a filename (listing single_ended fasta/q)");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->se_list,optarg);
	  }
	else
	  {
	    errx(1,"[--se_list] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,R_OK)==-1){
	  errx(1,"[--se_list] filename [%s] cannot be accessed",optarg);
	}
	cmdline_ptr->input_seq = true;
	break; 
      }


    case 'd': //pe_list, comma separated - two files, listing matched pairs in the same order
      {
	if (optarg==NULL)
	  errx(1,"[--pe_list] option requires two filenames, comma separated, listing paired-end fasta/q files (in matching order)");
	
	if (strlen(optarg)<2*MAX_FILENAME_LEN)
	  {
	    char* filename=NULL;
	    char delims[] = ",";
	    char temp1[MAX_FILENAME_LEN];
	    temp1[0]='\0';
	    strcpy(temp1, optarg);
	    filename = strtok(temp1, delims );
	    if (filename==NULL)
	      {
		errx(1,"[--pe_list] option requires two filenames, comma separated, listing paired-end fasta/q files (in matching order)");
	      }
	    strcpy(cmdline_ptr->pe_list_lh_mates,filename);
	    filename = strtok( NULL, delims );
	    if (filename==NULL)
	      {
		errx(1,"[--pe_list] option requires two filenames, comma separated, listing paired-end fasta/q files (in matching order)");
	      }
	    strcpy(cmdline_ptr->pe_list_rh_mates,filename);

	  }
	else
	  {
	    errx(1,"[--pe_list] filename too long [%s]",optarg);
	  }
	
	if (access(cmdline_ptr->pe_list_lh_mates,R_OK)==-1){
	  errx(1,"[--pe_list] filename [%s] cannot be accessed",cmdline_ptr->pe_list_lh_mates);
	}
	if (access(cmdline_ptr->pe_list_rh_mates,R_OK)==-1){
	  errx(1,"[--pe_list] filename [%s] cannot be accessed",cmdline_ptr->pe_list_rh_mates);
	}
	cmdline_ptr->input_seq = true;
	break; 
      }

      case 'e': //kmer size
	{
	  if (optarg==NULL)
	    errx(1,"[--kmer_size] option requires int argument [kmer size]");	    
	  cmdline_ptr->kmer_size = atoi(optarg);
	  
	  if (cmdline_ptr->kmer_size == 0)
	    errx(1,"[--kmer_size] option requires int argument bigger than 0");

	  if (cmdline_ptr->kmer_size % 2 == 0)
	    errx(1,"[--kmer_size] option requires int argument which is not divisible by 2");
	  
	  break;
	}

    case 'f': //--mem_width  = bucket size
      {
	if (optarg==NULL)
	  errx(1,"[--mem_width] option requires int argument [hash table bucket size]");
	cmdline_ptr->bucket_size = atoi(optarg);
	
	if (cmdline_ptr->bucket_size == 0) 
	  errx(1,"[--mem_width] option requires argument bigger than 0");
	if (cmdline_ptr->bucket_size > 32000) 
	  errx(1,"[--mem_width] option requires argument less than 32000 - recommend increasing the number of buckets instead of making each one very deep");
	break;
      }

    case 'g': //mem_height = number of buckets
      {
	if (optarg==NULL)
	  errx(1,"[--mem_height] option requires int argument [hash table,  number of buckets in bits - ie 2^(this number) is the number of buckets]");
	cmdline_ptr->number_of_buckets_bits = atoi(optarg);      
	break;
      }

    case 'i': //ref_colour
      {
	if (optarg==NULL)
	  errx(1,"[--ref_colour] option requires int argument [which colour is the reference - i.e. what position in colour_list (count starts from 0)]");
	if (optarg<0)
	  errx(1,"[--ref_colour] option requires positive argument [which colour is the reference - i.e. what position in colour_list (count starts from 0)]");
	if (atoi(optarg)>NUMBER_OF_COLOURS)
	  errx(1,"[--ref_colour] requires you specify the colour of the reference, and this must be less than the maximum number of colours allowed by Cortex. This maximum number is fixed at compile-time. See the Manual to find how to reset this, and check the number you have entered : %s is correct.", optarg);

	cmdline_ptr->using_ref=true;
	cmdline_ptr->ref_colour = atoi(optarg);      
	break;
      }


    case 'j': //remove_pcr_duplicates
      {
	cmdline_ptr->remove_pcr_dups = true;
	break;
      }
    case 'k': //cut_homopolymers
      {

	if (optarg==NULL)
	  errx(1,"[--cut_homopolymers] option requires positive integer  argument [maximum allowed length of homopolymer in read]");
	if (optarg<0)
	  errx(1,"[--cut_homopolymers] option requires positive integer  argument [maximum allowed length of homopolymer in read]");

	cmdline_ptr->cut_homopolymers  = true;
	cmdline_ptr->homopolymer_limit = atoi(optarg);      
	break;
      }
    case 'l': //path divergence caller
      {

	if (optarg==NULL)
	  errx(1,"[--path_divergence_caller] option requires at least one colour. If >1, must be comma separated.");
	
	parse_commasep_list(cmdline_ptr, optarg, strlen(optarg), "[--path_divergence_caller]");
	cmdline_ptr->make_pd_calls = true;
	break; 


      }

    case 'm'://quality threshold
      {
	if (optarg==NULL)
	  errx(1,"[--quality_score_threshold] option requires int argument [quality score threshold]");
	cmdline_ptr->quality_score_threshold = atoi(optarg);

	
	if ( (cmdline_ptr->format_of_input_seq != FASTQ) && (cmdline_ptr->format_of_input_seq != UNSPECIFIED) )
	  {
	    errx(1,"[--quality_score_threshold] has been called, implying the input format is FASTQ, but another command-line input has specified the input format is NOT FASTQ. Inconsistent args");
	  }
	cmdline_ptr->format_of_input_seq = FASTQ;
	
	if (cmdline_ptr->quality_score_threshold == 0)
	  errx(1,"[--quality_score_threshold] option requires int argument bigger than 0");
	break;    
      }

    case 'n':
      {

	if (optarg==NULL)
          errx(1,"[--fastq_offset] option requires int argument");
        cmdline_ptr->quality_score_offset = atoi(optarg);

	break;    
      }


    case 'o': //use coverage and topology to remove sequencing errors. Not compatible with any other error correction
      {
	cmdline_ptr->remove_seq_errors = true;
	break;
      }

    case 'p'://output of binary ctx
      {
	cmdline_ptr->dump_binary=true;
	if (optarg==NULL)
	  errx(1,"[--dump_binary] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->output_binary_filename,optarg);
	    cmdline_ptr->dump_binary=true;
	  }
	else
	  {
	    errx(1,"[--dump_binary] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0){
	  errx(1,"[--dump_binary] filename [%s] exists!",optarg);
	}
	
	break; 
      }


    case 'q': //output supernode contigs
      {
	cmdline_ptr->print_supernode_fasta = true;
	if (optarg==NULL)
	  errx(1,"[--output_supernodes] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->output_supernodes,optarg);
	  }
	else
	  {
	    errx(1,"[--output_supernodes] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0)
	  {
	    errx(1,"[--output_supernodes] filename [%s] already exists. Exiting, to prevent overwriting.",optarg);
	  }
	if (optarg[0]=='-')
	  {
	    errx(1, "[--output_supernodes] requires a filename, but finds this: [%s] starts with -. Have you omitted the filename?\n", optarg);
	  }
	break; 
      }

      
    case 'r': //detect_bubbles1 - will give something of this form 1,2,3/5,1,7 to show the colours it wants to distinguish
      {
	if (optarg==NULL)
	  errx(1,"[--detect_bubbles1] option requires two sets of comma-separated integers, separated by a forward-slash. eg 1/1, or 1,2,3/3,4,5");
	
	int ret = parse_colourinfo_argument(cmdline_ptr, optarg, strlen(optarg), "[--detect_bubbles1] ", 1);
	if (ret==-1)
	  {
	    errx(1, "Problem with  cmd line argument for [--detect_bubbles1]");
	  }
	cmdline_ptr->detect_bubbles1=true;

	break; 
      }
    case 's': //output file for detect_bubbles1
      {
	if (optarg==NULL)
	  errx(1,"[--output_bubbles1] option requires a filename");

	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->output_detect_bubbles1,optarg);
	  }
	else
	  {
	    errx(1,"[--output_bubbles1] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0){
	  errx(1,"[--output_bubbles1] filename [%s] exists!",optarg);
	}
	break;
      }


    case 't': //detect_bubbles2 - will give something of this form 1,2,3/5,1,7 to show the colours it wants to distinguish
      {
	if (optarg==NULL)
	  errx(1,"[--detect_bubbles2] option requires two sets of comma-separated integers, separated by a forward-slash. eg 1/1, or 1,2,3/3,4,5");
	
	int ret = parse_colourinfo_argument(cmdline_ptr, optarg, strlen(optarg), "[-t | --detect_bubbles2] ", 2);
	if (ret==-1)
	  {
	    errx(1, "Problem with  cmd line argument for [-t | --detect_bubbles2]");
	  }
	cmdline_ptr->detect_bubbles2=true;
	break; 
      }
    case 'u': //output file for detect_bubbles2
      {
	if (optarg==NULL)
	  errx(1,"[--output_bubbles2] option requires a filename");

	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->output_detect_bubbles2,optarg);
	  }
	else
	  {
	    errx(1,"[--output_bubbles2] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0){
	  errx(1,"[--output_bubbles2] filename [%s] exists!",optarg);
	}
	break;        
      } 
     
    case 'v': //file format - either fasta, fastq or ctx.
      {

	if (optarg==NULL)
	  errx(1,"[--format] option requires argument FASTQ, FASTA or CTX");

	if ( (strcmp(optarg, "FASTA") !=0) && (strcmp(optarg, "FASTQ") !=0) && (strcmp(optarg, "CTX") !=0)  )
	  {
	    errx(1,"[--format] option requires argument FASTQ, FASTA or CTX");
	  }
	cmdline_ptr->seq_file_format_known=true;
	
	if (strcmp(optarg, "FASTA") ==0)
	  {
	    cmdline_ptr->format_of_input_seq=FASTA;
	  }
	else if (strcmp(optarg, "FASTQ") ==0)
	  {
	    cmdline_ptr->format_of_input_seq=FASTQ;
	  } 
	else if (strcmp(optarg, "CTX") ==0)
	  {
	    cmdline_ptr->format_of_input_seq=CTX;
	  }


	break ;
      }      
    case 'w'://max_read_length
      {
	if (optarg==NULL)
	  errx(1,"[--max_read_len] option requires (positive) integer argument");
	if (atoi(optarg)<0)
	  {
	    errx(1,"[--max_read_len] option requires (positive) integer argument");
	  }
	cmdline_ptr->max_read_length = atoi(optarg);

	break ;
      }
    case 'x':
      {
	cmdline_ptr->print_colour_coverages=true;
	break ;
      }
    case 'y':
      {
	if (optarg==NULL)
	  errx(1,"[--max_var_len] option requires (positive) integer argument");
	if (atoi(optarg)<0)
	  {
	    errx(1,"[--max_var_len] option requires (positive) integer argument");
	  }
	cmdline_ptr->max_var_len = atoi(optarg);

	break ;
      }

    case 'z':
      {

	if (optarg==NULL)
	  errx(1,"[--list_ref_fasta] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->ref_chrom_fasta_list,optarg);
	  }
	else
	  {
	    errx(1,"[--list_ref_fasta] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)!=0){
	  errx(1,"[--list_ref_fasta] filename [%s] not found or cannot be opened",optarg);
	}


	break ;
      }
    case 'A':
      {
	
	if (optarg==NULL)
	  errx(1,"[--dump_covg_distribution] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->covg_distrib_outfile,optarg);
	    cmdline_ptr->dump_covg_distrib=true;
	  }
	else
	  {
	    errx(1,"[--dump_covg_distribution] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0){
	  errx(1,"[--dump_covg_distribution] filename [%s] already exists",optarg);
	}
	

	break ;
       }


   case 'B':
      {

	if (optarg==NULL)
	  {
	    printf("[--remove_low_coverage_kmers] option requires int argument [node coverage cut off]");
	    exit(1);
	  }
	cmdline_ptr->node_coverage_threshold = atoi(optarg);
	cmdline_ptr->remove_low_coverage_nodes = true;
	
	if (cmdline_ptr->node_coverage_threshold <= 0)
	  {
	    printf("[--remove_low_coverage_kmers] option requires int argument bigger than 0");
	    exit(1);
	  }
	break;    
	
      }

    case 'C':
      {
	cmdline_ptr->health_check=true;
	break;    
      }
    case 'D':
      {

	if (optarg==NULL)
	  errx(1,"[--dump_filtered_readlen_distribution FILENAME] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->readlen_distrib_outfile,optarg);
	    cmdline_ptr->dump_readlen_distrib=true;
	  }
	else
	  {
	    errx(1,"[--dump_filtered_readlen_distribution] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0){
	  errx(1,"[--dump_filtered_readlen_distribution] filename [%s] already exists",optarg);
	}

	break ;
      }
    case 'E':
      {
        if (optarg==NULL)
	  {
            printf("[--load_colours_only_where_overlap_clean_colour] option requires int argument [which is the clean/preferred colour]");
            exit(1);
          }
	cmdline_ptr->load_colours_only_where_overlap_clean_colour=true;
	cmdline_ptr->clean_colour=atoi(optarg);

        if (cmdline_ptr->clean_colour<0)
          {
            printf("[--load_colour_list_only_where_overlap_clean_colour] option requires int argument bigger than 0");
            exit(1);
          }
	else if  (cmdline_ptr->clean_colour>=NUMBER_OF_COLOURS)
          {
            printf("[--load_colour_list_only_where_overlap_clean_colour] option requires int argument <=  max colour limit");
            exit(1);
          }
        break;

      }
    case 'F':
      {
	
	if (optarg==NULL)
	  errx(1,"[--successively_dump_cleaned_colours] option requires a suffix - one cleaned binary will be dumped per colour in your colourlist, and this suffix will be added to the filename");
	
	if (strlen(optarg)<MAX_SUFFIX_LEN)
	  {
	    strcpy(cmdline_ptr->successively_dump_cleaned_colours_suffix,optarg);
	    cmdline_ptr->successively_dump_cleaned_colours=true;
	  }
	else
	  {
	    errx(1,"[--successively_dump_cleaned_colours] suffix too long [%s]",optarg);
	  }
	break ;
       }

    case 'G'://align a list of fasta/q to the graph and print colour coverages
      {
	
	if (optarg==NULL)
	  errx(1,"[--align FILENAME] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->list_fastaq_to_align,optarg);
	    cmdline_ptr->align_given_list=true;
	  }
	else
	  {
	    errx(1,"[--align] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)!=0){
	  errx(1,"[--align_fasta] filename [%s] cannot be found/opened",optarg);
	}
	
	break ;
      }
    case 'H': //file format of files that we will align to the graph - either fasta or  fastq
      {
	
	if (optarg==NULL)
	  errx(1,"[--align_input_format] option requires argument LIST_OF_FASTA or LIST_OF_FASTQ");
	
	if ( (strcmp(optarg, "LIST_OF_FASTA") !=0) && (strcmp(optarg, "LIST_OF_FASTQ") !=0)  )
	  {
	    errx(1,"[--align_input_format] option requires argument LIST_OF_FASTA or LIST_OF_FASTQ");
	  }
	
	if (strcmp(optarg, "LIST_OF_FASTA") ==0)
	  {
	    cmdline_ptr->format_of_files_to_align=FASTA;
	  }
	else if (strcmp(optarg, "LIST_OF_FASTQ") ==0)
	  {
	    cmdline_ptr->format_of_files_to_align=FASTQ;
	  } 
	break ;
      }
    case 'I':
	{
	  if (optarg==NULL)
	    errx(1,"[--path_divergence_caller_output] option requires a filename");

	  if (strlen(optarg)<MAX_FILENAME_LEN-30)
	    {
	      strcpy(cmdline_ptr->path_divergence_caller_output_stub,optarg);
	    }
	  else
	    {
	      errx(1,"[--path_divergence_caller_output] filename too long [%s]",optarg);
	    }

	  break;
	}
    case 'J':
	{
	  if (optarg==NULL)
	    errx(1,"[--colour_overlaps] option requires two sets of comma-separated integers, separated by a forward-slash. eg 1/1, or 1,2,3/3,4,5");


	  int ret = parse_colourinfo_argument(cmdline_ptr, optarg, strlen(optarg), "[--colour_overlaps] ", 3);
	if (ret==-1)
	  {
	    errx(1, "Problem with  cmd line argument for [--colour_overlaps]");
	  }
	cmdline_ptr->print_colour_overlap_matrix=true;

	  break;
	}

    }      


    opt = getopt_long(argc, argv, "ha:b:c:d:e:f:g:i:jk:lm:n:opqr:s:t:u:v:w:xy:z:A:B:C:D:E:F:G:H:I:J:", long_options, &longopt_index);

  }
  return 0;
}


int check_cmdline(CmdLine* cmd_ptr, char* error_string)
{
  //must have specified either sequence input, or binary. Binary may be colour_list, or a multicolour graph, or both
  if ( (cmd_ptr->input_seq ==true) && ( (cmd_ptr->input_colours==true) || (cmd_ptr->input_multicol_bin==true) ) )
    {
      char tmp[] = "Cannot pass in both binary and sequence (fasta/fastq) data\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if (cmd_ptr->input_colours ==true)
    {
      check_colour_list(cmd_ptr->colour_list, cmd_ptr->kmer_size);
      
      if ( (cmd_ptr->successively_dump_cleaned_colours==false) && (cmd_ptr->num_colours_in_input_colour_list>NUMBER_OF_COLOURS) )
	{
	  printf("You are trying to load more colours than you have compiled cortex for. Your colour_list contains %d colours, but cortex_var is compiled for a maximum of %d\n",
		 cmd_ptr->num_colours_in_input_colour_list, NUMBER_OF_COLOURS);
	  exit(1);
	}

   }

  

  if ( (cmd_ptr->load_colours_only_where_overlap_clean_colour==true)
       &&
       (cmd_ptr->input_multicol_bin==false)
       )
    {
     	  char tmp[]="If you want to specify --load_colours_only_where_overlap_clean_colour, then you need to have a prespecified colour already in the graph before you load your colour list. This must therefore be loaded by --multicolour_bin, but you have not specified that.\n";
	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      printf("coding error - this string is too long:\n%s\n", tmp);
	      exit(1);
	    }
	  strcpy(error_string, tmp);
	  return -1;
    }
  

  if ( (cmd_ptr->successively_dump_cleaned_colours==true) && (cmd_ptr->load_colours_only_where_overlap_clean_colour==false) )
    {
     	  char tmp[]="If you specify --successively_dump_cleaned_colours, you must also specify --load_colours_only_where_overlap_clean_colour\n";
	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      printf("coding error - this string is too long:\n%s\n", tmp);
	      exit(1);
	    }
	  strcpy(error_string, tmp);
	  return -1;
      
    }

  if ( (cmd_ptr->load_colours_only_where_overlap_clean_colour==true)
       &&
       (cmd_ptr->input_multicol_bin==true)
       )
    {
      //check that the clean_colour specified is consistent with having to be loaded from multicolour_bin
      // ie if the multicolour_bin has n colours, then clean_colour must be in [0,n-1]
      FILE* fp = fopen(cmd_ptr->multicolour_bin, "r");
      if (fp==NULL)
	{
	  char tmp[]="Cannot open the specified multicolour binary\n";
	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      printf("coding error - this string is too long:\n%s\n", tmp);
	      exit(1);
	    }
	  strcpy(error_string, tmp);
	  return -1;
	}

      int binary_version_in_header;
      int kmer_in_header;
      int num_bf;
      int num_cols;
      boolean check = query_binary(fp, &binary_version_in_header, &kmer_in_header, &num_bf, &num_cols);
      if (cmd_ptr->clean_colour >= num_cols)
	{
	  char tmp[]="Specified clean colour number is > number of colours in the multicolour binary. Consult the manual.\n";
	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      printf("coding error - this string is too long:\n%s\n", tmp);
	      exit(1);
	    }
	  strcpy(error_string, tmp);
	  return -1;
	}
    }


  if (cmd_ptr->input_seq ==false) 
    {
      if ( (cmd_ptr->input_colours==false) && (cmd_ptr->input_multicol_bin==false) )
	{
	  char tmp[]="Must specify either binary (--colour_list or --multicol_bin), or sequence data (--se_list and/or --pe_list) as input \n";
	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      printf("coding error - this string is too long:\n%s\n", tmp);
	      exit(1);
	    }
	  strcpy(error_string, tmp);
	  return -1;
	}
    }


  
  //check kmer_size is odd
  if (cmd_ptr->kmer_size % 2 == 0)
    {
      char tmp[]="[-k/--kmer_size] must be odd";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }
  
  
  if ( (cmd_ptr->using_ref==true) && (cmd_ptr->remove_seq_errors==true) )

    {

      char tmp[]="Error-correction techniques should not be used on a graph including a reference genome.\nWe recommend error-correcting single-colour binaries of sequencing data, and then loading them into a multicolour graph including a reference\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }


  if ( (cmd_ptr->dump_readlen_distrib==true) && (cmd_ptr->format_of_input_seq!=FASTQ) )
    {
      char tmp[]="--dump_filtered_readlen_distribution is only supported for fastq input";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }


  /*
  //prevent removal of sequencing errors if we are loading incto >1 colour, or loading a binary with >1 colour

  if ( ( (more_than_one_colour_in_list(cmd_ptr->colour_list)==true)
	 ||
	 (more_than_one_colour_in_multicol_binary(cmd_ptr->multicolour_bin, cmd_ptr->kmer_size)==true)
	 )
       &&
       (cmd_ptr->remove_seq_errors==true)
       )
    {
      char tmp[]="Incompatible cmd line arguments. Error cleaning options only allowed on single-colour graphs\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;

    }
  */

  
  
  //check only call detect_bubbles2 if already have called detect_bubbles1
  if ( (cmd_ptr->detect_bubbles2==true) && (cmd_ptr->detect_bubbles1==false) )
    {

      char tmp[]="Do not specify --detect_bubbles2 unless you have already specified --detect_bubbles1. Consult manual for further information\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }


  //if specify detect_bubbles, must give output file
  if ( (cmd_ptr->detect_bubbles1==true) && (cmd_ptr->output_detect_bubbles1[0]=='\0') )
    {

      char tmp[]="If you specify --detect_bubbles1, you must also specify an output file with --output_bubbles1";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }
  if ( (cmd_ptr->detect_bubbles2==true) && (cmd_ptr->output_detect_bubbles2[0]=='\0') )
    {

      char tmp[]="If you specify --detect_bubbles2, you must also specify an output file with --output_bubbles2";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }


  
  if ( (cmd_ptr->input_seq==true) && (cmd_ptr->max_read_length==0) )
    {
      char tmp[]="Must specify max read-length if inputting fasta/fastq data\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }
  

  //if making path divergence calls, must have >=2 colours, and must specify the ref colour, and give a ref_fasta list
  if (cmd_ptr->make_pd_calls==true)
    {
      if ((cmd_ptr->using_ref==false)||(cmd_ptr->ref_colour==-1) || (strlen(cmd_ptr->ref_chrom_fasta_list)==0)  )
	{
	  char tmp[]="If --path_divergence_caller is specified, then must also specify --ref_colour and --list_ref_fasta";

	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      printf("coding error - this string is too long:\n%s\n", tmp);
	      exit(1);
	    }
	  strcpy(error_string, tmp);
	  return -1;
	}
      
      if (cmd_ptr->path_divergence_caller_output_stub[0]=='\0')
	{
	  char tmp[]="If --path_divergence_caller is specified, then must also specify --path_divergence_caller_output";
          if (strlen(tmp)>LEN_ERROR_STRING)
            {
              printf("coding error - this string is too long:\n%s\n", tmp);
              exit(1);
            }
          strcpy(error_string, tmp);
          return -1;
	}      
    }

  if ((cmd_ptr->path_divergence_caller_output_stub[0]!='\0') && (cmd_ptr->make_pd_calls==false) )
    {
      char tmp[]="If --path_divergence_caller_output is specified, then must also specify --path_divergence_caller";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  // check multicolour binary  

  if (cmd_ptr->input_multicol_bin==true)
    {
      int num_m_cols;
      FILE* fp = fopen(cmd_ptr->multicolour_bin, "r");
      if (fp==NULL)
	{
	  char tmp[LEN_ERROR_STRING];
	  sprintf(tmp,"Unable to open multicolour bin %s for initial sanity checks\n", cmd_ptr->multicolour_bin) ;
	  strcpy(error_string, tmp);
          return -1;
	}

      int mean_readlens[NUMBER_OF_COLOURS];
      long long total_seqs[NUMBER_OF_COLOURS];
      int* mean_readlens_ptrs[NUMBER_OF_COLOURS];
      long long* total_seqs_ptrs[NUMBER_OF_COLOURS];
      int i;
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  mean_readlens[i]=0;
	  mean_readlens_ptrs[i]=&(mean_readlens[i]);
	  total_seqs[i]=0;
	  total_seqs_ptrs[i]=&(total_seqs[i]);
	}
      boolean is_multicol_bin_ok = check_binary_signature(fp, cmd_ptr->kmer_size, BINVERSION, &num_m_cols, mean_readlens_ptrs, total_seqs_ptrs);

      
      fclose(fp);

      char tmp[LEN_ERROR_STRING];
      if (is_multicol_bin_ok==false)
	{
	  sprintf(tmp,"This binary %s is not compatible with the current de Bruijn graph parameters\n", cmd_ptr->multicolour_bin);
          strcpy(error_string, tmp);
          return -1;
	  
	}
      
      if (num_m_cols<=0)
	{
	  sprintf(tmp,"Corrupt binary %s - signatire claims to have <=0 colours within\n", cmd_ptr->multicolour_bin);
          strcpy(error_string, tmp);
          return -1;

	}
      if (num_m_cols>NUMBER_OF_COLOURS)
	{
	  sprintf(tmp,"Multicolour binary %s contains %d colours, but cortex_var is compiled to support a maximum of %d colours\n", 
		  cmd_ptr->multicolour_bin, num_m_cols, NUMBER_OF_COLOURS);
          strcpy(error_string, tmp);
          return -1;

	}
      else if ( (cmd_ptr->successively_dump_cleaned_colours==false) && (num_m_cols+cmd_ptr->num_colours_in_input_colour_list > NUMBER_OF_COLOURS) )
	{
	  sprintf(tmp,"Between %s (containing %d colours) and %s (containing %d colours), you have exceeded the compile-time limit on colours, %d\n",
		 cmd_ptr->multicolour_bin, num_m_cols, cmd_ptr->colour_list, cmd_ptr->num_colours_in_input_colour_list, NUMBER_OF_COLOURS);
          strcpy(error_string, tmp);
          return -1;
	}

      
      cmd_ptr->num_colours_in_multicol_bin=num_m_cols;
    }



  
  if (cmd_ptr->align_given_list==true)
    {
      char tmp[LEN_ERROR_STRING];

      if (cmd_ptr->format_of_files_to_align==UNSPECIFIED)
	{
	  sprintf(tmp, "If --align is specified, then --align_input_format must also be specified\n");
	  strcpy(error_string, tmp);
          return -1;
	}

      if (cmd_ptr->max_read_length==0)
	{
	  sprintf(tmp, "If --align is specified, then --max_read_len must also be specified\n");
          strcpy(error_string, tmp);
	  return -1;
	}
    }


  if (cmd_ptr->print_colour_overlap_matrix==true)
    {
      int n;

      for (n=0; n<cmd_ptr->num_colours_in_colour_overlap_first_list; n++)
	{
	  if (cmd_ptr->colour_overlap_first_colour_list[n]>NUMBER_OF_COLOURS-1)
	    {
	      char tmp[LEN_ERROR_STRING];
	      sprintf(tmp, "One of the colours you listed in --colour_overlaps is greater than the maximum allowed number of colours for this executable - you specified this at compile time. Either adjust your command-line (typo?) or recompile for more colours\n");
	      strcpy(error_string, tmp);
	      return -1;
	    }
	}

      for (n=0; n<cmd_ptr->num_colours_in_colour_overlap_second_list; n++)
	{
	  if (cmd_ptr->colour_overlap_second_colour_list[n]>NUMBER_OF_COLOURS-1)
	    {
	      char tmp[LEN_ERROR_STRING];
	      sprintf(tmp, "One of the colours you listed in --colour_overlaps is greater than the maximum allowed number of colours for this executable - you specified this at compile time. Either adjust your command-line (typo?) or recompile for more colours\n");
	      strcpy(error_string, tmp);
	      return -1;
	    }
	}


    }


  return 0;
}


CmdLine parse_cmdline( int argc, char* argv[], int unit_size) 
{	
  int i;
  printf("Command: ");
  for(i=0;i<argc;i++){
    printf("%s ",argv[i]);
  }
  printf("\n");
  // printf("Unit size:%i\n",unit_size);

  CmdLine cmd_line;
  default_opts(&cmd_line);

  char error_string[LEN_ERROR_STRING]="";
  int err = parse_cmdline_inner_loop(argc, argv, unit_size, &cmd_line, error_string);


  if (err==-1) 
    {
      printf("Error in cmd-line input: %s\n", error_string);
      exit(1);
    }


  char error_string2[LEN_ERROR_STRING]="";
    
  err = check_cmdline(&cmd_line, error_string2);

  if (err == -1)
    {
      printf("Error in cmd-line input: %s\n", error_string2);
      exit(1);
    }
  
  return cmd_line;
  
}








//RELIES on the list (arg1) ending in '\0'.
    //returns -1 on error
int get_numbers_from_comma_sep_list(char* list, int* return_list, int max_len_return_list)
{
  int number[max_len_return_list];

  int i;
  int current=0;
  number[current]=0;


  for(i=0;i<strlen(list);i++){
    if (list[i]==',')
      {
	current++;
	number[current]=0;
	if (current>=max_len_return_list)
	  {
	    printf("Too many colours in list\n");
	    return -1;
	  }
      }
    else{
      if ( (list[i]>='0' && list[i]<='9') ||  (i+1 == strlen(list))   ){
	number[current]=number[current]*10 + (list[i]-'0');
      }
      else{
	printf("This part of cmd-line argument is badly formatted: %s\n", list);
	return -1;
      }
    }
  }

  
  for(i=0;i<=current;i++){
    return_list[i]=number[i];
  }
  return current+1;
}






  //allows user to enter something like 1,2,3/562,87 for --detect_bubbles, and will then detect bubbles where one branch is in union of 1,2,3 and other in union of 562,87,
  //or vice-versa.
  //returns -1 on error, 0 otherwise. Parses argument2 and gets the colours into arg1
int parse_colourinfo_argument(CmdLine* cmd, char* arg, int len_arg, char* text_for_error_describing_which_option_this_is, int which_detect_bubbles)
{

  if ( (which_detect_bubbles!=1) && (which_detect_bubbles!=2) && (which_detect_bubbles!=3) )
    {
      printf("Called parse_colourinfo_argument with last argument %d - should be 1 or 2\n", which_detect_bubbles);
      return -1;
    }

  
  if (len_arg<MAX_LEN_DETECT_BUB_COLOURINFO)
	  {
	    int k;
	    int num_slashes=0;
	    int posn_slash=0;
	    for (k=0; k<len_arg; k++)
	      {
		if (arg[k]=='/') 
		  {
		    num_slashes++;
		    posn_slash=k; //will only keep position of "last" slash - will check that there is in fact only 1
		  }
	      }

	    if (num_slashes != 1)
	      {
		printf("%s option requires two sets of comma-separated integers, separated by a forward-slash. eg 1/1, or 1,2,3/3,4,5", 
		     text_for_error_describing_which_option_this_is);		
		return -1;
	      }


	    char left[MAX_LEN_DETECT_BUB_COLOURINFO];
	    //left[0]='\0';
	    set_string_to_null(left, MAX_LEN_DETECT_BUB_COLOURINFO);
	    char right[MAX_LEN_DETECT_BUB_COLOURINFO];
	    //right[0]='\0';
	    set_string_to_null(right, MAX_LEN_DETECT_BUB_COLOURINFO);

	    strncpy(left, arg,posn_slash);
	    strncpy(right, arg+posn_slash+1, len_arg-posn_slash-1);
	    //printf("Split %s into %s and %s", arg, left, right);
	    
	    int left_colours[MAX_COLOURS_ALLOWED_TO_MERGE];
	    int right_colours[MAX_COLOURS_ALLOWED_TO_MERGE];
	    int num_left_colours = get_numbers_from_comma_sep_list(left, left_colours, MAX_COLOURS_ALLOWED_TO_MERGE);
	    int num_right_colours = get_numbers_from_comma_sep_list(right, right_colours, MAX_COLOURS_ALLOWED_TO_MERGE);


	    if ( (num_left_colours==-1) || (num_right_colours==-1) )
	      {
		printf("Num of left or right colours is -1. They are %d and %d\n", num_left_colours, num_right_colours);
		return -1; //error.
	      }
	    //printf("Num left colours is %d and num right is %d\n", num_left_colours, num_right_colours);
	    if (which_detect_bubbles==1)
	      {
		cmd->num_colours_in_detect_bubbles1_first_colour_list=num_left_colours;
		cmd->num_colours_in_detect_bubbles1_second_colour_list=num_right_colours;
		copy_list(left_colours, num_left_colours, cmd->detect_bubbles1_first_colour_list, cmd->num_colours_in_detect_bubbles1_first_colour_list);
		copy_list(right_colours, num_right_colours, cmd->detect_bubbles1_second_colour_list, cmd->num_colours_in_detect_bubbles1_second_colour_list);
	      }
	    else if (which_detect_bubbles==2)
	      {
		cmd->num_colours_in_detect_bubbles2_first_colour_list=num_left_colours;
		cmd->num_colours_in_detect_bubbles2_second_colour_list=num_right_colours;
		copy_list(left_colours, num_left_colours, cmd->detect_bubbles2_first_colour_list, cmd->num_colours_in_detect_bubbles2_first_colour_list);
		copy_list(right_colours, num_right_colours, cmd->detect_bubbles2_second_colour_list, cmd->num_colours_in_detect_bubbles2_second_colour_list);
	      }
	    else if (which_detect_bubbles==3)
	      {
		cmd->num_colours_in_colour_overlap_first_list=num_left_colours;
		cmd->num_colours_in_colour_overlap_second_list=num_right_colours;
		copy_list(left_colours, num_left_colours, cmd->colour_overlap_first_colour_list, cmd->num_colours_in_colour_overlap_first_list);
		copy_list(right_colours, num_right_colours, cmd->colour_overlap_second_colour_list, cmd->num_colours_in_colour_overlap_second_list);
	      }
	    else
	      {
		printf("Unexpected argument to parse_colourinfo_argument, which_detect_bubbles is not 1,2 or 3");
		exit(1);
	      }
	    
	  }
  else
    {
      printf("Coding error. Have colour info argument longer than %d, which should not happen - but should have been caught before now\n", MAX_LEN_DETECT_BUB_COLOURINFO);
      exit(1);
    }

  return 0;
}





  //returns -1 on error, 0 otherwise. Parses argument2 and gets the colours into arg1
int parse_commasep_list(CmdLine* cmd, char* arg, int len_arg, char* text_for_error_describing_which_option_this_is)
{

  if (len_arg<MAX_LEN_DETECT_BUB_COLOURINFO)
	  {
	    int k;

	    for (k=0; k<len_arg; k++)
	      {
		if (arg[k]=='/') 
		  {
		    printf("%s option requires just one comma-separated list of integers, with no forward-slash /. Perhaps you are confusing with --bubbles.", 
		     text_for_error_describing_which_option_this_is);		
		    return -1;
		  }
	      }

	    char list[MAX_LEN_DETECT_BUB_COLOURINFO];
	    set_string_to_null(list, MAX_LEN_DETECT_BUB_COLOURINFO);
	    strncpy(list, arg, strlen(arg));
	    
	    int list_colours[MAX_COLOURS_ALLOWED_TO_MERGE];
	    int num_list_colours = get_numbers_from_comma_sep_list(list, list_colours, MAX_COLOURS_ALLOWED_TO_MERGE);

	    if (num_list_colours==-1) 
	      {
		return -1; //error.
	      }

	    cmd->num_colours_in_pd_colour_list=num_list_colours;
	    copy_list(list_colours, num_list_colours, cmd->pd_colour_list, cmd->num_colours_in_pd_colour_list);
	    
	  }
  else
    {
      printf("Coding error. Have colour info argument longer than %d, which should not happen - but should have been caught before now\n", MAX_LEN_DETECT_BUB_COLOURINFO);
      exit(1);
    }

  return 0;
}


