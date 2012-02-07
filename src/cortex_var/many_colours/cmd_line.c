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


int isNumeric (const char * s)
{
  if (s == NULL || *s == '\0')
    return 0;
  char * p;
  strtod (s, &p);
  return *p == '\0';
}


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
  int binversion_in_header;
  check_binary_signature(fp, kmer_size, BINVERSION, &num_cols, mean_readlens_ptrs, total_seqs_ptrs, &binversion_in_header);
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


"   [--help] \t\t\t\t\t\t\t=\t This help screen.\n\n" \
" \n ** DATA LOADING ** \n\n"\
  // -v
"   [--format TYPE] \t\t\t\t\t\t=\t File format for input in se_list and pe_list. All files assumed to be of the same format.\n\t\t\t\t\t\t\t\t\t Type must be FASTQ or FASTA\n" \
"   [--colour_list FILENAME] \t\t\t\t\t=\t File of filenames, one per colour. n-th file is a list of\n\t\t\t\t\t\t\t\t\t single-colour binaries to be loaded into colour n.\n\t\t\t\t\t\t\t\t\t Cannot be used with --se_list or --pe_list \n" \
"   [--multicolour_bin FILENAME] \t\t\t\t=\t Filename of a multicolour binary, will be loaded first, into colours 0..n.\n\t\t\t\t\t\t\t\t\t If using --colour_list also, those will be loaded into subsequent colours, after this.\n" \
"   [--se_list FILENAME] \t\t\t\t\t=\t List of single-end fasta/q to be loaded into a single-colour graph.\n\t\t\t\t\t\t\t\t\t Cannot be used with --colour_list\n" \
"   [--pe_list FILENAME] \t\t\t\t\t=\t Two filenames, comma-separated: each is a list of paired-end fasta/q to be \n\t\t\t\t\t\t\t\t\t loaded into a single-colour graph. Lists are assumed to ordered so that \n\t\t\t\t\t\t\t\t\t corresponding paired-end fasta/q files are at the same positions in their lists.\n\t\t\t\t\t\t\t\t\t Currently Cortex only use paired-end information to remove\n\t\t\t\t\t\t\t\t\t PCR duplicate reads (if that flag is set).\n\t\t\t\t\t\t\t\t\t Cannot be used with --colour_list\n" \
"   [--kmer_size INT] \t\t\t\t\t\t=\t Kmer size. Must be an odd number.\n" \
"   [--mem_width INT] \t\t\t\t\t\t=\t Size of hash table buckets (default 100).\n" \
  //-g 
"   [--mem_height INT] \t\t\t\t\t\t=\t Number of buckets in hash table in bits (default 10). \n\t\t\t\t\t\t\t\t\t Actual number of buckets will be 2^(the number you enter)\n" \
  // -n 
"   [--fastq_offset INT] \t\t\t\t\t=\t Default 33, for standard fastq.\n\t\t\t\t\t\t\t\t\t Some fastq directly from different versions of Illumina machines require different offsets.\n" \
  // -p
"   [--dump_binary FILENAME] \t\t\t\t\t=\t Dump a binary file, with this name (after applying error-cleaning, if specified).\n" \
  // -w
"   [--max_read_len] \t\t\t\t\t\t=\t For fastq, this is the Maximum read length over all input files.\n\t\t\t\t\t\t\t\t\t For fasta it is the size of chunk in which the reads are read (if reading a whole chromosome\n\t\t\t\t\t\t\t\t\t a typical value to use is 10000. Values above 20000 forbidden.\n\t\t\t\t\t\t\t\t\t (Mandatory if fastq or fasta files are input.)\n" \
" \n**** FILTERING AND ERROR CLEANING OPTIONS ****\n\n"\
  // -m
"   [--quality_score_threshold INT] \t\t\t\t=\t Filter for quality scores in the input file (default 0).\n" \
  // -j
"   [--remove_pcr_duplicates] \t\t\t\t\t=\t Removes PCR duplicate reads by ignoring read pairs if both \n\t\t\t\t\t\t\t\t\t reads start at the same k-mer as a previous read,\n\t\t\t\t\t\t\t\t\t and single-ended reads if they start at the same k-mer as a previous read\n" \
  // -k
"   [--cut_homopolymers INT] \t\t\t\t\t=\t Breaks reads at homopolymers of length >= this threshold.\n\t\t\t\t\t\t\t\t\t (i.e. max homopolymer in filtered read==threshold-1, and New read starts after homopolymer)\n" \
  // -O
"   [--remove_low_coverage_supernodes INT]\t\t\t\t=\t Remove all supernodes where max coverage is <= the limit you set. Overrides --remove_seq_errors. Recommended method.\n" \
   // -o
"   [--remove_seq_errors] \t\t\t\t\t=\t Remove tips + remove supernodes with coverage everywhere==1\n" \
  // -B
"   [--remove_low_coverage_kmers INT] \t\t\t\t=\t Filter for kmers with coverage less than or equal to  threshold. Not recommended. See manual and our paper for why\n"  \
  // -E
"   [--load_colours_only_where_overlap_clean_colour INT] \t=\t Only load nodes from binary files in the colour-list when they overlap a\n\t\t\t\t\t\t\t\t\t specific colour (e.g. that contains a cleaned pooled graph);\n\t\t\t\t\t\t\t\t\t requires you to specify this particular colour. You must have loaded that colour beforehand, using --multicolour_bin\n"  \
  // -F
"   [--successively_dump_cleaned_colours SUFFIX] \t\t=\t Only to be used when also using --load_colours_only_where_overlap_clean_colour and --multicolour_bin\n\t\t\t\t\t\t\t\t\t Used to allow error-correction of low-coverage data on large numbers of individuals with large genomes.\n\t\t\t\t\t\t\t\t\t Requires the user specify a suffix which will be added to the names of cleaned binaries. See manual for details.\n"  \
" \n\n**** OUTPUT STATISTICS AND SUPERNODES****\n\n" \
  // -A
"   [--dump_covg_distribution FILENAME] \t\t\t\t=\t Print k-mer coverage distribution to the file specified\n" \

  // -D
"   [--dump_filtered_readlen_distribution FILENAME] \t\t=\t Dump to file the distribution of \"effective\" read lengths after quality/homopolymer/PCR dup filters \n"  \
  // -q
"   [--output_supernodes FILENAME] \t\t\t\t=\t Dump a fasta file of all the supernodes (after applying all specified actions on graph).\n" \
"\n\n** VARIANT CALLING **\n\n"\
  // -y
"   [--max_var_len INT] \t\t\t\t\t\t=\t Maximum variant size searched for. Default 10kb. \n" \
  // -r
"   [--detect_bubbles1 COMMA_SEP_COLOURS/COMMA_SEP_COLOURS] \t=\t Find all the bubbles in the graph where the two branches lie in the specified colours\n\t\t\t\t\t\t\t\t\t (after applying all specified actions on graph).\n\t\t\t\t\t\t\t\t\t Typical use would be --detect_bubbles1 1/1 to find hets in colour 1,\n\t\t\t\t\t\t\t\t\t or --detect_bubbles1 0/1 to find homozygous non-reference bubbles where one branch is in colour 0 (and not colour1)\n\t\t\t\t\t\t\t\t\t and the other branch is in colour1 (but not colour 0).\n\t\t\t\t\t\t\t\t\t However, one can do more complex things:\n\t\t\t\t\t\t\t\t\t e.g.  --detect_bubbles1 1,2,3/4,5,6 to find bubbles where one branch is in 1,2 or 3 (and not 4,5 or 6)\n\t\t\t\t\t\t\t\t\t and the other branch in colour 4,5 or 6 (but not 1,2, or 3).\n\t\t\t\t\t\t\t\t\t See the manual for a detailed description of the syntax. It is possible to\n\t\t\t\t\t\t\t\t\t specify \"all colour\" rather than enumerate them eplicitly, and also to exclude colours\n" \
  // -s
"   [--output_bubbles1 FILENAME]\t\t\t\t\t=\t Bubbles called in detect_bubbles1 are dumped to this file.\n" \
  // -x
"   [--print_colour_coverages]\t\t\t\t\t=\t Print coverages in all colours for supernodes and variants.\n" \
  // -M
"   [--exclude_ref_bubbles]\t\t\t\t\t=\t If you have specified --ref_colour, this will exclude any bubble in that colour from being called by the Bubble Caller.\n" \
  // -t
"   [--detect_bubbles2 COMMA_SEP_COLOURS/COMMA_SEP_COLOURS] \t=\t Exactly the same as detect_bubbles1, but allows you to make\n\t\t\t\t\t\t\t\t\t a second set of bubble calls immediately afterwards.\n\t\t\t\t\t\t\t\t\t This is to accomodate the common use-case where one loads a reference\n\t\t\t\t\t\t\t\t\t and an individual, and then wants to call homs, and hets.\n" \
  // -u
"   [--output_bubbles2 FILENAME]\t\t\t\t\t=\t Bubbles called in detect_bubbles2 are dumped to this file.\n" \

  // -l
"   [--path_divergence_caller [COMMA_SEP_COLOURS|SQUARE_OPEN_BRACKET_PRECEDED_AND_SEP_COLOURS]] \t\t=\t Make Path Divergence variant calls.\n\t\t\t\t\t\t\t\t\t Must specify colour of sample in which you want to find\n\t\t\t\t\t\t\t\t\t variants compared with the reference.\n\t\t\t\t\t\t\t\t\t This sample colour can be a union of colours (comma-separated list) \n Or, given a square open bracket [ PRECEDED AND SEPARATED list(example [2[3[10 ) the caller will call against each colour in turn\n\t\t\t\t\t\t\t\t\t Must also specify --ref_colour and --list_ref_fasta\n" \
  // -I
"   [--path_divergence_caller_output PATH_STUB]\t\t\t=\t Specifies the path and beginning of filenames of Path Divergence caller output files.\n\t\t\t\t\t\t\t\t\t One output file will be created per  reference fasta listed in --list_ref_fasta\n" \
  // -i
"   [--ref_colour INT] \t\t\t\t\t\t=\t Colour of reference genome.\n" \
 // -z
"   [--list_ref_fasta FILENAME] \t\t\t\t\t=\t File listing reference chromosome fasta file(s); needed for path-divergence calls. \n" \
"\n\n **** ADVANCED OPTIONS **** \n\n"\
  // -P  
"   [--experiment_type]\t\t\t\t\t\t=\t The statistical models for determining genotype likelihoods, and for deciding if bubbles are repeat or variants,\n\t\t\t\t\t\t\t\t\t require knowledge of whether each sample is a separate diploid/haploid individual. \n\t\t\t\t\t\t\t\t\t Enter type of experiment (EachColourADiploidSample, EachColourADiploidSampleExceptTheRefColour, \n\t\t\t\t\t\t\t\t\t EachColourAHaploidSample,EachColourAHaploidSampleExceptTheRefColour). \n\t\t\t\t\t\t\t\t\t This is only needed for determining likelihoods, so ignore this is you are pooling samples within a colour (support to be added for this later).\n" \
  // Q
"   [--estimated_error_rate]\t\t\t\t\t=\t If you have some idea of the sequencing error rate (per base-pair), enter it here. eg 0.01. Currently used in calculating likelihoods\n"\
  // -L
"   [--genome_size]\t\t\t\t\t\t=\t If you specify --experiment_type, and therefore want to calculate likelihoods, you must also specify the (estimated) genome size in bp.\n" \

  // -G
"   [--align FILENAME,{output binary name|no}] \t\t\t\t\t\t=\t Aligns a list of fasta/q files to the graph, and prints coverage of each kmer in each read in each colour.\n\t\t\t\t\t\t\t\t\t Takes two arguments. First, a list of fasta/q. Second, either an output filename (if you want it to dump a binary of the part of the graph touched by the alignment) OR just \"no\" \n\t\t\t\t\t\t\t\t\t Must also specify --align_input_format, and --max_read_len\n"  \
  // -H
"   [--align_input_format TYPE] \t\t\t\t\t=\t --align requires a list of fasta or fastq. This option specifies the input format as LIST_OF_FASTQ or LIST_OF_FASTA\n"  \

  // T
"   [--estimate_genome_complexity FILENAME]\t\t\t=\t Print estimated genome complexity by reading up to 10,000 reads from the given file\n"\
"\n\n **** OPTIONS THAT I MAY REMOVE - NOT CLEAR ANYONE HAS A GOOD USE FOR THEM - EMAIL ME IF YOU LIKE IT\n\n"\
  // -J
"   [--colour_overlaps COMMA_SEP_COLOURS/COMMA_SEP_COLOURS]\t=\t Compares each coloured subgraph in the first list with all of the \n\t\t\t\t\t\t\t\t\t coloured subgraphs in the second list. Outputs a matrix to stdout;\n\t\t\t\t\t\t\t\t\t (i,j)-element is the number of nodes in both colour-i (on first list)\n\t\t\t\t\t\t\t\t\t and colour-j (on second list).\n" \

"\n\n**** EARLY ACCESS/BETA OPTIONS **** \n\n"\
  // R
"   [--genotype_site]\t\t\t\t\t\t=\t (Beta code!) Genotype a single (typically multiallelic) site. Syntax is slightly complex. \n\t\t\t\t\t\t\t\t\t Requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN.\n\t\t\t\t\t\t\t\t\t x,y is a comma-sep list of colours to genotype.\n\t\t\t\t\t\t\t\t\t z is the reference-minus-site colour.\n\t\t\t\t\t\t\t\t\t N is the number of alleles for this site (which cortex assumes are loaded\n\t\t\t\t\t\t\t\t\t in a multicolour_bin containing those alleles first, one per colour).\n\t\t\t\t\t\t\t\t\t Cortex will genotype combinations A through B of the N choose 2 possible genotypes (allows parallelisation);\n\t\t\t\t\t\t\t\t\t fasta is the file listing one read per allele.\n\t\t\t\t\t\t\t\t\t CLEANED or UNCLEANED allows Cortex to tailor its genotyping model.\n\t\t\t\t\t\t\t\t\t p,q are two free/unused colours that Cortex will use internally. \n\t\t\t\t\t\t\t\t\t yes/no specifies whether to use the more sophisticated error model, which is still in development. \n\t\t\t\t\t\t\t\t\t I recommend you stick with \"no\" for now.\n\t\t\t\t\t\t\t\t\t The final argument, MIN, is optional and allows performance speedup\n\t\t\t\t\t\t\t\t\t by disarding any genotype with log-likelihood<MIN.\n\t\t\t\t\t\t\t\t\t See manual for details.Must also specify --max_var_len to give the length of the longest allele\n"\
  // -V
"   [--print_novel_contigs]\t\t\t\t\t=\t Allows printing of novel sequence absent from a reference (or more generally, absent from a set of colours)\n\t\t\t\t\t\t\t\t\t Takes arguments in this format a,b,../c,d,../x/y/<output filename>\n\t\t\t\t\t\t\t\t\t Cortex will find supernodes in the union graph of colours a,b,..\n\t\t\t\t\t\t\t\t\t Typically the list c,d,.. of colours is just one colour  - that of the reference.\n\t\t\t\t\t\t\t\t\t Cortex will print contigs (supernodes) to the output file which satisfy the following criteria\n\t\t\t\t\t\t\t\t\t Contigs must be at least x base pairs long\n\t\t\t\t\t\t\t\t\t The percentage (as integer) of kmers in the contig which are present in ANY of the colours c,d,... must be at most 1-y. \n\t\t\t\t\t\t\t\t\t i.e. y is the minimum proportion of novel kmers in a contig. Typically this is 100.\n\t\t\t\t\t\t\t\t\t We ignore the first and last kmer of the contig, as these will typically connect to the reference\n  " \


  // -K
  //"   [--require_hw]\t\t\t\t\t\t=\t For each bubble found, calculate likelihood of observed coverage \n\t\t\t\t\t\t\t\t\t under 3 models (repeat, error, variation obeying Hardy-Weinberg)\n\t\t\t\t\t\t\t\t\t Only call variants where the bubble is more likely (according to these models) to be a variant.\n" \

;






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
  //novelsseq stuff
  c->print_novel_contigs=false;
  c->novelseq_contig_min_len_bp=100;
  c->novelseq_min_percentage_novel=100;
  initialise_int_list(c->novelseq_colours_search, MAX_COLOURS_ALLOWED_TO_MERGE);
  c->numcols_novelseq_colours_search=0;
  initialise_int_list(c->novelseq_colours_avoid, MAX_COLOURS_ALLOWED_TO_MERGE);
  c->numcols_novelseq_colours_avoid=0;
  set_string_to_null(c->novelseq_outfile, MAX_FILENAME_LEN);

  c->working_colour1 = -1;
  c->working_colour2 = -1;
  // c->working_colour3_for_1net=-1;
  //c->working_colour4_for_2net=-1;
  c->manually_entered_seq_error_rate=-1;
  c->manually_override_error_rate=false;
  c->expt_type = Unspecified;
  c->genome_size=0;
  c->kmer_size = -1;
  c->bucket_size = 100;
  c->number_of_buckets_bits = 10;
  c->ref_colour=-1;//there are places where I specifically check to see if this is -1, and if so, assume there is no reference
  c->homopolymer_limit=-1;
  c->quality_score_threshold=0;
  c->node_coverage_threshold=0;
  c->quality_score_offset = 33;//standard fastq, not illumina v-whatever fastq  
  c->max_read_length = 0;
  c->max_var_len = 10000;
  c->specified_max_var_len = false;
  c->remv_low_covg_sups_threshold=-1;
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
  set_string_to_null(c->fasta_alleles_for_complex_genotyping, MAX_FILENAME_LEN);
  set_string_to_null(c->filelist_1net_binaries_for_alleles, MAX_FILENAME_LEN);
  set_string_to_null(c->fastaq_for_estimating_genome_complexity, MAX_FILENAME_LEN);
  // set_string_to_null(c->filelist_2net_binaries_for_alleles, MAX_FILENAME_LEN);

  //booleans
  c->exclude_ref_bubbles=false;
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
  c->pd_calls_against_each_listed_colour_consecutively=false;
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
  c->dump_aligned_overlap_binary=false;
  c->print_colour_overlap_matrix=false;
  c->format_of_files_to_align=UNSPECIFIED;
  c->apply_model_selection_at_bubbles=false;
  c->estimate_genome_complexity=false;

  c->num_colours_to_genotype=0;
  initialise_int_list(c->list_colours_to_genotype, NUMBER_OF_COLOURS);
  c->colour_of_reference_with_site_excised=-1;
  c->num_alleles_of_site=0;
  c->first_genotype_to_calc_likelihoods_for=0;
  c->last_genotype_to_calc_likelihoods_for=0;
  c->genotype_complex_site=false;
  c->using_1net=false;
  c->using_2net=false;
  c->assump_for_genotyping=AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice;
  c->min_acceptable_llk=-999999999;
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
    {"require_hw",no_argument,NULL,'K'},
    {"genome_size",required_argument,NULL,'L'},
    {"exclude_ref_bubbles",no_argument,NULL,'M'},
    {"remove_low_coverage_supernodes",required_argument,NULL,'O'},
    {"experiment_type", required_argument, NULL, 'P'},
    {"estimated_error_rate", required_argument, NULL, 'Q'},
    {"genotype_site", required_argument, NULL, 'R'},
    //  {"detect_alleles", required_argument, NULL, 'S'},
    {"estimate_genome_complexity", required_argument, NULL, 'T'},
    {"print_novel_contigs", required_argument, NULL, 'V'},
    {0,0,0,0}	
  };
  

  //do not change this! Only really matters for testing, but getopt_long uses 
  //variables which are not local to this function, so when we run tests, and call this function
  // repeatedly, those variable values are carried across. The following line resets this.
  optind=1;
  
 
  opt = getopt_long(argc, argv, "ha:b:c:d:e:f:g:i:jk:l:m:n:op:q:r:s:t:u:v:w:xy:z:A:B:C:D:E:F:G:H:I:J:KL:MO:P:Q:R:T:V:", long_options, &longopt_index);

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
	  errx(1,"[--path_divergence_caller] option requires at least one colour. If >1, must EITHER be comma separated (e.g. 1,2,3) OR [- preceded and separated  - eg [1[4[76.\n");

	//we need to decide whether we call variants ONCE against the union of the listed colours
	// or we call repeatedly, once against each colour listed (but genotyping all colours obviously)
	boolean commasep_so_call_against_union=true;
	if (optarg[0]=='[')
	  {
	    commasep_so_call_against_union=false;
	    cmdline_ptr->pd_calls_against_each_listed_colour_consecutively=true;
	  }
	parse_commasep_or_open_square_brack_sep_list(cmdline_ptr, //->num_colours_in_pd_colour_list, cmdline_ptr->pd_colour_list, 
					    optarg, strlen(optarg), "[--path_divergence_caller]", commasep_so_call_against_union);

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
	cmdline_ptr->remv_low_covg_sups_threshold = 1;
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
	    errx(1,"[--max_read_len] option requires (positive) integer argument. Either you have entered 0 or such an enormous number i has overflowed\n");
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
	cmdline_ptr->specified_max_var_len = true;
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

    case 'G'://align a list of fasta/q to the graph and print colour coverages, and optionally dump binary of subgraph which is touched by alignment
      {
	if (optarg==NULL)
	  errx(1,"[--align FILENAME,[output binary name|no]] option requires a filename, a comma, then either an output binary name or \"no\", depending on whether you want to dump a single-colour binary of the part of the graph that is aligned to.\n");

	if (strlen(optarg)<MAX_FILENAME_LEN*2)
	  {
	
	    char* filename = NULL;
	    char delims[] = ",";
	    char temp[MAX_FILENAME_LEN*2];
	    temp[0]='\0';
	    strcpy(temp, optarg);
	    filename = strtok(temp, delims );
	    if (filename==NULL)
	      {
		errx(1,"[--align] option requires a filename, a comma, then either an output binaryname or \"no\", depending on whether you want to dump a binary of part ofthe graph that is aligned to - i.e. the overlap of your sequences and the graph\n");
	      }
	    else if (access(filename,F_OK)!=0){
	      errx(1,"[--align]   - Cannot open this filelist of files to align  [%s] . Abort.\n", filename);
	    }

	    strcpy(cmdline_ptr->list_fastaq_to_align,filename);
	    cmdline_ptr->align_given_list=true;

	    char* dump_bin=NULL;
	    dump_bin = strtok(NULL, delims);
	    if (dump_bin==NULL)
	      {
		//do nothing - default is dump no binary
		printf("No output file specified, so will not dump binary of overlap of alignment and graph\n");
	      }
	    else if (strcmp(dump_bin, "no")==0)
	      {
		//do nothing - default is dump no binary
		printf("User has specified NOT to dump binary of overlap of alignment and graph\n");

	      }
	    else if (access(dump_bin,F_OK)==0){
	      errx(1,"[--align] output binary name [%s] already exists - will not overwrite. Abort.\n",dump_bin);
	    }
	    else
	      {
		cmdline_ptr->dump_aligned_overlap_binary=true;
		strcpy(cmdline_ptr->output_aligned_overlap_binname, dump_bin);
		printf("Will dump binary of overlap of sequences with graph in file %s\n", cmdline_ptr->output_aligned_overlap_binname);
	      }

	  }
	else
	  {
	    errx(1,"[--align] filename(s) too long [%s]",optarg);
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
    case 'K':
	{
	  //cmdline_ptr->apply_model_selection_at_bubbles=true;
	  break;
	}
    case 'L':
	{
	  if (optarg==NULL)
	    errx(1,"[--genome_size] option requires an integer argument, the genome size in bp");
	  cmdline_ptr->genome_size = atoll(optarg);

	  break;
	}
    case 'M':
      {
	cmdline_ptr->exclude_ref_bubbles=true;
	break;
      }
    case 'O':
	{
	  if (optarg==NULL)
	    errx(1,"[--remove_low_coverage_supernodes] option requires an integer argument; supernodes with covg<= this limit will be cleaned/removed");

	  cmdline_ptr->remv_low_covg_sups_threshold = atoi(optarg);

	  break;
	}
    case 'P':
	{
	  if (optarg==NULL)
	    errx(1,"[--experiment_type] option requires a string argument - the type of experiment. Acceptable options are EachColourADiploidSample, EachColourADiploidSampleExceptTheRefColour, EachColourAHaploidSample, EachColourAHaploidSampleExceptTheRefColour.");

	    if (strcmp(optarg, "EachColourADiploidSample")==0)
	      {
		cmdline_ptr->expt_type = EachColourADiploidSample;
	      }
	    else if (strcmp(optarg, "EachColourADiploidSampleExceptTheRefColour")==0)
	      {
		cmdline_ptr->expt_type = EachColourADiploidSampleExceptTheRefColour;
	      } 
	    else if (strcmp(optarg, "EachColourAHaploidSample")==0)
	      {
		cmdline_ptr->expt_type = EachColourAHaploidSample;
	      }
	    else if (strcmp(optarg, "EachColourAHaploidSampleExceptTheRefColour")==0)
	      {
		cmdline_ptr->expt_type = EachColourAHaploidSampleExceptTheRefColour;
	      }
	    else
	      {
		errx(1,"[--experiment_type] option requires a string argument - the type of experiment. Acceptable options are EachColourADiploidSample, EachColourADiploidSampleExceptTheRefColour, EachColourAHaploidSample, EachColourAHaploidSampleExceptTheRefColour. Your entry was not one of these options");
		
	      }
	    break;
	}
    case 'Q':
      {
	if (optarg==NULL)
	  errx(1,"[--estimated_error_rate] option requires an argument, the per base-pair error rate.");
	
	if (isNumeric(optarg))
	  {
	    cmdline_ptr->manually_override_error_rate=true;
	    cmdline_ptr->manually_entered_seq_error_rate = strtod(optarg,NULL);
	  }
	else
	  {
	    errx(1,"[--estimated_error_rate] option requires an argument, the per base-pair error rate. You have entered a non-numeric value");
	  }
	break;
      }

    case 'R':
      {
	if (optarg==NULL)
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q.  x,y is a comma-sep list of colours to genotype. z is the reference-minus-site colour. N is the number of alleles for this site (which cortex assumes are loaded in a multicolour_bin containing precisely and only those alleles, one per colour). Cortex will genotype combinations A through B of the N choose 2 possible genotypes (allows parallelisation); fasta is the file listing one read per allele. CLEANED or UNCLEANED allow Cortex to tailor its genotyping model. p,q are two free/unused colours that Cortex will use internally. See manual for details.");

	parse_genotype_site_argument(optarg, cmdline_ptr->list_colours_to_genotype, &(cmdline_ptr->num_colours_to_genotype),
				     &(cmdline_ptr->colour_of_reference_with_site_excised), &(cmdline_ptr->num_alleles_of_site),
				     &(cmdline_ptr->first_genotype_to_calc_likelihoods_for),
				     &(cmdline_ptr->last_genotype_to_calc_likelihoods_for),
				     cmdline_ptr->fasta_alleles_for_complex_genotyping, &(cmdline_ptr->assump_for_genotyping),
				     &(cmdline_ptr->working_colour1), &(cmdline_ptr->working_colour2),
				     &(cmdline_ptr->using_1net), &(cmdline_ptr->using_2net),
				     cmdline_ptr->filelist_1net_binaries_for_alleles, &(cmdline_ptr->min_acceptable_llk) );

	cmdline_ptr->genotype_complex_site=true;
	break;
      }
    case 'T'://estimate_genome_complexity
      {
	if (optarg==NULL)
	  errx(1,"[--estimate_genome_complexity] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    if (access(optarg,R_OK)==-1)
	      {
		errx(1,"[--estimate_genome_complexity] filename [%s] cannot be accessed",optarg);
	      }
	    
	    cmdline_ptr->estimate_genome_complexity=true;
	    strcpy(cmdline_ptr->fastaq_for_estimating_genome_complexity,optarg);
	  }
	else
	  {
	    errx(1,"[--estimate_genome_complexity] filename too long [%s]",optarg);
	  }
	
	break;
      }
    case 'V'://print_novel_contigs
	{
	  if (optarg==NULL)
	    {
	    errx(1,"[--print_novel_contigs] option requires a two sets of comma-separated numbers, separated by a slash, followed by an integer (minium contig length), followed by a slash, followed by a floating point number (minimum proportion of kmers in contig which are novel), followed by a slash, followed by an output filename. eg --print_novel_contigs 1,2/3/100/0.9/novel.txt will find contigs in union of colours 1 and 2, which are of minimum length 100bp, where 90 percent of the kmers in the contig do not overlap colour 3, and print it to file novel.txt");
	    }
          if (strlen(optarg)>MAX_FILENAME_LEN)
	    {
	      errx(1,"[--print_novel_contigs] argument too long [%s]",optarg);
	    }

	  parse_novelseq_args(optarg, 
			      &(cmdline_ptr->novelseq_colours_search), &(cmdline_ptr->numcols_novelseq_colours_search),
			      &(cmdline_ptr->novelseq_colours_avoid),  &(cmdline_ptr->numcols_novelseq_colours_avoid), 
			      &(cmdline_ptr->novelseq_contig_min_len_bp), &(cmdline_ptr->novelseq_min_percentage_novel),
			      &(cmdline_ptr->novelseq_outfile));
	  cmdline_ptr->print_novel_contigs=true;

	  break;
	}
      

    }
    opt = getopt_long(argc, argv, "ha:b:c:d:e:f:g:i:jk:lm:n:opqr:s:t:u:v:w:xy:z:A:B:C:D:E:F:G:H:I:J:KL:MO:P:Q:R:T:V:", long_options, &longopt_index);
  }   
  
  return 0;
}
  

int check_cmdline(CmdLine* cmd_ptr, char* error_string)
{
  
  if (cmd_ptr->kmer_size==-1)
    {
      char tmp[] = "You must specify kmer_size\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }
  //all colours in detect_bubbles must be less than NUMBER_OF_COLOURS!
  if (cmd_ptr->detect_bubbles1==true)
    {
      int i;
      for (i=0; i<cmd_ptr->num_colours_in_detect_bubbles1_first_colour_list; i++)
	{
	  if (cmd_ptr->detect_bubbles1_first_colour_list[i]>=NUMBER_OF_COLOURS)
	    {
	      char tmp[] = "in detect_bubbles1 Cannot specify a number > %d, which you specified at compile time as the number of colours supported. Recompile/see manual.\n";
	      if (strlen(tmp)>LEN_ERROR_STRING)
		{
		  printf("coding error - this string is too long:\n%s\n", tmp);
		  exit(1);
		}
	      strcpy(error_string, tmp);
	      return -1;
	    }
	}


      for (i=0; i<cmd_ptr->num_colours_in_detect_bubbles1_second_colour_list; i++)
	{
	  if (cmd_ptr->detect_bubbles1_second_colour_list[i]>=NUMBER_OF_COLOURS)
	    {
	      char tmp[] = "In detect_bubbles1 , you cannot specify a number > %d, which you specified at compile time as the number of colours supported. Recompile/see manual.\n";
	      if (strlen(tmp)>LEN_ERROR_STRING)
		{
		  printf("coding error - this string is too long:\n%s\n", tmp);
		  exit(1);
		}
	      strcpy(error_string, tmp);
	      return -1;
	    }
	}


    }


    if (cmd_ptr->detect_bubbles2==true)
    {
      int i;
      for (i=0; i<cmd_ptr->num_colours_in_detect_bubbles2_first_colour_list; i++)
	{
	  if (cmd_ptr->detect_bubbles2_first_colour_list[i]>=NUMBER_OF_COLOURS)
	    {
	      char tmp[] = "in detect_bubbles2, you cannot specify a number > %d, which you specified at compile time as the number of colours supported. Recompile/see manual.\n";
	      if (strlen(tmp)>LEN_ERROR_STRING)
		{
		  printf("coding error - this string is too long:\n%s\n", tmp);
		  exit(1);
		}
	      strcpy(error_string, tmp);
	      return -1;
	    }
	}




      for (i=0; i<cmd_ptr->num_colours_in_detect_bubbles2_second_colour_list; i++)
	{
	  if (cmd_ptr->detect_bubbles2_second_colour_list[i]>=NUMBER_OF_COLOURS)
	    {
	      char tmp[] = "In detect_bubbles2, you cannot specify a number > %d, which you specified at compile time as the number of colours supported. Recompile/see manual.\n";
	      if (strlen(tmp)>LEN_ERROR_STRING)
		{
		  printf("coding error - this string is too long:\n%s\n", tmp);
		  exit(1);
		}
	      strcpy(error_string, tmp);
	      return -1;
	    }
	}


    }

    if ( (cmd_ptr->expt_type!=Unspecified) && (cmd_ptr->genome_size==0) )
      {
      char tmp[] = "Experiment type is specified to allow genotyping, but to do that we also need genome length. Please also specify --genome_size\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
	
      }

  if ( (cmd_ptr->expt_type==EachColourADiploidSampleExceptTheRefColour) && (cmd_ptr->ref_colour==-1) )
    {
      char tmp[] = "If you specify the experiment as EachColourADiploidSampleExceptTheRefColour, then you need to specify the ref colour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if ( (cmd_ptr->expt_type==EachColourAHaploidSampleExceptTheRefColour) && (cmd_ptr->ref_colour==-1) )
    {
      char tmp[] = "If you specify the experiment as EachColourAHaploidSampleExceptTheRefColour, then you need to specify the ref colour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if ( (cmd_ptr->expt_type==EachColourADiploidSample) && (cmd_ptr->ref_colour!=-1) )
    {
      char tmp[] = "You have specified a reference colour, and yet also set experiment type to EachColourADiploidSample. Surely you meant EachColourADiploidSampleExceptTheRefColour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if ( (cmd_ptr->expt_type==EachColourAHaploidSample) && (cmd_ptr->ref_colour!=-1) )
    {
      char tmp[] = "You have specified a reference colour, and yet also set experiment type to EachColourAHaploidSample. Surely you meant EachColourAHaploidSampleExceptTheRefColour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }





  if ( (cmd_ptr->exclude_ref_bubbles==true) && (cmd_ptr->ref_colour==-1) )
    {
      char tmp[] = "If you specify --exclude_ref_bubbles, then you must specify a reference colour with --ref_colour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
    }

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
	  char tmp[]="Must specify either binary (--colour_list or --multicolour_bin), or sequence data (--se_list and/or --pe_list) as input \n";
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

  //if specify apply model selection at bubbles, must call bubbles!
  if ( (cmd_ptr->detect_bubbles1==false) && (cmd_ptr->apply_model_selection_at_bubbles==true) )
    {
      char tmp[]="If you specify --require_hw, you must also specify --detect_bubbles1";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }


  if ( (cmd_ptr->apply_model_selection_at_bubbles==true) && (cmd_ptr->genome_size==0) )
    {
      char tmp[]="If you specify --require_hw, you must also specify --genome_size";
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

  if (  (cmd_ptr->max_read_length>20000) && (cmd_ptr->align_given_list==false) )
    {
      char tmp[]="You are not allowed to set maximum read length >20000 when loading fasta/q into a graph. (However it IS allowed if you are aligning fasta/q to  preexisting graph). Since you have not specified --align, Cortexhas exited. Either reduce your max read length, or decide you are doing alignment\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  printf("coding error - this string is too long:\n%s\n", tmp);
	  exit(1);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }


  if (  (cmd_ptr->max_read_length>300000000) && (cmd_ptr->align_given_list==true) )
    {
      char tmp[]="If doing alignment, max read length is 300Mb (for whole chromosomes)\n";
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
      int binversion_in_header;
      boolean is_multicol_bin_ok = check_binary_signature(fp, cmd_ptr->kmer_size, BINVERSION, &num_m_cols, mean_readlens_ptrs, total_seqs_ptrs, &binversion_in_header);

      
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

  if (cmd_ptr->genotype_complex_site==true)
    {
      if (cmd_ptr->specified_max_var_len==false)
	{
	  char tmp[LEN_ERROR_STRING];
	  sprintf(tmp, "If you use --genotype_site, you must also specify --max_var_len and give the length of the longest allele in the (presumably multillelic) site you are genotyping\n");
	  strcpy(error_string, tmp);
	  return -1;
	  
	}
      if (cmd_ptr->num_colours_to_genotype + cmd_ptr->num_alleles_of_site >NUMBER_OF_COLOURS)
	{
	  printf("There is an inconsistency in your compile + command line specifications.\n");
	  printf("You want to genotype %d colours/samples, at a site with %d alleles, but cortex_var is only compiled for %d colours.\n", 
		 cmd_ptr->num_colours_to_genotype, cmd_ptr->num_alleles_of_site, NUMBER_OF_COLOURS
		 );
	  printf("cortex_var expects colours 0...%d-1 should each contain one of these alleles, and then you need one colour per sample and then perhaps one\n", cmd_ptr->num_colours_to_genotype);
	  printf(" more if you have a reference-minus-site colour. \n");
	  char tmp[LEN_ERROR_STRING];
	  sprintf(tmp, "These samples and alleles you have specified won't all fit in - recompile with more colours\n");
	  strcpy(error_string, tmp);
	  return -1;
	  
	}
      if (cmd_ptr->genome_size==0)
	{
	  char tmp[LEN_ERROR_STRING];
	  sprintf(tmp, "If you specify --genotype-site then you must also specify --genome_size\n");
	  strcpy(error_string, tmp);
	  return -1;
	  
	}

    }
  if ( (cmd_ptr->estimate_genome_complexity==true) && (cmd_ptr->format_of_input_seq==UNSPECIFIED) )
    {
	  char tmp[LEN_ERROR_STRING];
	  sprintf(tmp, "If you specify --estimate_genome_complexity then you must also specify --format\n");
	  strcpy(error_string, tmp);
	  return -1;
    }
  if ( (cmd_ptr->estimate_genome_complexity==true) && (cmd_ptr->max_read_length==0)  )
    {
	  char tmp[LEN_ERROR_STRING];
	  sprintf(tmp, "If you specify --estimate_genome_complexity then you must also specify --max_read_len\n");
	  strcpy(error_string, tmp);
	  return -1;
    }

  if ( (cmd_ptr->estimate_genome_complexity==true) && (NUMBER_OF_COLOURS<2)  )
    {
	  char tmp[LEN_ERROR_STRING];
	  sprintf(tmp, "If you specify --estimate_genome_complexity then you must also compile Cortex for at least 2 colours. Colour 0 should have your graph, and colour 1 is used internally by cortex\n");
	  strcpy(error_string, tmp);
	  return -1;
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
    
  //  extra_cmd_line_settings(&cmd_line);//some cmd_line settings depend on other options, and we need to wait til all parsing is done

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
// the list -1 is short hand for ALL colours, and a list prefixed by * means ALL colours except those on this list
int get_numbers_from_comma_sep_list(char* list, int* return_list, int max_len_return_list)
{
  int number[max_len_return_list];


  //exceptional case - if the list is just "-1" - no commas, just "-1", then this means ALL colours
  if (strcmp(list, "-1")==0)
    {
      int j;
      if (max_len_return_list<NUMBER_OF_COLOURS)
	{
	  printf("Passing array limit of %d int get_numbers_from_comma_sep_list which is less than the number of colours\n", max_len_return_list);
	  exit(1);
	}
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  return_list[j]=j;
	}
      
      return NUMBER_OF_COLOURS;
    }

  //exceptional case - if the list starts with "*," then all the numbers are to be EXCLUDED,
  if (list[0]=='*')
    {
      int colours_to_ignore[MAX_COLOURS_ALLOWED_TO_MERGE];

      int i;
      int current=0;
      colours_to_ignore[current]=0;

      for(i=1;i<strlen(list);i++){//start at 1!
	if (list[i]==',')
	  {
	    current++;
	    colours_to_ignore[current]=0;
	    if (current>=max_len_return_list)
	      {
		printf("Too many colours in list\n");
		return -1;
	      }
	  }
	else{
	  if ( (list[i]>='0' && list[i]<='9') ||  (i+1 == strlen(list))   ){
	    colours_to_ignore[current]=colours_to_ignore[current]*10 + (list[i]-'0');
	  }
	  else{
	    printf("This part of cmd-line argument is badly formatted: %s\n", list);
	    return -1;
	  }
	}
      }
      
      int index_in_returnlist=0;
      //this isn't very efficient but you only doing it once when you parse the cmd-line
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  int k;
	  boolean should_we_ignore=false;
	  for(k=0;k<=current;k++)
	    {
	      if (i==colours_to_ignore[k])
		{
		  should_we_ignore=true;
		}
	    }
	  if (should_we_ignore==false)
	    {
	      return_list[index_in_returnlist]=i;
	      index_in_returnlist++;
	    }

	}
      return index_in_returnlist;
    }

  else
    {
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
}


//RELIES on the list (arg1) ending in '\0'.
//expect the list to be in the form ;3;7;8 - PRECEDED and separated by semicolon
//returns -1 on error
// the list ;-1 is short hand for ALL colours, and a list prefixed by ;* means ALL colours except those on this list
int get_numbers_from_open_square_brack_sep_list(char* list, int* return_list, 
					int max_len_return_list)
{
  int number[max_len_return_list];


  //exceptional case - if the list is just "[-1" , 
  //                 then this means ALL colours
  if (strcmp(list, "[-1")==0)
    {
      int j;
      if (max_len_return_list<NUMBER_OF_COLOURS)
	{
	  printf("Passing array limit of %d int get_numbers_from_comma_sep_list which is less than the number of colours\n", max_len_return_list);
	  exit(1);
	}
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  return_list[j]=j;
	}
      
      return NUMBER_OF_COLOURS;
    }

  //exceptional case - if the list starts with ";*" then all the numbers are to be EXCLUDED,
  if (list[1]=='*')//we already know list[0] is [
    {
      int colours_to_ignore[MAX_COLOURS_ALLOWED_TO_MERGE];

      int i;
      int current=0;
      colours_to_ignore[current]=0;

      for(i=2;i<strlen(list);i++){//start at 2!
	if (list[i]=='[')
	  {
	    current++;
	    colours_to_ignore[current]=0;
	    if (current>=max_len_return_list)
	      {
		printf("Too many colours in list\n");
		return -1;
	      }
	  }
	else{
	  if ( (list[i]>='0' && list[i]<='9') ||  (i+1 == strlen(list))   ){
	    colours_to_ignore[current]=colours_to_ignore[current]*10 + (list[i]-'0');
	  }
	  else{
	    printf("This part of cmd-line argument is badly formatted: %s\n", list);
	    return -1;
	  }
	}
      }
      
      int index_in_returnlist=0;
      //this isn't very efficient but you only doing it once when you parse the cmd-line
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  int k;
	  boolean should_we_ignore=false;
	  for(k=0;k<=current;k++)
	    {
	      if (i==colours_to_ignore[k])
		{
		  should_we_ignore=true;
		}
	    }
	  if (should_we_ignore==false)
	    {
	      return_list[index_in_returnlist]=i;
	      index_in_returnlist++;
	    }

	}
      return index_in_returnlist;
    }

  else
    {
      int i;
      int current=0;
      number[current]=0;
      
      
      for(i=1;i<strlen(list);i++){//start from 1, as list[0] is a semicolon
	if (list[i]=='[')
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
}





//RELIES on the arg ending in '\0'.
//returns -1 on error

//note we implicitly assume colours 0...num_alleles are going to be one colour for each allele. 
//So the ref-minus-site colour must be > this, etc
int parse_genotype_site_argument(char* arg, int* colours_to_genotype_list, int* num_colours_to_genotype , int* ref_minus_site_colour, int* num_alleles,
				 int* start_gt_combin_num, int* end_gt_combin_num, char* fasta_file, AssumptionsOnGraphCleaning* assump,
				 int* wk_col1, int* wk_col2, boolean* using_1net, boolean* using_2net, char* file_1net_bins, double* min_llk)
{

  // we expect arg to be of this format: x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN   where z and A,B may be -1
  // x,y is a comma sep list of colours to genotype. z is the ref-minus-site colour. N is the number of alleles at the site (must be same as number of reads in the fasta)
  // A..B means genotype combinations A to B of the N choose 2 possible genotypes, p,q are spare colours, and yes/no specifies whether or not to use the 1net/2net form of likelihood,
  // and MIN is an optional argument to say ignore any genotype with log likelihood below MIN



  char delims[] = "[";
  char temp1[MAX_FILENAME_LEN];
  temp1[0]='\0';
  strcpy(temp1, arg);
  char* commaseplist = strtok(temp1, delims );
  if (commaseplist==NULL)
    {
      errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN, using the \"[\" character as a delimiter. you do not appear to have any of these delimiters - consult the manual");
    }
  
  if (strcmp(commaseplist, "-1")==0)
    {
      errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN where x,y denotes a comma-separated list of colours (eg 0,2,45) which are to be genotyped. In other parts of the cortex_var cmdline syntax, a \"-1\" is allowed to denote ALL colours. That makes no sense here. --genotype_site is designed for complex multiallelic sites, where the known alleles are loaded into colours 0..N-1 with --multicolour_bin, and then (optionally) another colour is used to refer to the reference-minus-site (see manual). So by definition, not all the colours in the graph are samples to be genotyped. In short, don't use -1 for A,B\n");
    }
  else
    {

      *num_colours_to_genotype  = get_numbers_from_comma_sep_list(commaseplist, colours_to_genotype_list, NUMBER_OF_COLOURS);
      if (*num_colours_to_genotype==-1)
	{
	  printf("cmdline error - Failed to parse out the colours which we should genotype in --genotype_site");
	  exit(1);
	}

    }
  
  char* refminussite_as_char  = strtok( NULL, delims );
  if (refminussite_as_char==NULL)
    {
      errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - you do not appear to have a \"z\" - consult the manual");
    }
  else
    {
      *ref_minus_site_colour = atoi(refminussite_as_char);
    }

  char* numalleles_as_char = strtok( NULL, delims );
  if (numalleles_as_char==NULL)
    {
      errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN. You have omitted the N, which is not permitted - consult the manual");
    }
  else
    {
      int num_a = atoi(numalleles_as_char);
      if (num_a<=0)
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - you have entered a negative value for N (the number of alleles), which is not permitted\n");
	}
      int num_cols_required_by_this_cmdline = *num_colours_to_genotype + num_a;
      if (*ref_minus_site_colour != -1)
	{
	  num_cols_required_by_this_cmdline++;
	}
      if (num_cols_required_by_this_cmdline>NUMBER_OF_COLOURS)
	{
	  printf("You have compiled for %d colours, but your commandline arg --genotype_site has arguments that together require at least %d\n", NUMBER_OF_COLOURS, num_cols_required_by_this_cmdline);
	}
      *num_alleles = atoi(numalleles_as_char);
    }

  char* startend_as_char = strtok( NULL, delims );

  if (startend_as_char==NULL)
    {
      errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - cortex cannot parse-find A,B\n");
    }
  else
    {
      // I am now ready to split out A,B but I'd rather finish strtokking my main string first.
      char* fa  = strtok( NULL, delims );
      if (fa==NULL)
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN  - unable to find a fasta in that string\n");
	}
      else if (test_file_existence(fa)==false)
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN  - unable to find/open the fasta file you mention");
	}
      else
	{
	  fasta_file[0]='\0';
	  strcpy(fasta_file, fa);
	}
      char* assump_as_char = strtok(NULL, delims);
      if (strcmp(assump_as_char,"CLEANED")==0)
	{
	  *assump = AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice;
	}
      else if (strcmp(assump_as_char, "UNCLEANED")==0)
	{
	  *assump = AssumeUncleaned;
	}
      else
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - after the fasta name, you appear to have writting something that is neither \"CLEANED\" nor \"UNCLEANED\" \n");	  
	}

      char* working_col1_as_char = strtok(NULL, delims);
      if (working_col1_as_char==NULL)
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - p and q should be colours in the graph (ie less than the max allowed number of colours that you compiled cortex for) - ypu appear to have left off p and q from your commandline\n");
	}
      int working_col1 = atoi(working_col1_as_char);
      if ( (working_col1>NUMBER_OF_COLOURS) || (working_col1<0) )
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - p and q should be colours in the graph (ie less than the max allowed number of colours that you compiled cortex for, but also should be unused - ie you have not loaded binaries into these colours\n");
	}
      char* working_col2_as_char = strtok(NULL, delims);
      if (working_col2_as_char==NULL)
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - p and q should be colours in the graph (ie less than the max allowed number of colours that you compiled cortex for) - ypu appear to have left off q from your commandline\n");
	}
      int working_col2 = atoi(working_col2_as_char);


      if ( (working_col2>NUMBER_OF_COLOURS) || (working_col2<0) )
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - p and q should be colours in the graph (ie less than the max allowed number of colours that you compiled cortex for, but also should be unused - ie you have not loaded binaries into these colours\n");
	}
      *wk_col1 = working_col1;
      *wk_col2 = working_col2;

      char* using_1net_as_char = strtok(NULL, delims);
      if (using_1net_as_char==NULL)
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - the final yes/no determines whether cortex uses  the full error model using edit distances of kmers from the expected genotype (yes) or not (no). You do not seem to have specified either\n");
	}
      else if (strcmp(using_1net_as_char, "yes")==0)
	{
	  *using_1net=true;
	  *using_2net=true;
	}
      else if (strcmp(using_1net_as_char, "no")==0)
	{
	  *using_1net=false;
	  *using_2net=false;
	}
      else
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN - the yes/no determines whether cortex uses  the full error model using edit distances of kmers from the expected genotype (yes) or not (no). You do not seem to have specified either\n");	  
	}


      char* min_llk_as_char = strtok(NULL, delims);
      if (min_llk_as_char==NULL)
	{
	  //do nothing, leave the min log likelihood as default
	}
      else
	{
	  int min_llk_as_int = atoi(min_llk_as_char);
	  if (min_llk_as_int<0)
	    {
	      *min_llk = (double) min_llk_as_int; 
	    }
	}


      /*
      //now get the file-list of 1-net binaries
      char* filelist_1net  = strtok( NULL, delims );
      if (filelist_1net==NULL)
	{
	  printf("[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN  - unable to find a file_1net in that string, so cortex interprets this as an active choice to use the simple error model (without the 1net etc). This is a valid choice to make! Just printing this to notify you, in case you have a typo.\n");
	}
      else if (test_file_existence(filelist_1net)==false)
	{
	  errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[file_1net  - unable to find/open the filelist_1net file you mention");
	}
      else
	{
	  file_1net_bins[0]='\0';
	  strcpy(file_1net_bins, filelist_1net);
	  *using_1net=true;
	  *using_2net=true;
	}
      */
      
      if (strcmp(startend_as_char, "-1")==0)
	{
	  *start_gt_combin_num = 1;
	  *end_gt_combin_num = (*num_alleles) * (*num_alleles -1) /2;  // N choose 2
	}
      else
	{
	  char delims2[] = ",";
	  char temp2[MAX_FILENAME_LEN];
	  temp2[0]='\0';
	  strcpy(temp2, startend_as_char);
	  char* startaschar = strtok(temp2, delims2);
	  if (startaschar==NULL)
	    {
	      errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN  - unable to split A,B out\n");
	    }
	  else
	    {
	      *start_gt_combin_num = atoi(startaschar);
	    }
	  char* endaschar = strtok( NULL, delims2 );
	  if (endaschar==NULL)
	    {
	      errx(1,"[--genotype_site] option requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN  - unable to split A,B out to get B\n");
	    }
	  else
	    {
	      *end_gt_combin_num = atoi(endaschar);
	    }
	}
    }
  
  
  return 0;
}



int parse_novelseq_args(char* arg, int* array_colours_to_look_in, int* num_cols_in_look_list,
			int* array_colours_to_avoid,  int* num_cols_in_avoid_list,
			int* min_contig_len, int* min_percentage_novel, char* outfile)
{
  //we expect arg tyo be of this format  a,b/c,d/x/y/filename
  // a,b,c,d are colours. a,b are the colours in which we will look for supernodes. c,d are "the reference" - ie the colours to avoid. 
  // you can use -1 to denote all colours instead of a,b. But that's not allowed for c,d - waste of time trying to find stuff not in the entire DBG
  // x is the min contig len
  // y is the min proportion of kmers in a contig that must be novel. Typically this will be 1. 
  // filename - output filename



  char delims[] = "/";
  char temp1[MAX_FILENAME_LEN];
  temp1[0]='\0';
  strcpy(temp1, arg);//we checked before callling this that it fits
  char* commaseplist1 = strtok(temp1, delims );
  if (commaseplist1==NULL)
    {
      errx(1,"[--print_novel_conrigs] option requires an argument of the form a,b/c,d/x/y/filename  you do not appear to have any of these / delimiters - consult the manual");
    }

  *num_cols_in_look_list = get_numbers_from_comma_sep_list(commaseplist1, array_colours_to_look_in, NUMBER_OF_COLOURS);
  if (*num_cols_in_look_list==-1)
    {
      printf("Badly formatted/chosen arguments for --print_novel_contigs. Consult the manual\n");
      exit(1);
    }



  char* commaseplist2 = strtok(NULL, delims );
  if (commaseplist2==NULL)
    {
      errx(1,"[--print_novel_conrigs] option requires an argument of the form a,b/c,d/x/y/filename  you do not appear to have anything after the first /\n");
    }

  *num_cols_in_avoid_list = get_numbers_from_comma_sep_list(commaseplist2, array_colours_to_avoid, NUMBER_OF_COLOURS);
  if (*num_cols_in_avoid_list==-1)
    {
      printf("Badly formatted/chosen arguments for --print_novel_contigs. Consult the manual\n");
      exit(1);
    }
  else if (*num_cols_in_avoid_list>=NUMBER_OF_COLOURS)
    {
      printf("[--print_novel_contigs] option requires an argument of the form a,b,../c,d,../x/y/filename. You have listed ALL the colours in the graph in c,d (either by enumerating them all, or by entering -1). this is meaningless - the idea is to look for stuff in the graph that is NOT in colours c,d,etc\n");
      exit(1);
    }
  char* min_contig_len_as_char = strtok(NULL, delims ); 
  if (min_contig_len_as_char==NULL)
    {
      errx(1,"[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, but you do not seem to have anything after c,d\n");
    }
  else if (isNumeric(min_contig_len_as_char))
    {
	  *min_contig_len = atoi(min_contig_len_as_char);
    }
  else
    {
      errx(1,"[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, where x is the minimum contig length (in bp) - you have entered a non-numeric thing for x\n");
    }
  
  char* percentage_novel_as_char = strtok(NULL, delims );
  if (percentage_novel_as_char==NULL)
    {
      errx(1,"[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, but you do not seem to have anything after x\n");
    }
  else if (isNumeric(percentage_novel_as_char))
    {
      int tmp = atoi(percentage_novel_as_char);
      if ( (tmp>0) && (tmp<=100) )
	{
	  *min_percentage_novel =  tmp;
	}
      else if (tmp<=0)
	{
	  printf("[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, where y is the minimum percentage of the kmers in a contig which we require to NOT be in any colours c,d,.... ie if we want to be very stringent, we set y=100, and demand ALL kmers be novel (ie outside colours c,d..). You have chosen a value for y which is <0, which we do not allow\n");
	  exit(1);
	}
      else if (tmp>100)
	{
	  printf("[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, where y is the minimum percentage of the kmers in a contig which we require to NOT be in any colours c,d,.... ie if we want to be very stringent, we set y=100, and demand ALL kmers be novel (ie outside colours c,d..). You have chosen a value for y which is >100, which is meaningless\n");
	  exit(1);
	  
	}
    }
  else
    {
	  printf("[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, where y is the minimum percentage of the kmers in a contig which we require to NOT be in any colours c,d,.... ie if we want to be very stringent, we set y=100, and demand ALL kmers be novel (ie outside colours c,d..). You have chosen a NON_NUMERIC value for y\n");
	  exit(1);
      
    }

  char* output_filename = strtok( NULL, delims );
  if (output_filename==NULL)
    {
      errx(1,"[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename - you have left off the filename\n");
    }
  else if (access(output_filename,F_OK)==0){
    errx(1,"You cannot specify filename [%s] as output file for --print_novel_contigs, - file exists!",output_filename);
  }
  outfile[0]='\0';
  strcpy(outfile, output_filename);
  
  return 0;
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

	    /*
	    if (num_left_colours==-2)//this means you have to merge ALL colours;
	      {
		int j;
		num_left_colours=NUMBER_OF_COLOURS;
		for (j=0; j<NUMBER_OF_COLOURS; j++)
		  {
		    left_colours[j]=j;
		  }
	      }
	    if (num_right_colours==-2)//this means you have to merge ALL colours;
	      {
		num_right_colours=NUMBER_OF_COLOURS;
		int j;
		for (j=0; j<NUMBER_OF_COLOURS; j++)
		  {
		    right_colours[j]=j;
		  }
	      }
	    */
	    

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




  //returns -1 on error, 0 otherwise. Parses argument3 and gets the colours into arg1 and arg2 
int parse_commasep_or_open_square_brack_sep_list(CmdLine* cmd, char* arg, int len_arg, char* text_for_error_describing_which_option_this_is, boolean commasep)
{

  if (len_arg<MAX_LEN_DETECT_BUB_COLOURINFO)
	  {
	    int k;

	    for (k=0; k<len_arg; k++)
	      {
		if (arg[k]=='/') 
		  {
		    
		    printf("%s option accepts a one comma-separated list of integers, or a open square-bracket-[ preceded-and-separated list, but in either case, with no forward-slash /. Perhaps you are confusing with --detect_bubbles1.", 
			   text_for_error_describing_which_option_this_is);		
		    return -1;
		  }
	      }

	    char list[MAX_LEN_DETECT_BUB_COLOURINFO];
	    set_string_to_null(list, MAX_LEN_DETECT_BUB_COLOURINFO);
	    strncpy(list, arg, strlen(arg));
	    
	    int list_colours[MAX_COLOURS_ALLOWED_TO_MERGE];
	    int num_list_colours;
	    if (commasep==true)
	      {
		num_list_colours= get_numbers_from_comma_sep_list(list, list_colours, MAX_COLOURS_ALLOWED_TO_MERGE);
	      }
	    else
	      {
		num_list_colours= get_numbers_from_open_square_brack_sep_list(list, list_colours, MAX_COLOURS_ALLOWED_TO_MERGE);
	      }

	    if (num_list_colours==-1) 
	      {
		return -1; //error.
	      }
	    int j;
	    for (j=0; j<num_list_colours; j++)
	      {
		if (list_colours[j]<0)
		  {
		    printf("Not allowed negative numbers in list, except for special case where there is just one number, and it is -1, which means \"all colours\". ");
		  }
	      }

	    //*num_cols= num_list_colours;
	    //copy_list(list_colours, num_list_colours, list_of_colours, num_list_colours);
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



