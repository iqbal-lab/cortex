
/*
 * Copyright 2009-2013 Zamin Iqbal and Mario Caccamo 
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
  cmd_line.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <err.h>

// cortex_var headers
#include "cmd_line.h"
#include "file_reader.h"
#include "file_format.h"

#define MAX_NUM_DIGITS_IN_COLOUR_ENTERED_ON_CMDLINE_COMMASEP 4
//#define LEN_ERROR_STRING 400
#define MAX_LINE 500

// Supress a warning from the compiler about an unused variable
#define UNUSED(x) ((void)(x))

int isNumeric (const char * s)
{
  if (s == NULL || *s == '\0')
    return 0;
  char * p;
  double d = strtod(s, &p);
  // to supress warning from compiler
  UNUSED(d);
  return *p == '\0';
}

boolean check_files_in_list_exist_ignoring_ref_line(char* filelist, char* error_msg, int ref_colour)
{
  error_msg[0]='\0';
  /*
  if (ref_colour==-1)
    {
      strcat(error_msg,"You should not call --estim_e_with_snps unless you also specify --ref_colour\n");
      return false;
      }*/
    if (filelist[0]=='\0')
    {
      strcat(error_msg,"Empty string in check_files_in_list_exist_ignoring_ref_line - call Zam - I only put this error msg in out of paranoia\n");
      return false;
    }

  FILE* fp = fopen(filelist, "r");
  if (fp==NULL)
    {
      die("Cannot open %s\n", filelist);
    }

  char filename[MAX_LINE];
  int count=-1;
  while (feof(fp) ==0)
    {
      if (fgets(filename,MAX_LINE, fp) !=NULL)
        {
	  count++;
	  if (count != ref_colour)
	    {
	      char* p;
	      if ((p = strchr(filename, '\n')) != NULL)
		*p = '\0';
	      
	      if (access(filename,R_OK)==-1) 
		{
		  sprintf(error_msg,
              "File Z%sZ on line %d (numbering lines from 0) (corresponding to "
              "colour %d) of file %s cannot be accessed - bad path?\n",
              filename, count, count, filelist);
		  return false;
		}
	    }
	}
    }
  fclose(fp);

  if (count!=NUMBER_OF_COLOURS-1)
    {
      sprintf(error_msg,
              "File %s is supposed to contain one line per colour in your graph "
              "(the assumption is your graph contains a ref-colour, plus a set "
              "of colours, all of which contain data from samples for who you "
              "want to estimate the sequencing error rate). However in fact it "
              "contains %d lines\n", filelist, count);
      return false;
      
    }
  return true;


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
      die("Cannot open %s\n", file);
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
      die("Cannot open %s\n", file);
    }


  GraphInfo* ginfo=graph_info_alloc_and_init();
  BinaryHeaderErrorCode ecode=EValid;
  BinaryHeaderInfo binfo;
  initialise_binary_header_info(&binfo, ginfo);

  query_binary_NEW(fp, &binfo, &ecode,0);
  fclose(fp);
  graph_info_free(ginfo);
  if (binfo.number_of_colours>1)
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
"\nusage: to build a single-colour binary:\ncortex_var --se_list <filename> --pe_list <filename> --format FASTQ --quality_score_threshold 5 --remove_pcr_duplicates --remove_low_coverage_supernodes 1 --dump_binary some_name.ctx\n" \
"\nusage: to build a multicolour graph from single-colour graphs and call variants between colours 1 and 2:\ncortex_var --colour_list <filename> --detect_bubbles1 1/2 --output_bubbles1 vars_between_cols1_and_2\n" \
"\nusage: to load a multicolour graph from single-colour graphs and call heterozygous variants in colour 0:\ncortex_var --colour_list <filename> --detect_bubbles1 0/0 --output_bubbles1 hets_in_colour_0\n" \
"\n" \


"   [--help] \t\t\t\t\t\t\t=\t This help screen.\n\n" \
" \n ** DATA LOADING ** \n\n"\
  // -v
  //"   [--format TYPE] \t\t\t\t\t\t=\t File format for input in se_list and pe_list. All files assumed to be of the same format.\n\t\t\t\t\t\t\t\t\t Type must be FASTQ or FASTA\n"
"   [--colour_list FILENAME] \t\t\t\t\t=\t File of filenames, one per colour. n-th file is a list of\n\t\t\t\t\t\t\t\t\t single-colour binaries to be loaded into colour n.\n\t\t\t\t\t\t\t\t\t Cannot be used with --se_list or --pe_list \n\t\t\t\t\t\t\t\t\t Optionally, this can contain a second column, containing a sample identifier/name for each colour\n" \
"   [--multicolour_bin FILENAME] \t\t\t\t=\t Filename of a multicolour binary, will be loaded first, into colours 0..n.\n\t\t\t\t\t\t\t\t\t If using --colour_list also, those will be loaded into subsequent colours, after this.\n" \
"   [--se_list FILENAME] \t\t\t\t\t=\t List of single-end fasta/q to be loaded into a single-colour graph.\n\t\t\t\t\t\t\t\t\t Cannot be used with --colour_list\n\t\t\t\t\t\t\t\t\t Optionally, the first line is allowed to have, after the filename, a tab, and then a sample identifier\n" \
"   [--pe_list FILENAME] \t\t\t\t\t=\t Two filenames, comma-separated: each is a list of paired-end fasta/q to be \n\t\t\t\t\t\t\t\t\t loaded into a single-colour graph. Lists are assumed to ordered so that \n\t\t\t\t\t\t\t\t\t corresponding paired-end fasta/q files are at the same positions in their lists.\n\t\t\t\t\t\t\t\t\t Currently Cortex only use paired-end information to remove\n\t\t\t\t\t\t\t\t\t PCR duplicate reads (if that flag is set).\n\t\t\t\t\t\t\t\t\t Cannot be used with --colour_list\n\t\t\t\t\t\t\t\t\t  Optionally, the first line is allowed to have, after the filename, a tab, and then a sample identifier\n" \
"   [--kmer_size INT] \t\t\t\t\t\t=\t Kmer size. Must be an odd number.\n" \
"   [--mem_width INT] \t\t\t\t\t\t=\t Size of hash table buckets (default 100).\n" \
  //-g 
"   [--mem_height INT] \t\t\t\t\t\t=\t Number of buckets in hash table in bits (default 10). \n\t\t\t\t\t\t\t\t\t Actual number of buckets will be 2^(the number you enter)\n" \
  // -n 
"   [--fastq_offset INT] \t\t\t\t\t=\t Default 33, for standard fastq.\n\t\t\t\t\t\t\t\t\t Some fastq directly from different versions of Illumina machines require different offsets.\n" \
  // -o
"   [--sample_id STRING] \t\t\t\t\t=\t (Only) if losding fasta/q, you can use this option to set the sample-identifier.\n\t\t\t\t\t\t\t\t\t This will be saved in any binary file you dump.\n" \
  // -p
"   [--dump_binary FILENAME] \t\t\t\t\t=\t Dump a binary file, with this name (after applying error-cleaning, if specified).\n" \
  // -T
"   [--subsample FRAC] \t\t\t\t\t=\t Subsample input data, taking fraction FRAC of data. If you want to dump a binary after having done this, use --dump_binary\n" \
  // -w
"   [--max_read_len] \t\t\t\t\t\t=\t (Unlike previous versions of Cortex) now required only if using --gt or --dump_filtered_readlen_distribution.\n" \
" \n**** FILTERING AND ERROR CORRECTION/CLEANING OPTIONS ****\n\n"\
  // -m
"   [--quality_score_threshold INT] \t\t\t\t=\t Filter for quality scores in the input file (default 0).\n" \
  // -j
"   [--remove_pcr_duplicates] \t\t\t\t\t=\t Removes PCR duplicate reads by ignoring read pairs if both \n\t\t\t\t\t\t\t\t\t reads start at the same k-mer as a previous read,\n" \
  // -k
"   [--cut_homopolymers INT] \t\t\t\t\t=\t Breaks reads at homopolymers of length >= this threshold.\n\t\t\t\t\t\t\t\t\t (i.e. max homopolymer in filtered read==threshold-1, and New read starts after homopolymer)\n" \
  // -O
"   [--remove_low_coverage_supernodes INT]\t\t\t\t=\t Remove all supernodes where max coverage is <= the limit you set. Recommended method.\n" \
  // -B
"   [--remove_low_coverage_kmers INT] \t\t\t\t=\t Filter for kmers with coverage less than or equal to  threshold. Not recommended. See manual and our paper for why\n"  \
  // -E
"   [--load_colours_only_where_overlap_clean_colour INT] \t=\t Only load nodes from binary files in the colour-list when they overlap a\n\t\t\t\t\t\t\t\t\t specific colour (e.g. that contains a cleaned pooled graph);\n\t\t\t\t\t\t\t\t\t requires you to specify this particular colour. You must have loaded that colour beforehand, using --multicolour_bin\n"  \
  // -F
"   [--successively_dump_cleaned_colours SUFFIX] \t\t=\t Only to be used when also using --load_colours_only_where_overlap_clean_colour and --multicolour_bin\n\t\t\t\t\t\t\t\t\t Used to allow error-correction of low-coverage data on large numbers of individuals with large genomes.\n\t\t\t\t\t\t\t\t\t Requires the user specify a suffix which will be added to the names of cleaned binaries. See manual for details.\n"  \
  // -S
  "   [--err_correct]\t\t\t\t\t\t=\t FIVE comma-separated arguments:\n\t\t\t\t\t\t\t\t\t(1) a filelist (listing FASTQ for correction)\n\t\t\t\t\t\t\t\t\t(2) a string/suffix less than 50 characters long \n\t\t\t\t\t\t\t\t\t  (/path/to/zam.fastq--->correction process --> outdir/zam.fastq.suffix)\n\t\t\t\t\t\t\t\t\t(3) an output directory\n\t\t\t\t\t\t\t\t\t(4) either 0 or 1. \n\t\t\t\t\t\t\t\t\t  O means discard a read if it has a low quality uncorrectable base, \n\t\t\t\t\t\t\t\t\t  1 means correct as many bases as you can, but if there is a low quality uncorrectable base, leave it uncorrected.\n\t\t\t\t\t\t\t\t\t(5) Experimental option - A number N of greedy extra bases to pad the end of the read with. \n\t\t\t\t\t\t\t\t\t  I suggest you use 0 (zero) for this, unless you are Jared Simpson\n\t\t\t\t\t\t\t\t\t   or Shane McCarthy, when you would enter 20.\n\t\t\t\t\t\t\t\t\t   The 1000 Genomes project wanted to be able to use a \n\t\t\t\t\t\t\t\t\t   particular BWT compression scheme, which required \n\t\t\t\t\t\t\t\t\t   the addition of N bases of sequence following the read, \n\t\t\t\t\t\t\t\t\t   from where it mapped in the reference. \n\t\t\t\t\t\t\t\t\t   Ignore this unless you are me, Jared or Shane for now. \n\t\t\t\t\t\t\t\t\t   In addition you must use --quality_score_threshold and --max_read_len.\n\t\t\t\t\t\t\t\t\tFinally, if you use --list_ref_fasta then Cortex will mark the + strand in the graph,\n\t\t\t\t\t\t\t\t\tand will print all reads in the + orientation \n"
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
"   [--detect_bubbles1 COMMA_SEP_COLOURS/COMMA_SEP_COLOURS] \t=\t Find all the bubbles in the graph where the two branches lie in the specified colours\n\t\t\t\t\t\t\t\t\t (after applying all specified actions on graph).\n\t\t\t\t\t\t\t\t\t Typical use would be --detect_bubbles1 1/1 to find hets in colour 1,\n\t\t\t\t\t\t\t\t\t or --detect_bubbles1 0/1 to find homozygous non-reference bubbles where one branch is in colour 0 (and not colour1)\n\t\t\t\t\t\t\t\t\t and the other branch is in colour1 (but not colour 0).\n\t\t\t\t\t\t\t\t\t However, one can do more complex things:\n\t\t\t\t\t\t\t\t\t e.g.  --detect_bubbles1 1,2,3/4,5,6 to find bubbles where one branch is in 1,2 or 3 (and not 4,5 or 6)\n\t\t\t\t\t\t\t\t\t and the other branch in colour 4,5 or 6 (but not 1,2, or 3).\n\t\t\t\t\t\t\t\t\t See the manual for a detailed description of the syntax. It is possible to\n\t\t\t\t\t\t\t\t\t specify \"all colour\" rather than enumerate them explicitly, and also to exclude colours\n" \
  // -s
"   [--output_bubbles1 FILENAME]\t\t\t\t\t=\t Bubbles called in detect_bubbles1 are dumped to this file.\n" \
  // -v
"   [--high_diff FLOAT]\t\t\t\t\t\t=\t Only call bubbles that (user-specified) minimum  allele-balance difference \n\t\t\t\t\t\t\t\t\t between at least 2 colours (ignores reference colour)\n\t\t\t\t\t\t\t\t\t Designed to be used when comparing pooled populations\n\t\t\t\t\t\t\t\t\t e.g --detect_bubbles -1/-1 --high_diff 0.9 will call bubbles where \n\t\t\t\t\t\t\t\t\t for example allele_balance(colour0)=0.01 and allele_balance(colour1)=0.95\n\t\t\t\t\t\t\t\t\t If you use this option, you must also set --genome_size" \

  // -x
"   [--print_colour_coverages]\t\t\t\t\t=\t Print coverages in all colours for supernodes and variants.\n" \
  // -M
"   [--exclude_ref_bubbles]\t\t\t\t\t=\t If you have specified --ref_colour, this will exclude any bubble in that colour from being called by the Bubble Caller.\n" \
  // -u
  //hidden from public
  //  "   [--estim_e_with_snps FILENAME]\t\t\t\t\t=\t Use SNPs alleles known (eg from SNP-chip genotyping) to estimate\n\t\t\t\t\t 
  //the sequencing error rate for each colour. Give a list of fasta, one per colour. For //ref colour, line is ignored.\n"

  // -l
"   [--path_divergence_caller [args]] \t\t\t\t\t= Make Path Divergence variant calls. Arguments can be specified in 2 ways.\n\t\t\t\t\t\t\t\t\t Option 1. Calls once, comparing reference and one colour (or union)\n\t\t\t\t\t\t\t\t\t e.g. --path_divergence_caller 1,2 --ref_colour 0 will look for differences\n\t\t\t\t\t\t\t\t\t between the union of colours 1,2 and the reference in colour 0\n\t\t\t\t\t\t\t\t\t Option2. Make several successive independent runs of the PD caller, each time against a different colour\n\t\t\t\t\t\t\t\t\tTo do this, use a square open bracket [ PRECEDED AND SEPARATED list\n\t\t\t\t\t\t\t\t\t For example --path_divergence_caller [2[3[10 --ref_colour 0 will make calls on samples 2 then 3 then 10)\n\t\t\t\t\t\t\t\t\t all output to the same file, with globally unique variant names. The caller will call against each colour in turn\n\t\t\t\t\t\t\t\t\t You must also specify --ref_colour and --list_ref_fasta\n" \
  // -I
"   [--path_divergence_caller_output PATH_STUB]\t\t\t=\t Specifies the path and beginning of filename of Path Divergence caller output file.\n" \
  // -i
"   [--ref_colour INT] \t\t\t\t\t\t=\t Colour of reference genome.\n" \
 // -z
"   [--list_ref_fasta FILENAME] \t\t\t\t\t=\t File listing reference chromosome fasta file(s); needed for path-divergence calls. \n" \
  // -t
"   [--gt INPUT,OUTPUT,{BC|PD}] \t\t\t\t\t=\t Given a file of calls in Cortex output format (5p, br1, br2, 3p), genotype all colours in the graph\n\t\t\t\t\t\t\t\t\t and output to specified filename. All calls must be either from the BubbleCaller or PathDivergence (not a mixture)\n\t\t\t\t\t\t\t\t\t and you specify this with either BC or PD. eg --gt infile,outfile,BC\n" \

"\n\n **** ADVANCED OPTIONS **** \n\n"\
  // -P  
"   [--experiment_type]\t\t\t\t\t\t=\t The statistical models for determining genotype likelihoods, and for deciding if bubbles are repeat or variants,\n\t\t\t\t\t\t\t\t\t require knowledge of whether each sample is a separate diploid/haploid individual. \n\t\t\t\t\t\t\t\t\t Enter type of experiment (EachColourADiploidSample, EachColourADiploidSampleExceptTheRefColour, \n\t\t\t\t\t\t\t\t\t EachColourAHaploidSample,EachColourAHaploidSampleExceptTheRefColour). \n\t\t\t\t\t\t\t\t\t This is only needed for determining likelihoods, so ignore this is you are pooling samples within a colour (support to be added for this later).\n" \
  // Q
"   [--estimated_error_rate]\t\t\t\t\t=\t If you have some idea of the sequencing error rate (per base-pair), enter it here. eg 0.01. Currently used in calculating likelihoods\n"\
  // -L
"   [--genome_size]\t\t\t\t\t\t=\t If you specify --experiment_type, and therefore want to calculate likelihoods, you must also specify the (estimated) genome size in bp.\n" \

  // -G
"   [--align FILENAME,{output binary name|no}] \t\t\t\t=\t Aligns a list of fasta/q files to the graph, \n\t\t\t\t\t\t\t and prints coverage of each kmer in each read in each colour.\n\t\t\t\t\t\t\t Takes two arguments. First, a LIST of fasta/q. \n\t\t\t\t\t\t\t Second, either an output filename (if you want it to dump a binary of the part of the graph touched by the alignment) OR just \"no\" \n\t\t\t\t\t\t\t Must also specify --align_input_format, and --max_read_len\n"  \
  // -H
"   [--align_input_format TYPE] \t\t\t\t\t=\t --align requires a list of fasta or fastq. This option specifies the input format as LIST_OF_FASTQ or LIST_OF_FASTA\n"  \

  // T
  //hidden from public
  //"   [--estimate_genome_complexity FILENAME]\t\t\t=\t Print estimated genome complexity by reading up to 10,000 reads from the given file\n"
  //"\n\n **** OPTIONS THAT I MAY REMOVE - NOT CLEAR ANYONE HAS A GOOD USE FOR THEM - EMAIL ME IF YOU LIKE IT\n\n"
  // -J
"   [--colour_overlaps COMMA_SEP_COLOURS/COMMA_SEP_COLOURS]\t=\t Compares each coloured subgraph in the first list with all of the \n\t\t\t\t\t\t\t\t\t coloured subgraphs in the second list. Outputs a matrix to stdout;\n\t\t\t\t\t\t\t\t\t (i,j)-element is the number of nodes in both colour-i (on first list)\n\t\t\t\t\t\t\t\t\t and colour-j (on second list).\n" \

"\n\n**** EARLY ACCESS/BETA OPTIONS **** \n\n"\
  // R
"   [--genotype_site]\t\t\t\t\t\t=\t (Beta code!) Genotype a single (typically multiallelic) site. Syntax is slightly complex. \n\t\t\t\t\t\t\t\t\t Requires an argument of the form x,y[z[N[A,B[fasta[<CLEANED|UNCLEANED>[p[q[<yes|no>[MIN.\n\t\t\t\t\t\t\t\t\t x,y is a comma-sep list of colours to genotype.\n\t\t\t\t\t\t\t\t\t z is the reference-minus-site colour.\n\t\t\t\t\t\t\t\t\t N is the number of alleles for this site (which cortex assumes are loaded\n\t\t\t\t\t\t\t\t\t in a multicolour_bin containing those alleles first, one per colour).\n\t\t\t\t\t\t\t\t\t Cortex will genotype combinations A through B of the N choose 2 possible genotypes (allows parallelisation);\n\t\t\t\t\t\t\t\t\t fasta is the file listing one read per allele.\n\t\t\t\t\t\t\t\t\t CLEANED or UNCLEANED allows Cortex to tailor its genotyping model.\n\t\t\t\t\t\t\t\t\t p,q are two free/unused colours that Cortex will use internally. \n\t\t\t\t\t\t\t\t\t yes/no specifies whether to use the more sophisticated error model, which is still in development. \n\t\t\t\t\t\t\t\t\t I recommend you stick with \"no\" for now.\n\t\t\t\t\t\t\t\t\t The final argument, MIN, is optional and allows performance speedup\n\t\t\t\t\t\t\t\t\t by disarding any genotype with log-likelihood<MIN.\n\t\t\t\t\t\t\t\t\t See manual for details.Must also specify --max_var_len to give the length of the longest allele\n"\
  // -V
"   [--print_novel_contigs]\t\t\t\t\t=\t Allows printing of novel sequence absent from a reference (or more generally, absent from a set of colours)\n\t\t\t\t\t\t\t\t\t Takes arguments in this format a,b,../c,d,../x/y/<output filename>\n\t\t\t\t\t\t\t\t\t Cortex will find supernodes in the union graph of colours a,b,..\n\t\t\t\t\t\t\t\t\t Typically the list c,d,.. of colours is just one colour  - that of the reference.\n\t\t\t\t\t\t\t\t\t Cortex will print contigs (supernodes) to the output file which satisfy the following criteria\n\t\t\t\t\t\t\t\t\t Contigs must be at least x base pairs long\n\t\t\t\t\t\t\t\t\t The percentage (as integer) of kmers in the contig which are present in ANY of the colours c,d,... must be at most 1-y. \n\t\t\t\t\t\t\t\t\t i.e. y is the minimum proportion of novel kmers in a contig. Typically this is 100.\n\t\t\t\t\t\t\t\t\t We ignore the first and last kmer of the contig, as these will typically connect to the reference\n  " \


  // -K
  "   [--pan_genome_matrix]\t\t\t\t\t\t=\t Pass in a gene fasta, and get percent overlap with all colours\n"

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
  c->high_diff=false;
  c->min_allele_balance_diff = 0.9; //when looking for high diff sites.
  c->subsample=false;
  c->subsample_propn=(float)1.0;//by default, take all reads

  c->get_pan_genome_matrix=false;
  set_string_to_null(c->pan_genome_genes_fasta, MAX_FILENAME_LEN);

  c->do_err_correction=false;
  c->err_correction_policy=DontWorryAboutLowQualBaseUnCorrectable;
  c->do_greedy_padding=false;
  c->greedy_pad=0;

  c->estimate_copy_num=false;
  set_string_to_null(c->copy_num_fasta, MAX_FILENAME_LEN);
  set_string_to_null(c->copy_num_output, MAX_FILENAME_LEN);

  c->for_each_colour_load_union_of_binaries=false;
  //novelsseq stuff
  c->use_snp_alleles_to_estim_seq_err_rate=false;
  c->entered_sampleid_as_cmdline_arg=false;
  c->loaded_sample_names=false;
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
  set_string_to_null(c->manually_entered_seq_error_rates_file, MAX_FILENAME_LEN);
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
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      set_string_to_null(c->colour_sample_ids[i], MAX_LEN_SAMPLE_NAME);
      strcpy(c->colour_sample_ids[i], "undefined");
    }
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
  set_string_to_null(c->file_of_calls_to_be_genotyped, MAX_FILENAME_LEN);
  set_string_to_null(c->output_genotyping, MAX_FILENAME_LEN);

  c->which_caller_was_used_for_calls_to_be_genotyped=BubbleCaller;
  // set_string_to_null(c->filelist_2net_binaries_for_alleles, MAX_FILENAME_LEN);

  //booleans
  c->exclude_ref_bubbles=false;
  c->cut_homopolymers = false;
  c->remove_pcr_dups=false;
  //c->clip_tips=false;
  c->print_colour_coverages=false;
  c->print_median_covg_only=false;
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
  c->do_genotyping_of_file_of_sites=false;

  c->num_colours_to_genotype=0;
  initialise_int_list(c->list_colours_to_genotype, NUMBER_OF_COLOURS);
  c->colour_of_reference_with_site_excised=-1;
  c->num_alleles_of_site=0;
  c->first_genotype_to_calc_likelihoods_for=0;
  c->last_genotype_to_calc_likelihoods_for=0;
  c->genotype_complex_site=false;
  c->using_1net=false;
  c->using_2net=false;
  c->assump_for_genotyping=AssumeUncleaned;
  c->min_acceptable_llk=-999999999;
  return 1;
}


CmdLine* cmd_line_alloc()
{
  CmdLine* cmd = (CmdLine*) malloc(sizeof(CmdLine));
  if (cmd==NULL)
    {
      die("Out of memory before we even start Cortex, cannot alloc space to hold the commandline variables! Abort\n");
    }
  cmd->colour_sample_ids=(char**) malloc(sizeof(char*) * NUMBER_OF_COLOURS);
  if (cmd->colour_sample_ids==NULL)
    {
      die("Out of memory before we even start Cortex, cannot alloc space to hold the commandline variables! Abort\n");      
    }
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      cmd->colour_sample_ids[i]=(char*) malloc(sizeof(char)*(MAX_LEN_SAMPLE_NAME+1));
      if (cmd->colour_sample_ids[i]==NULL)
	{
	  die("Out of memory before we even start Cortex, cannot alloc space to hold the commandline variables! Abort\n");      
	}
      
    }
 //string buffers
 cmd->err_correction_filelist = strbuf_new();
 cmd->err_correction_outdir   = strbuf_new();
 cmd->err_correction_suffix   = strbuf_new();
 return cmd;
}

void cmd_line_free(CmdLine* cmd)
{
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      free(cmd->colour_sample_ids[i]);
    }
  free(cmd->colour_sample_ids);
  strbuf_free(cmd->err_correction_filelist);
  strbuf_free(cmd->err_correction_outdir);
  strbuf_free(cmd->err_correction_suffix);
  free(cmd);
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
    {"sample_id",required_argument,NULL,'o'},
    {"dump_binary",required_argument,NULL,'p'},
    {"output_supernodes",required_argument,NULL,'q'},
    {"detect_bubbles1",required_argument, NULL, 'r'},
    {"output_bubbles1",required_argument, NULL, 's'},    
    {"gt",required_argument, NULL, 't'},
    {"estim_e_with_snps",required_argument, NULL, 'u'},    
    {"high_diff",required_argument,NULL,'v'},
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
    {"pan_genome_matrix",required_argument,NULL,'K'},
    {"genome_size",required_argument,NULL,'L'},
    {"exclude_ref_bubbles",no_argument,NULL,'M'},
    {"estimate_copy_num",required_argument,NULL,'N'},
    {"remove_low_coverage_supernodes",required_argument,NULL,'O'},
    {"experiment_type", required_argument, NULL, 'P'},
    {"estimated_error_rate", required_argument, NULL, 'Q'},
    {"genotype_site", required_argument, NULL, 'R'},
    {"err_correct", required_argument, NULL, 'S'},
    {"subsample", required_argument, NULL, 'T'},
    {"print_median_covg_only", no_argument, NULL, 'U'},
    {"print_novel_contigs", required_argument, NULL, 'V'},
    {0,0,0,0}	
  };
  

  //do not change this! Only really matters for testing, but getopt_long uses 
  //variables which are not local to this function, so when we run tests, and call this function
  // repeatedly, those variable values are carried across. The following line resets this.
  optind=1;
  
 
  opt = getopt_long(argc, argv, "ha:b:c:d:e:f:g:i:jk:l:m:n:o:p:q:r:s:t:u:v:w:xy:z:A:B:CD:E:F:G:H:I:J:K:L:MN:O:P:Q:R:S:T:UV:", long_options, &longopt_index);

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
	  errx(1,"[--colour_list] option requires a filename [file of filenames, one for each colour. Optional second column for sample names]");


	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    //special case. If optarg starts with [, then that's a flag to say "load this colourlist, and for each colour, load the UNION of these binaries (for shared kmers, dont increment covg)
	    char temp1[MAX_FILENAME_LEN];
	    temp1[0]='\0';
	    strcpy(temp1, optarg);
	    if (temp1[0]=='[')
	      {
		if (NUMBER_OF_COLOURS>1)
		  {
		    errx(1, "Sorry - I won't let you load a union of binaries when there is >1 colour. Gets too messy when you start loading colour 2, and find kmers in colour 1\n");
		  }
		cmdline_ptr->for_each_colour_load_union_of_binaries=true;
		printf("You have enabled a Cortex secret option! Congratulations.\n");
		printf("Since you put a [ before the colourlist name, Cortex will load the UNION of these binaries\n");
		printf("(Normally, it cumulates binaries, adding up coverages)\n");

	      }
	    char delims[] = "\t";
	    char* col_list = strtok(temp1, delims);
	    if (col_list==NULL)
	      {
		errx(1,"[--colour_list] option requires a filename [file of filenames, one for each colour, with optional second column for sample id's]");
	      }
	      
	    //check - is this a one-column file, just of filelists, or a 2-column file (col1=filelist, col2=sample name)
	    boolean colour_list_has_samples = check_if_colourlist_contains_samplenames(col_list);
	    int num_cols_in_input_list;
	    if (colour_list_has_samples==false)
	      {
		//how many colours in this filelist?
		num_cols_in_input_list = load_paths_from_filelist(col_list, NULL);
	      }
	    else //colour list has two columns and includes samples
	      {
		num_cols_in_input_list=get_number_of_files_and_check_existence_and_get_samplenames_from_col_list(col_list, cmdline_ptr); 
		cmdline_ptr->loaded_sample_names=true;//can be used to test if colour_list is 1 or 2 column
	      }
	    if (num_cols_in_input_list==0)
	      {
		errx(1, "[--colour_list] filename %s contains nothing", col_list);
	      }
	    cmdline_ptr ->num_colours_in_input_colour_list=num_cols_in_input_list;
	    strcpy(cmdline_ptr->colour_list,col_list);
	    printf("colour list is set to %s\n", cmdline_ptr->colour_list);
	  }
	else
	  {
	    errx(1,"[--colour_list] filename too long [%s]",optarg);
	  }
	
	if ( (cmdline_ptr->for_each_colour_load_union_of_binaries==false) && (access(optarg,R_OK)==-1))
	  {
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
	  errx(1,"[--se_list] option requires a filename (listing single_ended fasta/q)\nOptionally, the first line of the list is allowed to have (aftert the filename) a tab, and then a sample-identifier string, which will be loaded up, and stored in any binary graph you dump.");
	
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
	if (atoi(optarg)<0)
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
	if (atoi(optarg)<0)
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


    case 'o': //--sample_id - only to be used when loading sequence data
      {
	cmdline_ptr->entered_sampleid_as_cmdline_arg = true;
	strcpy(cmdline_ptr->colour_sample_ids[0], optarg);
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


    case 't': // --gt - genotype a list of sites, given in Cortex output format (ie as if called by Cortex)
      {

	if (optarg==NULL)
	  errx(1,"[--gt] option requires an input filename (should be a file of Cortex calls, just 5p, branches and 3p, no colour-coverage output), an output filename, and either BC or PD to specify which caller was used.\n");
	
	if (strlen(optarg)<2*MAX_FILENAME_LEN)
	  {
	    char msg[300];
	    int err= parse_arguments_for_genotyping(cmdline_ptr, optarg, msg);
	    if (err==0)
	      {
		cmdline_ptr->do_genotyping_of_file_of_sites=true;
	      }
	    else
	      {
		die("[--gt] error - %s\n", msg);
	      }
	  }
	else
	  {
	    die("[--gt] argument too long [%s]",optarg);
	  }
	break; 
      }
    case 'u': // --estim_e_with_snps
              //use SNP alleles known not to be present in each colour
              // to estimate sequencing error rates
      {
	// must be a list of fasta files (one per colour)  of SNP alleles which we expect to have zero coverage
	if (access(optarg,R_OK)==-1)
	  {
	    errx(1,"[--estim_e_with_snps] filename [%s] cannot be accessed",optarg);
	  }
	else
	  {
	    strcpy(cmdline_ptr->colourlist_snp_alleles,optarg);
            cmdline_ptr->use_snp_alleles_to_estim_seq_err_rate=true;
	  }
      	break;        
      } 
     

    case 'v': //high_diff - 
      {
	if (optarg==NULL)
	  {
	    errx(1,"[--high_diff requires a decimal argument between 0 and 1 - the minimum difference in allele balance required between two pools\n");
	  }
	cmdline_ptr->min_allele_balance_diff = atof(optarg);
	if ( (cmdline_ptr->min_allele_balance_diff<=0)
	     ||
	     (cmdline_ptr->min_allele_balance_diff>=1) )
	  {
	    errx(1,"[--high_diff requires a decimal argument between 0 and 1 - the minimum difference in allele balance required between two pools\n");
	  }
	cmdline_ptr->high_diff=true;
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
	    die("[--remove_low_coverage_kmers] option requires int argument [node coverage cut off]");
	  }
	cmdline_ptr->node_coverage_threshold = atoi(optarg);
	cmdline_ptr->remove_low_coverage_nodes = true;
	
	if (cmdline_ptr->node_coverage_threshold <= 0)
	  {
	    die("[--remove_low_coverage_kmers] option requires int argument bigger than 0");
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
            die("[--load_colours_only_where_overlap_clean_colour] option requires int argument [which is the clean/preferred colour]");
          }
	cmdline_ptr->load_colours_only_where_overlap_clean_colour=true;
	cmdline_ptr->clean_colour=atoi(optarg);

        if (cmdline_ptr->clean_colour<0)
          {
            die("[--load_colour_list_only_where_overlap_clean_colour] option requires int argument bigger than 0");
          }
	else if  (cmdline_ptr->clean_colour>=NUMBER_OF_COLOURS)
          {
            die("[--load_colour_list_only_where_overlap_clean_colour] option requires int argument <=  max colour limit");
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
    case 'K'://--pan_genome_matrix
	{
	  cmdline_ptr->get_pan_genome_matrix=true;
	  if (strlen(optarg)<MAX_FILENAME_LEN)
	    {
	      strcpy(cmdline_ptr->pan_genome_genes_fasta,optarg);
	    }
	  else
	    {
	      errx(1,"[--pan_genome_matrix] filename too long [%s]",optarg);
	    }
	  if (access(cmdline_ptr->pan_genome_genes_fasta,F_OK)!=0)
	    {
	      errx(1,"[--pan_genome_matrix]   - Cannot open this fasta  [%s] . Abort.\n", 
		   cmdline_ptr->pan_genome_genes_fasta);
	    }
	  
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
    case 'N'://--estimate_copy_num takes input,output
      {
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    //   if (access(optarg,R_OK)==-1)
	    //{
	    //	errx(1,"[--estimate_copy_num] filename [%s] cannot be accessed",optarg);
	    //}
	    char* filename = NULL;
	    char delims[] = ",";
	    char temp[MAX_FILENAME_LEN*2];
	    temp[0]='\0';
	    strcpy(temp, optarg);
	    filename = strtok(temp, delims);
	    if (filename==NULL)
	      {
		errx(1,"[--estimate_copy_num] option requires an input fasta filename, a comma, then an output filename\n");
	      }
	    else if (access(filename,F_OK)!=0){
	      errx(1,"[--estimate_copy_num]   - Cannot open this fasta  [%s] . Abort.\n", filename);
	    }

	    cmdline_ptr->estimate_copy_num = true;
	    strcpy(cmdline_ptr->copy_num_fasta,filename);

	    char* filename2=NULL;
	    filename2 = strtok(NULL, delims);
	    if (filename2==NULL)
	      {
		errx(1,"[--estimate_copy_num] option requires an input fasta filename, a comma, then an output filename - you haven't entered an output filename\n");
	      }
	    else if (access(filename2,F_OK)==0)
	      {
		errx(1,"[--estimate_copy_num] output filename [%s] already exists - will not overwrite. Abort.\n",filename2);
	      }
	    else
	      {
		strcpy(cmdline_ptr->copy_num_output, filename2);
		printf("Will dump copy number estimates in file %s\n", cmdline_ptr->copy_num_output);
	      }
	  }
	else
	  {
	    errx(1,"[--estimate_copy_num] filename too long [%s]",optarg);
	  }
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
	    cmdline_ptr->manually_entered_seq_error_rate = (long double) strtod(optarg,NULL);
	  }
	else if (access(optarg,R_OK)==0)
	  {
	    cmdline_ptr->manually_override_error_rate=true;
	    strcpy(cmdline_ptr->manually_entered_seq_error_rates_file,optarg);
	  }
	else
	  {
	    errx(1,"[--estimated_error_rate] option requires an argument. This can either be the per base-pair error rate (which is then assumed to be true for all colours), OR a file, containing one estimated error rate for each colour. You have neither entered a numeric argument, nor a file which Cortex can find/open\n");
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
    case 'S':
	{
	  parse_err_correction_args(cmdline_ptr, optarg);
	  cmdline_ptr->do_err_correction=true;
	  break;
	}

    case 'T'://subsample
      {
	if (optarg==NULL)
	  errx(1,"[--subsample] option requires a decimal number between 0 and 1 as argunemt\n");
	
	cmdline_ptr->subsample_propn = atof(optarg);
	cmdline_ptr->subsample=true;
	break;
      }
    case 'U'://print_median_covg_only
      {
	cmdline_ptr->print_median_covg_only=true;
	break ;
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
			      (cmdline_ptr->novelseq_colours_search), &(cmdline_ptr->numcols_novelseq_colours_search),
			      (cmdline_ptr->novelseq_colours_avoid),  &(cmdline_ptr->numcols_novelseq_colours_avoid), 
			      &(cmdline_ptr->novelseq_contig_min_len_bp), &(cmdline_ptr->novelseq_min_percentage_novel),
			      (cmdline_ptr->novelseq_outfile));
	  cmdline_ptr->print_novel_contigs=true;

	  break;
	}
    default:
      {
	die("Unknown option %c", opt);
      }      

    }
    opt = getopt_long(argc, argv, "ha:b:c:d:e:f:g:i:jk:l:m:n:o:p:q:r:s:t:u:v:w:xy:z:A:B:CD:E:F:G:H:I:J:K:L:MN:O:P:Q:R:S:T:UV:", long_options, &longopt_index);
    
  }   
  
  return 0;
}
  

int check_cmdline(CmdLine* cmd_ptr, char* error_string)
{
  if ( (cmd_ptr->print_median_covg_only==true) && (cmd_ptr->print_colour_coverages==true) )
    {
      char tmp[] = "--print_median_covg_only and --print_colour_coverages are mutually exclusive alternatives - choose one.\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }
  if  ((cmd_ptr->high_diff==true) &&  (cmd_ptr->genome_size==0))
    {
      char tmp[] = "If you specify --high_diff, you must also set --genome_size\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }
  if ( (cmd_ptr->subsample==true) && (cmd_ptr->input_seq==false) )
    {
      die("If you specify --subsample, you must be entering sequence data as fasta or fastq or bam\n");
    }

  if ( (cmd_ptr->subsample==true) &&   ( (cmd_ptr->input_colours==true) || (cmd_ptr->input_multicol_bin==true))    )
    {
      die("If you specify --subsample, you must be entering sequence data as fasta or fastq or bam - you are passing in binary files (--multicolour_bin or --colour_list)\n");
    }

  
  if ( (cmd_ptr->do_err_correction==true) &&  (cmd_ptr->quality_score_threshold==0) )
    {
      die("If you specify --err_correct then you must also specify a quality threshold with --quality_score_threshold\n");
    }
  if ( (cmd_ptr->do_err_correction==true) &&  (cmd_ptr->max_read_length==0) )
    {
      die("If you specify --err_correct then you must also specify --max_read_len\n");
    }

  if ( (cmd_ptr->get_pan_genome_matrix==true) && (cmd_ptr->max_read_length==0) )
    {
      char tmp[] = "If you want to get a pangenome matrix, you must specify --max_read_len in your file of genes\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }
  if ( (cmd_ptr->dump_readlen_distrib==true) && (cmd_ptr->max_read_length==0) )
    {
      char tmp[] = "If you want to dump the distribution of (filtered) read lengths, you must specify --max_read_len\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }
  
  if ( (cmd_ptr->for_each_colour_load_union_of_binaries==true) && (cmd_ptr->successively_dump_cleaned_colours==true) )
    {
      char tmp[] = "You can't use successively_dump_cleaned_colours and load a colour_list in \"union\" mode\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }
  if (cmd_ptr->use_snp_alleles_to_estim_seq_err_rate==true)
    {
      char error_msg[2*MAX_FILENAME_LEN];
      if (check_files_in_list_exist_ignoring_ref_line(cmd_ptr->colourlist_snp_alleles, error_msg, cmd_ptr->ref_colour)==false)
	{
	  die("%s\n",error_msg);
	}
    }
  if ( (cmd_ptr->do_genotyping_of_file_of_sites==true) && (cmd_ptr->max_read_length==0) )
    {
      char tmp[] = "If you use --gt, then you must specify --max_read_len (and it must be >= any of the reads (including flanks) in your input file).\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }

  if ( (cmd_ptr->do_genotyping_of_file_of_sites==true) && (cmd_ptr->which_caller_was_used_for_calls_to_be_genotyped==SimplePathDivergenceCaller) 
       &&  (cmd_ptr->ref_colour==-1)  )
    {
      char tmp[] = "If you use --gt, and specify that the calls to genotype were called with the Path Divergence caller, then you must specify --ref_colour.\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }


  if ( (cmd_ptr->do_genotyping_of_file_of_sites==true) && (cmd_ptr->which_caller_was_used_for_calls_to_be_genotyped==SimplePathDivergenceCaller) 
       &&  (cmd_ptr->expt_type!=EachColourADiploidSampleExceptTheRefColour) && (cmd_ptr->expt_type!=EachColourAHaploidSampleExceptTheRefColour) )

    {
      char tmp[] = "If you use --gt, and specify that the calls to genotype were called with the Path Divergence caller, then you must specify --experiment_type to be EachColourADiploidSampleExceptTheRefColour or EachColourAHaploidSampleExceptTheRefColour.\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }



  if ( (cmd_ptr->do_genotyping_of_file_of_sites==true) && (cmd_ptr->max_read_length<cmd_ptr->kmer_size) )
    {
      char tmp[] = "You have specified a max_read_length < kmer_size\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }



  if ( (cmd_ptr->do_genotyping_of_file_of_sites==true) && ( (cmd_ptr->expt_type==Unspecified)|| (cmd_ptr->genome_size==0) )     )
    {
      char tmp[] = "If you specify --gt, you must also set --genome_size and --experiment_type\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }


  if (cmd_ptr->kmer_size==-1)
    {
      char tmp[] = "You must specify kmer_size\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
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
		  die("coding error - this string is too long:\n%s\n", tmp);
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
		  die("coding error - this string is too long:\n%s\n", tmp);
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
		  die("coding error - this string is too long:\n%s\n", tmp);
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
		  die("coding error - this string is too long:\n%s\n", tmp);
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
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
	
      }

  if ( (cmd_ptr->expt_type==EachColourADiploidSampleExceptTheRefColour) && (cmd_ptr->ref_colour==-1) )
    {
      char tmp[] = "If you specify the experiment as EachColourADiploidSampleExceptTheRefColour, then you need to specify the ref colour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if ( (cmd_ptr->expt_type==EachColourAHaploidSampleExceptTheRefColour) && (cmd_ptr->ref_colour==-1) )
    {
      char tmp[] = "If you specify the experiment as EachColourAHaploidSampleExceptTheRefColour, then you need to specify the ref colour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if ( (cmd_ptr->expt_type==EachColourADiploidSample) && (cmd_ptr->ref_colour!=-1) )
    {
      char tmp[] = "You have specified a reference colour, and yet also set experiment type to EachColourADiploidSample. Surely you meant EachColourADiploidSampleExceptTheRefColour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if ( (cmd_ptr->expt_type==EachColourAHaploidSample) && (cmd_ptr->ref_colour!=-1) )
    {
      char tmp[] = "You have specified a reference colour, and yet also set experiment type to EachColourAHaploidSample. Surely you meant EachColourAHaploidSampleExceptTheRefColour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }





  if ( (cmd_ptr->exclude_ref_bubbles==true) && (cmd_ptr->ref_colour==-1) )
    {
      char tmp[] = "If you specify --exclude_ref_bubbles, then you must specify a reference colour with --ref_colour\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if ( (cmd_ptr->input_seq ==true) && ( (cmd_ptr->input_colours==true) || (cmd_ptr->input_multicol_bin==true) ) )
    {
      char tmp[] = "Cannot pass in both binary and sequence (fasta/fastq) data\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  if ( (cmd_ptr->entered_sampleid_as_cmdline_arg==true) && (cmd_ptr->input_seq ==false) )
    {
      char tmp[] = "You can only use --sample_id when loading sequence data (fasta/q), not loading binaries\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }

  if (cmd_ptr->input_colours ==true)
    {
      check_colour_list(cmd_ptr->colour_list, cmd_ptr->kmer_size);
      
      if ( (cmd_ptr->successively_dump_cleaned_colours==false) && (cmd_ptr->num_colours_in_input_colour_list>NUMBER_OF_COLOURS) )
	{
	  die("You are trying to load more colours than you have compiled cortex for. \n"
        "Your colour_list contains %d colours, but cortex_var is compiled for a "
        "maximum of %d\n",
		    cmd_ptr->num_colours_in_input_colour_list, NUMBER_OF_COLOURS);
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
	      die("coding error - this string is too long:\n%s\n", tmp);
	    }
	  strcpy(error_string, tmp);
	  return -1;
    }
  

  if ( (cmd_ptr->successively_dump_cleaned_colours==true) && (cmd_ptr->load_colours_only_where_overlap_clean_colour==false) )
    {
     	  char tmp[]="If you specify --successively_dump_cleaned_colours, you must also specify --load_colours_only_where_overlap_clean_colour\n";
	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      die("coding error - this string is too long:\n%s\n", tmp);
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
	      die("coding error - this string is too long:\n%s\n", tmp);
	    }
	  strcpy(error_string, tmp);
	  return -1;
	}

      GraphInfo* temp_ginfo=graph_info_alloc_and_init();
      BinaryHeaderErrorCode ecode=EValid;
      BinaryHeaderInfo binfo;
      initialise_binary_header_info(&binfo, temp_ginfo);
      query_binary_NEW(fp, &binfo, &ecode,0);
      fclose(fp);
      graph_info_free(temp_ginfo);

      if (cmd_ptr->clean_colour >= binfo.number_of_colours)
	{
	  char tmp[]="Specified clean colour number is > number of colours in the multicolour binary. Consult the manual.\n";
	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      die("coding error - this string is too long:\n%s\n", tmp);
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
	      die("coding error - this string is too long:\n%s\n", tmp);
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
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }
  
  
  if ( (cmd_ptr->using_ref==true) && ((cmd_ptr->remv_low_covg_sups_threshold!=-1)|| (cmd_ptr->remove_low_coverage_nodes==true)) )

    {

      char tmp[]="Error-correction techniques should not be used on a graph including a reference genome.\nWe recommend error-correcting single-colour binaries of sequencing data, and then loading them into a multicolour graph including a reference\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
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
	  die("coding error - this string is too long:\n%s\n", tmp);
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
	  die("coding error - this string is too long:\n%s\n", tmp);
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
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }


  if ( (cmd_ptr->apply_model_selection_at_bubbles==true) && (cmd_ptr->genome_size==0) )
    {
      char tmp[]="If you specify --require_hw, you must also specify --genome_size";
      if (strlen(tmp)>LEN_ERROR_STRING)
        {
          die("coding error - this string is too long:\n%s\n", tmp);
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
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }
  if ( (cmd_ptr->detect_bubbles2==true) && (cmd_ptr->output_detect_bubbles2[0]=='\0') )
    {

      char tmp[]="If you specify --detect_bubbles2, you must also specify an output file with --output_bubbles2";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }
  
  if (  (cmd_ptr->max_read_length>100000) && (cmd_ptr->align_given_list==false) && (cmd_ptr->do_genotyping_of_file_of_sites==false) )
    {
      char tmp[]="You are not allowed to set maximum read length >100000 when loading fasta/q into a graph. (However it IS allowed if you are aligning fasta/q to  preexisting graph, or genotyping a file of calls (which might be very long)). Since you have not specified --align, Cortex has exited.\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }
  

  if (  (cmd_ptr->max_read_length>300000000) && (cmd_ptr->align_given_list==true) )
    {
      char tmp[]="If doing alignment, max read length is 300Mb (for whole chromosomes)\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
      
    }

  /*  DEV - to remove - think this is no longer necessary
  if ( (cmd_ptr->input_seq==true) && (cmd_ptr->max_read_length==0) )
    {
      char tmp[]="Must specify max read-length if inputting fasta/fastq data\n";
      if (strlen(tmp)>LEN_ERROR_STRING)
	{
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }
  */
  //if making path divergence calls, must have >=2 colours, and must specify the ref colour, and give a ref_fasta list
  if (cmd_ptr->make_pd_calls==true)
    {
      if ((cmd_ptr->using_ref==false)||(cmd_ptr->ref_colour==-1) || (strlen(cmd_ptr->ref_chrom_fasta_list)==0)  )
	{
	  char tmp[]="If --path_divergence_caller is specified, then must also specify --ref_colour and --list_ref_fasta";

	  if (strlen(tmp)>LEN_ERROR_STRING)
	    {
	      die("coding error - this string is too long:\n%s\n", tmp);
	    }
	  strcpy(error_string, tmp);
	  return -1;
	}
      
      if (cmd_ptr->path_divergence_caller_output_stub[0]=='\0')
	{
	  char tmp[]="If --path_divergence_caller is specified, then must also specify --path_divergence_caller_output";
          if (strlen(tmp)>LEN_ERROR_STRING)
            {
              die("coding error - this string is too long:\n%s\n", tmp);
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
	  die("coding error - this string is too long:\n%s\n", tmp);
	}
      strcpy(error_string, tmp);
      return -1;
    }

  // check multicolour binary  

  if (cmd_ptr->input_multicol_bin==true)
    {
      FILE* fp = fopen(cmd_ptr->multicolour_bin, "r");
      if (fp==NULL)
	{
	  char tmp[LEN_ERROR_STRING];
	  sprintf(tmp,"Unable to open multicolour bin %s for initial sanity checks\n", cmd_ptr->multicolour_bin) ;
	  strcpy(error_string, tmp);
          return -1;
	}


      
      GraphInfo* ginfo=graph_info_alloc_and_init();
      BinaryHeaderErrorCode ecode=EValid;      
      BinaryHeaderInfo binfo;
      initialise_binary_header_info(&binfo, ginfo);
      boolean is_multicol_bin_ok = check_binary_signature_NEW(fp, cmd_ptr->kmer_size, &binfo, &ecode,0);
      fclose(fp);
      graph_info_free(ginfo);

      char tmp[LEN_ERROR_STRING];

      
       if (is_multicol_bin_ok==false)
	{
	  sprintf(tmp,"This binary %s is not compatible with the current de Bruijn graph parameters, error code %d\n", 
		  cmd_ptr->multicolour_bin, ecode);
          strcpy(error_string, tmp);
          return -1;
	  
	}

      if (binfo.number_of_colours<=0)
	{
	  sprintf(tmp, "Corrupt binary %s - signature claims to have <=0 colours within\n", cmd_ptr->multicolour_bin);
          strcpy(error_string, tmp);
          return -1;
	}
      else if (binfo.number_of_colours >NUMBER_OF_COLOURS)
	{
	  sprintf(tmp, "Multicolour binary %s contains %d colours, but cortex_var is compiled to support a maximum of %d colours\n", 
		  cmd_ptr->multicolour_bin, binfo.number_of_colours, NUMBER_OF_COLOURS);
          strcpy(error_string, tmp);
          return -1;


	}
      else if ( (cmd_ptr->successively_dump_cleaned_colours==false) && (binfo.number_of_colours+cmd_ptr->num_colours_in_input_colour_list > NUMBER_OF_COLOURS) )
	{
	  sprintf(tmp, "Between %s (containing %d colours) and %s (containing %d colours), you have exceeded the compile-time limit on colours, %d\n",
		 cmd_ptr->multicolour_bin, binfo.number_of_colours, cmd_ptr->colour_list, cmd_ptr->num_colours_in_input_colour_list, NUMBER_OF_COLOURS);
          strcpy(error_string, tmp);
          return -1;

	}

       if (is_multicol_bin_ok==false)
	{
	  sprintf(tmp,"This binary %s is not compatible with the current de Bruijn graph parameters, error code %d\n", 
		  cmd_ptr->multicolour_bin, ecode);
          strcpy(error_string, tmp);
          return -1;
	  
	}
      
      cmd_ptr->num_colours_in_multicol_bin=binfo.number_of_colours;
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
  /*
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
  */





  return 0;
}


void parse_cmdline(CmdLine* cmd_line, int argc, char* argv[], int unit_size) 
{	
  int i;
  printf("Command: ");
  for(i=0;i<argc;i++){
    printf("%s ",argv[i]);
  }
  printf("\n");
  // printf("Unit size:%i\n",unit_size);

  //CmdLine cmd_line;
  default_opts(cmd_line);

  char error_string[LEN_ERROR_STRING]="";
  int err = parse_cmdline_inner_loop(argc, argv, unit_size, cmd_line, error_string);


  if (err==-1) 
    {
      die("Error in cmd-line input: %s\n", error_string);
    }


  char error_string2[LEN_ERROR_STRING]="";
    
  //  extra_cmd_line_settings(&cmd_line);//some cmd_line settings depend on other options, and we need to wait til all parsing is done

  err = check_cmdline(cmd_line, error_string2);

  if (err == -1)
    {
      die("Error in cmd-line input: %s\n", error_string2);
    }
  
  
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
	  die("Passing array limit of %d int get_numbers_from_comma_sep_list which is less than the number of colours\n", max_len_return_list);
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
	  die("Passing array limit of %d int get_numbers_from_comma_sep_list which is less than the number of colours\n", max_len_return_list);
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
	  die("cmdline error - Failed to parse out the colours which we should genotype in --genotype_site");
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
      die("Badly formatted/chosen arguments for --print_novel_contigs. Consult the manual\n");
    }



  char* commaseplist2 = strtok(NULL, delims );
  if (commaseplist2==NULL)
    {
      errx(1,"[--print_novel_conrigs] option requires an argument of the form a,b/c,d/x/y/filename  you do not appear to have anything after the first /\n");
    }

  *num_cols_in_avoid_list = get_numbers_from_comma_sep_list(commaseplist2, array_colours_to_avoid, NUMBER_OF_COLOURS);
  if (*num_cols_in_avoid_list==-1)
    {
      die("Badly formatted/chosen arguments for --print_novel_contigs. Consult the manual\n");
    }
  else if (*num_cols_in_avoid_list>=NUMBER_OF_COLOURS)
    {
      die("[--print_novel_contigs] option requires an argument of the form a,b,../c,d,../x/y/filename. You have listed ALL the colours in the graph in c,d (either by enumerating them all, or by entering -1). this is meaningless - the idea is to look for stuff in the graph that is NOT in colours c,d,etc\n");
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
	  die("[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, where y is the minimum percentage of the kmers in a contig which we require to NOT be in any colours c,d,.... ie if we want to be very stringent, we set y=100, and demand ALL kmers be novel (ie outside colours c,d..). You have chosen a value for y which is <0, which we do not allow\n");
	}
      else if (tmp>100)
	{
	  die("[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, where y is the minimum percentage of the kmers in a contig which we require to NOT be in any colours c,d,.... ie if we want to be very stringent, we set y=100, and demand ALL kmers be novel (ie outside colours c,d..). You have chosen a value for y which is >100, which is meaningless\n");
	  
	}
    }
  else
    {
	  die("[--print_novel_contigs] option requires an argument of the form a,b/c,d/x/y/filename, where y is the minimum percentage of the kmers in a contig which we require to NOT be in any colours c,d,.... ie if we want to be very stringent, we set y=100, and demand ALL kmers be novel (ie outside colours c,d..). You have chosen a NON_NUMERIC value for y\n");
      
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
		die("Unexpected argument to parse_colourinfo_argument, which_detect_bubbles is not 1,2 or 3");
	      }
	    
	  }
  else
    {
      die("Coding error. Have colour info argument longer than %d, which should not happen - but should have been caught before now\n", MAX_LEN_DETECT_BUB_COLOURINFO);
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
      die("Coding error. Have colour info argument longer than %d, which should not happen - but should have been caught before now\n", MAX_LEN_DETECT_BUB_COLOURINFO);
    }

  return 0;
}


int parse_arguments_for_genotyping(CmdLine* cmdline, char* argmt, char* msg)
{

  //argument should be of form inputfilename,outputfilename,<CALLER>  where CALLER is "BC" or "PD".
  msg[0]='\0';
  char* filename1=NULL;
  char* filename2=NULL;
  char* caller=NULL;
  char delims[] = ",";
  char temp1[2*MAX_FILENAME_LEN];
  temp1[0]='\0';
  strcpy(temp1, argmt);
  filename1 = strtok(temp1, delims );
  if (filename1==NULL)
    {
      strcat(msg, "[--gt] option requires two filenames, plus either BC or PD, comma separated. The filenames are to be an input (file of cortex calls) and output filename.");
      return 1;
    }
  else if (access(filename1,R_OK)==-1)
    {
      strcat(msg, "[--gt] option requires two filenames, plus either BC or PD, comma separated. The filenames are to be an input (file of cortex calls) and output filename. In this case the input filename does not exist - unable to open it");
      return 1;

    }
  strcpy(cmdline->file_of_calls_to_be_genotyped,filename1);
  filename2 = strtok( NULL, delims );
  if (filename2==NULL)
    {
      strcat(msg, "[--gt] option requires two filenames, plus either BC or PD, comma separated. The filenames are to be an input (file of cortex calls) and output filename. Cannot find the output filename in your cmdline input\n");
      return 1;
    }
  strcpy(cmdline->output_genotyping, filename2);
  caller= strtok( NULL, delims );
  if (caller==NULL)
    {
      strcat(msg, "[--gt] option requires two filenames, plus either BC or PD, comma separated. You have entered the filenames, but not specified which caller was used to make these calls\n");
      return 1;
    }
  if (strcmp(caller, "BC")==0) 
    {
      cmdline->which_caller_was_used_for_calls_to_be_genotyped=BubbleCaller;
    }
  else if (strcmp(caller, "PD")==0 )
    {
      cmdline->which_caller_was_used_for_calls_to_be_genotyped=SimplePathDivergenceCaller;
    }
  else
    {
      strcat(msg, "[-gt] option - third part of argument MUST be either BC or PD.");
      return 1;
    }
  return 0;
  
  
}

/*
//returns true if it did find a sample identifier on the top line of the se/pe list.
//will exit if it finds an identifier, but finds the id is already set (means an identifier was
// set in some other file also) and it does not match. ASSUMPTION is fasta/q only goes into colour 0,
// and all has same sample id.
boolean get_sample_id_from_se_pe_list(char* cmdline_sampleid, char* se_pe_list)
{
  boolean found_an_id=false;
  FILE* fp = fopen(se_pe_list, "r");
  if (fp==NULL)
    {
      die("Unable to open se/pe filelist %s\n", se_pe_list);
    }
  else
    {
      char line[MAX_FILENAME_LENGTH+MAX_LEN_SAMPLE_NAME];
      line[0]='\0';
      if (fgets(line, MAX_FILENAME_LENGTH+MAX_LEN_SAMPLE_NAME, fp) != NULL)
	{
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(line, '\n')) != NULL)
	    {
	      *p = '\0';
	    }
	  
	  //now this line is either just a filename, or a filename and sample_idshould be tab separated with two columns.
	  char* file;
	  char* sample_name;
	  char delims[]="\t";
	  file = strtok(line, delims);
	  if (file==NULL)
	    {
	      //empty file
	    }
	  else if (access(file, R_OK)==-1)
	    {
	      die("Cannot access file %s listed in %s. Cortex thinks this is a filename. Either you are trying to list sample_identifiers in your se_list or pe_list file, (in which case - make sure you have a tab between the filename (first column) and sample-id (second column)), or you have listed a file without the right path, as it does not seem to exist/have the right permissions\n", file, se_pe_list);
	    }
	  else
	    {
	      sample_name = strtok(NULL, delims);
	      if (sample_name !=NULL)
		{
		  if (strcmp(cmdline_sampleid, "undefined")==0)
		    {
		      strcpy(cmdline_sampleid, sample_name);
		      found_an_id=true;
		    }
		  else if (strcmp(cmdline_sampleid, sample_name)!=0)
		    {
		      //you seem to have got a second sample id for the same data/sample
		      die("Your se/pe lists seem to contain TWO sample id's, %s and %s. One identifier only allowed when loading FASTA/Q, as it all goes into colour 0, and is treated like one sample\n", cmdline_sampleid, sample_name);
		    }

		}
	    }
	}
      fclose(fp);
    }
  return found_an_id;
}
*/

//returns number of filenames/colours in this file
int get_number_of_files_and_check_existence_and_get_samplenames_from_col_list(char* colour_list, CmdLine* cmd)
{
  int count=0;

  FILE* fp = fopen(colour_list, "r");
  if (fp==NULL)
    {
      die("Running initial sanity check: Cannot open the colourlist %s. Abort.\n", colour_list);
    }

  char line[2*MAX_FILENAME_LENGTH];
  line[0]='\0';
  
  while (feof(fp)==0)
    {
      if (fgets(line, 2*MAX_FILENAME_LENGTH, fp) != NULL)
	{
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(line, '\n')) != NULL)
	    {
	      *p = '\0';
	    }

	  //now this line should be tab separated with two columns.
	  char* file;
	  char* sample_name;
	  char delims[]="\t";
	  file = strtok(line, delims);
	  if (file==NULL)
	    {
	      die("Format issue in %s. Started off looking like two columns (tab sep), but one of the lines looks bad\n", line);
	    }
	  else if (access(file, R_OK)==-1)
	    {
	      die("Cannot access file %s listed in %s. Cortex thinks this is a filename. If you are trying to list sample_identifiers in your colour_list file, make sure you have a tab between the filename (first column) and sample-id (second column)\n", file, colour_list);
	    }
	  else
	    {
	      sample_name = strtok(NULL, delims);
	      strcpy(cmd->colour_sample_ids[count], sample_name);
	      count++;
	    }

	}
    }
  fclose(fp);

  return count;
}


boolean check_if_colourlist_contains_samplenames(char* filename)
{
  boolean contains_samplenames=false;

  FILE* fp = fopen(filename, "r");
  if (fp==NULL)
    {
      die("Running initial sanity check: Cannot open the colourlist %s. Abort.\n", filename);
    }

  char line[2*MAX_FILENAME_LENGTH];
  line[0]='\0';
  
  while (feof(fp)==0)
    {
      if (fgets(line, 2*MAX_FILENAME_LENGTH, fp) != NULL)
	{
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(line, '\n')) != NULL)
	    {
	      *p = '\0';
	    }

	  //Is this one column, or two?
	  char* file;
	  char* sample_name;
	  char delims[]="\t";
	  file = strtok(line, delims);
	  if (file==NULL)
	    {
	      die("Format issue in %s. Do you have an empty line in your colourlist?\n",line );
	    }
	  else
	    {
	      sample_name = strtok(NULL, delims);
	      if (sample_name !=NULL)
		{
		  contains_samplenames=true;
		}
	    }
	}
    }
  fclose(fp);
  return contains_samplenames;
}



void parse_err_correction_args(CmdLine* cmd, char* arg)
{
  StrBuf* arg_strbuf  = strbuf_create(arg);
  if (arg_strbuf==NULL)
    {
      die("Unable even to parse the cmd-lien without oom. Your server is out of memory\n");
    }
  
  
  char delims[] = ",";
  char* fastqlist = strtok(arg_strbuf->buff, delims );
  if (fastqlist==NULL)
    {
      errx(1,"--err_correct needs a comma-sep list of 5 arguments. A filelist (file, which contains a list of FASTQ). A suffix. An output directory. A 0 (discard a read if it has a low quality base that is uncorrectable) or 1 (print it correcting what you can). Finally, a number. If >0 this is the number of padding bases to add greedily walking through graph - otherwise put zero to do no padding. You don't seem to have given an argument\n");
    }
  
  //check file exists
  if (access(fastqlist, R_OK)==-1)
    {
      die("Cannot access file %s\n", fastqlist);
    }

  char* suffix = strtok(NULL, delims );
  if (suffix==NULL)
    {
      errx(1,"--err_correct needs a comma-sep list of 5 arguments. A filelist (file, which contains a list of FASTQ). A suffix. An output directory. A 0 (discard a read if it has a low quality base that is uncorrectable) or 1 (print it correcting what you can). Finally, a number. If >0 this is the number of padding bases to add greedily walking through graph - otherwise put zero to do no padding. You don't seem to have given a suffix or outdir or 0/1\n");
    }
  char* outdir = strtok(NULL, delims );
  if (outdir==NULL)
    {
            errx(1,"--err_correct needs a comma-sep list of 5 arguments. A filelist (file, which contains a list of FASTQ). A suffix. An output directory. A 0 (discard a read if it has a low quality base that is uncorrectable) or 1 (print it correcting what you can). Finally, a number. If >0 this is the number of padding bases to add greedily walking through graph - otherwise put zero to do no padding.You don't seem to have given an outdir or 0/1\n");
    }

  if (dir_exists(outdir)==false)
    {
      die("You have specified an output directory, %s, which does not exist", outdir);
    }



  char* num_as_str = strtok(NULL, delims );
  if (num_as_str==NULL)
    {
      errx(1,"--err_correct needs a comma-sep list of 5 arguments. A filelist (file, which contains a list of FASTQ). A suffix. An output directory. A 0 (discard a read if it has a low quality base that is uncorrectable) or 1 (print it correcting what you can). Finally, a number. If >0 this is the number of padding bases to add greedily walking through graph - otherwise put zero to do no padding.You don't seem to have givena 0/1 4th arg\n");

    }
  char* pad_as_str = strtok(NULL, delims );
  if (pad_as_str==NULL)
    {
      errx(1,"--err_correct needs a comma-sep list of 5 arguments. A filelist (file, which contains a list of FASTQ). A suffix. An output directory. A 0 (discard a read if it has a low quality base that is uncorrectable) or 1 (print it correcting what you can). Finally, a number. If >0 this is the number of padding bases to add greedily walking through graph - otherwise put zero to do no padding.You don't seem to have given a 5th arg\n");

    }


  //put values into the cmdline object
  strbuf_append_str(cmd->err_correction_filelist, fastqlist);
  strbuf_append_str(cmd->err_correction_outdir, outdir);
  strbuf_add_slash_on_end(cmd->err_correction_outdir);
  strbuf_append_str(cmd->err_correction_suffix, suffix);
  if (strcmp(num_as_str, "0")==0)
    {
      cmd->err_correction_policy=DiscardReadIfLowQualBaseUnCorrectable;
    }
  else if (strcmp(num_as_str, "1")==0)
    {
      cmd->err_correction_policy = DontWorryAboutLowQualBaseUnCorrectable;
    }
  else
    {
      errx(1,"--err_correct needs a comma-sep list of 5 arguments. A filelist (file, which contains a list of FASTQ). A suffix. An output directory. A 0 (discard a read if it has a low quality base that is uncorrectable) or 1 (print it correcting what you can). Finally, a number. If >0 this is the number of padding bases to add greedily walking through graph - otherwise put zero to do no padding.Your 4th arg is neither 0 nor 1\n");

    }


  if (strcmp(pad_as_str, "0")==0)
    {
    }
  else
    {
      int n = atoi(pad_as_str);
      if (n>0)
	{
	  cmd->do_greedy_padding=true;
	  cmd->greedy_pad=n;
	}
      else
	{
	  errx(1,"--err_correct needs a comma-sep list of 5 arguments. A filelist (file, which contains a list of FASTQ). A suffix. An output directory. A 0 (discard a read if it has a low quality base that is uncorrectable) or 1 (print it correcting what you can). Finally, a number. If >0 this is the number of padding bases to add greedily walking through graph - otherwise put zero to do no padding.Your 5th arg is neither 0 nor a positive integer\n");	  
	}
    }
  cmd->do_err_correction=true;

  strbuf_free(arg_strbuf);
}

