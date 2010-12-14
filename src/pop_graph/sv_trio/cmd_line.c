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
#include <pop_globals.h>

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
  check_binary_signature(fp, kmer_size, BINVERSION, &num_cols);
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
"Cortex, multicoloured target, cortex_var by Z. Iqbal (primary contact for cortex_var: zam@well.ox.ac.uk) and M. Caccamo (mario.caccamo@bbsrc.ac.uk)\n" \
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
"   [--cut_homopolymers INT] \t\t\t\t\t=\t Breaks reads at homopolymers of length > this threshold. (New read starts after homopolymer)\n" \
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
"   [--output_contigs FILENAME] \t\t\t\t\t=\t Dump a fasta file of all the supernodes (after applying all specified actions on graph).\n" \
  // -r
"   [--detect_bubbles1 COMMA_SEP_COLOURS/COMMA_SEP_COLOURS] \t=\t Output to standard output all the bubbles in the graph where the two branches lie in the specified colours\n\t\t\t\t\t\t\t\t\t (after applying all specified actions on graph).\n\t\t\t\t\t\t\t\t\t Typical use would be --detect_bubbles1 1/1 to find hets in colour 1,\n\t\t\t\t\t\t\t\t\t or --detect_bubbles1 0/1 to find homozygous non-reference bubbles where one branch is in colour 0 (and not colour1)\n\t\t\t\t\t\t\t\t\t and the other branch is in colour1 (but not colour 0).\n\t\t\t\t\t\t\t\t\t However, one can do more complex things:\n\t\t\t\t\t\t\t\t\t e.g.  --detect_het_bubbles 1,2,3/4,5,6 to find bubbles where one branch is in 1,2 or 3 (and not 4,5 or 6)\n\t\t\t\t\t\t\t\t\t and the other branch in colour 4,5 or 6 (but not 1,2, or 3).\n" \
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





int default_opts(CmdLine * c)
{

  c->kmer_size = 21;
  c->bucket_size = 100;
  c->number_of_buckets_bits = 10;
  c->ref_colour=-1;
  c->homopolymer_limit=-1;
  c->quality_score_threshold=0;
  //c->node_coverage_threshold=0;
  c->quality_score_offset = 33;//standard fastq, not illumina v-whatever fastq  
  c->max_read_length = 0;
  c->max_var_len = 10000;
  c->num_colours_in_detect_bubbles1_first_colour_list=0;
  initialise_int_list(c->detect_bubbles1_first_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);
  c->num_colours_in_detect_bubbles1_second_colour_list=0;
  initialise_int_list(c->detect_bubbles1_second_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);
  
  c->num_colours_in_detect_bubbles2_first_colour_list=0;
  initialise_int_list(c->detect_bubbles2_first_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);  
  c->num_colours_in_detect_bubbles2_second_colour_list=0;
  initialise_int_list(c->detect_bubbles2_second_colour_list, MAX_COLOURS_ALLOWED_TO_MERGE);

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
  set_string_to_null(c->ref_chrom_fasta_list,MAX_FILENAME_LEN);
  set_string_to_null(c->config,MAX_FILENAME_LEN);


  //booleans
  c->cut_homopolymers = false;
  c->remove_pcr_dups=false;
  //c->clip_tips=false;
  c->remove_seq_errors=false;
  c->print_colour_coverages=false;
  c->dump_binary=false;
  c->print_contig_fasta=false;
  //c->remove_low_coverage_nodes=false;
  c->detect_bubbles1=false;
  c->detect_bubbles2=false;
  c->make_pd_calls=false;
  c->using_ref=false;
  c->seq_file_format_known=false;
  c->input_colours=false;
  c->input_multicol_bin = false;
  c->input_seq = false;
  c->format_of_input_seq=UNSPECIFIED;


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
    {"output_contigs",required_argument,NULL,'q'},
    {"detect_bubbles1",required_argument, NULL, 'r'},
    {"output_bubbles1",required_argument, NULL, 's'},    
    {"detect_bubbles2",required_argument, NULL, 't'},
    {"output_bubbles2",required_argument, NULL, 'u'},    
    {"format",required_argument,NULL,'v'},
    {"max_read_len",required_argument,NULL,'w'},
    {"print_colour_coverages",no_argument,NULL,'x'},
    {"max_var_len",required_argument,NULL,'y'},
    {"list_ref_fasta",required_argument,NULL,'z'},
    {0,0,0,0}	
  };
  

  //do not change this! Only really matters for testing, but getopt_long uses 
  //variables which are not local to this function, so when we run tests, and call this function
  // repeatedly, those variable values are carried across. The following line resets this.
  optind=1;
  
 
  opt = getopt_long(argc, argv, "ha:b:c:d:e:f:g:i:jk:l:m:n:op:q:r:s:t:u:v:w:xy:z:", long_options, &longopt_index);

  while ((opt) > 0) {
	       
    //Parse the default options
    switch(opt) {

    case 'h':
      {
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
	    else if (num_cols_in_input_list>NUMBER_OF_COLOURS)
	      {
		errx(1, "[--colour_list] filename %s contains %d files, but cortex_ar is currently compiled to suppoer a maximum of %d colours\n", 
		     optarg, num_cols_in_input_list, NUMBER_OF_COLOURS);
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
	  errx(1,"[-b | --multicolour_bin] option requires a filename (of a multicolour binary)");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->multicolour_bin,optarg);
	  }
	else
	  {
	    errx(1,"[-a | --multicolour_bin] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,R_OK)==-1){
	  errx(1,"[-a | --multicolour_bin] filename [%s] cannot be accessed",optarg);
	}
	cmdline_ptr->input_multicol_bin = true;
	//we set num_colours_in_multicol_bin later on - need to open the file and check signature, and cant do that til we know what kmer, and cant be sure in this case 
	//that we have already parsed the --kmer_size case already
	break; 
      }





    case 'c': //se_list
      {
	if (optarg==NULL)
	  errx(1,"[-c | --se_list] option requires a filename (listing single_ended fasta/q)");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->se_list,optarg);
	  }
	else
	  {
	    errx(1,"[-c | --se_list] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,R_OK)==-1){
	  errx(1,"[-c | --se_list] filename [%s] cannot be accessed",optarg);
	}
	cmdline_ptr->input_seq = true;
	break; 
      }


    case 'd': //pe_list, comma separated - two files, listing matched pairs in the same order
      {
	if (optarg==NULL)
	  errx(1,"[-d | --pe_list] option requires two filenames, comma separated, listing paired-end fasta/q files (in matching order)");
	
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
		errx(1,"[-d | --pe_list] option requires two filenames, comma separated, listing paired-end fasta/q files (in matching order)");
	      }
	    strcpy(cmdline_ptr->pe_list_lh_mates,filename);
	    filename = strtok( NULL, delims );
	    if (filename==NULL)
	      {
		errx(1,"[-d | --pe_list] option requires two filenames, comma separated, listing paired-end fasta/q files (in matching order)");
	      }
	    strcpy(cmdline_ptr->pe_list_rh_mates,filename);

	  }
	else
	  {
	    errx(1,"[-d | --pe_list] filename too long [%s]",optarg);
	  }
	
	if (access(cmdline_ptr->pe_list_lh_mates,R_OK)==-1){
	  errx(1,"[-d | --pe_list] filename [%s] cannot be accessed",cmdline_ptr->pe_list_lh_mates);
	}
	if (access(cmdline_ptr->pe_list_rh_mates,R_OK)==-1){
	  errx(1,"[-d | --pe_list] filename [%s] cannot be accessed",cmdline_ptr->pe_list_rh_mates);
	}
	cmdline_ptr->input_seq = true;
	break; 
      }

      case 'e': //kmer size
	{
	  if (optarg==NULL)
	    errx(1,"[-e | --kmer_size] option requires int argument [kmer size]");	    
	  cmdline_ptr->kmer_size = atoi(optarg);
	  
	  if (cmdline_ptr->kmer_size == 0)
	    errx(1,"[-e | --kmer_size] option requires int argument bigger than 0");

	  if (cmdline_ptr->kmer_size % 2 == 0)
	    errx(1,"[-e | --kmer_size] option requires int argument which is not divisible by 2");
	  
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
	  errx(1,"[-g | --mem_height] option requires int argument [hash table,  number of buckets in bits - ie 2^(this number) is the number of buckets]");
	cmdline_ptr->number_of_buckets_bits = atoi(optarg);      
	break;
      }

    case 'i': //ref_colour
      {
	if (optarg==NULL)
	  errx(1,"[-i | --ref_colour] option requires int argument [which colour is the reference - i.e. what position in colour_list (count starts from 0)]");
	if (optarg<0)
	  errx(1,"[-i | --ref_colour] option requires positive argument [which colour is the reference - i.e. what position in colour_list (count starts from 0)]");
	if (atoi(optarg)>NUMBER_OF_COLOURS)
	  errx(1,"[-i | --ref_colour] requires you specify the colour of the reference, and this must be less than the maximum number of colours allowed by Cortex. This maximum number is fixed at compile-time. See the Manual to find how to reset this, and check the number you have entered : %s is correct.", optarg);

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
	  errx(1,"[-k | --cut_homopolymers] option requires positive integer  argument [maximum allowed length of homopolymer in read]");
	if (optarg<0)
	  errx(1,"[-k | --cut_homopolymers] option requires positive integer  argument [maximum allowed length of homopolymer in read]");

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
	  errx(1,"[-m  | --quality_score_threshold] option requires int argument [quality score threshold]");
	cmdline_ptr->quality_score_threshold = atoi(optarg);

	
	if ( (cmdline_ptr->format_of_input_seq != FASTQ) && (cmdline_ptr->format_of_input_seq != UNSPECIFIED) )
	  {
	    errx(1,"[-m  | --quality_score_threshold] has been called, implying the input format is FASTQ, but another command-line input has specified the input format is NOT FASTQ. Inconsistent args");
	  }
	cmdline_ptr->format_of_input_seq = FASTQ;
	
	if (cmdline_ptr->quality_score_threshold == 0)
	  errx(1,"[-m | --quality_score_threshold] option requires int argument bigger than 0");
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
	  errx(1,"[-p/--dump_binary] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->output_binary_filename,optarg);
	    cmdline_ptr->dump_binary=true;
	  }
	else
	  {
	    errx(1,"[-o/--dump_binary] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0){
	  errx(1,"[-o/--dump_binary] filename [%s] exists!",optarg);
	}
	
	break; 
      }


    case 'q': //output supernode contigs
      {
	cmdline_ptr->print_contig_fasta = true;
	if (optarg==NULL)
	  errx(1,"[-q/--output_contigs] option requires a filename");
	
	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->output_supernodes,optarg);
	  }
	else
	  {
	    errx(1,"[-q/--output_contigs] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0)
	  {
	    errx(1,"[-q/--output_contigs] filename [%s] already exists. Exiting, to prevent overwriting.",optarg);
	  }
	if (optarg[0]=='-')
	  {
	    errx(1, "[-q/--output_contigs] requires a filename, but finds this: [%s] starts with -. Have you omitted the filename?\n", optarg);
	  }
	break; 
      }

      
    case 'r': //detect_bubbles1 - will give something of this form 1,2,3/5,1,7 to show the colours it wants to distinguish
      {
	if (optarg==NULL)
	  errx(1,"[-r | --detect_bubbles1] option requires two sets of comma-separated integers, separated by a forward-slash. eg 1/1, or 1,2,3/3,4,5");
	
	int ret = parse_colourinfo_argument(cmdline_ptr, optarg, strlen(optarg), "[-r | --detect_bubbles1] ", 1);
	if (ret==-1)
	  {
	    errx(1, "Problem with  cmd line argument for [-r | --detect_bubbles1]");
	  }
	cmdline_ptr->detect_bubbles1=true;

	break; 
      }
    case 's': //output file for detect_bubbles1
      {
	if (optarg==NULL)
	  errx(1,"[-s | --output_bubbles1] option requires a filename");

	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->output_detect_bubbles1,optarg);
	  }
	else
	  {
	    errx(1,"[-s | --output_bubbles1] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0){
	  errx(1,"[-s | --output_bubbles1] filename [%s] exists!",optarg);
	}
	break;
      }


    case 't': //detect_bubbles2 - will give something of this form 1,2,3/5,1,7 to show the colours it wants to distinguish
      {
	if (optarg==NULL)
	  errx(1,"[-t | --detect_bubbles2] option requires two sets of comma-separated integers, separated by a forward-slash. eg 1/1, or 1,2,3/3,4,5");
	
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
	  errx(1,"[-u | --output_bubbles2] option requires a filename");

	if (strlen(optarg)<MAX_FILENAME_LEN)
	  {
	    strcpy(cmdline_ptr->output_detect_bubbles2,optarg);
	  }
	else
	  {
	    errx(1,"[-u | --output_bubbles2] filename too long [%s]",optarg);
	  }
	
	if (access(optarg,F_OK)==0){
	  errx(1,"[-s | --output_bubbles2] filename [%s] exists!",optarg);
	}
	break;        
      } 
     
    case 'v': //file format - either fasta, fastq or ctx.
      {

	if (optarg==NULL)
	  errx(1,"[-v | --format] option requires argument FASTQ, FASTA or CTX");

	if ( (strcmp(optarg, "FASTA") !=0) && (strcmp(optarg, "FASTQ") !=0) && (strcmp(optarg, "CTX") !=0)  )
	  {
	    errx(1,"[-v | --format] option requires argument FASTQ, FASTA or CTX");
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
	  errx(1,"[-w | --max_read_len] option requires (positive) integer argument");
	if (atoi(optarg)<0)
	  {
	    errx(1,"[-w | --max_read_len] option requires (positive) integer argument");
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
	  errx(1,"[-y | --max_var_len] option requires (positive) integer argument");
	if (atoi(optarg)<0)
	  {
	    errx(1,"[-y | --max_var_len] option requires (positive) integer argument");
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


      
    }
    opt = getopt_long(argc, argv, "ha:b:c:d:e:f:g:i:jk:lm:n:opqr:s:t:u:v:w:xy:z:", long_options, &longopt_index);
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
      if ((cmd_ptr->using_ref==false)||(cmd_ptr->ref_colour==-1) )
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
    }


  // check multicolour binary  

  if (cmd_ptr->input_multicol_bin==true)
    {
      int num_m_cols;
      FILE* fp = fopen(cmd_ptr->multicolour_bin, "r");
      if (fp==NULL)
	{
	  printf("Unable to open multicolour bin %s for checks\n", cmd_ptr->multicolour_bin);
	  exit(1);
	}
      boolean is_multicol_bin_ok = check_binary_signature(fp, cmd_ptr->kmer_size, BINVERSION, &num_m_cols);
      fclose(fp);
      
      if (is_multicol_bin_ok==false)
	{
	  printf("This binary %s is not compatible with the current de Bruijn graph parameters\n", cmd_ptr->multicolour_bin);
	  exit(1);
	}
      
      if (num_m_cols<=0)
	{
	  printf("Corrupt binary %s - signatire claims to have <=0 colours within\n", cmd_ptr->multicolour_bin);
	  exit(1);
	}
      if (num_m_cols>NUMBER_OF_COLOURS)
	{
	  printf("Multicolour binary %s contains %d colours, but cortex_var is compiled to support a maximum of %d colours\n", 
		 cmd_ptr->multicolour_bin, num_m_cols, NUMBER_OF_COLOURS);
	  exit(1);
	}
      else if (num_m_cols+cmd_ptr->num_colours_in_input_colour_list > NUMBER_OF_COLOURS)
	{
	  printf("Between %s (containing %d colours) and %s (containing %d colours), you have exceeded the compile-time limit on colours, %d\n",
		 cmd_ptr->multicolour_bin, num_m_cols, cmd_ptr->colour_list, cmd_ptr->num_colours_in_pd_colour_list, NUMBER_OF_COLOURS);
	  exit(1);
	}
      
      cmd_ptr->num_colours_in_multicol_bin=num_m_cols;
    }



  
  
  //check mem usage is reasonable
  
  //print what  mem usage will be, and output cmd line and what it means

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

  if ( (which_detect_bubbles!=1) && (which_detect_bubbles!=2) )
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
	    else
	      {
		cmd->num_colours_in_detect_bubbles2_first_colour_list=num_left_colours;
		cmd->num_colours_in_detect_bubbles2_second_colour_list=num_right_colours;
		copy_list(left_colours, num_left_colours, cmd->detect_bubbles2_first_colour_list, cmd->num_colours_in_detect_bubbles2_first_colour_list);
		copy_list(right_colours, num_right_colours, cmd->detect_bubbles2_second_colour_list, cmd->num_colours_in_detect_bubbles2_second_colour_list);
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
  

