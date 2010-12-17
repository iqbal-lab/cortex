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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <cmd_line.h>
#include <binary_kmer.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>

const char* usage=
"\nusage: cortex_con [-h] [--input file_of_files] [--mem_height 14] [--dump_binary bin_output] [--input_format fastq|fasta|binary] [--output_contigs contigs.fa] \n" \
"by M. Caccamo (mario.caccamo@bbsrc.ac.uk) (primary contact for cortex_con), and Z. Iqbal (zam@well.ox.ac.uk)\n" \
"\n" \
"   [-h | --help] = This help screen.\n" \
"   [--input FILENAME] = File of filenames to be processed (start and end read is optional, format <filename>  <start read index>  <end read index> ).\n" \
"   [--kmer_size INT] = Kmer size (default 21), it has to be an odd number.\n" \
"   [--mem_width INT] = Size of hash table buckets (default 100).\n" \
"   [--mem_height INT] = Number of buckets in hash table in bits (default 10, this is a power of 2, ie 2^mem_height).\n" \
"   [--tip_clip INT] = Clips the tips in the graph, the argument defines the max length for the tips.\n" \
"   [--quality_score_threshold INT] = Filter for quality scores in the input file, any k-mer wiht a base wiht quality in the threshold or smaller is not considered (default 0).\n" \
"   [--remove_low_coverage_kmers INT] = Filter for kmers with coverage in the threshold or smaller.\n" \
"   [--dump_binary FILENAME] = Dump binary for graph in file (after applying all specified actions on graph).\n" \
"   [--ouput_supernodes FILENAME] = Fasta file with all the supernodes (after applying all specified actions on graph).\n" \
"   [--ouput_contigs FILENAME] = Fasta file with all the contigs (after applying all specified actions on graph).\n" \
"   [--input_format FORMAT] = File format for input (binary|fasta|fastq).\n" \
"   [--print_coverages] = Print coverages for contigs/supernodes in a different file with _cov suffix.\n" \
"   [--remove_seq_errors] = remove sequence of kmers induced by errors.\n"\
"   [--remove_bubbles] = Removes the bubbles in the graph.\n"\
"   [--max_read_len] = Maximum read length over all input files.\n"\

"\n";


int default_opts(CmdLine * c)
{

  //core parameters
  c->verbose=false;
  c->kmer_size = 21;
  c->bucket_size = 100;
  c->number_of_buckets_bits = 10;
  c->input_file=false;

  //actions on/off
  c->tip_clip=false;
  c->low_coverage_node_clip=false;
  c->dump_binary=false;
  c->print_supernodes_fasta=false;
  c->print_paths_fasta=false;
  c->low_coverage_node_clip=false;
  c->detect_bubbles=false;
  c->detect_dirty_bubbles=false;
  c->detect_forks=false;
  c->print_coverages=false;
  c->input_file_format_known=false;
  c->remove_low_coverage_edges=false;
  c->remove_low_coverage_supernodes=false;
  c->remove_bubbles=false;
  c->check_binary = true;
  c->dump_hash_table = false;
  c->health_check_binary=false;
  c->remove_weak_edges=false;

  //parameters
  c->detect_vars_branch_length=0;
  c->quality_score_threshold=0;
  c->quality_score_offset = 33;
  c->coverage_threshold = 30;
  c->remove_low_coverage_supernodes_threshold=1;
  c->max_length_low_coverage_supernode=0;
  c->min_coverage_threshold_remove_weak_edges=30;
  c->max_read_len=500;
  c->binary_version=0;
  return 1;
}


CmdLine parse_cmdline( int argc, char* argv[], int unit_size) 
{	
  int i;
  printf("command: ");
  for(i=0;i<argc;i++){
    printf("%s ",argv[i]);
  }
  printf("\n");
  printf("Unit size:%i\n",unit_size);

  CmdLine cmd_line;
  default_opts(&cmd_line);

  int opt;
  int longopt_index;
 
  static struct option long_options[] = {
    {"remove_bubbles", no_argument, NULL, 'a'},
    {"mem_width", required_argument, NULL, 'b'},
    {"tip_clip", required_argument, NULL, 'c'},
    {"output_contigs", required_argument, NULL, 'd'},
    {"print_coverages", no_argument, NULL, 'e'},
    {"output_supernodes", required_argument, NULL, 'f'},
    {"dump_hash_table", required_argument, NULL, 'g'},
    {"help", no_argument, NULL, 'h'},
    {"input", required_argument, NULL, 'i'},
    {"health_check",no_argument,NULL,'j'},
    {"kmer_size", required_argument, NULL, 'k'},
    {"quality_score_offset", required_argument, NULL, 'l'},
    {"max_length_low_coverage_supernode", required_argument, NULL, 'm'},
    {"mem_height", required_argument, NULL, 'n'},
    {"dump_binary", required_argument, NULL, 'o'},
    {"remove_weak_edges", required_argument, NULL, 'p'},
    {"quality_score_threshold", required_argument, NULL, 'q'},
    {"no_check", no_argument, NULL, 'r'},
    {"remove_low_coverage_supernodes", required_argument, NULL,'s'},
    {"input_format", required_argument, NULL, 't'},
    {"remove_seq_errors", no_argument, NULL, 'u'}, 
    {"verbose", no_argument, NULL, 'v'},
    {"detect_bubbles", required_argument, NULL, 'w'},
    {"coverage_threshold",required_argument,NULL,'x'},
    {"remove_low_coverage_kmers", required_argument, NULL, 'z'},
    {"max_read_len",required_argument,NULL,'A'},
    {"version_binary",required_argument,NULL,'B'},
    {0, 0, 0, 0}
  };
  
  

  while (( opt = 
	   getopt_long(argc, argv,
		       "ab:c:d:ef:g:hi:jk:l:m:n:o:p:q:rs:t:uvw:x:z:A:B:",
		       long_options, &longopt_index)) > 0){

	       
    //Parse the default options
    switch(opt) {
    case 'a':
      cmd_line.remove_bubbles = true;
      break;
      
    case 'b': //bucket size
      if (optarg==NULL)
	errx(1,"[-b | --bsize] option requires int argument [hash table bucket size]");
      cmd_line.bucket_size = atoi(optarg);
      
      if (cmd_line.bucket_size == 0) //check that -b is not bigger than max_short 
	errx(1,"[-b | --bsize] option requires 'short' argument bigger than 0");
      break;

    case 'c'://clip tip
      if (optarg==NULL)
	errx(1,"[-c | --clip_tip] option requires int argument [max length of tips]");
      cmd_line.tip_length = atoi(optarg);
      
      if (cmd_line.tip_length <= 0)
	errx(1,"[-c | --clip_tip] option requires int argument bigger than 0");
      cmd_line.tip_clip = true;
      break ;

    case 'd': //dump paths from the tree decomposion of the graph
      cmd_line.print_paths_fasta = true;
      if (optarg==NULL)
	errx(1,"[-d/--output_paths] option requires a filename");
     
      if (strlen(optarg)<LENGTH_FILENAME)
	{
	  strcpy(cmd_line.output_paths_fasta_filename,optarg);
	}
      else
	{
	  errx(1,"[-d/--output_paths] filename too long [%s]",optarg);
	}
      
      if (access(optarg,F_OK)==0){
	errx(1,"[-d/--output_paths] filename [%s] exist!",optarg);
      }
      break; 
      
    case 'e'://print coverages
      cmd_line.print_coverages = true;
      break ;
 
    case 'f'://output of supernodes
      cmd_line.print_supernodes_fasta = true;
      if (optarg==NULL)
	errx(1,"[-f/--output_supernodes] option requires a filename");
      
      if (strlen(optarg)<LENGTH_FILENAME)
	{
	  strcpy(cmd_line.output_supernodes_fasta_filename,optarg);
	}
      else
	{
	  errx(1,"[-f/--output_supernodes] filename too long [%s]",optarg);
	}
      
      if (access(optarg,F_OK)==0){
	errx(1,"[-f/--output_supernodes] filename [%s] exist!",optarg);
      }
      break; 
      
    case 'g'://dump hash table
      cmd_line.dump_hash_table = true;
      if (optarg==NULL)
	errx(1,"[-g/--dump_hash_table] option requires a filename");
      
      if (strlen(optarg)<LENGTH_FILENAME)
	{
	  strcpy(cmd_line.output_hash_filename,optarg);
	}
      else
	{
	  errx(1,"[-g/--dump_hash_table] filename too long [%s]",optarg);
	}
      
      if (access(optarg,F_OK)==0){
	errx(1,"[-g/--dump_hash_table] filename [%s] exist!",optarg);
      }
      break; 
      
    case 'h':
      printf("***********************\n");
      printf("Cortex version %d.%d.%d\n", VERSION, SUBVERSION, SUBSUBVERSION);
      printf("***********************\n");
      printf("%s",usage);
      exit(0);
      break;
     
    case 'i'://file of filenames
      if (optarg==NULL)
	errx(1,"[-i | --input] option requires a filename [file of filenames]");
      
      if (strlen(optarg)<LENGTH_FILENAME)
	{
	  strcpy(cmd_line.input_filename,optarg);
	}
      else
	{
	  errx(1,"[-i | --input] filename too long [%s]",optarg);
	}
      
      if (access(optarg,R_OK)==-1){
	errx(1,"[-i | --input] filename [%s] cannot be accessed",optarg);
      }
      cmd_line.input_file = true;
      break;
      
    case 'j'://print coverages
      cmd_line.health_check_binary = true;
      break ;
 
    case 'k': //kmer size
      if (optarg==NULL)
	errx(1,"[-k | --kmer_size] option requires int argument [kmer size]");	    
      cmd_line.kmer_size = atoi(optarg);
      
      if (cmd_line.kmer_size == 0)
	errx(1,"[-k | --kmer_size] option requires int argument bigger than 0");
      
      cmd_line.detect_vars_branch_length=cmd_line.kmer_size;
      break;

    case 'l': //quality_score_offset
      if (optarg==NULL)
	errx(1,"[-l | --quality_score_offset] option requires int argument");	    
      cmd_line.quality_score_offset = atoi(optarg);
      
      if (cmd_line.quality_score_offset == 0)
	errx(1,"[-l | --quality_score_offset] option requires int argument bigger than 0");
      break;

    case 'm': //max length of supernode to remove when is below low coverage threshold
      if (optarg==NULL)
	errx(1,"[-m | --max_length_low_coverage_supernode] option requires int argument");	    
      cmd_line.max_length_low_coverage_supernode = atoi(optarg);
      
      if (cmd_line.max_length_low_coverage_supernode == 0)
	errx(1,"[-m | --max_length_low_coverage_supernode] option requires int argument bigger than 0");
      break;

    case 'n': //number of buckets
      if (optarg==NULL)
	errx(1,"[-n | --hsize] option requires int argument [hash table number of buckets in bits]");
      cmd_line.number_of_buckets_bits = atoi(optarg);      
      break;  

    case 'o'://output of binary ctx
      cmd_line.dump_binary=true;
      if (optarg==NULL)
	errx(1,"[-o/--output_ctx] option requires a filename");
      
      if (strlen(optarg)<LENGTH_FILENAME){
	  strcpy(cmd_line.output_ctx_filename,optarg);
	}
      else{
	  errx(1,"[-o/--output_ctx] filename too long [%s]",optarg);
	}
      
      if (access(optarg,F_OK)==0){
	errx(1,"[-o/--output_ctx] filename [%s] exists!",optarg);
      }
      break; 
      
    case 'p': //remove weak edges
      if (optarg==NULL)
	errx(1,"[-p | --remove_weak_edges] option requires int argument [node coverage]");
      cmd_line.min_coverage_threshold_remove_weak_edges= atoi(optarg);
      
      if (cmd_line.min_coverage_threshold_remove_weak_edges <= 0)
	errx(1,"[-p | --remove_weak_edges] option requires int argument bigger than 0");

      cmd_line.remove_weak_edges=true;
      break;    
    
    case 'q'://quality threshold
      if (optarg==NULL)
	errx(1,"[-q | --quality_score_threshold] option requires int argument [quality score threshold]");
      cmd_line.quality_score_threshold = atoi(optarg);
      
      cmd_line.input_file_format = FASTQ;
      
      if (cmd_line.quality_score_threshold == 0)
	errx(1,"[-q | --quality_score_threshold] option requires int argument bigger than 0");
      break;    

    case 'r': //no check for binaries (for old binaries deprecated)
      cmd_line.check_binary=false;
      break;
      
    case 's': //remove_low_coverage_supernodes
      cmd_line.remove_low_coverage_supernodes=true;
      if (optarg==NULL)
	errx(1,"[-s | --remove_low_coverage_supernodes] option requires int argument [coverage threshold]");
      cmd_line.remove_low_coverage_supernodes_threshold= atoi(optarg);
      
      if (cmd_line.remove_low_coverage_supernodes_threshold <= 0)
	errx(1,"[-s | --remove_low_coverage_supernodes] option requires intint argument bigger than 0");
      break;    
      
      cmd_line.remove_low_coverage_supernodes=true;
      cmd_line.remove_low_coverage_supernodes_threshold=1;
      break;  
    
    case 't':
      if (optarg == NULL){
	errx(1,
	     "[-t  | --input_format binary|fasta|fastq ] option requires a file type");
	exit(-1);
      }

      cmd_line.input_file_format_known=true;
      
      if(strcmp(optarg, "binary")==0){
	cmd_line.input_file_format = CTX;
      }else if(strcmp(optarg, "fasta")==0){
	cmd_line.input_file_format = FASTA;
      }else if(strcmp(optarg, "fastq")==0){
	cmd_line.input_file_format = FASTQ;
       }else if(strcmp(optarg, "hash")==0){
	cmd_line.input_file_format = HASH;
       }else  {
	fprintf(stderr,
		"[-t  | --input_format  binary|fasta|fastq ] invalid option %s ]", optarg);
	exit(-1);
      }

      break;

    case 'u': //remove_seq_errors --> this is consistent with an option on cortex_var -- it removes low coverate supernodes fixing the length tp k+1 and the coverage at 1
      cmd_line.remove_low_coverage_supernodes=true;
      cmd_line.remove_low_coverage_supernodes_threshold=1;
      break;  
      
    case 'v':
      cmd_line.verbose = true;
      break ;
 
   case 'w':
     cmd_line.detect_bubbles = true;
     break ;
     
    case 'x': //coverage threshold use as an indication of unique regions 
      if (optarg==NULL)
	errx(1,"[-x | --coverage_threshold] option requires int argument");
      cmd_line.coverage_threshold = atoi(optarg); 
      
      if (cmd_line.coverage_threshold == 0)     
	errx(1,"[-x | --coverage_threshold] option requires int argument bigger than 0"); 
      break;
      

    case 'z'://node coverage threshold
    
      if (optarg==NULL)
	errx(1,"[-z | --remove_low_coverage_kmers] option requires int argument [node coverage cut off]");

      cmd_line.node_coverage_threshold = atoi(optarg);
      cmd_line.low_coverage_node_clip = true;
      
      if (cmd_line.node_coverage_threshold == 0)
	errx(1,"[-z | --remove_low_coverage_kmers] option requires int argument bigger than 0");
      break;    
      
    case 'A'://max read length --> use by the parser to read fasta/fastq files. 
 
      if (optarg==NULL)
	errx(1,"[-A | --max_read_len] option requires int argument [maximum read length in input]");

      cmd_line.max_read_len = atoi(optarg);
      
      if (cmd_line.max_read_len <= 0)
	errx(1,"[-A | --max_read_len] option requires int argument bigger than 0");
      break;    

    case 'B'://version 
    
    if (optarg==NULL)
	errx(1,"[-B | --binary_version] option requires int argument [version of binary format]");
      
      cmd_line.binary_version = atoi(optarg);
      
      if (cmd_line.binary_version <= 0)
	errx(1,"[-B | --binary_version] option requires int argument bigger than 0");
      break;    
      

    }
  }
  
  //check input file
  if (! cmd_line.input_file){
    errx(1,"input file required [option -i | --input_file]");
  }

  if (cmd_line.number_of_buckets_bits<0){
    printf("args error -b %i -n %i\n",cmd_line.bucket_size,cmd_line.number_of_buckets_bits);
    errx(1,"memory configuration erorr - revise options -b -h and/or -m");
  }

  //check if input file format is known
  if (cmd_line.input_file_format_known == false)
    {
      errx(1,"file format not defined [--input_format fasta | --input_format fastq | --input_format binary]");
    }

  //check kmer_size is odd
  if (cmd_line.kmer_size % 2 == 0){
    errx(1,"[-k/--kmer_size] is even [%i]!",cmd_line.kmer_size);
  }

  //fix length of max length of low coverage supernodes
  if (cmd_line.remove_low_coverage_supernodes==true &&    cmd_line.max_length_low_coverage_supernode==0){
    cmd_line.max_length_low_coverage_supernode=cmd_line.kmer_size+1;
  }

  if (cmd_line.binary_version==0){
    cmd_line.binary_version=BINVERSION;
  }
  return cmd_line;
}



