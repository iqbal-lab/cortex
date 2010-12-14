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
#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <cmd_line.h>
#include <time.h>

void timestamp();

int main(int argc, char **argv){

  printf("Starting Cortex\n");
  timestamp();
  CmdLine cmd_line = parse_cmdline(argc,argv,sizeof(Element));

  FILE *fp;
  char filename[1000];
  int hash_key_bits;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int bucket_size;
  int ctg_length = 200000;
  
   //command line arguments 
  fp = fopen(cmd_line.input_filename, "r");    //open file of file names/hash table dump
  kmer_size        = cmd_line.kmer_size; 
  hash_key_bits    = cmd_line.number_of_buckets_bits;  //number of buckets: 2^hash_key_bits
  bucket_size      = cmd_line.bucket_size;
  DEBUG            = cmd_line.verbose;


  int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
  int max_kmer_size = (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1;
  int min_kmer_size = ((NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1) * sizeof(bitfield_of_64bits) * 4) + 1;

  if (cmd_line.check_binary && cmd_line.input_file_format != HASH && number_of_bitfields != NUMBER_OF_BITFIELDS_IN_BINARY_KMER){
    printf("K-mer %i  is not in current range of kmers [%i - %i]!\n",kmer_size,min_kmer_size,max_kmer_size);
    exit(0); 
  }
  
  
  printf("Maximum k-mer size: %i\n", max_kmer_size);

  if (cmd_line.kmer_size> (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1){
    errx(1,"k-mer size is too big [%i]!",cmd_line.kmer_size);
  }
  printf("Input file of filenames: %s\n",cmd_line.input_filename);

  //Create the de Bruijn graph/hash table
  if (cmd_line.input_file_format == HASH){
    timestamp();
    printf("Loading hash table from file -- parameters obtained from file (ignore command line parameters)\n");
    db_graph =  hash_table_load_from_dump(fp,20);
    
    printf("Kmer size: %d hash_table_size: %qd - bucket size: %d - total size: %qd\n",db_graph->kmer_size,db_graph->number_buckets,db_graph->bucket_size, (long long) db_graph->number_buckets * db_graph->bucket_size);

    printf("table created/loaded: %qd\n",db_graph->number_buckets * db_graph->bucket_size);
    
    
    if (db_graph == NULL){
      timestamp();
      exit(1);
    }
  }
  else{
    printf("Kmer size: %d hash_table_size (%d bits): %d - bucket size: %d - total size: %qd\n",cmd_line.kmer_size,cmd_line.number_of_buckets_bits,1 << cmd_line.number_of_buckets_bits, cmd_line.bucket_size, ((long long) 1<<cmd_line.number_of_buckets_bits)*cmd_line.bucket_size);
    
    db_graph = hash_table_new(cmd_line.number_of_buckets_bits,cmd_line.bucket_size, 20,cmd_line.kmer_size);
    
    if (db_graph == NULL){
      timestamp();
      exit(1);
    }
    
    printf("table created: %d\n",1 << cmd_line.number_of_buckets_bits);
    
    int count_file   = 0;
    long long total_length = 0; //total sequence length
    
    //Go through all the files, loading data into the graph
    
    boolean all_entries_are_unique = false;
    
    timestamp();
    printf("Start loading files\n");
    
    if (cmd_line.input_file_format == CTX){
      all_entries_are_unique = true;
    }
    
    long long bad_reads = 0;
    long long reads = 0;
    long long start=0;
    long long end=0;

    while (!feof(fp)){
      
      int match = fscanf(fp,"%s\t%qd\t%qd\n", filename,&start,&end);
      if (match == 1){
	start=0;
	end=0;
      }
      else{

	if (start<=0 || end<=0 || start>end){
	  printf("inconsistency start [%qd] end[%qd]\n",start,end);
	  exit(1);
	}
	
      }
      
      long long seq_length = 0;
      count_file++;
      if (start!=0){
	printf("reading file %s - start:%qd - end:%qd\n",filename,start,end);
      }
      else{
	printf("reading file %s - full file\n",filename);
      }
      fflush(stdout);

      switch (cmd_line.input_file_format) {
	
      case CTX:
	seq_length += load_binary_from_filename_into_graph(filename,cmd_line.binary_version,db_graph,all_entries_are_unique,cmd_line.check_binary);
	all_entries_are_unique = false;
	break;
	
      case FASTQ:
	seq_length += load_fastq_from_filename_into_graph(filename,&reads,&bad_reads, cmd_line.quality_score_threshold, cmd_line.quality_score_offset, cmd_line.max_read_len, start,end,db_graph);
      break;
      
      case FASTA:
	seq_length += load_fasta_from_filename_into_graph(filename,&reads,&bad_reads, cmd_line.max_read_len,start,end,db_graph);
	break;
	
      case HASH:
	printf("inconsistency\n");
	exit(1);
	break;  
      }
    
      total_length += seq_length;
      timestamp();
      printf("%i kmers: %qd file name:%s reads: %qd bad reads: %qd seq:%qd total seq:%qd\n\n",count_file,hash_table_get_unique_kmers(db_graph),filename,reads,bad_reads,seq_length, total_length);
      hash_table_print_stats(db_graph);
    }
  }

  if (cmd_line.low_coverage_node_clip){
    timestamp();
    printf("remove low coverage nodes (<=%i) \n",cmd_line.node_coverage_threshold);
    fflush(stdout);
    db_graph_remove_low_coverage_nodes(cmd_line.node_coverage_threshold,db_graph);
  }

  if (cmd_line.tip_clip){
    timestamp();
    printf("clip tips\n");
    fflush(stdout);
    printf("%i nodes clipped\n",db_graph_clip_tips(cmd_line.tip_length,db_graph));
  }
  

  if (cmd_line.remove_weak_edges){
    timestamp();
    printf("remove weak edges [threshold %i]\n",cmd_line.min_coverage_threshold_remove_weak_edges);
    fflush(stdout);
    printf("%qd edges removed\n",db_graph_remove_weak_edges(cmd_line.min_coverage_threshold_remove_weak_edges,db_graph));
  }


  if(cmd_line.remove_low_coverage_supernodes){
    timestamp();
    printf("removing low coverage supernodes - threshold %i - max length: %i\n",cmd_line.remove_low_coverage_supernodes_threshold,cmd_line.max_length_low_coverage_supernode);
    fflush(stdout);
    db_graph_clip_low_coverage_supernodes(cmd_line.remove_low_coverage_supernodes_threshold,cmd_line.max_length_low_coverage_supernode,ctg_length,db_graph);
  }
  
  if(cmd_line.remove_bubbles){
    timestamp();
    printf("removing bubbles\n");
    fflush(stdout);
    int kmers_removed_1 = db_graph_remove_bubbles(db_graph->kmer_size*10+1, db_graph);
    printf("%i kmers removed\n",kmers_removed_1);
  }
  

  if (cmd_line.health_check_binary){
    timestamp();
    printf("health checking the graph\n");
    db_graph_health_check(false,db_graph);
  }
  if (cmd_line.print_paths_fasta){
    timestamp();
    printf("check double Ys\n");
    fflush(stdout);
    int count_double_Ys = db_graph_find_double_Ys(ctg_length, db_graph);
    printf("%i double Ys found\n",count_double_Ys);
    timestamp();
    printf("dumping contigs %s -- coverage threshold %i\n",cmd_line.output_paths_fasta_filename,cmd_line.coverage_threshold);
    fflush(stdout);
    db_graph_print_paths(cmd_line.output_paths_fasta_filename,ctg_length,cmd_line.print_coverages,cmd_line.coverage_threshold,db_graph);    
  }

 if (cmd_line.print_supernodes_fasta){
   timestamp();
   printf("dumping supernodes %s\n",cmd_line.output_supernodes_fasta_filename);
   fflush(stdout);
   db_graph_print_supernodes(cmd_line.output_supernodes_fasta_filename,ctg_length,cmd_line.print_coverages,db_graph);
 }
 

  if (cmd_line.detect_bubbles){
    timestamp();
    printf("detect bubbles\n");
    fflush(stdout);
    db_graph_detect_vars_clean_bubbles(cmd_line.detect_vars_branch_length,db_graph);
  }
  
  if (cmd_line.detect_dirty_bubbles){
    timestamp();
    printf("detect dirty bubbles\n");
    db_graph_detect_vars_dirty_bubbles(cmd_line.detect_vars_branch_length,db_graph);
  }

 if (cmd_line.dump_binary){
   timestamp();
   printf("dumping graph %s\n",cmd_line.output_ctx_filename);
   fflush(stdout);
   db_graph_dump_binary(cmd_line.output_ctx_filename,&db_node_check_flag_not_pruned,db_graph);
 }

if (cmd_line.dump_hash_table){
   timestamp();
   printf("dumping hash table %s\n",cmd_line.output_hash_filename);
   fflush(stdout);
   db_graph_dump_hash_table(cmd_line.output_hash_filename,db_graph);
 }

 timestamp(); 
 printf("DONE\n");
 return 0;
  }
  

void timestamp(){
 time_t ltime;
 ltime = time(NULL);
 printf("\n-----\n%s",asctime(localtime(&ltime)));
 fflush(stdout);
}
