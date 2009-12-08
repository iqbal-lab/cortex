#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>

int main(int argc, char **argv){

  FILE *fp_fnames;
  
  char filename[1000];
  char* pair_fnames1;
  char* pair_fnames2;

  int hash_key_bits;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int quality_cut_off; //if 0 entry is fasta otherwise is quality cut-off
  int bucket_size;
  boolean rmdups_single;
  boolean rmdups_paired;
  boolean we_are_loading_some_single_ended_data;
  boolean we_are_loading_some_paired_end_data;
  int max_read_length; 

  long long bad_reads     = 0;
  long long dup_reads_in_single_ended     = 0;
  long long dup_reads_in_paired     = 0;
  long long single_seq_length = 0;
  long long paired_seq_length = 0;

  FILE * fout;

  void print_node_binary(dBNode * node){
    db_node_print_binary(fout,node);
  }


  //command line arguments 
  we_are_loading_some_single_ended_data = (boolean) atoi(argv[1]);
  we_are_loading_some_paired_end_data = (boolean) atoi(argv[2]);
  fp_fnames        = fopen(argv[3], "r");    //open file of single_ended fastq file names
  kmer_size        = atoi(argv[4]); 
  hash_key_bits    = atoi(argv[5]); //number of buckets: 2^hash_key_bits
  bucket_size      = atoi(argv[6]);
  quality_cut_off  = atoi(argv[7]);
  fout             = fopen(argv[8],"w");
  DEBUG            = atoi(argv[9]);
  pair_fnames1  = argv[10];//list of paired end fastq, left-hand mates
  pair_fnames2  = argv[11];//list of right hand mates in same order
  rmdups_single    = (boolean) atoi(argv[12]);// should we remove PCR duplicates from single ended fastq by discarding any that start at a previous read-start position?
  rmdups_paired    = (boolean) atoi(argv[13]);// should we remove PCR duplicates from PAIRED end fastq by discarding any pairs that BOTH start at a previous read-start position?
  max_read_length  = atoi(argv[14]);

  fprintf(stdout,"Input file of single ended fastq filenames: %s, and of paired end files: %s, and %s, \n",argv[3], pair_fnames1, pair_fnames2); //argv[10], argv[11]);
  fprintf(stdout,"Output bin file: %s\n",argv[8]);
  fprintf(stdout,"Kmer size: %d hash_table_size (%d bits): %d - bucket size: %d - total size: %qd\n",kmer_size,hash_key_bits,1 << hash_key_bits, bucket_size, ((long long) 1<<hash_key_bits)*bucket_size);
  if (quality_cut_off>0){
    fprintf(stdout,"quality cut-off: %i\n",quality_cut_off);
  }

  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(hash_key_bits,bucket_size, 10,kmer_size);
  fprintf(stdout,"table created: %d\n",1 << hash_key_bits);


  int num_single_ended_files_loaded   = 0;
  int num_file_pairs_loaded   = 0;
  long long total_length = 0; //total sequence length
  
  //Go through all the files, loading data into the graph


  if (we_are_loading_some_single_ended_data==true)
    {
      printf("Load files being treated as single-ended fastq\n");
      while (!feof(fp_fnames)){
	fscanf(fp_fnames, "%s\n", filename);
	num_single_ended_files_loaded++;
	
	if (quality_cut_off>0){
	  
	  if (rmdups_single==false)
	    {
	      single_seq_length += load_fastq_from_filename_into_graph(filename,&bad_reads, quality_cut_off, &dup_reads_in_single_ended, max_read_length, false,db_graph);
	    }
	  else
	    {
	      single_seq_length += load_fastq_from_filename_into_graph(filename,&bad_reads, quality_cut_off, &dup_reads_in_single_ended, max_read_length, true,db_graph);
	    }
	}
	else{
	  if (quality_cut_off==0){
	    if (rmdups_single==false)
	      {
		single_seq_length += load_fasta_from_filename_into_graph(filename,&bad_reads, &dup_reads_in_single_ended, max_read_length, false, db_graph);
	      }
	    else
	      {
		single_seq_length += load_fasta_from_filename_into_graph(filename,&bad_reads, &dup_reads_in_single_ended, max_read_length, true, db_graph);
	      }
	  }
	}
	
	total_length += single_seq_length;
	
	printf("\nNum SE files loaded:%i kmers: %qd cumulative bad reads: %qd total SE seq: %qd duplicates removed: %qd\n\n",
	       num_single_ended_files_loaded,hash_table_get_unique_kmers(db_graph),bad_reads,total_length, dup_reads_in_single_ended);
      }

    }

  if (we_are_loading_some_paired_end_data==true)
    {
      printf("Load the paired-end files\n");
      paired_seq_length = load_list_of_paired_end_fastq_into_graph(pair_fnames1, pair_fnames2, quality_cut_off, max_read_length, 
								   &bad_reads, &dup_reads_in_paired, &num_file_pairs_loaded, rmdups_paired, db_graph); 
      
      printf("\nNum PE files loaded:%i kmers: %qd cumulative bad reads: %qd total PE seq: %qd duplicates removed:%qd\n\n",
	     num_file_pairs_loaded,hash_table_get_unique_kmers(db_graph),bad_reads,paired_seq_length, dup_reads_in_paired);
      
      total_length += paired_seq_length;
    }



    hash_table_print_stats(db_graph);

    //print mem status
    //FILE* fmem=fopen("/proc/self/status", "r");
    //char line[500];
    //while (fgets(line,500,fmem) !=NULL){
    // if (line[0] == 'V' && line[1] == 'm'){
    //	fprintf(stderr,"%s",line);
    // }
    //}
    //fclose(fmem);
    //fprintf(stderr,"************\n");
    
  printf("print nodes binary\n");
  hash_table_traverse(&print_node_binary,db_graph);
  
  fclose(fout);
  return 0;
}
