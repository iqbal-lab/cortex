#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>

int main(int argc, char **argv){

  FILE *fp_fnames;
  char filename[1000];
  int hash_key_bits;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int fastq; //if 0 entry is fasta otherwise is quality cut-off
  int bucket_size;

  long long bad_reads     = 0;

  FILE * fout;

  void print_node_binary(dBNode * node){
    db_node_print_binary(fout,node);
  }


  //command line arguments 
  fp_fnames= fopen(argv[1], "r");    //open file of file names
  kmer_size        = atoi(argv[2]); 
  hash_key_bits    = atoi(argv[3]); //number of buckets: 2^hash_key_bits
  bucket_size      = atoi(argv[4]);
  fastq            = atoi(argv[5]);
  fout             = fopen(argv[6],"w");
  DEBUG            = atoi(argv[7]);
  

  fprintf(stdout,"Input file of filenames: %s\n",argv[1]);
  fprintf(stdout,"Output bin file: %s\n",argv[6]);
  fprintf(stdout,"Kmer size: %d hash_table_size (%d bits): %d - bucket size: %d - total size: %qd\n",kmer_size,hash_key_bits,1 << hash_key_bits, bucket_size, ((long long) 1<<hash_key_bits)*bucket_size);
  if (fastq>0){
    fprintf(stdout,"quality cut-off: %i\n",fastq);
  }

  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(hash_key_bits,bucket_size, 10,kmer_size);
  fprintf(stdout,"table created: %d\n",1 << hash_key_bits);


  int count_file   = 0;
  long long total_length = 0; //total sequence length
  
  //Go through all the files, loading data into the graph

  
  while (!feof(fp_fnames)){
  
    
    fscanf(fp_fnames, "%s\n", filename);
       
    long long seq_length = 0;
    count_file++;

    if (fastq>0){
      seq_length += load_fastq_from_filename_into_graph(filename,&bad_reads, fastq, 5000, db_graph);
    }
    else{
      if (fastq==0){
	seq_length += load_fasta_from_filename_into_graph(filename,&bad_reads, 500000, db_graph);
      }
    }

    total_length += seq_length;
    
    fprintf(stdout,"\n%i kmers: %qd file name:%s bad reads: %qd seq:%qd total seq:%qd\n\n",count_file,hash_table_get_unique_kmers(db_graph),filename,bad_reads,seq_length, total_length);

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
  }  
    
  printf("print nodes binary\n");
  hash_table_traverse(&print_node_binary,db_graph);
  
  fclose(fout);
  return 0;
}
