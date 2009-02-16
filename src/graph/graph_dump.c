
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
  int clip_limit;
  int fastq; //if 0 entry is fastq otherwise is quality cut-off
  long long kmers_clipped = 0;
  long long count_kmers   = 0;
  long long bad_reads     = 0;

  FILE * fout;

  void print_node_binary(dBNode * node){
    db_node_print_binary(fout,node);
  }


  void tip_clipping(dBNode * node){

    kmers_clipped += db_graph_clip_tip(node,clip_limit,db_graph);
    
  }


  //command line arguments 
  fp_fnames= fopen(argv[1], "r");    //open file of file names
  kmer_size        = atoi(argv[2]); 
  hash_key_bits    = atoi(argv[3]); //number of buckets: 2^hash_key_bits
  clip_limit       = atoi(argv[4]);
  fastq            = atoi(argv[5]);
  fout             = fopen(argv[6],"w");
  DEBUG            = atoi(argv[7]);
  
  fprintf(stderr,"Kmer size: %d hash_table_size (%d bits): %d\n",kmer_size,hash_key_bits,1 << hash_key_bits);
  if (fastq>0){
    fprintf(stderr,"quality cut-off: %i\n",fastq);
  }

  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(1 << hash_key_bits,kmer_size);
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits);


  int count_file   = 0;
  long long total_length = 0; //total sequence length

  //Go through all the files, loading data into the graph
  while (!feof(fp_fnames)){

    //int count_bad_reads = 0;
    fscanf(fp_fnames, "%s\n", filename);
    
    int seq_length = 0;
    count_file++;

    if (fastq>0){
      seq_length += load_fastq_data_from_filename_into_graph(filename, &count_kmers, &bad_reads, fastq, 5000, db_graph);
    }
    else{
      seq_length += load_fasta_data_from_filename_into_graph(filename, &count_kmers, &bad_reads, 5000, db_graph);
    }

    total_length += seq_length;
    
    fprintf(stderr,"\n%i kmers: %qd file name:%s bad reads: %qd seq:%i total seq:%qd\n\n",count_file,count_kmers,filename,bad_reads,seq_length, total_length);

    //print mem status
    FILE* fmem=fopen("/proc/self/status", "r");
    char line[500];
    while (fgets(line,500,fmem) !=NULL){
      if (line[0] == 'V' && line[1] == 'm'){
	fprintf(stderr,"%s",line);
      }
    }
    fclose(fmem);
    fprintf(stderr,"************\n");
  }  
    
  if (clip_limit>0){
    printf("clip tips\n");
    hash_table_traverse(&tip_clipping,db_graph);
    printf("\nkmers clipped %qd\n\n",kmers_clipped);
  }
  
  printf("print nodes binary\n");
  hash_table_traverse(&print_node_binary,db_graph);
  
  fclose(fout);
  return 1;
}
