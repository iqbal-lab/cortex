#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <dB_graph_population.h>

int main(int argc, char **argv){

  char* filename;
  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int action;


  //command line arguments 
  filename         = argv[1];        //open file that lists one file per individual in the trio (population), and each of those gives list of files.
  kmer_size        = atoi(argv[2]);  //global variable defined in element.h
  hash_key_bits    = atoi(argv[3]);  //number of buckets: 2^hash_key_bits
  bucket_size      = atoi(argv[4]);
  action           = atoi(argv[5]);
  DEBUG            = atoi(argv[6]);
  
  int max_retries=10;

  fprintf(stdout,"Kmer size: %d hash_table_size (%d bits): %d\n",kmer_size,hash_key_bits,1 << hash_key_bits);

  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits);

  load_population_as_binaries_from_graph(filename, db_graph);


  switch (action)
    {
    case 0:
      printf("Compare with chromosome 1\n");

      FILE* chrom_fptr = fopen("/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/Homo_sapiens.NCBI36.52.dna.chromosome.1.fa", "r");
      if (chrom_fptr==NULL)
	{
	  printf("Cannot open /nfs/1000g-work/G1K/work/zi/projects/marzam/humref/Homo_sapiens.NCBI36.52.dna.chromosome.1.fa \n");
	  exit(1);
	}
      
      int min_fiveprime_flank_anchor = 2;
      int min_threeprime_flank_anchor= 3;
      int max_anchor_span = 20000;
      int length_of_arrays=40000;
      int min_covg =1;
      int max_covg = 10000000;
      int max_expected_size_of_supernode=20000;

      int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 0, 
							    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							    max_expected_size_of_supernode, length_of_arrays, db_graph, stdout,
							    0, NULL, NULL, NULL, NULL, NULL);

      

    }


  printf("Finished\n");
  return 1;
}
