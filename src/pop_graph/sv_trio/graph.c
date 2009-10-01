#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <string.h>

int main(int argc, char **argv){

  char* filename;
  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int action;
  char* list_of_fasta;

  //command line arguments 
  filename         = argv[1];        //open file that lists one file per individual in the trio (population), and each of those gives list of files.
  kmer_size        = atoi(argv[2]);  //global variable defined in element.h
  hash_key_bits    = atoi(argv[3]);  //number of buckets: 2^hash_key_bits
  bucket_size      = atoi(argv[4]);
  action           = atoi(argv[5]);
  DEBUG            = atoi(argv[6]);
  list_of_fasta   = argv[7];


  int max_retries=10;

  fprintf(stdout,"Kmer size: %d hash_table_size (%d bits): %d\n",kmer_size,hash_key_bits,1 << hash_key_bits);


  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits);

 
  fprintf(stderr, "Start loading population\n");
  load_population_as_binaries_from_graph(filename, db_graph);
  fprintf(stderr, "Finished loading population\n");

  // ***************************************
  //locations of chromosome reference files:
  // ***************************************

  char** ref_chroms = malloc( sizeof(char*) * 25); //one for each chromosome, ignoring haplotypes like MHC
  if (ref_chroms==NULL)
    {
      printf("OOM. Give up can't even allocate space for the names of the ref chromosome files\n");
      exit(1);
    }
  int i;
  for (i=0; i< 25; i++)
    {
      ref_chroms[i] = malloc(sizeof(char)*150); //filenames including path are about 100 characters. 50 characters of leeway
      if (ref_chroms[i]==NULL)
	{
	  printf("OOM. Giveup can't even allocate space for the names of the ref chromosome file i = %d\n",i);
	  exit(1);
	}
    }

  ref_chroms[0] = "/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/Homo_sapiens.NCBI36.52.dna.chromosome.MT.fa";
  
  for (i=1; i<23; i++)
    {
      sprintf(ref_chroms[i], "/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/Homo_sapiens.NCBI36.52.dna.chromosome.%i.fa", i);
    }
  ref_chroms[23]="/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/Homo_sapiens.NCBI36.52.dna.chromosome.X.fa";
  ref_chroms[24]="/nfs/1000g-work/G1K/work/zi/projects/marzam/humref/Homo_sapiens.NCBI36.52.dna.chromosome.Y.fa";



  // **** set up one output file per chromosome ***** //
  
  char** output_files = malloc( sizeof(char*) * 25); //one for each chromosome, ignoring haplotypes like MHC
  if (output_files==NULL)
    {
      printf("OOM. Give up can't even allocate space for the names of the output  files \n");
      exit(1);
    }

  for (i=0; i< 25; i++)
    {
      output_files[i] = malloc(sizeof(char)*100); 
      if (output_files[i]==NULL)
	{
	  printf("OOM. Giveup can't even allocate space for the names of the ref chromosome file i = %d\n",i);
	  exit(1);
	}
    }

  output_files[0] = "sv_called_in_MT";
  
  for (i=1; i<23; i++)
    {
      sprintf(output_files[i], "sv_called_in_chrom_%i", i);
    }
  output_files[23]="sv_called_in_chrom_X";
  output_files[24]="sv_called_in_chrom_Y";





  switch (action)
    {
    case 0:
      {
	printf("Make SV calls based on the trusted-path/supernode algorithm, against the whole genome\n");
 
	int min_fiveprime_flank_anchor = 2;
	int min_threeprime_flank_anchor= 21;
	int max_anchor_span =  20000;
	int length_of_arrays = 40000;
	int min_covg =1;
	int max_covg = 10000000;
	int max_expected_size_of_supernode=20000;
      
      
      //ignore mitochondrion for now, so start with i=1
      for (i=1; i<25; i++) 
	{
	  printf("Call SV comparing individual with chromosome %s\n", ref_chroms[i]);
	  
	  FILE* chrom_fptr = fopen(ref_chroms[i], "r");
	  if (chrom_fptr==NULL)
	    {
	      printf("Cannot open %s \n", ref_chroms[i]);
	      exit(1);
	    }
      
	  FILE* out_fptr = fopen(output_files[i], "w");
	  if (out_fptr==NULL)
	    {
	      printf("Cannot open %s for output\n", output_files[i]);
	      exit(1);
	    }
	  
	  //Note we assume person 0 is the reference, and person 1 is the person we are interested in
	  int ret = db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
								individual_edge_array, 0,
								min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
								max_expected_size_of_supernode, length_of_arrays, db_graph, out_fptr,
								0, NULL, NULL, NULL, NULL, NULL);

	  
	  fclose(chrom_fptr);
	  fclose(out_fptr);
	}
      break;

      }
    case 1:
      {
	printf("For each reference chromosome, print a fasta file which has a 1 if the 31-mer starting at that base exists precisely once in the reference, and 0 otherwise\n");

	
	char** uniq_files = malloc( sizeof(char*) * 25); //one for each chromosome, ignoring haplotypes like MHC
	if (uniq_files==NULL)
	      {
		printf("OOM. Give up can't even allocate space for the names of the uniq  files \n");
		exit(1);
	      }
	for (i=0; i< 25; i++)
	      {
		uniq_files[i] = malloc(sizeof(char)*100); 
		if (uniq_files[i]==NULL)
		  {
		    printf("OOM. Giveup can't even allocate space for the names of the uniqueness files for  i = %d\n",i);
		    exit(1);
		  }
	      }
	    
	uniq_files[0] = "kmer_uniqueness_in_chrom_MT";
	for (i=1; i<23; i++)
	  {
	    sprintf(uniq_files[i], "kmer_uniqueness_in_chrom_%i", i);
	  }
	uniq_files[23]="kmer_uniqueness_in_chrom_X";
	uniq_files[24]="kmer_uniqueness_in_chrom_Y";
	    
	    
	for (i=1; i<25; i++) 
	  {
	    printf("Print kmer uniqueness-in-entire-ref file for  %s\n", ref_chroms[i]);
	    
	    FILE* chrom_fptr = fopen(ref_chroms[i], "r");
	    if (chrom_fptr==NULL)
	      {
		printf("Cannot open %s \n", ref_chroms[i]);
		exit(1);
	      }
	    
	    FILE* out_fptr = fopen(uniq_files[i], "w");
	    if (out_fptr==NULL)
	      {
		printf("Cannot open %s for output\n", uniq_files[i]);
		exit(1);
	      }
	    
	    void print_uniqueness(dBNode* node)
	      {
		if ((node==NULL) || (node->kmer == ~0) )
		  {
		    fprintf(out_fptr, "0");
		  }
		else
		  {
		    // we assume that the reference is individual at index 0.
		    if (db_node_get_coverage(node, individual_edge_array,0)==1)
		      {
			//this kmer exists precisely once in the ref
			fprintf(out_fptr, "1");
		      }
		    else
		      {
			fprintf(out_fptr,"0");
		      }
		  }
	      }

	    apply_to_all_nodes_in_path_defined_by_fasta(&print_uniqueness, chrom_fptr, 10000, db_graph);//work through fasta in chunks of size 10kb
	    fclose(out_fptr);
	    fclose(chrom_fptr);


	  }
	
	break;
      }
    case 2:
      {
	printf("Take the file %s, which was input as argument 7, and for each fasta file listed within, print out its own output file of coverage/number_of_occurrences_in_reference for each edge, tab separated. We expect no N's in the fasta", list_of_fasta);

	FILE* fptr = fopen(list_of_fasta, "r");
	if (fptr==NULL)
	  {
	    printf("Cannot open %s\n", list_of_fasta);
	    exit(1);
	  }

	//file contains a list of fasta file names
	char line[MAX_FILENAME_LENGTH+1];
	
	while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
	  {
	    
	    //remove newline from endof line- replace with \0
	    char* p;
	    if ((p = strchr(line, '\n')) != NULL)
	      *p = '\0';
	    
	    //open fasta file
	    FILE* path_fptr = fopen(line, "r");
	    if (path_fptr==NULL)
	      {
		printf("Cannot open %s\n", line);
		exit(1);
	      }

	    //create output file 
	    char name[300];
	    sprintf(name,"%s.covg",line);

	    FILE* out_fptr = fopen(name, "w");
	    if (out_fptr==NULL)
	      {
		printf("Cannot open %s\n", name);
		exit(1);
	      }
	    
	    void print_covg_and_ref_covg(dBNode* node)
	      {
		if ((node==NULL) || (node->kmer == ~0) )
		  {
		    fprintf(out_fptr, "N in fasta - unexpected\t");
		  }
		else
		  {
		      // we assume that the reference is individual at index 0.
		    int ref_covg =  db_node_get_coverage(node, individual_edge_array,0);
		    int covg     =  db_node_get_coverage(node, individual_edge_array,1);
		    fprintf(out_fptr, "%d/%d\t", covg, ref_covg);
		  }
	      }
	    
	    apply_to_all_nodes_in_path_defined_by_fasta(&print_covg_and_ref_covg, path_fptr, 50, db_graph);//work through fasta in chunks of size 50 bases - fasta not too big
	    fclose(out_fptr);
	    fclose(path_fptr);
	    
	    
	  }
	
	fclose(fptr);
	break;
      }

      

    }


  //cleanup

  /*
  for(i=0; i<25; i++)
    {
      free(ref_chroms[i]);
      free(output_files[i]);
    }

  free(ref_chroms);
  free(output_files);
  */

  printf("Finished getting SV from all chromosomes\n");
  return 0;
}
