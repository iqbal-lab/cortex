#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <string.h>
#include <internal_oxford.h>
//#include <internal_tgac.h>

int main(int argc, char **argv){

  char* filename;
  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int action;
  char* dumped_binary;
  char* list_of_fastq;
  //command line arguments 
  filename         = argv[1];        //open file that lists one file per individual in the trio (population), and each of those gives list of files.
  kmer_size        = atoi(argv[2]);  //global variable defined in element.h
  hash_key_bits    = atoi(argv[3]);  //number of buckets: 2^hash_key_bits
  bucket_size      = atoi(argv[4]);
  action           = atoi(argv[5]);
  DEBUG            = atoi(argv[6]);
  dumped_binary   = argv[7];
  list_of_fastq = argv[8];

  int max_retries=10;

  fprintf(stdout,"Kmer size: %d hash_table_size (%d bits): %d\n",kmer_size,hash_key_bits,1 << hash_key_bits);


  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits);

 
  fprintf(stderr, "Start loading population\n");
  load_population_as_binaries_from_graph(filename, db_graph);
  fprintf(stderr, "Finished loading population\n");








  switch (action)
    {
    case 0:
      {
	printf("Dump binary %s\n", dumped_binary);
	break;
      } 
    case 1:
      {
	printf("remove low coverage nodes (<=1) \n");

	void apply_reset_to_all_edges(dBNode* node, Orientation or, Nucleotide nuc)
	{
	  int j;
	  for (j=0; j<NUMBER_OF_BITFIELDS_IN_BINARY_KMER; j++)
	    {
	      reset_one_edge(node, or, nuc, individual_edge_array, j);
	    }
	}
	void apply_reset_to_all_edges_2(dBNode* node )
	{
	  int j;
	  for (j=0; j<NUMBER_OF_BITFIELDS_IN_BINARY_KMER; j++)
	    {
	      db_node_reset_edges(node, individual_edge_array, j);
	    }
	}

	db_graph_remove_low_coverage_nodes(1,db_graph, &element_get_covg_union_of_all_covgs, &element_get_colour_union_of_all_colours,
					   &apply_reset_to_all_edges, &apply_reset_to_all_edges_2);
	    
	printf("dumping graph %s\n",dumped_binary);
	db_graph_dump_binary(dumped_binary,&db_node_check_status_not_pruned,db_graph);

	break;
      }
    case 2:
      {
	//detect variants looking for bubbles in the union graph of all colours.
	//apply no condition to whether one branch or other shoud be one colour or another
	int max_allowed_branch_len=500000; //500kb
	db_graph_detect_vars(max_allowed_branch_len,db_graph, &detect_vars_condition_always_true, 
			     &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);
	break;
      }
    case 3:
      {
	printf("Print supernodes containing entirely novel(non-reference) sequence\n");
	read_all_ref_chromosomes_and_mark_graph(db_graph);
	int min_covg=2; //demand coverage of at least 2
	
	//todo - add condition to make sure ignore pruned
	printf("Start printing..\n");
	db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(db_graph, &db_node_check_status_is_not_exists_in_reference,min_covg, stdout, false, NULL, 0,
										     individual_edge_array,0);//assumes interested in person in colour 0
	
	//cleanup - the printing process has set all printed nodes to visited
	//        - or to visited_and_exists_in_reference
	hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph); 
	
      break;
      
      }      
    case 4:
      {
	printf("Print supernodes that match reference at start and end, but not the middle\n");
	read_all_ref_chromosomes_and_mark_graph(db_graph);
	int minimum_covg=5; //demand covg of at least 5
	printf("Start printing..\n");
	db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, minimum_covg,
													  2,25,40, stdout, false,NULL,0, individual_edge_array,0);//assumed interested in colour 0
	//2,25 are required numbers of overlapping nodes with ref at start and end
	//40 is min number of nodes which are NOT in reference
	//cleanup - the printing process has set all printed nodes to visited
	//           - or to visited_and_exists_in_reference
	hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph); 
	break;
	
      }
    case 5:
      {
	printf("Given this list of fastq files %s, read them each and dump a fasta file for each one containing only those reads lying inside our (cleaned) graph\n", argv[8]);
	
	
	int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
	  * full_entry = true;
	  
	  if (new_entry!= true){
	    puts("new_entry has to be true for fastq\n");
	    exit(1);
	  }
	  
	  return read_sequence_from_fastq(fp,seq,max_read_length);
	}
	
	
	FILE* list_fptr = fopen(argv[8], "r");
	if (list_fptr==NULL)
	  {
	    printf("Cannot open %s\n", argv[8]);
	    exit(1);
	  }
	
	char line[MAX_FILENAME_LENGTH+1];
	
	
	
	while(fgets(line,MAX_FILENAME_LENGTH, list_fptr) !=NULL)
	  {
	    
	    //remove newline from end of line - replace with \0
	    char* p;
	    if ((p = strchr(line, '\n')) != NULL)
	      {
		*p = '\0';
	      }
	    
	    printf("Looking at %s\n", line);
	    
	    char outputfile[200];
	    sprintf(outputfile,"%s.cleaned.fasta",line);
	    FILE* out = fopen(outputfile, "w");
	    if (out ==NULL)
	      {
		printf("Cannot open %s, exiting", outputfile);
		exit(1);
	      }
	    printf("We will output to %s\n", outputfile);
	    
	    
	    FILE* read_fptr = fopen(line, "r");
	    if (read_fptr==NULL)
	      {
		printf("Cannot open %s", line);
		exit(1);
	      }
	    
	    long long b_reads=0; //will ignore base reads, but needed by API
	    int max_read_length=10000;
	    
	    read_fastq_and_print_reads_that_lie_in_graph(read_fptr, out, &file_reader, 
							 &b_reads, max_read_length, db_graph,
							 false, NULL, NULL);//last 3 arguments show is not for testing
	    fclose(out);
	    fclose(read_fptr);
	    
	  }
	fclose(list_fptr);
	break;
	
      }      
    case 6:
      {
      printf("Given this list of fastq files %s, read them each and dump a fasta file of sub_reads lying inside our (cleaned) graph. Break when a kmer is not in the graph or an edge is not\n", argv[8]);


      int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
	* full_entry = true;
	
	if (new_entry!= true){
	  puts("new_entry has to be true for fastq\n");
	  exit(1);
	}
	
	return read_sequence_from_fastq(fp,seq,max_read_length);
      }
      

      FILE* list_fptr = fopen(argv[8], "r");
      if (list_fptr==NULL)
	{
	  printf("Cannot open %s\n", argv[8]);
	  exit(1);
	}
      
      char line[MAX_FILENAME_LENGTH+1];
      

      
      while(fgets(line,MAX_FILENAME_LENGTH, list_fptr) !=NULL)
	{
	  
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(line, '\n')) != NULL)
	    {
	      *p = '\0';
	    }

	  printf("Looking at %s\n", line);

	  char outputfile[200];
	  sprintf(outputfile,"%s.cleaned.fasta",line);
	  FILE* out = fopen(outputfile, "w");
	  if (out ==NULL)
	    {
	      printf("Cannot open %s, exiting", outputfile);
	      exit(1);
	    }
	  printf("We will output to %s\n", outputfile);


	  FILE* read_fptr = fopen(line, "r");
	  if (read_fptr==NULL)
	    {
	      printf("Cannot open %s", line);
	      exit(1);
	    }

	  long long b_reads=0; //will ignore base reads, but needed by API
	  int max_read_length=10000;
	  


	  read_fastq_and_print_subreads_that_lie_in_graph_breaking_at_edges_or_kmers_not_in_graph(read_fptr, out, &file_reader, 
												  &b_reads, max_read_length, db_graph,
												  individual_edge_array, 0, //assumes interested in colour 0
												  false, NULL, NULL);//last 3 arguments show is not for testing
	  fclose(out);
	  fclose(read_fptr);
	  
	}
      fclose(list_fptr);

      break;

      }
    case 7:
      {
	printf("Make SV calls based on the trusted-path algorithm, against the whole human genome\n");

	
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
	
	set_ref_chromosome_file_pointers(ref_chroms, 25);
	
	
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
								  0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true);
	    
	    
	    fclose(chrom_fptr);
	    fclose(out_fptr);
	  }
	//cleanup
	for(i=0; i<25; i++)
	  {
	    free(ref_chroms[i]);
	    free(output_files[i]);
	  }
	free(ref_chroms);
	free(output_files);
	break;
	
      }
    case 8:
      {
	//detect homozygous bubbles in an an individual by comparison with the reference.
	//the graph is a union of colours 0 (ref) and 1 (individual)
	//Apply condition that one branch lies in the reference, but has no coverage from the individual.


	boolean condition_is_hom_nonref


	break;
      }

      

    }



  return 0;
}
