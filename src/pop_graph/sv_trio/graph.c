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
  boolean allow_for_ref=false;
  
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
  printf("table created: %d\n",1 << hash_key_bits);

 
  printf("Start loading population\n");
  load_population_as_binaries_from_graph(filename, db_graph);
  printf("Finished loading population\n");
  printf("Total kmers in table: %qd\n", hash_table_get_unique_kmers(db_graph));







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
	db_graph_detect_vars(stdout, max_allowed_branch_len,db_graph, &detect_vars_condition_flanks_at_least_3,
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

	/*
	//Assume we have loaded the following into the dbg
	// colour 0 : reference
	// colour 1 : CEPH
	// colour 2 : YRI

	int num_ceph=49; //49 ceph people have slx data
	int num_yri =57; //57 yri have slx data
	double avg_ceph_covg_per_person = 0;
	double avg_yri_covg_per_person = 0; 
	//************************ MUST fix the above - is not 0!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	  

	//Call bubbles with high differentiation between populations
	//We don't knwo which branch is the ref branch in detect vars, but this does not matter
	boolean condition_AF_of_branch1_gtr_90_in_yri_and_less_10_in_ceph(VariantBranchesAndFlanks* var)
	{
	  //return true if the allele frequency that maximises the likelihood of seeing branch 1 > 90% in YRI, and <10% in CEPH
	  // **AND** the allele frequency of branch2  < 10% in YRI, and >10% in CEPH 

	  int br1_multiplicity_in_br1[var->len_one_allele];
	  int br1_multiplicity_in_br2[var->len_one_allele];
	  int br2_multiplicity_in_br2[var->len_other_allele];
	  int br2_multiplicity_in_br1[var->len_other_allele];
	  int* br1_multiplicity_in_br1_ptr[var->len_one_allele];
	  int* br1_multiplicity_in_br2_ptr[var->len_one_allele];
	  int* br2_multiplicity_in_br2_ptr[var->len_other_allele];
	  int* br2_multiplicity_in_br1_ptr[var->len_other_allele];
	  int br1_normalised_covg_ceph[var->len_one_allele];
	  int br2_normalised_covg_ceph[var->len_other_allele];
	  int* br1_normalised_covg_ptr_ceph[var->len_one_allele];
	  int* br2_normalised_covg_ptr_ceph[var->len_other_allele];
	  int br1_normalised_covg_yri[var->len_one_allele];
	  int br2_normalised_covg_yri[var->len_other_allele];
	  int* br1_normalised_covg_ptr_yri[var->len_one_allele];
	  int* br2_normalised_covg_ptr_yri[var->len_other_allele];
	  
	  int i;
	  for (i=0; i<var->len_one_allele; i++)
	    {
	      br1_multiplicity_in_br1[i]=0;
	      br1_multiplicity_in_br2[i]=0;
	      br1_multiplicity_in_br1_ptr[i]=&br1_multiplicity_in_br1[i];
	      br1_multiplicity_in_br2_ptr[i]=&br1_multiplicity_in_br2[i];
	      br1_normalised_covg_ceph[i]=0;
	      br1_normalised_covg_ptr_ceph[i] = &br1_normalised_covg_ceph[i];
	      br1_normalised_covg_yri[i]=0;
	      br1_normalised_covg_ptr_yri[i] = &br1_normalised_covg_yri[i];
	    }
	  for (i=0; i<var->len_other_allele; i++)
	    {
	      br2_multiplicity_in_br1[i]=0;
	      br2_multiplicity_in_br2[i]=0;
	      br2_multiplicity_in_br1_ptr[i]=&br2_multiplicity_in_br1[i];
	      br2_multiplicity_in_br2_ptr[i]=&br2_multiplicity_in_br2[i];
	      br2_normalised_covg_ceph[i]=0;
	      br2_normalised_covg_ptr_ceph[i] = &br1_normalised_covg_ceph[i];
	      br2_normalised_covg_yri[i]=0;
	      br2_normalised_covg_ptr_yri[i] = &br1_normalised_covg_yri[i];

	    }
	  
	  get_node_multiplicities(var->one_allele, var->len_one_allele, var->other_allele, var->len_other_allele, 
				  br1_multiplicity_in_br1_ptr,  br2_multiplicity_in_br2_ptr,
				  br1_multiplicity_in_br2_ptr,  br2_multiplicity_in_br1_ptr,
				  false,NULL, NULL);
	  get_rough_normalised_coverage(num_ceph, avg_ceph_covg_per_person, 36,
					br1_normalised_covg_ceph, br2_normalised_covg_ceph,
					br1_multiplicity_in_br1, br1_multiplicity_in_br2,
					var, db_graph, 1, 0);
	  get_rough_normalised_coverage(num_yri, avg_yri_covg_per_person, 36,
					br1_normalised_covg_yri, br2_normalised_covg_yri,
					br1_multiplicity_in_br1, br1_multiplicity_in_br2,
					var, db_graph, 2, 0);

	  double likelihoods_yri[2*num_yri];
	  double likelihoods_ceph[2*num_ceph];
	  double max_likelihood=LARGE_NEGATIVE_NUMBER;
	  int mle_af=-1;
	  
	  for (i=0; i<2*num_yri; i++)
	    {
	      likelihoods_yri[i]=calc_log_likelihood_of_data_seen_on_one_allele_excluding_nodes_on_both_alleles(br1_normalised_covg_yri, ref_multiplicity_in_ref, ref_multiplicity_in_alt, 
														ref_normalised_covg_yri,
														var->len_one_allele, i, avg_depth_of_covg_pop1/2, 
													    avg_read_len, db_graph->kmer_size)
		+ calc_log_likelihood_of_data_seen_on_one_allele_excluding_nodes_on_both_alleles(alt_normalised_covg, alt_multiplicity_in_alt, alt_multiplicity_in_ref, alt_normalised_covg,
												 var->len_other_allele, i, avg_depth_of_covg_pop1/2, 
											   avg_read_len, db_graph->kmer_size);
	if (likelihoods[i]>mle_af)
	  {
	    max_likelihood = likelihoods[i];
	    mle_af = i;
	  }
      }
    *mle_allele_freq = mle_af;
    return max_likelihood;;
			

	  
	  
	}
	boolean condition_AF_of_branch1_gtr_90_in_ceph_and_less_10_in_yri(VariantBranchesAndFlanks* var)
	{
	  
	}
	*/
	break;
      }
    case 9:
      {
	// Traverse graph marking any node that is not highly differentiated between pops
	// as visited. Then print supernodes for what remains. Then re-mark those nodes, and find bubbles.

	//Assume we have loaded the following into the dbg
	// colour 0 : reference
	// colour 1 : CEPH
	// colour 2 : YRI

	double avg_ceph_haploid_covg_per_person =1.075;
	double avg_yri_haploid_covg_per_person = 1.4; 
	int num_yri_haploids=114;// 2*pop size
	int num_ceph_haploids=98;//2* pop size

	void mark_nodes_gtr90_yri_less_10_ceph(dBNode* e)
	{

	  int multiplier=1; //to allow for multiples in ref

	  if (allow_for_ref==true)
	    {
	      if (e->coverage[0] > 0) //exists in reference
		{
		  multiplier=e->coverage[0];
		}
	    }
	  //so now, if the kmer is novel, multiplier is 1, otherwise multiplier is the number of times we expect to see it in one haploid.
	  printf("Multiplier is %d\n", multiplier);

	  int occurrences_of_kmer_in_yri = (int) (e->coverage[2])/avg_yri_haploid_covg_per_person;
	  int occurrences_of_kmer_in_ceph= (int) (e->coverage[1])/avg_ceph_haploid_covg_per_person;
	  int num_haploids_with_kmer_in_yri = occurrences_of_kmer_in_yri/multiplier;
	  int num_haploids_with_kmer_in_ceph = occurrences_of_kmer_in_ceph/multiplier;

	  if ( (num_haploids_with_kmer_in_yri> 0.9*num_yri_haploids)
	       && (num_haploids_with_kmer_in_ceph < 0.1*num_ceph_haploids) )
	    {
	      //do nothing
	    }
	  else
	    {
	      db_node_set_status(e, visited);
	    }
	}
	void mark_nodes_gtr90_ceph_less_10_yri(dBNode* e)
	{
	  int multiplier=1; //to allow for multiples in ref
	  if (allow_for_ref==true)
	    {
	      if (e->coverage[0] > 0) //exists in reference
		{
		  multiplier=e->coverage[0];
		}
	    }
	  //so now, if the kmer is novel, multiplier is 1, otherwise multiplier is the number of times we expect to see it in one haploid.

	  int occurrences_of_kmer_in_yri =  ( (e->coverage[2])/avg_yri_haploid_covg_per_person);
	  int occurrences_of_kmer_in_ceph= ( (e->coverage[1])/avg_ceph_haploid_covg_per_person);
	  int num_haploids_with_kmer_in_yri = occurrences_of_kmer_in_yri/multiplier;
	  int num_haploids_with_kmer_in_ceph = occurrences_of_kmer_in_ceph/multiplier;

	  if ( (num_haploids_with_kmer_in_ceph> 0.9*num_ceph_haploids)
	       && (num_haploids_with_kmer_in_yri < 0.1*num_yri_haploids) )
	    {
	      //do nothing
	    }
	  else
	    {
	      db_node_set_status(e, visited);
	    }
	}

	Edges get_colour_ceph(const dBNode* node)
	{
	  return node->individual_edges[1];
	}
	Edges get_colour_yri(const dBNode* node)
	{
	  return node->individual_edges[2];
	}
	int get_covg_ceph(const dBNode* node)
	{
	  return node->coverage[1];
	}
	int get_covg_yri(const dBNode* node)
	{
	  return node->coverage[2];
	}
	

	int max_length = 20000; //longest allowed branch in a bubble


	void traverse_graph_and_get_results(boolean do_we_allow_for_multiplicity_of_node_in_reference_genome)
	{

	  allow_for_ref = do_we_allow_for_multiplicity_of_node_in_reference_genome;

	  char allow_for_ref_char[25];
	  if (do_we_allow_for_multiplicity_of_node_in_reference_genome==true)
	    {
	      sprintf(allow_for_ref_char,"allow_for_ref");
	    }
	  else
	    {
	      sprintf(allow_for_ref_char,"not_allow_for_ref");
	    }

	  char filename_high_yri_low_ceph_sups_by_yri_edges[100];
	  char filename_high_yri_low_ceph_sings_by_yri_edges[100];
	  char filename_high_yri_low_ceph_bubbles[100];
	  sprintf(filename_high_yri_low_ceph_sups_by_yri_edges, "%s_prefilter_nodes_supernodes_yri_gtr_90_and_ceph_less_10",allow_for_ref_char);
	  sprintf(filename_high_yri_low_ceph_sings_by_yri_edges, "%s_prefilter_nodes_hub_kmers_yri_gtr_90_and_ceph_less_10",allow_for_ref_char);
	  sprintf (filename_high_yri_low_ceph_bubbles, "%s_prefilter_nodes_bubbles_yri_gtr_90_and_ceph_less_10",allow_for_ref_char);
	  char filename_high_ceph_low_yri_sups_by_ceph_edges[100];
	  char filename_high_ceph_low_yri_sings_by_ceph_edges[100];
	  char filename_high_ceph_low_yri_bubbles[100];
	  sprintf(filename_high_ceph_low_yri_sups_by_ceph_edges, "%s_prefilter_nodes_supernodes_ceph_gtr_90_and_yri_less_10",allow_for_ref_char);
	  sprintf(filename_high_ceph_low_yri_sings_by_ceph_edges, "%s_prefilter_nodes_hub_kmers_ceph_gtr_90_and_yri_less_10",allow_for_ref_char);
	  sprintf(filename_high_ceph_low_yri_bubbles, "%s_prefilter_nodes_bubbles_ceph_gtr_90_and_yri_less_10",allow_for_ref_char);
	  


	  //traverse graph, and mark everything as visited, EXCEPT stuff that is >90% freq in YRI and <10% in CEPH
	  hash_table_traverse(&mark_nodes_gtr90_yri_less_10_ceph, db_graph);
	  db_graph_print_supernodes_for_specific_person_or_pop(filename_high_yri_low_ceph_sups_by_yri_edges, filename_high_yri_low_ceph_sings_by_yri_edges,
							       max_length, db_graph, individual_edge_array,2);

	  //clean up node status marking
	  hash_table_traverse(&db_node_set_status_to_none, db_graph);
	  
	  
	  FILE* fptr_yri_high = fopen(filename_high_yri_low_ceph_bubbles, "w");
	  if (fptr_yri_high==NULL)
	    {
	      printf("Cannot open %s", filename_high_yri_low_ceph_bubbles);
	      exit(1);
	    }
	  //AGAIN - traverse graph, and mark everything as visited, EXCEPT stuff that is >90% freq in YRI and <10% in CEPH
	  hash_table_traverse(&mark_nodes_gtr90_yri_less_10_ceph, db_graph);
	  db_graph_detect_vars(fptr_yri_high, max_length, db_graph,
			       &detect_vars_condition_always_true,
			       get_colour_yri, get_covg_yri);
	  fclose(fptr_yri_high);

	  //clean up node status markings
	  hash_table_traverse(&db_node_set_status_to_none, db_graph);
	  

	  //traverse graph, and mark everything as visited, EXCEPT stuff that is >90% freq in CEPH and <10% in YRI
	  hash_table_traverse(&mark_nodes_gtr90_ceph_less_10_yri, db_graph);
	  db_graph_print_supernodes_for_specific_person_or_pop(filename_high_ceph_low_yri_sups_by_ceph_edges, filename_high_ceph_low_yri_sings_by_ceph_edges,
							       max_length, db_graph, individual_edge_array,1);

	//clean up node status markings
	  hash_table_traverse(&db_node_set_status_to_none, db_graph);


	  //AGAIN - traverse graph, and mark everything as visited, EXCEPT stuff that is >90% freq in CEPH and <10% in YRI
	  hash_table_traverse(&mark_nodes_gtr90_ceph_less_10_yri, db_graph);
	  
	  FILE* fptr_ceph_high = fopen(filename_high_ceph_low_yri_bubbles, "w");
	  if (fptr_ceph_high==NULL)
	    {
	      printf("Cannot open %s", filename_high_ceph_low_yri_bubbles);
	      exit(1);
	    }
	  
	  db_graph_detect_vars(fptr_ceph_high, max_length, db_graph,
			       &detect_vars_condition_always_true,
			       get_colour_ceph, get_covg_ceph);
	  fclose(fptr_ceph_high);
	  //clean up node status markings
	  hash_table_traverse(&db_node_set_status_to_none, db_graph);



	}


	printf("Get supernodes and bubbles that are highly differentiated between CEPH and YRI - with normalisation of reference:\n");
	traverse_graph_and_get_results(true);
	printf("Get supernodes and bubbles that are highly differentiated between CEPH and YRI - without making allowance for whether nodes are known in reference:\n");
	traverse_graph_and_get_results(false);

	break;
      }
      
      

    }

  hash_table_free(&db_graph);

  return 0;
}
