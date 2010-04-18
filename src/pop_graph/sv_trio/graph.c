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

  boolean loading_colours_separately=true;//either we load a population, with a file listing individuals, and for each individual you have a list of binaries (enter 1)
                                          // or we just load one single multicolour binary (in this case enter 0)
  int number_colours_in_multicolour_binary=0;
  char* filename;
  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int action;
  char* dumped_binary;
  char* list_of_fastq;
  int remove_low_covg_nodes;
  char* detectvars_filename;
  char* detectvars_after_remv_ref_bubble_filename;
  char* detectvars_hom_nonref_filename;
  char* ref_assisted_filename;
  char* ref_fasta;
  int low_cov_thresh;
  int D_slx_na12878;//depth of covg
  int D_slx_na19240;
  int D_454_na12878;//depth of covg
  int D_454_na19240;
  int R_slx;//mean read length
  int R_454;//mean read length

  //command line arguments 
  loading_colours_separately = (boolean) atoi(argv[1]);
  number_colours_in_multicolour_binary = atoi(argv[2]);
  filename         = argv[3];        //open file that lists one file per individual in the trio (population), and each of those gives list of files.
  kmer_size        = atoi(argv[4]);  //global variable defined in element.h
  hash_key_bits    = atoi(argv[5]);  //number of buckets: 2^hash_key_bits
  bucket_size      = atoi(argv[6]);
  action           = atoi(argv[7]);
  DEBUG            = atoi(argv[8]);
  dumped_binary   = argv[9];
  list_of_fastq = argv[10];
  remove_low_covg_nodes = atoi(argv[11]);
  detectvars_filename = argv[12];
  detectvars_after_remv_ref_bubble_filename = argv[13];
  detectvars_hom_nonref_filename=argv[14];
  ref_assisted_filename=argv[15];
  ref_fasta = argv[16];
  low_cov_thresh = atoi(argv[17]);
  D_slx_na12878 = atoi(argv[18]);
  D_slx_na19240 = atoi(argv[19]);
  D_454_na12878 = atoi(argv[20]);
  D_454_na19240 = atoi(argv[21]);
  R_slx = atoi(argv[22]);
  R_454 = atoi(argv[23]);

  int max_retries=10;

  fprintf(stdout,"Kmer size: %d hash_table_size (%d bits): %d\n",kmer_size,hash_key_bits,1 << hash_key_bits);


  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  printf("table created: %d\n",1 << hash_key_bits);

  if (loading_colours_separately==true)
    {
      printf("Start loading population\n");
      load_population_as_binaries_from_graph(filename, db_graph);
      printf("Finished loading population\n");
      printf("Total kmers in table: %qd\n", hash_table_get_unique_kmers(db_graph));
    }
  else
    {
      if (number_colours_in_multicolour_binary< NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  printf("Need to port this code from the _wth_genotyping repo\n");
	  exit(1);
	  // int kmers_loaded = load_multicolour_binary_with_strictly_less_colours_from_filename_into_graph(filename, db_graph, number_colours_in_multicolour_binary);
	  //printf("Loaded the single multicolour binary %s, and got %d kmers\n", filename, kmers_loaded);
	}
      else if (number_colours_in_multicolour_binary==NUMBER_OF_INDIVIDUALS_PER_POPULATION)
	{
	  int kmers_loaded = load_multicolour_binary_data_from_filename_into_graph(filename, db_graph);
	  printf("Loaded the single multicolour binary %s, and got %d kmers\n", filename, kmers_loaded);
	}
      else
	{
	  printf("Trying to load a binary with too many colours");
	  exit(1);
	}
    }





  switch (action)
    {
    case 0:
      {
	printf("Dump binary %s\n", dumped_binary);
	db_graph_dump_binary(dumped_binary,&db_node_check_status_not_pruned,db_graph);

	break;
      } 
    case 1:
      {
	printf("remove low coverage nodes (<= %d ) \n", low_cov_thresh);

	void apply_reset_to_all_edges(dBNode* node, Orientation or, Nucleotide nuc)
	{
	  int j;
	  for (j=0; j<NUMBER_OF_INDIVIDUALS_PER_POPULATION; j++)
	    {
	      reset_one_edge(node, or, nuc, individual_edge_array, j);
	    }
	}
	void apply_reset_to_all_edges_2(dBNode* node )
	{
	  int j;
	  for (j=0; j<NUMBER_OF_INDIVIDUALS_PER_POPULATION; j++)
	    {
	      db_node_reset_edges(node, individual_edge_array, j);
	    }
	}

	db_graph_remove_low_coverage_nodes(low_cov_thresh,db_graph, &element_get_covg_union_of_all_covgs, &element_get_colour_union_of_all_colours,
					   &apply_reset_to_all_edges, &apply_reset_to_all_edges_2);
	    
	printf("dumping graph %s\n",dumped_binary);
	db_graph_dump_binary(dumped_binary,&db_node_check_status_not_pruned,db_graph);

	break;
      }
    case 2:
      {

	//detect variants looking for bubbles in the union graph of all colours.
	//apply no condition to whether one branch or other shoud be one colour or another
	int max_allowed_branch_len=50000; //5000
	db_graph_detect_vars(stdout, max_allowed_branch_len,db_graph, 
			     &detect_vars_condition_flanks_at_least_3,
			     &db_node_action_set_status_visited,
			     &db_node_action_set_status_visited,
			     &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs,
			     &print_no_extra_info);
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
	//Make bubble calls for homozygous non-ref variants, using two colours (ref=0, individual=1)
	
	boolean condition_is_hom_nonref(VariantBranchesAndFlanks* var)
	{
	  //Assumes the reference is colour 0 and the individual is colour 1
	  int covg_threshold = 1;
	  int i;
	  int count_how_many_nodes_in_one_allele_have_covg_by_indiv=0;
	  int count_how_many_nodes_in_other_allele_have_covg_by_indiv=0;
	  int count_how_many_nodes_in_one_allele_have_covg_by_ref=0;
	  int count_how_many_nodes_in_other_allele_have_covg_by_ref=0;

	  for (i=0; i< var->len_one_allele; i++)
	    {
	      if ( (var->one_allele)[i]->coverage[1] >= covg_threshold )
		{
		  count_how_many_nodes_in_one_allele_have_covg_by_indiv++;
		}
	      if ( (var->one_allele)[i]->coverage[0] > 0 )
		{
		  count_how_many_nodes_in_one_allele_have_covg_by_ref++;
		}
	    }
	  for (i=0; i< var->len_other_allele; i++)
	    {
	      if ( (var->other_allele)[i]->coverage[1] >= covg_threshold )
		{
		  count_how_many_nodes_in_other_allele_have_covg_by_indiv++;
		}
	      if ( (var->other_allele)[i]->coverage[0] > 0 )
		{
		  count_how_many_nodes_in_other_allele_have_covg_by_ref++;
		}
	    }
	  if (//individual has branch1 but not branch2 
	      (count_how_many_nodes_in_one_allele_have_covg_by_indiv==var->len_one_allele)
	       &&
	      (count_how_many_nodes_in_other_allele_have_covg_by_indiv<=1)//last node of the two branches is same
	       &&
	      //reference has branch2 only
	      (count_how_many_nodes_in_one_allele_have_covg_by_ref<=1)//Mario - do you agree?
	      &&
	      (count_how_many_nodes_in_other_allele_have_covg_by_ref==var->len_other_allele)
	      )
	    {
	      return true;
	    }
	  else if (//individual has branch2 but not branch1
		   (count_how_many_nodes_in_one_allele_have_covg_by_indiv<=1)
		   &&
		   (count_how_many_nodes_in_other_allele_have_covg_by_indiv==var->len_one_allele)
		   &&
		   //reference has branch1 only
		   (count_how_many_nodes_in_one_allele_have_covg_by_ref==var->len_other_allele)
		   &&
		   (count_how_many_nodes_in_other_allele_have_covg_by_ref<=1)
		   )
	    {
	      return true;
	    }
	  else
	    {
	      return false;
	    }
		   
	}

	int max_allowed_branch_len=500; 
	db_graph_detect_vars(stdout, max_allowed_branch_len,db_graph, &condition_is_hom_nonref,
			     &db_node_action_set_status_visited,
			     &db_node_action_set_status_visited,
			     &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, &print_no_extra_info);

	break;
      }
    case 8:
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
	    db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, 1, 
							individual_edge_array, 0,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, db_graph, out_fptr,
							0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &db_variant_action_do_nothing,
							&print_no_extra_info);
	    
	    
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
    case 9:
      {


	int get_covg_ref(const dBNode* e)
	{
	  return e->coverage[0];
	}
	int get_covg_indiv(const dBNode* e)
	{
	  return e->coverage[1];
	}


	boolean detect_vars_condition_covg_gtr_1(VariantBranchesAndFlanks* var)
	{
	  int i;
	  for (i=0; i<var->len_one_allele; i++)
	    {
	      if (db_node_get_coverage((var->one_allele)[i], individual_edge_array, 1) <=1 )
		{
		  return false;
		}
	    }
	  
	  for (i=0; i<var->len_other_allele; i++)
	    {
	      if (db_node_get_coverage((var->other_allele)[i], individual_edge_array, 1) <=1)
		{
		  return false;
		}
	    }
	  
	  return true;
	}
	

	boolean detect_vars_condition_is_hom_nonref(VariantBranchesAndFlanks* var)
	{
	  return detect_vars_condition_is_hom_nonref_given_colour_funcs_for_ref_and_indiv(var, &get_covg_ref, &get_covg_indiv);
	}

	boolean detect_vars_condition_is_hom_nonref_with_covg_gtr_1(VariantBranchesAndFlanks* var)
	{
	  return detect_vars_condition_is_hom_nonref_with_covg_gtr_1_given_colour_funcs_for_ref_and_indiv(var, &get_covg_ref, &get_covg_indiv);
	}


	//I don't want to affect the nodes that are in the reference, but not the individual
	//Also, if it has covg 1 in the individual, but also covg 1 in the reference, I dont want to touch the reference edge
	void apply_reset_to_all_edges(dBNode* node, Orientation or, Nucleotide nuc)
	{
	  // int j;
	  //for (j=0; j<NUMBER_OF_INDIVIDUALS_PER_POPULATION; j++)
	  //  {
	      reset_one_edge(node, or, nuc, individual_edge_array, 1);
	      //  }
	}
	void apply_reset_to_all_edges_2(dBNode* node )
	{
	  int j;
	  //for (j=0; j<NUMBER_OF_INDIVIDUALS_PER_POPULATION; j++)
	  //  {
	      db_node_reset_edges(node, individual_edge_array, 1);
	      //  }
	}

	if (remove_low_covg_nodes>0)
	  {
	    printf("remove low coverage nodes (<= %d ) \n", low_cov_thresh);
	    db_graph_remove_low_coverage_nodes(low_cov_thresh,db_graph, &element_get_covg_colour1, &element_get_colour1,
					       &apply_reset_to_all_edges, &apply_reset_to_all_edges_2);
	    
	    printf("dumping graph %s\n",dumped_binary);
	    db_graph_dump_binary(dumped_binary,&db_node_check_status_not_pruned,db_graph);
	  }


	// STEP1 - detect vars, calls hets, with no conditions

	int max_allowed_branch_len=50000; 
	FILE* detect_vars_fptr = fopen(detectvars_filename, "w");
	if (detect_vars_fptr==NULL)
	  {
	    printf("Cannot open %s, so exit\n", detectvars_filename);
	    exit(1);
	  }
	printf("Going to output ref free hets to %s\n", detectvars_filename);
	db_graph_detect_vars(detect_vars_fptr, max_allowed_branch_len,db_graph, &detect_vars_condition_covg_gtr_1,
			     &db_node_action_set_status_visited,
			     &db_node_action_set_status_visited,
			     &element_get_colour1, &element_get_covg_colour1, &print_no_extra_info);
	fclose(detect_vars_fptr);

	//cleanu
	hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	


	//STEP2 - detect vars in the reference colour, and mark these branches as visited, so they are ignored. Then call vars in colour1
	FILE* detect_vars_after_remv_ref_bub_fptr = fopen(detectvars_after_remv_ref_bubble_filename, "w");
	if (detect_vars_after_remv_ref_bub_fptr==NULL)
	  {
	    printf("Cannot open %s so exit\n", detectvars_after_remv_ref_bubble_filename);
	    exit(1);
	  }
	printf("Call het variants after marking off the bubbles in the ref\n");
	db_graph_detect_vars_after_marking_vars_in_reference_to_be_ignored(detect_vars_after_remv_ref_bub_fptr, max_allowed_branch_len,db_graph, 
									   &detect_vars_condition_covg_gtr_1,
									   &element_get_colour0, &element_get_covg_colour0,
									   &element_get_colour1, &element_get_covg_colour1, &print_no_extra_info);
	fclose(detect_vars_after_remv_ref_bub_fptr);
	//no need to traverse and do cleanup, as db_graph_detect_vars_after_marking_vars_in_reference_to_be_ignored does it at the end

	//STEP3 - detect hom non ref variants
	FILE* detect_vars_hom_nonref_fptr = fopen(detectvars_hom_nonref_filename, "w");
	if (detect_vars_hom_nonref_fptr==NULL)
	  {
	    printf("Cannot open %s so exit\n", detectvars_hom_nonref_filename);
	    exit(1);
	  }

	printf("About to print hom nonref calls to %s\n", detectvars_hom_nonref_filename);
	db_graph_detect_vars( detect_vars_hom_nonref_fptr, max_allowed_branch_len,db_graph, &detect_vars_condition_is_hom_nonref_with_covg_gtr_1,
			      &db_node_action_set_status_visited,  &db_node_action_set_status_visited,
			      &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, &print_no_extra_info);

	//cleanup
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	


	//STEP 4 - detect homovariants using ref-assisted trusted-path algorithm

	int min_fiveprime_flank_anchor = 3;
	int min_threeprime_flank_anchor= 3;
	int max_anchor_span =  70000;
	int length_of_arrays = 140000;
	int min_covg =2;
	int max_covg = 100000000;
	int max_expected_size_of_supernode=70000;
	

	//needs a filepointer to traverse the reference as it walks the graph
	printf("Detect polymorphic sites using reference assisted caller, using this ref file %s\n", ref_fasta);

	FILE* ref_ass_fptr = fopen(ref_assisted_filename, "w");
	if (ref_ass_fptr==NULL)
	  {
	    printf("Cnnot open %s, so exit", ref_assisted_filename);
	    exit(1);
	  }

	FILE* ref_fptr = fopen(ref_fasta, "r");
	if (ref_fptr==NULL)
	  {
	    printf("Cannot open %s, so exit", ref_fasta);
	    exit(1);
	  }

	db_graph_make_reference_path_based_sv_calls(ref_fptr, individual_edge_array, 1, individual_edge_array, 0,
						    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
						    max_expected_size_of_supernode, length_of_arrays, db_graph, ref_ass_fptr,
						    0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &db_variant_action_do_nothing,
						    &print_no_extra_info);
	fclose(ref_ass_fptr);


	//cleanup
	hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	


	printf("Finished making all calls\n");

	break;
      }
    case 10:
	{

	  //Call bubbles on NA12878 and NA19240 separately, but for each, when you find a bubble, annotate how it lies in both.
	  // Then call joint hom non-ref sites. 
	  // Then use Ref Assisted caller on each individual separately. 
	  // Then print supernodes.

	  int colour_human_ref=0;
	  int colour_na12878=1;
	  int colour_na19240=2;
	  int colour_chimp =3;
	  int colour_gorilla=4;
	  //int colour_macaca=5;
	  
	  Edges get_colour_human_ref(const dBNode* e)
	  {
	    return get_edge_copy(*e, individual_edge_array, colour_human_ref);
	  }

	  Edges get_colour_na12878(const dBNode* e)
	  {
	    return get_edge_copy(*e, individual_edge_array, colour_na12878);
	  }

	  Edges get_colour_na19240(const dBNode* e)
	  {
	    return get_edge_copy(*e, individual_edge_array, colour_na19240);
	  }

	  Edges get_colour_chimp(const dBNode* e)
	  {
	    return get_edge_copy(*e, individual_edge_array, colour_chimp);
	  }

	  Edges get_colour_gorilla(const dBNode* e)
	  {
	    return get_edge_copy(*e, individual_edge_array, colour_gorilla);
	  }

	  //  Edges get_colour_macaca(const dBNode* e)
	  //{
	  //  return get_edge_copy(*e, individual_edge_array, colour_macaca);
	  // }

	  int get_covg_human_ref(const dBNode* e)
	  {
	    return e->coverage[colour_human_ref];
	  }
	  int get_covg_na12878(const dBNode* e)
	  {
	    return e->coverage[colour_na12878];
	  }
	  int get_covg_na19240(const dBNode* e)
	  {
	    return e->coverage[colour_na19240];
	  }
	  int get_covg_union_human_ref_and_na12878(const dBNode* e)
	  {
	    return  e->coverage[colour_human_ref] + e->coverage[colour_na12878];
	  }
	  int get_covg_union_human_ref_and_na19240(const dBNode* e)
	  {
	    return  e->coverage[colour_human_ref] + e->coverage[colour_na19240];
	  }
	  
	  int get_covg_chimp(const dBNode* e)
	  {
	    return e->coverage[colour_chimp];
	  }
	  int get_covg_gorilla(const dBNode* e)
	  {
	    return e->coverage[colour_gorilla];
	  }
	  //  int get_covg_macaca(const dBNode* e)
	  //{
	  //  return e->coverage[5];
	  // }

	  int get_covg_union_ancestral_species(const dBNode* e)
	  {
	    return e->coverage[colour_chimp] + e->coverage[colour_gorilla];
	  }
	  
	  //Edges element_get_colour_human_ref(const Element* e)
	  //{
	  //  return get_edge_copy(*e, individual_edge_array,0);
	  //}

	  Edges element_get_union_human_colours(const Element* e)
	  {
	    Edges edges=0;
	    Edges edges1=get_edge_copy(*e, individual_edge_array,colour_na12878);
	    edges |= edges1;
	    Edges edges2=get_edge_copy(*e, individual_edge_array,colour_na19240);
	    edges |= edges2;
	    return edges;
	  }

	  Edges element_get_union_human_indiv_and_ref_colours(const Element* e)
	  {
	    Edges edges=0;
	    Edges edges1=get_edge_copy(*e, individual_edge_array,colour_na12878);
	    edges |= edges1;
	    Edges edges2=get_edge_copy(*e, individual_edge_array,colour_na19240);
	    edges |= edges2;
	    Edges edges3=get_edge_copy(*e, individual_edge_array,colour_human_ref);
	    edges |= edges3;

	    return edges;
	  }

	  Edges element_get_union_human_ref_and_na12878_colours(const Element* e)
	  {
	    Edges edges=0;
	    Edges edges1=get_edge_copy(*e, individual_edge_array,colour_na12878);
	    edges |= edges1;
	    Edges edges3=get_edge_copy(*e, individual_edge_array,colour_human_ref);
	    edges |= edges3;

	    return edges;
	  }
	  Edges element_get_union_human_ref_and_na19240_colours(const Element* e)
	  {
	    Edges edges=0;
	    Edges edges1=get_edge_copy(*e, individual_edge_array,colour_na19240);
	    edges |= edges1;
	    Edges edges3=get_edge_copy(*e, individual_edge_array,colour_human_ref);
	    edges |= edges3;

	    return edges;
	  }
	  
	  int get_joint_human_covg(const dBNode* e)
	  {
	    return e->coverage[colour_na12878] + e->coverage[colour_na19240];
	  }

	  int get_joint_human_indiv_and_ref_covg(const dBNode* e)
	  {
	    return e->coverage[colour_human_ref] + e->coverage[colour_na12878] + e->coverage[colour_na19240];
	  }
	  
	  boolean detect_vars_condition_is_hom_nonref_in_both_individuals(VariantBranchesAndFlanks* var)
	  {
	    return detect_vars_condition_is_hom_nonref_given_colour_funcs_for_ref_and_indiv(var, &get_covg_human_ref, &get_joint_human_covg );
	  }

	  boolean detect_vars_condition_is_hom_nonref_in_na12878(VariantBranchesAndFlanks* var)
	  {
	    return detect_vars_condition_is_hom_nonref_given_colour_funcs_for_ref_and_indiv(var, &get_covg_human_ref, &get_covg_na12878 );
	  }

	  boolean detect_vars_condition_is_hom_nonref_in_na19240(VariantBranchesAndFlanks* var)
	  {
	    return detect_vars_condition_is_hom_nonref_given_colour_funcs_for_ref_and_indiv(var, &get_covg_human_ref, &get_covg_na19240 );
	  }
	  
	  int max(int a, int b)
	  {
	    if (a>b)
	      {
		return a;
	      }
	    else
	      {
		return b;
	      }
	  }


	  void print_zygosity_string(zygosity z, FILE* fptr)
	  {
	    if (z==het)
	      {
		fprintf(fptr, "HET");
	      }
	    else if (z==hom_one)
	      {
		fprintf(fptr, "HOM branch1");
	      }
	    else if (z==hom_other)
	      {
		fprintf(fptr, "HOM branch2");
	      }
	    else
	      {
		fprintf(fptr, "BOTH_ALLELES_ABSENT");
	      }
	  }
	  


	  void print_extra_supernode_info(dBNode** node_array, Orientation* or_array, int len, FILE* fout)
	  {
	    fprintf(fout, "Mult in  hum ref:\n");
	    int i;
	    for (i=0; i<len; i++)
	      {
		if (node_array[i]!=NULL)
		  {
		    fprintf(fout, "%d ", node_array[i]->coverage[colour_human_ref]);
		  }
		else
		  {
		    fprintf(fout, "0 ");
		  }
	      }
	    fprintf(fout, "\n");
	    fprintf(fout, "Covg in NA12878:\n");
	    for (i=0; i<len; i++)
	      {
		if (node_array[i]!=NULL)
		  {
		    fprintf(fout, "%d ", node_array[i]->coverage[colour_na12878]);
		  }
		else
		  {
		    fprintf(fout, "0 ");
		  }
	      }
	    fprintf(fout, "\n");
	    fprintf(fout, "Covg in NA19240:\n");
	    for (i=0; i<len; i++)
	      {
		if (node_array[i]!=NULL)
		  {
		    fprintf(fout, "%d ", node_array[i]->coverage[colour_na19240]);
		  }
		else
		  {
		    fprintf(fout, "0 ");
		  }
	      }
	    fprintf(fout, "\n");
	    fprintf(fout, "Mult in chimp:\n");
	    for (i=0; i<len; i++)
	      {
		if (node_array[i]!=NULL)
		  {		    
		    fprintf(fout, "%d ", node_array[i]->coverage[colour_chimp]);
		  }
		else
		  {
		    fprintf(fout, "0 ");
		  }
	      }
	    fprintf(fout, "\n");
	    fprintf(fout, "Mult in gorilla:\n");
	    for (i=0; i<len; i++)
	      {
		if (node_array[i]!=NULL)
		  {
		    fprintf(fout, "%d ", node_array[i]->coverage[colour_gorilla]);
		  }
		else
		  {
		    fprintf(fout, "0 ");
		  }
	      }
	    fprintf(fout, "\n\n");
		    
	  }


	  
	  void print_extra_info(VariantBranchesAndFlanks* var, FILE* fout)
	  {
	    // determine ancestral allele by comparing with chimp, gorilla
	    fprintf(fout, "\n");
	    
	    WhichAllele which_allele_matches_chimp;
	    WhichAllele which_allele_matches_gorilla;
	    boolean precisely_one_allele_matches_chimp=db_variant_precisely_one_allele_is_in_given_func_of_colours(var, &get_colour_chimp, db_graph, &which_allele_matches_chimp);
	    boolean precisely_one_allele_matches_gorilla=db_variant_precisely_one_allele_is_in_given_func_of_colours(var, &get_colour_gorilla, db_graph, &which_allele_matches_gorilla); ;

	    if (precisely_one_allele_matches_chimp==true)
	      {
		fprintf(fout, "ANCESTRAL ALLELE: ");
		if (which_allele_matches_chimp==allele_one)
		  {
		    fprintf(fout, "branch1 matches chimp\n");//and branch2 does not
		  }
		else
		  {
		    fprintf(fout, "branch2 matches chimp\n");//and branch1 does not
		  }
	      }
	    else if (precisely_one_allele_matches_gorilla==true)
	      {
		fprintf(fout, "ANCESTRAL ALLELE: ");
		if (which_allele_matches_gorilla==allele_one)
		  {
		    fprintf(fout, "branch1 matches gorilla\n");//and branch2 does not
		  }
		else
		  {
		    fprintf(fout, "branch2 matches gorilla\n");//and branch1 does not
		  }
	      }
	    else
	      {
		fprintf(fout, "ANCESTRAL ALLELE: unknown\n");
	      }
	    
	    // determine which individual this variant is on, and print
	    zygosity zygo_in_na12878 = db_variant_get_zygosity_in_given_func_of_colours(var, get_colour_na12878, db_graph);
	    zygosity zygo_in_na19240 = db_variant_get_zygosity_in_given_func_of_colours(var, get_colour_na19240, db_graph);
	    fprintf(fout, "NA12878:");
	    print_zygosity_string(zygo_in_na12878, fout);
	    fprintf(fout,"\n");
	    fprintf(fout, "NA19240:");
	    print_zygosity_string(zygo_in_na19240, fout);
	    fprintf(fout,"\n");
	    
	    //determine which allele matches reference
	    if (does_this_path_exist_in_this_colour(var->one_allele, var->one_allele_or, var->len_one_allele, get_colour_human_ref, db_graph)==true)
	      {
		fprintf(fout, "REF_ALLELE:branch1\n");
	      }
	    else if (does_this_path_exist_in_this_colour(var->other_allele, var->other_allele_or, var->len_other_allele, get_colour_human_ref, db_graph)==true)
	      {
		fprintf(fout, "REF_ALLELE:branch2\n");
	      }
	    else
	      {
		fprintf(fout, "BOTH_ALLELES_NOVEL\n");
	      }
	    //print coverages:
	    fprintf(fout, "branch1 coverages\n");
	    print_extra_supernode_info(var->one_allele, var->one_allele_or, var->len_one_allele, fout);
	    fprintf(fout, "branch2 coverages\n");
	    print_extra_supernode_info(var->other_allele, var->other_allele_or, var->len_other_allele, fout);
	    fprintf(fout, "\n\n");
	  }

	  
	  
	  int max_allowed_branch_len=5000;


	  //initialise output files and chrom ref ptrs for the trusted path caller:

	char** ref_chroms = malloc( sizeof(char*) * 24); //one for each chromosome, ignoring haplotypes like MHC, and the Y chrom
	if (ref_chroms==NULL)
	  {
	    printf("OOM. Give up can't even allocate space for the names of the ref chromosome files\n");
	    exit(1);
	  }
	int i;
	for (i=0; i< 24; i++)
	  {
	    ref_chroms[i] = malloc(sizeof(char)*150); //filenames including path are about 100 characters. 50 characters of leeway
	    if (ref_chroms[i]==NULL)
	      {
		printf("OOM. Giveup can't even allocate space for the names of the ref chromosome file i = %d\n",i);
		exit(1);
	      }
	  }
	
	set_ref_chromosome_file_pointers(ref_chroms, 24);
	
	
	// **** set up one output file per chromosome ***** //
	
	char** na12878_output_files = malloc( sizeof(char*) * 24); 
	char** na19240_output_files = malloc( sizeof(char*) * 24); 
	if ((na12878_output_files==NULL)||(na19240_output_files==NULL))
	  {
	    printf("OOM. Give up can't even allocate space for the names of the output  files \n");
	    exit(1);
	  }
	
	for (i=0; i< 24; i++)
	  {
	    na12878_output_files[i] = malloc(sizeof(char)*200); 
	    na19240_output_files[i] = malloc(sizeof(char)*200); 
	    if ((na12878_output_files[i]==NULL)||(na19240_output_files[i]==NULL))
	      {
		printf("OOM. Giveup can't even allocate space for the names of the ref chromosome file i = %d\n",i);
		exit(1);
	      }
	  }
	
	sprintf(na12878_output_files[0],"na12878_ref_assisted_variants_called_in_MT");
	sprintf(na19240_output_files[0],"na19240_ref_assisted_variants_called_in_MT");
	
	for (i=1; i<23; i++)
	  {
	    sprintf(na12878_output_files[i], "na12878_ref_assisted_variants_called_in_chrom_%i", i);
	    sprintf(na19240_output_files[i], "na19240_ref_assisted_variants_called_in_chrom_%i", i);
	  }
	sprintf(na12878_output_files[23],"na12878_ref_assisted_variants_called_in_chrom_X");
	sprintf(na19240_output_files[23],"na19240_ref_assisted_variants_called_in_chrom_X");
	//na12878 and na19240 are women,so don't look for variants in the Y chromosome!!!

	//***************************
	//end of initialisation
	// **************************
	  
	
	//comment out the bubble calls temporarily - want to re-run js the ref-assisted calls and supernodes


	/*
	//CALL het bubbles in na12878

	char na12878_bubbles[100];
	sprintf(na12878_bubbles, "na12878_bubbles");
	FILE* detect_vars_after_remv_ref_bub_fptr = fopen(na12878_bubbles, "w");
	if (detect_vars_after_remv_ref_bub_fptr==NULL)
	  {
	    printf("Cannot open %s so exit\n", na12878_bubbles);
	    exit(1);
	  }

	printf("Call het variants on na12878 alone, after marking off the bubbles in the ref");
	db_graph_detect_vars_after_marking_vars_in_reference_to_be_ignored(detect_vars_after_remv_ref_bub_fptr, max_allowed_branch_len,db_graph, 
									   &detect_vars_condition_always_true,
									   &get_colour_human_ref, &get_covg_human_ref,
									   &get_colour_na12878, &get_covg_na12878,
									   &print_extra_info);
	fclose(detect_vars_after_remv_ref_bub_fptr);
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	
       

	//do it again, just na19240
	char na19240_bubbles[100];
	sprintf(na19240_bubbles, "na19240_bubbles");
	detect_vars_after_remv_ref_bub_fptr = fopen(na19240_bubbles, "w");
	if (detect_vars_after_remv_ref_bub_fptr==NULL)
	  {
	    printf("Cannot open %s so exit\n", na19240_bubbles);
	    exit(1);
	  }
	printf("Call het variants on na19240 alone\n");
	db_graph_detect_vars_after_marking_vars_in_reference_to_be_ignored(detect_vars_after_remv_ref_bub_fptr, max_allowed_branch_len,db_graph, 
									   &detect_vars_condition_always_true,
									   &get_colour_human_ref, &get_covg_human_ref,
									   &get_colour_na19240, &get_covg_na19240,
									   &print_extra_info);
	fclose(detect_vars_after_remv_ref_bub_fptr);
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	



	//STEP3 - detect hom non ref variants jointly - ie these are hom-non-ref in both individuals
	char joint_hom_nonref_bubbles[100];
	sprintf(joint_hom_nonref_bubbles, "joint_hom_nonref_bubbles");
	FILE* detect_vars_hom_nonref_fptr = fopen(joint_hom_nonref_bubbles, "w");
	if (detect_vars_hom_nonref_fptr==NULL)
	  {
	    printf("Cannot open %s so exit\n", joint_hom_nonref_bubbles);
	    exit(1);
	  }

	printf("About to print hom nonref calls to %s\n", joint_hom_nonref_bubbles);
	db_graph_detect_vars( detect_vars_hom_nonref_fptr, max_allowed_branch_len,db_graph, &detect_vars_condition_is_hom_nonref_in_both_individuals,
			      &db_node_action_set_status_visited,  &db_node_action_set_status_visited,
			      &element_get_union_human_indiv_and_ref_colours, &get_joint_human_indiv_and_ref_covg, &print_extra_info);
	fclose(detect_vars_hom_nonref_fptr);
	//cleanup
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	
	*/


	//Detect hom non ref in each individual separately
	char hom_nonref_na12878[100];
	sprintf(hom_nonref_na12878, "hom_nonref_na12878");
	FILE* homnonref_fptr = fopen(hom_nonref_na12878, "w");

	if (homnonref_fptr==NULL)
	  {
	    printf("Cannot open %s so exit\n", hom_nonref_na12878);
	    exit(1);
	  }

	db_graph_detect_vars( homnonref_fptr, max_allowed_branch_len,db_graph, &detect_vars_condition_is_hom_nonref_in_na12878,
			      &db_node_action_set_status_visited,  &db_node_action_set_status_visited,
			      &element_get_union_human_ref_and_na12878_colours, &get_covg_union_human_ref_and_na12878, &print_extra_info);
	fclose(homnonref_fptr);
	//cleanup
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	


	char hom_nonref_na19240[100];
	sprintf(hom_nonref_na19240, "hom_nonref_na19240");
	homnonref_fptr = fopen(hom_nonref_na19240, "w");

	if (homnonref_fptr==NULL)
	  {
	    printf("Cannot open %s so exit\n", hom_nonref_na19240);
	    exit(1);
	  }

	db_graph_detect_vars( homnonref_fptr, max_allowed_branch_len,db_graph, &detect_vars_condition_is_hom_nonref_in_na19240,
			      &db_node_action_set_status_visited,  &db_node_action_set_status_visited,
			      &element_get_union_human_ref_and_na19240_colours, &get_covg_union_human_ref_and_na19240, &print_extra_info);

	fclose(homnonref_fptr);
	//cleanup
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	

	/*


	//STEP 4 - detect variants using ref-assisted trusted-path algorithm

	int min_fiveprime_flank_anchor = 3;
	int min_threeprime_flank_anchor= 3;
	int max_anchor_span =  70000;
	int length_of_arrays = 140000;
	int min_covg =1;
	int max_covg = 100000000;
	int max_expected_size_of_supernode=70000;
	
	

	printf("Detect polymorphic sites in NA12878 using reference assisted caller, putting no conditions on variant\n");
	
	//ignore mitochondrion for now, so start with i=1
	for (i=1; i<24; i++) 
	  {
	    printf("Call SV comparing NA12878 with chromosome %s\n", ref_chroms[i]);
	    
	    FILE* chrom_fptr = fopen(ref_chroms[i], "r");
	    if (chrom_fptr==NULL)
	      {
		printf("Cannot open %s \n", ref_chroms[i]);
		exit(1);
	      }
	    
	    FILE* out_fptr = fopen(na12878_output_files[i], "w");
	    if (out_fptr==NULL)
	      {
		printf("Cannot open %s for output\n", na12878_output_files[i]);
		exit(1);
	      }
	    
	    db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, colour_na12878, 
							individual_edge_array, colour_human_ref,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, db_graph, out_fptr,
							0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &db_variant_action_do_nothing,
							&print_extra_info);
	    fclose(chrom_fptr);
	    fclose(out_fptr);
	  }

	//cleanup
	hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
	
	
	//Now do the same for NA19240
	
	printf("Detect polymorphic sites in NA19240 using reference assisted caller, putting no conditions on variant\n");
	
	for (i=1; i<24; i++) 
	  {
	    printf("Call SV comparing NA19240 with chromosome %s\n", ref_chroms[i]);
	    
	    FILE* chrom_fptr = fopen(ref_chroms[i], "r");
	    if (chrom_fptr==NULL)
	      {
		printf("Cannot open %s \n", ref_chroms[i]);
		exit(1);
	      }
	    
	    FILE* out_fptr = fopen(na19240_output_files[i], "w");
	    if (out_fptr==NULL)
	      {
		printf("Cannot open %s for output\n", na19240_output_files[i]);
		exit(1);
	      }
	    
	    db_graph_make_reference_path_based_sv_calls(chrom_fptr, individual_edge_array, colour_na19240, 
							individual_edge_array, colour_human_ref,
							min_fiveprime_flank_anchor, min_threeprime_flank_anchor, max_anchor_span, min_covg, max_covg, 
							max_expected_size_of_supernode, length_of_arrays, db_graph, out_fptr,
							0, NULL, NULL, NULL, NULL, NULL, &make_reference_path_based_sv_calls_condition_always_true, &db_variant_action_do_nothing,
							&print_extra_info);
	    fclose(chrom_fptr);
	    fclose(out_fptr);
	  }


	
	//cleanup
	hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
	*/	

	printf("Finished making all calls\n");


	/*
	printf("Now print annotated supernodes of  NA12878\n");
	db_graph_print_supernodes_for_specific_person_or_pop("na12878_sups", "na12878_sings", 3000, db_graph, individual_edge_array,colour_na12878, &print_extra_supernode_info);
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	

	printf("Now print annotated supernodes of  NA19240\n");
	db_graph_print_supernodes_for_specific_person_or_pop("na19240_sups", "na19240_sings", 3000, db_graph, individual_edge_array,colour_na19240, &print_extra_supernode_info);
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	
	printf("Now print annotated SHARED supernodes of  NA12878 and NA19240\n");
	db_graph_print_supernodes_defined_by_func_of_colours("shared_sups", "shared_sings", 3000, 
							     db_graph, &element_get_union_human_colours, &get_joint_human_covg,
							     &print_extra_supernode_info);
	hash_table_traverse(&db_node_action_set_status_none, db_graph);	
	*/
	printf("FINISHED ALL\n");
	break;
	}
    case 11:
      {
	printf("Get covg distribution in colour 0\n");
	db_graph_get_covg_distribution(db_graph, individual_edge_array, 0);

	break;
      }
      
      

    }

  hash_table_free(&db_graph);

  return 0;
}
