#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <string.h>
#include <internal_oxford.h>
#include <cmd_line.h>

//#include <internal_tgac.h>


//TODO

//3. trusted-path calls
//4. handling illumina offset




void run_bubble_calls(CmdLine* cmd_line, int which, dBGraph* db_graph, 
		      void (*print_appropriate_extra_var_info)(VariantBranchesAndFlanks* var, FILE* fp),
		      Edges(*get_col_ref) (const dBNode* e),
		      int (*get_cov_ref)(const dBNode* e)
		      )
{
  
      printf("Detecting bubbles between first colour list: ");
      int k;
      if (which==1)
	{
	  for (k=0; k<cmd_line->num_colours_in_detect_bubbles1_first_colour_list; k++)
	    {
	      printf("%d, ", cmd_line->detect_bubbles1_first_colour_list[k]);
	    }
	}
      else
	{
	  for (k=0; k<cmd_line->num_colours_in_detect_bubbles2_first_colour_list; k++)
	    {
	      printf("%d, ", cmd_line->detect_bubbles2_first_colour_list[k]);
	    }
	}
      printf(" and second colour list:");
      if (which==1)
	{
	  for (k=0; k<cmd_line->num_colours_in_detect_bubbles1_second_colour_list; k++)
	    {
	      printf("%d, ", cmd_line->detect_bubbles1_second_colour_list[k]);
	    }
	}
      else
	{
	  for (k=0; k<cmd_line->num_colours_in_detect_bubbles2_second_colour_list; k++)
	    {
	      printf("%d, ", cmd_line->detect_bubbles2_second_colour_list[k]);
	    }
	  printf("\n");
	}


      FILE* fp;

      if (which==1)
	{
	  fp = fopen(cmd_line->output_detect_bubbles1, "w");
	}
      else
	{
	  fp = fopen(cmd_line->output_detect_bubbles2, "w");
	}

      if (fp==NULL)
	{
	  if (which==1)
	    {
	      printf("Cannot open %s. Exit.", cmd_line->output_detect_bubbles1);
	    }
	  else
	    {
	      printf("Cannot open %s. Exit.", cmd_line->output_detect_bubbles2);
	    }
	  exit(1);
	}

      if (cmd_line->using_ref==false)
	{
	  if (which==1)
	    {
	      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_supernode,db_graph, 
							  cmd_line->detect_bubbles1_first_colour_list, 
							  cmd_line->num_colours_in_detect_bubbles1_first_colour_list,
							  cmd_line->detect_bubbles1_second_colour_list, 
							  cmd_line->num_colours_in_detect_bubbles1_second_colour_list,
							  &detect_vars_condition_always_true, print_appropriate_extra_var_info,
							  false, NULL, NULL);
	    }
	  else
	    {
	      printf("Zam here - about to call vars\n");
	      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_supernode,db_graph, 
							  cmd_line->detect_bubbles2_first_colour_list, 
							  cmd_line->num_colours_in_detect_bubbles2_first_colour_list,
							  cmd_line->detect_bubbles2_second_colour_list, 
							  cmd_line->num_colours_in_detect_bubbles2_second_colour_list,
							  &detect_vars_condition_always_true, print_appropriate_extra_var_info,
							  false, NULL, NULL);
	    }
	}
      else
	{
	  printf("(First exclude bubbles from ref colour %d) \n", cmd_line->ref_colour);
	  if (which==1)
	    {
	      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_supernode,db_graph, 
							  cmd_line->detect_bubbles1_first_colour_list, 
							  cmd_line->num_colours_in_detect_bubbles1_first_colour_list,
							  cmd_line->detect_bubbles1_second_colour_list, 
							  cmd_line->num_colours_in_detect_bubbles1_second_colour_list,
							  &detect_vars_condition_always_true,
							  print_appropriate_extra_var_info,
							  true, get_col_ref, get_cov_ref);
	    }
	  else
	    {
	      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_supernode,db_graph, 
							  cmd_line->detect_bubbles2_first_colour_list, 
							  cmd_line->num_colours_in_detect_bubbles2_first_colour_list,
							  cmd_line->detect_bubbles2_second_colour_list, 
							  cmd_line->num_colours_in_detect_bubbles2_second_colour_list,
							  &detect_vars_condition_always_true,
							  print_appropriate_extra_var_info,
							  true, get_col_ref, get_cov_ref);
	    }
	}
}






int main(int argc, char **argv){

  CmdLine cmd_line = parse_cmdline(argc,argv,sizeof(Element));

  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL;
  short kmer_size;


  //***************************************************************************
  //define local functions
  //***************************************************************************

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

  void print_appropriate_extra_variant_info(VariantBranchesAndFlanks* var, FILE* fp)
  {
    if (cmd_line.print_colour_coverages==true)
      {
	print_standard_extra_info(var, fp);
      }
    else
      {
	print_no_extra_info(var, fp);
      }
  }

  void print_appropriate_extra_supernode_info(dBNode** node_array, Orientation* or_array, int len, FILE* fout)
  {
    if (cmd_line.print_colour_coverages==true)
      {
	print_standard_extra_supernode_info(node_array, or_array, len, fout);
      }
    else
      {
	print_no_extra_supernode_info(node_array, or_array, len, fout);
      }
  }


  Edges get_colour_ref(const dBNode* e)
  {

    if (cmd_line.using_ref==false)
      {
	printf("Do not call get_colour_ref when --using_ref was not specified. Exiting - coding error\n");
	exit(1);
      }

    if ( (cmd_line.ref_colour<0) || (cmd_line.ref_colour>NUMBER_OF_INDIVIDUALS_PER_POPULATION-1) )
      {
	printf("Calling get_colour_ref, but the reference colour %d has not been specified, or is > than the compile-time limit, of %d\n", 
	       cmd_line.ref_colour, NUMBER_OF_INDIVIDUALS_PER_POPULATION-1);
	exit(1);
      }
    Edges ed = get_edge_copy(*e, individual_edge_array, cmd_line.ref_colour);
    return ed;
  }

  int get_covg_ref(const dBNode* e)
  {
    if (cmd_line.using_ref==false)
      {
	printf("Do not call get_coverage_ref when --using_ref was not specified. Exiting - coding error\n");
      }
    
    return e->coverage[cmd_line.ref_colour];
  }

  //***************************************************************************
  // end local functions
  //***************************************************************************



  
  //set hash table variables:
  kmer_size        = cmd_line.kmer_size;
  hash_key_bits    = cmd_line.number_of_buckets_bits; //number of buckets: 2^hash_key_bits
  bucket_size      = cmd_line.bucket_size;

  //zam - these are from Mario, want to check
  //int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
  
  printf("Maximum k-mer size (compile-time setting): %ld\n", (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1);
  
  if (cmd_line.kmer_size> (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1){
    printf("k-mer size is too big [%i]!",cmd_line.kmer_size);
    exit(1);
  }
  printf("Actual K-mer size: %d\n", cmd_line.kmer_size);

  
  //Create the de Bruijn graph/hash table
  int max_retries=15;
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  printf("Hash table created, number of buckets: %d\n",1 << hash_key_bits);


  // input data:
  if (cmd_line.input_seq==true)
    {
      if (strcmp(cmd_line.se_list, "")!=0)
	{
	  fprintf(stdout,"Input file of single ended data filenames: %s\n",
		  cmd_line.se_list); 
	}
      else
	{
	  printf("No SE data\n");
	}
      if (strcmp(cmd_line.pe_list_lh_mates, "") !=0)
	{
	  fprintf(stdout,"Input file of paired end data: %s, and %s, \n",
		  cmd_line.pe_list_lh_mates, cmd_line.pe_list_rh_mates); 
	}
      else
	{
	  printf("No paired-end data\n");
	}
      if (cmd_line.quality_score_threshold>0)
	{
	  fprintf(stdout,"quality cut-off: %i\n",cmd_line.quality_score_threshold);
	}
      else
	{
	  fprintf(stdout, "No quality filtering\n");
	}
      if (cmd_line.remove_pcr_dups==true)
	{
	  printf("Removing duplicates from single-ended files single-endedly - ie if any read starts with same k-mer as a previous read\n");
	  printf("Removing duplicates from the paired end files when both mates start with the same kmer\n");
	}
      else
	{
	  printf("No PCR duplicate removal\n");
	}
      if (cmd_line.cut_homopolymers==true)
	{
	  printf("Breaking reads at homopolymer runs of length %d or greater. Read restarts at first base after the run\n", cmd_line.homopolymer_limit);
	}

      boolean there_is_se_data=false;
      boolean there_is_pe_data=false;
      if (strcmp(cmd_line.se_list, "")!=0)
	{
	  there_is_se_data=true;
	}
      if (strcmp(cmd_line.pe_list_lh_mates, "") != 0)
	{
	  there_is_pe_data=true;
	}
      
      //note, if we load fasta/fastq data it always goes into colour 0;
      load_se_and_pe_filelists_into_graph_of_specific_person_or_pop(there_is_se_data, there_is_pe_data, 
								    cmd_line.se_list, cmd_line.pe_list_lh_mates, cmd_line.pe_list_rh_mates,
								    cmd_line.quality_score_threshold, cmd_line.remove_pcr_dups, cmd_line.remove_pcr_dups,
								    cmd_line.cut_homopolymers, cmd_line.homopolymer_limit, cmd_line.format_of_input_seq,
								    cmd_line.max_read_length, 0, db_graph);
      
    }
  else
    {
      //if there is a multicolour binary, load that in first
      
      int first_colour_data_starts_going_into=0;

      if (cmd_line.input_multicol_bin==true)
	{
	  int kmers_loaded = load_multicolour_binary_from_filename_into_graph(cmd_line.multicolour_bin,db_graph, &first_colour_data_starts_going_into);
	  printf("Loaded the multicolour binary %s, and got %d kmers\n", cmd_line.multicolour_bin, kmers_loaded);
	}

      if (cmd_line.input_colours==true)
	{
	  printf("List of colours: %s (contains one filelist per colour). Load data into consecutive colours starting at %d\n", 
		 cmd_line.colour_list, first_colour_data_starts_going_into);
	  load_population_as_binaries_from_graph(cmd_line.colour_list, db_graph);
	  //printf("Finished loading single_colour binaries\n");
	}


    }
      
  printf("Total kmers in table: %qd\n", hash_table_get_unique_kmers(db_graph));	  


  
  // Error Correction actions
  if (cmd_line.remove_seq_errors==true)
    {
      printf("Remove nodes that look like sequencing errors. Clip tips first\n");
      db_graph_clip_tips_for_specific_person_or_pop(db_graph, individual_edge_array, 0);
      
      printf("Then remove low coverage nodes (<= 1) based on topology as well as covg  -  must look like induced by a single base error) \n");
      db_graph_remove_errors_considering_covg_and_topology(1,db_graph, &element_get_covg_union_of_all_covgs, &element_get_colour_union_of_all_colours,
							   &apply_reset_to_all_edges, &apply_reset_to_all_edges_2,0);//last argument is colour of ref, to be protected.

    }
  else
    {
      if (cmd_line.clip_tips==true)
	{
	  printf("Clipping tips from colour0 graph.\n");
	  db_graph_clip_tips_for_specific_person_or_pop(db_graph, individual_edge_array, 0);
	}
      if (cmd_line.remove_low_coverage_nodes==true)
	{
	  printf("Removing nodes with coverage <= %d (from union graph of all colours)\n", cmd_line.node_coverage_threshold);
	  db_graph_remove_low_coverage_nodes(cmd_line.node_coverage_threshold,db_graph, &element_get_covg_union_of_all_covgs, &element_get_colour_union_of_all_colours,
					     &apply_reset_to_all_edges, &apply_reset_to_all_edges_2);
	}

    }

  if (cmd_line.dump_binary==true)
    {
      if (cmd_line.input_seq==true)
	{
	  //dump single colour
	  printf("Input data was fasta/q, so dump single colour binary file: %s\n", cmd_line.output_binary_filename);
	  db_graph_dump_single_colour_binary_of_colour0(cmd_line.output_binary_filename, &db_node_check_status_not_pruned,db_graph);
	}
      else
	{
	  //printf("Dump multicolour binary with %d colours (compile-time setting), of which input data .......\n");
	  //db_graph_dump_binary(cmd_line.output_binary_filename, &db_node_check_status_not_pruned,db_graph);
	}
    }

  if (cmd_line.print_contig_fasta==true)
    {
      printf("Print contigs(supernodes) in the graph created by the union of all colours.\n");
      
      db_graph_print_supernodes_defined_by_func_of_colours(cmd_line.output_supernodes, "", cmd_line.max_supernode,
							   db_graph, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, 
							   &print_appropriate_extra_supernode_info);

    }


  // DETECT BUBBLES

  if (cmd_line.detect_bubbles1==true)
    {
      run_bubble_calls(&cmd_line, 1, db_graph, &print_appropriate_extra_variant_info, &get_colour_ref, &get_covg_ref);

      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
    }

  //second detect bubbles
  if (cmd_line.detect_bubbles2==true)
    {
      run_bubble_calls(&cmd_line, 2, db_graph, &print_appropriate_extra_variant_info, &get_colour_ref, &get_covg_ref);

      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      

    }

  
  hash_table_free(&db_graph);

  return 0;
}
