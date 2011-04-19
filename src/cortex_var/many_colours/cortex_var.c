/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo 
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
#include <dB_graph_population.h>
#include <string.h>
#include <cmd_line.h>
#include <time.h>
#include <graph_info.h>
#include <db_differentiation.h>

void timestamp();



long long calculate_mean(long long* array, long long len)
{
  long long sum=0;
  long long num=0;
  long long i;
  for (i=0; i<len; i++)
    {
      sum += i*array[i];
      num += array[i];
    }
  return  (sum/num);
}

void run_pd_calls(CmdLine* cmd_line, dBGraph* db_graph, 
		  void (*print_some_extra_var_info)(VariantBranchesAndFlanks* var, FILE* fp))
{
  printf("Calling variants using Path Divergence Caller.\n");
  printf("Calls are made between the reference path as specified by the fasta in %s\n", cmd_line->ref_chrom_fasta_list);
  printf(" and the sample, which is the union of the following colour(s): ");
  int i;
  for (i=0; i< cmd_line->num_colours_in_pd_colour_list-1; i++)
    {
      printf("%d,", cmd_line->pd_colour_list[i]);
    }
  printf("%d\n", cmd_line->pd_colour_list[cmd_line->num_colours_in_pd_colour_list-1]);
  printf("The reference colour is %d\n", cmd_line->ref_colour);

  //this will also check that all the ref chrom fasta files exist
  int num_ref_chroms = get_number_of_files_and_check_existence_from_filelist(cmd_line->ref_chrom_fasta_list);

  char** ref_chroms = malloc( sizeof(char*) * num_ref_chroms);
  if (ref_chroms==NULL)
    {
      printf("OOM. Give up can't even allocate space for the names of the ref chromosome files\n");
      exit(1);
    }

  for (i=0; i< num_ref_chroms; i++)
    {
      ref_chroms[i] = malloc(sizeof(char)*500);
      if (ref_chroms[i]==NULL)
	{
	  printf("OOM. Giveup can't even allocate space for the names of the ref chromosome file i = %d\n",i);
	  exit(1);
	}
    }

  get_filenames_from_list(cmd_line->ref_chrom_fasta_list, ref_chroms, num_ref_chroms);

  //now set up output files

  char** output_files = malloc( sizeof(char*) * num_ref_chroms); //one for each chromosome
  if (output_files==NULL)
    {
      printf("OOM. Give up can't even allocate space for the names of the output  files \n");
      exit(1);
    }
	
  for (i=0; i< num_ref_chroms; i++)
    {
      output_files[i] = malloc(sizeof(char)*1000); 
      if (output_files[i]==NULL)
	{
	  printf("OOM. Giveup can't even allocate space for the names of the ref chromosome file i = %d\n",i);
	  exit(1);
	}
    }
  
  
  for (i=0; i<num_ref_chroms; i++)
    {
      sprintf(output_files[i], "%s_pd_chr_%i", cmd_line->path_divergence_caller_output_stub, i+1);
    }


  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= cmd_line->kmer_size;
  int max_anchor_span =  cmd_line -> max_var_len;
  int length_of_arrays = 2*max_anchor_span;
  int min_covg =1;
  int max_covg = 10000000;//this is ignored. will be changing API
  int max_expected_size_of_supernode=cmd_line -> max_var_len;
	
      
  for (i=0; i<num_ref_chroms; i++) 
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
      
      
      db_graph_make_reference_path_based_sv_calls_given_list_of_colours_for_indiv(cmd_line->pd_colour_list, cmd_line->num_colours_in_pd_colour_list,
										  chrom_fptr, cmd_line->ref_colour,
										  min_fiveprime_flank_anchor, min_threeprime_flank_anchor, 
										  max_anchor_span, min_covg, max_covg, 
										  max_expected_size_of_supernode, length_of_arrays, db_graph, out_fptr,
										  0, NULL, NULL, NULL, NULL, NULL, 
										  &make_reference_path_based_sv_calls_condition_always_true_in_subgraph_defined_by_func_of_colours, 
										  &db_variant_action_do_nothing,
										  print_some_extra_var_info);
      
      
      fclose(chrom_fptr);
      fclose(out_fptr);
    }
  //cleanup
  for(i=0; i<num_ref_chroms; i++)
    {
      free(ref_chroms[i]);
      free(output_files[i]);
    }
  free(ref_chroms);
  free(output_files);

  

}

void run_bubble_calls(CmdLine* cmd_line, int which, dBGraph* db_graph, 
		      void (*print_appropriate_extra_var_info)(VariantBranchesAndFlanks* var, FILE* fp),
		      Edges(*get_col_ref) (const dBNode* e),
		      int (*get_cov_ref)(const dBNode* e)
		      )
{
  
      printf("Detecting bubbles between the union of this set of colours: ");
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
      printf(" and the union of this set of colours: ");
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
	      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_var_len,db_graph, 
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
	      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_var_len,db_graph, 
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
	      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_var_len,db_graph, 
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
	      db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_var_len,db_graph, 
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

  timestamp();
  printf("Starting Cortex\n");

  CmdLine cmd_line = parse_cmdline(argc,argv,sizeof(Element));

  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL;
  short kmer_size;


  //***************************************************************************
  //define local functions
  //***************************************************************************
  /*
  void apply_reset_to_all_edges(dBNode* node, Orientation or, Nucleotide nuc)
  {
    int j;
    for (j=0; j<NUMBER_OF_COLOURS; j++)
      {
	reset_one_edge(node, or, nuc, individual_edge_array, j);
      }
  }
  void apply_reset_to_all_edges_2(dBNode* node )
  {
    int j;
    for (j=0; j<NUMBER_OF_COLOURS; j++)
      {
	      db_node_reset_edges(node, individual_edge_array, j);
      }
  }
  */

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

    if ( (cmd_line.ref_colour<0) || (cmd_line.ref_colour>NUMBER_OF_COLOURS-1) )
      {
	printf("Calling get_colour_ref, but the reference colour %d has not been specified, or is > than the compile-time limit, of %d\n", 
	       cmd_line.ref_colour, NUMBER_OF_COLOURS-1);
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

  int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
  int max_kmer_size = (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1;
  int min_kmer_size = ((NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1) * sizeof(bitfield_of_64bits) * 4) + 1;

  if (number_of_bitfields != NUMBER_OF_BITFIELDS_IN_BINARY_KMER){
    printf("K-mer %i  is not in current range of kmers [%i - %i] required for this executable.!\n",kmer_size,min_kmer_size,max_kmer_size);
    exit(0);
  }
  
  printf("Maximum k-mer size (compile-time setting): %ld\n", (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1);
  
  if (cmd_line.kmer_size> (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1){
    printf("k-mer size is too big [%i]!",cmd_line.kmer_size);
    exit(1);
  }
  printf("Actual K-mer size: %d\n", cmd_line.kmer_size);

  
  //Create the de Bruijn graph/hash table
  int max_retries=15;
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  if (db_graph==NULL)
    {
      printf("Giving up - unable to allocate memory for the hash table\n");
      exit(1);
    }
  printf("Hash table created, number of buckets: %d\n",1 << hash_key_bits);
  GraphInfo db_graph_info;
  graph_info_initialise(&db_graph_info);



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
      long long bases_parsed=0;
      long long bases_pass_filters_and_loaded=0;
      //get read-len distribution (after filters):
      long long* readlen_distrib=(long long*) malloc(sizeof(long long) * (cmd_line.max_read_length+1));
      long long** readlen_distrib_ptrs = (long long**) malloc(sizeof(long long*) * (cmd_line.max_read_length+1));

      if (readlen_distrib==NULL)
	{
	  printf("Unable to malloc array to hold readlen distirbution!Exit.\n");
	  exit(1);
	}
      int i;
      for (i=0; i<=cmd_line.max_read_length; i++)
	{
	  readlen_distrib[i]=0;
	  readlen_distrib_ptrs[i]=&readlen_distrib[i];
	}

      timestamp();

      load_se_and_pe_filelists_into_graph_of_specific_person_or_pop(there_is_se_data, there_is_pe_data, 
								    cmd_line.se_list, cmd_line.pe_list_lh_mates, cmd_line.pe_list_rh_mates,
								    &bases_parsed, &bases_pass_filters_and_loaded,readlen_distrib_ptrs,
								    cmd_line.quality_score_threshold, false, cmd_line.remove_pcr_dups,
								    cmd_line.cut_homopolymers, cmd_line.homopolymer_limit, cmd_line.quality_score_offset,
								    cmd_line.format_of_input_seq,
								    cmd_line.max_read_length, 0, db_graph);

      //update the graph info object
      if (cmd_line.format_of_input_seq==FASTQ)
	{
	  graph_info_update_mean_readlen_and_total_seq(&db_graph_info, 0, calculate_mean(readlen_distrib, (long long) (cmd_line.max_read_length+1)), bases_pass_filters_and_loaded);
	}
      else//for FASTA we do not get read length distribution
	{
	  graph_info_increment_seq(&db_graph_info, 0, bases_pass_filters_and_loaded);
	}



      //cleanup marks left on nodes by loading process (for PE reads)
      hash_table_traverse(&db_node_set_status_to_none, db_graph);


      timestamp();
      if (cmd_line.format_of_input_seq==FASTQ)
	{
	  printf("Fastq data loaded\nTotal bases parsed:%qd\nTotal bases passing filters and loaded into graph:%qd\nMean read length after filters applied:%d\n", 
		 bases_parsed, bases_pass_filters_and_loaded, db_graph_info.mean_read_length[0]);
	}
      else
	{
	  printf("Fasta data loaded\nTotal bases parsed:%qd\nTotal bases passing filters and loaded into graph:%qd\n", 
		 bases_parsed, bases_pass_filters_and_loaded);
	}

      if (cmd_line.dump_readlen_distrib==true)
	{
	  FILE* rd_distrib_fptr = fopen(cmd_line.readlen_distrib_outfile, "w");
	  if (rd_distrib_fptr==NULL)
	    {
	      printf("Cannot open %s, so will dumping distribution of filtered read-lengths to stdout\n", cmd_line.readlen_distrib_outfile);
	      for (i=db_graph->kmer_size; i<=cmd_line.max_read_length; i++)
		{
		  printf("%d\t%qd\n", i, readlen_distrib[i]);
		}
	    }
	  else
	    {
	      printf("Dumping distribution of effective read lengths (ie after quality, homopolymer and/or PCR duplicate filters to file %s.\n", cmd_line.readlen_distrib_outfile);
	      for (i=db_graph->kmer_size; i<=cmd_line.max_read_length; i++)
		{
		  fprintf(rd_distrib_fptr, "%d\t%qd\n", i, readlen_distrib[i]);
		}
	      fclose(rd_distrib_fptr);
	    }
	  
	}


      free(readlen_distrib_ptrs);
      free(readlen_distrib);
      
    }
  else
    {
      //if there is a multicolour binary, load that in first
      timestamp();      
      int first_colour_data_starts_going_into=0;
      boolean graph_has_had_no_other_binaries_loaded=true;

      if (cmd_line.input_multicol_bin==true)
	{
	  int mean_readlens[NUMBER_OF_COLOURS];
	  int* mean_readlens_ptrs[NUMBER_OF_COLOURS];
	  long long total_seq_in_that_colour[NUMBER_OF_COLOURS];
	  long long* total_seq_in_that_colour_ptrs[NUMBER_OF_COLOURS];
	  int j;
	  for (j=0; j<NUMBER_OF_COLOURS; j++)
	    {
	      mean_readlens[j]=0;
	      mean_readlens_ptrs[j]=&(mean_readlens[j]);
	      total_seq_in_that_colour[j]=0;
	      total_seq_in_that_colour_ptrs[j]=&(total_seq_in_that_colour[j]);
	    }
	  long long  bp_loaded = load_multicolour_binary_from_filename_into_graph(cmd_line.multicolour_bin,db_graph, &first_colour_data_starts_going_into,
										  mean_readlens_ptrs, total_seq_in_that_colour_ptrs);
	  //update graph_info object
	  for (j=0; j<first_colour_data_starts_going_into; j++)
            {
	      graph_info_update_mean_readlen_and_total_seq(&db_graph_info, j, mean_readlens[j], total_seq_in_that_colour[j]);
	    }
	  timestamp();
	  printf("Loaded the multicolour binary %s, and got %qd kmers\n", cmd_line.multicolour_bin, bp_loaded/db_graph->kmer_size);
	  graph_has_had_no_other_binaries_loaded=false;
	}

      if (cmd_line.input_colours==true)
	{
	  timestamp();

	  //normal use
	  if (cmd_line.successively_dump_cleaned_colours==false)
	    {
	      
	      printf("List of colours: %s (contains one filelist per colour). Load data into consecutive colours starting at %d\n", 
		     cmd_line.colour_list, first_colour_data_starts_going_into);
	      if (cmd_line.load_colours_only_where_overlap_clean_colour==true)
		{
		  printf("When loading the binaries specified in %s, we only load nodes that are already in colour %d\n", 
			 cmd_line.colour_list, cmd_line.clean_colour);
		}
	      load_population_as_binaries_from_graph(cmd_line.colour_list, first_colour_data_starts_going_into, 
						     graph_has_had_no_other_binaries_loaded, db_graph, &db_graph_info,
						     cmd_line.load_colours_only_where_overlap_clean_colour, cmd_line.clean_colour);
	      timestamp();
	      printf("Finished loading single_colour binaries\n");
	    }
	  else
	    {
	      //we have loaded a multicolour binary, and we have checked that the clean_colour is one of the colours in that binary
	      if (cmd_line.load_colours_only_where_overlap_clean_colour==false)
		{
		  printf("If you specify --successively_dump_cleaned_colours, you must also specify --load_colours_only_where_overlap_clean_colour\n");
		  printf("That should fix your problem, however, this should have been caught as soon as Cortex parsed your command-line. Please inform Zam Iqbal (zam@well.ox.ac.uk) so he can fix that bug\n");
		  exit(1);
		}
	      printf("For each colour in %s, load data into graph, cleaning by comparison with colour %d, then dump a single-colour binary\n",
		     cmd_line.colour_list,cmd_line.clean_colour);
	      dump_successive_cleaned_binaries(cmd_line.colour_list, first_colour_data_starts_going_into,cmd_line.clean_colour,
					       cmd_line.successively_dump_cleaned_colours_suffix, db_graph, &db_graph_info);
	      printf("Completed dumping of clean binaries\n");
	    }



	}
    }
      
  printf("Total kmers in table: %qd\n", hash_table_get_unique_kmers(db_graph));	  
  printf("****************************************\n");
  printf("SUMMARY:\nColour:\tMeanReadLen\tTotalSeq\n");
  int j;
  for (j=0; j<NUMBER_OF_COLOURS; j++)
    {
      printf("%d\t%d\t%qd\n", j, db_graph_info.mean_read_length[j], db_graph_info.total_sequence[j]);
    }
  printf("****************************************\n");

  if (cmd_line.health_check==true)
    {
      timestamp();
      printf("Run health check on loaded graph\n");
      db_graph_health_check(false, db_graph);
      printf("End of health check\n");
      timestamp();
    }


  
  // Error Correction actions
  if (cmd_line.remove_seq_errors==true)
    {
      timestamp();
      printf("Remove nodes that look like sequencing errors. Clip tips first\n");
      db_graph_clip_tips_in_union_of_all_colours(db_graph);
      
      printf("Then remove low coverage nodes (<= 1) based on topology as well as covg  -  must be more likely to be single base error than due to sampling\n");
      db_graph_remove_errors_considering_covg_and_topology(1,db_graph, &element_get_covg_union_of_all_covgs, &element_get_colour_union_of_all_colours,
							   &apply_reset_to_specific_edge_in_union_of_all_colours, &apply_reset_to_all_edges_in_union_of_all_colours,
							   cmd_line.max_var_len);
      timestamp();
      printf("Error correction done\n");

    }
  else if (cmd_line.remove_low_coverage_nodes==true)
    {
      timestamp();
      printf("Start to to remove nodes with covg (in union of all colours)  <= %d\n", cmd_line.node_coverage_threshold);
      db_graph_remove_low_coverage_nodes_ignoring_colours(cmd_line.node_coverage_threshold, db_graph);
      timestamp();
      printf("Error correction done\n");
      
    }

  if (cmd_line.health_check==true)
    {
      timestamp();
      printf("Run health check on cleaned graph\n");
      db_graph_health_check(false, db_graph);
      printf("Health check done\n");
      timestamp();
    }

  if (cmd_line.dump_binary==true)
    {
      if (cmd_line.input_seq==true)
	{
	  //dump single colour
	  timestamp();
	  printf("Input data was fasta/q, so dump single colour binary file: %s\n", cmd_line.output_binary_filename);
	  db_graph_dump_single_colour_binary_of_colour0(cmd_line.output_binary_filename, &db_node_check_status_not_pruned,db_graph, &db_graph_info);
	  timestamp();
	  printf("Binary dumped\n");


	}
      else
	{
	  timestamp();
	  printf("Dump multicolour binary with %d colours (compile-time setting)\n", NUMBER_OF_COLOURS);
	  db_graph_dump_binary(cmd_line.output_binary_filename, &db_node_check_status_not_pruned,db_graph, &db_graph_info);
	  timestamp();
	  printf("Binary dumped\n");
	}
    }

  if (cmd_line.print_supernode_fasta==true)
    {
      timestamp();
      printf("Print contigs(supernodes) in the graph created by the union of all colours.\n");
      
      db_graph_print_supernodes_defined_by_func_of_colours(cmd_line.output_supernodes, "", cmd_line.max_var_len,// max_var_len is the public face of maximum expected supernode size
							   db_graph, &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs, 
							   &print_appropriate_extra_supernode_info);
      timestamp();
      printf("Supernodes dumped\n");


    }


  if (cmd_line.dump_covg_distrib==true)
    {
      timestamp();
      printf("Dump kmer coverage distribution for colour 0 to file %s\n", cmd_line.covg_distrib_outfile);
      db_graph_get_covg_distribution(cmd_line.covg_distrib_outfile, db_graph, individual_edge_array, 0, &db_node_check_status_not_pruned);
      timestamp();
      printf("Covg distribution dumped\n");
    }

  // DETECT BUBBLES

  if (cmd_line.detect_bubbles1==true)
    {
      timestamp();
      printf("Start first set of bubble calls\n");
      run_bubble_calls(&cmd_line, 1, db_graph, &print_appropriate_extra_variant_info, &get_colour_ref, &get_covg_ref);

      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Detect Bubbles 1, completed\n");
    }

  //second detect bubbles
  if (cmd_line.detect_bubbles2==true)
    {
      timestamp();
      printf("Start second set of bubble calls\n");
      run_bubble_calls(&cmd_line, 2, db_graph, &print_appropriate_extra_variant_info, &get_colour_ref, &get_covg_ref);
      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Detect Bubbles 2, completed\n");
    }

  if (cmd_line.make_pd_calls==true)
    {
      timestamp();
      printf("Run Path-Divergence Calls\n");
      run_pd_calls(&cmd_line, db_graph, &print_appropriate_extra_variant_info);
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Finished Path Divergence calls\n");
    }
  if (cmd_line.align_given_list==true)
    {
      timestamp();
      printf("Start aligning the fasta/q listed in this file: %s\n", cmd_line.list_fastaq_to_align);
      int array_of_colours[NUMBER_OF_COLOURS];
      int j;
      char* array_of_colournames[NUMBER_OF_COLOURS];
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  array_of_colours[j]=j;
	  array_of_colournames[j]=(char*)malloc(sizeof(char) * 50);
	  if (array_of_colournames[j]==NULL)
	    {
	      printf("Severe lack of memory. Cannot even allocate 50 chars. Give up\n");
	      exit(1);
	    }
	  sprintf(array_of_colournames[j], "colour_%d", j);
	}
      align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours(cmd_line.format_of_files_to_align, cmd_line.list_fastaq_to_align,
								       cmd_line.max_read_length, array_of_colours, array_of_colournames,
								       NUMBER_OF_COLOURS,db_graph,cmd_line.quality_score_offset,
								       false, NULL, NULL);
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  free(array_of_colournames[j]);
	}
      printf("Completed alignment of fasta/q to graph to print coverages in all colours\n");
    }
  if (cmd_line.print_colour_overlap_matrix==true)
    {
      timestamp();
      printf("Do a direct graph comparison of each colour x,y,z.. listed in --colour_subgraph_overlap_matrix a,b,c../x,y,z... against the colours a,b,c.. (process is symmetrical)\n");
      
      db_graph_print_colour_overlap_matrix(cmd_line.colour_overlap_first_colour_list,
					   cmd_line.num_colours_in_colour_overlap_first_list,
					   cmd_line.colour_overlap_second_colour_list,
					   cmd_line.num_colours_in_colour_overlap_second_list,
					   db_graph);
      printf("\nCompleted graph overlap matrix printing\n");
    }

  
  hash_table_free(&db_graph);
  timestamp();
  printf("Cortex completed\n");
  return 0;
}



void timestamp(){
 time_t ltime;
 ltime = time(NULL);
 printf("\n-----\n%s",asctime(localtime(&ltime)));
 fflush(stdout);
}
