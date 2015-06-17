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
/*
  cortex_var.c - cortex_var executable source file (contains 'main()' function)
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <string_buffer.h>

// cortex_var headers
#include "element.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "cmd_line.h"
#include "graph_info.h"
#include "db_differentiation.h"
#include "db_complex_genotyping.h"
#include "model_selection.h"
#include "experiment.h"
#include "genome_complexity.h"
#include "maths.h"
#include "seq_error_rate_estimation.h"

void timestamp();





int main(int argc, char **argv)
{
  CmdLine* cmd_line = cmd_line_alloc();
  if (cmd_line==NULL)
    {
      die("Out of memory!! Cannot even malloc the space to store your command-line arguments. Check who else is using your server, you seem to have severe problems\n");
    }
  
  parse_cmdline(cmd_line, argc,argv,sizeof(Element));

  int hash_key_bits, bucket_size;
  dBGraph * db_graph = NULL;
  short kmer_size;

  //next two needed much later
  //int num_kmers_dumped_after_alignment=0;
  char tmp_dump[300];




  CovgArray* working_ca_for_median=NULL;
  int lim = 2* cmd_line->max_var_len;//this 2* is deliberate, for the PD.
  if (cmd_line->max_read_length>lim)
    {
      lim = cmd_line->max_read_length;
    }
  if (cmd_line->print_median_covg_only==true)
    {
      working_ca_for_median = alloc_and_init_covg_array(lim);//will die if fails to alloc
    }

  void print_appropriate_extra_variant_info(AnnotatedPutativeVariant* annovar, FILE* fp)
  {
    if (cmd_line->print_colour_coverages==true)
      {
	print_standard_extra_info(annovar, fp);
      }
    else if (cmd_line->print_median_covg_only==true)
      {
	if (cmd_line->which_caller_was_used_for_calls_to_be_genotyped==BubbleCaller)
	  {
	    print_median_covg_extra_info(annovar, working_ca_for_median, fp);
	  }
	else
	  {
	    print_median_covg_on_informative_kmers_extra_info(annovar, working_ca_for_median, fp);
	  }
      }
    else
      {
	print_no_extra_info(annovar, fp);
      }
  }

  void print_appropriate_extra_supernode_info(dBNode** node_array, Orientation* or_array, int len, FILE* fout)
  {
    if (cmd_line->print_colour_coverages==true)
      {
	print_standard_extra_supernode_info(node_array, or_array, len, fout);
      }
    else if (cmd_line->print_median_covg_only==true)
      {
	print_median_extra_supernode_info(node_array, or_array, 
					  len, working_ca_for_median, fout);
      }
    else
      {
	print_no_extra_supernode_info(node_array, or_array, len, fout);
      }
  }


  Edges get_colour_ref(const dBNode* e)
  {

    if (cmd_line->using_ref==false)
      {
	die("Do not call get_colour_ref when --using_ref was not specified. Exiting - coding error\n");
      }

    if ( (cmd_line->ref_colour<0) || (cmd_line->ref_colour>NUMBER_OF_COLOURS-1) )
      {
	die("Calling get_colour_ref, but the reference colour %d has not been specified, "
      "or is > than the compile-time limit, of %d\n", 
	       cmd_line->ref_colour, NUMBER_OF_COLOURS-1);
      }
    Edges ed = get_edge_copy(*e, cmd_line->ref_colour);
    return ed;
  }

  Covg get_covg_ref(const dBNode* e)
  {
    if (cmd_line->using_ref==false)
      {
	die("Do not call get_coverage_ref when --using_ref was not specified. Exiting - coding error\n");
      }
    
    return e->coverage[cmd_line->ref_colour];
  }


  Covg get_covg_in_union_all_colours_except_ref(const dBNode* e)
  {
    Covg cov = element_get_covg_union_of_all_covgs(e);
    if (cmd_line->ref_colour!=-1)
    {
      cov -= e->coverage[cmd_line->ref_colour];
    }
    return cov;
  }

  //***************************************************************************
  // end local functions
  //***************************************************************************



  
  //set hash table variables:
  kmer_size        = cmd_line->kmer_size;
  hash_key_bits    = cmd_line->number_of_buckets_bits; //number of buckets: 2^hash_key_bits
  bucket_size      = cmd_line->bucket_size;

  int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
  int max_kmer_size = (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1;
  int min_kmer_size = ((NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1) * sizeof(bitfield_of_64bits) * 4) + 1;

  if (number_of_bitfields != NUMBER_OF_BITFIELDS_IN_BINARY_KMER)
  {
    die("K-mer %i  is not in current range of kmers [%i - %i] required for this executable!\n",
        kmer_size, min_kmer_size, max_kmer_size);
  }
  


  if(cmd_line->kmer_size > max_kmer_size)
  {
    die("k-mer size is too big [%i]!",cmd_line->kmer_size);
  }


  //Create the de Bruijn graph/hash table
  int max_retries=15;
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  if (db_graph==NULL)
    {
      die("Giving up - unable to allocate memory for the hash table\n");
    }



  GraphInfo* db_graph_info=graph_info_alloc_and_init();//will exit it fails to alloc.

  // input data:
  if (cmd_line->input_seq==true)
    {
      if (cmd_line->subsample==true)
	{
	  //	  printf("User has specified subsampling %f of the reads, so as we parse the fasta/q or BAM  we will\n load each read with probability %f.\n",
	  //	 cmd_line->subsample_propn, cmd_line->subsample_propn);
	}

      if (strcmp(cmd_line->se_list, "")!=0)
	{
	  //fprintf(stdout,"Input file of single ended data filenames: %s\n",
	  //	  cmd_line->se_list);
	  //cmd_line->loaded_sample_names = get_sample_id_from_se_pe_list(cmd_line->colour_sample_ids[0], cmd_line->se_list);
	}
      else
	{
	  printf("No SE data\n");
	}
      if (strcmp(cmd_line->pe_list_lh_mates, "") !=0)
	{
	  fprintf(stdout,"Input file of paired end data: %s, and %s, \n",
		  cmd_line->pe_list_lh_mates, cmd_line->pe_list_rh_mates); 
	}
      else
	{
	  //printf("No paired-end data\n");
	}
      if (cmd_line->quality_score_threshold>0)
	{
	  fprintf(stdout,"quality cut-off: %i\n",cmd_line->quality_score_threshold);
	}
      else
	{
	  fprintf(stdout, "No quality filtering\n");
	}
      if (cmd_line->remove_pcr_dups==true)
	{
	  //printf("Removing duplicates from the paired end files when both mates "
	  //           "start with the same kmer\n");
	}
      else
	{
	  //printf("No PCR duplicate removal\n");
	}
      if (cmd_line->cut_homopolymers==true)
	{
	  printf("Breaking reads at homopolymer runs of length %d or greater. "
           "Read restarts at first base after the run\n",
           cmd_line->homopolymer_limit);
	}

    // note, if we load fasta/fastq data it always goes into colour 0
    int into_colour = 0;

    unsigned int num_files_loaded = 0;
    unsigned long long num_bad_reads = 0;
    unsigned long long num_dup_reads = 0;
    unsigned long long num_bases_parsed = 0;
    unsigned long long num_bases_loaded = 0;

    // Get read-len distribution (after filters):
    int max_expected_read_len = cmd_line->max_read_length;
    if (max_expected_read_len==0)
      {
	max_expected_read_len=20000;
      }
    unsigned long readlen_distrib_size = max_expected_read_len + 1;
    unsigned long *readlen_distrib
      = (unsigned long*) malloc(sizeof(unsigned long) * readlen_distrib_size);
    
    if(readlen_distrib == NULL)
    {
      die("Unable to malloc array to hold readlen distirbution!Exit.\n");
    }

    unsigned long i;
    for(i = 0; i < readlen_distrib_size; i++)
    {
      readlen_distrib[i] = 0;
    }

    timestamp();

    int homopolymer_cutoff
      = cmd_line->cut_homopolymers ? cmd_line->homopolymer_limit : 0;

    boolean (*subsample_function)();

    //local func
    boolean subsample_as_specified()
    {
      double ran = drand48();
      if (ran <= cmd_line->subsample_propn)
	{
	  return true;
	}
      return false;
    }
    //end of local func

    if (cmd_line->subsample==true)
      {
	subsample_function = &subsample_as_specified;
      }
    else
      {
	subsample_function = &subsample_null;
      }

    if(strcmp(cmd_line->se_list, "") != 0)
      {
	load_se_filelist_into_graph_colour(cmd_line->se_list,
					   cmd_line->quality_score_threshold, homopolymer_cutoff, false,
					   cmd_line->quality_score_offset,
					   into_colour, db_graph, 0, // 0 => filelist not colourlist
					   &num_files_loaded, &num_bad_reads, &num_dup_reads,
					   &num_bases_parsed, &num_bases_loaded,
					   readlen_distrib, readlen_distrib_size,
					   subsample_function);
      }

    if(strcmp(cmd_line->pe_list_lh_mates, "") != 0)
      {
	load_pe_filelists_into_graph_colour(
					    cmd_line->pe_list_lh_mates, cmd_line->pe_list_rh_mates,
					    cmd_line->quality_score_threshold, homopolymer_cutoff, 
					    cmd_line->remove_pcr_dups,
					    cmd_line->quality_score_offset,
					    into_colour, db_graph, 0, // 0 => filelist not colourlist
					    &num_files_loaded, &num_bad_reads, &num_dup_reads,
					    &num_bases_parsed, &num_bases_loaded,
					    readlen_distrib, readlen_distrib_size,
					    subsample_function);
    }

    // Update the graph info object
    unsigned long mean_contig_length = calculate_mean_ulong(readlen_distrib,
                                                            readlen_distrib_size);
    
    graph_info_update_mean_readlen_and_total_seq(db_graph_info, 0,
                                                 mean_contig_length,
                                                 num_bases_loaded);
    
    if(cmd_line->entered_sampleid_as_cmdline_arg == true)
      {
	graph_info_set_sample_ids(cmd_line->colour_sample_ids, 1, db_graph_info, 0);
      }
    
    // Cleanup marks left on nodes by loading process (for PE reads)
    hash_table_traverse(&db_node_set_status_to_none, db_graph);
    
    hash_table_print_stats(db_graph);
    
    timestamp();
    
    //printf("Sequence data loaded\n");
    //printf("Total bases parsed:%llu\n", num_bases_parsed);
    //printf("Total bases passing filters and loaded into graph:%llu\n",
    //           num_bases_loaded);
    //printf("Mean read length after filters applied:%lu\n", mean_contig_length);
    
    if(cmd_line->dump_readlen_distrib == true)
      {
	FILE* rd_distrib_fptr = fopen(cmd_line->readlen_distrib_outfile, "w");
	
	if(rd_distrib_fptr == NULL)
	  {
	    printf("Cannot open %s, so dump distribution of filtered "
		   "read-lengths to stdout (ie to the screen)\n", cmd_line->readlen_distrib_outfile);
	    
	    rd_distrib_fptr = fdopen(fileno(stdout), "w");
	  }
	else
	  {
	    printf("Dumping distribution of effective read lengths (ie after "
		   "quality, homopolymer and/or PCR duplicate filters) to file %s.\n",
		   cmd_line->readlen_distrib_outfile);
	  }
	
	for(i = db_graph->kmer_size; i < readlen_distrib_size; i++)
	  {
	    fprintf(rd_distrib_fptr, "%lu\t%lu\n", i, readlen_distrib[i]);
	  }
	
	fclose(rd_distrib_fptr);
      }
    
    free(readlen_distrib);
    }
  else
    {
      //if there is a multicolour binary, load that in first
      timestamp();      
      int first_colour_data_starts_going_into=0;
      boolean graph_has_had_no_other_binaries_loaded=true;
      
      if (cmd_line->input_multicol_bin==true)
	{
	  long long  bp_loaded = load_multicolour_binary_from_filename_into_graph(cmd_line->multicolour_bin,db_graph, 
										  db_graph_info, &first_colour_data_starts_going_into);
	  timestamp();
	  printf("Loaded the multicolour binary %s, and got %qd kmers\n", cmd_line->multicolour_bin, bp_loaded/db_graph->kmer_size);
	  graph_has_had_no_other_binaries_loaded=false;
	  timestamp();
	  
	}
      
      if (cmd_line->input_colours==true)
	{
	  timestamp();
	  
	  //normal use
	  
	  if (cmd_line->successively_dump_cleaned_colours==false)
	    {
	      
	      printf("List of colours: %s (contains one filelist per colour). Load data into consecutive colours starting at %d\n", 
		     cmd_line->colour_list, first_colour_data_starts_going_into);
	      if (cmd_line->load_colours_only_where_overlap_clean_colour==true)
		{
		  printf("When loading the binaries specified in %s, we only load nodes that are already in colour %d\n", 
			 cmd_line->colour_list, cmd_line->clean_colour);
		}
	      
	      load_population_as_binaries_from_graph(cmd_line->colour_list, first_colour_data_starts_going_into, 
						     graph_has_had_no_other_binaries_loaded, db_graph, db_graph_info,
						     cmd_line->load_colours_only_where_overlap_clean_colour, cmd_line->clean_colour,
						     cmd_line->for_each_colour_load_union_of_binaries);
	      
	      //if the colour_list contained sample_ids, add them to the GraphInfo object
	      // these will override the sample-id in the binary
	      if (cmd_line->loaded_sample_names==true)
		{
		  graph_info_set_sample_ids(cmd_line->colour_sample_ids, cmd_line->num_colours_in_input_colour_list,
					    db_graph_info, first_colour_data_starts_going_into);
		}
	      timestamp();
	      printf("Finished loading single_colour binaries\n");
	    }
	  else//we are going to clean a list of binaries against one of the colours in the multicolir bin
	    {
	      //we have loaded a multicolour binary, and we have checked that the clean_colour is one of the colours in that binary
	      if (cmd_line->load_colours_only_where_overlap_clean_colour==false)
		{
		  die(
		      "If you specify --successively_dump_cleaned_colours, you must also specify\n"
		      "--load_colours_only_where_overlap_clean_colour That should fix your problem,\n"
		      "however, this should have been caught as soon as Cortex parsed your command-line.\n"
		      "Please inform Zam Iqbal (zam@well.ox.ac.uk) so he can fix that UI bug\n");
		}
	      
	      
	      printf("For each colour in %s, load data into graph, cleaning by comparison with colour %d, then dump a single-colour binary\n",
		     cmd_line->colour_list,cmd_line->clean_colour);
	      graph_info_set_specific_colour_to_cleaned_against_pool(db_graph_info,  first_colour_data_starts_going_into, 
								     cmd_line->multicolour_bin, cmd_line->clean_colour);
	      
	      db_graph_info->cleaning[first_colour_data_starts_going_into]->cleaned_against_another_graph=true;
	      dump_successive_cleaned_binaries(cmd_line->colour_list, first_colour_data_starts_going_into,cmd_line->clean_colour,
					       cmd_line->successively_dump_cleaned_colours_suffix, db_graph, db_graph_info);
	      
	      graph_info_unset_specific_colour_from_cleaned_against_pool(db_graph_info, first_colour_data_starts_going_into);
	      db_graph_info->cleaning[first_colour_data_starts_going_into]->cleaned_against_another_graph=false;
	      printf("Completed dumping of clean binaries\n");
	      
	      timestamp();
	      
	    }
	  
	  
	  
	}
    }
  
  
  
  
  GraphAndModelInfo model_info;
  float repeat_geometric_param_mu = 0.8;
  // float seq_err_rate_per_base;
  int i;
  // need estimated sequencing error rates for genotyping
  if ((cmd_line->manually_override_error_rate==false) && (cmd_line->use_snp_alleles_to_estim_seq_err_rate==false))
    {
      //if you manually override, you set the same seq error rate for all colours 
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  if (i != cmd_line->ref_colour)
	    {
	      db_graph_info->seq_err[i]=0.01;
	    }
	}
    }
  else if ((cmd_line->manually_override_error_rate==true) && (cmd_line->use_snp_alleles_to_estim_seq_err_rate==false))
    {
      if (strcmp(cmd_line->manually_entered_seq_error_rates_file, "")==0)
	{
	  for (i=0; i<NUMBER_OF_COLOURS; i++)
	    {
	      if (i != cmd_line->ref_colour)
		{
		  db_graph_info->seq_err[i]=cmd_line->manually_entered_seq_error_rate;
		}
	    }
	}
      else
	{
	  FILE* fp_err = fopen(cmd_line->manually_entered_seq_error_rates_file, "r");
	  if (fp_err==NULL)
	    {
	      die("Unable to open file %s\n", cmd_line->manually_entered_seq_error_rates_file);
	    }
	  read_estimated_seq_errors_from_file(db_graph_info, fp_err);
	  fclose(fp_err);
	}
    }
  else if (cmd_line->use_snp_alleles_to_estim_seq_err_rate==true)
    {
      long double default_seq_err = 0.01;
      estimate_seq_error_rate_from_snps_for_each_colour(cmd_line->colourlist_snp_alleles, 
							db_graph_info, db_graph, 
							cmd_line->ref_colour, 
							//cmd_line->genome_size, //unused
							default_seq_err,
							"seq_err_estimation.txt");
    }
  else
    {
      die("Should never reach here %s:%i", __FILE__, __LINE__);
    }

  
  print_seq_err_rates_to_screen(db_graph_info);


  int num_chroms_in_expt=NUMBER_OF_COLOURS;
  if (cmd_line->expt_type==EachColourADiploidSample)
    {
      num_chroms_in_expt=2*NUMBER_OF_COLOURS;
    }
  else if (cmd_line->expt_type==EachColourADiploidSampleExceptTheRefColour)
    {
      num_chroms_in_expt=2*NUMBER_OF_COLOURS-2;
    }
  else if (cmd_line->expt_type==EachColourAHaploidSample)
    {
      num_chroms_in_expt=NUMBER_OF_COLOURS;
    }
  else if (cmd_line->expt_type==EachColourAHaploidSampleExceptTheRefColour)
    {
      num_chroms_in_expt=NUMBER_OF_COLOURS-1;
    }
  else if (cmd_line->expt_type==Unspecified)
    {
      num_chroms_in_expt=NUMBER_OF_COLOURS;
    }
  else
    {
      die("Coding error - expt type not specified - not even as Unspecified");
    }

  AssumptionsOnGraphCleaning assump = AssumeUncleaned;
  //  AssumptionsOnGraphCleaning assump = AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice;
  initialise_model_info(&model_info, db_graph_info, cmd_line->genome_size, 
			repeat_geometric_param_mu, //seq_err_rate_per_base, 
			cmd_line->ref_colour, num_chroms_in_expt, 
			cmd_line->expt_type, assump);


  int j;







  if (cmd_line->phelim==true)
    {

      timestamp();
      printf("This is a demo for Phelim. I'm going to assume we have loaded a graph with k=31, and print the supernode containing the kmer TTGGTGAATAGTAATTTAATTCTATCTCATA, which is in the middle of blaZ\n");
      
      BinaryKmer tmp_kmer1;
      BinaryKmer tmp_kmer2;

      //first get the node which corresponds to our kmer string
      dBNode* query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TTGGTGAATAGTAATTTAATTCTATCTCATA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);


      //some setting up of variables and memory

      //a supernode consists of 4 arrays: the nodes, the orientations (since a node encodes both forwards and rev comp), the edges (labels) and a string which gives the sequence of the edges
      dBNode * *    path_nodes;
      Orientation * path_orientations;
      Nucleotide *  path_labels;
      char * seq;
      boolean is_cycle;//we will call a function below to get the supernode containing query_node, and it will tell us if it is a loop/cycle
      double avg_coverage;//that function will also give us average/min/max coverage on the supernode. legacy.
      Covg min_covg, max_covg;
  
      int max_length=10000;//whatever
      path_nodes        = calloc(max_length,sizeof(dBNode*));
      path_orientations = calloc(max_length,sizeof(Orientation));
      path_labels       = calloc(max_length,sizeof(Nucleotide));
      seq               = calloc(max_length+1+db_graph->kmer_size,sizeof(char));
      //end of setup
      

      //sorry about the function name. "for_specific_person_or_pop" really means for some colour
      int colour=0;
      int supernode_length = db_graph_supernode_for_specific_person_or_pop(query_node,max_length,
									   &db_node_action_set_status_visited,//you can decide what to do to nodes you have seen. Here mark as visited
									   path_nodes,path_orientations,path_labels,seq,
									   &avg_coverage, &min_covg, &max_covg, &is_cycle,
									   db_graph, colour);
      

      //fyi there are other functions that allow you to get a supernode in the graph determined not by a single colour, but any combination of colours.
      //there are also functions for getting supernodes which satisfy conditions (where you pass in the condition function)

      //OK< we now have the supernode in path_nodes, path_orientations, path_labels and seq
      printf("Supernode is %s, of length %d\n", seq, supernode_length);


      //now what happens if we try the reverse complement of that kmer?
      query_node = hash_table_find(element_get_key(seq_to_binary_kmer("TATGAGATAGAATTAAATTACTATTCACCAA", kmer_size, &tmp_kmer1), kmer_size, &tmp_kmer2) ,db_graph);
      supernode_length = db_graph_supernode_for_specific_person_or_pop(query_node,max_length,
								       &db_node_action_set_status_visited,//you can decide what to do to nodes you have seen. Here mark as visited
								       path_nodes,path_orientations,path_labels,seq,
								       &avg_coverage, &min_covg, &max_covg, &is_cycle,
								       db_graph, colour);

      printf("\n\nIf we pass in the reverse complement node: TATGAGATAGAATTAAATTACTATTCACCAA, we get this Supernode : %s, of length %d\n", seq, supernode_length);

      
      //cleanup
      free(path_nodes);
      free(path_orientations);
      free(path_labels);
      free(seq);
    }


  if (cmd_line->dump_covg_distrib==true)
    {
      timestamp();
      printf("Dump kmer coverage distribution for colour 0 to file %s\n", 
	     cmd_line->covg_distrib_outfile);
      db_graph_get_covg_distribution(cmd_line->covg_distrib_outfile, 
				     db_graph, 
				     0, 
				     &db_node_check_status_not_pruned);
      timestamp();
      printf("Covg distribution dumped\n");
    }


  
  hash_table_free(&db_graph);



  cmd_line_free(cmd_line);
  graph_info_free(db_graph_info);
  printf("Cortex completed - y'all have a nice day!\n");
  return 0;
}



void timestamp(){
 time_t ltime;
 ltime = time(NULL);
 printf("\n-----\n%s",asctime(localtime(&ltime)));
 fflush(stdout);
}
