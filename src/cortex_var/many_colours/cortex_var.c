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



void run_genotyping(CmdLine* cmd_line, dBGraph* db_graph,
                    void (*print_whatever_extra_variant_info)(VariantBranchesAndFlanks*, FILE*), 
                    GraphAndModelInfo* model_info)
                    //Edges(*get_col_ref)(const dBNode* e)
                    //Covg (*get_cov_ref)(const dBNode* e)
                    //GraphInfo* db_graph_info
{
  FILE* fp = fopen(cmd_line->file_of_calls_to_be_genotyped, "r");
  if (fp==NULL)
    {
      die("Cannot open file %s\n", cmd_line->file_of_calls_to_be_genotyped);
    }
  else
    {
      int gt_file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
	long long ret;
	int offset = 0;
	if (new_entry!= true){
	  die("new_entry has to be true for read_next_variant_from_full_flank_file\n");
	}	
	ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
	return ret;
      }


      //----------------------------------
      // allocate the memory used to read the sequences
      //----------------------------------
      Sequence * seq = malloc(sizeof(Sequence));
      if (seq == NULL){
	die("Out of memory trying to allocate Sequence");
      }
      alloc_sequence(seq,cmd_line->max_read_length,LINE_MAX);

      Sequence * seq_inc_prev_kmer = malloc(sizeof(Sequence));
      if (seq == NULL){
	die("Out of memory trying to allocate Sequence");
      }
      alloc_sequence(seq_inc_prev_kmer,cmd_line->max_read_length+db_graph->kmer_size,LINE_MAX);

      
      
      //We are going to load all the bases into a single sliding window 
      KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
      if (kmer_window==NULL)
	{
	  die("Failed to malloc kmer sliding window for genotyping. Exit.\n");
	}
      kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(cmd_line->max_read_length-db_graph->kmer_size-1));
      if (kmer_window->kmer==NULL)
	{
	  die("Failed to malloc kmer_window->kmer for genotyping. Tried to alloc\n"
        "an array of %d binary kmers. Max read len:%d, kmer size %d. Exit.\n", 
		    cmd_line->max_read_length-db_graph->kmer_size-1,
        cmd_line->max_read_length,
        db_graph->kmer_size);
	}
      kmer_window->nkmers=0;
      
      CovgArray* working_ca = alloc_and_init_covg_array(cmd_line->max_read_length+1);
      if (working_ca==NULL)
	{
	  die("Unable to alloc an array of coverages for genotyping - abort. Your server must have run out of memory.\n");
	}
      VariantBranchesAndFlanks* var = alloc_VariantBranchesAndFlanks_object(cmd_line->max_read_length+1, cmd_line->max_read_length+1, 
									    cmd_line->max_read_length+1, cmd_line->max_read_length+1, db_graph->kmer_size);
      if (var==NULL)
	{
	  die("Abort - unable to allocate memory for buffers for reading callfile - \n"
        "either sever oom conditions or you have specified very very large max_read_len\n");
	}
      GenotypingWorkingPackage* gwp=NULL;
      LittleHashTable* little_dbg=NULL;
      if (cmd_line->which_caller_was_used_for_calls_to_be_genotyped==SimplePathDivergenceCaller)
	{
	  int max = cmd_line->max_read_length;
	  if (cmd_line->max_var_len > max)
	    {
	      max = cmd_line->max_var_len;
	    }
	  gwp = alloc_genotyping_work_package(max,max,
					      NUMBER_OF_COLOURS, 
					      NUMBER_OF_COLOURS+1);
	  if (gwp==NULL)
	    {
	      die("Unable to alloc resources for genotyping prior to starting. Abort\n");
	    }
	  int little_width=100;
	  int little_height=(int) log((cmd_line->max_read_length) *4)+1;
	  int little_retries=20;
	  printf("Building little hash for genotyping: height %d and width %d\n", little_height, little_width);
	  little_dbg = little_hash_table_new(little_height, little_width, little_retries, db_graph->kmer_size);
	  if (little_dbg==NULL)
	    {
	      die("Out of memory - failed to alloc tiny auxiliary hash table\n");
	    }
	  var->which=first;
	}
      //end of initialisation 
      
      FILE* fout = fopen(cmd_line->output_genotyping, "w");
      if (fout==NULL)
	{
	  die("Unable to open output file %s - abort.\n", cmd_line->output_genotyping);
	}
      
      int ret=1;

      while (ret)
	{

	  //note you are reading a bunch of variants that may have been called on another sample,
	  //and potentially are genotyping them on a different sample. So these reads can contain kmers
	  // that are not in our graph. These become NULL points in our array of dBNode* 's.
	  ret = read_next_variant_from_full_flank_file(fp, cmd_line->max_read_length,
						       var, db_graph, &gt_file_reader, seq, seq_inc_prev_kmer, kmer_window);
	  if (ret==1)
	    {
	      AssumptionsOnGraphCleaning assump=AssumeUncleaned; 
	      print_call_given_var_and_modelinfo(var, fout, model_info, cmd_line->which_caller_was_used_for_calls_to_be_genotyped, db_graph,
						 print_whatever_extra_variant_info,
						 assump, gwp, little_dbg, working_ca);
	    }
	  else
	    {

	    }
	}

      //cleanup
      fclose(fout);
      fclose(fp);
      free_VariantBranchesAndFlanks_object(var);
      free_covg_array(working_ca);
      if (cmd_line->which_caller_was_used_for_calls_to_be_genotyped==SimplePathDivergenceCaller)
	{
	  free_genotyping_work_package(gwp); 
	  little_hash_table_free(&little_dbg);
	}
      free_sequence(&seq);
      free_sequence(&seq_inc_prev_kmer);
      free(kmer_window->kmer);
      free(kmer_window);
    }
}


void run_pd_calls(CmdLine* cmd_line, dBGraph* db_graph, 
		  void (*print_some_extra_var_info)(VariantBranchesAndFlanks* var, FILE* fp),
		  GraphAndModelInfo* model_info)
{
  if (cmd_line->pd_calls_against_each_listed_colour_consecutively==false)
    {
      printf("Calling variants using Path Divergence Caller.\n");
      printf("Calls are made between the reference path as specified by the fasta in %s\n", cmd_line->ref_chrom_fasta_list);
      printf(" and the union of the following colour(s): ");
      int i;
      for (i=0; i< cmd_line->num_colours_in_pd_colour_list-1; i++)
	{
	  printf("%d,", cmd_line->pd_colour_list[i]);
	}
      printf("%d\n", cmd_line->pd_colour_list[cmd_line->num_colours_in_pd_colour_list-1]);
    }
  else
    {
      printf("Calling variants using Path Divergence Caller. Since you used semicolons preceding and separating your list of colours,\n");
      printf("we treat each of those colours separately. i.e we do DISCOVERY on each colour in the list independently,\n");
      printf(" but for all variants called, however they were discovered, (if you specify --experiment_type and --genome_size, then) we will GENOTYPE all samples/colours.\n");
      printf("Calls are made between the reference path as specified by the fasta in %s\n", cmd_line->ref_chrom_fasta_list);
      printf(" and each of the following samples in turn. ");
      int i;
      for (i=0; i< cmd_line->num_colours_in_pd_colour_list-1; i++)
	{
	  printf("%d,", cmd_line->pd_colour_list[i]);
	}
      printf("%d\n", cmd_line->pd_colour_list[cmd_line->num_colours_in_pd_colour_list-1]);
      
    }
  printf("The reference colour is %d\n", cmd_line->ref_colour);
  if (model_info->expt_type==Unspecified)
    {
      printf("Since you did not use --experiment_type, Cortex does not know how if each colour is a diploid/haploid sample or pool, so\nwill not calculate genotypes or likelihoods\n");
    }
  else if (cmd_line->genome_size==0)
    {
      printf("Since you did not specify the genome size/length, Cortex cannot calculate genotype likelihoods or call genotypes\n");
    }

  // This will also check that all the ref chrom files exist
  int num_ref_chroms = load_paths_from_filelist(cmd_line->ref_chrom_fasta_list,
                                                NULL);

  char** ref_chroms = malloc(sizeof(char*) * num_ref_chroms);

  if(ref_chroms==NULL)
  {
    die("Can't allocate memory for paths to ref chromosome files\n");
  }

  load_paths_from_filelist(cmd_line->ref_chrom_fasta_list, ref_chroms);

  // Now set up output file names
  char* output_file = malloc(sizeof(char)*1000);
  sprintf(output_file, "%s_pd_calls", cmd_line->path_divergence_caller_output_stub);

  FILE* out_fptr = fopen(output_file, "w");
  if (out_fptr==NULL)
    {
      die("Cannot open %s for output\n", output_file);
    }



  int min_fiveprime_flank_anchor = 2;
  int min_threeprime_flank_anchor= cmd_line->kmer_size;
  int max_anchor_span =  cmd_line -> max_var_len;
  int length_of_arrays = 2*max_anchor_span;
  int min_covg =1;
  int max_covg = 10000000;//this is ignored. will be changing API
  int max_expected_size_of_supernode=cmd_line -> max_var_len;
	int i;
      
  if (cmd_line->pd_calls_against_each_listed_colour_consecutively==false)
    {
      int global_var_counter=0;
      for (i=0; i<num_ref_chroms; i++) 
	{
	  printf("Call SV comparing individual/sample  with chromosome %s\n", ref_chroms[i]);
	  
	  FILE* chrom_fptr = fopen(ref_chroms[i], "r");
	  if (chrom_fptr==NULL)
	    {
	      die("Cannot open %s \n", ref_chroms[i]);
	    }
	  
	  
	  
	  global_var_counter += db_graph_make_reference_path_based_sv_calls_given_list_of_colours_for_indiv(cmd_line->pd_colour_list, cmd_line->num_colours_in_pd_colour_list,
													    chrom_fptr, cmd_line->ref_colour,
													    min_fiveprime_flank_anchor, min_threeprime_flank_anchor, 
													    max_anchor_span, min_covg, max_covg, 
													    max_expected_size_of_supernode, length_of_arrays, db_graph, out_fptr,
													    0, NULL, NULL, NULL, NULL, NULL, 
													    &make_reference_path_based_sv_calls_condition_always_true_in_subgraph_defined_by_func_of_colours, 
													    &db_variant_action_do_nothing,
													    //print_some_extra_var_info, model_info,  AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice, global_var_counter+1);
													    print_some_extra_var_info, model_info,  AssumeUncleaned, global_var_counter+1);
	  
	  
	  fclose(chrom_fptr);
	}      
    }
  else
    {
      int global_var_counter=0;
      for (i=0; i<num_ref_chroms; i++) 
	{
	  printf("Call SVs by comparing each colour in turn  with chromosome %s\n", ref_chroms[i]);

	  int p;
	  for (p=0; p<cmd_line->num_colours_in_pd_colour_list; p++)
	    {
	      FILE* chrom_fptr = fopen(ref_chroms[i], "r");
	      if (chrom_fptr==NULL)
		{
		  die("Cannot open %s \n", ref_chroms[i]);
		}
	      
	      
	      int list_one_col[1];
	      list_one_col[0]=cmd_line->pd_colour_list[p];
	      
	      global_var_counter += db_graph_make_reference_path_based_sv_calls_given_list_of_colours_for_indiv(list_one_col, 1,
														chrom_fptr, cmd_line->ref_colour,
														min_fiveprime_flank_anchor, min_threeprime_flank_anchor, 
														max_anchor_span, min_covg, max_covg, 
														max_expected_size_of_supernode, length_of_arrays, db_graph, out_fptr,
														0, NULL, NULL, NULL, NULL, NULL, 
														&make_reference_path_based_sv_calls_condition_always_true_in_subgraph_defined_by_func_of_colours, //just always returns true
														&db_variant_action_do_nothing,
														//print_some_extra_var_info, model_info, AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice,
														print_some_extra_var_info, model_info, AssumeUncleaned,
														global_var_counter+1);
	  
	  
	      fclose(chrom_fptr);
	      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
	    }

	}
      
    }
  //cleanup
  fclose(out_fptr);

  for(i=0; i<num_ref_chroms; i++)
    {
      free(ref_chroms[i]);
      //free(output_files[i]);
    }
  free(ref_chroms);
  //free(output_files);

  

}

void run_bubble_calls(CmdLine* cmd_line, int which, dBGraph* db_graph, 
                      void (*print_appropriate_extra_var_info)(VariantBranchesAndFlanks* var, FILE* fp),
                      Edges(*get_col_ref) (const dBNode* e),
                      Covg (*get_cov_ref)(const dBNode* e),
                      //GraphInfo* db_graph_info,
                      GraphAndModelInfo* model_info)
{


  printf("Detecting bubbles between the union of this set of colours: ");
  int k;

  for (k=0; k<cmd_line->num_colours_in_detect_bubbles1_first_colour_list; k++)
    {
      printf("%d, ", cmd_line->detect_bubbles1_first_colour_list[k]);
    }

  printf("\nand the union of this set of colours: ");
    
  for (k=0; k<cmd_line->num_colours_in_detect_bubbles1_second_colour_list; k++)
    {
      printf("%d, ", cmd_line->detect_bubbles1_second_colour_list[k]);
    }
  

  if (cmd_line->exclude_ref_bubbles==true)
    {
      printf("Will remove bubbles in the reference colour %d before doing any calling\n", cmd_line->ref_colour);
    }
  
  long double* allele_balances=NULL;
  CovgArray* working_ca=NULL;
  if (cmd_line->high_diff==true)
    {
      printf("User specified calling only if highly-differentiated sites - will only call sites where allele balance differs by at least %f between at least 2 colours (excluding nay reference\n", cmd_line->min_allele_balance_diff);

      allele_balances=calloc(NUMBER_OF_COLOURS, sizeof(long double));
      working_ca = alloc_and_init_covg_array(cmd_line->max_var_len);//will die if fails to alloc
      if (allele_balances==NULL)
	{
	  die("Completely out of memory - can't even alloc a small array\n");
	}
    }

  if (cmd_line->apply_model_selection_at_bubbles==true)
    {
      printf("Will compare likelihoods of two models (Repeat and  Variation (Hardy-Weinberg)) at all bubbles,\nand mark those more likely to be repeats for filtering\n");
    }
  
  if (model_info->expt_type==Unspecified)
    {
      printf("Since you did not use --experiment_type, Cortex does not know how if each colour is a diploid/haploid sample or pool, so\nwill not calculate genotypes or likelihoods\n");
    }
  else if (cmd_line->genome_size==0)
    {
      printf("Since you did not specify the genome size/length, Cortex cannot calculate genotype likelihoods\n");
    }
  FILE* fp;
  
  fp = fopen(cmd_line->output_detect_bubbles1, "w");

  //local function
  boolean bubble_caller_condition(VariantBranchesAndFlanks* var)
  {
    if (cmd_line->high_diff==false)
      {
	return true;
      }
    else//we want only high differentiation sites
      {
	int i;
	int count=0;
	
	for (i=0; i<NUMBER_OF_COLOURS; i++)//for each colour
	  {
	    boolean too_short=false;
	    if ( (model_info->ginfo->total_sequence[i]>0) //ignore colours with no data
	      &&
		 (model_info->ginfo->mean_read_length[i]>0) )
	      {
		if (i==cmd_line->ref_colour)
		  {
		    float median_covg_allele1_in_ref = median_covg_on_allele_in_specific_colour(var->one_allele, var->len_one_allele, 
												working_ca, cmd_line->ref_colour, &too_short);
		    float median_covg_allele2_in_ref = median_covg_on_allele_in_specific_colour(var->other_allele, var->len_other_allele, 
												working_ca, cmd_line->ref_colour, &too_short);
		    if ( (median_covg_allele1_in_ref==0) && (median_covg_allele2_in_ref==0) )
		      {
			return false; //I insist on covg on ref allele as one of these two
		      }
		  }
		else//this is not the reference colour
		  {
		    uint64_t       R = (uint64_t) (model_info->ginfo->mean_read_length[i]);
		    //uint64_t factor1 = (R -db_graph->kmer_size+1);
		    //float factor2 = factor1/R;
		    //uint64_t seq     = (uint64_t) (model_info->ginfo->total_sequence[i]);
		    //uint64_t depth   = seq/( (uint64_t)(model_info->genome_len) );
		    //float exp_covg = factor2 * ((float)depth);
		    
		    Covg median_covg_allele1 = median_covg_on_allele_in_specific_colour(var->one_allele, var->len_one_allele,
											working_ca, i, &too_short);
		    Covg median_covg_allele2 = median_covg_on_allele_in_specific_colour(var->other_allele, var->len_other_allele,
											working_ca, i, &too_short);
		    
		    Covg sum =sum_covgs(median_covg_allele1,median_covg_allele2);
		    // exp depth = depth*factor1/R = seq*factor1/(R*genome_size)
		    //if ( (sum>0) && (2*R*sum*model_info->genome_len >=seq * factor1) )  //right hand expression means sum>= exp covg/2, but avoids division
		    if (sum>0)
		      {
			long double ab = (long double)median_covg_allele1/(long double)sum;
			long double rounded_ab = floorl(ab * 100 + 0.5) / 100;
			allele_balances[count]=rounded_ab;
			count++;//now, count is the number of stored allele_balances
		      }
		    else
		      {
			//one of the populations has no covg on BOTH alleles. Discard this site
			return false;
		      }
		    
		  }
	      }//end of "if there is some data in all colours

	  }//end of for each colour

	qsort(allele_balances, count, sizeof(long double), long_double_cmp);
	if (abs(allele_balances[0]-allele_balances[count-1])>= cmd_line->min_allele_balance_diff)
	  {
	    fprintf(fp,"top and bottom allele balances are %Lf and %Lf and count is %d\n", allele_balances[0], allele_balances[count-1], count);
	    return true;
	  }
	return false;

      }//end of calling high-diff sites
  }
  //end local function

  
  if (fp==NULL)
    {
      die("Cannot open %s. Exit.", cmd_line->output_detect_bubbles1);
    }
  
  boolean (*mod_sel_criterion)(AnnotatedPutativeVariant* annovar,  GraphAndModelInfo* model_info)=NULL;
  
  if (cmd_line->apply_model_selection_at_bubbles==true)
    {
      mod_sel_criterion = &basic_model_selection_condition; 
    }
  
  
  if (cmd_line->exclude_ref_bubbles==true)
    {
      printf("(First exclude bubbles from ref colour %d) \n", cmd_line->ref_colour);
    }
  db_graph_detect_vars_given_lists_of_colours(fp,cmd_line->max_var_len,db_graph, 
					      cmd_line->detect_bubbles1_first_colour_list, 
					      cmd_line->num_colours_in_detect_bubbles1_first_colour_list,
					      cmd_line->detect_bubbles1_second_colour_list, 
					      cmd_line->num_colours_in_detect_bubbles1_second_colour_list,
					      &bubble_caller_condition, print_appropriate_extra_var_info,
					      cmd_line->exclude_ref_bubbles, get_col_ref, get_cov_ref, 
					      cmd_line->apply_model_selection_at_bubbles, mod_sel_criterion, 
					      model_info);

  fclose(fp);
  if (cmd_line->high_diff==true)
    {
      free(allele_balances);
      free_covg_array(working_ca);
    }
}






int main(int argc, char **argv)
{
  timestamp();

  // VERSION_STR is passed from the makefile -- usually last commit hash
  printf("Starting Cortex, version %d.%d.%d.%d"VERSION_STR"\n",
         VERSION, SUBVERSION, SUBSUBVERSION, SUBSUBSUBVERSION);

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


  //***************************************************************************
  //define local functions
  //***************************************************************************
  /*
  void apply_reset_to_all_edges(dBNode* node, Orientation or, Nucleotide nuc)
  {
    int j;
    for (j=0; j<NUMBER_OF_COLOURS; j++)
      {
	reset_one_edge(node, or, nuc, j);
      }
  }
  void apply_reset_to_all_edges_2(dBNode* node )
  {
    int j;
    for (j=0; j<NUMBER_OF_COLOURS; j++)
      {
	      db_node_reset_edges(node, j);
      }
  }
  */


  CovgArray* working_ca_for_median=NULL;
  if (cmd_line->print_median_covg_only==true)
    {
      working_ca_for_median = alloc_and_init_covg_array(cmd_line->max_var_len);//will die if fails to alloc
    }

  void print_appropriate_extra_variant_info(VariantBranchesAndFlanks* var, FILE* fp)
  {
    if (cmd_line->print_colour_coverages==true)
      {
	print_standard_extra_info(var, fp);
      }
    else if (cmd_line->print_median_covg_only==true)
      {
	print_median_covg_extra_info(var, working_ca_for_median, fp);
      }
    else
      {
	print_no_extra_info(var, fp);
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
  
  printf("Maximum k-mer size (compile-time setting): %i\n", max_kmer_size);

  if(cmd_line->kmer_size > max_kmer_size)
  {
    die("k-mer size is too big [%i]!",cmd_line->kmer_size);
  }

  printf("Actual K-mer size: %d\n", cmd_line->kmer_size);

  //Create the de Bruijn graph/hash table
  int max_retries=15;
  db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
  if (db_graph==NULL)
    {
      die("Giving up - unable to allocate memory for the hash table\n");
    }
  printf("Hash table created, number of buckets: %d\n",1 << hash_key_bits);



  GraphInfo* db_graph_info=graph_info_alloc_and_init();//will exit it fails to alloc.

  // input data:
  if (cmd_line->input_seq==true)
    {
      if (cmd_line->subsample==true)
	{
	  printf("User has specified subsampling %f of the reads, so as we parse the fasta/q or BAM  we will\n load each read with probability %f.\n",
		 cmd_line->subsample_propn, cmd_line->subsample_propn);
	}

      if (strcmp(cmd_line->se_list, "")!=0)
	{
	  fprintf(stdout,"Input file of single ended data filenames: %s\n",
		  cmd_line->se_list);
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
	  printf("No paired-end data\n");
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
	  printf("Removing duplicates from the paired end files when both mates "
           "start with the same kmer\n");
	}
      else
	{
	  printf("No PCR duplicate removal\n");
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
    
    printf("Sequence data loaded\n");
    printf("Total bases parsed:%llu\n", num_bases_parsed);
    printf("Total bases passing filters and loaded into graph:%llu\n",
           num_bases_loaded);
    printf("Mean read length after filters applied:%lu\n", mean_contig_length);
    
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
  /*
  printf("The following mod is for the sims for the paper only: i need the graphinfo to contain read lengths for the fasta based sim binaries\n");

  for (j=0; j<NUMBER_OF_COLOURS; j++)
    {
      if (db_graph_info->mean_read_length[j]==0)
	{
	  graph_info_set_mean_readlen(db_graph_info, j, cmd_line->max_read_length);
	}
    }
  */


  if (cmd_line->successively_dump_cleaned_colours==false)
    {      
      printf("Total kmers in table: %qd\n", hash_table_get_unique_kmers(db_graph));	  
      printf("The following is a summary of the data that has been loaded, immediately after loading (prior to any error cleaning, calling etc)\n");
      printf("****************************************\n");
      printf("SUMMARY:\nColour\tSampleID\tMeanReadLen\tTotalSeq\tErrorCleaning\tLowCovSupsThresh\tLowCovNodesThresh\tPoolagainstWhichCleaned\n");
      
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  //colour sampleid readlen totalseq 
	  printf("%d\t%s\t%d\t%" PRIu64 "\t", 
		 j, 
		 db_graph_info->sample_ids[j],
		 db_graph_info->mean_read_length[j], 
		 db_graph_info->total_sequence[j]);
	  if ( (db_graph_info->cleaning[j]->tip_clipping==false)
	       &&
	       (db_graph_info->cleaning[j]->remv_low_cov_sups==false)
	       &&
	       (db_graph_info->cleaning[j]->remv_low_cov_nodes==false)
	       &&
	       (db_graph_info->cleaning[j]->cleaned_against_another_graph==false)
	       )
	    {
	      printf("UNCLEANED\t");
	    }
	  else
	    {
	      printf("CLEANED\t");
	    }
	  printf("%d\t%d\t%s\n",
		 db_graph_info->cleaning[j]->remv_low_cov_sups_thresh,
		 db_graph_info->cleaning[j]->remv_low_cov_nodes_thresh,
		 db_graph_info->cleaning[j]->name_of_graph_against_which_was_cleaned);
	}
      printf("****************************************\n");
    }



  if (cmd_line->health_check==true)
    {
      timestamp();
      printf("Run health check\n");
      db_graph_health_check(false, db_graph);
      printf("End of health check\n");
      timestamp();
    }


  if (cmd_line->remv_low_covg_sups_threshold!=-1)
    {
      printf("Clip tips first\n");
      db_graph_clip_tips_in_union_of_all_colours(db_graph);

      printf("Remove low coverage supernodes covg (<= %d) \n", cmd_line->remv_low_covg_sups_threshold);
      db_graph_remove_errors_considering_covg_and_topology(cmd_line->remv_low_covg_sups_threshold,
							   db_graph, 
							   &element_get_covg_union_of_all_covgs, 
							   &element_get_colour_union_of_all_colours,
							   &apply_reset_to_specific_edge_in_union_of_all_colours, 
							   &apply_reset_to_all_edges_in_union_of_all_colours,
							   cmd_line->max_var_len);
      timestamp();
      printf("Error correction done\n");
      int z;
      for (z=0; z<NUMBER_OF_COLOURS; z++)
	{
	  graph_info_set_remv_low_cov_sups(db_graph_info, z, cmd_line->remv_low_covg_sups_threshold);
	  graph_info_set_tip_clipping(db_graph_info, z);
	}
    }
  else if (cmd_line->remove_low_coverage_nodes==true)
    {
      timestamp();
      printf("Start to to remove nodes with covg (in union of all colours)  <= %d\n", cmd_line->node_coverage_threshold);
      db_graph_remove_low_coverage_nodes_ignoring_colours(cmd_line->node_coverage_threshold, db_graph);

      timestamp();
      printf("Error correction done\n");
      int z;
      for (z=0; z<NUMBER_OF_COLOURS; z++)
	{
	  graph_info_set_remv_low_cov_nodes(db_graph_info, z, cmd_line->node_coverage_threshold);
	}
      
    }

  if (cmd_line->health_check==true)
    {
      timestamp();
      printf("Run health check on cleaned graph\n");
      db_graph_health_check(false, db_graph);
      printf("Health check done\n");
      timestamp();
    }

  if ( (cmd_line->dump_binary==true)     
       ||
       (cmd_line->subsample==true) )
    {
      if (cmd_line->input_seq==true)
	{
	  //dump single colour
	  timestamp();
	  printf("Input data was fasta/q, so dump single colour binary file: %s\n", cmd_line->output_binary_filename);
	  db_graph_dump_single_colour_binary_of_colour0(cmd_line->output_binary_filename, &db_node_check_status_not_pruned,
							db_graph, db_graph_info, BINVERSION);
	  timestamp();
	  printf("Binary dumped\n");
	  
	}
      else
	{
	  timestamp();
	  printf("Dump multicolour binary with %d colours (compile-time setting)\n", NUMBER_OF_COLOURS);
	  db_graph_dump_binary(cmd_line->output_binary_filename, &db_node_check_status_not_pruned,db_graph, db_graph_info, BINVERSION);
	  timestamp();
	  printf("Binary dumped\n");
	}
    }

  if (cmd_line->do_err_correction==true)
    {
      timestamp();
      printf("Error correct against a population graph\n");
      boolean reverse_comp_reads_to_match_strand=false;
      if (strcmp(cmd_line->ref_chrom_fasta_list, "")!=0)
	{
	  printf("Since you gave a reference fasta list, will mark plus strand on graph and print reads in that strand\n");
	  FILE* list_fp = fopen(cmd_line->ref_chrom_fasta_list, "r");
	  if (list_fp==NULL)
	    {
	      printf("Can't open your list of reference fasta: %s\n", cmd_line->ref_chrom_fasta_list);
	    }
	  StrBuf *line = strbuf_new();
	  while(strbuf_reset_readline(line, list_fp))
	    {
	      strbuf_chomp(line);
	      if(strbuf_len(line) > 0)
		{
		  read_ref_fasta_and_mark_strand(line->buff, db_graph);
		}
	    }
	  strbuf_free(line);
	  fclose(list_fp);
	  reverse_comp_reads_to_match_strand=true;
	  printf("Finished marking strand. Start correcting\n");
	  timestamp();
	}
      error_correct_list_of_files(cmd_line->err_correction_filelist, cmd_line->quality_score_threshold, cmd_line->quality_score_offset,
				  db_graph, cmd_line->err_correction_policy,
				  cmd_line->max_read_length, cmd_line->err_correction_suffix,
				  cmd_line->err_correction_outdir->buff,
				  cmd_line->do_greedy_padding, 
				  cmd_line->greedy_pad,
				  reverse_comp_reads_to_match_strand);
      timestamp();
      printf("Error correction done\n");
    }
  if (cmd_line->print_supernode_fasta==true)
    {

      timestamp();
      printf("Print contigs(supernodes) in the graph created by the union of all colours.\n");
      
      db_graph_print_supernodes_defined_by_func_of_colours(cmd_line->output_supernodes, "", 
							   cmd_line->max_var_len,// max_var_len is the public face of maximum expected supernode size
							   db_graph, 
							   &element_get_colour_union_of_all_colours, 
							   &element_get_covg_union_of_all_covgs, 
							   &print_appropriate_extra_supernode_info);


      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Supernodes dumped\n");
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

  // DETECT BUBBLES

  if (cmd_line->detect_bubbles1==true)
    {
      timestamp();
      printf("Start first set of bubble calls\n");
      run_bubble_calls(cmd_line, 1, db_graph, &print_appropriate_extra_variant_info,
                       &get_colour_ref, &get_covg_ref, &model_info);

      //unset the nodes marked as visited, but not those marked as to be ignored
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Detect Bubbles 1, completed\n");
    }


  if (cmd_line->do_genotyping_of_file_of_sites==true)
    {
      timestamp();

      printf("Genotype the calls in this file %s\n", cmd_line->file_of_calls_to_be_genotyped);

      run_genotyping(cmd_line, db_graph, &print_appropriate_extra_variant_info,
                     //&get_colour_ref, &get_covg_ref, db_graph_info,
                     &model_info);

      //unset the nodes marked as visited, but not those marked as to be ignored
      timestamp();
      printf("Genotyping completed\n");
    }

  if (cmd_line->estimate_copy_num==true)
    {
      // add code to do copy number estimation here!!!
    }

  if (cmd_line->make_pd_calls==true)
    {
      timestamp();
      printf("Run Path-Divergence Calls\n");
      run_pd_calls(cmd_line, db_graph, &print_appropriate_extra_variant_info, &model_info);
      //hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph);	
      timestamp();
      printf("Finished Path Divergence calls\n");
    }
  if (cmd_line->align_given_list==true)
    {
      timestamp();
      printf("Start aligning the fasta/q listed in this file: %s\n", cmd_line->list_fastaq_to_align);
      int array_of_colours[NUMBER_OF_COLOURS];
      int j;
      char* array_of_colournames[NUMBER_OF_COLOURS];
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  array_of_colours[j]=j;
	  array_of_colournames[j]=(char*)malloc(sizeof(char) * 50);
	  if (array_of_colournames[j]==NULL)
	    {
	      die("Severe lack of memory. Cannot even allocate 50 chars. Give up\n");
	    }
	  sprintf(array_of_colournames[j], "colour_%d", j);
	}
      align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours(cmd_line->format_of_files_to_align, cmd_line->list_fastaq_to_align,
								       cmd_line->max_read_length, array_of_colours, array_of_colournames,
								       NUMBER_OF_COLOURS,db_graph,cmd_line->quality_score_offset,
								       false, NULL, NULL, cmd_line->dump_aligned_overlap_binary);
      for (j=0; j<NUMBER_OF_COLOURS; j++)
	{
	  free(array_of_colournames[j]);
	}
      printf("Completed alignment of fasta/q to graph to print coverages in all colours\n");
      printf("Dumping a binary, %s,  of all the nodes which were hit by the alignment process\n", cmd_line->output_aligned_overlap_binname);
      sprintf(tmp_dump, "%s.temporary_delete_me", cmd_line->output_aligned_overlap_binname);
      printf("In the process we have to create a temporary file, %s, which you can/should delete when cortex has completed\n", tmp_dump);

      //num_kmers_dumped_after_alignment = db_graph_dump_binary(tmp_dump, &db_node_check_status_to_be_dumped, db_graph, db_graph_info, BINVERSION);
      db_graph_dump_binary(tmp_dump, &db_node_check_status_to_be_dumped,
                           db_graph, db_graph_info, BINVERSION);

      hash_table_traverse(&db_node_action_set_status_of_unpruned_to_none, db_graph);
    }
  if (cmd_line->get_pan_genome_matrix==true)
    {
      timestamp();
      printf("Print pangenome matrix. For each read(gene) in this fasta (%s)\nprint the percent of kmers which are present in each sample\n",
	     cmd_line->pan_genome_genes_fasta);
      print_percent_agreement_for_each_colour_for_each_read(cmd_line->pan_genome_genes_fasta, 
							    cmd_line->max_read_length, db_graph,
							    db_graph_info->sample_ids);
      timestamp();
    }
  if (cmd_line->print_colour_overlap_matrix==true)
    {
      timestamp();
      printf("Do a direct graph comparison of each colour x,y,z.. listed in --colour_subgraph_overlap_matrix a,b,c../x,y,z... against the colours a,b,c.. (process is symmetrical)\n");
      
      db_graph_print_colour_overlap_matrix(cmd_line->colour_overlap_first_colour_list,
					   cmd_line->num_colours_in_colour_overlap_first_list,
					   cmd_line->colour_overlap_second_colour_list,
					   cmd_line->num_colours_in_colour_overlap_second_list,
					   db_graph);
      printf("\nCompleted graph overlap matrix printing\n");
    }


  if (cmd_line->genotype_complex_site==true)
    {
      printf("Genotyping a complex site with %d known alleles\n", cmd_line->num_alleles_of_site);
      printf("Of all possible genotypes, we will calculate likelihoods for those numbered %d to %d\n",
	     cmd_line->first_genotype_to_calc_likelihoods_for, cmd_line->last_genotype_to_calc_likelihoods_for);
      printf("We do this for these samples (which we assume are diploid): colours ");
      int k;
      for (k=0; k<cmd_line->num_colours_to_genotype; k++)
	{
	  printf("%d,", cmd_line->list_colours_to_genotype[k]);
	}
      printf("\n");

      //sanity check before we start
      //check that we have read lengths for the colours we want to genotype
      for (k=0; k<cmd_line->num_colours_to_genotype; k++)
	{
	  if (db_graph_info->mean_read_length[cmd_line->list_colours_to_genotype[k]] < (uint32_t) db_graph->kmer_size )
	    {
	      die(
		  "This will not work. If you scroll up to the summary of read-lengths and covgs\n"
		  "in your colours, you will see that at least one of the colours you want to\n"
		  "genotype has mean read length < kmer_size. This should be impossible - this is just a paranoid error check.\n");
	    }

	  if (db_graph_info->total_sequence[cmd_line->list_colours_to_genotype[k]]==0)
	    {
	      die(
		  "This will not work. If you scroll up to the summary of read-lengths and covgs\n"
		  "in your colours, you will see that at least one of the colours you want to\n"
		  "genotype has total sequence 0. \n");
	    }
	}

      double* current_max_lik_array         = alloc_ML_results_array(cmd_line->num_colours_to_genotype);
      double* current_max_but_one_lik_array = alloc_ML_results_array(cmd_line->num_colours_to_genotype);
      char** name_current_max_lik_array             = alloc_ML_results_names_array(cmd_line->num_colours_to_genotype);
      char** name_current_max_but_one_lik_array     = alloc_ML_results_names_array(cmd_line->num_colours_to_genotype);

      calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site(cmd_line->list_colours_to_genotype, cmd_line->num_colours_to_genotype,
										       cmd_line->colour_of_reference_with_site_excised,
										       cmd_line->num_alleles_of_site,
										       cmd_line->first_genotype_to_calc_likelihoods_for,
										       cmd_line->last_genotype_to_calc_likelihoods_for,
										       cmd_line->max_var_len, cmd_line->fasta_alleles_for_complex_genotyping,
										       cmd_line->assump_for_genotyping,
										       current_max_lik_array, current_max_but_one_lik_array,
										       name_current_max_lik_array, name_current_max_but_one_lik_array,
										       true, &model_info, db_graph, cmd_line->working_colour1, cmd_line->working_colour2,
										       cmd_line->using_1net, cmd_line->using_2net,
										       cmd_line->min_acceptable_llk
										       );


    }
  if (cmd_line->estimate_genome_complexity==true)
    {
      int num_reads_used_in_estimate=0;
      double g = estimate_genome_complexity(db_graph, cmd_line->fastaq_for_estimating_genome_complexity,
					    true, 1,cmd_line->max_read_length, cmd_line->format_of_input_seq,
					    cmd_line->quality_score_offset, &num_reads_used_in_estimate);
      printf("We estimate genome complexity at k=%d (for SNPs) as %f\n", db_graph->kmer_size, g);
      printf("This estimate used a sample of %d high-quality reads\n", num_reads_used_in_estimate);

    }

  if (cmd_line->print_novel_contigs==true)
    {
      timestamp();
      printf("Start to search for and print novel contigs\n");
      printf("Definition of novel: contig must lie in union of these colours: ");
      int j;
      for (j=0; j<cmd_line->numcols_novelseq_colours_search; j++)
	{
	  printf("%d,", cmd_line->novelseq_colours_search[j]);
	}
      printf("\n and at least %d percent of the kmers in this contig (excluding first and last) must have zero coverage in the union of these colours: ", cmd_line->novelseq_min_percentage_novel);
      for (j=0; j<cmd_line->numcols_novelseq_colours_avoid; j++)
	{
	  printf("%d,", cmd_line->novelseq_colours_avoid[j]);
	}
      printf("\nAlso contig must be at least %d bp long\n", cmd_line->novelseq_contig_min_len_bp);
      db_graph_print_novel_supernodes(cmd_line->novelseq_outfile, cmd_line->max_var_len, db_graph, 
				      cmd_line->novelseq_colours_search, cmd_line->numcols_novelseq_colours_search,
				      cmd_line->novelseq_colours_avoid, cmd_line->numcols_novelseq_colours_avoid,
				      cmd_line->novelseq_contig_min_len_bp, cmd_line->novelseq_min_percentage_novel,
				      &print_appropriate_extra_supernode_info);
      timestamp();
      printf("Finished printing novel contigs\n");
      
    }

  
  hash_table_free(&db_graph);
  if (cmd_line->print_median_covg_only==true)
    {
      free_covg_array(working_ca_for_median);
    }

  timestamp();

  if (cmd_line->dump_aligned_overlap_binary==true)
    {

      //reload the binary you dumped, clean off the edges, and then dump again.
      //malloc a new hash table. Only needs to be as large as you need.
      //float s = log(num_kmers_dumped_after_alignment/bucket_size)/log(2); 
      // now 2^s = num_kmers_dumped_after_alignment/bucket_size, so s is my height
      dBGraph* db_graph2;
      printf("Reload and fix dangling edges in the temporary binary we have created\n");
      db_graph2 = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
      if (db_graph2==NULL)
	{
	  die("Cortex has nearly finished. It's done everything you asked it to do, and has dumped a binary of the overlap of your alignment with the graph. However, by \"ripping out\" nodes from the main graph, that dumped binary now has edges pointing out to nodes that are not in the binary. So the idea is that we have deallocated the main graph now, and we were going to load the dumped binary, clean it up and re-dump it. However that has failed, becaause we could not malloc the memory to do it - the most likely reason is that someone else is sharing your server and their mempory use has gone up.\n");
	}
      
      //this is all for the API - we wont use this info
      GraphInfo* temp_info = graph_info_alloc_and_init();

      int num_c;//number of colours in binary      
      load_multicolour_binary_from_filename_into_graph(tmp_dump,db_graph2, temp_info, &num_c);


      //this is why we are going to all this bother - cleaning edges
      db_graph_clean_orphan_edges(db_graph2);
      db_graph_dump_binary(cmd_line->output_aligned_overlap_binname, 
			   &db_node_condition_always_true,
			   db_graph2, db_graph_info, BINVERSION);//deliberately using original graph info - we want the same info
      hash_table_free(&db_graph2);
      graph_info_free(temp_info);
    }
  


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
