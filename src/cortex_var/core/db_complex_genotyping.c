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

#include <math.h>
#include <db_variants.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <maths.h>
#include <gsl_sf_gamma.h>
#include <db_complex_genotyping.h>
#include <string.h>
#include <model_selection.h>
#include <stdio.h>
#include <stdlib.h>
#include <global.h>
#include <file_reader.h>

// Lowest log likelihood, at which we just cutoff
int MIN_LLK = -999999999;

MultiplicitiesAndOverlapsOfBiallelicVariant*  alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(int len_allele1, int len_allele2)
{
  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv = (MultiplicitiesAndOverlapsOfBiallelicVariant*) malloc(sizeof(MultiplicitiesAndOverlapsOfBiallelicVariant));

  if (mobv==NULL)
    {
      printf("Ridiculous - failed to alloc an array of 4 ints\n");
      exit(1);
    }
  mobv->mult11=(int*) malloc( (sizeof(int))*len_allele1);
  mobv->mult12=(int*) malloc( (sizeof(int))*len_allele1);
  mobv->mult21=(int*) malloc( (sizeof(int))*len_allele2);
  mobv->mult22=(int*) malloc( (sizeof(int))*len_allele2);

  if ( (mobv->mult11==NULL) || (mobv->mult12==NULL) || (mobv->mult21==NULL) || (mobv->mult22==NULL) )
    {
      return NULL;
    }
  mobv->len1 = len_allele1;
  mobv->len2 = len_allele2;

  return mobv;

}


void dealloc_MultiplicitiesAndOverlapsOfBiallelicVariant(MultiplicitiesAndOverlapsOfBiallelicVariant* mobv)
{
  free(mobv->mult11);
  free(mobv->mult12);
  free(mobv->mult21);
  free(mobv->mult22);
  free(mobv);
}

//does not change the lengths - these are set when you malloc
void reset_MultiplicitiesAndOverlapsOfBiallelicVariant(MultiplicitiesAndOverlapsOfBiallelicVariant* mobv)
{
  if (mobv==NULL)
    {
      return;
    }
  else
    {
      int i;
      for (i=0; i< mobv->len1; i++)
	{
	  mobv->mult11[i]=0;
	  mobv->mult12[i]=0;
	}
      for (i=0; i< mobv->len2; i++)
	{
	  mobv->mult21[i]=0;
	  mobv->mult22[i]=0;
	}
    }
  return;
}

//Utility function - only exported so I can test it.
void initialise_multiplicities_of_allele_nodes_wrt_both_alleles(VariantBranchesAndFlanks* var, MultiplicitiesAndOverlapsOfBiallelicVariant* mult,
								boolean only_count_nodes_with_edge_in_specified_colour_func,
								Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*) )
{

  void get_mult(dBNode** br_src, int len_br_src, dBNode** br_target, int len_br_target, int* mult_array)
  {
    int i,j;
    int count_occurrences=0; //will be number of things we have put in this array

    for (i=0; i<len_br_src ; i++)
      {
	count_occurrences=0;
	
	for (j=0 ; j<len_br_target; j++)
	  {
	    //if i-th and j-th elements are the same, AND they exist in the colour (or function of colours) we are interested in
	    if (db_node_addr_cmp(&br_src[i], &br_target[j])==0 )
	      {
		if ( (only_count_nodes_with_edge_in_specified_colour_func==true) && 
		     (!db_node_is_this_node_in_subgraph_defined_by_func_of_colours(br_src[i], get_colour)) )
		  {
		    //does not count if node does not exist in the specified subgraph
		  }
		else
		  {
		    count_occurrences++;
		  }
	      }
	  }
	
	mult_array[i]=count_occurrences;
      }
  }
  get_mult(var->one_allele, var->len_one_allele, var->one_allele, var->len_one_allele,  mult->mult11);
  get_mult(var->other_allele, var->len_other_allele, var->other_allele, var->len_other_allele,  mult->mult22);
  get_mult(var->one_allele, var->len_one_allele, var->other_allele, var->len_other_allele,  mult->mult12);
  get_mult(var->other_allele, var->len_other_allele, var->one_allele, var->len_one_allele,  mult->mult21);
}

void utility_set_one_array_equal_to_another(dBNode** src, int len, dBNode** target)
{
  int i;
  for (i=0; i<len; i++)
    {
      target[i]=src[i];
    }
}

//Utility function - only exported so I can test it.
void improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(VariantBranchesAndFlanks* var, MultiplicitiesAndOverlapsOfBiallelicVariant* mult,
									 boolean only_count_nodes_with_edge_in_specified_colour_func,
									 Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*),
									 int working_colour1, int working_colour2, dBGraph* db_graph)
{
  //wipe the working colours:
  db_graph_wipe_colour(working_colour1, db_graph);
  db_graph_wipe_colour(working_colour2, db_graph);
  //walk through allele1, and as you do so, add multicplicity counts to working_colour1
  int i;
  for (i=0; i<var->len_one_allele; i++)
    {
      db_node_increment_coverage(var->one_allele[i], individual_edge_array, working_colour1);
    }
  //same for allele2
  for (i=0; i<var->len_other_allele; i++)
    {
      db_node_increment_coverage(var->other_allele[i], individual_edge_array, working_colour2);
    }
  //now, as you walk through allele1, for each node, you can see how many times it occurs in allele1 and in allele2
  for (i=0; i<var->len_one_allele; i++)
    {
      mult->mult11[i] += db_node_get_coverage(var->one_allele[i], individual_edge_array, working_colour1);
      mult->mult12[i] += db_node_get_coverage(var->one_allele[i], individual_edge_array, working_colour2);
    }

  for (i=0; i<var->len_other_allele; i++)
    {
      mult->mult21[i] += db_node_get_coverage(var->other_allele[i], individual_edge_array, working_colour1);
      mult->mult22[i] += db_node_get_coverage(var->other_allele[i], individual_edge_array, working_colour2);
    }

}



// WARNING - the implicit assumption here is that the graph/hash table ONLY contains nodes in the union of all known alleles for 
//           the site in question.
// Pass in the current_max and current_max_but_one log likelihoods, and as soon as this drops below both, we abort mission.
// and return -99999999
// var_mults must be pre-allocated ***and pre-initialised***

double calc_log_likelihood_of_genotype_with_complex_alleles(VariantBranchesAndFlanks* var,
							    char* name_of_this_genotype,
							    MultiplicitiesAndOverlapsOfBiallelicVariant* var_mults,
							    GraphAndModelInfo* model_info,
							    int colour_indiv, 
							    int colour_ref_minus_our_site, dBGraph* db_graph,
							    int* working_array_self, int* working_array_shared,
							    double* current_max_lik, double* current_max_but_one_lik,
							    char* current_max_lik_name, char* current_max_but_one_lik_name,
							    AssumptionsOnGraphCleaning assump)
//assume the  2 working arrays have length = max read length
{
    
  // *********************************************************************************************************************************************************
  // Likelihood of this genotype is given by prob(see number_errors of errors in sequence with length=sum of lengths of two alleles)
  //                                         * Product ( dpois(observed covg, expected covg) )     (product over subchunks, splitting where alleles overlap)
  // *********************************************************************************************************************************************************

  Edges element_get_colour_indiv(const Element* e)
  {
    return e->individual_edges[colour_indiv];
  }
  int element_get_covg_indiv(const dBNode* e)
  {
    return db_node_get_coverage(e, individual_edge_array, colour_indiv);
  }


  // 1. Count number of errors:

  int check_covg_in_ref_with_site_excised(dBNode* n)
  {
    if (colour_ref_minus_our_site==-1)
      {
	return 0;
      }
    else
      {
	return db_node_get_coverage(n, individual_edge_array, colour_ref_minus_our_site);
      }
  }
  
  void get_number_nodes_outside_our_expected_two_alleles_with_covg( dBNode* e, int* total)
  {
    if (db_node_check_status(e, in_desired_genotype)==true)
      {
      }
    else if ( (db_node_get_coverage(e, individual_edge_array,colour_indiv)>0) && (check_covg_in_ref_with_site_excised(e)==0) )
      {
	*total +=1;
      }
    return;
  }



  int number_errors=0;//get number of nodes that are (in the union of all known alleles but) outside our genotype alleles, and are not in the rest of the genome
  int total_len_alleles = (var->len_one_allele + var->len_other_allele);
  hash_table_traverse_passing_int(&get_number_nodes_outside_our_expected_two_alleles_with_covg  ,db_graph, &number_errors);

  // Errors on union of two alleles occur as a Poisson process with rate = sequencing_error_rate_per_base * kmer * length of two alleles
  double poisson_error_rate;
  if (assump==AssumeUncleaned)
    {
      poisson_error_rate = (model_info->seq_error_rate_per_base) * (db_graph->kmer_size) *total_len_alleles;
    }
  else if (assump==AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice)
    {
      poisson_error_rate = (model_info->seq_error_rate_per_base)*(model_info->seq_error_rate_per_base) * (db_graph->kmer_size) *total_len_alleles;
    }
  double log_prob_error = 0;


  //log of (exp(-poisson_error_rate) (poisson_error_rate^number_errors)/(number_errors!) ) = -poisson_error_rate + number_errors*log(poisson_error_rate) - log(number_errors!)
  if (number_errors>0)
    {
	  log_prob_error = -poisson_error_rate + number_errors*log(poisson_error_rate) - gsl_sf_lnfact(number_errors);
    }

  if (log_prob_error < *current_max_but_one_lik)
    {
      return MIN_LLK;//abort - this genotype cannot be the max likelihood genotype
    }


  // 2.  Get sum of log likelihoods of all the subchunks of the two alleles

  // walk along allele 1. For sections unique to allele1, work out log likelihood with eff haploid covg
  // for sections shared with allele2, use eff diploid covg
  // and finally, walk along allele2, and consider only sections unique to allele2



  int max;
  if (var->len_one_allele > var->len_other_allele)
    {
      max = var->len_one_allele;
    }
  else
    {
      max = var->len_other_allele;
    }

  int working_array_self_count=0;
  int working_array_shared_count=0;
  double log_prob_data=0;
  double hap_D  = (double) (model_info->ginfo->total_sequence[colour_indiv])/(double) (2*model_info->genome_len) ;
  double hap_D_over_R = hap_D/(model_info->ginfo->mean_read_length[colour_indiv]);

  //break allele into intervals that only lie on allele1, and those shared with allele2.
  int k=0;
  while (k < var->len_one_allele)
      {
	if (var_mults->mult12[k]>0) //this node occurs  >0 times in the other allele
	  {
	    if (working_array_self_count>0)
	      {
		boolean too_short = 0;
		int covg_on_self_in_this_chunk = count_reads_on_allele_in_specific_colour_given_array_of_cvgs(working_array_self, working_array_self_count, 
													      &too_short);
		if (too_short==false)
		  {
 		    // add log dpois ( hap_D_over_R * working_array_self_count , covg_on_self_in_this_chunk)
		    log_prob_data += -hap_D_over_R * working_array_self_count + covg_on_self_in_this_chunk * log(hap_D_over_R * working_array_self_count) - gsl_sf_lnfact(covg_on_self_in_this_chunk);
		    if (log_prob_error + log_prob_data < *current_max_but_one_lik)
		      {
			return MIN_LLK;
		      }
		  }
	      }

	    working_array_self_count=0;
	    working_array_shared_count=0;
	    while ( (var_mults->mult12[k] > 0) && 
		    (k < var->len_one_allele) && 
		    (check_covg_in_ref_with_site_excised(var->one_allele[k])==0 )//automatically handles case when there is no colour for ref-minus-site     
		  )
	      {
		// this node might happen 2 times on allele 1 and 3 times on allele2. in total 5 times, annd here it counts for 2 of them
		//working_array_shared[working_array_shared_count]=covg_allele1[k]* 2/(var_mults->mult11[k] + var_mults->mult12[k])  ;<<<<<<<<<<<< this is what I had before (when I did NA19240)
		working_array_shared[working_array_shared_count]= db_node_get_coverage(var->one_allele[k], individual_edge_array, colour_indiv)/(var_mults->mult11[k] + var_mults->mult12[k])  ;

		k++;
		working_array_shared_count++;
	      }
	    while ( (k<var->len_one_allele) && ( db_node_get_coverage( (var->one_allele)[k], individual_edge_array, colour_ref_minus_our_site)>0 ) )
	      {
		k++;
	      }

	  }
	if (var_mults->mult12[k]==0)
	  {
	    //start of contiguous chunk which occurs only in second allele, not the fisrt
	    if (working_array_shared_count>0)
              {
		boolean too_short = 0;
		int covg_on_shared_in_this_chunk = count_reads_on_allele_in_specific_colour_given_array_of_cvgs(working_array_shared, working_array_shared_count, 
														&too_short);
		if (too_short==false)
		  {
		    //            NOTE  2*  - diploid covg
 		    // add log dpois ( 2*hap_D_over_R * working_array_shared_count , covg_on_shared_in_this_chunk)
		    log_prob_data += -2*hap_D_over_R * working_array_shared_count + covg_on_shared_in_this_chunk * log(2*hap_D_over_R * working_array_shared_count) - gsl_sf_lnfact(covg_on_shared_in_this_chunk);
		    if (log_prob_error + log_prob_data < *current_max_but_one_lik)
		      {
			return MIN_LLK;
		      }
		  }
              }

	    working_array_shared_count=0;
	    working_array_self_count=0;
	    
	    while ( (var_mults->mult12[k]==0) && 
		    (k < var->len_one_allele) &&
		    (check_covg_in_ref_with_site_excised(var->one_allele[k])==0)
		    )
	      {
		if (var_mults->mult11[k]==0)
		  {
		    printf("At k= %d get zero\n", k);
		  }
		working_array_self[working_array_self_count]=db_node_get_coverage(var->one_allele[k], individual_edge_array, colour_indiv)/var_mults->mult11[k];
		k++;
		working_array_self_count++;
	      }

	    while ( (k < var->len_one_allele) &&( check_covg_in_ref_with_site_excised(var->one_allele[k])>0 ) )

	      { 
		k++;
              }

	  }
      }
  //may have exited loop with some left-over data in arrays, so:  
  if (working_array_self_count>0)
    {
      boolean too_short = 0;
      int covg_on_self_in_this_chunk = count_reads_on_allele_in_specific_colour_given_array_of_cvgs(working_array_self, working_array_self_count,
												    &too_short);
      if (too_short==false)
	{
	  // add log dpois ( hap_D_over_R * working_array_self_count , covg_on_self_in_this_chunk)
	  log_prob_data += -hap_D_over_R * working_array_self_count + covg_on_self_in_this_chunk * log(hap_D_over_R * working_array_self_count) - gsl_sf_lnfact(covg_on_self_in_this_chunk);
	  if (log_prob_error + log_prob_data < *current_max_but_one_lik)
	    {
	      return MIN_LLK;
	    }
	}

    }

  if (working_array_shared_count>0)
    {
      boolean too_short = 0;
      int covg_on_shared_in_this_chunk = count_reads_on_allele_in_specific_colour_given_array_of_cvgs(working_array_shared, working_array_shared_count, 
												      &too_short);
      if (too_short==false)
	{
	  //            NOTE  2*  - diploid covg
	  // add log dpois ( 2*hap_D_over_R * working_array_shared_count , covg_on_shared_in_this_chunk)
	  log_prob_data += -2*hap_D_over_R * working_array_shared_count + covg_on_shared_in_this_chunk * log(2*hap_D_over_R * working_array_shared_count) - gsl_sf_lnfact(covg_on_shared_in_this_chunk);
	  if (log_prob_error + log_prob_data < *current_max_but_one_lik)
	    {
	      return MIN_LLK;
	    }
	}
    }



  //now walk along the other allele, only considering the intervals that are NOT shared
  k=0;
  working_array_self_count=0;
  working_array_shared_count=0;
  
  while (k < var->len_other_allele)
    {
      if (var_mults->mult21[k]>0) //this node occurs  >0 times in the other allele
	{
	  if (working_array_self_count>0)
	    {
	      boolean too_short = 0;
	      int covg_on_self_in_this_chunk = count_reads_on_allele_in_specific_colour_given_array_of_cvgs(working_array_self, working_array_self_count, 
													    &too_short);
	      if (too_short==false)
		{
		  // add log dpois ( hap_D_over_R * working_array_self_count , covg_on_self_in_this_chunk)
		  log_prob_data += -hap_D_over_R * working_array_self_count 
		    + covg_on_self_in_this_chunk * log(hap_D_over_R * working_array_self_count) - gsl_sf_lnfact(covg_on_self_in_this_chunk);
		  
		  if (log_prob_error + log_prob_data < *current_max_but_one_lik)
		    {
		      return MIN_LLK;
		    }
		}
	    }
	  
	  working_array_self_count=0;
	  k++;
	}
      else
	{
	  while ( (check_covg_in_ref_with_site_excised(var->other_allele[k])==0) &&
		  (k < var->len_other_allele) &&
		  (var_mults->mult21[k]==0)
		  )
	    {
	      working_array_self[working_array_self_count]=db_node_get_coverage(var->other_allele[k], individual_edge_array, colour_indiv)/var_mults->mult22[k];
	      k++;
	      working_array_self_count++;
	    }
	  
	  while (  (k<var->len_other_allele) &&
		   (check_covg_in_ref_with_site_excised(var->other_allele[k])>0)
		   )
	    {
	      k++;
	    } 
	}
      
    }

  if (log_prob_error + log_prob_data < *current_max_but_one_lik)
    {
      return MIN_LLK;
    }
  else
    {
      if (log_prob_data + log_prob_error > (*current_max_lik) )
	{
	  //this is the new ML genotype
	  //sort out new names for top and next best genotype:
	  current_max_but_one_lik_name[0]='\0';
	  strcat(current_max_but_one_lik_name, current_max_lik_name);
	  current_max_lik_name[0]='\0';
	  strcat(current_max_lik_name, name_of_this_genotype);
	  //and do the likelihoods themselves
	  *current_max_but_one_lik = *current_max_lik;
	  *current_max_lik= log_prob_data + log_prob_error;
	}
      else if (log_prob_data + log_prob_error > *current_max_but_one_lik)
	{
	  //this is the new second best genotype
	  //sort out name
	  current_max_but_one_lik_name[0]='\0';
          strcat(current_max_but_one_lik_name, name_of_this_genotype);
	  //..and the likelihood
	  *current_max_but_one_lik = log_prob_data + log_prob_error;
	}

      return log_prob_data + log_prob_error;
    }
  
    
}


//we ASSUME colours 0 to number_alleles are the various alternate alleles, loading in multicolour_bin
void calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site(int* colours_to_genotype, int num_colours_to_genotype,
										      int colour_ref_minus_site, int number_alleles,
										      int first_gt, int last_gt, // of all the possible gt's
										      int max_allele_length, char* fasta,//one read per allele
										      AssumptionsOnGraphCleaning assump,
										      double* current_max_lik_array, double* current_max_but_one_lik_array,
										      char** name_current_max_lik_array, char** name_current_max_but_one_lik_array,
										      boolean print_all_liks_calculated,//not just the top two
										      GraphAndModelInfo* model_info, dBGraph* db_graph,
										      int working_colour1, int working_colour2
										      )
{

  printf("Start calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site\n");
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_allele_length,MAX_READ_NAME_LEN);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  
  
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_allele_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in db_graph_make_reference_path_based_sv_calls. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
  int j;
  int i;
  
  //we are going to need to hold ALL the allele paths in memory at the same time
  
  dBNode*** array_of_node_arrays = (dBNode***) malloc( sizeof(dBNode**) * number_alleles );
  Orientation** array_of_or_arrays = (Orientation**) malloc( sizeof(Orientation*) * number_alleles );
  int* lengths_of_alleles = (int*) malloc(sizeof(int) * number_alleles);
  char** array_of_allele_names = (char**) malloc( sizeof(char*) * number_alleles );
  
  if ( (array_of_node_arrays==NULL) || (array_of_or_arrays==NULL) || (lengths_of_alleles==NULL) || (array_of_allele_names==NULL) )
    {
      printf("Cannot alloc arrays of arrays in print_log_liks_of_specified_set_of_genotypes_of_complex_site\n");
      exit(1);
    }
  
  for (i=0; i<number_alleles; i++)
    {
      array_of_node_arrays[i] = (dBNode**) malloc(sizeof(dBNode*) * max_allele_length);
      array_of_or_arrays[i]   = (Orientation*) malloc( sizeof(Orientation) * max_allele_length);
      array_of_allele_names[i]= (char*) malloc(sizeof(char) * 200 );
      
      if ( (array_of_node_arrays[i]==NULL) || (array_of_or_arrays[i]==NULL) || (array_of_allele_names==NULL) )
	{
	  printf("Cannot alloc the %d -th node and or array in print_log_liks_of_specified_set_of_genotypes_of_complex_site", i);
	  exit(1);
	}
      lengths_of_alleles[i]=0;
      array_of_allele_names[i][0]='\0';
    }



  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv = alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(max_allele_length, max_allele_length);
  
  int* working_array_self = (int*) malloc(sizeof(int) * max_allele_length);
  int* working_array_shared = (int*) malloc(sizeof(int) * max_allele_length);
  
  if ( (mobv==NULL)||(working_array_self==NULL) || (working_array_shared==NULL))
    {
      printf("Cannot alloc all the arrays in calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site.  Give up and exit.");
      exit(1);
    }
  

  //this will give the indices WITHIN array_of_node_arrays of the alleles - there will be an offset compared with the 
  // colour numbers
  int list_alleles_with_covg[number_alleles];
  for (j=0; j<number_alleles; j++)
    {
      list_alleles_with_covg[j]=-1;
    }
  
  //end of initialisation 

  //create file reader
  int file_reader(FILE * fp, Sequence * seq, int max_allele_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      printf("new_entry must be true in hsi test function");
      exit(1);
    }
    ret =  read_sequence_from_fasta(fp,seq,max_allele_length,new_entry,full_entry,offset);
    
    return ret;
  }


  //Get list of alleles with coverage 
  FILE* fp = fopen(fasta, "r");
  if (fp==NULL)
    {
      printf("UNable to open %s. Exit", fasta);
      exit(1);
    }

  int index_alleles_with_nonzero_covg=0;

  for (j=0; j<= number_alleles-1; j++)
    {
      int num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_allele_length, 
								     array_of_node_arrays[j],//assumed alleles are colours 0..num_alleles-1 
								     array_of_or_arrays[j], 
								     false, file_reader, seq, kmer_window, db_graph, 0);
      
      strcat(array_of_allele_names[j], seq->name);
      lengths_of_alleles[j]=num_kmers;
      
      //Do any of these nodes have covg in the colours/samples we want to genotype
      int p;
      boolean there_is_no_covg=true; //so far, as far as we know, there is no covg in any of our samples
      for (p=0; (p< num_kmers) && (there_is_no_covg==true); p++)
	{
	  int q;
	  for (q=0; q<num_colours_to_genotype; q++)
	    {
	      if (db_node_get_coverage(array_of_node_arrays[j][p], individual_edge_array, colours_to_genotype[q])>0)
		{
		  there_is_no_covg=false;
		  break;
		}
	    }
	}
      if (there_is_no_covg==false)
	{
	  //add this to list of alleles for further consideration		
	  list_alleles_with_covg[index_alleles_with_nonzero_covg]=j;
	  index_alleles_with_nonzero_covg++;
	}
      
    }
  fclose(fp);

  int num_alleles_with_nonzero_covg = index_alleles_with_nonzero_covg;

  // For each pair in that list, keep incrementing a counter, and if the counter is betweenfirst_gt and last_gt
  // then go ahead and check it out

  int genotype_count=0;

  for (i=0; i<num_alleles_with_nonzero_covg; i++)
    {
      for (j=i; j<num_alleles_with_nonzero_covg; j++) 
	{
	  genotype_count++;
	  
	  if (genotype_count<first_gt)
	    {
	      continue;
	    }
	  if (genotype_count>last_gt)
	    {
	      continue;
	    }


	  VariantBranchesAndFlanks var;
	  set_variant_branches_but_flanks_to_null(&var, 
						  array_of_node_arrays[list_alleles_with_covg[i]], 
						  array_of_or_arrays[list_alleles_with_covg[i]], 
						  lengths_of_alleles[list_alleles_with_covg[i]],
						  array_of_node_arrays[list_alleles_with_covg[j]], 
						  array_of_or_arrays[list_alleles_with_covg[j]],
						  lengths_of_alleles[list_alleles_with_covg[j]],
						  unknown);

	  set_status_of_nodes_in_branches(&var, in_desired_genotype);
	  reset_MultiplicitiesAndOverlapsOfBiallelicVariant(mobv);

	  printf("Call initialise_multiplicities_of_allele_nodes_wrt_both_alleles\n");
	  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, mobv, false, NULL, NULL, working_colour1, working_colour2, db_graph);
	  printf("returns from initialise_multiplicities_of_allele_nodes_wrt_both_alleles\n");

	  int z;
	  for (z=0; z<num_colours_to_genotype; z++)
	    {
	      char name[400];
	      sprintf(name, "%s/%s", array_of_allele_names[i], array_of_allele_names[j]); 

	      printf("Call the actual calc\n");
	      double llk= calc_log_likelihood_of_genotype_with_complex_alleles(&var, name,
									       mobv, model_info, colours_to_genotype[z],
									       colour_ref_minus_site, db_graph, 
									       working_array_self, working_array_shared,
									       &(current_max_lik_array[z]), &(current_max_but_one_lik_array[z]),
									       name_current_max_lik_array[z], name_current_max_but_one_lik_array[z],
									       assump);
	      printf("Got the actual calc\n");


	      if (print_all_liks_calculated==true)
		{
		  printf("Colour %d, GENOTYPE %s : LLK=%f\n", colours_to_genotype[z], name, llk);
		}
	      set_status_of_nodes_in_branches(&var, none);
	      
	    }
	}
      



    }

  //finished this set
  printf("Finished calculating likelihoods of genotypes %d to %d\n", first_gt, last_gt);
  int z;
  for (z=0; z<num_colours_to_genotype; z++)
    {
      printf("Colour %d, MAX_LIKELIHOOD GENOTYPE %s : LLK=%f\n", colours_to_genotype[z], name_current_max_lik_array[z], current_max_lik_array[z]);
      printf("Colour %d, NEXT BEST GENOTYPE %s : LLK=%f\n", colours_to_genotype[z], name_current_max_but_one_lik_array[z], current_max_but_one_lik_array[z]);
    }

  dealloc_MultiplicitiesAndOverlapsOfBiallelicVariant(mobv);
  //many other things you should free here.

}

double* alloc_ML_results_array(int num_samples_to_genotype)
{
  double* retarray = (double*) malloc(sizeof(double) * num_samples_to_genotype);
  if (retarray==NULL)
    {
      printf("UNable to malloc %d colours in alloc_ML_results_array\n",num_samples_to_genotype);
      exit(1);
    }
  int i;
  for (i=0; i<num_samples_to_genotype; i++)
    {
      retarray[i]=MIN_LLK;
    }
  return retarray;
}

char** alloc_ML_results_names_array(int num_samples_to_genotype)
{
  char** retarray = (char**) malloc(sizeof(char*) * num_samples_to_genotype);
  if (retarray==NULL)
    {
      printf("Ridiculous - unanle to malloc array of names in alloc_ML_results_names_array");
      exit(1);
    }
  int i;
  for (i=0; i<num_samples_to_genotype; i++)
    {
      retarray[i] = (char*) malloc(sizeof(char) * MAX_READ_NAME_LEN);
      if (retarray[i]==NULL)
	{
	  printf("Unable to malloc name array in alloc_ML_results_names_array\n");
	  exit(1);
	}
      else
	{
	  retarray[i][0]='\0';
	}
    }
  return retarray;
}
