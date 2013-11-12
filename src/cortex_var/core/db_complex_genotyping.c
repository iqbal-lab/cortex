/*
 * Copyright 2009-2013 Zamin Iqbal and Mario Caccamo  
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
  db_complex_genotyping.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl_sf_gamma.h>

#include "db_variants.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "maths.h"
#include "db_complex_genotyping.h"
#include "model_selection.h"
#include "file_reader.h"
#include "genotyping_element.h"

/*
//assumes 3 colour graph, colour 0 is the union of all alleles, colour 1 is the ref-minus site. colour 2 will be which ever allele has just been loaded
//for each allele, create all 86000*3 possible SNPs and see what % of them touch the union/reference intersect union/allele.
// repeat for 2bp changes.
// once you have
void calculate_edit_probabilities(dBGraph* db_graph, char** array_allele_names, dBNode** array_alleles, int* array_allele_lens, FILE* fout, char* list_allele_binaries, int colour_union)
{
  //wipe colour 2

  //load j-th allele binary
  int i;

  for (i=0; i<array_allele_lens[j]; i+=50)
    {

  char orig_seq[db_graph->kmer_size+1];
  binary_kmer_to_seq(&(n->kmer), db_graph->kmer_size, orig_seq);

  for (i=0; i<db_graph->kmer_size; i++)
    {
      //there are 3 possible kmers created by changing position i 
      int p;
      for (p=0; p<=2; p++)
	{
	  char copy_orig[db_graph->kmer_size+1];
	  strcpy(copy_orig, orig_seq);
	  modify_character(copy_orig, i,p);//modify the i-th character and change it to the p-th of the other 3. eg if the base is C, and p=2 change it to the second one of A,G,T.
	  BinaryKmer tmp_kmer1, tmp_kmer2;
	  dBNode* found_node = hash_table_find(element_get_key(seq_to_binary_kmer(copy_orig, db_graph->kmer_size, &tmp_kmer1), db_graph->kmer_size, &tmp_kmer2), db_graph);
	  if (found_node==NULL)
	    {
	      //this particular edit of 1 base does not take you back to the desired genotype, but maybe some other edit does
	    }
	  else if (db_node_check_status(found_node, in_desired_genotype)==true)
	    {
	      return true;
	    }
	  else
	    {
	      //this particular edit of 1 base does not take you back to the desired genotype, but maybe some other edit does
	    }
	}

    }
  //all possible edits have failed to hit the genotype
  return false;




}

*/

// Lowest log likelihood, at which we just cutoff
int MIN_LLK = -9999999;
    
GenotypingWorkingPackage* alloc_genotyping_work_package(int max_allele_len,
                                                        int max_sup_len, 
                                                        int working_colour1,
                                                        int working_colour2)
{
  GenotypingWorkingPackage* gwp =(GenotypingWorkingPackage*)  malloc(sizeof(GenotypingWorkingPackage));
  if (gwp==NULL)
    {
      die("Unable to malloc GenotypingWorkingPackage 1");
    }
  gwp->max_allele_len = max_allele_len;
  gwp->max_sup_len    = max_sup_len;
  
  gwp->working_g_e_one      =(GenotypingElement**) malloc(sizeof(GenotypingElement*)*max_allele_len);
  gwp->working_o_one        =(Orientation*)        malloc(sizeof(Orientation)*max_allele_len);
  gwp->working_g_e_other    =(GenotypingElement**) malloc(sizeof(GenotypingElement*)*max_allele_len);
  gwp->working_o_other      =(Orientation*)        malloc(sizeof(Orientation)*max_allele_len);
  gwp->working_array_self   =(Covg*)               malloc(sizeof(Covg)*max_allele_len);
  gwp->working_array_shared =(Covg*)               malloc(sizeof(Covg)*max_allele_len);
  gwp->mobv                 =                      alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(max_allele_len, max_allele_len);
  gwp->working_colour1      =                      working_colour1;
  gwp->working_colour2      =                      working_colour2;
  gwp->path_nodes           =(dBNode**)            malloc(sizeof(dBNode*)*max_sup_len);
  gwp->path_orientations    =(Orientation*)        malloc(sizeof(Orientation) *max_sup_len);
  gwp->path_labels          =(Nucleotide*)         malloc(sizeof(Nucleotide) *max_sup_len);
  gwp->path_string          =(char*)               malloc(sizeof(char) *max_sup_len);

  if ( (gwp->working_g_e_one == NULL) ||
       (gwp->working_o_one  == NULL) ||
       (gwp->working_g_e_other == NULL) ||
       (gwp->working_o_other == NULL) ||
       (gwp->working_array_self == NULL) ||
       (gwp->working_array_shared == NULL) ||
       (gwp->path_nodes == NULL) ||
       (gwp->path_orientations == NULL) ||
       (gwp->path_labels== NULL) ||
       (gwp->path_string  == NULL) )
    {
      die("Unable to alloc GenotypingWorkingPackage members");
    }

  memset(gwp->working_array_self,   0, sizeof(Covg) * gwp->max_allele_len);
  memset(gwp->working_array_shared, 0, sizeof(Covg) * gwp->max_allele_len);
  
  return gwp;
}
 
 void free_genotyping_work_package(GenotypingWorkingPackage* gwp)
 {
   free(gwp->working_g_e_one);
   free(gwp->working_o_one);
   free(gwp->working_g_e_other);
   free(gwp->working_o_other);
   free(gwp->working_array_self);
   free(gwp->working_array_shared);
   free(gwp->path_nodes);
   free(gwp->path_orientations);
   free(gwp->path_labels);
   free(gwp->path_string);
   dealloc_MultiplicitiesAndOverlapsOfBiallelicVariant(gwp->mobv);
   free(gwp);
 }






 MultiplicitiesAndOverlapsOfBiallelicVariant*  alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(int len_allele1, int len_allele2)
{
  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv = (MultiplicitiesAndOverlapsOfBiallelicVariant*) malloc(sizeof(MultiplicitiesAndOverlapsOfBiallelicVariant));

  if (mobv==NULL)
    {
      die("Ridiculous - failed to alloc an array of 4 ints");
    }
  mobv->mult11=(Covg*) malloc( (sizeof(Covg))*len_allele1);
  mobv->mult12=(Covg*) malloc( (sizeof(Covg))*len_allele1);
  mobv->mult21=(Covg*) malloc( (sizeof(Covg))*len_allele2);
  mobv->mult22=(Covg*) malloc( (sizeof(Covg))*len_allele2);

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
void initialise_multiplicities_of_allele_nodes_wrt_both_alleles(
  VariantBranchesAndFlanks* var, MultiplicitiesAndOverlapsOfBiallelicVariant* mult,
  boolean only_count_nodes_with_edge_in_specified_colour_func,
  Edges (*get_colour)(const dBNode*))
  //, Covg (*get_covg)(const dBNode*) //unsued
{

  void get_mult(dBNode** br_src, int len_br_src, dBNode** br_target, int len_br_target, Covg* mult_array)
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
void improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(
  VariantBranchesAndFlanks* var, MultiplicitiesAndOverlapsOfBiallelicVariant* mult,
  int working_colour1, int working_colour2)
  // unused paramters
  //boolean only_count_nodes_with_edge_in_specified_colour_func,
  //Edges (*get_colour)(const dBNode*), Covg (*get_covg)(const dBNode*),
{
  // *** ASSUME THE WORKING COLOURS ARE CLEAN ***


  //walk through allele1, and as you do so, add multicplicity counts to working_colour1
  int i;
  for (i=0; i<var->len_one_allele; i++)
    {
      db_node_increment_coverage(var->one_allele[i], working_colour1);
    }
  //same for allele2
  for (i=0; i<var->len_other_allele; i++)
    {
      db_node_increment_coverage(var->other_allele[i], working_colour2);
    }
  //now, as you walk through allele1, for each node, you can see how many times it occurs in allele1 and in allele2
  for (i=0; i<var->len_one_allele; i++)
    {
      mult->mult11[i] += db_node_get_coverage(var->one_allele[i], working_colour1);
      mult->mult12[i] += db_node_get_coverage(var->one_allele[i], working_colour2);
    }

  for (i=0; i<var->len_other_allele; i++)
    {
      mult->mult21[i] += db_node_get_coverage(var->other_allele[i], working_colour1);
      mult->mult22[i] += db_node_get_coverage(var->other_allele[i], working_colour2);
    }

  //cleanup
  for (i=0; i<var->len_one_allele; i++)
    {
      db_node_set_coverage(var->one_allele[i], working_colour1, 0);
      db_node_set_coverage(var->one_allele[i], working_colour2, 0);
    }
  for (i=0; i<var->len_other_allele; i++)
    {
      db_node_set_coverage(var->other_allele[i], working_colour1,0);
      db_node_set_coverage(var->other_allele[i], working_colour2,0);
    }
  

}


void improved_initialise_multiplicities_of_allele_genotyping_nodes_wrt_both_alleles(
  GenotypingVariantBranchesAndFlanks* var, 
  MultiplicitiesAndOverlapsOfBiallelicVariant* mult,
  int working_colour1, int working_colour2)
{
  // *** ASSUME THE WORKING COLOURS ARE CLEAN ***


  //walk through allele1, and as you do so, add multicplicity counts to working_colour1
  int i;
  for (i=0; i<var->len_one_allele; i++)
    {
      db_genotyping_node_increment_coverage(var->one_allele[i], working_colour1);
    }
  //same for allele2
  for (i=0; i<var->len_other_allele; i++)
    {
      db_genotyping_node_increment_coverage(var->other_allele[i], working_colour2);
    }
  //now, as you walk through allele1, for each node, you can see how many times it occurs in allele1 and in allele2
  for (i=0; i<var->len_one_allele; i++)
    {
      mult->mult11[i] += db_genotyping_node_get_coverage(var->one_allele[i], working_colour1);
      mult->mult12[i] += db_genotyping_node_get_coverage(var->one_allele[i], working_colour2);
    }

  for (i=0; i<var->len_other_allele; i++)
    {
      mult->mult21[i] += db_genotyping_node_get_coverage(var->other_allele[i], working_colour1);
      mult->mult22[i] += db_genotyping_node_get_coverage(var->other_allele[i], working_colour2);
    }

  //cleanup
  for (i=0; i<var->len_one_allele; i++)
    {
      db_genotyping_node_set_coverage(var->one_allele[i], working_colour1, 0);
      db_genotyping_node_set_coverage(var->one_allele[i], working_colour2, 0);
    }
  for (i=0; i<var->len_other_allele; i++)
    {
      db_genotyping_node_set_coverage(var->other_allele[i], working_colour1,0);
      db_genotyping_node_set_coverage(var->other_allele[i], working_colour2,0);
    }
  



}





// WARNING - the implicit assumption here is that the graph/hash table ONLY
//           contains nodes in the union of all known alleles for 
//           the site in question.
// Pass in the current_max and current_max_but_one log likelihoods, and as soon
// as this drops below both, we abort mission.
// and return -99999999
// var_mults must be pre-allocated ***and pre-initialised***

double calc_log_likelihood_of_genotype_with_complex_alleles(
  VariantBranchesAndFlanks* var,
  char* name_of_this_genotype,
  MultiplicitiesAndOverlapsOfBiallelicVariant* var_mults,
  GraphAndModelInfo* model_info,
  int colour_indiv, 
  int colour_ref_minus_our_site, dBGraph* db_graph,
  Covg* working_array_self, Covg* working_array_shared,
  double* current_max_lik, double* current_max_but_one_lik,
  char* current_max_lik_name, char* current_max_but_one_lik_name,
  AssumptionsOnGraphCleaning assump,
  dBNode** p_nodes, Orientation* p_orientations,
  Nucleotide* p_labels, char* p_string, int max_allele_length,
  boolean using_1net, Covg (*get_covg_in_1net_of_genotype)(dBNode*), 
  boolean using_2net, Covg (*get_covg_in_2net_of_genotype)(dBNode*),
  double min_acceptable_llk)
//assume the  2 working arrays have length = max read length
{
    
  // *********************************************************************************************************************************************************
  // Likelihood of this genotype is given by prob(see number_errors of errors in sequence with length=sum of lengths of two alleles)
  //                                         * Product ( dpois(observed covg, expected covg) )     (product over subchunks, splitting where alleles overlap)
  // *********************************************************************************************************************************************************

  if ( (using_1net==false) && (using_2net==true))
    {
      die("you cannot use the 2net woithout using the 1net - code error");
    }

  Edges element_get_colour_indiv(const Element* e)
  {
    return e->individual_edges[colour_indiv];
  }

  Covg element_get_covg_indiv(const dBNode* e)
  {
    return db_node_get_coverage(e, colour_indiv);
  }


  // 1. Count number of errors:

  Covg check_covg_in_ref_with_site_excised(dBNode* n)
  {
    if (colour_ref_minus_our_site==-1)
      {
	return 0;
      }
    else
      {
	return db_node_get_coverage(n, colour_ref_minus_our_site);
      }
  }
  
  void get_number_nodes_outside_our_expected_two_alleles_with_covg( dBNode* e, Covg* total) 
  { 
    if (db_node_check_status(e, in_desired_genotype)==true) 
      {
      } 
    else if ( (db_node_get_coverage(e,colour_indiv)>0) && (check_covg_in_ref_with_site_excised(e)==0) ) 
      {
	//*total = (*total) + 1;
	//char zam[db_graph->kmer_size +1];
	//printf("Error node is %s, with covg %d\n", binary_kmer_to_seq(&(e->kmer), db_graph->kmer_size, zam), db_node_get_coverage(e, colour_indiv) );
	(*total) = (*total) + db_node_get_coverage(e,colour_indiv);
      }
  }
  

   // the three ints you return are basically levels of badness of error. 
  void count_reads_in_1net_2net_and_beyond(dBNode* e, int* total_1net, int* total_2net, int* total_3net,
					   dBNode** p_nodes, Orientation* p_or, Nucleotide* p_lab, char* p_str, int max_len)
  {

    if ( (db_node_check_status(e, in_desired_genotype)==true) || (db_node_check_status(e, visited)==true) )
      {
      }
    else if ( (db_node_get_coverage(e,colour_indiv)>0) && (check_covg_in_ref_with_site_excised(e)==0) )
      {
	double avg_coverage=0;
	Covg min = 0, max=0;
	boolean is_cycle=false;
	//get the whole supernode. If any of them are in_desired, then is 1bp away, else not
	int length = db_graph_supernode_for_specific_person_or_pop(e,max_len,&db_node_action_set_status_visited,
								   p_nodes,p_or,p_lab, p_str,
								   &avg_coverage,&min,&max,&is_cycle,
								   db_graph, colour_indiv);
	int i;
	typedef enum
	{
	  worst_is_1net=0,
	  worst_is_2net=1,
	  worst_is_beyond_2net=2,
	} WorstNodeInSup;

	WorstNodeInSup wor =worst_is_1net;
	for (i=1; (i<length) && (using_1net==true); i++)//if not using_1net, then don't bother - assume everything in this supernode is bad, there is only 1 category of bad, 
	  {
	    if (get_covg_in_1net_of_genotype(p_nodes[i])>0) 
	      {
	      }
	    else if ( (using_2net==true) && (get_covg_in_2net_of_genotype(p_nodes[i])>0) )
	      {
		if (wor !=worst_is_beyond_2net )
		  {
		    wor=worst_is_2net;
		  }
	      }
	    else if (using_2net==true)
	      {
		wor=worst_is_beyond_2net;
		i=length+1;
	      }
	    else 
	      {
		wor=worst_is_2net;
	      }
	  }

	boolean too_short=false;
	Covg ct = count_reads_on_allele_in_specific_colour(p_nodes, length, colour_indiv, &too_short);
	int extra=0;
	if (too_short==false)
	  {
	    if (length>db_graph->kmer_size)
	      {
		extra=length-db_graph->kmer_size;
	      }//tells us if >1 SNP
	    if (wor==worst_is_1net)
	      {
		/*
		  if (too_short==false)
		  {
		  *total_1net = (*total_1net) + ct;
		  }
		  else
		  {
		  *total_1net = (*total_1net) + db_node_get_coverage(e,colour_indiv);
		  }
		*/
		//*total_1net = (*total_1net) + 1;
		//*total_1net = (*total_1net) + db_node_get_coverage(e,colour_indiv) ;
		*total_1net = (*total_1net) + ct + extra+1;
		//*total_1net = (*total_1net) + extra+1;
		//*total_1net = (*total_1net) + length;
	      }
	    else if (wor==worst_is_2net)
	      {
		/*
		  if (too_short==false)
		  {
		  *total_2net = (*total_2net) + ct;	
		  }
		  else
		  {
		  *total_2net = (*total_2net) + db_node_get_coverage(e,colour_indiv) * (extra + 1);
		  }
		*/
		*total_2net = (*total_2net) + ct+extra+1;
		//*total_2net = (*total_2net) + 1;
		//*total_2net = (*total_2net) + db_node_get_coverage(e,colour_indiv);
		//*total_2net = (*total_2net) + length;
	      }
	    else
	      {
		/*
		  if (too_short==false)
		  {
		  *total_3net = (*total_3net) + ct;	
		  }
		  else
		  {
		  *total_3net = (*total_3net) + db_node_get_coverage(e,colour_indiv) * (extra + 1);
		  }
		*/
		*total_3net = (*total_3net) + ct+extra+1;
		//*total_3net = (*total_3net) + 1;
		//*total_3net = (*total_3net) + db_node_get_coverage(e,colour_indiv);
		//*total_3net = (*total_3net) + length;
	      }
	  }
      }
    
    return;
  }

  /*
  void get_number_nodes_outside_our_expected_two_alleles_with_covg_and_outside_1net( dBNode* e, int* total_1net, int* total_beyond_1net) 
  {
    if (db_node_check_status(e, in_desired_genotype)==true) 
      {
      }
    else if (get_covg_in_colours_for_errors_from_genotype(e)>0)
      {
	*total_1net = (*total_1net) +1;
      }
    else
      {
	*total_beyond_1net = (*total_beyond_1net)+1;
      }
  }
  */

  void count_errors_1bp_and_further_from_desired_genotype( dBNode* e, int* total_1bp_away, int* total_further,
							   dBNode** p_nodes, Orientation* p_or, Nucleotide* p_lab, char* p_str, int max_len)
  {
    if ( (db_node_check_status(e, in_desired_genotype)==true) || (db_node_check_status(e, visited)==true) )
      {
      }
    else if ( (db_node_get_coverage(e,colour_indiv)>0) && (check_covg_in_ref_with_site_excised(e)==0) )
      {
	double avg_coverage=0;
	Covg min = 0, max=0;
	boolean is_cycle=false;
	//get the whole supernode. If any of them are in_desired, then is 1bp away, else not
	int length = db_graph_supernode_for_specific_person_or_pop(e,max_len,&db_node_action_set_status_visited,
								   p_nodes,p_or,p_lab, p_str,
								   &avg_coverage,&min,&max,&is_cycle,
								   db_graph, colour_indiv);
	int i;
	boolean at_least_one_node_in_sup_is_in_desired=false;
	for (i=0; i<length; i++)
	  {
	    if (db_node_check_status(p_nodes[i], in_desired_genotype)==true)
	      {
		at_least_one_node_in_sup_is_in_desired=true;
		i=length+1;
	      }
	  }
	if (at_least_one_node_in_sup_is_in_desired==true)
	  {
	    /*//the following seems like  a good idea, but it fails my unit tests on carefuly checked test cases
	      // should examine why
	    boolean too_short=false;
	    Covg ct = count_reads_on_allele_in_specific_colour(p_nodes, length, colour_indiv, &too_short);
	    if (too_short==false)
	      {
		*total_1bp_away = (*total_1bp_away)+ct;
	      }
	    */
	
	     *total_1bp_away = (*total_1bp_away)+1;
	  }
	else
	  {
	    /*
	    //the following seems like  a good idea, but it fails my unit tests on carefuly checked test cases
	      // should examine why

	    boolean too_short=false;
	    Covg ct = count_reads_on_allele_in_specific_colour(p_nodes, length, colour_indiv, &too_short);
	    if (too_short==false)
	      {
		*total_further = (*total_further)+ct;
	      }
	    */

	     *total_further = (*total_further)+1;
	  }
      }
    
    return;
  }


  
  boolean condition_is_error_supernode(dBNode** path, int len,  int* count_error_nodes, Covg* count_error_covg)
  {
    boolean sup_is_error=false;
    if (len>0)
      {
	int i;
	
	*count_error_nodes=0;
	*count_error_covg=0;
	Covg total_covg=0;
	for (i=0; i<len; i++)
	  {
	    total_covg += db_node_get_coverage(path[i], colour_indiv);
	    
	    if (db_node_check_status(path[i],in_desired_genotype)==false) 
	      {		
		if (check_covg_in_ref_with_site_excised(path[i])==0)
		  {
		    //this node is not in our genotype AND not in rest of genome ===> error
		    sup_is_error=true;
		    (*count_error_nodes) = (*count_error_nodes) +1;
		    *count_error_covg = *count_error_covg + db_node_get_coverage(path[i], colour_indiv);//wont' really use this, but interesting - what you really want is number of error reads
		    
		  }
		
	      }
	  }
      }
    else if ( (db_node_check_status(path[0],in_desired_genotype)==false) && (check_covg_in_ref_with_site_excised(path[0])==0) )
      {
	(*count_error_nodes) = (*count_error_nodes) +1;
	*count_error_covg = *count_error_covg + db_node_get_coverage(path[0], colour_indiv);
	sup_is_error=true;
      }
    return sup_is_error;
  }
  
  
  
  int num_1net_errors=0;
  int num_2net_errors=0;
  int num_3net_errors=0;//everything worse than2net is called 3net
  
  //  int number_errors=0;//get number of nodes that are (in the union of all known alleles but) outside our genotype alleles, and are not in the rest of the genome, and are a SNP away from the genotype
  int total_len_alleles = (var->len_one_allele + var->len_other_allele);
  //removed next line - trying better error counting
  //hash_table_traverse_passing_int(&get_number_nodes_outside_our_expected_two_alleles_with_covg  ,db_graph, &number_errors);


  //new option - look at errors more closely
  //int number_bad_errors=0; //errors > 1 SNP from genotype

  //  hash_table_traverse_passing_ints_and_path(&count_errors_1bp_and_further_from_desired_genotype  ,db_graph, &number_errors, &number_bad_errors,
  //					    p_nodes, p_orientations, p_labels, p_string, max_allele_length);
  // hash_table_traverse_passing_ints_and_path(&count_errors_in_1net_of_desired_genotype_or_worse  ,db_graph, &number_errors, &number_bad_errors,
  //					    p_nodes, p_orientations, p_labels, p_string, max_allele_length);

  hash_table_traverse_passing_3ints_and_path(&count_reads_in_1net_2net_and_beyond, db_graph, &num_1net_errors, &num_2net_errors, &num_3net_errors,
					     p_nodes, p_orientations, p_labels, p_string, max_allele_length);
  hash_table_traverse(&db_node_action_set_status_none, db_graph); 
  set_status_of_nodes_in_branches(var, in_desired_genotype);
  
  //number_errors=number_errors/db_graph->kmer_size;
  
  /*
    long long better_number_errors = db_graph_count_error_supernodes(max_allele_length, db_graph, colour_indiv, 
								   p_nodes, p_orientations, p_labels, p_string,
								   &db_node_check_status_is_not_visited, &condition_is_error_supernode,&db_node_action_set_status_visited);
  number_errors=(int)better_number_errors;
  
  //cleanup after that traversal << sorry, you cannot remove the following two lines to improve performance - screws things up if >1 sample
  hash_table_traverse(&db_node_action_set_status_none, db_graph); 
  set_status_of_nodes_in_branches(var, in_desired_genotype);
  */

  printf("Genotype %s - number errors in 1net %d, 2net  %d and 3net %d\n",
         name_of_this_genotype, num_1net_errors, num_2net_errors, num_3net_errors);

  // Errors on union of two alleles occur as a Poisson process with rate = sequencing_error_rate_per_base * kmer * length of two alleles
  //double poisson_error_rate;
  //double poisson_bad_error_rate;
  double poisson_1net_err_rate=0;
  double poisson_2net_err_rate=0;
  double poisson_3net_err_rate=0;

  //double m=(model_info->seq_error_rate_per_base) ;
  double m=model_info->ginfo->seq_err[colour_indiv];
  //int kmer=(db_graph->kmer_size);
  //double prop = max_allele_length/(model_info->genome_len);
  if (assump==AssumeUncleaned)
    {
      //commented out
      // poisson_error_rate = m * kmer *total_len_alleles;
      //poisson_bad_error_rate = m*m*m*kmer*kmer*kmer*total_len_alleles; 

      poisson_1net_err_rate= m*total_len_alleles;
      poisson_2net_err_rate= m*m*total_len_alleles;
      poisson_3net_err_rate= m*m*m*total_len_alleles;

    }
  //  else if (assump==AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice)
  else
    {

      //next line is the one I had
      //poisson_error_rate = m*m* kmer*kmer* total_len_alleles;
      
      poisson_1net_err_rate= m*m*total_len_alleles;
      poisson_2net_err_rate= m*m*m*m*total_len_alleles;
      poisson_3net_err_rate= m*m*m*m*m*m*total_len_alleles;
    }

  double log_prob_error = 0;
  //printf("Poisson error rate is %.20f\n", poisson_error_rate);

  /*
  //log of (exp(-poisson_error_rate) (poisson_error_rate^number_errors)/(number_errors!) ) = -poisson_error_rate + number_errors*log(poisson_error_rate) - log(number_errors!)
  if (number_errors>0)
    {
      log_prob_error += -poisson_error_rate + number_errors*log(poisson_error_rate) - gsl_sf_lnfact(number_errors);
    }
  if (number_bad_errors>0)
    {
      //log_prob_error = -poisson_error_rate*poisson_error_rate + number_bad_errors*log(poisson_error_rate*poisson_error_rate) - gsl_sf_lnfact(number_bad_errors);
      log_prob_error += -poisson_bad_error_rate + number_bad_errors*log(poisson_bad_error_rate) - gsl_sf_lnfact(number_bad_errors);
    }

  if (log_prob_error < *current_max_but_one_lik)
    {
      return MIN_LLK;//abort - this genotype cannot be the max likelihood genotype
    }
  */


  if (num_1net_errors>0)
    {
      log_prob_error += -poisson_1net_err_rate + num_1net_errors*log(poisson_1net_err_rate) - gsl_sf_lnfact(num_1net_errors);
    }
  if (num_2net_errors>0)
    {
      log_prob_error += -poisson_2net_err_rate + num_2net_errors*log(poisson_2net_err_rate) - gsl_sf_lnfact(num_2net_errors);
    }
  if (num_3net_errors>0)
    {
      log_prob_error += -poisson_3net_err_rate + num_3net_errors*log(poisson_3net_err_rate) - gsl_sf_lnfact(num_3net_errors);
    }
  



  // 2.  Get sum of log likelihoods of all the subchunks of the two alleles

  // walk along allele 1. For sections unique to allele1, work out log likelihood with eff haploid covg
  // for sections shared with allele2, use eff diploid covg
  // and finally, walk along allele2, and consider only sections unique to allele2


  //int working_array_self_count=0;
  //int working_array_shared_count=0;
  double log_prob_data=0;
  double hap_D  = (double) (model_info->ginfo->total_sequence[colour_indiv])/(double) (2*model_info->genome_len) ;
  double hap_D_over_R;
  if (model_info->ginfo->mean_read_length[colour_indiv]>0)
    {
      hap_D_over_R = hap_D/(model_info->ginfo->mean_read_length[colour_indiv]);
    }
  else
    {
      hap_D_over_R=0;
    }
  //int eff_r_plus_one = model_info->ginfo->mean_read_length[colour_indiv] - db_graph->kmer_size;

  //printf("total sequence is %qd and genome len is %qd , and hapD is %f and hapD over R is %f\n", model_info->ginfo->total_sequence[colour_indiv], model_info->genome_len, hap_D, hap_D_over_R);

  //Option 1
  //break allele1 into intervals that only lie on allele1, and those shared with allele2.
  //traverse allele 1 in chunks. If a chunk exists only in allele 1, then divide covg by multiplicity in 1, and calc prob
  //                             If a chunk exists in both alleles, divide covg by multiplicty in 1 plus mult in 2, and calc prob.
  //then do allele 2 in identical manner

  //Option 2 - Ignore for now
  //traverse allele 1 abd then allele 2. For a given chunk, take full covg, but use combined mult in alleles 1 and 2 to determne poisson mean, and then mark chunk visited and do not return to visited chunks.

  //Option 3 - Ignore for now  - dont use chunks at all, use supernodes <<< will come with sup hash


  //Implement Option1  below:


  //walk the first allele
  log_prob_data
    += get_log_probability_of_covg_on_one_allele_given_second_allele_and_multiplicities(
         hap_D_over_R, var->one_allele, var->len_one_allele,
         var_mults->mult11, var_mults->mult12,
         working_array_self, working_array_shared,
         check_covg_in_ref_with_site_excised, colour_indiv);

  if(log_prob_error + log_prob_data < *current_max_but_one_lik ||
     log_prob_error + log_prob_data < min_acceptable_llk)
    {
      return MIN_LLK;
    }

  //now walk the second (which will be the same one again if this is hom)
  log_prob_data
    += get_log_probability_of_covg_on_one_allele_given_second_allele_and_multiplicities(
         hap_D_over_R, var->other_allele, var->len_other_allele,
         var_mults->mult22, var_mults->mult21,
         working_array_self, working_array_shared,
         check_covg_in_ref_with_site_excised, colour_indiv);
  
  if(log_prob_error + log_prob_data < *current_max_but_one_lik ||
     log_prob_error + log_prob_data < min_acceptable_llk)
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


//eff_D_over_R is either hap_D_over_R or 2*hap_D_over_R
double calc_log_prob_of_covg_on_chunk(double eff_D_over_R, Covg* working_array,
                                      int working_array_len)
{
  boolean too_short = false;

  Covg covg_in_this_chunk
    = count_reads_on_allele_in_specific_colour_given_array_of_cvgs(working_array,
                                                                   working_array_len,
                                                                   &too_short);
  
  double log_prob=0;
  if (too_short==false)
  {
    // add log dpois ( eff_D_over_R * working_array_self_count , covg_on_self_in_this_chunk)
    log_prob += -eff_D_over_R * working_array_len +
                (double)covg_in_this_chunk * log(eff_D_over_R * working_array_len) -
                gsl_sf_lnfact(covg_in_this_chunk);
  }

  return log_prob;
}


//traverse an allele, with a second allele in mind, given multiplicities of each node in either allele
// break allele into unique and shared chunks, and (given haplid covg), return the log probabiity
// of seeing this covg on this allele, given that we know the truth consists of this allele AND the other one
// thus if our allele consists of these chunks ABCABAB and the other allele consists of DEFAGBB
// then we are going to account for 3/4 of the covg on A adnd 3/5 of the covg on B while traversing this allele, and all of C
// in the knowldge that we will come back and traverse the second allele later and account for the remainin 1/4 of A covg and 3/5 of B covg, plus D,E,F,G.

double get_log_probability_of_covg_on_one_allele_given_second_allele_and_multiplicities(
  double hap_D_over_R, dBNode** allele, int len_allele,
  Covg* mult_this_allele_in_self, Covg* mult_this_allele_in_other,
  Covg* working_array_self, Covg* working_array_shared,
  Covg (*check_covg_in_ref_with_site_excised)(dBNode*), int colour_indiv)
{
  
  double log_prob_data=0;
  int k=0;
  int working_array_self_count  =0;
  int working_array_shared_count=0;
  while (k < len_allele)
  {
    //this node occurs  >0 times in the other allele
    if (mult_this_allele_in_other[k]>0)
    {
      //contribution from previous chunk, which was unique to this allele
      if (working_array_self_count>0)
      {
        log_prob_data += calc_log_prob_of_covg_on_chunk(hap_D_over_R,
                                                        working_array_self,
                                                        working_array_self_count);
        working_array_self_count = 0;
      }
    }

    // automatically handles case when there is no colour for ref-minus-site
    while(k < len_allele && mult_this_allele_in_other[k] > 0 && 
          check_covg_in_ref_with_site_excised(allele[k]) == 0)
    {
      // this node might happen 2 times on allele 1 and 3 times on allele2. in
      // total 5 times, and here it counts for 2 of them
      working_array_shared[working_array_shared_count]
        = db_node_get_coverage(allele[k], colour_indiv) /
          (mult_this_allele_in_self[k] + mult_this_allele_in_other[k]);

      //printf("Shared segment: Dividing covg by %d\n",
      //       mult_this_allele_in_self[k] + mult_this_allele_in_other[k]);

      k++;
      working_array_shared_count++;
    }

    while(k < len_allele && check_covg_in_ref_with_site_excised(allele[k]) > 0)
    {
      k++;
    }
	
	
	if (mult_this_allele_in_other[k]==0)
	  {
	    //start of contiguous chunk which occurs only in first allele
	    if (working_array_shared_count>0)
	      {
		//have broken up allle and dividec covg by multiplicity in both alleles, so just hap_D below
		log_prob_data += calc_log_prob_of_covg_on_chunk(hap_D_over_R, working_array_shared, working_array_shared_count);
		working_array_shared_count=0;
	      }
	  }
	
	    
	while ((k < len_allele) && 
	       (mult_this_allele_in_other[k]==0) && 
	       (check_covg_in_ref_with_site_excised(allele[k])==0)
	       )
	  {
	    if (mult_this_allele_in_self[k]>0)
	      {
		working_array_self[working_array_self_count]=db_node_get_coverage(allele[k], colour_indiv)/mult_this_allele_in_self[k];
		working_array_self_count++;
	      }
	    //printf("Self segment: Dividing covg by %d, k is %d and len allele is %d\n", mult_this_allele_in_self[k], k, len_allele);
	    k++;


	  }

	while ( (k < len_allele) &&( check_covg_in_ref_with_site_excised(allele[k])>0 ) )
	  { 
	    k++;
	  }
	
  }

  //may have exited loop with some left-over data in arrays, so:  
  if (working_array_self_count>0)
    {
      log_prob_data += calc_log_prob_of_covg_on_chunk(hap_D_over_R, working_array_self, working_array_self_count);
      working_array_self_count=0;
    }
  
  if (working_array_shared_count>0)
    {
      log_prob_data += calc_log_prob_of_covg_on_chunk(2*hap_D_over_R, working_array_shared, working_array_shared_count);
      working_array_shared_count=0;
    }


  return log_prob_data;

}




//traverse an allele, with a second allele in mind, given multiplicities of each node in either allele
// break allele into unique and shared chunks, and (given haplid covg), return the log probabiity
// of seeing this covg on this allele, given that we know the truth consists of this allele AND the other one
// thus if our allele consists of these chunks ABCABAB and the other allele consists of DEFAGBB
// then we are going to account for 3/4 of the covg on A adnd 3/5 of the covg on B while traversing this allele, and all of C
// in the knowldge that we will come back and traverse the second allele later and account for the remainin 1/4 of A covg and 3/5 of B covg, plus D,E,F,G.

double get_log_probability_of_covg_on_one_allele_given_second_allele_and_multiplicities_using_little_hash(
  double hap_D_over_R, GenotypingElement** allele, int len_allele,
  Covg* mult_this_allele_in_self, Covg* mult_this_allele_in_other,
  Covg* working_array_self, Covg* working_array_shared,
  Covg (*check_covg_in_ref_with_site_excised)(GenotypingElement*),
  int colour_indiv)
{
  double log_prob_data=0;
  int k=0;
  int working_array_self_count  =0;
  int working_array_shared_count=0;
  while (k < len_allele)
      {

	if (mult_this_allele_in_other[k]>0) //this node occurs  >0 times in the other allele
	  {
	    if (working_array_self_count>0) //contribution from previous chunk, which was unique to this allele
	      {
		log_prob_data += calc_log_prob_of_covg_on_chunk(hap_D_over_R, working_array_self, working_array_self_count);
		working_array_self_count=0;
	      }
	  }

	while ((k < len_allele) &&  
	       (mult_this_allele_in_other[k] > 0) && 
	       (check_covg_in_ref_with_site_excised(allele[k])==0 )//automatically handles case when there is no colour for ref-minus-site     
	       )
	  {
	    // this node might happen 2 times on allele 1 and 3 times on allele2. in total 5 times, annd here it counts for 2 of them
	    working_array_shared[working_array_shared_count]= 
	      db_genotyping_node_get_coverage(allele[k], colour_indiv)/(mult_this_allele_in_self[k] + mult_this_allele_in_other[k])  ;
	    //printf("Shared segment: Dividing covg by %d\n", mult_this_allele_in_self[k] + mult_this_allele_in_other[k]);
	    k++;
	    working_array_shared_count++;

	  }

	while ( (k<len_allele) && (check_covg_in_ref_with_site_excised(allele[k])>0 ) )
	  {
	    k++;
	  }
	
	
	if (mult_this_allele_in_other[k]==0)
	  {
	    //start of contiguous chunk which occurs only in first allele
	    if (working_array_shared_count>0)
	      {
		//have broken up allle and dividec covg by multiplicity in both alleles, so just hap_D below
		log_prob_data += calc_log_prob_of_covg_on_chunk(hap_D_over_R, working_array_shared, working_array_shared_count);
		working_array_shared_count=0;
	      }
	  }
	
	    
	while ((k < len_allele) && 
	       (mult_this_allele_in_other[k]==0) && 
	       (check_covg_in_ref_with_site_excised(allele[k])==0)
	       )
	  {
	    if (mult_this_allele_in_self[k]>0)
	      {//in PD calls, can have N in the ref allele
		working_array_self[working_array_self_count]=db_genotyping_node_get_coverage(allele[k], colour_indiv)/mult_this_allele_in_self[k];
		working_array_self_count++;
	      }
	    //printf("Self segment: Dividing covg by %d, k is %d and len allele is %d\n", mult_this_allele_in_self[k], k, len_allele);
	    k++;


	  }

	while ( (k < len_allele) &&( check_covg_in_ref_with_site_excised(allele[k])>0 ) )
	  { 
	    k++;
	  }
	
      }

  //may have exited loop with some left-over data in arrays, so:  
  if (working_array_self_count>0)
    {
      log_prob_data += calc_log_prob_of_covg_on_chunk(hap_D_over_R, working_array_self, working_array_self_count);
      working_array_self_count=0;
    }
  
  if (working_array_shared_count>0)
    {
      log_prob_data += calc_log_prob_of_covg_on_chunk(2*hap_D_over_R, working_array_shared, working_array_shared_count);
      working_array_shared_count=0;
    }


  return log_prob_data;

}









//for now, does not support 1net or 2net
//caller ensures that in the LittleHash the colour specified as model_info->ref_colour actually has covg in ref-minus-site.

double calc_log_likelihood_of_genotype_with_complex_alleles_using_little_hash(GenotypingVariantBranchesAndFlanks* var,
									      MultiplicitiesAndOverlapsOfBiallelicVariant* var_mults,
									      GraphAndModelInfo* model_info,
									      int colour_indiv, 
									      LittleHashTable* little_db_graph, dBGraph* db_graph,
									      Covg* working_array_self, Covg* working_array_shared,
									      AssumptionsOnGraphCleaning assump,
									      dBNode** p_nodes, Orientation* p_orientations, 
									      Nucleotide* p_labels, char* p_string, 
									      int max_sup_len
									      )
//assume the  2 working arrays have length = max read length
{
    
  // *********************************************************************************************************************************************************
  // Likelihood of this genotype is given by prob(see number_errors of errors in sequence with length=sum of lengths of two alleles)
  //                                         * Product ( dpois(observed covg, expected covg) )     (product over subchunks, splitting where alleles overlap)
  // *********************************************************************************************************************************************************

  // 1. Count number of errors:
  int colour_ref_minus_site = model_info->ref_colour; //caller of function has ensured this now contains covg in ref-minus-site
  Covg check_covg_in_ref_with_site_excised(GenotypingElement* n)
  {
    if (colour_ref_minus_site==-1)
      {
	return 0;
      }
    else
      {
	return db_genotyping_node_get_coverage(n, colour_ref_minus_site);
      }
  }
  
  

  // the three ints you return are basically levels of badness of error. 
  void count_reads_in_1net(GenotypingElement* e, int* total_1net, int* total_2net, int* total_3net,
			   dBNode** p_nod, Orientation* p_or, Nucleotide* p_lab, char* p_str, int max_len)
  {
    //zam - fixed tis in debugging if ( (db_genotyping_node_check_status(e, in_desired_genotype)==true) || (db_genotyping_node_check_status(e, visited)==true) )
    if (db_genotyping_node_check_status(e, in_desired_genotype)==true)
      {
	//	printf("Zam2 - traversing little graph - ignore this is visited/in genotype %s\n", binary_kmer_to_seq(&(e->kmer), db_graph->kmer_size, zamzam));
      }
    else if (db_genotyping_node_get_coverage(e,colour_indiv) > 0 &&
             check_covg_in_ref_with_site_excised(e) == 0)
      {
	double avg_coverage=0;
	Covg min = 0, max = 0;
	boolean is_cycle=false;
	//get the whole supernode IN THE MAIN GRAPH. If any of them are in_desired, then is 1bp away, else not
	dBNode* node_corresponding_to_e = hash_table_find(&(e->kmer), db_graph);
	if (node_corresponding_to_e==NULL)
	  {
	    char tmp[db_graph->kmer_size+1];
	    die("During counting, cannot find node %s (which IS in the little graph) "
          "in the main graph. Abort.",
          binary_kmer_to_seq(&(e->kmer), db_graph->kmer_size, tmp));
	  }
	else if (db_node_check_status_special(node_corresponding_to_e)==true)
	  {
	    //already done this
	    return;
	  }
	int length = db_graph_supernode_for_specific_person_or_pop(node_corresponding_to_e,max_len,
								   &db_node_action_specialise_status,
								   p_nod,p_or,p_lab, p_str,
								   &avg_coverage,&min,&max,&is_cycle,
								   db_graph, colour_indiv);
	//assume everything in supernode is bad. Not using 1net2net etc

	boolean too_short=false;
	Covg ct = count_reads_on_allele_in_specific_colour(p_nod, length, colour_indiv, &too_short);
	if (too_short==false)
	  {
	    int extra=0;
	    if (length>db_graph->kmer_size)
	      {
		extra=length-db_graph->kmer_size;
	      }//tells us if >1 SNP
	    
	    *total_1net = (*total_1net) + ct + extra+1;
	  }
      }
    
    return;
  }


  void reset_nodes_in_main_graph_to_visited(GenotypingElement* e, dBNode** p_nod, Orientation* p_or, Nucleotide* p_lab, char* p_str, int max_len)
  {

    if (db_genotyping_node_check_status(e, in_desired_genotype)==true)
      {
      }
    else if ( (db_genotyping_node_get_coverage(e,colour_indiv)>0) && (check_covg_in_ref_with_site_excised(e)==0) )
      {
	double avg_coverage=0;
	Covg min = 0, max = 0;
	boolean is_cycle=false;
	//get the corresponding node IN THE MAIN GRAPH. 
	dBNode* node_corresponding_to_e = hash_table_find(&(e->kmer), db_graph);
	if (node_corresponding_to_e==NULL)
	  {
	    char tmp[db_graph->kmer_size+1];
	    die("During genotyping reset, cannot find node %s "
          "(which IS in the little graph) in the main graph. Abort.",
          binary_kmer_to_seq(&(e->kmer), db_graph->kmer_size, tmp));
	  }
	else 
	  {
	    db_graph_supernode_for_specific_person_or_pop(node_corresponding_to_e,max_len,
							  &db_node_action_unspecialise_status,
							  p_nod,p_or,p_lab, p_str,
							  &avg_coverage,&min,&max,&is_cycle,
							  db_graph, colour_indiv);
	    
	    
	  }
	
      }
    
    return;
  }


  
  int num_1net_errors=0;
  //not using the next two for now - preserving due to api, and because may well start using them
  int num_2net_errors=0;
  int num_3net_errors=0;
  
  int total_len_alleles = (var->len_one_allele + var->len_other_allele);

  //new option - look at errors more closely
  // int number_bad_errors=0; //errors > 1 SNP from genotype

  //printf("Zam start traverse to count\n");
  little_hash_table_traverse_passing_3ints_and_big_graph_path(&count_reads_in_1net, little_db_graph, 
							      &num_1net_errors, 
							      &num_2net_errors, 
							      &num_3net_errors,
							      p_nodes, p_orientations, 
							      p_labels, p_string, max_sup_len);

  //printf("Zam end traverse to count\n");
  //printf("Zam start traverse to reset\n");
  //now go back one more time and un-set the visited flags you put on all those nodes in the main graph
  little_hash_table_traverse_passing_big_graph_path(&reset_nodes_in_main_graph_to_visited, little_db_graph, p_nodes, p_orientations, p_labels, p_string, max_sup_len);
  //printf("Zam end reset\n");

  // Errors on union of two alleles occur as a Poisson process with rate = sequencing_error_rate_per_base * kmer * length of two alleles

  double poisson_1net_err_rate=0;
  //double poisson_2net_err_rate=0;
  //double poisson_3net_err_rate=0;

  //double m=(model_info->seq_error_rate_per_base) ;
  double m=model_info->ginfo->seq_err[colour_indiv];
  //int kmer=(db_graph->kmer_size);
  //double prop = max_allele_length/(model_info->genome_len);

  if (assump==AssumeUncleaned)
    {
      poisson_1net_err_rate= m*total_len_alleles;
    }
  else if (assump==AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice)
    {
      poisson_1net_err_rate= m*m*total_len_alleles;
    }

  double log_prob_error = 0;


  if (num_1net_errors>0)
    {
      log_prob_error += -poisson_1net_err_rate + 
	num_1net_errors*log(poisson_1net_err_rate) - gsl_sf_lnfact(num_1net_errors);
    }


  // 2.  Get sum of log likelihoods of all the subchunks of the two alleles


  double log_prob_data=0;
  double hap_D  = (double) (model_info->ginfo->total_sequence[colour_indiv])/(double) (2*model_info->genome_len) ;
  double hap_D_over_R = hap_D/(model_info->ginfo->mean_read_length[colour_indiv]);
  //int eff_r_plus_one = model_info->ginfo->mean_read_length[colour_indiv] - db_graph->kmer_size;


  //walk the first allele

  log_prob_data
    += get_log_probability_of_covg_on_one_allele_given_second_allele_and_multiplicities_using_little_hash(
         hap_D_over_R, var->one_allele, var->len_one_allele,
         var_mults->mult11, var_mults->mult12,
         working_array_self, working_array_shared,
         check_covg_in_ref_with_site_excised, colour_indiv);

  //now walk the second (which will be the same one again if this is hom)
  log_prob_data
    += get_log_probability_of_covg_on_one_allele_given_second_allele_and_multiplicities_using_little_hash(
         hap_D_over_R, var->other_allele, var->len_other_allele,
         var_mults->mult22, var_mults->mult21,
         working_array_self, working_array_shared,
         check_covg_in_ref_with_site_excised, colour_indiv);

  
  return log_prob_data + log_prob_error;

  
    
}



void dealloc_array_of_files(char** array_files, int num_files_in_list)
{
  int i;
  for (i=0; i<num_files_in_list; i++)
    {
      free(array_files[i]);
    }
  free(array_files);
}

//alloc and make a list of binaries, one per allele
char** alloc_array_and_get_files_from_list(char* filelist, int num_files_in_list)
{
  char** array_files=(char**) malloc(num_files_in_list*sizeof(char*));
  if (array_files==NULL)
    {
      die("Unable to malloc a filelist array! Either your machine is out-of-mem, "
          "or you have trying to use an obscene number of colours accidentally");
    }
  int i;
  for (i=0; i<num_files_in_list; i++)
    {
      array_files[i]=(char*)malloc(MAX_FILENAME_LENGTH*sizeof(char));
      if (array_files[i]==NULL)
	{
	  die("Unable to malloc the %d -th filename in an array, max filename length "
        "set to %d", i, MAX_FILENAME_LENGTH);
	}
      array_files[i][0]='\0';
    }

  //now fill the array
  FILE* fp = fopen(filelist, "r");
  if (fp==NULL)
    {
      die("Unable to open the filelist of 1net or2net binaries %s - surprised "
          "this did not get caught by checks of the commandline", filelist);
    }

  char filename1[MAX_FILENAME_LENGTH+1];
  filename1[0]='\0';

  for (i=0; i<num_files_in_list; i++)
    {
      if (fgets(filename1,MAX_FILENAME_LENGTH, fp) !=NULL)
	{
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(filename1, '\n')) != NULL)
	    {
	      *p = '\0';
	    }
	  strcpy(array_files[i], filename1);
	}
      else
	{
	  die("Expected %d files in this list:%s, but failed to read the %d-th one from it",
        num_files_in_list, filelist, i);
	}
    }
  fclose(fp);
  return array_files;
}


/* 
void wipe_colour_and_load_binaries(dBGraph* db_graph, int colour, char* bin1, char* bin2)
{
  db_graph_wipe_colour(colour, db_graph);
  int mean_readlen=0;
  long long seq=0;
  char dummy_sample_id[MAX_LEN_SAMPLE_NAME];
  dummy_sample_id[0]='\0';
  if (bin1 != NULL)
    {
      load_single_colour_binary_data_from_filename_into_graph(bin1, db_graph, &mean_readlen, &seq,
							      false, colour,
							      false, 0, dummy_sample_id, false);
    }
  if ( (bin2 !=NULL) && (strcmp(bin2, bin1) !=0) )
    {
      load_single_colour_binary_data_from_filename_into_graph(bin2, db_graph, &mean_readlen, &seq,
							      false, colour,
							      false, 0, dummy_sample_id, false);
    }
}

void wipe_two_colours_and_load_two_binaries(dBGraph* db_graph, int colour1, int colour2,
					    char* binary11, char* binary12, char* binary21, char* binary22)
{
  char dummy_sample_id[MAX_LEN_SAMPLE_NAME];
  dummy_sample_id[0]='\0';

  db_graph_wipe_two_colours_in_one_traversal(colour1, colour2, db_graph);
  int mean_readlen=0;
  long long seq=0;
  //normal use case - load two binaryies into th 1net colour,  for the 1net of each allele in the genotype
  if (binary11 !=NULL)
    {
      load_single_colour_binary_data_from_filename_into_graph(binary11, db_graph, &mean_readlen, &seq,
							      false, colour1,
							      false, 0,dummy_sample_id,  false);
    }
  if ( (binary12 != NULL) && (strcmp(binary12, binary11)!=0) )
    {
      load_single_colour_binary_data_from_filename_into_graph(binary12, db_graph, &mean_readlen, &seq,
							      false, colour1,
							      false, 0, dummy_sample_id, false);
    }
  //same for 2net
  if (binary21 != NULL)
    {
      load_single_colour_binary_data_from_filename_into_graph(binary21, db_graph, &mean_readlen, &seq,
							      false, colour2,
							      false, 0, dummy_sample_id, false);
    }
  if ( (binary22 != NULL) && (strcmp(binary22, binary21) !=0) )
    {
      load_single_colour_binary_data_from_filename_into_graph(binary22, db_graph, &mean_readlen, &seq,
							      false, colour2,
							      false, 0, dummy_sample_id, false);
    }


}
*/



void modify_character(char* str, int which_base, int which_mutant)
{
  int mutants[3];
  if (str[which_base]=='A')
    {
      mutants[0]='C';
      mutants[1]='G';
      mutants[2]='T';
    }
  else if (str[which_base]=='C')
    {
      mutants[0]='A';
      mutants[1]='G';
      mutants[2]='T';
    }
  else if (str[which_base]=='G')
    {
      mutants[0]='A';
      mutants[1]='C';
      mutants[2]='T';
    }
  else if (str[which_base]=='T')
    {
      mutants[0]='A';
      mutants[1]='C';
      mutants[2]='G';
    }
  else
    {
      die("Bug - str[which_base] is %c, whole string is %s, which_base is %d",
          str[which_base], str, which_base);
    }

  str[which_base]=mutants[which_mutant];

}

//idea is: first check if is in 1net, and if not, call this. It is now either in the 2net or beyond and this return from this tells you which
boolean is_this_kmer_beyond_the_2net(dBNode* n, dBGraph* db_graph, Covg (*get_covg_1net)(dBNode* e) )
{
  int i;

  char orig_seq[db_graph->kmer_size+1];
  binary_kmer_to_seq(&(n->kmer), db_graph->kmer_size, orig_seq);

  for (i=0; i<db_graph->kmer_size; i++)
    {
      //there are 3 possible kmers created by changing position i 
      int p;
      for (p=0; p<=2; p++)
	{
	  char copy_orig[db_graph->kmer_size+1];
	  strcpy(copy_orig, orig_seq);
	  modify_character(copy_orig, i,p);//modify the i-th character and change it to the p-th of the other 3. eg if the base is C, and p=2 change it to the second one of A,G,T.
	  BinaryKmer tmp_kmer1, tmp_kmer2;
	  dBNode* found_node = hash_table_find(element_get_key(seq_to_binary_kmer(copy_orig, db_graph->kmer_size, &tmp_kmer1), db_graph->kmer_size, &tmp_kmer2), db_graph);
	  if (found_node==NULL)
	    {
	      //this particular edit of 1 base does not take you back to the desired genotype or to the 1net, but maybe some other edit does
	    }
	  else if (db_node_check_status(found_node, in_desired_genotype)==true)
	    {
	      //ignore
	      //return false;
	    }
	  else if (get_covg_1net(found_node)>0)
	    {
	      return false;
	    }
	  else
	    {

	    }
	}

    }
  //all possible edits have failed to hit the genotype or 1net
  return true;
}


/*
// Currently noone appears to call this function
// Removed unused parameter Covg (*get_covg_1net)(dBNode* e)
boolean is_this_kmer_in_the_1net(dBNode* n, dBGraph* db_graph)
{
  int i;

  char orig_seq[db_graph->kmer_size+1];
  binary_kmer_to_seq(&(n->kmer), db_graph->kmer_size, orig_seq);

  for (i=0; i<db_graph->kmer_size; i++)
    {
      //there are 3 possible kmers created by changing position i 
      int p;
      for (p=0; p<=2; p++)
	{
	  char copy_orig[db_graph->kmer_size+1];
	  strcpy(copy_orig, orig_seq);
	  modify_character(copy_orig, i,p);//modify the i-th character and change it to the p-th of the other 3. eg if the base is C, and p=2 change it to the second one of A,G,T.
	  BinaryKmer tmp_kmer1, tmp_kmer2;
	  
	  seq_to_binary_kmer(copy_orig, db_graph->kmer_size, &tmp_kmer1);
	  element_get_key(tmp_kmer1, db_graph->kmer_size, &tmp_kmer2);
	  dBNode* found_node = hash_table_find(&tmp_kmer2, db_graph);

	  if (found_node==NULL)
	    {
	      //this particular edit of 1 base does not take you back to the desired genotype, but maybe some other edit does
	    }
	  else if (db_node_check_status(found_node, in_desired_genotype)==true)
	    {
	      return true;
	    }
	  else
	    {
	      //this particular edit of 1 base does not take you back to the desired genotype, but maybe some other edit does
	    }
	}

    }
  //all possible edits have failed to hit the genotype
  return false;
}
*/

 //we ASSUME colours 0 to number_alleles are the various alternate alleles, loading in multicolour_bin
void calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site(
  int* colours_to_genotype, int num_colours_to_genotype,
  int colour_ref_minus_site, int number_alleles,
  int first_gt, int last_gt, // of all the possible gt's
  int max_allele_length, char* fasta,//one read per allele
  AssumptionsOnGraphCleaning assump,
  double* current_max_lik_array, double* current_max_but_one_lik_array,
  char** name_current_max_lik_array, char** name_current_max_but_one_lik_array,
  boolean print_all_liks_calculated,//not just the top two
  GraphAndModelInfo* model_info, dBGraph* db_graph,
  int working_colour1, int working_colour2,
  boolean using_1net, boolean using_2net,
  double min_acceptable_llk)
//										      char* filelist_1net_binaries )
										      
{
  Covg get_covg_in_union_of_colours_to_genotype(dBNode* e)
  {
    return element_get_covg_for_colourlist(e, colours_to_genotype,
                                           num_colours_to_genotype);
  }

  //wipe the working colours:
  //db_graph_wipe_colour(working_colour1, db_graph);
  //db_graph_wipe_colour(working_colour2, db_graph);


  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    die("Out of memory trying to allocate Sequence");
  }
  alloc_sequence(seq,max_allele_length,MAX_READ_NAME_LEN);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      die("Failed to malloc kmer sliding window in "
          "db_graph_make_reference_path_based_sv_calls. Exit.");
    }
  
  
  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_allele_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      die("Failed to malloc kmer_window->kmer in "
          "db_graph_make_reference_path_based_sv_calls. Exit.");
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
      die("Cannot alloc arrays of arrays in "
          "print_log_liks_of_specified_set_of_genotypes_of_complex_site");
    }
  
  for (i=0; i<number_alleles; i++)
    {
      array_of_node_arrays[i] = (dBNode**) malloc(sizeof(dBNode*) * max_allele_length);
      array_of_or_arrays[i]   = (Orientation*) malloc( sizeof(Orientation) * max_allele_length);
      array_of_allele_names[i]= (char*) malloc(sizeof(char) * 200 );
      
      if ( (array_of_node_arrays[i]==NULL) || (array_of_or_arrays[i]==NULL) || (array_of_allele_names==NULL) )
	{
	  die("Cannot alloc the %d -th node and or array in "
        "print_log_liks_of_specified_set_of_genotypes_of_complex_site", i);
	}
      lengths_of_alleles[i]=0;
      array_of_allele_names[i][0]='\0';
    }


  //we also need to get supernodes to count up errors
  dBNode**     path_nodes        = (dBNode**) malloc(sizeof(dBNode*)*max_allele_length); 
  Orientation* path_orientations = (Orientation*) malloc(sizeof(Orientation)*max_allele_length); 
  Nucleotide*  path_labels       = (Nucleotide*) malloc(sizeof(Nucleotide)*max_allele_length);
  char*        path_string  = (char*) malloc(sizeof(char)*max_allele_length+1); //+1 for \0

  if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL) )
    {
      die("Cannot malloc arrays for db_graph_remove_errors_considering_covg_and_topology");
    }

  


  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv
    = alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(max_allele_length,
                                                        max_allele_length);
  
  // DEV: what are working_arrays ?
  //      Assuming temporary arrays to store coverages for calc genotype llk
  Covg* working_array_self = (Covg*) malloc(sizeof(Covg) * max_allele_length);
  Covg* working_array_shared = (Covg*) malloc(sizeof(Covg) * max_allele_length);
  
  if ( (mobv==NULL)||(working_array_self==NULL) || (working_array_shared==NULL))
    {
      die("Cannot alloc all the arrays in "
          "calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site. "
          "Give up and exit.");
    }
  


  /*
  char** array_files_1net_binaries=NULL;
  if (using_1net==true)
    {
      array_files_1net_binaries = alloc_array_and_get_files_from_list(filelist_1net_binaries, number_alleles);
    }
  */







  
  //end of initialisation 

  //create file reader
  int file_reader(FILE * fp, Sequence * seq, int max_allele_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      die("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_allele_length,new_entry,full_entry,offset);
    
    return ret;
  }


  //Get list of alleles with coverage 
  FILE* fp = fopen(fasta, "r");
  if (fp==NULL)
    {
      die("UNable to open %s. Exit", fasta);
    }


  for (j=0; j<= number_alleles-1; j++)
    {
      boolean f_entry=true;
      int num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_allele_length, 
								     array_of_node_arrays[j],//assumed alleles are colours 0..num_alleles-1 
								     array_of_or_arrays[j], 
								     false, &f_entry, file_reader, seq, kmer_window, db_graph, 0);
      if (!f_entry)
	{
	  die("One of these alleles is longer than the specified max_read_len");
	}
      strcat(array_of_allele_names[j], seq->name);
      lengths_of_alleles[j]=num_kmers;
    }
  fclose(fp);


      /*      
      //Do any of these nodes have covg in the colours/samples we want to genotype
      int p;
      boolean there_is_no_covg=true; //so far, as far as we know, there is no covg in any of our samples
      for (p=0; (p< num_kmers) && (there_is_no_covg==true); p++)
	{

	  if (get_covg_in_union_of_colours_to_genotype(array_of_node_arrays[j][p])>0)
	    {
		  there_is_no_covg=false;
		  break;	      
	    }
	
	//  int q;
	 // for (q=0; q<num_colours_to_genotype; q++)
	  //  {
	   //   if (db_node_get_coverage(array_of_node_arrays[j][p], colours_to_genotype[q])>0)  //<<<< MUCH MUCH FASTER TO DO THIS ONCE IN UNION OF COLOURS TO GENOTYPE
//		{
		 // there_is_no_covg=false;
		  //break;
	//	}
	 //   }
	
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
*/


  // For each pair in that list, keep incrementing a counter, and if the counter is betweenfirst_gt and last_gt
  // then go ahead and check it out

  int genotype_count=0;


  for (i=0; i<number_alleles; i++)
    {
      for (j=i; j<number_alleles; j++) // **** now that i and j are specified, we have specified a genotype
	{

	  //printf("ZAM genotype is i=%d, j=%d\n", i,j);
	  genotype_count++;

	  //printf("Start genotype number %d. i is %d and j is %d\n", genotype_count, i,j);
	  //fflush(stdout);
	  
	  if (genotype_count<first_gt)
	    {
	      continue;
	    }
	  if (genotype_count>last_gt)
	    {
	      continue;
	    }
	  
	  /*
	  //it's a bit messier having these if's and elses but it's more efficient to ensure you only traverse the hash once since we're inside
	  // a loop that will run many times
	  if ((using_1net==true)&& (using_2net==false))
	    {
	      wipe_colour_and_load_binaries(db_graph, working_colour_1net, array_files_1net_binaries[i], array_files_1net_binaries[j]);
	    }
	  else if ((using_2net==true)&&(using_1net==false))
	    {
	      die("Using 2net and not 1net? SHould have caight this earlier");
	      //wipe_colour_and_load_binaries(db_graph, working_colour_2net, array_files_2net_binaries[i], array_files_2net_binaries[j]);
	    }
	  else if ((using_1net==true)&& (using_2net==true) )
	    {
	      printf("Going to load these 2 binaries: %s and %s\n", array_files_1net_binaries[i], array_files_1net_binaries[j] );
	      wipe_two_colours_and_load_two_binaries(db_graph, working_colour_1net, working_colour_2net, 
						     array_files_1net_binaries[i], array_files_1net_binaries[j],
						     array_files_2net_binaries[i], array_files_2net_binaries[j]);
						     }*/
	

	  Covg get_covg_in_1net_errors_from_genotype(dBNode* n)
	  {
	    //return db_node_get_coverage(n, working_colour_1net);
      int arr[] = {i+number_alleles, j+number_alleles};
	    return element_get_covg_for_colourlist(n, arr, 2);
	  }

	  Covg get_covg_in_2net_errors_from_genotype(dBNode* n)
	  {
	    if (is_this_kmer_beyond_the_2net(n, db_graph, &get_covg_in_1net_errors_from_genotype)==false)
	      {
		return 1;
	      }
	    else
	      {
		return 0;
	      }
	  }
	  
	  
	  VariantBranchesAndFlanks var;
	  set_variant_branches_but_flanks_to_null(&var, 
						  array_of_node_arrays[i], 
						  array_of_or_arrays[i], 
						  lengths_of_alleles[i],
						  array_of_node_arrays[j], 
						  array_of_or_arrays[j],
						  lengths_of_alleles[j],
						  unknown);
	  
	  set_status_of_nodes_in_branches(&var, in_desired_genotype);
	  reset_MultiplicitiesAndOverlapsOfBiallelicVariant(mobv);
	  improved_initialise_multiplicities_of_allele_nodes_wrt_both_alleles(&var, mobv, working_colour1, working_colour2);
	  int z;
	  for (z=0; z<num_colours_to_genotype; z++)
	    {
	      //printf("Start next sample - this time z is %d", z);
	      //fflush(stdout);

	      char name[300];
	      if (strlen(array_of_allele_names[i]) + strlen(array_of_allele_names[j])>=300 )
		{
		  die("Names of alleles %s and %s are too long (%d) - concatenated, Cortex "
          "requires them to be less than 300 characters", 
			 array_of_allele_names[i], array_of_allele_names[j],(int)( strlen(array_of_allele_names[i]) + strlen(array_of_allele_names[j])) );
		}
	      sprintf(name, "%s/%s", array_of_allele_names[i], array_of_allele_names[j]); 
	      
	      double llk= calc_log_likelihood_of_genotype_with_complex_alleles(&var, name,
									       mobv, model_info, colours_to_genotype[z],
									       colour_ref_minus_site, db_graph, 
									       working_array_self, working_array_shared,
									       &(current_max_lik_array[z]), &(current_max_but_one_lik_array[z]),
									       name_current_max_lik_array[z], name_current_max_but_one_lik_array[z],
									       assump, path_nodes, path_orientations, path_labels, path_string, max_allele_length,
									       using_1net, get_covg_in_1net_errors_from_genotype, using_2net, get_covg_in_2net_errors_from_genotype, min_acceptable_llk);
	      
	      

	      if (print_all_liks_calculated==true)
		{
		  printf("Colour %d, GENOTYPE %s : LLK=%f\n", colours_to_genotype[z], name, llk);
		}
	      
	    }
	  set_status_of_nodes_in_branches(&var, none);
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
  free_sequence(&seq);
  free(kmer_window->kmer);
  free(kmer_window);
  for (i=0; i<number_alleles; i++)
    {
      free(array_of_node_arrays[i]);
      free(array_of_or_arrays[i]);
      free(array_of_allele_names[i]);
    }
  
  free(array_of_node_arrays);
  free(array_of_or_arrays);
  free(lengths_of_alleles);
  free(array_of_allele_names);
  free(working_array_self);
  free(working_array_shared);
  free(path_nodes);
  free(path_orientations);
  free(path_labels);
  free(path_string);
  /*
  if (using_1net==true)
    {
      dealloc_array_of_files(array_files_1net_binaries, num_colours_to_genotype);
    }
  */
  
  
}
  
double* alloc_ML_results_array(int num_samples_to_genotype)
{
  double* retarray = (double*) malloc(sizeof(double) * num_samples_to_genotype);
  if (retarray==NULL)
    {
      die("Unable to malloc %d colours in alloc_ML_results_array",
          num_samples_to_genotype);
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
      die("Ridiculous - unanle to malloc array of names in alloc_ML_results_names_array");
    }
  int i;
  for (i=0; i<num_samples_to_genotype; i++)
    {
      retarray[i] = (char*) malloc(sizeof(char) * MAX_READ_NAME_LEN);
      if (retarray[i]==NULL)
	{
	  die("Unable to malloc name array in alloc_ML_results_names_array");
	}
      else
	{
	  retarray[i][0]='\0';
	}
    }
  return retarray;
}


void wipe_little_graph(LittleHashTable* little_graph)
{
  void wipe_node(GenotypingElement* node)
  {
    /*
    int colour;
    for (colour=0; colour<MAX_ALLELES_SUPPORTED_FOR_STANDARD_GENOTYPING+NUMBER_OF_COLOURS+2 ; colour++)
      {
	node->individual_edges[colour]=0;
	node->coverage[colour]=0;
	node->status = (char)unassigned;
	}
    */
    db_genotyping_node_set_status(node, unassigned);
  }
  little_hash_table_traverse(&wipe_node, little_graph);
}



// this requires that we KNOW which is the ref allele
// since there are only 2 alleles, there is no chance of having a node in the graph with covg, but not in either
// of the two alleles being genotyperd. so no need for the 1net and 2net
// ASSUME you are being passed in zeroed arrays for working_array_self_and_shared
void calculate_llks_for_biallelic_site_using_full_model_for_one_colour_with_known_ref_allele(AnnotatedPutativeVariant* annovar,
											     AssumptionsOnGraphCleaning assump,
											     GraphAndModelInfo* model_info, 
											     LittleHashTable* little_db_graph,
											     dBGraph* db_graph, 
											     GenotypingElement** working_g_e_one, 
											     Orientation* working_o_one,
											     GenotypingElement** working_g_e_other, 
											     Orientation* working_o_other,
											     Covg* working_array_self,
											     Covg* working_array_shared,
											     int len_working_arrays,
											     MultiplicitiesAndOverlapsOfBiallelicVariant* mobv,
											     int colour_to_genotype, 
											     int working_colour1, int working_colour2,
											     dBNode** path_nodes, Orientation* path_orientations, 
											     Nucleotide* path_labels, char* path_string , int max_sup_len
											     )
{
  ExperimentType expt = model_info->expt_type;
  boolean first_allele_is_ref;
  int len_ref;

  if (annovar->var->which==unknown)
    {
      die("Cannot call calculate_llks_for_biallelic_site_using_full_model_for_one_colour_with_known_ref_allele "
          "unless the annovar specifies which allele is ref");
    }
  else if (annovar->var->which==first)//this should always be true for PD
    {
      first_allele_is_ref=true;
      len_ref=annovar->var->len_one_allele;
    }
  else
    {
      first_allele_is_ref=false;
      len_ref=annovar->var->len_other_allele;
    }

  //initialise the little graph - set all nodes to unassigned
  wipe_little_graph(little_db_graph);


  //load in the two alleles
  int i;
  for (i=0; i<annovar->var->len_one_allele; i++)
    {
      //the ref allele for PD may contain Ns
      if (annovar->var->one_allele[i]==NULL)
	{
	  continue;
	}
      //for non-N's, add to little hash
      boolean found=false;    
      GenotypingElement* ge = little_hash_table_find_or_insert(&(annovar->var->one_allele[i]->kmer), &found, little_db_graph);
      if (ge==NULL)
	{
	  die("Could neither find nor insert into the little hash - coding error, "
        "should be large enough");
	}
      if (found==false)
	{
	  //we have just inserted a new node into the little hash. Fill in its data from the main graph:
	  genotyping_element_initialise_from_normal_element(ge, annovar->var->one_allele[i], true);
	  //char tmpzam[db_graph->kmer_size+1];
	  //printf("Zam added %s to little hash\n", binary_kmer_to_seq(&(ge->kmer), db_graph->kmer_size, tmpzam));
	}
      else
	{
	  //is already in the graph, no need to update the edges or covg
	}
      working_g_e_one[i]=ge;
      working_o_one[i]=annovar->var->one_allele_or[i];

    }


  // I don't expect ever to see an N in the alt allele. But make robust to it and print a warning
  for (i=0; i<annovar->var->len_other_allele; i++)
    {
      if (annovar->var->other_allele[i]==NULL)
	{
	  warn("alt allele to be genotyped contains an N (or a kmer missing from the graph)-  unexpected, but will try "
         "to continue. Worth telling Zam (zam@well.ox.ac.uk).\n");
	  continue;
	}
      
      boolean found=false;
      GenotypingElement* ge = little_hash_table_find_or_insert(&(annovar->var->other_allele[i]->kmer), &found, little_db_graph);
      if (ge==NULL)
	{
	  die("Could neither find nor insert into the little hash - coding error, "
        "should be large enough");
	}
      if (found==false)
	{
	  //we have just inserted a new node into the little hash. Fill in its data from the main graph:
	  genotyping_element_initialise_from_normal_element(ge, annovar->var->other_allele[i], true);
	  //char tmpzam[db_graph->kmer_size+1];
	  //printf("Zam added %s to little hash\n", binary_kmer_to_seq(&(ge->kmer), db_graph->kmer_size, tmpzam));

	}
      else
	{
	  //is already in the graph, no need to update the edges or covg
	}
      working_g_e_other[i]=ge;
      working_o_other[i]=annovar->var->other_allele_or[i];

    }


  //now comes the clever bit. The Little Hash now contains only nodes corresponding to the two alleles.
  //the ref colour for each node contains the count for that node in the ENTIRE ref. 
  //so walk the REF allele, and decrement the covg in the ref colour once each time you hit it.
  // at the end, the ref colour has become the ref-minus-site colour :-)
  
  int ref_minus_site_col=model_info->ref_colour;
  for (i=0; i<=len_ref; i++)//zam changed this from < to <=
    {
      if (first_allele_is_ref)
	{
	  db_genotyping_node_decrement_coverage(working_g_e_one[i], ref_minus_site_col);
	}
      else
	{
	  db_genotyping_node_decrement_coverage(working_g_e_other[i], ref_minus_site_col);
	}
    }





  //  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv= alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(annovar->var->len_one_allele, 
  //												       annovar->var->len_other_allele);
  
  if (mobv==NULL)
    {
      die("Cannot alloc all the arrays in calculate_llks_for_biallelic_site_using_full_model_for_one_colour_with_known_ref_allele.  Give up and exit.");
    }
  


  //end of initialisation 


  //hom1
  GenotypingVariantBranchesAndFlanks var_test;
  set_genotyping_variant_branches_but_flanks_to_null(&var_test, 
						     working_g_e_one,
						     working_o_one,
						     annovar->var->len_one_allele,
						     working_g_e_one,
						     working_o_one,
						     annovar->var->len_one_allele,
						     unknown);//not going to use the WhichIsRef in here


  set_status_of_genotyping_nodes_in_branches(&var_test, in_desired_genotype);
  little_graph_wipe_colour(working_colour1, little_db_graph);
  little_graph_wipe_colour(working_colour2, little_db_graph);
  reset_MultiplicitiesAndOverlapsOfBiallelicVariant(mobv);
  improved_initialise_multiplicities_of_allele_genotyping_nodes_wrt_both_alleles(&var_test, mobv, working_colour1, working_colour2);

   annovar->gen_log_lh[colour_to_genotype].log_lh[hom_one] = 
     calc_log_likelihood_of_genotype_with_complex_alleles_using_little_hash(&var_test,
									    mobv, model_info, colour_to_genotype,
									    little_db_graph, db_graph, 
									    working_array_self, working_array_shared,
									    assump, 
									    path_nodes, path_orientations, path_labels, path_string, 
									    max_sup_len);

   set_status_of_genotyping_nodes_in_branches(&var_test, none);

   //memset just the amount you used back to zero
   memset(working_array_self, 0, 
   	  max_of_ints(annovar->var->len_one_allele, annovar->var->len_other_allele)*sizeof(int));
   memset(working_array_shared, 0, 
   	  max_of_ints(annovar->var->len_one_allele, annovar->var->len_other_allele)*sizeof(int));


   if ( (expt==EachColourADiploidSample) || (expt==EachColourADiploidSampleExceptTheRefColour) )
     {
       //het
       set_genotyping_variant_branches_but_flanks_to_null(&var_test, 
							  working_g_e_one,
							  working_o_one,
							  annovar->var->len_one_allele,
							  working_g_e_other,
							  working_o_other,
							  annovar->var->len_other_allele,
							  unknown);//not going to use the WhichIsRef in here


       set_status_of_genotyping_nodes_in_branches(&var_test, in_desired_genotype);
       little_graph_wipe_colour(working_colour1, little_db_graph);
       little_graph_wipe_colour(working_colour2, little_db_graph);
       reset_MultiplicitiesAndOverlapsOfBiallelicVariant(mobv);
       improved_initialise_multiplicities_of_allele_genotyping_nodes_wrt_both_alleles(&var_test, mobv, working_colour1, working_colour2);
       

       annovar->gen_log_lh[colour_to_genotype].log_lh[het]= 
	 calc_log_likelihood_of_genotype_with_complex_alleles_using_little_hash(&var_test,
										mobv, model_info, colour_to_genotype,
										little_db_graph, db_graph, 
										working_array_self, working_array_shared,
										assump, 
										path_nodes, path_orientations, path_labels, path_string, 
										max_sup_len);
       set_status_of_genotyping_nodes_in_branches(&var_test, none);
       memset(working_array_self, 0, 
	      max_of_ints(annovar->var->len_one_allele, annovar->var->len_other_allele)*sizeof(int));
       memset(working_array_shared, 0, 
   	  max_of_ints(annovar->var->len_one_allele, annovar->var->len_other_allele)*sizeof(int));




     }
   else//haploid
     {
       annovar->gen_log_lh[colour_to_genotype].log_lh[het]=-9999999;
     }

  //hom other

  set_genotyping_variant_branches_but_flanks_to_null(&var_test, 
						     working_g_e_other,
						     working_o_other,
						     annovar->var->len_other_allele,
						     working_g_e_other,
						     working_o_other,
						     annovar->var->len_other_allele,
						     unknown);//not going to use the WhichIsRef in here
	  
  set_status_of_genotyping_nodes_in_branches(&var_test, in_desired_genotype);
  little_graph_wipe_colour(working_colour1, little_db_graph);
  little_graph_wipe_colour(working_colour2, little_db_graph);
  reset_MultiplicitiesAndOverlapsOfBiallelicVariant(mobv);
  improved_initialise_multiplicities_of_allele_genotyping_nodes_wrt_both_alleles(&var_test, mobv, working_colour1, working_colour2);


  annovar->gen_log_lh[colour_to_genotype].log_lh[hom_other]= 
    calc_log_likelihood_of_genotype_with_complex_alleles_using_little_hash(&var_test,
									   mobv, model_info, colour_to_genotype,
									   little_db_graph, db_graph, 
									   working_array_self, working_array_shared,
									   assump, 
									   path_nodes, path_orientations, path_labels, path_string, 
									   max_sup_len);
  set_status_of_genotyping_nodes_in_branches(&var_test, none);

  memset(working_array_self, 0, 
	 max_of_ints(annovar->var->len_one_allele, annovar->var->len_other_allele)*sizeof(int));
  memset(working_array_shared, 0, 
	 max_of_ints(annovar->var->len_one_allele, annovar->var->len_other_allele)*sizeof(int));
  
   
  
}


void get_all_full_model_genotype_log_likelihoods_at_PD_call_for_one_colour(
  AnnotatedPutativeVariant* annovar, AssumptionsOnGraphCleaning assump,
  GraphAndModelInfo* model_info, LittleHashTable* little_db_graph, 
  dBGraph* db_graph, GenotypingWorkingPackage* gwp, int colour_to_genotype)
{
  GenotypingElement** working_g_e_one   = gwp->working_g_e_one;
  Orientation* working_o_one            = gwp->working_o_one;
  GenotypingElement** working_g_e_other = gwp->working_g_e_other;
  Orientation* working_o_other          = gwp->working_o_other;
  Covg* working_array_self              = gwp->working_array_self;
  Covg* working_array_shared            = gwp->working_array_shared;
  int len_working_arrays                = gwp->max_allele_len;

  // These are zeroed by the gwp allocation function, so provided we clean up at
  // the end of this call, everything is fine
  int working_colour1                   = gwp->working_colour1;
  int working_colour2                   = gwp->working_colour2;

  dBNode** path_nodes                   = gwp->path_nodes;
  Orientation* path_orientations        = gwp->path_orientations;
  Nucleotide* path_labels               = gwp->path_labels; 
  char* path_string                     = gwp->path_string;
  int max_sup_len                       = gwp->max_sup_len;

  MultiplicitiesAndOverlapsOfBiallelicVariant* mobv=gwp->mobv;
    
  //  if ( (working_colour1<0) || (working_colour1>=NUMBER_OF_COLOURS+1) || 
  //     (working_colour2<0) || (working_colour2>=NUMBER_OF_COLOURS+1)  )
  //  {
  //    die("Calling get_all_full_model_genotype_log_likelihoods_at_PD_call_for_one_colour with bad working colours %d, %d\n", working_colour1, working_colour2);
  //  }
  if ( (annovar->var->len_one_allele> gwp->max_allele_len)
       ||
       (annovar->var->len_other_allele> gwp->max_allele_len)
       )
    {
      die("Trying to genotype a variant where one of the alleles (lengths %d, %d) is longer than the arrays we have malloced (length %d)",
	     annovar->var->len_one_allele, annovar->var->len_other_allele, gwp->max_allele_len);
    }
  calculate_llks_for_biallelic_site_using_full_model_for_one_colour_with_known_ref_allele(annovar, assump,
											  model_info, little_db_graph, db_graph,
											  working_g_e_one, working_o_one, 
											  working_g_e_other, working_o_other,
											  working_array_self, working_array_shared, len_working_arrays,
											  mobv, colour_to_genotype, working_colour1, working_colour2,
											  path_nodes, path_orientations, path_labels, path_string,
											  max_sup_len);
  
}




void get_all_genotype_log_likelihoods_at_non_SNP_PD_call_for_one_colour(AnnotatedPutativeVariant* annovar, double seq_error_rate_per_base, double sequencing_depth_of_coverage, int read_length, int colour)
{
  //boolean too_short = false;
  //next two lines are an expenseive way of setting too_short!!
  //count_reads_on_allele_in_specific_colour(annovar->var->one_allele, annovar->var->len_one_allele, colour, &too_short);
  //count_reads_on_allele_in_specific_colour(annovar->var->other_allele, annovar->var->len_other_allele, colour, &too_short);

  boolean too_short = (annovar->var->len_one_allele <= 1 ||
                       annovar->var->len_other_allele <= 1);

  if (too_short==true)
  {
    annovar->too_short=true;
    int j;
    //set all the log likelihoods to zero.
    for (j=0; j<NUMBER_OF_COLOURS; j++)
  	{
  	  initialise_genotype_log_likelihoods(&(annovar->gen_log_lh[j]));
  	}
  }
  else //assume caller has checked is not a SNP
  {
    //set to hom_other (the variant allele)
    annovar->gen_log_lh[colour].log_lh[hom_one]   =-999;
    annovar->gen_log_lh[colour].log_lh[het]       =-999;
    annovar->gen_log_lh[colour].log_lh[hom_other] = 0;
  }
}


//returns true if was able to initialise
// if you have no GraphInfo, enter NULL
//if you do not know what sequencing_error_per_base is, enter -1, will use default of 0.01
//if you do not know what genome length is, enter -1, will use default of 3 billion (human)
// NOTE THIS GENOTYPES THE SITE
// if you enter ref_colour=-1, there is no ref. Otherwise, the ref colour will be ignored for calculations like theta, where
// we are aggregating covg.
// AssumptionsOnGraphCleaning assump, and GenotypingWorkingPackage gwp is only used for PD call genotyping, otherwise you can pass a NULL
// db_graph and little_db_graph only used by PD caller genotyping, so again you can pass NULL there for bubble calls
boolean initialise_putative_variant(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info,
				    VariantBranchesAndFlanks* var, DiscoveryMethod caller,  int kmer,
				    AssumptionsOnGraphCleaning assump, GenotypingWorkingPackage* gwp,
				    dBGraph* db_graph, LittleHashTable* little_db_graph, CovgArray* working_ca,
				    boolean do_genotyping, boolean use_median)
				    
{
  GraphInfo* ginfo = NULL;
  //double seq_error_rate_per_base= -1;
  long long genome_length       = -1;
  ExperimentType expt           = Unspecified;
  int ref_colour                = -1;

  if (do_genotyping==true)
    {
      if (model_info !=NULL)
	{
	  ginfo                  = model_info->ginfo;
	  //seq_error_rate_per_base= model_info->seq_error_rate_per_base;
	  genome_length          = model_info->genome_len;
	  expt                   = model_info->expt_type;
	  ref_colour             = model_info->ref_colour;
	}
      else
	{
	  die("Abort - asking to genotype a putative site without specifying the model");
	}
    }

  int number_individals=NUMBER_OF_COLOURS;
  if ( (ref_colour !=-1) && ( (ref_colour<0) || (ref_colour>=NUMBER_OF_COLOURS) ) )
    {
      die("Called initialise_putative_variant with a ref_colour which is not -1, "
          "nor one of the colours this executable is compiled for. \n"
          "Coding error - this should have been caught earlier - call Zam");
    }
  if (ref_colour != -1)
    {
      number_individals--;
    }

  annovar->kmer = kmer;
  annovar->caller = caller;
  annovar->var=var;

  boolean flag1=false;
  boolean flag2=false;

  flag1=get_num_effective_reads_on_branch(annovar->br1_covg, var->one_allele, var->len_one_allele, use_median, working_ca, ginfo, kmer);
  flag2=get_num_effective_reads_on_branch(annovar->br2_covg, var->other_allele, var->len_other_allele, use_median, working_ca, ginfo, kmer);

  annovar->too_short = false;
  if ( (flag1==true)||(flag2==true) )
    {
      annovar->too_short = true;
    }
  
  if (var->len_one_allele < var->len_other_allele)
    {
      annovar->len_start=var->len_one_allele ;
    }
  else
    {
      annovar->len_start=var->len_other_allele;
    }

  if (annovar->too_short==true)
    {
      //deliberately not setting values to zero, for performance reasons.
      //this annovar should be discarded as soon as this function returns.
    }
  else//is a valid putative site
    {
      get_num_effective_reads_on_branch(annovar->theta1, var->one_allele, annovar->len_start, use_median, working_ca, ginfo, kmer);
      get_num_effective_reads_on_branch(annovar->theta2, var->other_allele, annovar->len_start, use_median, working_ca, ginfo, kmer);

      int i;
      annovar->BigTheta = 0;
      annovar->BigThetaStart = 0;
      
      for (i=0; i<NUMBER_OF_COLOURS; i++)

	{
	  if (i==ref_colour)
	    {
	      initialise_genotype_log_likelihoods(&(annovar->gen_log_lh[i]));
	      continue;
	    }
	  annovar->BigTheta      += annovar->br1_covg[i] + annovar->br2_covg[i];
	  annovar->BigThetaStart += annovar->theta1[i] + annovar->theta2[i];
	  initialise_genotype_log_likelihoods(&(annovar->gen_log_lh[i]));
	}
      

      if (do_genotyping==true)
	{
	  for (i=0; i<NUMBER_OF_COLOURS; i++)
	    {
	      if (i==ref_colour)
		{
		  continue;
		}
	      double sequencing_depth_of_coverage=0;
	      //if (seq_error_rate_per_base==-1)
	      //	{
	      //  seq_error_rate_per_base=0.01;//default
	      //}
	      if (genome_length==-1)
		{
		  warn("Genome length not specified, so assuming is human (3,000,000,000bp)\n");
		  genome_length=3000000000;//default
		}
	      int mean_read_len=100;
	      if (ginfo !=NULL)
		{
		  sequencing_depth_of_coverage=(double) ginfo->total_sequence[i]/genome_length;
		  mean_read_len = ginfo->mean_read_length[i];
		}
	      if ( (sequencing_depth_of_coverage==0)|| (annovar->too_short==true) )
		{
		  annovar->genotype[i]=absent;
		}
	      else if ( (expt==EachColourADiploidSample) || (expt==EachColourADiploidSampleExceptTheRefColour) )
		{
		  
		  //the full genotyping model reduces to a simple formula for bubbles
		  //do standard genotyping according to our model for bubble calls and SNP calls from PD
		  if ( (caller==BubbleCaller) 
		       //|| 
		       //( (caller==SimplePathDivergenceCaller) && (annovar->var->len_one_allele==annovar->var->len_other_allele) && (annovar->var->len_one_allele<=annovar->kmer+1)    ) 
		       )
		    {
		      get_all_genotype_log_likelihoods_at_bubble_call_for_one_colour(annovar,  
										     model_info->ginfo->seq_err[i],
										     sequencing_depth_of_coverage,
										     mean_read_len,i);
		    }
		  else
		    {
		      get_all_full_model_genotype_log_likelihoods_at_PD_call_for_one_colour(annovar, assump,
											    model_info, little_db_graph, db_graph,gwp,i);
		    }
		  
		  if (annovar->gen_log_lh[i].log_lh[hom_one]>= annovar->gen_log_lh[i].log_lh[het])
		    {
		      if (annovar->gen_log_lh[i].log_lh[hom_one]>=annovar->gen_log_lh[i].log_lh[hom_other])
			{
			  annovar->genotype[i]=hom_one;
			}
		      else
			{
			  annovar->genotype[i]=hom_other;
			}
		    }
		  else if (annovar->gen_log_lh[i].log_lh[het]>=annovar->gen_log_lh[i].log_lh[hom_other])
		    {
		      annovar->genotype[i]=het;
		    }
		  else
		    {
		      annovar->genotype[i]=hom_other;
		    }
		  
		  
		}
	      else if ( (expt==EachColourAHaploidSample) || (expt==EachColourAHaploidSampleExceptTheRefColour) )
		{
		  
		  if (caller==BubbleCaller)
		    {
		      get_all_haploid_genotype_log_likelihoods_at_bubble_call_for_one_colour(annovar,  
											     model_info->ginfo->seq_err[i], 
											     sequencing_depth_of_coverage,
											     mean_read_len,
											     i);
		    }
		  else
		    {
		      get_all_full_model_genotype_log_likelihoods_at_PD_call_for_one_colour(annovar, 
											    assump,
											    model_info, 
											    little_db_graph, 
											    db_graph,
											    gwp,
											    i);

		    }
		  
		  
		  if (annovar->gen_log_lh[i].log_lh[hom_one]> annovar->gen_log_lh[i].log_lh[hom_other])
		    {
		      annovar->genotype[i]=hom_one;
		    }
		  else
		    {
		      annovar->genotype[i]=hom_other;
		    }
		  
		}
	      else
		{
		  annovar->genotype[i]=hom_one;
		}
	      
	      

	    }
	}
    }

  return true;
}
