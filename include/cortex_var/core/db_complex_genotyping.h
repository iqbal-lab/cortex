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



#ifndef DB_COMPLEX_GENOTYPING_H_
#define DB_COMPLEX_GENOTYPING_H_

#include <element.h>
#include <db_variants.h>
#include <graph_info.h>
#include <model_selection.h>

typedef struct {
  int* mult11; //multiplicity of allele 1 nodes in allele 1
  int* mult12; //multiplicity of allele 1 nodes in allele 2,   etc
  int* mult21;
  int* mult22;
  int len1;//length of array of nodes for allele1 - ie mult11 and mult12
  int len2;

} MultiplicitiesAndOverlapsOfBiallelicVariant;

typedef enum {
  AssumeUncleaned                                = 1,
  AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice = 2,//because we cleaned off stuff that occurred only once
}AssumptionsOnGraphCleaning;


MultiplicitiesAndOverlapsOfBiallelicVariant* alloc_MultiplicitiesAndOverlapsOfBiallelicVariant(int len_allele1, int len_allele2);
void dealloc_MultiplicitiesAndOverlapsOfBiallelicVariant(MultiplicitiesAndOverlapsOfBiallelicVariant* mobv);
void reset_MultiplicitiesAndOverlapsOfBiallelicVariant(MultiplicitiesAndOverlapsOfBiallelicVariant* mobv);

void initialise_multiplicities_of_allele_nodes_wrt_both_alleles(VariantBranchesAndFlanks* var, MultiplicitiesAndOverlapsOfBiallelicVariant* mult,
								boolean only_count_nodes_with_edge_in_specified_colour_func,
								Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*) );

double calc_log_likelihood_of_genotype_with_complex_alleles(VariantBranchesAndFlanks* var,
							    char* name_of_this_genotype,
							    MultiplicitiesAndOverlapsOfBiallelicVariant* var_mults,
							    GraphAndModelInfo* model_info,
							    int colour_indiv, 
							    int colour_ref_minus_our_site, dBGraph* db_graph,
							    int* working_array_self, int* working_array_shared,
							    double* current_max_lik, double* current_max_but_one_lik,
							    char* current_max_lik_name, char* current_max_but_one_lik_name,
							    AssumptionsOnGraphCleaning assump);

//we ASSUME colours 0 to number_alleles are the various alternate alleles, loading in multicolour_bin
void calculate_max_and_max_but_one_llks_of_specified_set_of_genotypes_of_complex_site(int* colours_to_genotype, int num_colours_to_genotype,
										      int colour_ref_minus_site, int number_alleles,
										      int first_gt, int last_gt, // of the number_alleles choose 2
										      int max_allele_length, char* fasta,//one read per allele
										      AssumptionsOnGraphCleaning assump,
										      double* current_max_lik_array, double* current_max_but_one_lik_array,
										      char** name_current_max_lik_array, char** name_current_max_but_one_lik_array,
										      boolean print_all_liks_calculated, //not just the top two
										      GraphAndModelInfo* model_info, dBGraph* db_graph
										      );

double* alloc_ML_results_array(int num_samples_to_genotype);
char** alloc_ML_results_names_array(int num_samples_to_genotype);
#endif
