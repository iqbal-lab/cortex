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
  db_variants.h
*/

#ifndef DB_VARIANTS_H_
#define DB_VARIANTS_H_

// system libraries                                                                                                                                                                      
#include <stdlib.h>
#include <math.h>

#include "global.h"
#include "graph_info.h"
#include "model_info.h"
#include "genotyping_element.h"
#include "element.h"
#include "experiment.h"
#include "dB_graph.h"
#include "experiment.h"
#include "dB_graph_supernode.h"

#define MAX_VARNAME_LEN 200

typedef struct {
  double llk_var;
  double llk_rep;
  double llk_err;
  //  double log_bayes_factor_var_over_rep;
  // double log_bayes_factor_var_over_error;
  // double log_bayes_factor_rep_over_error;
} ModelLogLikelihoodsAndBayesFactors;// this is for deciding if a site is variant or repeat or error

void initialise_stats(ModelLogLikelihoodsAndBayesFactors* stats);




//a variant can be hom for oen allele or the other, het or absent in an individual
typedef enum
  {
    hom_one = 0,
    het = 1,
    hom_other = 2,
    absent =3
  }zygosity;

typedef enum
 {
    full_flank_format   = 0,
    glf  = 1,
 } Variant_File_Format ;

typedef enum
  {
    allele_one = 0,
    allele_other = 1,
  } WhichAllele;

typedef enum
  {
    first = 0,
    second = 1,
    unknown = 2,
  } WhichAlleleIsRef ;

typedef struct{
  dBNode** flank5p;
  Orientation* flank5p_or;
  int len_flank5p;
  dBNode** one_allele; 
  Orientation* one_allele_or; 
  int len_one_allele;
  dBNode** other_allele; 
  Orientation* other_allele_or;
  int len_other_allele;
  dBNode** flank3p;    
  Orientation* flank3p_or;
  int len_flank3p;
  WhichAlleleIsRef which;
  char* var_name;//make sure this is malloced if you use it
  char* seq5p;
  char* seq_one;
  char* seq_other;
  char* seq3p;
  int num_alleles;
} VariantBranchesAndFlanks;


typedef struct{
  GenotypingElement** flank5p;
  Orientation* flank5p_or;
  int len_flank5p;
  GenotypingElement** one_allele; 
  Orientation* one_allele_or; 
  int len_one_allele;
  GenotypingElement** other_allele; 
  Orientation* other_allele_or;
  int len_other_allele;
  GenotypingElement** flank3p;    
  Orientation* flank3p_or;
  int len_flank3p;
  WhichAlleleIsRef which;
} GenotypingVariantBranchesAndFlanks;

typedef enum{
  BubbleCaller=0,
  SimplePathDivergenceCaller=1,
} DiscoveryMethod;

typedef struct {
  double log_lh[3];// hom_one, het and hom_other. Use zygosity as index to ensure this
} GenotypeLogLikelihoods;


// AnnotatedPutativeVariant contains a VariantBranchesAndFlanks object, plus info that we do not want to keep recalculating:
//covg on both alleles, and also, is sometimes useful only to look at the "breakpoint"/start - since
//one branch can be longer than the other, define the breakpoint to be the whole of the shorter
//branch, and the same length, of the longer branch. ie take min(len(br1), len(br2)) length of both branches
typedef struct {
  int kmer;
  boolean too_short; //if one or both alleles are too short (containing no interior nodes), this is set to true and a call is not made
  DiscoveryMethod caller;
  VariantBranchesAndFlanks* var;
  Covg br1_covg[NUMBER_OF_COLOURS];//covg on first interior kmer + read jumps
  Covg br2_covg[NUMBER_OF_COLOURS];
  int len_start;// = min of length(br1) and length(br2)
  Covg theta1[NUMBER_OF_COLOURS];//coverage on start of branch1 
  Covg theta2[NUMBER_OF_COLOURS];//coverage on start of branch2
  Covg BigTheta; // total coverage of all colours over both branches (sum of all elements of br1_covg and br2_covg)
  Covg BigThetaStart; // total coverage of all colours over start of
                           //  both branches (sum of all elements of theta1 and theta2)

  //under the model that this IS a variant, compare likelihoods of being hom-br1, het, hom-br2
  zygosity genotype[NUMBER_OF_COLOURS]; //for each colour, genotype calls. determine this by comparing Log Likelihoods of standard model as used for HLA etc
  GenotypeLogLikelihoods gen_log_lh[NUMBER_OF_COLOURS];
  ModelLogLikelihoodsAndBayesFactors model_llks;
} AnnotatedPutativeVariant;




//functions for VariantBranchesAndFlanks

VariantBranchesAndFlanks* alloc_VariantBranchesAndFlanks_object(
  int len_5p, int len_br1, int len_br2, int len_3p, int kmer_size);

void free_VariantBranchesAndFlanks_object(VariantBranchesAndFlanks* var);


void set_alloced_variant_branches_and_flanks_allowing_inputargs_in_either_order(VariantBranchesAndFlanks* var, 
										dBNode** flank5p,    Orientation* flank5p_or,    int len_flank5p,  Orientation arraydir_5p, 
										dBNode** one_allele, Orientation* one_allele_or, int len_one_allele, Orientation arraydir_one,
										dBNode** other_allele, Orientation* other_allele_or, int len_other_allele, Orientation arraydir_other,
										dBNode** flank3p,    Orientation* flank3p_or,    int len_flank3p, Orientation arraydir_3p,
										WhichAlleleIsRef which);

void set_variant_branches_and_flanks(VariantBranchesAndFlanks* var, 
				     dBNode** flank5p,    Orientation* flank5p_or,    int len_flank5p,
				     dBNode** one_allele, Orientation* one_allele_or, int len_one_allele, 
				     dBNode** other_allele, Orientation* other_allele_or, int len_other_allele, 
				     dBNode** flank3p,    Orientation* flank3p_or,    int len_flank3p, WhichAlleleIsRef which);

void set_variant_branches_but_flanks_to_null(VariantBranchesAndFlanks* var, 
					     dBNode** one_allele, Orientation* one_allele_or, int len_one_allele, 
					     dBNode** other_allele, Orientation* other_allele_or, int len_other_allele, 
					     WhichAlleleIsRef which);


void set_genotyping_variant_branches_and_flanks(GenotypingVariantBranchesAndFlanks* var, 
						GenotypingElement** flank5p,      Orientation* flank5p_or,      int len_flank5p,
						GenotypingElement** one_allele,   Orientation* one_allele_or,   int len_one_allele, 
						GenotypingElement** other_allele, Orientation* other_allele_or, int len_other_allele, 
						GenotypingElement** flank3p,      Orientation* flank3p_or,      int len_flank3p, WhichAlleleIsRef which);

void set_genotyping_variant_branches_but_flanks_to_null(GenotypingVariantBranchesAndFlanks* var, 
							GenotypingElement** one_allele, Orientation* one_allele_or, int len_one_allele, 
							GenotypingElement** other_allele, Orientation* other_allele_or, int len_other_allele, 
							WhichAlleleIsRef which);

void set_status_of_nodes_in_branches(VariantBranchesAndFlanks* var, NodeStatus status);
void set_status_of_genotyping_nodes_in_branches(GenotypingVariantBranchesAndFlanks* var, NodeStatus status);

void exact_copy_variant_branches_and_flanks(VariantBranchesAndFlanks copy_to,
                                            const VariantBranchesAndFlanks copy_from);

void copy_variant_branches_and_flanks_switching_branches(
  VariantBranchesAndFlanks copy_to, const VariantBranchesAndFlanks copy_from);

void action_set_flanks_and_branches_to_be_ignored(VariantBranchesAndFlanks* var);


void db_variant_action_do_nothing(VariantBranchesAndFlanks* var);
boolean db_variant_precisely_one_allele_is_in_given_func_of_colours(
  VariantBranchesAndFlanks* var, Edges (*get_colour)(const dBNode*),
  dBGraph* db_graph, WhichAllele* which);

//toy function - do not use
zygosity db_variant_get_zygosity_in_given_func_of_colours(
  VariantBranchesAndFlanks* var, Edges (*get_colour)(const dBNode*),
  dBGraph* db_graph);


//genotyping of a site known to be a variant
void initialise_genotype_log_likelihoods(GenotypeLogLikelihoods* gl);

void get_all_haploid_genotype_log_likelihoods_at_bubble_call_for_one_colour(
  AnnotatedPutativeVariant* annovar, double seq_error_rate_per_base,
  double sequencing_depth_of_coverage, int read_length, int colour);

void get_all_genotype_log_likelihoods_at_bubble_call_for_one_colour(
  AnnotatedPutativeVariant* annovar, double seq_error_rate_per_base,
  double sequencing_depth_of_coverage, int read_length, int colour);


double get_log_likelihood_of_genotype_on_variant_called_by_bubblecaller(
  zygosity genotype, double error_rate_per_base,
  Covg covg_branch_1, Covg covg_branch_2, 
  double theta_one, double theta_other); //int kmer was an unused param


long long get_big_theta(AnnotatedPutativeVariant* annovar);


//utility functions
boolean get_num_effective_reads_on_branch(Covg* array, dBNode** allele, int how_many_nodes, 
					  boolean use_median, CovgArray* working_ca, GraphInfo* ginfo, int kmer);

Covg count_reads_on_allele_in_specific_colour(dBNode** allele, int len, int colour, boolean* too_short);
Covg count_reads_on_allele_in_specific_colour_given_array_of_cvgs(Covg* covgs, int len, boolean* too_short);
Covg count_reads_on_allele_in_specific_func_of_colours(
  dBNode** allele, int len,
  Covg (*sum_of_covgs_in_desired_colours)(const Element *),
  boolean* too_short);

Covg median_covg_on_allele_in_specific_colour(dBNode** allele, int len, CovgArray* working_ca,
					      int colour, boolean* too_short);

Covg median_of_CovgArray(CovgArray* array, CovgArray* working_array);

#endif /* DB_VARIANTS_H_ */
