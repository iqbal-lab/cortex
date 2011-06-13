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


// model_selection.h

#ifndef MODEL_SELECTION_H_
#define MODEL_SELECTION_H_


#define MAX_EXPECTED_REPEAT_COPIES 10
#define NUM_STEPS  48

typedef struct {
  double log_bayes_factor_var_over_rep;
  double log_bayes_factor_var_over_error;
  double log_bayes_factor_rep_over_error;
} LogLikelihoodsAndBayesFactors;// this is for a single colour; only ever used locally during variant discovery, at a putativ variant site

typedef struct{
  GraphInfo* ginfo;
  long long genome_len;
  double mu; //parameter of geometric distirb describing prior for repeat copy number
  double seq_error_rate_per_base;
  int ref_colour;
} GraphAndModelInfo;

void initialise_stats(LogLikelihoodsAndBayesFactors* stats);
void initialise_model_info(GraphAndModelInfo* model_info, GraphInfo* ginfo, long long genome_len, double mu, double seq_err_rate_per_base);

void set_BF_var_over_rep(LogLikelihoodsAndBayesFactors* stats, double val);



boolean basic_model_selection_condition(AnnotatedPutativeVariant* annovar, LogLikelihoodsAndBayesFactors* stats, GraphAndModelInfo* model_info);
double get_log_bayesfactor_varmodel_over_repeatmodel(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info);
double calculate_integrated_loglikelihood_of_repeat_model_given_data(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info, double allele_balance_prior_coeff);
double calculate_integrated_loglikelihood_of_snp_model_given_data(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info);

#endif /* MODEL_SELECTION_H_ */
