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
  model_selection.h
*/

#ifndef MODEL_SELECTION_H_
#define MODEL_SELECTION_H_

#include "global.h"
#include "model_info.h"
#include "graph_info.h"
#include "db_variants.h"
#include "experiment.h"

#define MAX_EXPECTED_REPEAT_COPIES 10
#define NUM_STEPS  48


boolean basic_model_selection_condition(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info);
double calculate_integrated_loglikelihood_of_repeat_model_given_data(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info, double allele_balance_prior_coeff);
double calculate_integrated_loglikelihood_of_snp_model_given_data(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info);

// Just an empty funciton:
double get_log_bayesfactor_varmodel_over_repeatmodel(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info);

#endif /* MODEL_SELECTION_H_ */
