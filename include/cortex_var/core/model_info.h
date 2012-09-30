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
  model_info.h
*/

#ifndef MODEL_INFO_H_
#define MODEL_INFO_H_

#include "global.h"
#include "graph_info.h"
#include "experiment.h"

typedef struct{
  GraphInfo* ginfo;
  long long genome_len;
  double mu; //parameter of geometric distirb describing prior for repeat copy number
  //double seq_error_rate_per_base;
  int ref_colour;
  int num_haploid_chromosomes;
  ExperimentType expt_type;
  AssumptionsOnGraphCleaning assump;
} GraphAndModelInfo;

void initialise_model_info(GraphAndModelInfo* model_info, GraphInfo* ginfo, long long genome_len, double mu, //double seq_err_rate_per_base, 
			   int ref_colour, int num_chroms, ExperimentType expt_type, AssumptionsOnGraphCleaning assump);

#endif /* MODEL_INFO_H_ */
