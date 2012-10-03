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
  model_info.c
*/

#include "model_info.h"

void initialise_model_info(GraphAndModelInfo* model_info, GraphInfo* ginfo, 
			   long long genome_len, double mu, //double seq_err_rate_per_base,
			   int ref_colour, int num_chroms, ExperimentType type, AssumptionsOnGraphCleaning assump)
{
  model_info->ginfo = ginfo;
  model_info->genome_len = genome_len;
  model_info->mu = mu;
  // model_info->seq_error_rate_per_base = seq_err_rate_per_base;
  model_info->ref_colour=ref_colour; //if -1, that means no ref colour
  model_info->num_haploid_chromosomes = num_chroms;
  model_info->expt_type=type;
  model_info->assump = assump;
}
