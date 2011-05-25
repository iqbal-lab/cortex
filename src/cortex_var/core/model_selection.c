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

#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <string.h>
#include <cmd_line.h>
#include <time.h>
#include <graph_info.h>
#include <db_differentiation.h>
#include <maths.h>
#include <math.h> //we need both!

void initialise_stats(LogLikelihoodsAndBayesFactors* stats)
{
  stats->log_bayes_factor_var_over_rep = 0;
  stats->log_bayes_factor_var_over_error = 0;
  stats->log_bayes_factor_rep_over_error = 0;
}

void initialise_model_info(GraphAndModelInfo* model_info, GraphInfo* ginfo, long long genome_len, double mu, double seq_err_rate_per_base)
{
  model_info->ginfo = ginfo;
  model_info->genome_len = genome_len;
  model_info->mu = mu;
  model_info->seq_error_rate_per_base = seq_err_rate_per_base;
}

void set_BF_var_over_rep(LogLikelihoodsAndBayesFactors* stats, double val)
{
  stats->log_bayes_factor_var_over_rep = val;
}

boolean basic_model_selection_condition(AnnotatedPutativeVariant* annovar, LogLikelihoodsAndBayesFactors* stats, GraphAndModelInfo* model_info)
{

  double log_bf_var_over_rep = get_log_bayesfactor_varmodel_over_repeatmodel(annovar, model_info); 

  if  (log_bf_var_over_rep>0) 
    {
      set_BF_var_over_rep(stats, log_bf_var_over_rep);
      return true;
    }
  else
    {
      return false;
    }
}

//mu is the paramter of the geometric distribution describing the prior
//for number of copies in a repeat. lies between 0 and 1.
double get_log_bayesfactor_varmodel_over_repeatmodel(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info)
{
  
  
  uint64_t total_covg_across_colours = get_total_coverage_across_colours(model_info->ginfo, model_info->genome_len);
  double mu = model_info->mu;

  //lambda = c_tot * len /2R
  double lambda = ((double) total_covg_across_colours) * ((double) annovar->len_start) /((double) (2*get_mean_readlen_across_colours(model_info->ginfo)));

  double alpha = NUMBER_OF_COLOURS*lambda -log(1-mu);
  
  double ret=log( alpha*(1-mu)/mu ) - 2*lambda + annovar->BigThetaStart * log(2*alpha) -log_factorial(annovar->BigThetaStart);
  
  double p;
  double local_sum=0;

  double sfs_times_prior(zygosity genotype, double allele_freq)
  {
    double zeta = 0.128; //1 over  integral 0.02 to 0.98 of 1/x(1-x)
    if (genotype==hom_one)
      {
	return p* zeta /(1-p);  // p^2 * zeta /( p(1-p) )
      }
    else if (genotype==het)
      {
	return 2 * zeta; // 2*p*q * zeta/(p(1-p))
      }
    else if (genotype==hom_other)
      {
	return (1-p)* zeta /p;// (1-p)^2 * zeta / (p(1-p))
      }
    else
      {
	printf("Do not call sfs_times_prior with zygosity=absent arg1.\n");
	exit(1);
	return 0;
      }
  }
  
  double sum =0;
  for (p=0.02 ; p<=0.98 ; p += 0.02)
    {
      int colour;
      double product = 1;
      for (colour=0; colour<NUMBER_OF_COLOURS; colour++)
	{
	  product = product * sfs_times_prior(annovar->genotype[colour], p);
	}
      sum +=product;
    }
  
  ret +=log(sum);
  return ret;
}
