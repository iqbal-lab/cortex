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
  
  //uint64_t total_covg_across_colours = get_total_coverage_across_colours(model_info->ginfo, model_info->genome_len);
  double mu = model_info->mu;

  //lambda = c_tot * len /2R
  double len = (double) annovar->len_start;
  //double lambda = ((double) total_covg_across_colours) * len /((double) ());
  long long  lambda_s[NUMBER_OF_COLOURS];
  long long  total_covg_across_samples=0;
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      //lambda_s = c_s * len /2R - c_s is DEPTH of covg  so total_seq/genome_len
      long long depth = model_info->ginfo->total_sequence[i] / model_info->genome_len;
      int  read_len = model_info->ginfo->mean_read_length[i];
      lambda_s[i]=  (depth * len /(2*read_len)) ;
      total_covg_across_samples += depth;
    }
  long long  lambda =( total_covg_across_samples * len /(  2*get_mean_readlen_across_colours(model_info->ginfo) ));
  double alpha = (double) lambda -log(1-mu);
  
  //double ret=log( alpha*(1-mu)/mu ) - 2*lambda + annovar->BigThetaStart * log(2*alpha) -log_factorial(annovar->BigThetaStart);

  double log_prob_repeat = log((1-mu)/mu);

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      int thet = annovar->theta1[i]+ annovar->theta2[i];
      log_prob_repeat +=  thet * log(lambda_s[i]) - log_factorial(thet);
    }
  
  log_prob_repeat += log_factorial(annovar->BigThetaStart);
  log_prob_repeat -= (annovar->BigThetaStart + 1) * log(alpha);
  
  //allow for alleilic ratio
  int theta1_bar = 0;
  int theta2_bar=0;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      theta1_bar += annovar->theta1[i];
      theta2_bar += annovar->theta2[i];
    }
  log_prob_repeat += log_factorial(theta1_bar);
  for (i=theta2_bar+1; i<=theta2_bar+theta1_bar+1; i++)
    {
      log_prob_repeat -= log(i);
    }
  printf("repeat component is %f\n", log_prob_repeat);

  double log_prob_var = 0;


  double p;
  double local_sum=0;

  /*
  double zeta = 0.128; //1 over  integral 0.02 to 0.98 of 1/x(1-x)
  for (p=0.02 ; p<=0.98 ; p += 0.02)
    {
      int colour;
      double product = 1;
      for (colour=0; colour<NUMBER_OF_COLOURS; colour++)
	{
	  product = product * (p*p*exp(annovar->gen_log_lh[colour].log_lh[hom_one]) + 2*p*(1-p)*exp(annovar->gen_log_lh[colour].log_lh[het])
			       + (1-p)*(1-p)*exp(annovar->gen_log_lh[colour].log_lh[hom_other]) );

	}
      printf("p is %f, add component %f\n", product*zeta/(p*(1-p)) );
      log_prob_var = log_prob_var + product*zeta/(p*(1-p));
    }
  */




  /*
  double K1=0.375;
  double K2=0.246;
  double K3=0.375;
  double zeta=0.128;
  */
  /*
  double sum=0;
  int colour;
  for (colour=0; colour<NUMBER_OF_COLOURS; colour++)
    {
      local_sum += K1 * exp(annovar->gen_log_lh[colour].log_lh[hom_one]) + K2 * exp(annovar->gen_log_lh[colour].log_lh[het]) + K3 * exp(annovar->gen_log_lh[colour].log_lh[hom_other]);
      sum = sum + log(local_sum);
    }
  */
  




  //let's trust the genotyping, since we are assuming the var model. So work out the allele frequency of branch1 and branch2
 int freq_br1 = 0;
 int freq_br2 = 0;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (annovar->genotype[i]==hom_one)
	{
	  freq_br1 +=2;
	}
      else if (annovar->genotype[i]==het)
	{
	  freq_br1 +=1;
	}
    }
 
 if (freq_br1>0.98*NUMBER_OF_COLOURS)
   {
     freq_br1=(int) (0.98 * NUMBER_OF_COLOURS);
   }


 double dec_freq_br1 = ((double)freq_br1)/NUMBER_OF_COLOURS;




  int colour;      
  double zeta = 0.128;

  double factor1 = zeta/(dec_freq_br1 * (1-dec_freq_br1));
  for (colour=0; colour<NUMBER_OF_COLOURS; colour++)                                                                                                                                            
    {    
      //log_prob_var += annovar->gen_log_lh[colour].log_lh[annovar->genotype[colour]];
      if (annovar->genotype[colour]==het)
	{
	  log_prob_var += log(2*dec_freq_br1*(1-dec_freq_br1) );
	}
      else if (annovar->genotype[colour]==hom_one)
	{
	  log_prob_var += log(dec_freq_br1*dec_freq_br1 );
	}
      else 
	{
	  log_prob_var += log((1-dec_freq_br1)*(1-dec_freq_br1) );
	}
    }
  log_prob_var += log_factorial(NUMBER_OF_COLOURS) - log_factorial(NUMBER_OF_COLOURS-freq_br1) - log_factorial(freq_br1);
  log_prob_var += log(factor1);




  printf("Var element is %f\n", log_prob_var);

  return log_prob_var - log_prob_repeat;
}
