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
  model_selection.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <maths.h>

#include <gsl_sf_gamma.h>

#include "element.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "cmd_line.h"
#include "graph_info.h"
#include "db_differentiation.h"
#include "math.h"
#include "db_variants.h"



/*
void set_BF_var_over_rep(ModelLogLikelihoodsAndBayesFactors* stats, double val)
{
  stats->log_bayes_factor_var_over_rep = val;
}
*/
boolean basic_model_selection_condition(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info)
{
  //demand variant is 100 times as likely as repeat
  int lthresh = log(100);
  //double allele_balance_prior = 5;

  if (annovar->len_start<=2)
    {
      //one of these branches is just too short. Fail this
      annovar->model_llks.llk_var=12345;
      annovar->model_llks.llk_rep=12345;
      return false;
    }

  //calling these functions sets the values of the likelihoods inside the annovar object/struct
  
  // Warning: commented out the next 2 lines when I added one seq error rate per
  // colour. uncomment them when you fix those functions to support this.
  // Also delete the line (void)model_info;
  //calculate_integrated_loglikelihood_of_snp_model_given_data(annovar, model_info);
  //calculate_integrated_loglikelihood_of_repeat_model_given_data(annovar, model_info, allele_balance_prior);
  
  // Let the compiler know that we know that we're not using model_info param
  (void)model_info;

  double log_bf_var_over_rep = annovar->model_llks.llk_var - annovar->model_llks.llk_rep;

  if  (log_bf_var_over_rep> lthresh) 
    {
      //set_BF_var_over_rep(stats, log_bf_var_over_rep);
      return true;
    }
  else
    {
      return false;
    }
}


/*
double calculate_integrated_loglikelihood_of_snp_model_given_data(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info)
{
   double step = (0.98-0.02)/NUM_STEPS; //grid/mesh/step  value for frequency
  
  //
  //sort out priors for frequency and genotype
  //
  double prior_freq[NUM_STEPS+1];
  int i;
  int j;
  double fsum=0;
  for (i=0; i<NUM_STEPS+1; i++)
    {
      double f = 0.02 + i*step;
      prior_freq[i] = 1/(f*(1-f));
      fsum += prior_freq[i];
    }
  for (i=0; i<NUM_STEPS+1; i++)
    {
      prior_freq[i] = prior_freq[i] /fsum;
    }

  double prior_gt[3][NUM_STEPS+1];//TODO - pass in malloced??
  for (i=0; i<NUM_STEPS+1; i++)
    {
      double f = 0.02 + i*step;
      prior_gt[0][i]=(1-f)*(1-f);
      prior_gt[1][i]=2*f*(1-f);
      prior_gt[2][i]=f*f;
    }
    
  //int mean_observed_cov = (annovar->BigThetaStart)/2;
  int total_covg_br1 = 0;
  int total_covg_br2 = 0;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (i==model_info->ref_colour)
	{
	  continue;
	}
      total_covg_br1 +=annovar->theta1[i];
      total_covg_br2 +=annovar->theta2[i];
    }

  //
  //Expected covg per allele on each sample
  //
  //len is the number of nodes on this branch we are cconsidering as informative
  double len = (double) annovar->len_start+ 1 - 2; // IMPORTANT   - +1 because number of nodes is var->len_one/other_allele +1, and subtract 2, so you ignore the first and last nodes
  //double len = (double) annovar->len_start - 2; // IMPORTANT   -  subtract 2, so you ignore the first and last nodes

  if (len<=0)
    {
      die("Should not call calculate_integrated_loglikelihood_of_snp_model_given_data when len_start <=2. Coding error somewhere");
    }

  double  lambda_s[NUMBER_OF_COLOURS];//expected Poisson parameter for depth of covg, for each sample, haploid (ie per allele)
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      //lambda_s = depth of covg * len /2R 
      double depth = model_info->ginfo->total_sequence[i] / model_info->genome_len;
      int  read_len = model_info->ginfo->mean_read_length[i];
      if ( (read_len<=0) && (i != model_info->ref_colour ) )
	{
	  die("Exiting - you seem to have read length %d in colour %d", read_len, i);
	}
      lambda_s[i]=  (depth * len /(2*read_len)) ;

    }

  
  //genotypes
  // int    gs[3] = {2,1,0};
  //expected counts for each genotype, allowing for errors
  //(can think of the following as this vector sum gs*(1-err)+(2-gs)*err)
  // es is (2(1-err) , 1, err)  <<< this is how many copies of allele1 you expect to see if is hom1, het, hom2
  // 2-es is ( 2err,1,2-err)    <<< this is how many copies of allele2 you expect to see if is hom1, het, hom2
  double es[3] = {2*(1-model_info->seq_error_rate_per_base) , 1, 2*model_info->seq_error_rate_per_base };


  //log likelihoods for each allele per sample, then combine a per sample per GT likelihood 
  double llk_samples_b1[3][NUMBER_OF_COLOURS];
  double llk_samples_b2[3][NUMBER_OF_COLOURS];
  double likelihoods[3][NUMBER_OF_COLOURS];//remember to ignore ref_colour
  double max[NUMBER_OF_COLOURS]; 
  double maxsum=0;

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      max[i]=1;
      if (i==model_info->ref_colour)
	{
	  //ignore
	  int j;
	  for (j=0; j<3; j++)
	    {
	      llk_samples_b1[j][i]=0;
	      llk_samples_b2[j][i]=0;
	    }
	}
      else
	{
	  for (j=0; j<3; j++)
	    {
	      //read count on allele 1 annovar->theta1[i]
	      //expected depth of covg on one copy is:lambda_s[i]
	      //genotype tells you how many copies, so you have 
	      // on br1 --> Poisson rate es[j]*lambda_s[i] taking value annovar->theta1[i]
	      // on br2 ->  Poisson rate (2-es[j])*lambda_s[i] taking value annovar->theta2[i]

	      llk_samples_b1[j][i] = -es[j]*lambda_s[i] + annovar->theta1[i] * log(es[j]*lambda_s[i]) - gsl_sf_lnfact(annovar->theta1[i]);
	      llk_samples_b2[j][i] = - (2-es[j])*lambda_s[i] + annovar->theta2[i] * log((2-es[j])*lambda_s[i]) - gsl_sf_lnfact(annovar->theta2[i]);
	      
	      
	      if ( (llk_samples_b1[j][i]+ llk_samples_b2[j][i] > max[i]) || (max[i]>0) )
		{
		  max[i] = llk_samples_b1[j][i] + llk_samples_b2[j][i];
		}
	    }
	}
    }

  // will offset by max
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      for (j=0; j<3; j++)
	{

	  likelihoods[j][i]=exp(llk_samples_b1[j][i] + llk_samples_b2[j][i] - max[i]);

	}
      maxsum +=max[i];
    }


  //log likelihoods for each frequency
  double llk_freq[NUM_STEPS+1];
  llk_freq[0]=0;
  int f;
  double maxf=1;
  for (f=1; f<NUM_STEPS+1; f++)
    {
      llk_freq[f]=0;
      for (i=0; i<NUMBER_OF_COLOURS; i++) //sum over samples
	{
	  if (i==model_info->ref_colour)
	    {
	      continue;
	    }
	  double sum=0;
	  for (j=0; j<3; j++)
	    {
	      sum += prior_gt[j][f] * likelihoods[j][i];  
	    }
	  llk_freq[f] += log(sum);
	}
      llk_freq[f] +=maxsum;
      if ( (llk_freq[f] > maxf) || (maxf>0) )
	{
	  maxf = llk_freq[f];
	}
    }
  
  //include prior and integrate
  
  double local_sum=0;
  for (f=1; f<NUM_STEPS+1; f++)
    {
      local_sum += prior_freq[f] * exp(llk_freq[f]-maxf); 
    }
  double llk = maxf + log(local_sum);

  //set the value in the annovar
  annovar->model_llks.llk_var=llk;
  //and return it also
  return llk;

}
*/



int factorial(int n)
{
  if (n<0)
    {
      die("Passing neg number %d into factorial", n);
    }
  if (n==0)
    {
      return 1;
    }
  else
    {
      return factorial(n-1);
    }
}

double beta_func(double alpha, double beta)
{
  double l1 = log_factorial(alpha);
  double l2 = log_factorial(beta);
  double l3 = log_factorial(alpha+beta);
  return exp(l1+l2-l3);
}

//bal_prior - we use a symetric Beta prior distribution on allele balance
// and this is that coefficient, so the prior is B(bal_prior,bal_prior)
double calculate_integrated_loglikelihood_of_repeat_model_given_data(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info, double bal_prior)
{

  int total_covg_br1 = 0;
  int total_covg_br2 = 0;
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (i==model_info->ref_colour)
	{
	  continue;
	}
      total_covg_br1 +=annovar->theta1[i];
      total_covg_br2 +=annovar->theta2[i];
    }
  

  //Log likelihood from allele balance

  double llk_bal= gsl_sf_lnbeta(bal_prior+total_covg_br1, bal_prior+total_covg_br2) - gsl_sf_lnbeta(bal_prior, bal_prior) ;
  //double llk_bal=beta_func(bal_prior+total_covg_br1, bal_prior+total_covg_br2)
  //                -beta_func(bal_prior, bal_prior);
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (i==model_info->ref_colour)
	{
	  continue;
	}
      //      llk_bal += log_factorial(annovar->theta1[i] + annovar->theta2[i]);
      //llk_bal -= log_factorial(annovar->theta1[i]);
      //llk_bal -= log_factorial(annovar->theta2[i]);
      llk_bal += gsl_sf_lnchoose(annovar->theta1[i] + annovar->theta2[i], annovar->theta1[i]);
    }


  //Log likelihood from coverage
  double prior_vec[MAX_EXPECTED_REPEAT_COPIES];
  double tmp_total=0;
  double mu = model_info->mu;
  for (i=0; i<MAX_EXPECTED_REPEAT_COPIES; i++)
    {
      prior_vec[i]=mu*pow(1-mu,i);
      tmp_total+=mu*pow(1-mu,i);
    }
  for (i=0; i<MAX_EXPECTED_REPEAT_COPIES; i++)
    {
      prior_vec[i]=prior_vec[i]/tmp_total;
    }

  double copy_number_arr[MAX_EXPECTED_REPEAT_COPIES+1];//index j corresponds to j copies
  double len = (double) annovar->len_start +1 - 2;//ignoring first and last nodes

  int r;
  int max=1;

  //int mean_observed_cov = (annovar->BigThetaStart)/2;
  double lambda_per_copy[NUMBER_OF_COLOURS];
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (i==model_info->ref_colour)
	{
	  lambda_per_copy[i]=0; //will ignore this
	}
      else
	{
	  double depth = model_info->ginfo->total_sequence[i] / model_info->genome_len;
	  int  read_len = model_info->ginfo->mean_read_length[i];
	  lambda_per_copy[i]=  (depth * len /(2*read_len)) ;
	}

    }


  for (r=1; r<=MAX_EXPECTED_REPEAT_COPIES; r++)
    {
      copy_number_arr[r]=0;//initialise

      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  if (i==model_info->ref_colour)
	    {
	      continue;
	    }
	  double pois_rate = lambda_per_copy[i]*2*r; // r copies on each haplotype

	  //add log probability of seeing what we see for idividual i, assuming r copies
	  copy_number_arr[r] += (annovar->theta1[i]+annovar->theta2[i])*log(pois_rate) - gsl_sf_lnfact(annovar->theta1[i]+annovar->theta2[i])  -pois_rate;
	}
      if (max>0)//will only happen the first time
	{
	  max = copy_number_arr[r];
	  if (r>1)
	    {
	      die("Should never happen - had positive maximum after the first step\n");
	    }
	}
      if (copy_number_arr[r]>max)
	{
	  max=copy_number_arr[r];
	}
    }
  
  double sum = 0;
  for (r=1; r<=MAX_EXPECTED_REPEAT_COPIES; r++)
    {
      sum += prior_vec[r-1] * exp(copy_number_arr[r]-max);
    }

  double llk_cov = log(sum)+max;

  //set the value inside the annovar
  annovar->model_llks.llk_rep = llk_bal + llk_cov;
  //but also return it
  return llk_bal + llk_cov;
  

  
}


// Can this be removed?
double get_log_bayesfactor_varmodel_over_repeatmodel(AnnotatedPutativeVariant* annovar,
                                                     GraphAndModelInfo* model_info)
{
  // These parameters are currently being ignored
  (void)annovar;
  (void)model_info;

  //  double allele_balance_prior=5;//goes into symmetric Beta
  return 0; //ZAM - commented this out when I moved to one seq error rate per colour,
  // as I have not modified calculate_integrated_loglikelihood_of_snp_model_given_data
  // to support that commenting out as I know I have to reimplement this anyway


  //return calculate_integrated_loglikelihood_of_snp_model_given_data(annovar, model_info)
  //- calculate_integrated_loglikelihood_of_repeat_model_given_data(annovar, model_info, allele_balance_prior);
}


//mu is the paramter of the geometric distribution describing the prior
//for number of copies in a repeat. lies between 0 and 1.
double ignore_get_log_bayesfactor_varmodel_over_repeatmodel(AnnotatedPutativeVariant* annovar, GraphAndModelInfo* model_info)
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
      //lambda_s = c_s * len /2R - c_s is DEPTH of covg  (i.e.  total_seq/genome_len)
      long long depth = model_info->ginfo->total_sequence[i] / model_info->genome_len;
      int  read_len = model_info->ginfo->mean_read_length[i];
      lambda_s[i]=  (depth * len /(2*read_len)) ;
      total_covg_across_samples += depth;
    }
  long long total_depth=total_covg_across_samples;
  long long  lambda =( total_depth * len /(  2*get_mean_readlen_across_colours(model_info->ginfo) ));
  double alpha = (double) lambda -log(1-mu);
  
  int theta1_bar = 0;
  int theta2_bar=0;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      theta1_bar += annovar->theta1[i];
      theta2_bar += annovar->theta2[i];
    }

  double log_prob_repeat = log(mu) + gsl_sf_lnfact(theta1_bar) + gsl_sf_lnfact(theta2_bar)
    - log(annovar->BigThetaStart + 1) - log(alpha) ; //(annovar->BigThetaStart + 1)*log(alpha);

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      int thet = annovar->theta1[i]+ annovar->theta2[i];
      log_prob_repeat +=  thet * log(lambda_s[i]/alpha) - gsl_sf_lnfact(annovar->theta1[i]) - gsl_sf_lnfact(annovar->theta2[i]);
    }
  

  printf("\nrepeat component is %f\n", log_prob_repeat);

  double log_prob_var = 0;

  /*
  // using exactly Gils model - covg as poisson, copy num2, plus allelic baance

  log_prob_var  = -2*lambda + (annovar->BigThetaStart) * log(2);

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      int thet = annovar->theta1[i]+ annovar->theta2[i];
      log_prob_var +=  thet * log(lambda_s[i]) - log_factorial(thet) ;
    }
  
  double varsum=0;
  double zeta = 0.128;
  double p;
  for (p=0.02; p<=0.98; p+= 0.02)
    {
      double prod=1;
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  double per_sample_sum=0;
	  if (annovar->theta1[i]==0)
	    {
	      per_sample_sum += p*p;
	    }
	  if (annovar->theta2[i]==0)
	    {
	      per_sample_sum += (1-p)*(1-p);
	    }
	  printf("Pow is %f, so add %f to the per sample sum \n", pow(0.5, annovar->theta1[i] + annovar->theta2[i]), 
		 2*p*(1-p)*pow(0.5, annovar->theta1[i] + annovar->theta2[i]));
	  per_sample_sum += 2*p*(1-p)*pow(0.5, annovar->theta1[i] + annovar->theta2[i]);
	  prod = prod*per_sample_sum;
	  printf("Prod becomes %f\n", prod);
	}
      varsum += prod*p/(1-p);
    }
  printf("Varsum is %f\n", varsum);
  log_prob_var += log(varsum*zeta);

  */


  // Gil's model, but choose Max Likelihood frequency, and condition on that, so not multiplying
  //   hundreds of small numbers


  log_prob_var  = -2*lambda + (annovar->BigThetaStart) * log(2);

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      int thet = annovar->theta1[i]+ annovar->theta2[i];
      log_prob_var +=  thet * log(lambda_s[i]) - log_factorial(thet) ;
    }
  
  //double varsum=0;
  //double zeta = 0.128;
  double p;
  double max_lik_p=0;
  double current_best=0;
  for (p=0.02; p<=0.98; p+= 0.02)
    {
      double prod=1;
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  double per_sample_sum=0;
	  if (annovar->theta1[i]==0)
	    {
	      per_sample_sum += p*p;
	    }
	  if (annovar->theta2[i]==0)
	    {
	      per_sample_sum += (1-p)*(1-p);
	    }
	  per_sample_sum += 2*p*(1-p)*pow(0.5, annovar->theta1[i] + annovar->theta2[i]);

	  prod = prod*per_sample_sum;
	  printf("Prod becomes %f\n", prod);
	}
      if (prod>current_best)
	{
	  current_best=prod;
	  max_lik_p=p;
	}
    }
  
  //then having chosen frequency, apply it
  double prod=1;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      double per_sample_sum=0;
      if (annovar->theta1[i]==0)
	{
	  per_sample_sum += max_lik_p*max_lik_p;
	}
      if (annovar->theta2[i]==0)
	{
	  per_sample_sum += (1-max_lik_p)*(1-max_lik_p);
	}
      per_sample_sum += 2*max_lik_p*(1-max_lik_p)*pow(0.5, annovar->theta1[i] + annovar->theta2[i]);
      prod=prod*per_sample_sum;
    }

  log_prob_var += log(prod);






  




  /*
  //this is the one we like
  
  double p;
  double prob_var=0;


  double zeta = 0.128; //1 over  integral 0.02 to 0.98 of 1/x(1-x)
  for (p=0.02 ; p<=0.98 ; p += 0.02)
    {
      int colour;
      double product = 1;
      for (colour=0; colour<NUMBER_OF_COLOURS; colour++)
	{
	  //	  printf("First %f\n", exp(annovar->gen_log_lh[colour].log_lh[hom_one]) );
	  // printf("Second %f\n", exp(annovar->gen_log_lh[colour].log_lh[het]) );
	  // printf("Third %f\n", exp(annovar->gen_log_lh[colour].log_lh[hom_other]) );
	  product = product * (      p*exp(annovar->gen_log_lh[colour].log_lh[hom_one])/(1-p) 
			       +     2*exp(annovar->gen_log_lh[colour].log_lh[het])
			       + (1-p)*exp(annovar->gen_log_lh[colour].log_lh[hom_other])/p );

	}

      prob_var += product*zeta;
    }
  printf("prob var is %.12f and log of it is %f\n", prob_var, log(prob_var));
  log_prob_var = log(prob_var); // if we want to multiply by permuations:  + log_factorial(NUMBER_OF_COLOURS);;

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
  



  /*

  //let's trust the genotyping, since we are assuming the var model. So work out the allele frequency of branch1 and branch2
 int max_lik_num_chroms_with_br1 = 0;
 int num_hets=0;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (annovar->genotype[i]==hom_one)
	{
	  max_lik_num_chroms_with_br1 +=2;
	}
      else if (annovar->genotype[i]==het)
	{
	  num_hets++;
	  max_lik_num_chroms_with_br1 +=1;
	}
    }

 
 if (max_lik_num_chroms_with_br1>0.98*2*NUMBER_OF_COLOURS)
   {
     max_lik_num_chroms_with_br1=(int) (0.98 * 2*NUMBER_OF_COLOURS);
   }


 double freq_br1 = ((double)max_lik_num_chroms_with_br1)/(2*NUMBER_OF_COLOURS);


 int colour;      
 double zeta = 0.128;

 
 double prob_of_this_allele_freq = zeta/(freq_br1 * (1-freq_br1));


 double prob_var_given_freq=1;
  for (colour=0; colour<NUMBER_OF_COLOURS; colour++)                         
    {    
  */
      /*
      double tmp;      
      if (annovar->genotype[colour]==het)
	{
	  tmp = 2*freq_br1*(1-freq_br1) ;
	}
      else if (annovar->genotype[colour]==hom_one)
	{
	  tmp  = freq_br1*freq_br1 ;
	}
      else 
	{
	  tmp = (1-freq_br1)*(1-freq_br1) ;
	}

      log_prob_var += log(tmp) + annovar->gen_log_lh[colour].log_lh[annovar->genotype[colour]];
      */
  /*
      prob_var_given_freq =prob_var_given_freq*
	((2*freq_br1*(1-freq_br1)) * exp(annovar->gen_log_lh[colour].log_lh[het])
	 +  freq_br1*freq_br1   * exp(annovar->gen_log_lh[colour].log_lh[hom_one])
	 +  (1-freq_br1)*(1-freq_br1)   * exp(annovar->gen_log_lh[colour].log_lh[hom_other]) ) ;
      
    }


  
  //printf("prob_var_given_freq is %f\n", prob_var_given_freq);
  log_prob_var += log(prob_of_this_allele_freq) + log(prob_var_given_freq) ;

  //printf("before Var element is %f\n", log_prob_var);
  //printf("log-factorial NUMBER_OF_COLOURS is %f\n", log_factorial(NUMBER_OF_COLOURS) );
  // log_prob_var += log_factorial(NUMBER_OF_COLOURS) ;
  //printf("1 after Var element is %f\n", log_prob_var);
  //printf("log_factorial(num_hets) is %f\n", log_factorial(num_hets) );
  //log_prob_var -= log_factorial(num_hets);
  //printf("2 after Var element is %f\n", log_prob_var);
  //printf("log_factorial(NUMBER_OF_COLOURS-num_hets) is %f\n", log_factorial(NUMBER_OF_COLOURS-num_hets));
  //log_prob_var -= log_factorial(NUMBER_OF_COLOURS-num_hets);
  //printf("3 after Var element is %f\n", log_prob_var);

  */




  /*
  //be blunt - look for excess hets only

  //just do a count o hets

 int max_lik_num_chroms_with_br1 = 0;
 int num_hets=0;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      if (annovar->genotype[i]==hom_one)
	{
	  max_lik_num_chroms_with_br1 +=2;
	}
      else if (annovar->genotype[i]==het)
	{
	  num_hets++;
	  max_lik_num_chroms_with_br1 +=1;
	}
    }

 double freq_br1 = ((double)max_lik_num_chroms_with_br1)/(2*NUMBER_OF_COLOURS);

 //binomial distribution for this many hets out of NUMBER_OF_COLOURS individuals, at this given allele frequency (chosen by Max Lik)

 //either use the MAx LIk freq only:
 // log_prob_var = log_factorial(NUMBER_OF_COLOURS) - log_factorial(NUMBER_OF_COLOURS-num_hets)
 //  -log_factorial(num_hets) 
 //  + num_hets*log(2*freq_br1*(1-freq_br1)) + (NUMBER_OF_COLOURS-num_hets)*log(1 - 2*freq_br1*(1-freq_br1) );

 //or sum over them all
 double zeta=0.128;
 log_prob_var = log_factorial(NUMBER_OF_COLOURS) - log_factorial(NUMBER_OF_COLOURS-num_hets)
   -log_factorial(num_hets);
 double p;
 for (p=0.02; p<=0.98; p+=0.02)
   {
     log_prob_var += num_hets*log(2*p*(1-p)) + (NUMBER_OF_COLOURS-num_hets)*log(1-2*p*(1-p)) -log(p)-log(1-p) + log(zeta);
   }
  */








 printf("Var bit is %f\n", log_prob_var);

  return log_prob_var - log_prob_repeat;



}
