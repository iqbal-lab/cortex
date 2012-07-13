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

#include <dB_graph.h>
#include <graph_info.h>
#include <file_reader.h>
#include <string.h>
#include <global.h>
#include <stdlib.h>

int MAX_LEN_SAMPLE_NAME=1000;

ErrorCleaning* error_cleaning_alloc_and_init()
{
  ErrorCleaning* ec = (ErrorCleaning*) malloc(sizeof(ErrorCleaning));
  if (ec==NULL)
    {
      printf("Cannot even alloc a tiny error-cleaning info object. Must be some OOM problem. Abort\n");
      exit(1);
    }
  memset(ec, 0, sizeof(ErrorCleaning));
  ec->name_of_graph_against_which_was_cleaned = (char*) malloc(sizeof(char)*MAX_FILENAME_LENGTH);
  if (ec->name_of_graph_against_which_was_cleaned==NULL)
    {
      printf("Cannot even alloc a tiny error-cleaning info object. Must be some OOM problem. Abort\n");
      exit(1);      
    }
  error_cleaning_initialise(ec);
  return ec;
}
void error_cleaning_free(ErrorCleaning* cl)
{
  free(cl->name_of_graph_against_which_was_cleaned);
  free(cl);
}

void error_cleaning_initialise(ErrorCleaning* cl)
{
  cl->tip_clipping=false;
  cl->remv_low_cov_sups=false;
  cl->remv_low_cov_nodes=false;
  cl->cleaned_against_another_graph=false;
  cl->remv_low_cov_sups_thresh=-1;
  cl->remv_low_cov_nodes_thresh=-1;
  cl->name_of_graph_against_which_was_cleaned[0]='\0';
  strcat(cl->name_of_graph_against_which_was_cleaned, "undefined");
  cl->len_name_of_graph_against_which_was_cleaned=
    strlen(cl->name_of_graph_against_which_was_cleaned);
}

GraphInfo* graph_info_alloc_and_init()
{
  GraphInfo* ginfo = (GraphInfo*) malloc(sizeof(GraphInfo));
  if (ginfo==NULL)
    {
      printf("Cannot even alloc a GraphInfo object. Your machine must be severely out of memory. Abort\n");
      exit(1);
    }
  else
    {
      memset(ginfo, 0, sizeof(GraphInfo));
      int i;
      for (i=0; i<NUMBER_OF_COLOURS; i++)
	{
	  ginfo->sample_ids[i] = (char*) malloc(sizeof(char)*MAX_LEN_SAMPLE_NAME) ;
	  if (ginfo->sample_ids[i]==NULL)
	    {
	      printf("Cannot even alloc a GraphInfo object. Your machine must be severely out of memory. Abort\n");
	      exit(1);
	    }
	  ginfo->cleaning[i]=error_cleaning_alloc_and_init(); //will abort if cannot alloc. Should never happen
	}
      
      graph_info_initialise(ginfo);
      return ginfo;
    }
  return NULL;
}

void graph_info_free(GraphInfo* ginfo)
{
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      free(ginfo->sample_ids[i]);
      error_cleaning_free(ginfo->cleaning[i]);
    }
  free(ginfo);
}

void graph_info_initialise(GraphInfo* ginfo)
{
  int i;

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      (ginfo->sample_ids[i])[0]='\0';
      strcat(ginfo->sample_ids[i], "undefined");
      ginfo->sample_id_lens[i]=strlen(ginfo->sample_ids[i]);

      graph_info_set_seq(ginfo, i, 0);
      graph_info_set_mean_readlen(ginfo, i, 0);
      ginfo->seq_err[i]=0.01;
      error_cleaning_initialise(ginfo->cleaning[i]);
    }
}



void graph_info_set_seq(GraphInfo* ginfo, int colour, long long num_bp)
{
  ginfo->total_sequence[colour]=num_bp;
}


long long graph_info_increment_seq(GraphInfo* ginfo, int colour, long long num_bp)
{
  ginfo->total_sequence[colour]+=num_bp;
  return ginfo->total_sequence[colour];
}

void graph_info_set_mean_readlen(GraphInfo* ginfo, int colour, int len)
{
  ginfo->mean_read_length[colour]=len;
}


//note if you are updating both mean read len AND total seq, 
// then do this one first (once you update the total seq, you no longer know 
//  what it used to be, ie the previous_seq)
int graph_info_update_mean_readlen(GraphInfo* ginfo, int colour, int previous_mean, long long previous_seq, 
			int mean_readlen_in_added_data, long long added_seq)
{
  if (added_seq==0)
    {
      return previous_mean;
    }
  else if (mean_readlen_in_added_data==0)
    {
      //printf("Warning - adding data with mean read-length 0\n");
      return previous_mean;
    }
  long long numerator = ((long long) previous_mean) * previous_seq + 
                         ((long long)mean_readlen_in_added_data) * added_seq;

  if (previous_seq+added_seq==0)
    {
      printf("WARNING - Updating graph_info object which contains no data with extra zero data!\n");
      return 0;
    }
  int new_mean = (int) (numerator/( previous_seq+added_seq));
  ginfo->mean_read_length[colour]=new_mean;
  return new_mean;
}

void graph_info_update_mean_readlen_and_total_seq(GraphInfo* ginfo, int colour,int mean_readlen_in_added_data, 
				       long long added_seq)
{
  graph_info_update_mean_readlen(ginfo, colour, 
		      ginfo->mean_read_length[colour], ginfo->total_sequence[colour], 
		      mean_readlen_in_added_data, added_seq);
  graph_info_increment_seq(ginfo, colour, added_seq);
}



double get_total_depth_of_coverage_across_colours(GraphInfo* ginfo, long long genome_length)
{
  int i;
  long long total_seq =0;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      total_seq += ginfo->total_sequence[i];
    }
  
  if (total_seq/genome_length<1)
    {
      printf("Warning - total sequence contained in union of all colours is less than 1x coverage of your genome. Is this really what you intend?\n");
    }

  return ((double) total_seq)/ ((double) genome_length);
}


int get_mean_readlen_across_colours(GraphInfo* ginfo)
{
  int colour;

  long long alpha=0;
  long long beta=0;
  for (colour=0; colour<NUMBER_OF_COLOURS; colour++)
    {
      alpha += (ginfo->total_sequence[colour]) * (ginfo->mean_read_length[colour]);
      beta += (ginfo->total_sequence[colour]);
    }
  return (int) (alpha/beta);
}

void read_estimated_seq_errors_from_file(GraphInfo* ginfo, FILE* fp)
{
  char line[MAX_FILENAME_LENGTH+1];
  int col=0;

  while (fgets(line, MAX_FILENAME_LENGTH, fp) !=NULL)
    {
            //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	*p = '\0';
      
      if (col<NUMBER_OF_COLOURS) //just for robustness - should always be true
	{
	  ginfo->seq_err[col]=  (long double) strtod(line ,NULL);
	  col++;
	}
	
    }
}

void print_seq_err_rates_to_screen(GraphInfo* ginfo)
{
  printf("Setting the following per-colour sequencing error rates (used only for genotyping):\nColour\tRate\n");
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      printf("%d\t%.3Lf\n", i, ginfo->seq_err[i]);
    }
}
