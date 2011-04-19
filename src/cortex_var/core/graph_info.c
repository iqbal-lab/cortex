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

void graph_info_initialise(GraphInfo* ginfo)
{
  int i;

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      graph_info_set_seq(ginfo, i, 0);
      graph_info_set_mean_readlen(ginfo, i, 0);
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
