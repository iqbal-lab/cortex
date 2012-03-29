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

#ifndef GRAPH_INFO_H_
#define GRAPH_INFO_H_

#include <stdint.h>
#include <global.h>
#include <stdio.h>

typedef struct{
  boolean tip_clipping;
  boolean remv_low_cov_sups;
  boolean remv_low_cov_nodes;
  int remv_low_cov_sups_thresh;
  int remv_low_cov_nodes_thresh;
}ErrorCleaning;

void error_cleaning_initialise(ErrorCleaning cl);


typedef struct{
  long long     total_sequence[NUMBER_OF_COLOURS];
  int           mean_read_length[NUMBER_OF_COLOURS];
  long double   seq_err[NUMBER_OF_COLOURS];
  ErrorCleaning cleaning[NUMBER_OF_COLOURS];
} GraphInfo;

void graph_info_initialise(GraphInfo* ginfo);

//set total amount of sequence in a colour
void graph_info_set_seq(GraphInfo* ginfo, int colour, long long num_bp);
//increment total amoutn of sequence in a colour (eg when you merge in a new binary)
long long graph_info_increment_seq(GraphInfo* ginfo, int colour, long long num_bp);
//set mean read length in a colour
void graph_info_set_mean_readlen(GraphInfo* ginfo, int colour, int len);
//update mean read length in a colour, eg when you merge a new binary
int graph_info_update_mean_readlen(GraphInfo* ginfo, int colour, int previous_mean, long long previous_seq, int mean_readlen_in_added_data, long long added_seq);
void graph_info_update_mean_readlen_and_total_seq(GraphInfo* ginfo, int colour,int mean_readlen_in_added_data, long long added_seq);
double get_total_depth_of_coverage_across_colours(GraphInfo* ginfo, long long genome_length);
int get_mean_readlen_across_colours(GraphInfo* ginfo);
void read_estimated_seq_errors_from_file(GraphInfo* ginfo, FILE* fp);
void print_seq_err_rates_to_screen(GraphInfo* ginfo);

#endif
