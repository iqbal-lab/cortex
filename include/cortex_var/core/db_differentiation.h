
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
#include <dB_graph.h>
#include <seq.h>
#include <limits.h>
#include <file_reader.h>
#include <db_variants.h>


void align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours(FileFormat format, char* list_of_fastaq, int max_read_length, int* array_of_colours, char** array_of_names_of_colours,
								      int num_of_colours, dBGraph* db_graph,int fastq_ascii_offset,
								      boolean is_for_testing, char** for_test_array_of_strings, int* for_test_index);
