
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
  internal_oxford.h
*/

#ifndef INTERNAL_OX_H_
#define INTERNAL_OX_H_

#include "global.h"

void set_ref_chromosome_file_pointers(char** reference_chromosome_file_ptrs, int num_chromosomes);
void create_uniqueness_file(char** ref_chroms, dBGraph* db_graph);
void print_coverage_and_ref_multiplicities_for_list_of_fasta(char* list_of_fasta, dBGraph* db_graph);

#endif /* INTERNAL_OX_H_ */
