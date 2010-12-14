/*
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
  db_graph.h 

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <open_hash/hash_table.h>
#include <stdio.h>

typedef HashTable dBGraph;

int db_graph_clip_tips(int threshold, dBGraph * db_graph);
int db_graph_clip_low_coverage_supernodes(int threshold, int, int, dBGraph * db_graph);

long long db_graph_health_check(boolean fix, dBGraph * db_graph);

int db_graph_clip_tip(int threshold, dBGraph * db_graph);

int db_graph_db_node_prune_edges_with_single_coverage(dBNode * node, 
						      int coverage,
						      void (*node_action)(dBNode * node), dBGraph * db_graph );

long long db_graph_prune_low_coverage_edges(int coverage_threshold, dBGraph * db_graph);

int db_graph_get_perfect_path(dBNode * node, Orientation orientation, int limit, void (*node_action)(dBNode * node),
			      dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
			      char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,
			      boolean * is_cycle, dBGraph * db_graph);





boolean db_graph_detect_bubble(dBNode * node,
			       Orientation orientation,
			       int limit, 
			       void (*node_action)(dBNode * node), 
			       int * length1,dBNode ** path_nodes1, Orientation * path_orientations1, Nucleotide * path_labels1,
			       char * seq1, double * avg_coverage1, int * min_coverage1, int * max_coverage1,
			       int * length2,dBNode ** path_nodes2, Orientation * path_orientations2, Nucleotide * path_labels2,
			       char * seq2, double * avg_coverage2, int * min_coverage2, int * max_coverage2,
			       dBGraph * db_graph);

boolean db_graph_db_node_smooth_bubble(dBNode * node, Orientation orientation, int limit,int coverage_limit,
				       void (*node_action)(dBNode * node),
				       dBGraph * db_graph);

boolean db_graph_db_node_prune_low_coverage(dBNode * node, int coverage,
					    void (*node_action)(dBNode * node),
					    dBGraph * db_graph);


// limit is the max length
// min_coverage, max_coverage and avg_coveragte refer to the internal nodes
int db_graph_supernode(dBNode * node,int limit,
		       void (*node_action)(dBNode * node), 
		       dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels,
		       char * seq, double * avg_coverage,int * min_coverage, int * max_coverage,boolean * is_cycle,
		       dBGraph * db_graph);


int db_graph_db_node_clip_tip(dBNode * node, int limit,
			      void (*node_action)(dBNode * node),
			      dBGraph * db_graph);



void db_graph_print_supernodes(char * filename,int max_length, boolean with_coverages, dBGraph * db_graph);


void db_graph_print_coverage(dBGraph * db_graph);

void db_graph_print_paths(char * filename, int max_length,boolean with_coverages, int coverage_threshold, dBGraph * db_graph);


long long db_graph_find_double_Ys(int max_length, dBGraph * db_graph);

int db_graph_remove_low_coverage_nodes(int coverage, dBGraph * db_graph);
void db_graph_dump_binary(char * filename, boolean (*condition)(dBNode * node), dBGraph * db_graph);



void db_graph_detect_vars_clean_bubbles(int max_length, dBGraph * db_graph);
void db_graph_detect_vars_dirty_bubbles(int max_length, dBGraph * db_graph);

void db_graph_traverse_with_array(void (*f)(HashTable*, Element *, int**, int),HashTable * hash_table, int** array, int length_of_array);
int db_graph_get_N50_of_supernodes(dBGraph* db_graph);


int int_cmp(const void *a, const void *b);


dBNode* db_graph_get_first_node_in_supernode_containing_given_node(dBNode* node,  dBGraph* db_graph);

dBNode* db_graph_get_next_node_in_supernode(dBNode* node, Orientation orientation, Orientation* next_orientation,  dBGraph* db_graph);

void db_graph_get_supernode_length_marking_it_as_visited(dBGraph* db_graph, Element* node, int** array_of_supernode_lengths, int length_of_array);

boolean db_graph_detect_X_with_orientation(dBNode * node, Orientation orientation, dBGraph * db_graph);

int db_graph_remove_bubbles(int limit, dBGraph * db_graph);

void db_graph_dump_hash_table(char * filename, dBGraph * db_graph);

long long db_graph_remove_weak_edges(int threshold,dBGraph * db_graph);
#endif /* DB_GRAPH_H_ */
