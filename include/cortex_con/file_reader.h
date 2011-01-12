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
#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>
#include <seq.h>

int MAX_READ_LENGTH;

extern int MAX_FILENAME_LENGTH;


//these routines return the length of the read sequence, for the binary file is all the kmers conctenated

//for short fasta entries - reads or similar 

long long load_fasta_from_filename_into_graph(char* filename, long long * reads, long long * bad_reads, int max_chunk_length, long long start, long long end, dBGraph* db_graph);

//for fastq
long long load_fastq_from_filename_into_graph(char* filename, long long * reads, long long * bad_reads,  char quality_cut_off, int ascii_fq_offset,
					      int max_read_length, long long start, long long end, dBGraph* db_graph);

//for reference

void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
                                                                            long long * bad_reads, char quality_cut_off, int max_read_length, dBGraph * db_graph);

void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph);

void read_all_ref_chromosomes_and_mark_graph(dBGraph* db_graph);



// This would be in seq.c, but that has no knowledge of dBGraph
int get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(char * seq,  char * qualities, int length, char quality_cut_off,
                                                                                  KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph);

// for dumping clean fasta files from fastq - ie only holding reads that lie entirely in a (presumably cleaned) graph
void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
						  long long * bad_reads, int max_read_length, dBGraph * db_graph,
						  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);


//for binary
long long  load_binary_from_filename_into_graph(char* filename, int binary_version,dBGraph* db_graph, boolean all_entries_are_unique, boolean check_binary);

void print_binary_signature(FILE * fp,int kmer_size, int num_cols, int mean_read_len, long long total_seq);
boolean check_binary_signature(FILE * fp,int kmer_size, int bin_version, int* number_of_colours_in_binary,int* mean_read_len, long long* total_seq);




#endif /* FILE_READER_H_ */
